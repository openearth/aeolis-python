/* Viewer tab: the single place to manage everything shown on the map.
 *
 * - One collapsible section per layer group (Model output /
 *   Interpolated .grd / Sample data / Grid / Objects / Background),
 *   each layer in its own draggable card - same look as the Domain tab.
 * - Stale .grd files (grid changed) cannot be displayed: the eye is
 *   disabled and an orange warning icon explains why on hover.
 * - Every colormapped layer uses a saved colormap style (preset) or a
 *   custom one; the legend bottom-left lists the visible styles.
 * - netCDF output: one card per output variable (own eye + style), with
 *   instant time scrubbing (client LRU + prefetch + GPU interpolation).
 */
"use strict";

const ViewerTab = (() => {

  let meta = null;              // /api/output/meta
  let srcInfo = null;           // {override, path} — where the output is read from
  let mesh = null;              // {x, y, n, s}
  let frameCache = new Map();   // "var|extra|t" -> Float32Array (LRU)
  let frameRange = new Map();
  let inflight = new Map();
  let els = {};

  const FRAME_CACHE_MAX = 60;

  /* Every output variable is its own registry layer "output-<var>" with
   * field id "field-output-<var>" - shown/hidden, styled and reordered
   * exactly like any other layer card. Per-variable runtime state
   * (extra-dim indices, current frame bracket, auto-range flag). */
  const OUT_PREFIX = "output-";
  const _isOut = (id) => id.startsWith(OUT_PREFIX);
  const _outVar = (id) => id.slice(OUT_PREFIX.length);
  const outState = new Map();   // var -> {dims: [], bracket: null, auto: true}

  function _outStateOf(varName) {
    let st = outState.get(varName);
    if (!st) { st = { dims: [], bracket: null, auto: true }; outState.set(varName, st); }
    return st;
  }

  function _extraIdxOf(varName) {
    const dims = _outStateOf(varName).dims;
    return dims.length ? dims.join(",") : "0";
  }

  // which output variables are shown - persisted per project
  function _outVisibility() {
    return App.state.ui.outputVisible || (App.state.ui.outputVisible = {});
  }

  const GROUPS = [
    ["output", "Model output"],
    ["custom", "Custom layers"],
    ["domain", "Model files"],
    ["rawdata", "Sample data"],
    ["grid", "Grid"],
    ["transect", "Transects"],
    ["objects", "Objects"],
    ["background", "Background"],
  ];

  const loadingLayers = new Set();
  const dataRanges = new Map();   // layer id -> [lo, hi] of the loaded data
  const pointLayers = new Set();  // layer ids rendered as circle layers
  const pointData = new Map();    // layer id -> {x,y,z,radius} for hover/transect sampling
  const customIds = new Set();    // live field ids of custom computed layers

  function init() {
    Tabs.register("viewer", { enter: _enter });
    App.on("project", () => { _reset(); _loadOutput(); });
    App.on("run-finished", () => _loadOutput(true));
    // the Run tab pointed the output API at another run folder (or back)
    App.on("output-source", () => _loadOutput(true));
    App.on("clock-tick", _onClock);
    App.on("layer-visibility", _onLayerVisibility);
    App.on("layers", () => _renderGroups());
    App.on("objects", () => _renderGroups());
    App.on("basemap", () => _renderGroups());
    App.on("layer-order", _applyLayerOrder);
    App.on("clock-tick", _onCustomClock);
    App.on("styles", () => { _applyAllStyles(); _buildPanel(); });
    MapView.setHoverSampler(_hoverSample);
    _buildPanel();
  }

  function _reset() {
    meta = null; mesh = null; srcInfo = null;
    frameCache.clear(); frameRange.clear(); inflight.clear();
    dataRanges.clear(); pointLayers.clear(); pointData.clear();
    outState.clear();
    for (const l of [...Layers.byGroup("output")]) {
      FieldLayer.remove(`field-${l.id}`);
      Layers.unregister(l.id);
    }
    for (const fid of customIds) FieldLayer.remove(fid);
    customIds.clear();
    Playbar.removeSource("output");
    Playbar.setIndexTimes(null);
    _syncLegend();
  }

  function _enter() {
    if (!meta) _loadOutput();
    _buildPanel();
  }

  function _setLoading(id, busy) {
    if (busy) loadingLayers.add(id); else loadingLayers.delete(id);
    App.emit("layer-loading", { id, busy });
    _renderGroups();
  }

  /* ================= layer groups (collapsible sections) ================= */

  function _collapsedKey(group) { return `viewer-${group}`; }

  /* Display order of the groups: user-draggable, persisted. The saved
   * order is deduplicated on read: state written by older versions could
   * contain repeated group names, which would render every section (and
   * its cards) many times over. */
  function _orderedGroups() {
    const keys = GROUPS.map(([g]) => g);
    const saved = [...new Set(App.state.ui.viewerGroupOrder || [])];
    const order = [...saved.filter((g) => keys.includes(g)),
      ...keys.filter((g) => !saved.includes(g))];
    return order.map((g) => GROUPS.find(([key]) => key === g));
  }

  /* Eye button in a section header: show/hide the whole group. The
   * objects and transect groups both live in the Objects store (split by
   * kind), so their eyes toggle store objects, not registry layers. */
  function _groupEye(group, layers) {
    const storeObjects = () => (group === "transect"
      ? App.state.objects.filter((o) => o.kind === "transect")
      : App.state.objects.filter((o) => o.kind !== "transect"));
    let anyVisible;
    if (group === "objects" || group === "transect") {
      anyVisible = storeObjects().some((o) => o.visible);
    } else {
      anyVisible = layers.some((l) => l.visible);
    }
    const eye = U.el("span", {
      class: `eye group-eye ${anyVisible ? "" : "off"}`,
      title: anyVisible ? "Hide the whole group" : "Show the whole group",
    }, "👁");
    eye.addEventListener("click", (ev) => {
      ev.stopPropagation();
      if (group === "objects" || group === "transect") {
        for (const obj of storeObjects()) {
          Objects.update(obj.id, { visible: !anyVisible });
        }
        return;
      }
      for (const layer of layers) {
        // never turn a stale layer on
        if (!anyVisible && layer.stale) continue;
        if (layer.visible === !anyVisible) continue;
        layer.visible = !anyVisible;
        App.emit("layer-visibility", layer);
      }
      _renderGroups();
    });
    return eye;
  }

  // re-entrancy guard: rendering triggers registry syncs whose events can
  // ask for another render; run that ONE follow-up after this pass instead
  // of recursing (unbounded recursion here duplicated every section and
  // froze the UI)
  let _renderingGroups = false;
  let _renderQueued = false;

  function _renderGroups() {
    if (_renderingGroups) { _renderQueued = true; return; }
    _renderingGroups = true;
    try {
      _renderGroupsNow();
    } finally {
      _renderingGroups = false;
      if (_renderQueued) {
        _renderQueued = false;
        setTimeout(_renderGroups, 0);
      }
    }
  }

  function _renderGroupsNow() {
    const box = els.groups;
    if (!box) return;
    U.clear(box);
    _syncCustomRegistry();

    for (const [group, title] of _orderedGroups()) {
      let contentBuilder = null;
      let count = 0;
      let groupLayers = [];
      if (group === "objects") {
        const polys = App.state.objects.filter((o) => o.kind !== "transect");
        if (!polys.length) continue;
        count = polys.length;
        contentBuilder = (body) => _renderObjectCards(body);
      } else if (group === "transect") {
        // always shown (even empty) so the draw button stays reachable
        count = Objects.byKind("transect").length;
        contentBuilder = (body) => _renderTransectCards(body);
      } else if (group === "background") {
        if (CRS.isLocal()) continue;
        contentBuilder = (body) => _renderBackgroundCards(body);
      } else if (group === "custom") {
        // always shown (even empty) so the "add" control stays reachable
        count = (App.state.ui.customLayers || []).length;
        contentBuilder = (body) => _renderCustomCards(body);
      } else if (group === "output") {
        // always shown (even before the first output step exists) so
        // the refresh button stays reachable while a run is writing
        groupLayers = Layers.byGroup(group);
        count = groupLayers.length;
        contentBuilder = (body) => {
          if (groupLayers.length) _renderLayerCards(body, group, groupLayers);
          _outputFooter(body);
        };
      } else {
        groupLayers = Layers.byGroup(group);
        if (!groupLayers.length) continue;
        count = groupLayers.length;
        contentBuilder = (body) => _renderLayerCards(body, group, groupLayers);
      }
      const collapsed = Boolean((App.state.ui.collapsed || {})[_collapsedKey(group)]);
      const section = U.section(title, { collapsed, count: count || null });
      section.wrap.dataset.group = group;
      section.head.addEventListener("click", () => {
        App.state.ui.collapsed = App.state.ui.collapsed || {};
        App.state.ui.collapsed[_collapsedKey(group)] =
          section.wrap.classList.contains("collapsed");
        App.touchUi();
      });
      // grip to drag whole groups + group show/hide. Only the grip is
      // draggable (exactly like the layer cards), so the drag feel — grab
      // point, ghost, insertion line — is identical to reordering cells.
      section.head.prepend(U.el("span", {
        class: "drag-grip", draggable: "true",
        title: "Drag to reorder groups (top = drawn on top)",
      }, "⠿"));
      if (group !== "background" && group !== "custom") {
        const countEl = section.head.querySelector(".count");
        section.head.insertBefore(_groupEye(group, groupLayers), countEl);
      }
      contentBuilder(section.body);
      box.append(section.wrap);
    }
    if (!box.children.length) {
      box.append(U.el("div", { class: "muted" }, "No layers yet."));
    }
    _wireGroupDrag(box);
  }

  /* drag & drop of whole sections to reorder groups (wired once per
   * container element - _renderGroups re-runs often). Drag starts from a
   * section header only, so nested layer-card drags don't hijack it. The
   * new order is read straight from the DOM sequence, which makes every
   * position reachable (incl. dropping a group above the first one). */
  function _wireGroupDrag(box) {
    if (box._groupDragWired) return;
    box._groupDragWired = true;
    U.wireSortable(box, (from, to) => {
      const order = [...box.querySelectorAll(".section[data-group]")]
        .map((s) => s.dataset.group);
      const [moved] = order.splice(from, 1);
      order.splice(to, 0, moved);
      App.state.ui.viewerGroupOrder = [...new Set(order)];
      App.touchUi();
      _renderGroups();
      _applyLayerOrder();
    }, { itemSel: ".section[data-group]", handleSel: ".section > header .drag-grip" });
  }

  function _renderLayerCards(body, group, layers) {
    const list = U.el("div", { class: "obj-list" });
    layers.forEach((layer, idx) => {
      list.append(_layerCard(layer, idx, layers.length));
    });
    _wireCardDrag(list, (from, to) => {
      const ids = Layers.byGroup(group).map((l) => l.id);
      const [moved] = ids.splice(from, 1);
      ids.splice(to, 0, moved);
      Layers.reorderGroup(group, ids);
    });
    body.append(list);
  }

  function _layerCard(layer, idx, count) {
    const loading = loadingLayers.has(layer.id);
    let eye;
    if (loading) {
      eye = U.el("span", { class: "spin", title: "Loading…" });
    } else {
      // stale layers cannot be turned ON (a visible one can be hidden)
      const stale = Boolean(layer.stale) && !layer.visible;
      eye = U.el("span", {
        class: `eye ${layer.visible ? "" : "off"} ${stale ? "disabled" : ""}`,
        title: stale
          ? "This file was made for an older grid — re-interpolate it in the Domain tab"
          : "Show/hide",
      }, "👁");
      if (!stale) {
        eye.addEventListener("click", () => {
          layer.visible = !layer.visible;
          App.emit("layer-visibility", layer);
          _renderGroups();
        });
      }
    }

    const card = U.el("div", {
      class: "obj-card", dataset: { idx },
    },
      count > 1 ? U.el("span", { class: "drag-grip", draggable: "true", title: "Drag to reorder (top = drawn on top)" }, "⠿") : null,
      eye,
      U.el("span", { class: "lp-name", title: layer.title }, layer.title),
      layer.subtitle ? U.el("span", { class: "lp-mini" }, layer.subtitle) : null);

    if (layer.stale) {
      card.append(U.el("span", {
        class: "warn-icon",
        title: "Grid changed after this file was interpolated — re-interpolate in the Domain tab",
      }, "⚠"));
    }

    // extra output dimensions (e.g. sediment fractions): cycle button
    // per dim, like the species selector on stacked vegetation grids
    if (_isOut(layer.id) && meta) {
      const varName = _outVar(layer.id);
      const info = meta.variables.find((v) => v.name === varName);
      for (const [di, dim] of ((info && info.extra_dims) || []).entries()) {
        const st = _outStateOf(varName);
        const cur = st.dims[di] || 0;
        const dimBtn = U.el("button", {
          class: "mini-btn", style: "width:auto;padding:0 6px;font-size:11px",
          title: `${dim.name} ${cur + 1} of ${dim.size} — click for next`,
        }, `${dim.name} ${cur + 1}/${dim.size}`);
        dimBtn.addEventListener("click", async (ev) => {
          ev.stopPropagation();
          st.dims[di] = ((st.dims[di] || 0) + 1) % dim.size;
          st.bracket = null;
          if (layer.visible !== false) {
            await _updateVarFrames(varName, App.state.clock.t).catch(() => {});
          }
          _renderGroups();
        });
        card.append(dimBtn);
      }
    }

    // species selector for stacked vegetation grids (hveg/Nt)
    if (layer.species > 1) {
      const k = layer.speciesIdx || 0;
      const spBtn = U.el("button", {
        class: "mini-btn", style: "width:auto;padding:0 6px;font-size:11px",
        title: `Showing species ${k + 1} of ${layer.species} — click for next`,
      }, `${k + 1}/${layer.species}`);
      spBtn.addEventListener("click", async () => {
        layer.speciesIdx = (k + 1) % layer.species;
        if (layer.visible) await _toggleDomainLayer(layer);
        _renderGroups();
      });
      card.append(spBtn);
    }

    // styling: every colormapped layer gets a style dropdown (saved
    // presets + custom) and a gear opening the full editor. Both work
    // whether or not the layer is currently shown — the assignment
    // persists and applies when the layer is (re)shown.
    const fieldId = _fieldIdFor(layer);
    const isColormapped = Boolean(fieldId) || pointLayers.has(layer.id);
    if (isColormapped) {
      card.append(_styleSelect(layer));
      card.append(U.miniBtn("gear", "Layer style (colormap, limits, opacity…)",
        () => _styleEditor(layer)));
    }
    return card;
  }

  /* Footer of the Model-output group: probe + zoom + refresh + step
   * count. Also rendered when no output is loaded yet, so the refresh
   * button is reachable while a run has not written its first step. */
  function _outputFooter(body) {
    const row = U.el("div", { class: "btn-row", style: "margin-top:6px" });
    if (meta) {
      const probeBtn = U.el("button", { class: "ghost", style: "font-size:12px" },
        "Probe cell (click map)");
      probeBtn.addEventListener("click", () => _armProbe(probeBtn));
      const zoomBtn = U.el("button", { class: "ghost", style: "font-size:12px" },
        "Zoom to output");
      zoomBtn.addEventListener("click", _zoomToOutput);
      row.append(probeBtn, zoomBtn);
    }
    const refreshBtn = U.el("button", {
      class: "ghost", style: "font-size:12px",
      title: "Re-read the output file — picks up steps written since it "
        + "was loaded (e.g. by a simulation that is still running)",
    }, "⟳ Refresh");
    refreshBtn.addEventListener("click", () => _refreshOutput(refreshBtn));
    row.append(refreshBtn);
    body.append(row,
      U.el("div", { class: "muted", style: "font-size:11.5px" },
        meta
          ? `${meta.times.length} output steps — scrub or play with the time bar below.`
          : "No output loaded — Refresh checks for a (new) output file."));
    // Output-source switch: shown as soon as this project has an HPC
    // run folder to read from (or an override is already active), so
    // switching between the local file and the P-drive run is one
    // click either way. The active side is highlighted and disabled.
    if (srcInfo && (srcInfo.override || srcInfo.known)) {
      const remoteDir = srcInfo.override
        ? srcInfo.dir
        : (srcInfo.known && srcInfo.known.dir);
      const remoteCfg = srcInfo.known ? srcInfo.known.config : null;
      const switchTo = async (dir, config) => {
        try {
          const res = await Api.post("/api/output/source",
            dir ? { dir, config: config || null } : { dir: null });
          U.toast(dir
            ? (res.exists
              ? "Viewer output: the HPC run folder on P:"
              : "Following the run folder (no output file written yet)")
            : "Viewer output: the project's own file", "ok");
          // the output-source listeners reload the layers AND the
          // data-availability strip in the graphs panel
          App.emit("output-source");
        } catch (err) { U.toast(err.message, "error"); }
      };
      const segBtn = (label, active, title, onClick) => {
        const b = U.el("button", {
          class: active ? "primary" : "ghost",
          style: "font-size:11px", title,
          disabled: active ? "" : null,
        }, label);
        if (!active) b.addEventListener("click", onClick);
        return b;
      };
      const jobBit = srcInfo.known && srcInfo.known.job_id
        ? ` (job ${srcInfo.known.job_id})` : "";
      body.append(
        U.el("div", {
          class: "muted",
          style: "font-size:11.5px;display:flex;gap:6px;align-items:center;flex-wrap:wrap;margin-top:4px",
        },
          U.el("span", {}, "Output source:"),
          segBtn("Project file", !srcInfo.override,
            "Read the output file in the project folder",
            () => switchTo(null)),
          remoteDir ? segBtn(`HPC run${jobBit}`, srcInfo.override,
            `Read the output written by the HPC run in ${remoteDir}`,
            () => switchTo(remoteDir, remoteCfg)) : null),
        U.el("div", {
          class: "muted",
          style: "font-size:11px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap",
          title: srcInfo.path || "",
        }, `reading: ${srcInfo.path || "?"}`));
    }
  }

  function _zoomToOutput() {
    if (!mesh) return;
    let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
    for (let i = 0; i < mesh.x.length; i += 1) {
      if (mesh.x[i] < minX) minX = mesh.x[i];
      if (mesh.x[i] > maxX) maxX = mesh.x[i];
      if (mesh.y[i] < minY) minY = mesh.y[i];
      if (mesh.y[i] > maxY) maxY = mesh.y[i];
    }
    MapView.fitModelBounds(minX, minY, maxX, maxY);
  }

  /* ================= per-layer styles (presets or custom) =================
   * Every colormapped layer resolves its style from ui.layerStyles:
   *   {preset: "<id>"} -> a saved Styles preset (shared, editable)
   *   {custom: {...}}  -> a one-off style for this layer
   *   absent           -> a sensible default inferred from the name
   * A style = {cmap, mode, dotSize, opacity, min, max}; min/max null
   * means "auto" (the layer's own data range). */

  function _styleKeyOf(layer) {
    return _isOut(layer.id) ? `output:${_outVar(layer.id)}` : layer.id;
  }

  function _assignments() {
    return App.state.ui.layerStyles || (App.state.ui.layerStyles = {});
  }

  const DEFAULT_CMAPS = {
    elevation: "topo_dutch", bed_change: "RdBu", veg_density: "Greens",
    veg_height: "Greens", mask: "gray", wind_speed: "viridis",
    shear_stress: "plasma", shear_velocity: "plasma",
    sed_conc: "sand", sed_transport: "sand",
  };

  function _defaultStyleFor(layer) {
    const name = _isOut(layer.id) ? _outVar(layer.id) : _quantityOf(layer);
    let cmap = DEFAULT_CMAPS[_categoryOfName(name)] || null;
    if (!cmap && layer.id.startsWith("custom-")) {
      const c = _customList().find((x) => `custom-${x.id}` === layer.id);
      cmap = (c && (c.op === "sub" || c.op === "offset")) ? "RdBu" : "viridis";
    }
    if (!cmap && layer.entry && (layer.entry.bands || 1) > 1) cmap = "gray";
    return { cmap: cmap || "viridis", min: null, max: null,
      opacity: _isOut(layer.id) ? 1 : 0.9, mode: "cells", dotSize: 6 };
  }

  /* The style DEFINITION for a layer (min/max may be null = auto). */
  function _styleDefOf(layer) {
    const a = _assignments()[_styleKeyOf(layer)];
    if (a && a.preset) {
      const p = Styles.get(a.preset);
      if (p) return { ...p, _preset: p.id, _name: p.name };
    }
    if (a && a.custom) return { ...a.custom };
    return _defaultStyleFor(layer);
  }

  /* Concrete draw style: auto limits filled from the layer's data range. */
  function _styleOf(layer) {
    const def = _styleDefOf(layer);
    const r = dataRanges.get(layer.id) || [0, 1];
    return {
      cmap: def.cmap || "viridis",
      mode: def.mode === "dots" ? "dots" : "cells",
      dotSize: Number(def.dotSize) || 6,
      opacity: def.opacity != null ? Number(def.opacity) : 0.9,
      min: def.min != null ? def.min : r[0],
      max: def.max != null ? def.max : r[1],
      _preset: def._preset || null, _name: def._name || null,
    };
  }

  /* Push a layer's resolved style to its live map object. */
  function _applyLayerStyle(layer) {
    const st = _styleOf(layer);
    const map = MapView.instance();
    if (_isOut(layer.id)) {
      const varName = _outVar(layer.id);
      const def = _styleDefOf(layer);
      const state = _outStateOf(varName);
      state.auto = def.min == null && def.max == null;
      const l = FieldLayer.get(`field-${layer.id}`);
      if (l) {
        const patch = { cmap: st.cmap, opacity: st.opacity,
          mode: st.mode, dotSize: st.dotSize };
        if (!state.auto) { patch.min = st.min; patch.max = st.max; }
        l.setStyle(patch);
        if (state.auto && state.bracket) _applyAutoRange(varName, state.bracket.k);
      }
    } else if (pointLayers.has(layer.id)) {
      if (map.getLayer(`pts-${layer.id}`)) {
        map.setPaintProperty(`pts-${layer.id}`, "circle-color",
          _pointColorExpr(st.cmap, st.min, st.max));
        map.setPaintProperty(`pts-${layer.id}`, "circle-radius", Math.max(1, st.dotSize / 2.5));
        map.setPaintProperty(`pts-${layer.id}`, "circle-opacity", st.opacity);
      }
    } else {
      const fl = FieldLayer.get(_fieldIdFor(layer));
      if (fl) {
        fl.setStyle({ cmap: st.cmap, min: st.min, max: st.max,
          opacity: st.opacity, mode: st.mode, dotSize: st.dotSize });
      }
    }
    _syncLegend();
  }

  function _assignPreset(layer, presetId) {
    _assignments()[_styleKeyOf(layer)] = { preset: presetId };
    App.touchUi();
    _applyLayerStyle(layer);
  }

  /* Patch a layer's custom style (detaches it from any preset). */
  function _updateCustom(layer, patch) {
    const cur = _styleDefOf(layer);
    const def = { cmap: cur.cmap, mode: cur.mode || "cells",
      dotSize: cur.dotSize != null ? cur.dotSize : 6,
      opacity: cur.opacity != null ? cur.opacity : 0.9,
      min: cur.min != null ? cur.min : null,
      max: cur.max != null ? cur.max : null, ...patch };
    _assignments()[_styleKeyOf(layer)] = { custom: def };
    App.touchUi();
    _applyLayerStyle(layer);
  }

  /* Re-apply every visible colormapped layer's style (after a preset
   * changed or was deleted). */
  function _applyAllStyles() {
    for (const layer of App.state.layers) {
      if (layer.visible === false) continue;
      if (_fieldIdFor(layer) || pointLayers.has(layer.id)) {
        _applyLayerStyle(layer);
      }
    }
    _syncLegend();
  }

  /* Compact style dropdown on a layer card: saved presets + custom. */
  function _styleSelect(layer) {
    const sel = U.el("select", { class: "cmap-mini", title: "Colormap style for this layer" });
    const a = _assignments()[_styleKeyOf(layer)];
    const cur = a && a.preset && Styles.get(a.preset) ? a.preset : "__custom";
    for (const p of Styles.all()) {
      sel.append(U.el("option", { value: p.id, selected: p.id === cur ? "" : null }, p.name));
    }
    sel.append(U.el("option", { value: "__custom", selected: cur === "__custom" ? "" : null }, "custom…"));
    sel.addEventListener("click", (ev) => ev.stopPropagation());
    sel.addEventListener("change", () => {
      if (sel.value === "__custom") _styleEditor(layer);
      else _assignPreset(layer, sel.value);
    });
    return sel;
  }

  /* drag & drop within one list (shared, de-lagged, insertion-line UI) */
  function _wireCardDrag(list, onDrop) {
    U.wireSortable(list, onDrop);
  }

  function _fieldIdFor(layer) {
    if (_isOut(layer.id) || layer.id.startsWith("domain-")
        || layer.id.startsWith("raw-") || layer.id.startsWith("custom-")) {
      return `field-${layer.id}`;
    }
    return null;
  }

  function _renderTransectCards(body) {
    const list = U.el("div", { class: "obj-list" });
    const transects = Objects.byKind("transect");
    for (const obj of transects) {
      const eye = U.el("span", { class: `eye ${obj.visible ? "" : "off"}`, title: "Show/hide" }, "👁");
      eye.addEventListener("click", () => Objects.update(obj.id, { visible: !obj.visible }));
      const name = U.el("span", { class: "lp-name", title: "Double-click to rename" }, obj.name);
      name.addEventListener("dblclick", () => {
        const input = U.el("input", { type: "text", value: obj.name, style: "flex:1;font-size:12px" });
        name.replaceWith(input); input.focus(); input.select();
        const commit = () => Objects.update(obj.id, { name: input.value.trim() || obj.name });
        input.addEventListener("blur", commit);
        input.addEventListener("keydown", (ev) => {
          if (ev.key === "Enter") input.blur();
          if (ev.key === "Escape") { input.value = obj.name; input.blur(); }
        });
      });
      const zoom = U.miniBtn("target", "Zoom to", () => {
        const xs = obj.coords.map((c) => c[0]), ys = obj.coords.map((c) => c[1]);
        MapView.fitModelBounds(Math.min(...xs), Math.min(...ys), Math.max(...xs), Math.max(...ys));
      });
      const del = U.miniBtn("trash", "Delete", () => {
        if (window.confirm(`Delete ${obj.name}?`)) Objects.remove(obj.id);
      });
      del.classList.add("danger-hover");
      list.append(U.el("div", { class: "obj-card" }, eye, name,
        U.el("span", { class: "lp-mini" }, "transect"), zoom, del));
    }
    if (!transects.length) {
      list.append(U.el("div", { class: "muted" }, "No transects — draw one, then plot it in the Transect tab."));
    }
    body.append(list);
    const drawBtn = U.el("button", { class: "add-optional" }, "+ Draw transect");
    drawBtn.addEventListener("click", _drawTransect);
    body.append(drawBtn);
  }

  async function _drawTransect() {
    try {
      const n = Objects.byKind("transect").length + 1;
      await Draw.transect({ name: `Transect ${n}`, color: "#111111" });
      _renderGroups();
    } catch (e) { /* draw cancelled */ }
  }

  function _renderObjectCards(body) {
    const list = U.el("div", { class: "obj-list" });
    for (const obj of App.state.objects.filter((o) => o.kind !== "transect")) {
      const eye = U.el("span", { class: `eye ${obj.visible ? "" : "off"}` }, "👁");
      eye.addEventListener("click", () => Objects.update(obj.id, { visible: !obj.visible }));

      const swatch = U.el("button", {
        class: "mini-btn", style: `color:${obj.color}`, title: "Change colour",
      }, "■");
      swatch.addEventListener("click", () => {
        const palette = ["#e6552f", "#2f7fe6", "#27a355", "#a034c6", "#e0a020", "#12a5b5", "#d1387f"];
        const next = palette[(palette.indexOf(obj.color) + 1) % palette.length];
        Objects.update(obj.id, { color: next });
      });

      const name = U.el("span", { class: "lp-name", title: "Double-click to rename" }, obj.name);
      name.addEventListener("dblclick", () => {
        const input = U.el("input", { type: "text", value: obj.name, style: "flex:1;font-size:12px" });
        name.replaceWith(input);
        input.focus(); input.select();
        const commit = () => Objects.update(obj.id, { name: input.value.trim() || obj.name });
        input.addEventListener("blur", commit);
        input.addEventListener("keydown", (ev) => {
          if (ev.key === "Enter") input.blur();
          if (ev.key === "Escape") { input.value = obj.name; input.blur(); }
        });
      });

      const zoom = U.miniBtn("target", "Zoom to", () => {
        const xs = obj.coords.map((c) => c[0]), ys = obj.coords.map((c) => c[1]);
        MapView.fitModelBounds(Math.min(...xs), Math.min(...ys), Math.max(...xs), Math.max(...ys));
      });

      const del = U.miniBtn("trash", "Delete", () => {
        if (window.confirm(`Delete ${obj.name}?`)) Objects.remove(obj.id);
      });
      del.classList.add("danger-hover");

      list.append(U.el("div", { class: "obj-card" }, eye, swatch, name,
        U.el("span", { class: "lp-mini" }, obj.kind), zoom, del));
    }
    body.append(list);
  }

  function _renderBackgroundCards(body) {
    const list = U.el("div", { class: "obj-list" });
    for (const [key, label] of [["gray", "Grey map"], ["sat", "Satellite"], ["none", "None"]]) {
      const active = App.state.ui.basemap === key;
      const card = U.el("div", { class: `obj-card ${active ? "selected" : ""}`, style: "cursor:pointer" },
        U.el("span", { class: "eye" }, active ? "●" : "○"),
        U.el("span", { class: "lp-name" }, label));
      card.addEventListener("click", () => MapView.setBasemap(key));
      list.append(card);
    }
    body.append(list);
  }

  /* ================= custom computed layers =================
   * A custom layer combines OUTPUT variable snapshots on the shared
   * output mesh (e.g. bed-level change zb@t2 − zb@t1, or a scalar
   * multiple). Definitions persist per project in ui.customLayers;
   * arithmetic runs client-side on the cached frame arrays. */

  const CUSTOM_OPS = [
    ["sub", "A − B (difference)"],
    ["add", "A + B"],
    ["mul", "A × B"],
    ["div", "A ÷ B"],
    ["scale", "k × A (scalar)"],
    ["offset", "A + k (scalar)"],
  ];
  const _binaryOp = (op) => op !== "scale" && op !== "offset";
  const SENTINEL = -1e30;
  const _bad = (v) => !(v > -1e29);   // NaN or nodata sentinel

  function _customList() { return App.state.ui.customLayers || (App.state.ui.customLayers = []); }

  const customBusy = new Set();   // custom ids currently (re)computing

  function _syncCustomRegistry() {
    const cl = _customList();
    const wanted = new Set(cl.map((c) => `custom-${c.id}`));
    for (const c of cl) {
      Layers.register({ id: `custom-${c.id}`, group: "custom", title: c.name, visible: !!c.visible });
    }
    for (const l of Layers.byGroup("custom")) {
      if (!wanted.has(l.id)) Layers.unregister(l.id);
    }
  }

  // resolve one operand: "current" interpolates the two bracketing frames
  async function _operandArray(varName, timeSel) {
    if (!meta || !varName) return null;
    if (timeSel === "current" || timeSel == null) {
      const { k, frac } = _bracket(App.state.clock.t);
      const k2 = Math.min(k + 1, meta.times.length - 1);
      const [A, B] = await Promise.all([
        _fetchFrame(varName, "0", k), _fetchFrame(varName, "0", k2)]);
      if (!frac || A === B) return A;
      const out = new Float32Array(A.length);
      for (let i = 0; i < A.length; i += 1) {
        out[i] = (_bad(A[i]) || _bad(B[i])) ? SENTINEL : A[i] + frac * (B[i] - A[i]);
      }
      return out;
    }
    const t = Math.max(0, Math.min(Number(timeSel) || 0, meta.times.length - 1));
    return _fetchFrame(varName, "0", t);
  }

  // ---- operands: any grid of consistent shape (output var@time OR a loaded grid) ----
  // grid (raster) sources for custom-layer operands (async: uses the catalog)
  async function _fieldSources() {
    return (await datasetCatalog()).filter((s) => s.kind === "raster");
  }

  function _srcLabel(key, time) {
    if (!key) return "?";
    const base = _srcLabels[key] || key.replace(/^(out|tgt|ent):/, "");
    if (key.startsWith("out:")) {
      return `${key.slice(4)}@${(time === "current" || time == null) ? "now" : "#" + time}`;
    }
    return base;
  }

  async function _resolveOperand(key, time) {
    if (key && key.startsWith("out:")) {
      if (!mesh) return null;
      const arr = await _operandArray(key.slice(4), time || "current");
      return arr ? { arr, mesh } : null;
    }
    const src = await _ensureSource(key);   // tgt: / ent:
    if (src && src.mesh && src.data) return { arr: src.data, mesh: src.mesh };
    return null;   // point datasets can't be used for grid arithmetic
  }

  function _sameMesh(m1, m2) {
    if (!m1 || !m2 || m1.x.length !== m2.x.length) return false;
    const n = m1.x.length;
    const mid = n >> 1;
    const eq = (a, b) => Math.abs(a - b) < 1e-6 * (1 + Math.abs(a));
    return eq(m1.x[0], m2.x[0]) && eq(m1.y[0], m2.y[0])
      && eq(m1.x[mid], m2.x[mid]) && eq(m1.y[mid], m2.y[mid])
      && eq(m1.x[n - 1], m2.x[n - 1]) && eq(m1.y[n - 1], m2.y[n - 1]);
  }

  // nearest-neighbour resample of srcData (on srcMesh) onto dstMesh, with a
  // cached index map (so per-clock recomputes stay cheap)
  const _resampleMaps = new WeakMap();
  function _nearestMap(dstMesh, srcMesh) {
    let byDst = _resampleMaps.get(dstMesh);
    if (!byDst) { byDst = new WeakMap(); _resampleMaps.set(dstMesh, byDst); }
    let idx = byDst.get(srcMesh);
    if (idx) return idx;
    idx = new Int32Array(dstMesh.x.length);
    for (let d = 0; d < dstMesh.x.length; d += 1) {
      idx[d] = _nearestIdx(srcMesh, dstMesh.x[d], dstMesh.y[d]).best;
    }
    byDst.set(srcMesh, idx);
    return idx;
  }
  function _resampleTo(dstMesh, src) {
    const idx = _nearestMap(dstMesh, src.mesh);
    const out = new Float32Array(idx.length);
    for (let d = 0; d < idx.length; d += 1) out[d] = src.arr[idx[d]];
    return { arr: out, mesh: dstMesh };
  }

  async function _computeCustom(c) {
    const A = await _resolveOperand(c.aKey, c.aTime);
    if (!A) return null;
    let B = null;
    if (_binaryOp(c.op)) {
      B = await _resolveOperand(c.bKey, c.bTime);
      if (!B) return null;
      // operands must share the grid; optionally resample B onto A's grid
      if (A.mesh !== B.mesh && !_sameMesh(A.mesh, B.mesh)) {
        if (c.resample) B = _resampleTo(A.mesh, B);
        else return { error: "shape-mismatch" };
      }
    }
    const S = Number(c.scalar) || 0;
    const a = A.arr, bArr = B && B.arr;
    const out = new Float32Array(a.length);
    for (let i = 0; i < a.length; i += 1) {
      const av = a[i];
      if (_bad(av)) { out[i] = SENTINEL; continue; }
      let v;
      if (c.op === "scale") v = S * av;
      else if (c.op === "offset") v = av + S;
      else {
        const bv = bArr[i];
        if (_bad(bv)) { out[i] = SENTINEL; continue; }
        v = c.op === "sub" ? av - bv : c.op === "add" ? av + bv
          : c.op === "mul" ? av * bv : (bv !== 0 ? av / bv : NaN);
      }
      out[i] = Number.isFinite(v) ? v : SENTINEL;
    }
    return { out, mesh: A.mesh };
  }

  async function _recomputeCustom(c) {
    const res = await _computeCustom(c);
    if (!res) return;
    const fid = `field-custom-${c.id}`;
    if (res.error) {                       // operands not ready or mismatched shapes
      if (res.error === "shape-mismatch" && c._warned !== true) {
        c._warned = true;
        U.toast(`"${c.name}": datasets have different grids — tick "resample to match" in the layer`, "error");
      }
      return;
    }
    c._warned = false;
    let layer = FieldLayer.get(fid);
    if (!layer || layer.mesh !== res.mesh) {
      FieldLayer.remove(fid);
      layer = FieldLayer.create(fid, res.mesh, { cmap: "RdBu", min: -1, max: 1, opacity: 1 });
      customIds.add(fid);
    }
    // data range for auto limits — symmetric around 0 for differences so
    // the diverging default colormap centres correctly
    let lo = Infinity, hi = -Infinity;
    for (let i = 0; i < res.out.length; i += 1) {
      const v = res.out[i];
      if (!_bad(v) && Number.isFinite(v)) { if (v < lo) lo = v; if (v > hi) hi = v; }
    }
    if (lo > hi) { lo = 0; hi = 1; }
    if (c.op === "sub" || c.op === "offset") {
      const m = Math.max(Math.abs(lo), Math.abs(hi)) || 1;
      lo = -m; hi = m;
    }
    dataRanges.set(`custom-${c.id}`, [lo, hi]);
    layer.setFrames(res.out);
    _applyLayerStyle({ id: `custom-${c.id}` });
  }

  async function _restoreCustomLayers() {
    for (const c of _customList()) if (c.visible) await _recomputeCustom(c);
    _syncCustomRegistry();
  }

  // recompute custom layers that depend on the current time as the clock moves
  function _onCustomClock() {
    for (const c of _customList()) {
      if (!c.visible) continue;
      const usesNow = (c.aKey && c.aKey.startsWith("out:") && (c.aTime === "current" || c.aTime == null)) ||
                      (_binaryOp(c.op) && c.bKey && c.bKey.startsWith("out:") && (c.bTime === "current" || c.bTime == null));
      if (usesNow) _recomputeCustom(c).catch(() => {});
    }
  }

  async function _toggleCustom(c) {
    c.visible = !c.visible;
    App.touchUi();
    if (c.visible) {
      customBusy.add(c.id);
      _renderGroups();
      try { await _recomputeCustom(c); }
      finally { customBusy.delete(c.id); }
    } else { FieldLayer.remove(`field-custom-${c.id}`); customIds.delete(`field-custom-${c.id}`); }
    _syncCustomRegistry();
    _applyLayerOrder();
    _renderGroups();
  }

  function _renderCustomCards(body) {
    const list = U.el("div", { class: "obj-list" });
    _customList().forEach((c) => list.append(_customCard(c)));
    body.append(list);
    const addBtn = U.el("button", { class: "add-optional" }, "+ Add custom layer");
    addBtn.addEventListener("click", () => _customDialog(null));
    addBtn.title = "Combine output variables / timesteps or any model-file & raster grids";
    body.append(addBtn);
  }

  function _customCard(c) {
    let eye;
    if (customBusy.has(c.id)) {
      eye = U.el("span", { class: "spin", title: "Computing…" });
    } else {
      eye = U.el("span", { class: `eye ${c.visible ? "" : "off"}`, title: "Show/hide" }, "\u{1F441}");
      eye.addEventListener("click", () => _toggleCustom(c));
    }
    const cmapSel = _styleSelect({ id: `custom-${c.id}` });
    const styleBtn = U.miniBtn("gear", "Layer style (colormap, limits, opacity…)",
      () => _styleEditor({ id: `custom-${c.id}`, title: c.name }));
    const edit = U.miniBtn("modify", "Edit…", () => _customDialog(c));
    const del = U.miniBtn("trash", "Delete", () => {
      FieldLayer.remove(`field-custom-${c.id}`);
      customIds.delete(`field-custom-${c.id}`);
      App.state.ui.customLayers = _customList().filter((x) => x.id !== c.id);
      App.touchUi();
      _syncCustomRegistry();
      _renderGroups();
    });
    del.classList.add("danger-hover");
    return U.el("div", { class: "obj-card" },
      eye,
      U.el("span", { class: "lp-name", title: _customFormula(c) }, c.name),
      U.el("span", { class: "lp-mini" }, _customFormula(c)),
      cmapSel, styleBtn, edit, del);
  }

  function _customFormula(c) {
    const a = _srcLabel(c.aKey, c.aTime);
    if (c.op === "scale") return `${c.scalar} × ${a}`;
    if (c.op === "offset") return `${a} + ${c.scalar}`;
    const b = _srcLabel(c.bKey, c.bTime);
    const sym = { sub: "−", add: "+", mul: "×", div: "÷" }[c.op] || "?";
    return `${a} ${sym} ${b}`;
  }

  async function _customDialog(existing) {
    const sources = await _fieldSources();
    if (!sources.length) { U.toast("Load a layer or run a model output first", "error"); return; }
    const keys = sources.map((s) => s.key);
    const c = existing || {
      id: `${Date.now()}`, name: "", op: "sub",
      aKey: keys[0], aTime: "current", bKey: keys[1] || keys[0], bTime: "current",
      scalar: 1,
    };
    const popup = Popup.open({ title: existing ? "Edit custom layer" : "New custom layer", width: 480 });

    const name = U.el("input", { type: "text", value: c.name, placeholder: "e.g. bed level change" });
    const opSel = U.el("select", {});
    for (const [k, lbl] of CUSTOM_OPS) opSel.append(U.el("option", { value: k, selected: k === c.op ? "" : null }, lbl));

    const srcSel = (val) => {
      const s = U.el("select", {});
      for (const src of sources) s.append(U.el("option", { value: src.key, selected: src.key === val ? "" : null }, src.label));
      return s;
    };
    const timeSel = (val) => {
      const s = U.el("select", { style: "max-width:150px" });
      if (meta) {
        s.append(U.el("option", { value: "current", selected: val === "current" ? "" : null }, "Current time"));
        for (let i = 0; i < meta.times.length; i += 1) {
          const lbl = meta.times_epoch ? `#${i} · ${U.fmtDate(meta.times_epoch[i])}` : `#${i}`;
          s.append(U.el("option", { value: String(i), selected: String(val) === String(i) ? "" : null }, lbl));
        }
      }
      return s;
    };
    const aSel = srcSel(c.aKey), aTime = timeSel(c.aTime);
    const bSel = srcSel(c.bKey), bTime = timeSel(c.bTime);
    const scalar = U.el("input", { type: "number", step: "any", value: c.scalar });

    // the time selector only applies when the operand is an output variable
    const syncTime = (sel, timeEl) => { timeEl.style.display = sel.value.startsWith("out:") ? "" : "none"; };
    aSel.addEventListener("change", () => syncTime(aSel, aTime));
    bSel.addEventListener("change", () => syncTime(bSel, bTime));

    const resampleCb = U.el("input", { type: "checkbox" });
    resampleCb.checked = !!c.resample;

    const rowA = U.el("div", { class: "form-row" }, U.el("label", {}, "A"), aSel, aTime);
    const rowB = U.el("div", { class: "form-row" }, U.el("label", {}, "B"), bSel, bTime);
    const rowK = U.el("div", { class: "form-row" }, U.el("label", {}, "Scalar k"), scalar);
    const rowR = U.el("div", { class: "form-row" }, U.el("label", {}, "Different grids"),
      U.el("label", { class: "choice-row", style: "font-size:12px" },
        resampleCb, U.el("span", {}, "resample B onto A's grid (nearest)")));
    const syncRows = () => {
      const bin = _binaryOp(opSel.value);
      rowB.style.display = bin ? "" : "none";
      rowR.style.display = bin ? "" : "none";
      rowK.style.display = bin ? "none" : "";
    };
    opSel.addEventListener("change", syncRows);
    syncRows(); syncTime(aSel, aTime); syncTime(bSel, bTime);

    const saveBtn = U.el("button", { class: "primary" }, "Save");
    saveBtn.addEventListener("click", async () => {
      const def = {
        id: c.id, op: opSel.value,
        aKey: aSel.value, aTime: aTime.value, bKey: bSel.value, bTime: bTime.value,
        scalar: Number(scalar.value) || 0,
        resample: resampleCb.checked,
        visible: existing ? c.visible : true,
      };
      def.name = name.value.trim() || _customFormula(def);
      const cl = _customList();
      const idx = cl.findIndex((x) => x.id === def.id);
      if (idx >= 0) cl[idx] = def; else cl.push(def);
      App.touchUi();
      FieldLayer.remove(`field-custom-${def.id}`);
      customIds.delete(`field-custom-${def.id}`);
      // computing can take a while (frame/grid fetches) — show it
      saveBtn.disabled = true;
      saveBtn.textContent = "Computing…";
      saveBtn.prepend(U.el("span", { class: "spin", style: "margin-right:6px" }));
      customBusy.add(def.id);
      _renderGroups();
      try {
        if (def.visible) await _recomputeCustom(def);
        popup.close();
      } catch (err) {
        U.toast(err.message, "error");
        saveBtn.disabled = false;
        saveBtn.textContent = "Save";
      } finally {
        customBusy.delete(def.id);
        _syncCustomRegistry();
        _applyLayerOrder();
        _renderGroups();
      }
    });
    const cancelBtn = U.el("button", { class: "ghost" }, "Cancel");
    cancelBtn.addEventListener("click", popup.close);

    popup.body.append(
      U.el("div", { class: "form-row" }, U.el("label", {}, "Name"), name),
      U.el("div", { class: "form-row" }, U.el("label", {}, "Compute"), opSel),
      rowA, rowB, rowK, rowR,
      U.el("div", { class: "muted", style: "font-size:11.5px" },
        "Combine output timesteps or any model-file / raster grids. If the two grids differ, tick "
        + "“resample”. Colours are set on the layer card (style dropdown / gear)."),
      U.el("div", { class: "btn-row", style: "justify-content:flex-end;margin-top:8px" }, cancelBtn, saveBtn),
    );
  }

  /* Apply the tree order to the map. Groups draw in their (draggable)
   * tree order - the first group/row of the tree ends up on top - so
   * e.g. the grid can be placed in front of or behind the objects. */
  const GROUP_MAP_IDS = {
    grid: ["grid-fill", "grid-lines", "grid-outline", "grid-shear", "grid-shear-inner"],
    objects: ["objects-fill", "objects-fill-outline"],
    // transect lines (objects-lines) are floated to the very top separately
  };

  function _applyLayerOrder() {
    const map = MapView.instance();
    if (!map || !map.isStyleLoaded()) return;
    const ordered = [];   // bottom -> top
    const groups = _orderedGroups().map(([g]) => g).reverse();
    for (const group of groups) {
      if (GROUP_MAP_IDS[group]) {
        ordered.push(...GROUP_MAP_IDS[group]);
        continue;
      }
      const layers = Layers.byGroup(group);
      for (let i = layers.length - 1; i >= 0; i -= 1) {
        const fieldId = _fieldIdFor(layers[i]);
        if (fieldId) ordered.push(fieldId);
        ordered.push(`pts-${layers[i].id}`);
      }
    }
    for (const id of ordered) {
      if (map.getLayer(id)) map.moveLayer(id);
    }
    // transect lines always float on top of everything
    if (map.getLayer("objects-lines")) map.moveLayer("objects-lines");
  }

  /* ================= quantities & default colormaps ================= */

  /* Which physical quantity a layer shows (bed, ne, veg, ...). Raw
   * downloads are elevation samples -> "bed"; layers duplicated from a
   * .grd keep that target's quantity (label "<target> (file) → samples"). */
  function _quantityOf(layer) {
    if (_isOut(layer.id)) return _outVar(layer.id);
    if (layer.id.startsWith("domain-")) return layer.id.slice("domain-".length);
    if (layer.id.startsWith("raw-")) {
      const entry = layer.entry || {};
      // multi-band rasters are imagery (CIR/RGB reflectance 0-255), not
      // elevation - linking them to the elevation colormap+limits would
      // clamp everything to one constant colour
      if ((entry.bands || 1) > 1) return null;
      if (entry.source === "converted") {
        const m = (entry.label || "").match(/^([A-Za-z_]\w*)\s*\(/);
        if (m) return m[1];
      }
      if (entry.source === "derived") {
        // band-math products (NDVI, masks, …): adopt the label only when it
        // clearly names a known quantity — never assume elevation
        return _categoryOfName(entry.label) ? entry.label : null;
      }
      return "bed";
    }
    return null;
  }

  /* Map a quantity/target/output-variable name to a category (used only
   * to pick a sensible DEFAULT colormap for unstyled layers). */
  function _categoryOfName(name) {
    const n = String(name || "").toLowerCase();
    if (/mask/.test(n)) return "mask";
    if (/dzb|dzbveg|sedero|erosion|deposition|bedchange|bed_change/.test(n)) return "bed_change";
    if (/hveg/.test(n)) return "veg_height";
    if (/rhoveg|^nt$|veget|^veg$/.test(n)) return "veg_density";
    if (/ustar/.test(n)) return "shear_velocity";
    if (/tau|shear/.test(n)) return "shear_stress";
    if (/^ct$|^cu$|conc/.test(n)) return "sed_conc";
    if (/^q[sxn]?$|transport|flux|pickup/.test(n)) return "sed_transport";
    if (/uw|wind|^u10|u_star|^u$/.test(n)) return "wind_speed";
    if (/zb|zne|^ne$|zsep|^zs$|bed|elev|topo/.test(n)) return "elevation";
    return null;
  }

  /* ================= colormap preset management ================= */

  /* Panel section listing the saved styles: preview, limits, edit, delete. */
  function _buildStylesPanel(container) {
    U.clear(container);
    container.append(U.el("div", { class: "muted", style: "font-size:11.5px;margin:0 0 6px" },
      "Saved colormap styles (colormap, limits, opacity, cells/dots). Pick one on "
      + "any layer card — editing a style restyles every layer that uses it."));
    for (const p of Styles.all()) {
      const preview = U.el("span", {
        class: "cmap-preview", style: `background:${Colormaps.cssGradient(p.cmap)}`,
      });
      const lim = (p.min != null || p.max != null)
        ? `${p.min != null ? U.fmtNum(p.min, 3) : "auto"} … ${p.max != null ? U.fmtNum(p.max, 3) : "auto"}`
        : "auto";
      const edit = U.miniBtn("modify", "Edit this colormap style…", () => _presetDialog(p));
      const del = U.miniBtn("trash", "Delete", async () => {
        if (!window.confirm(`Delete colormap "${p.name}"? Layers using it keep its look as a custom style.`)) return;
        // freeze the preset's look into every assignment that references it
        for (const [key, a] of Object.entries(_assignments())) {
          if (a && a.preset === p.id) {
            _assignments()[key] = { custom: { cmap: p.cmap, mode: p.mode,
              dotSize: p.dotSize, opacity: p.opacity, min: p.min, max: p.max } };
          }
        }
        App.touchUi();
        await Styles.remove(p.id);
      });
      del.classList.add("danger-hover");
      container.append(U.el("div", { class: "cmap-row" },
        U.el("span", { class: "cmap-cat", title: p.name }, p.name),
        preview,
        U.el("span", { class: "lp-mini" }, lim),
        edit, del));
    }
    const addBtn = U.el("button", { class: "add-optional" }, "+ New colormap style");
    addBtn.addEventListener("click", () => _presetDialog(null));
    container.append(addBtn);
  }

  /* Form with the full style controls; onChange(patch) fires live. */
  function _styleForm(def, onChange) {
    const box = U.el("div");
    const cmapSel = U.el("select", {});
    for (const nm of Colormaps.names()) {
      cmapSel.append(U.el("option", {
        value: nm, selected: nm === Colormaps.baseName(def.cmap) ? "" : null,
      }, nm));
    }
    const revCb = U.el("input", { type: "checkbox" });
    revCb.checked = Colormaps.isReversed(def.cmap);
    const preview = U.el("div", {
      style: `height:10px;border-radius:5px;margin:4px 0;background:${Colormaps.cssGradient(def.cmap)}`,
    });
    const pushCmap = () => {
      const cmap = Colormaps.withReverse(cmapSel.value, revCb.checked);
      preview.style.background = Colormaps.cssGradient(cmap);
      onChange({ cmap });
    };
    cmapSel.addEventListener("change", pushCmap);
    revCb.addEventListener("change", pushCmap);

    const minIn = U.el("input", { type: "text", value: def.min != null ? def.min : "",
      placeholder: "auto", style: "width:80px" });
    const maxIn = U.el("input", { type: "text", value: def.max != null ? def.max : "",
      placeholder: "auto", style: "width:80px" });
    const pushLim = () => {
      const mn = minIn.value.trim() === "" ? null : Number(minIn.value);
      const mx = maxIn.value.trim() === "" ? null : Number(maxIn.value);
      onChange({ min: Number.isFinite(mn) ? mn : null, max: Number.isFinite(mx) ? mx : null });
    };
    for (const inp of [minIn, maxIn]) {
      inp.addEventListener("blur", pushLim);
      inp.addEventListener("keydown", (ev) => { if (ev.key === "Enter") inp.blur(); });
    }

    const opac = U.el("input", { type: "range", min: 0, max: 1, step: 0.05,
      value: def.opacity != null ? def.opacity : 0.9 });
    opac.addEventListener("input", () => onChange({ opacity: Number(opac.value) }));

    const modeSel = U.el("select", {},
      U.el("option", { value: "cells" }, "cells (pcolormesh)"),
      U.el("option", { value: "dots" }, "dots (data points)"));
    modeSel.value = def.mode === "dots" ? "dots" : "cells";
    const dotSize = U.el("input", { type: "range", min: 2, max: 14, step: 1,
      value: def.dotSize || 6 });
    const sizeRow = U.el("div", { class: "form-row" }, U.el("label", {}, "Dot size"), dotSize);
    const syncSize = () => { sizeRow.style.display = modeSel.value === "dots" ? "" : "none"; };
    modeSel.addEventListener("change", () => { onChange({ mode: modeSel.value }); syncSize(); });
    dotSize.addEventListener("input", () => onChange({ dotSize: Number(dotSize.value) }));
    syncSize();

    box.append(
      U.el("div", { class: "form-row" }, U.el("label", {}, "Colormap"), cmapSel,
        U.el("label", { class: "choice-row", style: "font-size:12px" },
          revCb, U.el("span", {}, "reversed"))),
      preview,
      U.el("div", { class: "form-row" }, U.el("label", {}, "Min / max"), minIn, maxIn),
      U.el("div", { class: "muted", style: "font-size:11px;margin:-2px 0 4px" },
        "Leave empty for auto (each layer's own data range)."),
      U.el("div", { class: "form-row" }, U.el("label", {}, "Opacity"), opac),
      U.el("div", { class: "form-row" }, U.el("label", {}, "Display"), modeSel),
      sizeRow);
    return box;
  }

  /* Create/edit a saved preset; saving restyles every layer using it. */
  function _presetDialog(existing) {
    const def = existing
      ? { ...existing }
      : { id: Styles.newId(), name: "", cmap: "viridis", mode: "cells",
          dotSize: 6, opacity: 0.9, min: null, max: null };
    const popup = Popup.open({
      title: existing ? `Edit colormap — ${existing.name}` : "New colormap style",
      width: 400,
    });
    const nameIn = U.el("input", { type: "text", value: def.name, placeholder: "e.g. Elevation (NAP)" });
    const form = _styleForm(def, (patch) => Object.assign(def, patch));
    const saveBtn = U.el("button", { class: "primary" }, existing ? "Save" : "Create");
    saveBtn.addEventListener("click", async () => {
      def.name = nameIn.value.trim();
      if (!def.name) { U.toast("Give the colormap style a name", "error"); return; }
      await Styles.save({ id: def.id, name: def.name, cmap: def.cmap, mode: def.mode,
        dotSize: def.dotSize, opacity: def.opacity, min: def.min, max: def.max });
      _applyAllStyles();
      popup.close();
    });
    const cancel = U.el("button", { class: "ghost" }, "Cancel");
    cancel.addEventListener("click", popup.close);
    popup.body.append(
      U.el("div", { class: "form-row" }, U.el("label", {}, "Name"), nameIn),
      form,
      U.el("div", { class: "btn-row", style: "justify-content:flex-end;margin-top:8px" },
        cancel, saveBtn));
  }

  function _pointColorExpr(cmap, lo, hi) {
    // MapLibre requires strictly ascending stops: a constant-valued
    // dataset (lo == hi) or an empty one would otherwise produce an
    // invalid expression and the layer would silently fail
    if (!Number.isFinite(lo) || !Number.isFinite(hi)) { lo = 0; hi = 1; }
    if (!(hi > lo)) hi = lo + 1;
    const stops = Colormaps.stops(cmap);
    const expr = ["interpolate", ["linear"], ["get", "z"]];
    stops.forEach((c, i) => {
      expr.push(lo + (hi - lo) * (i / (stops.length - 1)),
        `rgb(${c[0]},${c[1]},${c[2]})`);
    });
    return expr;
  }

  /* ================= map legend (bottom-left) =================
   * One row per distinct visible style: the colormap ramp with its
   * limits, the style name (if a saved preset) and the layers that
   * are drawn with it. Sits just above the coordinates box. */

  function _legendLabel(layer) {
    if (_isOut(layer.id)) return _outVar(layer.id);
    if (layer.id.startsWith("domain-")) return layer.id.slice("domain-".length);
    if (layer.id.startsWith("custom-")) {
      const c = _customList().find((x) => `custom-${x.id}` === layer.id);
      return (c && c.name) || "custom";
    }
    const title = layer.title || layer.id;
    if (title.length <= 26) return title;
    return _quantityOf(layer) || `${title.slice(0, 24)}…`;
  }

  function _legendEntries() {
    const rows = new Map();   // signature -> {style, name, layers: []}
    const push = (style, name, label) => {
      const sig = `${style.cmap}|${U.fmtNum(style.min, 5)}|${U.fmtNum(style.max, 5)}`;
      if (!rows.has(sig)) rows.set(sig, { style, name: name || null, layers: [] });
      const row = rows.get(sig);
      if (!row.layers.includes(label)) row.layers.push(label);
    };
    for (const [group] of _orderedGroups()) {
      if (!["output", "custom", "domain", "rawdata"].includes(group)) continue;
      for (const layer of Layers.byGroup(group)) {
        if (layer.visible === false) continue;
        const fid = _fieldIdFor(layer);
        const live = fid && FieldLayer.get(fid);
        if (_isOut(layer.id)) {
          if (!live || !live.visible) continue;
          // the live style: output limits may come from per-frame auto range
          const st = _styleOf(layer);
          push({ ...st, cmap: live.style.cmap, min: live.style.min, max: live.style.max },
            st._name, _outVar(layer.id));
          continue;
        }
        if (!live && !pointLayers.has(layer.id)) continue;
        const st = _styleOf(layer);
        push(st, st._name, _legendLabel(layer));
      }
    }
    return [...rows.values()];
  }

  function _syncLegend() {
    const wrap = document.getElementById("map-wrap");
    if (!wrap) return;
    let box = document.getElementById("map-legend");
    if (!box) {
      box = U.el("div", { id: "map-legend" });
      wrap.append(box);
    }
    U.clear(box);
    const entries = _legendEntries();
    if (!entries.length) { box.style.display = "none"; return; }
    box.style.display = "";

    const hidden = App.state.ui.colorbars === false;
    const toggle = U.el("span", {
      id: "legend-toggle", class: `eye ${hidden ? "off" : ""}`,
      title: hidden ? "Show legend" : "Hide legend",
    }, "👁");
    toggle.addEventListener("click", () => {
      App.state.ui.colorbars = hidden;   // toggled
      App.touchUi();
      _syncLegend();
    });
    if (hidden) { box.append(toggle); return; }

    for (const e of entries) {
      const layersTxt = e.layers.join(", ");
      box.append(U.el("div", { class: "legend-row" },
        U.el("span", { class: "legend-min" }, U.fmtNum(e.style.min, 3)),
        U.el("span", {
          class: "legend-ramp",
          style: `background:${Colormaps.cssGradient(e.style.cmap)}`,
        }),
        U.el("span", { class: "legend-max" }, U.fmtNum(e.style.max, 3)),
        e.name ? U.el("span", { class: "legend-name" }, e.name) : null,
        U.el("span", { class: "legend-layers", title: layersTxt }, layersTxt)));
    }
    box.append(toggle);
  }

  /* Per-layer style editor: tune a custom style live, or save it as a
   * reusable preset. Works for shown AND hidden layers (the assignment
   * persists and applies when the layer is shown). */
  function _styleEditor(layerInfo) {
    const popup = Popup.open({ title: `Style — ${layerInfo.title || layerInfo.id}`, width: 400 });
    const a = _assignments()[_styleKeyOf(layerInfo)];
    const preset = a && a.preset ? Styles.get(a.preset) : null;
    popup.body.append(U.el("div", { class: "muted", style: "font-size:11.5px;margin:0 0 6px" },
      preset
        ? `Using saved style "${preset.name}" — changing anything below turns this layer custom.`
        : "Custom style for this layer only — save it below to reuse it on other layers."));
    popup.body.append(_styleForm(_styleDefOf(layerInfo),
      (patch) => _updateCustom(layerInfo, patch)));

    // save-as-preset: inline name row (window.prompt is unreliable in webview)
    const nameIn = U.el("input", { type: "text", placeholder: "name for the new style" });
    const confirmBtn = U.el("button", { class: "primary" }, "Save");
    const nameRow = U.el("div", { class: "form-row", style: "display:none" },
      U.el("label", {}, "Name"), nameIn, confirmBtn);
    const saveAs = U.el("button", { class: "ghost" }, "Save as colormap style…");
    saveAs.addEventListener("click", () => {
      nameRow.style.display = "";
      nameIn.focus();
    });
    confirmBtn.addEventListener("click", async () => {
      const name = nameIn.value.trim();
      if (!name) { U.toast("Give the colormap style a name", "error"); return; }
      const def = _styleDefOf(layerInfo);
      const newPreset = { id: Styles.newId(), name, cmap: def.cmap,
        mode: def.mode || "cells",
        dotSize: def.dotSize != null ? def.dotSize : 6,
        opacity: def.opacity != null ? def.opacity : 0.9,
        min: def.min != null ? def.min : null,
        max: def.max != null ? def.max : null };
      await Styles.save(newPreset);
      _assignPreset(layerInfo, newPreset.id);
      U.toast(`Saved "${name}" — now selectable on every layer`, "ok");
      popup.close();
    });
    nameIn.addEventListener("keydown", (ev) => { if (ev.key === "Enter") confirmBtn.click(); });
    popup.body.append(
      U.el("div", { class: "btn-row", style: "justify-content:flex-end;margin-top:8px" }, saveAs),
      nameRow);
  }

  /* ================= output loading ================= */

  async function _loadOutput(force = false) {
    if (!App.state.project) return;
    try {
      const m = await Api.get("/api/output/meta");
      // the full source info also carries the known HPC run folder this
      // project submitted to — the other half of the source switch
      const src = await Api.get("/api/output/source").catch(() => null);
      srcInfo = src
        ? { override: !!src.override, path: src.path || null, dir: src.dir || null,
            known: src.known || null }
        : { override: !!m.override, path: m.path || null, dir: null, known: null };
      // a run that crashed before the first write leaves a netCDF with
      // an empty time dimension - treat it as "no output"
      if (!m.exists || !m.times_epoch || !m.times_epoch.length) {
        meta = null;
        Playbar.setIndexTimes(null);
        _buildPanel();
        return;
      }
      if (meta && !force && m.file === meta.file &&
          m.times.length === meta.times.length) return;
      meta = m;
      frameCache.clear(); frameRange.clear(); inflight.clear();

      const bin = await Api.binary("/api/output/mesh");
      const head = new Uint32Array(bin.buffer, 0, 2);
      const n = head[0], s = head[1];
      mesh = {
        n, s,
        x: new Float32Array(bin.buffer, 8, n * s),
        y: new Float32Array(bin.buffer, 8 + n * s * 4, n * s),
      };

      // one card per output variable; visibility persists per project
      // (zb starts visible, the rest hidden)
      const vis = _outVisibility();
      const wanted = new Set();
      for (const v of meta.variables) {
        const id = OUT_PREFIX + v.name;
        wanted.add(id);
        const shown = vis[v.name] != null ? !!vis[v.name] : v.name === "zb";
        Layers.register({ id, group: "output", title: v.name,
          subtitle: v.units || "", visible: shown });
      }
      for (const l of [...Layers.byGroup("output")]) {
        if (!wanted.has(l.id)) {
          FieldLayer.remove(`field-${l.id}`);
          Layers.unregister(l.id);
        }
      }

      Playbar.setSource("output",
        meta.times_epoch[0], meta.times_epoch[meta.times_epoch.length - 1]);
      Playbar.setIndexTimes(meta.times_epoch);
      // keep the data-availability strip's Output lane in sync with the
      // (possibly overridden) output file; debounced, so the event path
      // triggering it twice is harmless
      if (typeof Graphs !== "undefined") Graphs.refreshAvailability();

      _buildPanel();
      for (const l of Layers.byGroup("output")) {
        if (l.visible !== false) _ensureOutputLayer(_outVar(l.id));
      }
      Playbar.setTime(meta.times_epoch[0]);
      await _restoreCustomLayers();
      _applyLayerOrder();
    } catch (err) {
      console.warn("output load failed", err);
    }
  }

  /* Manual refresh of the output catalog: during a run the netCDF file
   * grows, but _loadOutput only fires on project open / run end. When
   * the same file simply gained steps, extend the time axis in place —
   * the file is append-only, so cached frames stay valid and the
   * playbar keeps its position. Anything else (first output, another
   * file, changed variables) goes through the normal full reload. */
  async function _refreshOutput(btn) {
    if (!App.state.project) return;
    if (btn) btn.disabled = true;
    try {
      const m = await Api.get("/api/output/meta").catch(() => null);
      if (m) {
        srcInfo = { ...(srcInfo || {}),
          override: !!m.override, path: m.path || null };
      }
      if (!m || !m.exists || !m.times_epoch || !m.times_epoch.length) {
        U.toast("The output file has no time steps yet", "");
        return;
      }
      const known = meta ? meta.times.length : 0;
      // full-path compare: with the same filename in another folder (an
      // output-source switch) the frame indices mean different data
      const grew = meta && m.path === meta.path
        && m.times.length >= known
        && m.times_epoch[0] === meta.times_epoch[0]
        && m.times_epoch[known - 1] === meta.times_epoch[known - 1]
        && JSON.stringify(m.variables) === JSON.stringify(meta.variables);
      if (!grew) {
        await _loadOutput(true);
        U.toast(meta ? `Output loaded — ${meta.times.length} steps` : "No output found",
          meta ? "ok" : "");
        return;
      }
      const added = m.times.length - known;
      meta = m;
      // the previously-newest step may have been read while the model
      // was still writing it — drop it from the cache and redraw
      for (const key of [...frameCache.keys()]) {
        if (key.endsWith(`|${known - 1}`)) {
          frameCache.delete(key);
          frameRange.delete(key);
        }
      }
      Playbar.setSource("output",
        m.times_epoch[0], m.times_epoch[m.times_epoch.length - 1]);
      Playbar.setIndexTimes(m.times_epoch);
      if (typeof Graphs !== "undefined") Graphs.refreshAvailability();
      _buildPanel();
      for (const l of Layers.byGroup("output")) {
        if (l.visible === false) continue;
        _outStateOf(_outVar(l.id)).bracket = null;
        _updateVarFrames(_outVar(l.id), App.state.clock.t).catch(() => {});
      }
      U.toast(added
        ? `${added} new output step${added === 1 ? "" : "s"}`
        : "No new output steps", "ok");
    } finally {
      if (btn) btn.disabled = false;
    }
  }

  /* Create (or recreate after a mesh change) the field layer of one
   * output variable, styled from its saved/default style. */
  function _ensureOutputLayer(varName) {
    const fid = `field-${OUT_PREFIX}${varName}`;
    let layer = FieldLayer.get(fid);
    if (!layer || layer.mesh !== mesh) {
      FieldLayer.remove(fid);
      const st = _styleOf({ id: OUT_PREFIX + varName });
      layer = FieldLayer.create(fid, mesh, {
        cmap: st.cmap, min: st.min, max: st.max, opacity: st.opacity,
        mode: st.mode, dotSize: st.dotSize,
      });
      _outStateOf(varName).bracket = null;
    }
    _applyLayerStyle({ id: OUT_PREFIX + varName });
    return layer;
  }

  /* ================= frames & scrubbing ================= */

  function _frameKey(varName, extraIdx, t) { return `${varName}|${extraIdx}|${t}`; }

  async function _fetchFrame(varName, extraIdx, t) {
    const key = _frameKey(varName, extraIdx, t);
    if (frameCache.has(key)) return frameCache.get(key);
    if (inflight.has(key)) return inflight.get(key);
    const promise = Api.binary(`/api/output/field?var=${varName}&t=${t}&k=${extraIdx}`)
      .then(({ buffer, headers }) => {
        const arr = new Float32Array(buffer);
        frameCache.set(key, arr);
        const range = headers.get("X-Data-Range");
        if (range) frameRange.set(key, range.split(",").map(Number));
        while (frameCache.size > FRAME_CACHE_MAX) {
          const oldest = frameCache.keys().next().value;
          frameCache.delete(oldest);
          frameRange.delete(oldest);
        }
        inflight.delete(key);
        return arr;
      })
      .catch((err) => { inflight.delete(key); throw err; });
    inflight.set(key, promise);
    return promise;
  }

  function _bracket(epoch) {
    const times = meta.times_epoch;
    if (epoch <= times[0]) return { k: 0, frac: 0 };
    if (epoch >= times[times.length - 1]) return { k: times.length - 1, frac: 0 };
    let lo = 0, hi = times.length - 1;
    while (hi - lo > 1) {
      const mid = (lo + hi) >> 1;
      if (times[mid] <= epoch) lo = mid; else hi = mid;
    }
    const span = times[hi] - times[lo];
    return { k: lo, frac: span > 0 ? (epoch - times[lo]) / span : 0 };
  }

  function _onClock(epoch) {
    if (!meta || !mesh || !Number.isFinite(epoch)) return;
    for (const l of Layers.byGroup("output")) {
      if (l.visible === false) continue;
      _updateVarFrames(_outVar(l.id), epoch).catch((err) =>
        console.warn("frame fetch failed", err.message));
    }
  }

  /* Move one output variable's field layer to the given time. */
  async function _updateVarFrames(varName, epoch) {
    if (!meta || !mesh || !Number.isFinite(epoch)) return;
    const layer = FieldLayer.get(`field-${OUT_PREFIX}${varName}`);
    if (!layer || !layer.visible) return;
    const state = _outStateOf(varName);
    const extraIdx = _extraIdxOf(varName);
    const { k, frac } = _bracket(epoch);
    const k2 = Math.min(k + 1, meta.times.length - 1);

    const cachedA = frameCache.get(_frameKey(varName, extraIdx, k));
    const cachedB = frameCache.get(_frameKey(varName, extraIdx, k2));
    if (cachedA && cachedB) {
      if (!state.bracket || state.bracket.k !== k || state.bracket.extra !== extraIdx) {
        layer.setFrames(cachedA, cachedB, frac);
        state.bracket = { k, extra: extraIdx };
        if (state.auto) _applyAutoRange(varName, k);
      } else {
        layer.setFrac(frac);
      }
      if (k2 + 1 < meta.times.length) _fetchFrame(varName, extraIdx, k2 + 1).catch(() => {});
      return;
    }
    const [a, b] = await Promise.all([
      _fetchFrame(varName, extraIdx, k), _fetchFrame(varName, extraIdx, k2)]);
    layer.setFrames(a, b, frac);
    state.bracket = { k, extra: extraIdx };
    if (state.auto) _applyAutoRange(varName, k);
    if (k2 + 1 < meta.times.length) _fetchFrame(varName, extraIdx, k2 + 1).catch(() => {});
  }

  function _applyAutoRange(varName, k) {
    const id = OUT_PREFIX + varName;
    const layer = FieldLayer.get(`field-${id}`);
    const range = frameRange.get(_frameKey(varName, _extraIdxOf(varName), k));
    if (layer && range && Number.isFinite(range[0])) {
      layer.setStyle({ min: range[0], max: range[1] });
      // the editor's "auto" placeholder + legend reflect the real range
      dataRanges.set(id, range);
      _syncLegend();
    }
  }

  /* ================= panel ================= */

  function _buildPanel() {
    let panel = document.getElementById("viewer-panel");
    if (!panel) return;
    U.keepScroll(panel);
    U.clear(panel);
    els = {};

    if (!App.state.project) {
      panel.append(U.el("div", { class: "muted" }, "Open a project first."));
      return;
    }

    // --- saved colormap styles — own section at the top ---
    const cmapSection = U.section("Colormaps", { collapsed: true });
    cmapSection.wrap.id = "viewer-cmaps-section";
    _buildStylesPanel(cmapSection.body);
    panel.append(cmapSection.wrap);

    // --- layer groups (each its own collapsible section) ---
    els.groups = U.el("div");
    panel.append(els.groups);
    _renderGroups();

    // all output controls live on the per-variable cards in the Model
    // output group (plus the probe/zoom footer under the cards)
    if (!meta) {
      panel.append(U.el("div", { class: "muted", style: "margin-top:10px" },
        "No model output yet — run a simulation in the Run tab."));
    }
  }

  /* ================= hover value readout ================= */

  function _nearestIdx(m, px, py) {
    let best = 0, bestD = Infinity;
    const xs = m.x, ys = m.y;
    for (let idx = 0; idx < xs.length; idx += 1) {
      const dx = xs[idx] - px, dy = ys[idx] - py;
      const d = dx * dx + dy * dy;
      if (d < bestD) { bestD = d; best = idx; }
    }
    return { best, dist: Math.sqrt(bestD) };
  }

  function _cellSize(m) {
    const c = Math.hypot((m.x[1] || 0) - (m.x[0] || 0), (m.y[1] || 0) - (m.y[0] || 0));
    return c > 0 ? c : Infinity;
  }

  // interpolated value of a field layer at the cursor (output uses its
  // live A/B/frac; static grids use A). null if off-grid or no data.
  function _sampleFieldLayer(fl, px, py) {
    if (!fl || !fl.mesh || !fl.mesh.x || !fl.mesh.x.length || !fl._pending || !fl._pending.a) return null;
    const { best, dist } = _nearestIdx(fl.mesh, px, py);
    if (dist > 2.5 * _cellSize(fl.mesh)) return null;
    const p = fl._pending;
    const v = (p.b && p.frac) ? p.a[best] + p.frac * (p.b[best] - p.a[best]) : p.a[best];
    return _bad(v) ? null : v;
  }

  function _samplePoints(pd, px, py) {
    let best = -1, bestD = Infinity;
    for (let i = 0; i < pd.x.length; i += 1) {
      const dx = pd.x[i] - px, dy = pd.y[i] - py;
      const d = dx * dx + dy * dy;
      if (d < bestD) { bestD = d; best = i; }
    }
    if (best < 0 || Math.sqrt(bestD) > pd.radius) return null;
    return _bad(pd.z[best]) ? null : pd.z[best];
  }

  function _labelFor(layer, fid) {
    if (_isOut(layer.id)) {
      const varName = _outVar(layer.id);
      const info = (meta && meta.variables.find((v) => v.name === varName)) || {};
      return { label: varName, unit: info.units || "" };
    }
    return { label: "z", unit: "m" };
  }

  /* Value of the TOP visible layer under the cursor — walks the viewer's
   * group order (top first) and samples whichever kind the layer is:
   * output (interpolated), interpolated .grd, imported raster, or point
   * samples. Returns e.g. "zb = 5.230 m" or null when nothing is hit. */
  function _hoverSample(xy) {
    const px = xy[0], py = xy[1];
    for (const [group] of _orderedGroups()) {
      if (!["output", "custom", "domain", "rawdata"].includes(group)) continue;
      for (const layer of Layers.byGroup(group)) {
        if (layer.visible === false) continue;
        const fid = _fieldIdFor(layer);
        const fl = fid && FieldLayer.get(fid);
        let v = null;
        if (fl) v = _sampleFieldLayer(fl, px, py);
        else if (pointData.has(layer.id)) v = _samplePoints(pointData.get(layer.id), px, py);
        if (v != null) {
          const lu = _labelFor(layer, fid);
          const dec = Math.abs(v) < 100 ? 3 : 1;
          return `${lu.label} = ${v.toFixed(dec)}${lu.unit ? " " + lu.unit : ""}`;
        }
      }
    }
    return null;
  }

  /* ================= probe ================= */

  function _armProbe(button) {
    const map = MapView.instance();
    map.getCanvas().style.cursor = "crosshair";
    button.disabled = true;
    map.once("click", async (ev) => {
      map.getCanvas().style.cursor = "";
      button.disabled = false;
      if (!mesh || !meta) return;
      // probe every visible output variable at the clicked cell
      const vars = Layers.byGroup("output")
        .filter((l) => l.visible !== false).map((l) => _outVar(l.id));
      if (!vars.length) {
        U.toast("Show at least one output variable first", "error");
        return;
      }
      const [px, py] = CRS.fromLngLat(ev.lngLat);
      let best = 0, bestDist = Infinity;
      for (let idx = 0; idx < mesh.x.length; idx += 1) {
        const dx = mesh.x[idx] - px, dy = mesh.y[idx] - py;
        const d = dx * dx + dy * dy;
        if (d < bestDist) { bestDist = d; best = idx; }
      }
      const j = Math.floor(best / mesh.s);
      const i = best % mesh.s;
      try {
        for (const [vi, varName] of vars.entries()) {
          const res = await Api.get(
            `/api/output/series?var=${varName}&j=${j}&i=${i}&k=${_extraIdxOf(varName)}`);
          const info = meta.variables.find((v) => v.name === varName) || {};
          Graphs.registerSource(`probe-${varName}-${j}-${i}`, {
            group: "Output",
            label: `${varName} @ cell (${j},${i})`,
            unit: info.units || "",
            data: [res.t_epoch, res.values],
            range: res.t_epoch.length
              ? [res.t_epoch[0], res.t_epoch[res.t_epoch.length - 1]] : null,
          }, { select: vi === 0 });
        }
        MapView.setLabel(`probe-${j}-${i}`,
          [mesh.x[best], mesh.y[best]], `(${j},${i})`, "boundary-lateral");
      } catch (err) {
        U.toast(err.message, "error");
      }
    });
  }

  /* ================= domain / raw layer toggling ================= */

  async function _onLayerVisibility(layerInfo) {
    if (_isOut(layerInfo.id)) {
      const varName = _outVar(layerInfo.id);
      _outVisibility()[varName] = layerInfo.visible !== false;
      App.touchUi();
      if (layerInfo.visible === false) {
        FieldLayer.remove(`field-${layerInfo.id}`);
      } else if (mesh) {
        _ensureOutputLayer(varName);
        await _updateVarFrames(varName, App.state.clock.t).catch(() => {});
        _applyLayerOrder();
      }
      _syncLegend();
      return;
    }
    if (layerInfo.id.startsWith("raw-") && layerInfo.entry) {
      await _toggleRawLayer(layerInfo);
      _applyLayerOrder();
      _syncLegend();
    }
    if (layerInfo.id.startsWith("domain-")) {
      await _toggleDomainLayer(layerInfo);
      _applyLayerOrder();
      _syncLegend();
    }
  }

  async function _toggleRawLayer(layerInfo) {
    const id = `field-${layerInfo.id}`;
    if (!layerInfo.visible) {
      FieldLayer.remove(id);
      MapView.removeLayerAndSource(`pts-${layerInfo.id}`);
      pointLayers.delete(layerInfo.id);
      pointData.delete(layerInfo.id);
      return;
    }
    _setLoading(layerInfo.id, true);
    try {
      const entry = layerInfo.entry;
      if (entry.kind === "points") {
        const res = await Api.get(`/api/domain/rawfield?id=${entry.id}`);
        const features = res.x.map((x, i) => ({
          type: "Feature",
          properties: { z: res.z[i] },
          geometry: { type: "Point", coordinates: CRS.toLngLat([x, res.y[i]]) },
        }));
        let zmin = Infinity, zmax = -Infinity;
        for (const z of res.z) { if (z < zmin) zmin = z; if (z > zmax) zmax = z; }
        dataRanges.set(layerInfo.id, [zmin, zmax]);
        // cache the samples so hover/transect can read z at a location
        let minx = Infinity, miny = Infinity, maxx = -Infinity, maxy = -Infinity;
        for (let i = 0; i < res.x.length; i += 1) {
          if (res.x[i] < minx) minx = res.x[i]; if (res.x[i] > maxx) maxx = res.x[i];
          if (res.y[i] < miny) miny = res.y[i]; if (res.y[i] > maxy) maxy = res.y[i];
        }
        const spacing = Math.sqrt(Math.max(1, (maxx - minx) * (maxy - miny)) / Math.max(1, res.x.length));
        pointData.set(layerInfo.id, { x: res.x, y: res.y, z: res.z, radius: 2.5 * spacing });
        pointLayers.add(layerInfo.id);
        MapView.upsertGeojson(`pts-${layerInfo.id}`, { type: "FeatureCollection", features });
        const style = _styleOf(layerInfo);
        MapView.ensureLayer({
          id: `pts-${layerInfo.id}`, type: "circle", source: `pts-${layerInfo.id}`,
          paint: {
            "circle-radius": Math.max(1, style.dotSize / 2.5),
            "circle-color": _pointColorExpr(style.cmap, style.min, style.max),
            "circle-opacity": style.opacity,
          },
        });
      } else {
        const { buffer, headers } = await Api.binary(`/api/domain/rawfield?id=${entry.id}`);
        const parsed = FieldLayer.parseGridfield(buffer, headers);
        dataRanges.set(layerInfo.id, parsed.range);
        const st = _styleOf(layerInfo);
        const layer = FieldLayer.create(id, parsed.mesh, {
          cmap: st.cmap, min: st.min, max: st.max, opacity: st.opacity,
          mode: st.mode, dotSize: st.dotSize,
        });
        layer.setFrames(parsed.values);
      }
    } catch (err) {
      U.toast(`Layer failed: ${err.message}`, "error");
      layerInfo.visible = false;
    } finally {
      _setLoading(layerInfo.id, false);
    }
  }

  async function _toggleDomainLayer(layerInfo) {
    const id = `field-${layerInfo.id}`;
    const target = layerInfo.id.replace("domain-", "");
    if (!layerInfo.visible) {
      FieldLayer.remove(id);
      _syncLegend();
      return;
    }
    _setLoading(layerInfo.id, true);
    try {
      const k = layerInfo.speciesIdx || 0;
      const { buffer, headers } = await Api.binary(`/api/domain/gridfield?target=${target}&k=${k}`);
      const nSpecies = Number(headers.get("X-Species") || 1);
      Layers.register({ id: layerInfo.id, group: layerInfo.group, title: layerInfo.title,
        species: nSpecies, speciesIdx: k });
      const parsed = FieldLayer.parseGridfield(buffer, headers);
      dataRanges.set(layerInfo.id, parsed.range);
      const st = _styleOf(layerInfo);
      const layer = FieldLayer.create(id, parsed.mesh, {
        cmap: st.cmap, min: st.min, max: st.max, opacity: st.opacity,
        mode: st.mode, dotSize: st.dotSize,
      });
      layer.setFrames(parsed.values);
    } catch (err) {
      U.toast(`Layer failed: ${err.message}`, "error");
      // roll the eye back so the tree reflects reality
      layerInfo.visible = false;
    } finally {
      _setLoading(layerInfo.id, false);
    }
  }

  /* ================= dataset catalog + on-demand sampling ==============
   * Sources for the Transect tab and Custom layers, whether or not shown
   * on the map. Keys: "out:<var>" (output variable at the current time),
   * "tgt:<target>" (interpolated .grd), "ent:<entryId>" (raw/imported).
   * Field data is fetched and cached on demand. */
  function hasOutput() { return !!meta; }

  const sourceCache = new Map();   // key -> {mesh,data} | {points}
  const _entKind = {};             // entry id -> "points" | "raster"
  const _srcLabels = {};           // key -> friendly label (for formulas)

  async function datasetCatalog() {
    const out = [];
    if (meta) {
      for (const v of meta.variables) {
        const key = `out:${v.name}`, label = `output: ${v.name}`;
        _srcLabels[key] = label; out.push({ key, label, kind: "raster" });
      }
    }
    let ov = null;
    try { ov = await Api.get("/api/domain"); } catch (e) { ov = null; }
    if (ov) {
      for (const [t, info] of Object.entries(ov.targets || {})) {
        if (info && (info.exists || info.has_draft)) {
          const key = `tgt:${t}`, label = `model file: ${t}`;
          _srcLabels[key] = label; out.push({ key, label, kind: "raster" });
        }
      }
      for (const e of (ov.entries || [])) {
        _entKind[e.id] = e.kind === "points" ? "points" : "raster";
        const key = `ent:${e.id}`, label = e.label || e.id;
        _srcLabels[key] = label; out.push({ key, label, kind: _entKind[e.id] });
      }
    }
    return out;
  }

  function _pointsFrom(res) {
    let minx = Infinity, miny = Infinity, maxx = -Infinity, maxy = -Infinity;
    for (let i = 0; i < res.x.length; i += 1) {
      if (res.x[i] < minx) minx = res.x[i]; if (res.x[i] > maxx) maxx = res.x[i];
      if (res.y[i] < miny) miny = res.y[i]; if (res.y[i] > maxy) maxy = res.y[i];
    }
    const spacing = Math.sqrt(Math.max(1, (maxx - minx) * (maxy - miny)) / Math.max(1, res.x.length));
    return { x: res.x, y: res.y, z: res.z, radius: 2.5 * spacing };
  }

  // resolve a key to {mesh,data} or {points}, fetching + caching if needed
  async function _ensureSource(key) {
    if (key.startsWith("out:")) {
      if (!mesh) return null;
      const data = await _operandArray(key.slice(4), "current");
      return data ? { mesh, data } : null;
    }
    if (key.startsWith("tgt:")) {
      const target = key.slice(4);
      const live = FieldLayer.get(`field-domain-${target}`);
      if (live && live._pending && live._pending.a) return { mesh: live.mesh, data: live._pending.a };
      if (sourceCache.has(key)) return sourceCache.get(key);
      try {
        const { buffer, headers } = await Api.binary(`/api/domain/gridfield?target=${target}&k=0`);
        const parsed = FieldLayer.parseGridfield(buffer, headers);
        const r = { mesh: parsed.mesh, data: parsed.values };
        sourceCache.set(key, r); return r;
      } catch (e) { return null; }
    }
    if (key.startsWith("ent:")) {
      const id = key.slice(4);
      const live = FieldLayer.get(`field-raw-${id}`);
      if (live && live._pending && live._pending.a) return { mesh: live.mesh, data: live._pending.a };
      if (pointData.has(`raw-${id}`)) return { points: pointData.get(`raw-${id}`) };
      if (sourceCache.has(key)) return sourceCache.get(key);
      try {
        if (_entKind[id] === "points") {
          const res = await Api.get(`/api/domain/rawfield?id=${id}`);
          const r = { points: _pointsFrom(res) };
          sourceCache.set(key, r); return r;
        }
        const { buffer, headers } = await Api.binary(`/api/domain/rawfield?id=${id}`);
        const parsed = FieldLayer.parseGridfield(buffer, headers);
        const r = { mesh: parsed.mesh, data: parsed.values };
        sourceCache.set(key, r); return r;
      } catch (e) { return null; }
    }
    return null;
  }

  function _polylinePoints(coords, n) {
    const segs = [];
    let total = 0;
    for (let i = 1; i < coords.length; i += 1) {
      const dx = coords[i][0] - coords[i - 1][0], dy = coords[i][1] - coords[i - 1][1];
      const L = Math.hypot(dx, dy);
      segs.push({ x0: coords[i - 1][0], y0: coords[i - 1][1], dx, dy, L, acc: total });
      total += L;
    }
    if (total <= 0) return null;
    const dist = new Array(n), px = new Array(n), py = new Array(n);
    for (let k = 0; k < n; k += 1) {
      const d = total * k / (n - 1);
      dist[k] = d;
      let s = segs[segs.length - 1];
      for (const sg of segs) { if (d <= sg.acc + sg.L) { s = sg; break; } }
      const f = s.L > 0 ? (d - s.acc) / s.L : 0;
      px[k] = s.x0 + s.dx * f; py[k] = s.y0 + s.dy * f;
    }
    return { dist, px, py };
  }

  async function sampleTransect(coords, keys, nPoints = 240) {
    if (!coords || coords.length < 2 || !keys || !keys.length) return null;
    const line = _polylinePoints(coords, nPoints);
    if (!line) return null;
    const values = {};
    for (const key of keys) {
      const src = await _ensureSource(key);
      if (!src) continue;
      const out = new Array(nPoints);
      if (src.mesh) {
        for (let k = 0; k < nPoints; k += 1) {
          const best = _nearestIdx(src.mesh, line.px[k], line.py[k]).best;
          const v = src.data[best];
          out[k] = _bad(v) ? null : v;
        }
      } else {
        for (let k = 0; k < nPoints; k += 1) out[k] = _samplePoints(src.points, line.px[k], line.py[k]);
      }
      values[key] = out;
    }
    return { dist: line.dist, values };
  }

  return { init, sampleTransect, datasetCatalog, hasOutput };
})();
