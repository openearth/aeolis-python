/* Viewer tab: the single place to manage everything shown on the map.
 *
 * - One collapsible section per layer group (Model output /
 *   Interpolated .grd / Sample data / Grid / Objects / Background),
 *   each layer in its own draggable card - same look as the Domain tab.
 * - Stale .grd files (grid changed) cannot be displayed: the eye is
 *   disabled and an orange warning icon explains why on hover.
 * - Layers of the same quantity (e.g. several bathymetry sources)
 *   share one colormap and range; colorbars overlay the map.
 * - netCDF output: variable/statistic selection and instant time
 *   scrubbing (client LRU + prefetch + GPU frame interpolation).
 */
"use strict";

const ViewerTab = (() => {

  let meta = null;              // /api/output/meta
  let mesh = null;              // {x, y, n, s}
  let variable = null;
  let extraIdx = "0";
  let frameCache = new Map();
  let frameRange = new Map();
  let inflight = new Map();
  let currentBracket = null;
  let autoRange = true;
  let els = {};

  const FRAME_CACHE_MAX = 60;
  const OUTPUT_LAYER = "field-output";

  const GROUPS = [
    ["output", "Model output"],
    ["custom", "Custom layers"],
    ["domain", "Interpolated (.grd)"],
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
    App.on("clock-tick", _onClock);
    App.on("layer-visibility", _onLayerVisibility);
    App.on("layers", () => _renderGroups());
    App.on("objects", () => _renderGroups());
    App.on("basemap", () => _renderGroups());
    App.on("layer-order", _applyLayerOrder);
    App.on("clock-tick", _onCustomClock);
    MapView.setHoverSampler(_hoverSample);
    _buildPanel();
  }

  function _reset() {
    meta = null; mesh = null; variable = null;
    frameCache.clear(); frameRange.clear(); inflight.clear();
    dataRanges.clear(); pointLayers.clear(); pointData.clear();
    currentBracket = null;
    FieldLayer.remove(OUTPUT_LAYER);
    for (const fid of customIds) FieldLayer.remove(fid);
    customIds.clear();
    Playbar.removeSource("output");
    _syncColorbars();
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

  /* Display order of the groups: user-draggable, persisted. */
  function _orderedGroups() {
    const keys = GROUPS.map(([g]) => g);
    const saved = App.state.ui.viewerGroupOrder || [];
    const order = [...saved.filter((g) => keys.includes(g)),
      ...keys.filter((g) => !saved.includes(g))];
    return order.map((g) => GROUPS.find(([key]) => key === g));
  }

  /* Eye button in a section header: show/hide the whole group. */
  function _groupEye(group, layers) {
    let anyVisible;
    if (group === "objects") {
      anyVisible = App.state.objects.some((o) => o.visible);
    } else {
      anyVisible = layers.some((l) => l.visible);
    }
    const eye = U.el("span", {
      class: `eye group-eye ${anyVisible ? "" : "off"}`,
      title: anyVisible ? "Hide the whole group" : "Show the whole group",
    }, "👁");
    eye.addEventListener("click", (ev) => {
      ev.stopPropagation();
      if (group === "objects") {
        for (const obj of App.state.objects) {
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

  function _renderGroups() {
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
      App.state.ui.viewerGroupOrder = order;
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

    // styling for field layers and point (sample) layers. The category picker
    // is shown for every colormapped layer *type* (output / domain-* / raw-*),
    // even when the layer is currently hidden and its map object destroyed, so
    // the colour choice stays reachable; it applies when the layer is re-shown.
    const fieldId = _fieldIdFor(layer);
    const isColormapped = Boolean(fieldId) || pointLayers.has(layer.id);
    if (isColormapped) {
      // per-layer VARIABLE picker: choose which category's colormap this layer
      // uses (colours & limits are set once per category in the Colormaps group)
      const curCat = _categoryOfLayer(layer);
      const sel = U.el("select", { class: "cmap-mini", title: "Which variable colormap this layer uses" });
      for (const [ck, clabel] of CMAP_CATEGORIES) {
        sel.append(U.el("option", { value: ck, selected: ck === curCat ? "" : null }, clabel));
      }
      sel.addEventListener("click", (ev) => ev.stopPropagation());
      sel.addEventListener("change", () => _setLayerCategory(layer, sel.value));
      card.append(sel);
      // the full style editor needs the live map object, so keep it gated
      const live = (fieldId && FieldLayer.get(fieldId)) || pointLayers.has(layer.id);
      if (live) card.append(U.miniBtn("gear", "More layer options…", () => _styleEditor(layer)));
    }
    return card;
  }

  /* The colormap a layer is currently drawn with. */
  function _currentCmapOf(layer) {
    if (layer.id === OUTPUT_LAYER) {
      const l = FieldLayer.get(OUTPUT_LAYER);
      return l ? l.style.cmap : _defaultCmapForName(variable);
    }
    if (layer.id.startsWith("raw-")) return _effectiveRawStyle(layer).cmap;
    const fl = FieldLayer.get(_fieldIdFor(layer));
    return fl ? fl.style.cmap : _defaultCmapForName(_quantityOf(layer));
  }

  /* Apply a full style patch ({cmap?, min?, max?, opacity?}) to one layer. */
  function _applyStyleToLayer(layer, upd) {
    const map = MapView.instance();
    if (layer.id === OUTPUT_LAYER) {
      const l = FieldLayer.get(OUTPUT_LAYER);
      if (l) l.setStyle(upd);
    } else if (pointLayers.has(layer.id)) {
      const st = _ownPointStyle(layer);
      Object.assign(st, upd);
      _setStyleLink(layer, "own");
      if (map.getLayer(`pts-${layer.id}`)) {
        map.setPaintProperty(`pts-${layer.id}`, "circle-color",
          _pointColorExpr(st.cmap, st.min, st.max));
      }
    } else {
      const fl = FieldLayer.get(_fieldIdFor(layer));
      const target = _quantityOf(layer);
      if (fl) {
        fl.setStyle(upd);
        if (layer.id.startsWith("domain-") && target) _rememberTargetStyle(target, fl.style);
      } else if (layer.id.startsWith("domain-") && target) {
        // layer hidden (no live object): persist to the remembered target
        // style so the change is honoured when the layer is next shown
        _rememberTargetStyle(target, Object.assign({}, _targetStyles()[target], upd));
      }
    }
    _syncColorbars();
  }

  /* Point a layer at a variable category; it adopts that category's colormap
   * and z-limits so every source in the category shares one scale. */
  function _setLayerCategory(layer, cat) {
    _catOverride()[_layerCatKey(layer)] = cat;
    App.touchUi();
    const upd = { cmap: _cmapForCategory(cat) };
    const clim = _categoryClim(cat);
    if (clim.min != null) upd.min = clim.min;
    if (clim.max != null) upd.max = clim.max;
    _applyStyleToLayer(layer, upd);
  }

  /* Set a single layer's colormap (its own override). */
  function _setLayerCmap(layer, cmap) {
    const map = MapView.instance();
    if (layer.id === OUTPUT_LAYER) {
      const l = FieldLayer.get(OUTPUT_LAYER);
      if (l) l.setStyle({ cmap });
      _syncColorbars();
      return;
    }
    if (pointLayers.has(layer.id)) {
      const st = _ownPointStyle(layer);
      st.cmap = cmap;
      _setStyleLink(layer, "own");   // detach from any linked target
      if (map.getLayer(`pts-${layer.id}`)) {
        map.setPaintProperty(`pts-${layer.id}`, "circle-color",
          _pointColorExpr(cmap, st.min, st.max));
      }
      _syncColorbars();
      return;
    }
    const fl = FieldLayer.get(_fieldIdFor(layer));
    if (fl) {
      fl.setStyle({ cmap });
      const target = _quantityOf(layer);
      if (layer.id.startsWith("domain-") && target) _rememberTargetStyle(target, fl.style);
      _syncColorbars();
    }
  }

  /* drag & drop within one list (shared, de-lagged, insertion-line UI) */
  function _wireCardDrag(list, onDrop) {
    U.wireSortable(list, onDrop);
  }

  function _fieldIdFor(layer) {
    if (layer.id === OUTPUT_LAYER) return OUTPUT_LAYER;
    if (layer.id.startsWith("domain-") || layer.id.startsWith("raw-")) return `field-${layer.id}`;
    if (layer.id.startsWith("custom-")) return `field-${layer.id}`;
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
      const zoom = U.miniBtn("eye", "Zoom to", () => {
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

      const zoom = U.miniBtn("eye", "Zoom to", () => {
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

  async function _fetchFrameFor(varName, t) {
    const key = `${varName}|0|${t}`;
    if (frameCache.has(key)) return frameCache.get(key);
    if (inflight.has(key)) return inflight.get(key);
    const promise = Api.binary(`/api/output/field?var=${varName}&t=${t}&k=0`)
      .then(({ buffer }) => {
        const arr = new Float32Array(buffer);
        frameCache.set(key, arr);
        while (frameCache.size > FRAME_CACHE_MAX) {
          const oldest = frameCache.keys().next().value;
          frameCache.delete(oldest); frameRange.delete(oldest);
        }
        inflight.delete(key);
        return arr;
      })
      .catch((err) => { inflight.delete(key); throw err; });
    inflight.set(key, promise);
    return promise;
  }

  // resolve one operand: "current" interpolates the two bracketing frames
  async function _operandArray(varName, timeSel) {
    if (!meta || !varName) return null;
    if (timeSel === "current" || timeSel == null) {
      const { k, frac } = _bracket(App.state.clock.t);
      const k2 = Math.min(k + 1, meta.times.length - 1);
      const [A, B] = await Promise.all([_fetchFrameFor(varName, k), _fetchFrameFor(varName, k2)]);
      if (!frac || A === B) return A;
      const out = new Float32Array(A.length);
      for (let i = 0; i < A.length; i += 1) {
        out[i] = (_bad(A[i]) || _bad(B[i])) ? SENTINEL : A[i] + frac * (B[i] - A[i]);
      }
      return out;
    }
    const t = Math.max(0, Math.min(Number(timeSel) || 0, meta.times.length - 1));
    return _fetchFrameFor(varName, t);
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
    const eq = (a, b) => Math.abs(a - b) < 1e-6 * (1 + Math.abs(a));
    return eq(m1.x[0], m2.x[0]) && eq(m1.y[0], m2.y[0])
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
      layer = FieldLayer.create(fid, res.mesh, {
        cmap: c.cmap || "RdBu", min: c.min != null ? c.min : -1, max: c.max != null ? c.max : 1, opacity: 1,
      });
      customIds.add(fid);
    }
    const patch = { cmap: c.cmap || "RdBu" };
    if (c.min != null && c.max != null) { patch.min = c.min; patch.max = c.max; }
    else {
      let lo = Infinity, hi = -Infinity;
      for (let i = 0; i < res.out.length; i += 1) {
        const v = res.out[i];
        if (!_bad(v) && Number.isFinite(v)) { if (v < lo) lo = v; if (v > hi) hi = v; }
      }
      if (lo <= hi) {
        if (c.op === "sub" || c.op === "offset") { const m = Math.max(Math.abs(lo), Math.abs(hi)) || 1; patch.min = -m; patch.max = m; }
        else { patch.min = lo; patch.max = hi; }
      }
    }
    layer.setStyle(patch);
    layer.setFrames(res.out);
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
    if (c.visible) await _recomputeCustom(c);
    else { FieldLayer.remove(`field-custom-${c.id}`); customIds.delete(`field-custom-${c.id}`); }
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
    addBtn.title = "Combine output variables / timesteps or any interpolated & raster grids";
    body.append(addBtn);
  }

  function _customCard(c) {
    const eye = U.el("span", { class: `eye ${c.visible ? "" : "off"}`, title: "Show/hide" }, "\u{1F441}");
    eye.addEventListener("click", () => _toggleCustom(c));
    const cmapSel = U.el("select", { class: "cmap-mini", title: "Colormap" });
    for (const name of Colormaps.names()) {
      cmapSel.append(U.el("option", { value: name, selected: name === c.cmap ? "" : null }, name));
    }
    cmapSel.addEventListener("click", (ev) => ev.stopPropagation());
    cmapSel.addEventListener("change", () => { c.cmap = cmapSel.value; App.touchUi(); if (c.visible) _recomputeCustom(c); });
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
      cmapSel, edit, del);
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
      id: `${Date.now()}`, name: "", op: "sub", cmap: "RdBu",
      aKey: keys[0], aTime: "current", bKey: keys[1] || keys[0], bTime: "current",
      scalar: 1, min: null, max: null,
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
    const cmap = U.el("select", {});
    for (const n of Colormaps.names()) cmap.append(U.el("option", { value: n, selected: n === c.cmap ? "" : null }, n));
    const minI = U.el("input", { type: "number", step: "any", value: c.min != null ? c.min : "", placeholder: "auto" });
    const maxI = U.el("input", { type: "number", step: "any", value: c.max != null ? c.max : "", placeholder: "auto" });

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
        id: c.id, op: opSel.value, cmap: cmap.value,
        aKey: aSel.value, aTime: aTime.value, bKey: bSel.value, bTime: bTime.value,
        scalar: Number(scalar.value) || 0,
        resample: resampleCb.checked,
        min: minI.value === "" ? null : Number(minI.value),
        max: maxI.value === "" ? null : Number(maxI.value),
        visible: existing ? c.visible : true,
      };
      def.name = name.value.trim() || _customFormula(def);
      const cl = _customList();
      const idx = cl.findIndex((x) => x.id === def.id);
      if (idx >= 0) cl[idx] = def; else cl.push(def);
      App.touchUi();
      FieldLayer.remove(`field-custom-${def.id}`);
      customIds.delete(`field-custom-${def.id}`);
      if (def.visible) await _recomputeCustom(def);
      _syncCustomRegistry();
      _applyLayerOrder();
      _renderGroups();
      popup.close();
    });
    const cancelBtn = U.el("button", { class: "ghost" }, "Cancel");
    cancelBtn.addEventListener("click", popup.close);

    popup.body.append(
      U.el("div", { class: "form-row" }, U.el("label", {}, "Name"), name),
      U.el("div", { class: "form-row" }, U.el("label", {}, "Compute"), opSel),
      rowA, rowB, rowK, rowR,
      U.el("div", { class: "form-row" }, U.el("label", {}, "Colormap"), cmap),
      U.el("div", { class: "form-row" }, U.el("label", {}, "z-limits"), minI, maxI),
      U.el("div", { class: "muted", style: "font-size:11.5px" },
        "Combine output timesteps or any interpolated / raster grids. If the two grids differ, tick "
        + "“resample”. Leave z-limits blank to auto-scale."),
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

  /* ================= shared colormaps + colorbars ================= */

  /* Which physical quantity a layer shows (bed, ne, veg, ...). Raw
   * downloads are elevation samples -> "bed"; layers duplicated from a
   * .grd keep that target's quantity (label "<target> (file) → samples"). */
  function _quantityOf(layer) {
    if (layer.id === OUTPUT_LAYER) return null;   // output has its own controls
    if (layer.id.startsWith("domain-")) return layer.id.slice("domain-".length);
    if (layer.id.startsWith("raw-")) {
      const entry = layer.entry || {};
      if (entry.source === "converted") {
        const m = (entry.label || "").match(/^([A-Za-z_]\w*)\s*\(/);
        if (m) return m[1];
      }
      return "bed";
    }
    return null;
  }

  /* ---- colormap-per-variable-category ----
   * The user picks one colormap per physical quantity category; every
   * layer of that category defaults to it. Choices persist per project
   * (ui.varCmaps). A layer can still override its own via the dropdown. */
  const CMAP_CATEGORIES = [
    ["elevation", "Elevation", "topo_dutch"],
    ["bed_change", "Bed level change", "RdBu"],
    ["veg_density", "Vegetation density", "Greens"],
    ["veg_height", "Vegetation height", "Greens"],
    ["mask", "Masks", "gray"],
    ["wind_speed", "Wind speed", "viridis"],
    ["shear_stress", "Shear stress", "plasma"],
    ["shear_velocity", "Shear velocity", "plasma"],
    ["sed_conc", "Sediment concentration", "sand"],
    ["sed_transport", "Sediment transport", "sand"],
  ];
  const CMAP_DEFAULTS = Object.fromEntries(CMAP_CATEGORIES.map(([k, , c]) => [k, c]));

  /* Map a quantity/target/output-variable name to a category (or null). */
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

  /* Per-layer category overrides: a layer's colour is decided by the variable
   * CATEGORY it belongs to (Elevation, Wind speed, …), inferred from its name
   * but overridable per layer via the layer dropdown. All colours/limits then
   * live at the category level (the Colormaps group), so every source in a
   * category shares one scale. */
  function _catOverride() {
    App.state.ui.layerCategory = App.state.ui.layerCategory || {};
    return App.state.ui.layerCategory;
  }
  function _layerCatKey(layer) {
    if (layer.id === OUTPUT_LAYER) return `output:${variable}`;
    if (layer.id.startsWith("raw-")) return `raw:${layer.id}`;
    return `target:${_quantityOf(layer)}`;
  }
  function _layerName(layer) {
    return layer.id === OUTPUT_LAYER ? variable : _quantityOf(layer);
  }
  function _categoryOfKey(key, name) {
    return _catOverride()[key] || _categoryOfName(name);
  }
  function _categoryOfLayer(layer) {
    return _categoryOfKey(_layerCatKey(layer), _layerName(layer));
  }

  function _varCmaps() {
    App.state.ui.varCmaps = App.state.ui.varCmaps || {};
    return App.state.ui.varCmaps;
  }
  function _cmapForCategory(cat) {
    if (!cat) return null;
    return _varCmaps()[cat] || CMAP_DEFAULTS[cat] || "viridis";
  }
  /* Default colormap for a quantity/variable via its category. */
  function _defaultCmapForName(name, fallback = "viridis") {
    return _cmapForCategory(_categoryOfName(name)) || fallback;
  }

  function _categoryClim(cat) {
    return (App.state.ui.varClim && App.state.ui.varClim[cat]) || {};
  }

  /* Apply a style patch ({cmap?, min?, max?}) to every currently loaded
   * layer of a category (interpolated sets, sample layers, output). */
  function _applyCategoryStyle(cat, patch) {
    if (patch.cmap != null) _varCmaps()[cat] = patch.cmap;
    if (patch.min != null || patch.max != null) {
      App.state.ui.varClim = App.state.ui.varClim || {};
      App.state.ui.varClim[cat] = Object.assign({}, App.state.ui.varClim[cat],
        patch.min != null ? { min: patch.min } : {},
        patch.max != null ? { max: patch.max } : {});
    }
    App.touchUi();
    const map = MapView.instance();
    const upd = {};
    if (patch.cmap != null) upd.cmap = patch.cmap;
    if (patch.min != null) upd.min = patch.min;
    if (patch.max != null) upd.max = patch.max;

    for (const [target, style] of Object.entries(_targetStyles())) {
      if (_categoryOfKey(`target:${target}`, target) !== cat) continue;
      Object.assign(style, upd);
      const live = FieldLayer.get(`field-domain-${target}`);
      if (live) live.setStyle(upd);
    }
    if (variable && _categoryOfKey(`output:${variable}`, variable) === cat) {
      const l = FieldLayer.get(OUTPUT_LAYER);
      if (l) l.setStyle(upd);
    }
    for (const layer of App.state.layers) {
      if (!layer.id.startsWith("raw-")) continue;
      if (_categoryOfKey(`raw:${layer.id}`, _quantityOf(layer)) !== cat) continue;
      const own = _ownPointStyle(layer);
      Object.assign(own, upd);
      if (pointLayers.has(layer.id) && map.getLayer(`pts-${layer.id}`)) {
        map.setPaintProperty(`pts-${layer.id}`, "circle-color",
          _pointColorExpr(own.cmap, own.min, own.max));
      }
    }
    _syncSharedStyles();   // also refreshes colorbars
  }

  function _applyCategoryCmap(cat, cmap) { _applyCategoryStyle(cat, { cmap }); }

  /* Popup to set z-limits + reverse for a whole variable category. */
  function _categoryStyleDialog(cat, label, rebuild) {
    const full = _cmapForCategory(cat);
    const clim = _categoryClim(cat);
    const popup = Popup.open({ title: `${label} — colormap settings`, width: 360 });
    const minIn = U.el("input", { type: "text", value: clim.min != null ? clim.min : "",
      placeholder: "auto", style: "width:90px" });
    const maxIn = U.el("input", { type: "text", value: clim.max != null ? clim.max : "",
      placeholder: "auto", style: "width:90px" });
    const revCb = U.el("input", { type: "checkbox", id: "cmap-rev" });
    revCb.checked = Colormaps.isReversed(full);
    const applyBtn = U.el("button", { class: "primary" }, "Apply");
    applyBtn.addEventListener("click", () => {
      const patch = { cmap: Colormaps.withReverse(Colormaps.baseName(full), revCb.checked) };
      const mn = Number(minIn.value), mx = Number(maxIn.value);
      if (minIn.value.trim() !== "" && Number.isFinite(mn)) patch.min = mn;
      if (maxIn.value.trim() !== "" && Number.isFinite(mx)) patch.max = mx;
      _applyCategoryStyle(cat, patch);
      popup.close();
      if (rebuild) rebuild();
    });
    popup.body.append(
      U.el("div", { class: "form-row" }, U.el("label", {}, "min / max"), minIn, maxIn),
      U.el("div", { class: "choice-row" }, revCb, U.el("label", { for: "cmap-rev" }, "reverse colormap")),
      U.el("div", { class: "muted", style: "font-size:11.5px" },
        "Applies to every layer of this variable; leave min/max empty to keep the current range."),
      U.el("div", { class: "btn-row", style: "justify-content:flex-end" }, applyBtn));
  }

  /* Build the category → colormap controls into a container (used both in
   * the Viewer panel and in the popup opened from the map colorbars). */
  function _buildColormapControls(container) {
    U.clear(container);
    container.append(U.el("div", { class: "muted", style: "font-size:11.5px;margin:0 0 6px" },
      "Set the colormap per variable; the gear sets z-limits & reverse. Every "
      + "layer of that variable uses it (override an individual layer with its own dropdown)."));
    for (const [key, label] of CMAP_CATEGORIES) {
      const full = _cmapForCategory(key);
      const base = Colormaps.baseName(full);
      const sel = U.el("select", { style: "flex:0 0 38%" });
      for (const name of Colormaps.names()) {
        sel.append(U.el("option", { value: name, selected: name === base ? "" : null }, name));
      }
      const preview = U.el("span", {
        class: "cmap-preview",
        style: `background:${Colormaps.cssGradient(full)}`,
      });
      sel.addEventListener("change", () => {
        const cmap = Colormaps.withReverse(sel.value, Colormaps.isReversed(_cmapForCategory(key)));
        _applyCategoryStyle(key, { cmap });
        preview.style.background = Colormaps.cssGradient(cmap);
      });
      const gear = U.miniBtn("gear", "Z-limits & reverse…",
        () => _categoryStyleDialog(key, label, () => _buildColormapControls(container)));
      container.append(U.el("div", { class: "cmap-row" },
        U.el("span", { class: "cmap-cat" }, label), sel, preview, gear));
    }
  }

  /* Formatting of interpolated sets is remembered (persisted per
   * project); every sample layer picks WHICH interpolated set's
   * formatting to use - defaulting to its own quantity - so e.g. LiDAR
   * samples plot on the same colormap+range as z.grd. */

  function _targetStyles() {
    App.state.ui.targetStyles = App.state.ui.targetStyles || {};
    return App.state.ui.targetStyles;
  }

  function _rememberTargetStyle(target, style) {
    _targetStyles()[target] = { cmap: style.cmap, min: style.min, max: style.max };
    App.touchUi();
  }

  /* "own" or the target name whose formatting this sample layer uses */
  function _styleLinkOf(layer) {
    const links = App.state.ui.styleLink || {};
    if (links[layer.id]) return links[layer.id];
    const q = _quantityOf(layer);
    if (q && layer.id.startsWith("raw-") && _targetStyles()[q]) return q;
    return "own";
  }

  function _setStyleLink(layer, link) {
    App.state.ui.styleLink = App.state.ui.styleLink || {};
    App.state.ui.styleLink[layer.id] = link;
    App.touchUi();
    _syncSharedStyles();
  }

  const pointStyles = new Map();   // raw layer id -> own {cmap, min, max}

  function _ownPointStyle(layer) {
    if (!pointStyles.has(layer.id)) {
      const r = dataRanges.get(layer.id) || [0, 1];
      pointStyles.set(layer.id, {
        cmap: _defaultCmapForName(_quantityOf(layer), "topo_dutch"), min: r[0], max: r[1],
      });
    }
    return pointStyles.get(layer.id);
  }

  /* The style a raw layer should currently display with. */
  function _effectiveRawStyle(layer) {
    const link = _styleLinkOf(layer);
    if (link !== "own") {
      const live = FieldLayer.get(`field-domain-${link}`);
      if (live) return { ...live.style, _source: link };
      const stored = _targetStyles()[link];
      if (stored) return { ...stored, _source: link };
    }
    const fieldLayer = FieldLayer.get(_fieldIdFor(layer));
    if (fieldLayer && !pointLayers.has(layer.id)) return { ...fieldLayer.style, _source: "own" };
    return { ..._ownPointStyle(layer), _source: "own" };
  }

  function _syncSharedStyles() {
    const map = MapView.instance();
    for (const layer of App.state.layers) {
      if (!layer.visible || !layer.id.startsWith("raw-")) continue;
      const style = _effectiveRawStyle(layer);
      const fieldLayer = FieldLayer.get(_fieldIdFor(layer));
      if (fieldLayer && !pointLayers.has(layer.id)) {
        fieldLayer.setStyle({ cmap: style.cmap, min: style.min, max: style.max });
      }
      if (pointLayers.has(layer.id) && map.getLayer(`pts-${layer.id}`)) {
        map.setPaintProperty(`pts-${layer.id}`, "circle-color",
          _pointColorExpr(style.cmap, style.min, style.max));
      }
    }
    _syncColorbars();
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

  /* One colorbar per distinct visible style (interpolated sets + any
   * unlinked sample layers), top-center on the map, hideable. */
  function _syncColorbars() {
    const wrap = document.getElementById("map-wrap");
    if (!wrap) return;
    let box = document.getElementById("colorbars");
    if (!box) {
      box = U.el("div", { id: "colorbars" });
      wrap.append(box);
    }
    U.clear(box);

    const entries = new Map();   // key -> {title, style}
    for (const layer of App.state.layers) {
      if (!layer.visible) continue;
      if (layer.id.startsWith("domain-")) {
        const fieldLayer = FieldLayer.get(_fieldIdFor(layer));
        if (fieldLayer) {
          const target = layer.id.slice("domain-".length);
          entries.set(`target:${target}`, { title: target, style: fieldLayer.style });
        }
      } else if (layer.id.startsWith("raw-")) {
        const style = _effectiveRawStyle(layer);
        if (style._source !== "own") {
          if (!entries.has(`target:${style._source}`)) {
            entries.set(`target:${style._source}`, { title: style._source, style });
          }
        } else {
          entries.set(layer.id, { title: _quantityOf(layer) || layer.title, style });
        }
      }
    }

    if (!entries.size) { box.style.display = "none"; return; }
    box.style.display = "";

    // eye toggle to hide/show the legend (colormap settings live in the
    // Viewer tab's Colormaps group, so there's no gear here)
    const hidden = App.state.ui.colorbars === false;
    const toggle = U.el("span", {
      id: "colorbars-toggle", class: `eye ${hidden ? "off" : ""}`,
      title: hidden ? "Show legend" : "Hide legend",
    }, "👁");
    toggle.addEventListener("click", () => {
      App.state.ui.colorbars = hidden;   // toggled
      App.touchUi();
      _syncColorbars();
    });

    if (hidden) {
      box.append(toggle);
      return;
    }
    for (const { title, style } of entries.values()) {
      box.append(U.el("div", { class: "colorbar" },
        U.el("span", { class: "colorbar-title" }, title),
        U.el("span", { class: "colorbar-min" }, U.fmtNum(style.min, 3)),
        U.el("span", {
          class: "colorbar-ramp",
          style: `background:${Colormaps.cssGradient(style.cmap)}`,
        }),
        U.el("span", { class: "colorbar-max" }, U.fmtNum(style.max, 3))));
    }
    box.append(U.el("div", { class: "colorbars-tools" }, toggle));
  }

  /* Per-layer style editor popup. Sample (raw) layers additionally
   * choose WHOSE formatting to use (an interpolated set, or custom). */
  function _styleEditor(layerInfo) {
    const isRaw = layerInfo.id.startsWith("raw-");
    const isPoints = pointLayers.has(layerInfo.id);
    const fieldLayer = FieldLayer.get(_fieldIdFor(layerInfo));
    if (!fieldLayer && !isPoints) return;
    const target = layerInfo.id.startsWith("domain-")
      ? layerInfo.id.slice("domain-".length) : null;
    const popup = Popup.open({ title: `Style — ${layerInfo.title}`, width: 400 });
    const map = MapView.instance();

    const current = () => isRaw
      ? _effectiveRawStyle(layerInfo)
      : fieldLayer.style;

    const applyCustom = (patch) => {
      if (isRaw && isPoints) {
        const own = _ownPointStyle(layerInfo);
        Object.assign(own, patch);
        if (map.getLayer(`pts-${layerInfo.id}`)) {
          map.setPaintProperty(`pts-${layerInfo.id}`, "circle-color",
            _pointColorExpr(own.cmap, own.min, own.max));
        }
      } else if (fieldLayer) {
        fieldLayer.setStyle(patch);
        if (target) _rememberTargetStyle(target, fieldLayer.style);
      }
      _syncColorbars();
    };

    // ---- formatting source (sample layers only) ----
    const customBox = U.el("div");
    if (isRaw) {
      const linkSel = U.el("select", {});
      linkSel.append(U.el("option", { value: "own" }, "custom (this layer only)"));
      const targets = new Set([
        ...Layers.byGroup("domain").map((l) => l.id.slice("domain-".length)),
        ...Object.keys(_targetStyles()),
      ]);
      for (const t of targets) {
        linkSel.append(U.el("option", { value: t }, `like ${t} (.grd)`));
      }
      linkSel.value = _styleLinkOf(layerInfo);
      const syncEnabled = () => {
        customBox.style.opacity = linkSel.value === "own" ? "" : ".45";
        customBox.style.pointerEvents = linkSel.value === "own" ? "" : "none";
      };
      linkSel.addEventListener("change", () => {
        _setStyleLink(layerInfo, linkSel.value);
        syncEnabled();
      });
      popup.body.append(
        U.el("div", { class: "form-row" }, U.el("label", {}, "Formatting"), linkSel));
      setTimeout(syncEnabled);
    }

    // colours & z-limits are a property of the variable category, edited once
    // in the Colormaps group; here the layer only picks WHICH category it uses
    const catSel = U.el("select", {});
    const curCat = _categoryOfLayer(layerInfo);
    for (const [ck, clabel] of CMAP_CATEGORIES) {
      catSel.append(U.el("option", { value: ck, selected: ck === curCat ? "" : null }, clabel));
    }
    const preview = U.el("div", {
      style: `height:10px;border-radius:5px;margin:4px 0;background:${Colormaps.cssGradient(_cmapForCategory(curCat))}`,
    });
    catSel.addEventListener("change", () => {
      _setLayerCategory(layerInfo, catSel.value);
      preview.style.background = Colormaps.cssGradient(_cmapForCategory(catSel.value));
    });
    const gear = U.el("button", { class: "ghost", style: "font-size:11.5px" },
      "Edit this variable's colormap & limits…");
    gear.addEventListener("click", () => {
      const cat = catSel.value;
      const label = (CMAP_CATEGORIES.find(([k]) => k === cat) || [, cat])[1];
      _categoryStyleDialog(cat, label, () => {
        preview.style.background = Colormaps.cssGradient(_cmapForCategory(cat));
      });
    });

    customBox.append(
      U.el("div", { class: "form-row" }, U.el("label", {}, "Variable"), catSel),
      preview,
      gear,
      U.el("div", { class: "muted", style: "font-size:11.5px;margin-top:4px" },
        "Colormap & z-limits are shared by every source of this variable."),
    );
    if (fieldLayer && !isPoints) {
      const opacity = U.el("input", { type: "range", min: 0, max: 1, step: 0.05, value: fieldLayer.style.opacity });
      opacity.addEventListener("input", () => fieldLayer.setStyle({ opacity: Number(opacity.value) }));
      customBox.append(U.el("div", { class: "form-row" }, U.el("label", {}, "Opacity"), opacity));
    }
    popup.body.append(customBox);
  }

  /* ================= output loading ================= */

  async function _loadOutput(force = false) {
    if (!App.state.project) return;
    try {
      const m = await Api.get("/api/output/meta");
      // a run that crashed before the first write leaves a netCDF with
      // an empty time dimension - treat it as "no output"
      if (!m.exists || !m.times_epoch || !m.times_epoch.length) {
        meta = null; _buildPanel(); return;
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

      if (!variable || !meta.variables.find((v) => v.name === variable)) {
        const preferred = meta.variables.find((v) => v.name === "zb");
        variable = preferred ? "zb" : (meta.variables[0] || {}).name;
      }

      Layers.register({
        id: OUTPUT_LAYER, group: "output",
        title: `Output: ${meta.file}`,
        subtitle: `${meta.times.length} steps`,
      });
      Playbar.setSource("output",
        meta.times_epoch[0], meta.times_epoch[meta.times_epoch.length - 1]);

      _createOutputLayer();
      _buildPanel();
      await _showStep(0);
      await _restoreCustomLayers();
      _applyLayerOrder();
    } catch (err) {
      console.warn("output load failed", err);
    }
  }

  function _createOutputLayer() {
    const layer = FieldLayer.create(OUTPUT_LAYER, mesh, {
      cmap: "topo_dutch", min: -5, max: 15, opacity: 1,
    });
    const reg = Layers.get(OUTPUT_LAYER);
    if (reg && reg.visible === false) layer.setVisible(false);
  }

  /* ================= frames & scrubbing ================= */

  function _frameKey(t) { return `${variable}|${extraIdx}|${t}`; }

  async function _fetchFrame(t) {
    const key = _frameKey(t);
    if (frameCache.has(key)) return frameCache.get(key);
    if (inflight.has(key)) return inflight.get(key);
    const promise = Api.binary(`/api/output/field?var=${variable}&t=${t}&k=${extraIdx}`)
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

  async function _onClock(epoch) {
    if (!meta || !mesh || !Number.isFinite(epoch)) return;
    const layer = FieldLayer.get(OUTPUT_LAYER);
    if (!layer || !layer.visible) return;
    const { k, frac } = _bracket(epoch);
    const k2 = Math.min(k + 1, meta.times.length - 1);

    const cachedA = frameCache.get(_frameKey(k));
    const cachedB = frameCache.get(_frameKey(k2));
    if (cachedA && cachedB) {
      if (!currentBracket || currentBracket.k !== k || currentBracket.var !== variable) {
        layer.setFrames(cachedA, cachedB, frac);
        currentBracket = { k, var: variable };
        if (autoRange) _applyAutoRange(k);
      } else {
        layer.setFrac(frac);
      }
      if (k2 + 1 < meta.times.length) _fetchFrame(k2 + 1).catch(() => {});
      return;
    }
    try {
      const [a, b] = await Promise.all([_fetchFrame(k), _fetchFrame(k2)]);
      layer.setFrames(a, b, frac);
      currentBracket = { k, var: variable };
      if (autoRange) _applyAutoRange(k);
      if (k2 + 1 < meta.times.length) _fetchFrame(k2 + 1).catch(() => {});
    } catch (err) {
      console.warn("frame fetch failed", err.message);
    }
  }

  async function _showStep(t) {
    if (!meta) return;
    Playbar.setTime(meta.times_epoch[Math.min(t, meta.times_epoch.length - 1)]);
  }

  function _applyAutoRange(k) {
    const layer = FieldLayer.get(OUTPUT_LAYER);
    const range = frameRange.get(_frameKey(k));
    if (layer && range && Number.isFinite(range[0])) {
      layer.setStyle({ min: range[0], max: range[1] });
      _syncRangeInputs(range);
    }
  }

  function _syncRangeInputs(range) {
    if (els.min && document.activeElement !== els.min) els.min.value = U.fmtNum(range[0], 4);
    if (els.max && document.activeElement !== els.max) els.max.value = U.fmtNum(range[1], 4);
  }

  /* ================= panel ================= */

  function _buildPanel() {
    let panel = document.getElementById("viewer-panel");
    if (!panel) return;
    U.clear(panel);
    els = {};

    if (!App.state.project) {
      panel.append(U.el("div", { class: "muted" }, "Open a project first."));
      return;
    }

    // --- colormaps (per variable category) — own section at the top ---
    const cmapSection = U.section("Colormaps", { collapsed: true });
    cmapSection.wrap.id = "viewer-cmaps-section";
    _buildColormapControls(cmapSection.body);
    panel.append(cmapSection.wrap);

    // --- layer groups (each its own collapsible section) ---
    els.groups = U.el("div");
    panel.append(els.groups);
    _renderGroups();

    // --- output controls ---
    if (!meta) {
      panel.append(U.el("div", { class: "muted", style: "margin-top:10px" },
        "No model output yet — run a simulation in the Run tab."));
      return;
    }

    const outputSection = U.section("Output display");
    panel.append(outputSection.wrap);
    // the remaining controls all land inside the section body
    panel = outputSection.body;

    const varSel = U.el("select", {});
    for (const v of meta.variables) {
      varSel.append(U.el("option", {
        value: v.name, selected: v.name === variable ? "" : null,
        title: v.long_name,
      }, `${v.name}${v.units ? ` [${v.units}]` : ""}`));
    }
    varSel.addEventListener("change", async () => {
      variable = varSel.value;
      currentBracket = null;
      _renderExtraDims();
      await _onClock(App.state.clock.t);
    });
    panel.append(U.el("div", { class: "form-row" }, U.el("label", {}, "Variable"), varSel));

    els.extraBox = U.el("div");
    panel.append(els.extraBox);
    _renderExtraDims();

    const layer = FieldLayer.get(OUTPUT_LAYER);
    // pick the VARIABLE category (Elevation, Wind speed, …) — the colormap
    // itself and its z-limits are set once per category in the Colormaps group
    const outCat = _categoryOfKey(`output:${variable}`, variable);
    const cmapSel = U.el("select", {});
    for (const [ck, clabel] of CMAP_CATEGORIES) {
      cmapSel.append(U.el("option", { value: ck, selected: ck === outCat ? "" : null }, clabel));
    }
    const cmapPreview = U.el("div", {
      style: `height:10px;border-radius:5px;margin:4px 0;background:${Colormaps.cssGradient(_cmapForCategory(outCat))}`,
    });
    cmapSel.addEventListener("change", () => {
      _setLayerCategory({ id: OUTPUT_LAYER }, cmapSel.value);
      cmapPreview.style.background = Colormaps.cssGradient(_cmapForCategory(cmapSel.value));
    });
    panel.append(U.el("div", { class: "form-row" }, U.el("label", {}, "Colormap"), cmapSel),
      cmapPreview,
      U.el("div", { class: "muted", style: "font-size:11.5px;margin:-2px 0 4px" },
        "Edit the colormap & default limits in the Colormaps group above."));

    const autoCb = U.el("input", { type: "checkbox" });
    autoCb.checked = autoRange;
    autoCb.addEventListener("change", () => {
      autoRange = autoCb.checked;
      if (autoRange && currentBracket) _applyAutoRange(currentBracket.k);
    });
    els.min = U.el("input", { type: "text", style: "width:70px" });
    els.max = U.el("input", { type: "text", style: "width:70px" });
    if (layer) _syncRangeInputs([layer.style.min, layer.style.max]);
    const commitRange = () => {
      const lo = Number(els.min.value), hi = Number(els.max.value);
      if (Number.isFinite(lo) && Number.isFinite(hi) && hi > lo) {
        autoRange = false;
        autoCb.checked = false;
        const l = FieldLayer.get(OUTPUT_LAYER);
        if (l) l.setStyle({ min: lo, max: hi });
      }
    };
    for (const input of [els.min, els.max]) {
      input.addEventListener("blur", commitRange);
      input.addEventListener("keydown", (ev) => { if (ev.key === "Enter") input.blur(); });
    }
    panel.append(
      U.el("div", { class: "form-row" }, U.el("label", {}, "Auto range"), autoCb),
      U.el("div", { class: "form-row" }, U.el("label", {}, "Min / max"), els.min, els.max),
    );

    const opacity = U.el("input", { type: "range", min: 0, max: 1, step: 0.05,
      value: layer ? layer.style.opacity : 1 });
    opacity.addEventListener("input", () => {
      const l = FieldLayer.get(OUTPUT_LAYER);
      if (l) l.setStyle({ opacity: Number(opacity.value) });
    });
    panel.append(U.el("div", { class: "form-row" }, U.el("label", {}, "Opacity"), opacity));

    const probeBtn = U.el("button", { class: "ghost" }, "Probe cell (click map)");
    probeBtn.addEventListener("click", () => _armProbe(probeBtn));
    const zoomBtn = U.el("button", { class: "ghost" }, "Zoom to output");
    zoomBtn.addEventListener("click", () => {
      if (!mesh) return;
      let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
      for (let i = 0; i < mesh.x.length; i += 1) {
        if (mesh.x[i] < minX) minX = mesh.x[i];
        if (mesh.x[i] > maxX) maxX = mesh.x[i];
        if (mesh.y[i] < minY) minY = mesh.y[i];
        if (mesh.y[i] > maxY) maxY = mesh.y[i];
      }
      MapView.fitModelBounds(minX, minY, maxX, maxY);
    });
    panel.append(U.el("div", { class: "btn-row" }, probeBtn, zoomBtn));

    panel.append(U.el("div", { class: "muted", style: "font-size:11.5px" },
      `${meta.times.length} output steps — scrub or play with the time bar below.`));
  }

  function _renderExtraDims() {
    const box = els.extraBox;
    if (!box) return;
    U.clear(box);
    const info = meta.variables.find((v) => v.name === variable);
    if (!info || !info.extra_dims.length) { extraIdx = "0"; return; }
    const selects = [];
    for (const dim of info.extra_dims) {
      const select = U.el("select", {});
      for (let i = 0; i < dim.size; i += 1) {
        select.append(U.el("option", { value: i }, `${dim.name} ${i}`));
      }
      select.addEventListener("change", async () => {
        extraIdx = selects.map((sel) => sel.value).join(",");
        currentBracket = null;
        await _onClock(App.state.clock.t);
      });
      selects.push(select);
      box.append(U.el("div", { class: "form-row" }, U.el("label", {}, dim.name), select));
    }
    extraIdx = selects.map((sel) => sel.value).join(",");
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
    if (fid === OUTPUT_LAYER) {
      const info = (meta && meta.variables.find((v) => v.name === variable)) || {};
      return { label: variable || "z", unit: info.units || "" };
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
        const res = await Api.get(
          `/api/output/series?var=${variable}&j=${j}&i=${i}&k=${extraIdx}`);
        const info = meta.variables.find((v) => v.name === variable) || {};
        Graphs.registerSource(`probe-${variable}-${j}-${i}`, {
          group: "Output",
          label: `${variable} @ cell (${j},${i})`,
          unit: info.units || "",
          data: [res.t_epoch, res.values],
          range: res.t_epoch.length
            ? [res.t_epoch[0], res.t_epoch[res.t_epoch.length - 1]] : null,
        }, { select: true });
        MapView.setLabel(`probe-${j}-${i}`,
          [mesh.x[best], mesh.y[best]], `(${j},${i})`, "boundary-lateral");
      } catch (err) {
        U.toast(err.message, "error");
      }
    });
  }

  /* ================= domain / raw layer toggling ================= */

  async function _onLayerVisibility(layerInfo) {
    if (layerInfo.id === OUTPUT_LAYER) {
      const layer = FieldLayer.get(OUTPUT_LAYER);
      if (layer) layer.setVisible(layerInfo.visible);
      return;
    }
    if (layerInfo.id.startsWith("raw-") && layerInfo.entry) {
      await _toggleRawLayer(layerInfo);
      _applyLayerOrder();
      _syncSharedStyles();
    }
    if (layerInfo.id.startsWith("domain-")) {
      await _toggleDomainLayer(layerInfo);
      _applyLayerOrder();
      _syncSharedStyles();
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
        const style = _effectiveRawStyle(layerInfo);
        MapView.ensureLayer({
          id: `pts-${layerInfo.id}`, type: "circle", source: `pts-${layerInfo.id}`,
          paint: {
            "circle-radius": 2.4,
            "circle-color": _pointColorExpr(style.cmap, style.min, style.max),
          },
        });
      } else {
        const { buffer, headers } = await Api.binary(`/api/domain/rawfield?id=${entry.id}`);
        const parsed = FieldLayer.parseGridfield(buffer, headers);
        dataRanges.set(layerInfo.id, parsed.range);
        const layer = FieldLayer.create(id, parsed.mesh, {
          cmap: _defaultCmapForName(_quantityOf(layerInfo), "topo_dutch"),
          min: parsed.range[0], max: parsed.range[1], opacity: 0.85,
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
      _syncColorbars();
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
      // reuse the remembered formatting of this interpolated set (also
      // the source style for linked sample layers)
      const stored = _targetStyles()[target];
      const layer = FieldLayer.create(id, parsed.mesh, stored
        ? { cmap: stored.cmap, min: stored.min, max: stored.max, opacity: 0.9 }
        : {
          cmap: _defaultCmapForName(target, "viridis"),
          min: parsed.range[0], max: parsed.range[1], opacity: 0.9,
        });
      layer.setFrames(parsed.values);
      if (!stored) _rememberTargetStyle(target, layer.style);
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
          const key = `tgt:${t}`, label = `interpolated: ${t}`;
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
