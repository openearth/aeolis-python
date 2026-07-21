/* Viewer tab: the single place to manage everything shown on the map.
 *
 * - Layer tree (Output / Interpolated .grd / Sample data / Grid /
 *   Objects / Background) with visibility, ordering and per-layer
 *   styling.
 * - Objects management: rename, recolor, zoom, delete.
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
    ["domain", "Interpolated (.grd)"],
    ["rawdata", "Sample data"],
    ["grid", "Grid"],
    ["objects", "Objects"],
    ["background", "Background"],
  ];

  function init() {
    Tabs.register("viewer", { enter: _enter });
    App.on("project", () => { _reset(); _loadOutput(); });
    App.on("run-finished", () => _loadOutput(true));
    App.on("clock-tick", _onClock);
    App.on("layer-visibility", _onLayerVisibility);
    App.on("layers", () => _renderTree());
    App.on("objects", () => _renderTree());
    App.on("basemap", () => _renderTree());
    App.on("layer-order", _applyLayerOrder);
    _buildPanel();
  }

  function _reset() {
    meta = null; mesh = null; variable = null;
    frameCache.clear(); frameRange.clear(); inflight.clear();
    currentBracket = null;
    FieldLayer.remove(OUTPUT_LAYER);
    Playbar.removeSource("output");
  }

  function _enter() {
    if (!meta) _loadOutput();
    _buildPanel();
  }

  /* ================= layer tree ================= */

  function _renderTree() {
    const tree = els.tree;
    if (!tree) return;
    U.clear(tree);

    for (const [group, title] of GROUPS) {
      if (group === "objects") { _renderObjectsGroup(tree, title); continue; }
      if (group === "background") { _renderBackgroundGroup(tree, title); continue; }
      const layers = Layers.byGroup(group);
      if (!layers.length) continue;
      tree.append(U.el("div", { class: "lp-group-head" }, title));
      layers.forEach((layer, idx) => {
        tree.append(_layerRow(layer, idx, layers.length));
      });
    }
    if (!tree.children.length) {
      tree.append(U.el("div", { class: "muted" }, "No layers yet."));
    }
  }

  const loadingLayers = new Set();

  function _layerRow(layer, idx, count) {
    const loading = loadingLayers.has(layer.id);
    const eye = loading
      ? U.el("span", { class: "spin", title: "Loading…" })
      : U.el("span", { class: `eye ${layer.visible ? "" : "off"}`, title: "Show/hide" }, "👁");
    if (!loading) {
      eye.addEventListener("click", () => {
        layer.visible = !layer.visible;
        App.emit("layer-visibility", layer);
        _renderTree();
      });
    }

    const row = U.el("div", { class: "lp-row" }, eye,
      U.el("span", { class: "lp-name", title: layer.title }, layer.title),
      layer.subtitle ? U.el("span", { class: "lp-mini" }, layer.subtitle) : null);

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
        _renderTree();
      });
      row.append(spBtn);
    }

    // ordering within group (top row = drawn on top)
    if (count > 1) {
      const up = U.miniBtn("up", "Raise (drawn on top)", () => Layers.move(layer.id, -1));
      const down = U.miniBtn("down", "Lower", () => Layers.move(layer.id, +1));
      up.disabled = idx === 0;
      down.disabled = idx === count - 1;
      row.append(up, down);
    }

    // styling for field layers
    const fieldId = _fieldIdFor(layer);
    if (fieldId && FieldLayer.get(fieldId)) {
      row.append(U.miniBtn("gear", "Style…", () => _styleEditor(fieldId, layer.title)));
    }
    return row;
  }

  function _fieldIdFor(layer) {
    if (layer.id === OUTPUT_LAYER) return OUTPUT_LAYER;
    if (layer.id.startsWith("domain-") || layer.id.startsWith("raw-")) return `field-${layer.id}`;
    return null;
  }

  function _renderObjectsGroup(tree, title) {
    if (!App.state.objects.length) return;
    tree.append(U.el("div", { class: "lp-group-head" }, title));
    for (const obj of App.state.objects) {
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

      tree.append(U.el("div", { class: "lp-row" }, eye, swatch, name,
        U.el("span", { class: "lp-mini" }, obj.kind), zoom, del));
    }
  }

  function _renderBackgroundGroup(tree, title) {
    if (CRS.isLocal()) return;
    tree.append(U.el("div", { class: "lp-group-head" }, title));
    for (const [key, label] of [["gray", "Grey map"], ["sat", "Satellite"], ["none", "None"]]) {
      const active = App.state.ui.basemap === key;
      const row = U.el("div", { class: `lp-row ${active ? "selected" : ""}` },
        U.el("span", { class: "eye" }, active ? "●" : "○"),
        U.el("span", { class: "lp-name" }, label));
      row.addEventListener("click", () => MapView.setBasemap(key));
      tree.append(row);
    }
  }

  /* Apply the tree order to the map: iterate groups bottom-up so the
   * first row of the tree ends up on top; grid and objects stay above
   * the field layers. */
  function _applyLayerOrder() {
    const map = MapView.instance();
    if (!map || !map.isStyleLoaded()) return;
    const ordered = [];
    for (const group of ["rawdata", "domain", "output"]) {
      const layers = Layers.byGroup(group);
      for (let i = layers.length - 1; i >= 0; i -= 1) {
        const fieldId = _fieldIdFor(layers[i]);
        if (fieldId && map.getLayer(fieldId)) ordered.push(fieldId);
        const ptsId = `pts-${layers[i].id}`;
        if (map.getLayer(ptsId)) ordered.push(ptsId);
      }
    }
    // fields bottom-to-top, then grid, then objects on top
    for (const id of [...ordered,
      "grid-fill", "grid-lines", "grid-outline", "grid-shear", "grid-shear-inner",
      "objects-fill", "objects-fill-outline", "objects-lines"]) {
      if (map.getLayer(id)) map.moveLayer(id);
    }
  }

  /* Per-layer style editor popup (any field layer). */
  function _styleEditor(fieldId, title) {
    const layer = FieldLayer.get(fieldId);
    if (!layer) return;
    const popup = Popup.open({ title: `Style — ${title}`, width: 380 });

    const cmapSel = U.el("select", {});
    for (const name of Colormaps.names()) {
      cmapSel.append(U.el("option", { value: name, selected: name === layer.style.cmap ? "" : null }, name));
    }
    const preview = U.el("div", {
      style: `height:10px;border-radius:5px;margin:4px 0;background:${Colormaps.cssGradient(layer.style.cmap)}`,
    });
    cmapSel.addEventListener("change", () => {
      layer.setStyle({ cmap: cmapSel.value });
      preview.style.background = Colormaps.cssGradient(cmapSel.value);
    });

    const minIn = U.el("input", { type: "text", value: U.fmtNum(layer.style.min, 4), style: "width:80px" });
    const maxIn = U.el("input", { type: "text", value: U.fmtNum(layer.style.max, 4), style: "width:80px" });
    const commitRange = () => {
      const lo = Number(minIn.value), hi = Number(maxIn.value);
      if (Number.isFinite(lo) && Number.isFinite(hi) && hi > lo) layer.setStyle({ min: lo, max: hi });
    };
    for (const input of [minIn, maxIn]) {
      input.addEventListener("blur", commitRange);
      input.addEventListener("keydown", (ev) => { if (ev.key === "Enter") input.blur(); });
    }

    const opacity = U.el("input", { type: "range", min: 0, max: 1, step: 0.05, value: layer.style.opacity });
    opacity.addEventListener("input", () => layer.setStyle({ opacity: Number(opacity.value) }));

    popup.body.append(
      U.el("div", { class: "form-row" }, U.el("label", {}, "Colormap"), cmapSel),
      preview,
      U.el("div", { class: "form-row" }, U.el("label", {}, "Min / max"), minIn, maxIn),
      U.el("div", { class: "form-row" }, U.el("label", {}, "Opacity"), opacity),
    );
  }

  /* ================= output loading ================= */

  async function _loadOutput(force = false) {
    if (!App.state.project) return;
    try {
      const m = await Api.get("/api/output/meta");
      if (!m.exists) { meta = null; _buildPanel(); return; }
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

    // --- layer tree ---
    const layersSection = U.section("Layers");
    els.tree = U.el("div", { class: "layer-tree" });
    layersSection.body.append(els.tree);
    panel.append(layersSection.wrap);
    _renderTree();

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
    const cmapSel = U.el("select", {});
    for (const name of Colormaps.names()) {
      cmapSel.append(U.el("option", { value: name }, name));
    }
    if (layer) cmapSel.value = layer.style.cmap;
    const cmapPreview = U.el("div", {
      style: `height:10px;border-radius:5px;margin:4px 0;background:${Colormaps.cssGradient(cmapSel.value)}`,
    });
    cmapSel.addEventListener("change", () => {
      const l = FieldLayer.get(OUTPUT_LAYER);
      if (l) l.setStyle({ cmap: cmapSel.value });
      cmapPreview.style.background = Colormaps.cssGradient(cmapSel.value);
    });
    panel.append(U.el("div", { class: "form-row" }, U.el("label", {}, "Colormap"), cmapSel));
    panel.append(cmapPreview);

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
        Graphs.add(`probe-${variable}-${j}-${i}`, {
          title: `${variable} @ cell (${j},${i})`,
          data: [res.t_epoch, res.values],
          series: [{}, { label: variable, stroke: "#b4423b", width: 1.5 }],
          timeBased: true,
          axes: [
            { values: (u, ticks) => ticks.map((t) => U.fmtDate(t).slice(5, 16)) },
            { size: 55 },
          ],
          scales: { x: { time: false } },
        });
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
    }
    if (layerInfo.id.startsWith("domain-")) {
      await _toggleDomainLayer(layerInfo);
      _applyLayerOrder();
    }
  }

  async function _toggleRawLayer(layerInfo) {
    const id = `field-${layerInfo.id}`;
    if (!layerInfo.visible) {
      FieldLayer.remove(id);
      MapView.removeLayerAndSource(`pts-${layerInfo.id}`);
      return;
    }
    loadingLayers.add(layerInfo.id);
    _renderTree();
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
        MapView.upsertGeojson(`pts-${layerInfo.id}`, { type: "FeatureCollection", features });
        MapView.ensureLayer({
          id: `pts-${layerInfo.id}`, type: "circle", source: `pts-${layerInfo.id}`,
          paint: {
            "circle-radius": 2.4,
            "circle-color": ["interpolate", ["linear"], ["get", "z"],
              zmin, "#2c3e70", (zmin + zmax) / 2, "#e8d8a0", zmax, "#7a5230"],
          },
        });
      } else {
        const { buffer, headers } = await Api.binary(`/api/domain/rawfield?id=${entry.id}`);
        const parsed = FieldLayer.parseGridfield(buffer, headers);
        const layer = FieldLayer.create(id, parsed.mesh, {
          cmap: "topo_dutch", min: parsed.range[0], max: parsed.range[1], opacity: 0.85,
        });
        layer.setFrames(parsed.values);
      }
    } catch (err) {
      U.toast(`Layer failed: ${err.message}`, "error");
      layerInfo.visible = false;
    } finally {
      loadingLayers.delete(layerInfo.id);
      _renderTree();
    }
  }

  async function _toggleDomainLayer(layerInfo) {
    const id = `field-${layerInfo.id}`;
    const target = layerInfo.id.replace("domain-", "");
    if (!layerInfo.visible) {
      FieldLayer.remove(id);
      return;
    }
    loadingLayers.add(layerInfo.id);
    _renderTree();
    try {
      const k = layerInfo.speciesIdx || 0;
      const { buffer, headers } = await Api.binary(`/api/domain/gridfield?target=${target}&k=${k}`);
      const nSpecies = Number(headers.get("X-Species") || 1);
      Layers.register({ id: layerInfo.id, group: layerInfo.group, title: layerInfo.title,
        species: nSpecies, speciesIdx: k });
      const parsed = FieldLayer.parseGridfield(buffer, headers);
      const layer = FieldLayer.create(id, parsed.mesh, {
        cmap: target === "bed" || target === "ne" ? "topo_dutch" : "viridis",
        min: parsed.range[0], max: parsed.range[1], opacity: 0.9,
      });
      layer.setFrames(parsed.values);
    } catch (err) {
      U.toast(`Layer failed: ${err.message}`, "error");
      // roll the eye back so the tree reflects reality
      layerInfo.visible = false;
    } finally {
      loadingLayers.delete(layerInfo.id);
      _renderTree();
    }
  }

  return { init };
})();
