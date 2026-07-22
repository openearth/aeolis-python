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
    ["domain", "Interpolated (.grd)"],
    ["rawdata", "Sample data"],
    ["grid", "Grid"],
    ["objects", "Objects"],
    ["background", "Background"],
  ];

  const loadingLayers = new Set();
  const dataRanges = new Map();   // layer id -> [lo, hi] of the loaded data
  const pointLayers = new Set();  // layer ids rendered as circle layers

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
    _buildPanel();
  }

  function _reset() {
    meta = null; mesh = null; variable = null;
    frameCache.clear(); frameRange.clear(); inflight.clear();
    dataRanges.clear(); pointLayers.clear();
    currentBracket = null;
    FieldLayer.remove(OUTPUT_LAYER);
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

    for (const [group, title] of _orderedGroups()) {
      let contentBuilder = null;
      let count = 0;
      let groupLayers = [];
      if (group === "objects") {
        if (!App.state.objects.length) continue;
        count = App.state.objects.length;
        contentBuilder = (body) => _renderObjectCards(body);
      } else if (group === "background") {
        if (CRS.isLocal()) continue;
        contentBuilder = (body) => _renderBackgroundCards(body);
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
      // grip to drag whole groups + group show/hide
      section.head.prepend(U.el("span", {
        class: "drag-grip", title: "Drag to reorder groups (top = drawn on top)",
      }, "⠿"));
      section.head.draggable = true;
      if (group !== "background") {
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
   * container element - _renderGroups re-runs often) */
  function _wireGroupDrag(box) {
    if (box._groupDragWired) return;
    box._groupDragWired = true;
    let fromGroup = null;
    box.addEventListener("dragstart", (ev) => {
      const head = ev.target.closest && ev.target.closest(".section > header");
      const section = ev.target.closest && ev.target.closest(".section[data-group]");
      if (!head || !section) return;
      fromGroup = section.dataset.group;
      section.classList.add("dragging");
      ev.dataTransfer.effectAllowed = "move";
      ev.dataTransfer.setData("text/plain", "");
    });
    box.addEventListener("dragend", () => {
      fromGroup = null;
      for (const s of box.querySelectorAll(".section")) {
        s.classList.remove("dragging", "drop-above", "drop-below");
      }
    });
    box.addEventListener("dragover", (ev) => {
      if (!fromGroup) return;
      ev.preventDefault();
      const section = ev.target.closest && ev.target.closest(".section[data-group]");
      for (const s of box.querySelectorAll(".section")) {
        s.classList.remove("drop-above", "drop-below");
      }
      if (!section || section.dataset.group === fromGroup) return;
      const rect = section.getBoundingClientRect();
      section.classList.add(ev.clientY > rect.top + rect.height / 2 ? "drop-below" : "drop-above");
    });
    box.addEventListener("drop", (ev) => {
      if (!fromGroup) return;
      ev.preventDefault();
      const section = ev.target.closest && ev.target.closest(".section[data-group]");
      if (!section || section.dataset.group === fromGroup) return;
      const rect = section.getBoundingClientRect();
      const below = ev.clientY > rect.top + rect.height / 2;
      const order = _orderedGroups().map(([g]) => g);
      const from = order.indexOf(fromGroup);
      order.splice(from, 1);
      let to = order.indexOf(section.dataset.group) + (below ? 1 : 0);
      order.splice(to, 0, fromGroup);
      App.state.ui.viewerGroupOrder = order;
      App.touchUi();
      fromGroup = null;
      _renderGroups();
      _applyLayerOrder();
    });
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
      class: "obj-card", draggable: count > 1 ? "true" : null, dataset: { idx },
    },
      count > 1 ? U.el("span", { class: "drag-grip", title: "Drag to reorder (top = drawn on top)" }, "⠿") : null,
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

    // styling for field layers and point (sample) layers
    const fieldId = _fieldIdFor(layer);
    if ((fieldId && FieldLayer.get(fieldId)) || pointLayers.has(layer.id)) {
      // simple per-layer colormap dropdown (overrides the category default)
      const cur = _currentCmapOf(layer);
      const sel = U.el("select", { class: "cmap-mini", title: "Colormap for this layer" });
      for (const name of Colormaps.names()) {
        sel.append(U.el("option", { value: name, selected: name === cur ? "" : null }, name));
      }
      sel.addEventListener("click", (ev) => ev.stopPropagation());
      sel.addEventListener("change", () => _setLayerCmap(layer, sel.value));
      card.append(sel);
      card.append(U.miniBtn("gear", "More style options…", () => _styleEditor(layer)));
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

  /* drag & drop within one list (shared pattern with the Domain tab) */
  function _wireCardDrag(list, onDrop) {
    let fromIdx = null;
    list.addEventListener("dragstart", (ev) => {
      const card = ev.target.closest(".obj-card");
      if (!card) return;
      fromIdx = Number(card.dataset.idx);
      card.classList.add("dragging");
      ev.dataTransfer.effectAllowed = "move";
      ev.dataTransfer.setData("text/plain", "");   // Firefox needs data to drag
    });
    list.addEventListener("dragend", () => {
      fromIdx = null;
      for (const c of list.querySelectorAll(".obj-card")) {
        c.classList.remove("dragging", "drop-above", "drop-below");
      }
    });
    list.addEventListener("dragover", (ev) => {
      if (fromIdx === null) return;
      ev.preventDefault();
      const card = ev.target.closest(".obj-card");
      for (const c of list.querySelectorAll(".obj-card")) {
        c.classList.remove("drop-above", "drop-below");
      }
      if (!card) return;
      const rect = card.getBoundingClientRect();
      const below = ev.clientY > rect.top + rect.height / 2;
      card.classList.add(below ? "drop-below" : "drop-above");
    });
    list.addEventListener("drop", (ev) => {
      if (fromIdx === null) return;
      ev.preventDefault();
      const card = ev.target.closest(".obj-card");
      if (!card) return;
      const rect = card.getBoundingClientRect();
      const below = ev.clientY > rect.top + rect.height / 2;
      let toIdx = Number(card.dataset.idx) + (below ? 1 : 0);
      if (toIdx > fromIdx) toIdx -= 1;
      if (toIdx !== fromIdx) onDrop(fromIdx, toIdx);
      fromIdx = null;
    });
  }

  function _fieldIdFor(layer) {
    if (layer.id === OUTPUT_LAYER) return OUTPUT_LAYER;
    if (layer.id.startsWith("domain-") || layer.id.startsWith("raw-")) return `field-${layer.id}`;
    return null;
  }

  function _renderObjectCards(body) {
    const list = U.el("div", { class: "obj-list" });
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

  /* Apply the tree order to the map. Groups draw in their (draggable)
   * tree order - the first group/row of the tree ends up on top - so
   * e.g. the grid can be placed in front of or behind the objects. */
  const GROUP_MAP_IDS = {
    grid: ["grid-fill", "grid-lines", "grid-outline", "grid-shear", "grid-shear-inner"],
    objects: ["objects-fill", "objects-fill-outline", "objects-lines"],
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
    ["elevation", "Elevation (zb, zne, zsep)", "topo_dutch"],
    ["bed_change", "Bed level change", "RdBu"],
    ["veg_density", "Vegetation density", "viridis"],
    ["veg_height", "Vegetation height", "viridis"],
    ["mask", "Masks", "gray"],
    ["wind_speed", "Wind speed", "turbo"],
    ["shear_stress", "Shear stress", "plasma"],
    ["shear_velocity", "Shear velocity", "plasma"],
    ["sed_conc", "Sediment concentration", "turbo"],
    ["sed_transport", "Sediment transport", "turbo"],
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

  /* Apply a category's colormap to every currently loaded layer of that
   * category (interpolated sets, sample layers and the output layer). */
  function _applyCategoryCmap(cat, cmap) {
    _varCmaps()[cat] = cmap;
    App.touchUi();
    const map = MapView.instance();
    // interpolated sets + remembered target styles
    for (const [target, style] of Object.entries(_targetStyles())) {
      if (_categoryOfName(target) !== cat) continue;
      style.cmap = cmap;
      const live = FieldLayer.get(`field-domain-${target}`);
      if (live) live.setStyle({ cmap });
    }
    // output layer, if its variable belongs to this category
    if (variable && _categoryOfName(variable) === cat) {
      const l = FieldLayer.get(OUTPUT_LAYER);
      if (l) l.setStyle({ cmap });
    }
    // sample layers styled on their own
    for (const layer of App.state.layers) {
      if (!layer.id.startsWith("raw-")) continue;
      if (_categoryOfName(_quantityOf(layer)) !== cat) continue;
      const own = _ownPointStyle(layer);
      own.cmap = cmap;
      if (pointLayers.has(layer.id) && map.getLayer(`pts-${layer.id}`)) {
        map.setPaintProperty(`pts-${layer.id}`, "circle-color",
          _pointColorExpr(cmap, own.min, own.max));
      }
    }
    _syncSharedStyles();   // also refreshes colorbars
  }

  /* Build the category → colormap controls into a container (used both in
   * the Viewer panel and in the popup opened from the map colorbars). */
  function _buildColormapControls(container) {
    U.clear(container);
    container.append(U.el("div", { class: "muted", style: "font-size:11.5px;margin:0 0 6px" },
      "Set the colormap per variable; every layer of that variable uses it "
      + "(override an individual layer with its own dropdown)."));
    for (const [key, label] of CMAP_CATEGORIES) {
      const current = _cmapForCategory(key);
      const sel = U.el("select", { style: "flex:0 0 42%" });
      for (const name of Colormaps.names()) {
        sel.append(U.el("option", { value: name, selected: name === current ? "" : null }, name));
      }
      const preview = U.el("span", {
        class: "cmap-preview",
        style: `background:${Colormaps.cssGradient(current)}`,
      });
      sel.addEventListener("change", () => {
        _applyCategoryCmap(key, sel.value);
        preview.style.background = Colormaps.cssGradient(sel.value);
      });
      container.append(U.el("div", { class: "cmap-row" },
        U.el("span", { class: "cmap-cat" }, label), sel, preview));
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

    const hidden = App.state.ui.colorbars === false;
    const toggle = U.el("button", {
      id: "colorbars-toggle",
      title: hidden ? "Show colorbars" : "Hide colorbars",
    }, hidden ? "▤ legend" : "✕");
    toggle.addEventListener("click", () => {
      App.state.ui.colorbars = hidden;   // toggled
      App.touchUi();
      _syncColorbars();
    });

    // settings gear: opens the per-variable colormap controls (works from
    // any tab since the colorbars overlay the map everywhere)
    const gear = U.el("button", {
      id: "colorbars-settings", title: "Colormap settings",
    }, U.icon("palette", 13));
    gear.addEventListener("click", () => {
      const popup = Popup.open({ title: "Colormaps", width: 440 });
      _buildColormapControls(popup.body);
    });

    if (hidden) {
      box.append(gear, toggle);
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
    box.append(U.el("div", { class: "colorbars-tools" }, gear, toggle));
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

    const style = current();
    const cmapSel = U.el("select", {});
    for (const name of Colormaps.names()) {
      cmapSel.append(U.el("option", { value: name, selected: name === style.cmap ? "" : null }, name));
    }
    const preview = U.el("div", {
      style: `height:10px;border-radius:5px;margin:4px 0;background:${Colormaps.cssGradient(style.cmap)}`,
    });
    cmapSel.addEventListener("change", () => {
      applyCustom({ cmap: cmapSel.value });
      preview.style.background = Colormaps.cssGradient(cmapSel.value);
    });

    const minIn = U.el("input", { type: "text", value: U.fmtNum(style.min, 4), style: "width:80px" });
    const maxIn = U.el("input", { type: "text", value: U.fmtNum(style.max, 4), style: "width:80px" });
    const commitRange = () => {
      const lo = Number(minIn.value), hi = Number(maxIn.value);
      if (Number.isFinite(lo) && Number.isFinite(hi) && hi > lo) applyCustom({ min: lo, max: hi });
    };
    for (const input of [minIn, maxIn]) {
      input.addEventListener("blur", commitRange);
      input.addEventListener("keydown", (ev) => { if (ev.key === "Enter") input.blur(); });
    }

    customBox.append(
      U.el("div", { class: "form-row" }, U.el("label", {}, "Colormap"), cmapSel),
      preview,
      U.el("div", { class: "form-row" }, U.el("label", {}, "Min / max"), minIn, maxIn),
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

  return { init };
})();
