/* Transect viewer: a second tab in the shared time-series area.
 *
 * Plots datasets ALONG a drawn transect (distance on x) at the current
 * playbar time, reusing the shared TSPlot component so pan / zoom / hover
 * behave exactly like the time graph. Datasets are chosen from a catalog
 * of ALL output variables + interpolated + raw datasets (whether or not
 * they are shown on the map) and sampled on demand.
 *
 * Coastal look: bed level filled sand, water level filled sea-blue down
 * to the bed, separation bubble as a black dotted line.
 */
"use strict";

const Transect = (() => {

  let mode = "graphs";          // "graphs" | "transect"
  let host = null;              // #transect-panel
  let chartEl = null;
  let plot = null;              // TSPlot api
  let curSig = "";              // signature of the currently-built series
  let catalog = [];             // [{key,label,kind}]
  let redrawTimer = null;
  const tabBtns = {};
  const ctrl = {};

  const BED = ["zb"];
  const WATER = ["zs", "zw", "zsw", "water_level"];
  const SEP = ["zsep", "sep", "zsepbub", "z_separation"];
  const isBed = (v) => BED.includes(v);
  const isWater = (v) => WATER.includes(v);
  const isSep = (v) => SEP.includes(v);
  const _varOf = (key) => (key && key.startsWith("out:") ? key.slice(4) : null);
  const PALETTE = ["#2f7fe6", "#27a355", "#a034c6", "#e0a020", "#12a5b5", "#d1387f"];
  const SAND = "#d8c18a";
  const SEA = "rgba(70,130,200,.35)";

  function init() {
    const wrap = document.getElementById("graphs-wrap");
    if (!wrap) return;

    wrap.prepend(U.el("div", { class: "gt-tabs" },
      _tab("graphs", "Time series"), _tab("transect", "Transect")));

    ctrl.transectSel = U.el("select", { class: "gt-sel", title: "Transect to plot" });
    ctrl.transectSel.addEventListener("change", () => _scheduleRedraw(true));
    const drawBtn = U.el("button", { class: "gbtn" }, "Draw transect");
    drawBtn.addEventListener("click", _drawTransect);
    const varsBtn = U.el("button", { class: "gbtn" }, "Data…");
    varsBtn.addEventListener("click", _pickVars);
    ctrl.bar = U.el("div", { class: "gt-bar" },
      U.el("label", {}, "Transect"), ctrl.transectSel, drawBtn, varsBtn);
    chartEl = U.el("div", { class: "gt-chart" });
    host = U.el("div", { id: "transect-panel" }, ctrl.bar, chartEl);
    wrap.insertBefore(host, document.getElementById("graphs"));

    App.on("clock-tick", () => { if (mode === "transect") _scheduleRedraw(false); });
    App.on("objects", () => { if (mode === "transect") { _fillTransectSel(); _scheduleRedraw(false); } });
    App.on("run-finished", () => { catalog = []; if (mode === "transect") _refreshAndDraw(); });
    App.on("project", () => { catalog = []; _destroyPlot(); });

    if (window.ResizeObserver) {
      new ResizeObserver(() => {
        if (mode === "transect" && plot) plot.setSize({ width: chartEl.clientWidth, height: chartEl.clientHeight });
      }).observe(chartEl);
    }
    _setMode("graphs");
  }

  function _tab(key, label) {
    const btn = U.el("button", { class: "gt-tab" }, label);
    btn.addEventListener("click", () => _setMode(key));
    tabBtns[key] = btn;
    return btn;
  }

  function _setMode(m) {
    mode = m;
    for (const [k, b] of Object.entries(tabBtns)) b.classList.toggle("active", k === m);
    const graphs = document.getElementById("graphs");
    const empty = document.getElementById("graphs-empty");
    if (m === "transect") {
      if (graphs) graphs.style.display = "none";
      if (empty) empty.style.display = "none";
      host.style.display = "";
      _fillTransectSel();
      _refreshAndDraw();
    } else {
      host.style.display = "none";
      if (graphs) graphs.style.display = "";
    }
  }

  async function _refreshAndDraw() {
    try { catalog = await ViewerTab.datasetCatalog(); } catch (e) { catalog = []; }
    _scheduleRedraw(true);
  }

  async function _drawTransect() {
    try {
      const n = Objects.byKind("transect").length + 1;
      const obj = await Draw.transect({ name: `Transect ${n}`, color: "#111111" });
      _fillTransectSel();
      if (obj && obj.id) ctrl.transectSel.value = obj.id;
      _scheduleRedraw(true);
    } catch (e) { /* draw cancelled */ }
  }

  function _fillTransectSel() {
    const cur = ctrl.transectSel.value;
    U.clear(ctrl.transectSel);
    const transects = Objects.byKind("transect");
    for (const t of transects) {
      ctrl.transectSel.append(U.el("option", { value: t.id, selected: t.id === cur ? "" : null }, t.name));
    }
    if (!transects.length) {
      ctrl.transectSel.append(U.el("option", { value: "" }, "— draw a transect first —"));
    }
  }

  function _selectedKeys() {
    const avail = catalog.map((s) => s.key);
    const stored = App.state.ui.transectVars;
    if (stored && stored.length) {
      const keep = stored.filter((k) => avail.includes(k));
      if (keep.length) return keep;
    }
    const def = avail.filter((k) => { const v = _varOf(k); return v && (isBed(v) || isWater(v) || isSep(v)); });
    return def.length ? def : avail.slice(0, 1);
  }

  async function _pickVars() {
    if (!catalog.length) { try { catalog = await ViewerTab.datasetCatalog(); } catch (e) { /* */ } }
    if (!catalog.length) { U.toast("Nothing to plot yet — load a layer or run a simulation", "error"); return; }
    const popup = Popup.open({ title: "Transect data", width: 360 });
    const sel = new Set(_selectedKeys());
    const boxes = {};
    const list = U.el("div", { style: "max-height:340px;overflow:auto" });
    for (const s of catalog) {
      const cb = U.el("input", { type: "checkbox", checked: sel.has(s.key) ? "" : null });
      boxes[s.key] = cb;
      list.append(U.el("label", { class: "choice-row", style: "display:flex;gap:6px;align-items:center" },
        cb, U.el("span", {}, s.label)));
    }
    const ok = U.el("button", { class: "primary" }, "Apply");
    ok.addEventListener("click", () => {
      App.state.ui.transectVars = catalog.map((s) => s.key).filter((k) => boxes[k].checked);
      App.touchUi();
      popup.close();
      _scheduleRedraw(true);
    });
    popup.body.append(list,
      U.el("div", { class: "btn-row", style: "justify-content:flex-end;margin-top:8px" }, ok));
  }

  function _scheduleRedraw(rebuild) {
    clearTimeout(redrawTimer);
    redrawTimer = setTimeout(() => _redraw(rebuild).catch(() => {}), rebuild ? 0 : 60);
  }

  function _destroyPlot() { if (plot) { plot.destroy(); plot = null; curSig = ""; } }

  function _message(text) {
    _destroyPlot();
    U.clear(chartEl);
    chartEl.append(U.el("div", { class: "gt-msg" }, text));
  }

  function _color(key, i) {
    const v = _varOf(key);
    if (v && isWater(v)) return "#2a6fb0";
    if (v && isBed(v)) return "#8a6d3b";
    if (v && isSep(v)) return "#111";
    return PALETTE[i % PALETTE.length];
  }

  async function _redraw(rebuild) {
    if (mode !== "transect" || !host) return;
    if (!catalog.length) { _message("Load a layer (or run a simulation), then pick data."); return; }
    const obj = ctrl.transectSel.value && Objects.get(ctrl.transectSel.value);
    if (!obj) { _message("Draw a transect and select it above."); return; }
    const keys = _selectedKeys();
    if (!keys.length) { _message("Pick at least one dataset (Data…)."); return; }

    const prof = await ViewerTab.sampleTransect(obj.coords, keys);
    if (!prof || !prof.dist.length) { _message("Transect is outside the data extent."); return; }

    const data = [prof.dist];
    for (const k of keys) data.push(prof.values[k] || new Array(prof.dist.length).fill(null));

    const sig = keys.join("|");
    if (!rebuild && plot && sig === curSig) { plot.setData(data); return; }

    // (re)build the chart when the set of series changes
    _destroyPlot();
    U.clear(chartEl);
    const byKey = new Map(catalog.map((s) => [s.key, s]));
    let bedIdx = -1;
    keys.forEach((k, i) => { const v = _varOf(k); if (v && isBed(v)) bedIdx = i + 1; });
    const series = [];
    const bands = [];
    keys.forEach((k, i) => {
      const v = _varOf(k);
      const s = { label: (byKey.get(k) || {}).label || k, color: _color(k, i), unit: "m",
        width: (v && isSep(v)) ? 1.5 : 2, spanGaps: false };
      if (v && isBed(v)) { s.fill = SAND; s.fillTo = (u) => u.scales.y.min; }
      if (v && isSep(v)) s.dash = [4, 3];
      series.push(s);
      if (v && isWater(v) && bedIdx > 0) bands.push({ series: [i + 1, bedIdx], fill: SEA });
    });

    plot = TSPlot.create(chartEl, {
      timeBased: false, minSpan: 1, yLabel: "value (m)",
      xFormat: (d) => `${U.fmtNum(d, 1)} m`,
      data, series, bands,
      width: chartEl.clientWidth || 600, height: chartEl.clientHeight || 260,
    });
    curSig = sig;
  }

  return { init };
})();
