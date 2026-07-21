/* Graphs panel: a dynamic stack of uPlot charts in the bottom-right
 * area. Charts register with add(); the shared time cursor (from the
 * playbar clock) is drawn as a vertical line in every time-based chart. */
"use strict";

const Graphs = (() => {

  const charts = new Map();   // id -> {plot, opts, wrap, timeBased}

  function _container() { return document.getElementById("graphs"); }

  function _updateEmptyNote() {
    const note = document.getElementById("graphs-empty");
    if (note) note.style.display = charts.size ? "none" : "block";
  }

  /* Add (or replace) a chart.
   * spec: {title, height?, timeBased?, data: aligned uPlot data,
   *        series: uPlot series defs, axes?, scales?} */
  function add(id, spec) {
    remove(id);
    const wrap = U.el("div", { class: "graph-card", dataset: { graph: id } });
    const head = U.el("div", { class: "graph-head" },
      U.el("span", { class: "graph-title" }, spec.title || id),
      U.el("button", { class: "ghost graph-close", title: "Close", onclick: () => remove(id) }, "✕"),
    );
    const plotEl = U.el("div");
    wrap.append(head, plotEl);
    _container().append(wrap);

    const width = plotEl.clientWidth || _container().clientWidth - 20;
    const opts = {
      width,
      height: spec.height || 160,
      title: undefined,
      scales: spec.scales || {},
      series: spec.series,
      axes: spec.axes || [
        {}, { size: 55 },
      ],
      legend: { show: spec.legend !== false },
      cursor: { drag: { x: true, y: false } },
      hooks: {
        draw: [(u) => _drawTimeCursor(u, id)],
      },
      ...spec.uplot,
    };
    const plot = new uPlot(opts, spec.data, plotEl);
    charts.set(id, { plot, wrap, timeBased: spec.timeBased !== false });
    _updateEmptyNote();
    return plot;
  }

  function setData(id, data) {
    const entry = charts.get(id);
    if (entry) entry.plot.setData(data);
  }

  function remove(id) {
    const entry = charts.get(id);
    if (entry) {
      entry.plot.destroy();
      entry.wrap.remove();
      charts.delete(id);
      _updateEmptyNote();
    }
  }

  function clearAll() {
    for (const id of Array.from(charts.keys())) remove(id);
  }

  function has(id) { return charts.has(id); }

  /* ---- shared time cursor ---- */

  function _drawTimeCursor(u, id) {
    const entry = charts.get(id);
    const t = App.state.clock.t;
    if (!entry || !entry.timeBased || !Number.isFinite(t)) return;
    const scale = u.scales.x;
    if (t < scale.min || t > scale.max) return;
    const x = u.valToPos(t, "x", true);
    const ctx = u.ctx;
    ctx.save();
    ctx.strokeStyle = "#0f766e";
    ctx.lineWidth = 1.5;
    ctx.setLineDash([5, 4]);
    ctx.beginPath();
    ctx.moveTo(x, u.bbox.top);
    ctx.lineTo(x, u.bbox.top + u.bbox.height);
    ctx.stroke();
    ctx.restore();
  }

  function redrawCursors() {
    for (const entry of charts.values()) {
      if (entry.timeBased) entry.plot.redraw(false, true);
    }
  }

  function resizeAll() {
    const container = _container();
    if (!container) return;
    for (const entry of charts.values()) {
      const width = container.clientWidth - 20;
      if (width > 50) entry.plot.setSize({ width, height: entry.plot.height });
    }
  }

  function init() {
    App.on("clock-tick", redrawCursors);
    window.addEventListener("resize", U.debounce(resizeAll, 150));
    _updateEmptyNote();
  }

  return { init, add, setData, remove, clearAll, has, redrawCursors, resizeAll };
})();
