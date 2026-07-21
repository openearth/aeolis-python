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

  /* =================================================================
   * Data availability chart: simulation window, boundary-condition
   * series (with repeat shading), output steps and raw-data surveys on
   * one shared time axis, with the playbar cursor.
   * ================================================================= */

  let availBox = null;
  let availRange = null;   // [t0, t1] epoch

  const refreshAvailability = U.debounce(async () => {
    if (!App.state.project) return;
    let conditions = null, output = null, domain = null;
    try { conditions = await Api.get("/api/conditions"); } catch { /* no config */ }
    try { output = await Api.get("/api/output/meta"); } catch { /* none */ }
    try { domain = await Api.get("/api/domain"); } catch { /* none */ }
    if (!conditions) return;
    _renderAvailability(conditions, output, domain);
  }, 400);

  function _renderAvailability(conditions, output, domain) {
    const rows = [];
    // prefer live (possibly unsaved) config values so edits in the
    // Settings tab update the chart immediately
    let refEpoch = conditions.refdate_epoch;
    let tstart = conditions.tstart || 0;
    let tstop = conditions.tstop || 0;
    const cfg = App.state.config;
    if (cfg && cfg.refdate) {
      const raw = String(cfg.refdate);
      const parsed = Date.parse(raw.replace(" ", "T") + (raw.length <= 16 ? ":00Z" : "Z"));
      if (Number.isFinite(parsed)) refEpoch = parsed / 1000;
      if (Number.isFinite(cfg.tstart)) tstart = cfg.tstart;
      if (Number.isFinite(cfg.tstop)) tstop = cfg.tstop;
    }
    const simT0 = refEpoch + tstart;
    const simT1 = refEpoch + tstop;
    rows.push({ label: "Simulation", spans: [[simT0, simT1, "sim"]] });

    for (const [kind, title] of [["wind", "Wind"], ["tide", "Water levels"], ["wave", "Waves"]]) {
      const info = conditions.kinds[kind];
      if (!info || !info.series) continue;
      const s = info.series;
      const spans = [[s.t0_epoch, s.t1_epoch, "data"]];
      if (s.t1_epoch < simT1 - 1) spans.push([s.t1_epoch, simT1, "repeat"]);
      rows.push({ label: title, spans });
    }

    if (output && output.exists && output.times_epoch.length) {
      rows.push({
        label: "Output",
        spans: [[output.times_epoch[0], output.times_epoch[output.times_epoch.length - 1], "output"]],
        marks: output.times_epoch,
      });
    }

    if (domain && domain.entries && domain.entries.length) {
      const marks = [];
      for (const entry of domain.entries) {
        if (entry.date) marks.push(Date.parse(entry.date) / 1000);
        else if (entry.year) marks.push(Date.UTC(entry.year, 6, 1) / 1000);
      }
      if (marks.length) rows.push({ label: "Surveys (raw data)", marks });
    }

    // shared scale
    let t0 = simT0, t1 = simT1;
    for (const row of rows) {
      for (const [a, b] of row.spans || []) { t0 = Math.min(t0, a); t1 = Math.max(t1, b); }
      for (const m of row.marks || []) { t0 = Math.min(t0, m); t1 = Math.max(t1, m); }
    }
    const pad = (t1 - t0) * 0.03 || 1;
    t0 -= pad; t1 += pad;
    availRange = [t0, t1];
    const pct = (t) => `${(100 * (t - t0) / (t1 - t0)).toFixed(2)}%`;

    if (availBox) availBox.remove();
    availBox = U.el("div", { class: "graph-card", id: "avail-chart" });
    availBox.append(U.el("div", { class: "graph-head" },
      U.el("span", { class: "graph-title" }, "Data availability"),
      U.el("span", { class: "muted", style: "font-size:11px" },
        `${U.fmtDate(t0)} — ${U.fmtDate(t1)}`)));

    const labels = U.el("div", { class: "avail-labels" });
    const tracks = U.el("div", { class: "avail-tracks" });
    for (const row of rows) {
      labels.append(U.el("div", { class: "avail-label" }, row.label));
      const track = U.el("div", { class: "avail-track" });
      for (const [a, b, cls] of row.spans || []) {
        const left = pct(a);
        const width = `${(100 * (b - a) / (t1 - t0)).toFixed(2)}%`;
        track.append(U.el("div", {
          class: `avail-span ${cls}`,
          style: `left:${left};width:${width}`,
          title: `${U.fmtDate(a)} — ${U.fmtDate(b)}${cls === "repeat" ? " (repeated)" : ""}`,
        }));
      }
      for (const m of row.marks || []) {
        track.append(U.el("div", {
          class: "avail-mark", style: `left:${pct(m)}`, title: U.fmtDate(m),
        }));
      }
      tracks.append(track);
    }
    tracks.append(U.el("div", { class: "avail-cursor", id: "avail-cursor" }));
    const body = U.el("div", { class: "avail-body" }, labels, tracks);
    availBox.append(body);
    _container().prepend(availBox);
    _updateEmptyNote2();
    _updateAvailCursor();
  }

  function _updateAvailCursor() {
    const cursor = document.getElementById("avail-cursor");
    if (!cursor || !availRange) return;
    const t = App.state.clock.t;
    const [t0, t1] = availRange;
    if (!Number.isFinite(t) || t < t0 || t > t1) {
      cursor.style.display = "none";
      return;
    }
    cursor.style.display = "";
    cursor.style.left = `${(100 * (t - t0) / (t1 - t0)).toFixed(2)}%`;
  }

  function _updateEmptyNote2() {
    const note = document.getElementById("graphs-empty");
    if (note) note.style.display = (charts.size || availBox) ? "none" : "block";
  }

  function init() {
    App.on("clock-tick", () => { redrawCursors(); _updateAvailCursor(); });
    window.addEventListener("resize", U.debounce(resizeAll, 150));
    App.on("project", refreshAvailability);
    App.on("config-changed", refreshAvailability);
    App.on("run-finished", refreshAvailability);
    _updateEmptyNote();
  }

  return { init, add, setData, remove, clearAll, has, redrawCursors, resizeAll, refreshAvailability };
})();
