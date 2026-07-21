/* Graphs panel: full-height chart "slides" with scroll snapping and a
 * dot navigator. All time-based charts share one time window: zooming
 * any chart (drag) applies to all of them, to the availability chart
 * and highlights the window on the playbar; double-click or the ⤢
 * button resets. Tick positions align to clean hour/day/month
 * boundaries. The playbar clock draws a shared vertical cursor. */
"use strict";

const Graphs = (() => {

  const charts = new Map();   // id -> {plot, wrap, timeBased}
  let timeWindow = null;      // [min, max] epoch or null (= full range)
  let syncing = false;        // guard against setScale recursion

  const Y_AXIS_SIZE = 56;     // shared so all plots align horizontally

  function _container() { return document.getElementById("graphs"); }

  function _updateEmptyNote() {
    const note = document.getElementById("graphs-empty");
    if (note) note.style.display = (charts.size || document.getElementById("avail-chart"))
      ? "none" : "block";
  }

  function _cardHeight() {
    const wrapEl = document.getElementById("graphs-wrap");
    return Math.max(140, (wrapEl ? wrapEl.clientHeight : 240) - 12);
  }

  /* time ticks aligned to clean boundaries */
  const STEPS = [60, 600, 3600, 6 * 3600, 86400, 7 * 86400, 30.44 * 86400, 365.25 * 86400];
  function _timeSplits(u, axisIdx, min, max) {
    const span = max - min;
    let step = STEPS[STEPS.length - 1];
    for (const s of STEPS) {
      if (span / s <= 9) { step = s; break; }
    }
    const splits = [];
    if (step >= 365.25 * 86400 || step >= 30.44 * 86400) {
      // month/year boundaries via calendar arithmetic
      const d0 = new Date(min * 1000);
      let y = d0.getUTCFullYear(), m = step >= 365.25 * 86400 ? 0 : d0.getUTCMonth();
      const stepMonths = step >= 365.25 * 86400 ? 12 : (span / (30.44 * 86400) > 30 ? 6 : 1);
      let t = Date.UTC(y, m, 1) / 1000;
      while (t <= max) {
        if (t >= min) splits.push(t);
        m += stepMonths;
        y += Math.floor(m / 12); m %= 12;
        t = Date.UTC(y, m, 1) / 1000;
      }
    } else {
      const first = Math.ceil(min / step) * step;
      for (let t = first; t <= max; t += step) splits.push(t);
    }
    return splits.length >= 2 ? splits : null;
  }

  function _timeValues(u, ticks) {
    const span = (u.scales.x.max || 1) - (u.scales.x.min || 0);
    return ticks.map((t) => {
      const d = new Date(t * 1000);
      const pad = (v) => String(v).padStart(2, "0");
      if (span > 300 * 86400) {
        return d.getUTCMonth() === 0 && d.getUTCDate() === 1
          ? String(d.getUTCFullYear())
          : `${pad(d.getUTCDate())}-${pad(d.getUTCMonth() + 1)}-${String(d.getUTCFullYear()).slice(2)}`;
      }
      if (span > 3 * 86400) {
        return `${pad(d.getUTCDate())} ${["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul",
          "Aug", "Sep", "Oct", "Nov", "Dec"][d.getUTCMonth()]}`;
      }
      return `${pad(d.getUTCHours())}:${pad(d.getUTCMinutes())}`;
    });
  }

  /* Add (or replace) a chart.
   * spec: {title, timeBased?, data, series, axes?, scales?} */
  function add(id, spec) {
    remove(id);
    const timeBased = spec.timeBased !== false;
    const resetBtn = U.el("button", {
      class: "ghost graph-close", title: "Reset zoom", onclick: () => setTimeWindow(null),
    }, "⤢");
    const wrap = U.el("div", { class: "graph-card", dataset: { graph: id } });
    const head = U.el("div", { class: "graph-head" },
      U.el("span", { class: "graph-title" }, spec.title || id),
      U.el("span", {},
        timeBased ? resetBtn : null,
        U.el("button", { class: "ghost graph-close", title: "Close", onclick: () => remove(id) }, "✕")),
    );
    const plotEl = U.el("div");
    wrap.append(head, plotEl);
    wrap.style.minHeight = `${_cardHeight()}px`;
    _container().append(wrap);

    const width = _container().clientWidth - 24;
    const height = _cardHeight() - 58;

    const axes = spec.axes ? [...spec.axes] : [{}, { size: Y_AXIS_SIZE }];
    if (timeBased) {
      axes[0] = { ...axes[0], splits: _timeSplits, values: _timeValues };
    }
    for (let i = 1; i < axes.length; i += 1) {
      if (axes[i] && axes[i].side !== 1) axes[i] = { size: Y_AXIS_SIZE, ...axes[i] };
    }

    const opts = {
      width,
      height: Math.max(110, height),
      title: undefined,
      scales: spec.scales || {},
      series: spec.series,
      axes,
      legend: { show: spec.legend !== false },
      cursor: { drag: { x: true, y: false } },
      hooks: {
        draw: [(u) => _drawTimeCursor(u, id)],
        setScale: [(u, key) => {
          if (key !== "x" || syncing || !timeBased) return;
          const { min, max } = u.scales.x;
          _onUserZoom(id, min, max);
        }],
      },
      ...spec.uplot,
    };
    const plot = new uPlot(opts, spec.data, plotEl);
    plotEl.addEventListener("dblclick", () => setTimeWindow(null));
    charts.set(id, { plot, wrap, timeBased });
    if (timeBased && timeWindow) {
      syncing = true;
      plot.setScale("x", { min: timeWindow[0], max: timeWindow[1] });
      syncing = false;
    }
    _updateEmptyNote();
    _renderDots();
    return plot;
  }

  /* ---- shared time window ---- */

  function _fullRange(entry) {
    const xs = entry.plot.data[0];
    return xs && xs.length ? [xs[0], xs[xs.length - 1]] : null;
  }

  function _onUserZoom(sourceId, min, max) {
    const src = charts.get(sourceId);
    if (!src) return;
    const full = _fullRange(src);
    // treat a zoom equal to the full extent as "no window"
    if (full && Math.abs(min - full[0]) < 1 && Math.abs(max - full[1]) < 1) {
      timeWindow = null;
    } else {
      timeWindow = [min, max];
    }
    _applyTimeWindow(sourceId);
  }

  function setTimeWindow(window) {
    timeWindow = window;
    _applyTimeWindow(null);
  }

  function _applyTimeWindow(exceptId) {
    syncing = true;
    for (const [id, entry] of charts) {
      if (!entry.timeBased || id === exceptId) continue;
      const range = timeWindow || _fullRange(entry);
      if (range) entry.plot.setScale("x", { min: range[0], max: range[1] });
    }
    syncing = false;
    _updateAvailWindow();
    if (typeof Playbar !== "undefined") Playbar.setViewWindow(timeWindow);
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
      _renderDots();
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
    const cardHeight = _cardHeight();
    const height = Math.max(110, cardHeight - 58);
    for (const entry of charts.values()) {
      entry.wrap.style.minHeight = `${cardHeight}px`;
      const width = container.clientWidth - 24;
      if (width > 50) entry.plot.setSize({ width, height });
    }
    if (availBox) availBox.style.minHeight = `${cardHeight}px`;
  }

  /* ---- slide dot navigator ---- */

  let dotsEl = null;
  let dotObserver = null;

  function _renderDots() {
    const wrapEl = document.getElementById("graphs-wrap");
    if (!wrapEl) return;
    if (!dotsEl) {
      dotsEl = U.el("div", { id: "graph-dots" });
      wrapEl.append(dotsEl);
    }
    U.clear(dotsEl);
    const cards = _container().querySelectorAll(".graph-card");
    if (cards.length < 2) { dotsEl.style.display = "none"; return; }
    dotsEl.style.display = "";
    cards.forEach((card) => {
      const dot = U.el("span", {
        class: "graph-dot",
        title: (card.querySelector(".graph-title") || {}).textContent || "",
      });
      dot.addEventListener("click", () => card.scrollIntoView({ behavior: "smooth" }));
      dot._card = card;
      dotsEl.append(dot);
    });
    if (dotObserver) dotObserver.disconnect();
    dotObserver = new IntersectionObserver((entries) => {
      for (const entry of entries) {
        for (const dot of dotsEl.children) {
          if (dot._card === entry.target) dot.classList.toggle("on", entry.isIntersecting);
        }
      }
    }, { root: wrapEl, threshold: 0.55 });
    for (const card of cards) dotObserver.observe(card);
  }

  /* =================================================================
   * Data availability chart: simulation window, boundary-condition
   * series (with repeat shading), output steps and raw-data surveys on
   * one shared time axis, with the playbar cursor.
   * ================================================================= */

  let availBox = null;
  let availRange = null;   // [t0, t1] epoch
  let availLast = null;    // last payloads for cheap re-render on zoom

  const refreshAvailability = U.debounce(async () => {
    if (!App.state.project) return;
    let conditions = null, output = null, domain = null;
    try { conditions = await Api.get("/api/conditions"); } catch { /* no config */ }
    try { output = await Api.get("/api/output/meta"); } catch { /* none */ }
    try { domain = await Api.get("/api/domain"); } catch { /* none */ }
    if (!conditions) return;
    availLast = { conditions, output, domain };
    _renderAvailability(conditions, output, domain);
  }, 400);

  function _updateAvailWindow() {
    if (availLast) {
      _renderAvailability(availLast.conditions, availLast.output, availLast.domain);
    }
  }

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
    // follow the shared zoom window of the time-based graphs
    if (timeWindow) { [t0, t1] = timeWindow; }
    availRange = [t0, t1];
    const pct = (t) => `${(100 * (t - t0) / (t1 - t0)).toFixed(2)}%`;

    if (availBox) availBox.remove();
    availBox = U.el("div", { class: "graph-card", id: "avail-chart" });
    availBox.style.minHeight = `${_cardHeight()}px`;
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
    _renderDots();
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

  return { init, add, setData, remove, clearAll, has, redrawCursors, resizeAll,
    refreshAvailability, setTimeWindow, timeWindow: () => timeWindow };
})();
