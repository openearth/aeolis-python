/* TSPlot: the one shared time-series chart component (uPlot wrapper).
 *
 * Used by the graphs panel, the wizard previews and the probe charts so
 * every time series in the app looks and behaves the same:
 *
 *   drag        pan (x only, clamped to the data; y refits itself)
 *   wheel       zoom around the cursor
 *   shift+drag  select a range to zoom into
 *   dblclick    reset to the full range
 *
 * The header shows the title, one colored chip per series and - while
 * hovering - the values under the cursor. Ticks align to clean
 * hour/day/month boundaries. Direction-like series render as dots
 * (lines are meaningless across the 360-degree wrap). Gaps in the data
 * (dt >> median dt) are broken instead of bridged.
 */
"use strict";

const TSPlot = (() => {

  const Y_AXIS_SIZE = 56;      // shared so all plots align horizontally
  const MIN_SPAN = 60;         // never zoom below one minute
  const COLORS = ["#0f766e", "#b4423b", "#2f7fe6", "#a034c6", "#e0a020", "#12a5b5"];

  /* ================= calendar-aligned time axis ================= */

  const STEPS = [60, 600, 3600, 6 * 3600, 86400, 7 * 86400, 30.44 * 86400, 365.25 * 86400];

  function timeSplits(u, axisIdx, min, max) {
    const span = max - min;
    let step = STEPS[STEPS.length - 1];
    for (const s of STEPS) {
      if (span / s <= 9) { step = s; break; }
    }
    const splits = [];
    if (step >= 30.44 * 86400) {
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
    // never return null: uPlot feeds the result straight into values()
    return splits;
  }

  const MONTHS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"];

  function timeValues(u, ticks) {
    const span = (u.scales.x.max || 1) - (u.scales.x.min || 0);
    const pad = (v) => String(v).padStart(2, "0");
    return ticks.map((t) => {
      const d = new Date(t * 1000);
      if (span > 300 * 86400) {
        return d.getUTCMonth() === 0 && d.getUTCDate() === 1
          ? String(d.getUTCFullYear())
          : `${MONTHS[d.getUTCMonth()]} ${String(d.getUTCFullYear()).slice(2)}`;
      }
      if (span > 3 * 86400) {
        return `${pad(d.getUTCDate())} ${MONTHS[d.getUTCMonth()]}`;
      }
      return `${pad(d.getUTCHours())}:${pad(d.getUTCMinutes())}`;
    });
  }

  /* Insert nulls where the time step is much larger than the median so
   * gaps are not bridged by long straight lines. Returns new arrays. */
  function withGaps(data, factor = 5) {
    const t = data[0];
    if (!t || t.length < 3) return data;
    const dts = [];
    for (let i = 1; i < t.length; i += 1) dts.push(t[i] - t[i - 1]);
    const sorted = [...dts].sort((a, b) => a - b);
    const median = sorted[Math.floor(sorted.length / 2)] || 1;
    const cut = Math.max(median * factor, 1);
    const gapIdx = [];
    for (let i = 1; i < t.length; i += 1) {
      if (t[i] - t[i - 1] > cut) gapIdx.push(i);
    }
    if (!gapIdx.length) return data;
    const out = data.map(() => []);
    let g = 0;
    for (let i = 0; i < t.length; i += 1) {
      if (g < gapIdx.length && i === gapIdx[g]) {
        // splice a null sample just after the gap start
        out[0].push(t[i - 1] + cut / 2);
        for (let c = 1; c < data.length; c += 1) out[c].push(null);
        g += 1;
      }
      for (let c = 0; c < data.length; c += 1) out[c].push(data[c][i]);
    }
    return out;
  }

  /* ================= instance ================= */

  /* spec: {
   *   title, data: [t, col...],
   *   series: [{label, color?, unit?, points?, width?, scale?}],   // one per column
   *   timeBased?, yRange?: [lo,hi],        // fixed y range (e.g. 0..360)
   *   axes2?: bool,                        // second y-axis (series scale "y2")
   *   height, width,
   *   onWindow?: (win|null, self) => {},   // user pan/zoom notification
   *   drawExtra?: (u) => {},               // e.g. playbar time cursor
   *   actions?: [button nodes]             // extra header buttons (close etc.)
   * } */
  function create(host, spec) {
    const timeBased = spec.timeBased !== false;
    const seriesSpecs = (spec.series || []).map((s, i) => ({
      color: s.color || COLORS[i % COLORS.length], ...s,
    }));

    /* ---- header ---- */
    const chips = [];
    const chipBox = U.el("span", { class: "ts-chips" });
    seriesSpecs.forEach((s) => {
      const chip = U.el("span", { class: "ts-chip" },
        U.el("span", { class: "ts-chip-dot", style: `background:${s.color}` }),
        U.el("span", { class: "ts-chip-label" }, s.label || ""),
        U.el("span", { class: "ts-chip-val" }, ""));
      chips.push(chip);
      chipBox.append(chip);
    });
    const timeReadout = U.el("span", { class: "ts-time muted" }, "");
    const head = U.el("div", { class: "ts-head" },
      ...(spec.actionsLeft || []),
      U.el("span", { class: "ts-title", title: spec.title || "" }, spec.title || ""),
      chipBox,
      U.el("span", { class: "grow" }),
      timeReadout,
      ...(spec.actions || []));
    const plotEl = U.el("div", { class: "ts-plot" });
    const wrap = U.el("div", { class: "ts-wrap" }, head, plotEl);
    host.append(wrap);

    /* ---- ranges & window ---- */
    let window_ = null;         // [t0,t1] or null = full
    let fullOverride = spec.fullRange || null;   // e.g. include the sim window

    function fullRange() {
      if (fullOverride) return fullOverride;
      const xs = plot ? plot.data[0] : spec.data[0];
      if (!xs || !xs.length) return null;
      return [xs[0], xs[xs.length - 1]];
    }

    function clampWindow(win) {
      const full = fullRange();
      if (!full || !win) return win;
      let [lo, hi] = win;
      const span = Math.min(Math.max(hi - lo, MIN_SPAN), full[1] - full[0]);
      lo = U.clamp(lo, full[0], full[1] - span);
      return [lo, lo + span];
    }

    /* ---- y auto-fit (visible window) ---- */
    const yFit = U.debounce(() => {
      if (!plot || spec.yRange) return;
      const { min, max } = plot.scales.x;
      const t = plot.data[0];
      const perScale = {};
      seriesSpecs.forEach((s, i) => {
        const scale = s.scale || "y";
        const col = plot.data[i + 1];
        if (!col) return;
        let acc = perScale[scale];
        if (!acc) acc = perScale[scale] = { lo: Infinity, hi: -Infinity };
        for (let k = 0; k < t.length; k += 1) {
          if (t[k] < min || t[k] > max) continue;
          const v = col[k];
          if (v === null || v === undefined || Number.isNaN(v)) continue;
          if (v < acc.lo) acc.lo = v;
          if (v > acc.hi) acc.hi = v;
        }
      });
      plot.batch(() => {
        for (const [scale, acc] of Object.entries(perScale)) {
          if (!Number.isFinite(acc.lo)) continue;
          let pad = (acc.hi - acc.lo) * 0.08;
          if (pad === 0) pad = Math.abs(acc.hi) * 0.1 || 0.5;
          // yZeroFloor: anchor the axis at 0 for non-negative data (wind,
          // waves, magnitudes) instead of floating the baseline; series
          // that dip below 0 (water level) still show their full range.
          const min = (spec.yZeroFloor && acc.lo >= 0) ? 0 : acc.lo - pad;
          plot.setScale(scale, { min, max: acc.hi + pad });
        }
      });
    }, 120);

    /* ---- apply / report window ---- */
    let applying = false;

    function setWindow(win, silent = true) {
      if (!plot) return;
      window_ = win ? clampWindow(win) : null;
      const range = window_ || fullRange();
      if (!range) return;
      applying = true;
      try {
        plot.setScale("x", { min: range[0], max: range[1] });
      } finally {
        applying = false;
      }
      yFit();
      if (!silent && spec.onWindow) spec.onWindow(window_, api);
    }

    function reportUserWindow(win) {
      const full = fullRange();
      if (full && win && win[1] - win[0] >= full[1] - full[0] - 1) win = null;
      window_ = win ? clampWindow(win) : null;
      if (spec.onWindow) spec.onWindow(window_, api);
    }

    /* ---- interactions on the plot overlay ---- */

    function wireInteractions(u) {
      const over = u.over;
      over.style.cursor = "grab";

      // wheel = zoom around the cursor
      over.addEventListener("wheel", (ev) => {
        ev.preventDefault();
        const full = fullRange();
        if (!full) return;
        const { min, max } = u.scales.x;
        const span = max - min;
        const factor = ev.deltaY < 0 ? 0.8 : 1.25;
        let newSpan = U.clamp(span * factor, MIN_SPAN, full[1] - full[0]);
        const tAt = u.posToVal(ev.offsetX, "x");
        const frac = (tAt - min) / span;
        let lo = tAt - newSpan * frac;
        let hi = lo + newSpan;
        if (newSpan >= full[1] - full[0]) { setWindow(null, false); return; }
        if (lo < full[0]) { hi += full[0] - lo; lo = full[0]; }
        if (hi > full[1]) { lo -= hi - full[1]; hi = full[1]; }
        setWindow([lo, hi], false);
      }, { passive: false });

      // drag = pan; shift+drag = range select
      over.addEventListener("mousedown", (ev) => {
        if (ev.button !== 0) return;
        const full = fullRange();
        if (!full) return;
        const x0 = ev.clientX;
        const { min: m0, max: M0 } = u.scales.x;
        const pxPerSec = u.bbox.width / window.devicePixelRatio / (M0 - m0);

        if (ev.shiftKey) {
          // range-select zoom with a highlight band
          const band = U.el("div", { class: "ts-select" });
          over.append(band);
          const rect = over.getBoundingClientRect();
          const move = (mv) => {
            const a = Math.min(x0, mv.clientX) - rect.left;
            const b = Math.max(x0, mv.clientX) - rect.left;
            band.style.left = `${a}px`;
            band.style.width = `${b - a}px`;
          };
          const up = (mv) => {
            window.removeEventListener("mousemove", move);
            window.removeEventListener("mouseup", up);
            band.remove();
            const a = Math.min(x0, mv.clientX) - rect.left;
            const b = Math.max(x0, mv.clientX) - rect.left;
            if (b - a > 4) {
              const lo = u.posToVal(a, "x");
              const hi = u.posToVal(b, "x");
              setWindow([lo, hi], false);
            }
          };
          window.addEventListener("mousemove", move);
          window.addEventListener("mouseup", up);
          ev.preventDefault();
          return;
        }

        over.style.cursor = "grabbing";
        let moved = false;
        const move = (mv) => {
          const dt = (x0 - mv.clientX) / pxPerSec;
          if (Math.abs(mv.clientX - x0) > 2) moved = true;
          let lo = m0 + dt, hi = M0 + dt;
          if (lo < full[0]) { hi += full[0] - lo; lo = full[0]; }
          if (hi > full[1]) { lo -= hi - full[1]; hi = full[1]; }
          applying = true;
          try {
            u.setScale("x", { min: lo, max: hi });
          } finally {
            applying = false;
          }
        };
        const up = () => {
          window.removeEventListener("mousemove", move);
          window.removeEventListener("mouseup", up);
          over.style.cursor = "grab";
          if (moved) {
            const { min, max } = u.scales.x;
            reportUserWindow([min, max]);
            yFit();
          }
        };
        window.addEventListener("mousemove", move);
        window.addEventListener("mouseup", up);
        ev.preventDefault();
      });

      over.addEventListener("dblclick", () => setWindow(null, false));
    }

    /* ---- cursor readout ---- */

    function onCursor(u) {
      const idx = u.cursor.idx;
      if (idx === null || idx === undefined) {
        timeReadout.textContent = "";
        chips.forEach((chip) => { chip.querySelector(".ts-chip-val").textContent = ""; });
        return;
      }
      const t = u.data[0][idx];
      timeReadout.textContent = U.fmtDate(t);
      seriesSpecs.forEach((s, i) => {
        const v = u.data[i + 1] ? u.data[i + 1][idx] : null;
        chips[i].querySelector(".ts-chip-val").textContent =
          (v === null || v === undefined) ? "–" : U.fmtNum(v, 3) + (s.unit ? ` ${s.unit}` : "");
      });
    }

    /* ---- build uPlot ---- */

    const uSeries = [{}];
    seriesSpecs.forEach((s) => {
      const entry = {
        label: s.label || "",
        stroke: s.color,
        width: s.width ?? 1.5,
        scale: s.scale || "y",
        spanGaps: Boolean(s.spanGaps),
        points: { show: Boolean(s.points), size: 4, fill: s.color, stroke: s.color },
      };
      if (s.points) entry.paths = () => null;   // dots only (e.g. direction)
      uSeries.push(entry);
    });

    const axes = [
      timeBased
        ? { splits: timeSplits, values: timeValues, grid: { stroke: "rgba(28,36,48,.07)" } }
        : { grid: { stroke: "rgba(28,36,48,.07)" } },
      { size: Y_AXIS_SIZE, scale: "y", grid: { stroke: "rgba(28,36,48,.07)" },
        label: spec.yLabel || null, labelSize: spec.yLabel ? 14 : 0 },
    ];
    if (spec.axes2) {
      axes.push({ size: Y_AXIS_SIZE - 8, scale: "y2", side: 1, grid: { show: false } });
    }

    const scales = { x: { time: false } };
    scales.y = spec.yRange
      ? { range: () => spec.yRange }
      : { range: (u, lo, hi) => [lo, hi] };
    if (spec.axes2) scales.y2 = { range: (u, lo, hi) => [lo, hi] };

    const opts = {
      width: spec.width || 600,
      height: spec.height || 200,
      series: uSeries,
      axes,
      scales,
      legend: { show: false },
      cursor: {
        drag: { x: false, y: false, setScale: false },
        sync: timeBased ? { key: "aeolis-ts" } : undefined,
        points: { size: 6 },
      },
      hooks: {
        ready: [(u) => wireInteractions(u)],
        setCursor: [(u) => onCursor(u)],
        setScale: [(u, key) => { if (key === "x" && !applying) yFit(); }],
        draw: [(u) => { if (spec.drawExtra) spec.drawExtra(u); }],
      },
    };

    const plot = new uPlot(opts, spec.data, plotEl);
    yFit();

    /* ---- public api ---- */
    const api = {
      plot,
      wrap,
      head,
      setWindow,
      window: () => window_,
      fullRange,
      setData(data) {
        plot.setData(data);
        if (window_) setWindow(window_);
        else yFit();
      },
      setFullRange(range) {
        fullOverride = range || null;
        if (!window_) setWindow(null);
      },
      setSize({ width, height }) {
        plot.setSize({ width, height });
      },
      redraw() { plot.redraw(false, true); },
      destroy() {
        plot.destroy();
        wrap.remove();
      },
    };
    return api;
  }

  return { create, withGaps, timeSplits, timeValues,
    timeAxis: () => ({ splits: timeSplits, values: timeValues }),
    Y_AXIS_SIZE, COLORS };
})();
