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
  const MAX_SPAN = 100 * 365.25 * 86400;   // generous zoom-out wall (~100 years)
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
    const tooltip = U.el("div", { class: "ts-tooltip", style: "display:none" });
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

    // pan/zoom may roam a bit beyond the data (and past the simulation start)
    // so the user isn't hard-stopped at the edges. Zoom-out is capped at a
    // generous MAX_SPAN (~100 years) rather than snapping back to the data
    // range, so it just hits a wall instead of auto-readjusting.
    const OUTER_MARGIN = 0.4;
    function _maxSpan() {
      const full = fullRange();
      const dataSpan = full ? full[1] - full[0] : 0;
      // at least a ~100-year window for typical data; wider only if the data
      // itself already spans more than that
      return Math.max(dataSpan * (1 + 2 * OUTER_MARGIN), MAX_SPAN);
    }

    // Clamp both the span (MIN_SPAN … MAX_SPAN) and the position: when zoomed
    // in the window is kept overlapping the data (± a small margin) so the
    // data can never be panned entirely off-screen; when zoomed out wider than
    // the data it is centred on the data.
    function clampWindow(win) {
      const full = fullRange();
      if (!full || !win) return win;
      const dataSpan = full[1] - full[0];
      const marginAbs = dataSpan * OUTER_MARGIN;
      const span = U.clamp(win[1] - win[0], MIN_SPAN, _maxSpan());
      const loMin = full[0] - marginAbs;
      const loMax = full[1] + marginAbs - span;
      let lo;
      if (loMax >= loMin) lo = U.clamp(win[0], loMin, loMax);
      else lo = (full[0] + full[1]) / 2 - span / 2;   // wider than data: centre
      return [lo, lo + span];
    }

    /* ---- y auto-fit (visible window) ----
     * Fits once after the interaction settles (trailing debounce). A small
     * hysteresis skips near-identical re-fits, so the axis does not jitter
     * smaller/larger/smaller as the zoom-aware refetch swaps denser data in
     * for the same window. */
    let yApplied = {};   // scale -> {min, max} last set
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
          const nmin = (spec.yZeroFloor && acc.lo >= 0) ? 0 : acc.lo - pad;
          const nmax = acc.hi + pad;
          const prev = yApplied[scale];
          if (prev) {
            const tol = Math.max(1e-9, prev.max - prev.min) * 0.06;
            if (Math.abs(nmin - prev.min) < tol && Math.abs(nmax - prev.max) < tol) continue;
          }
          yApplied[scale] = { min: nmin, max: nmax };
          plot.setScale(scale, { min: nmin, max: nmax });
        }
      });
    }, 160);

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
      // zoom-out hits the MAX_SPAN wall (via clampWindow) rather than snapping
      // back to the full range — no auto-readjust. "Home"/dblclick still reset.
      window_ = win ? clampWindow(win) : null;
      if (spec.onWindow) spec.onWindow(window_, api);
    }

    /* ---- interactions on the plot overlay ---- */

    function wireInteractions(u) {
      // listen on the whole plot (u.root), not just the plotting overlay, so
      // zoom/pan work identically when the pointer is over the availability
      // band above the series. x-pixels are measured from the plot overlay's
      // left edge so posToVal stays correct wherever the event started.
      const root = u.root;
      root.style.cursor = "grab";
      const pxOf = (ev) => ev.clientX - u.over.getBoundingClientRect().left;

      // wheel = zoom around the cursor
      root.addEventListener("wheel", (ev) => {
        ev.preventDefault();
        if (!fullRange()) return;
        const { min, max } = u.scales.x;
        const span = max - min;
        const factor = ev.deltaY < 0 ? 0.8 : 1.25;
        const newSpan = U.clamp(span * factor, MIN_SPAN, _maxSpan());
        const tAt = u.posToVal(pxOf(ev), "x");
        const frac = (tAt - min) / span;
        const lo = tAt - newSpan * frac;
        // clampWindow keeps the position sane and enforces the 100-year wall
        setWindow([lo, lo + newSpan], false);
      }, { passive: false });

      // drag = pan; shift+drag = range select
      root.addEventListener("mousedown", (ev) => {
        if (ev.button !== 0) return;
        if (!fullRange()) return;
        const x0 = ev.clientX;
        const { min: m0, max: M0 } = u.scales.x;
        const pxPerSec = u.bbox.width / window.devicePixelRatio / (M0 - m0);

        if (ev.shiftKey) {
          // range-select zoom with a highlight band
          const band = U.el("div", { class: "ts-select" });
          u.over.append(band);
          const rect = u.over.getBoundingClientRect();
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

        root.style.cursor = "grabbing";
        let moved = false;
        const move = (mv) => {
          const dt = (x0 - mv.clientX) / pxPerSec;
          if (Math.abs(mv.clientX - x0) > 2) moved = true;
          // clampWindow keeps the data in view (can't pan it off-screen)
          const win = clampWindow([m0 + dt, M0 + dt]);
          applying = true;
          try {
            u.setScale("x", { min: win[0], max: win[1] });
          } finally {
            applying = false;
          }
        };
        const up = () => {
          window.removeEventListener("mousemove", move);
          window.removeEventListener("mouseup", up);
          root.style.cursor = "grab";
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

      root.addEventListener("dblclick", () => setWindow(null, false));
    }

    /* ---- cursor readout ---- */

    function onCursor(u) {
      const idx = u.cursor.idx;
      const left = u.cursor.left, top = u.cursor.top;
      if (idx === null || idx === undefined || left < 0 || top < 0) {
        tooltip.style.display = "none";
        chips.forEach((chip) => { chip.querySelector(".ts-chip-val").textContent = ""; });
        return;
      }
      const t = u.data[0][idx];
      // floating tooltip: time (with minutes) + each series' value
      let html = `<div class="tt-time">${U.fmtDate(t)}</div>`;
      seriesSpecs.forEach((s, i) => {
        const v = u.data[i + 1] ? u.data[i + 1][idx] : null;
        const val = (v === null || v === undefined) ? "–"
          : U.fmtNum(v, 3) + (s.unit ? ` ${s.unit}` : "");
        html += `<div class="tt-row"><span class="tt-dot" style="background:${s.color}"></span>`
          + `<span class="tt-lab">${s.label || ""}</span><b>${val}</b></div>`;
        chips[i].querySelector(".ts-chip-val").textContent = val;
      });
      tooltip.innerHTML = html;
      tooltip.style.display = "";
      // position near the cursor, flipping away from the edges
      const ow = u.over ? u.over.clientWidth : 0;
      const oh = u.over ? u.over.clientHeight : 0;
      const tw = tooltip.offsetWidth, th = tooltip.offsetHeight;
      let x = left + 14, y = top + 12;
      if (x + tw > ow) x = left - 14 - tw;
      if (y + th > oh) y = Math.max(2, oh - th - 2);
      tooltip.style.left = `${Math.max(2, x)}px`;
      tooltip.style.top = `${y}px`;
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

    // uPlot draws axis tick labels on the canvas, so CSS can't colour them:
    // read the theme tokens and pass explicit stroke/font (default is black,
    // invisible on the dark album themes). A theme-agnostic grid stroke keeps
    // gridlines faint on both light and dark backgrounds.
    const css = getComputedStyle(document.documentElement);
    const axisStroke = (css.getPropertyValue("--muted") || "#8a95a5").trim();
    const axisFont = `10.5px ${(css.getPropertyValue("--font-ui") || "system-ui").trim()}`;
    const gridStroke = "rgba(128,128,128,.18)";
    // no tick marks (they protrude below the baseline); the x-grid is drawn
    // full-height in drawExtra so it crosses the availability band too
    const axisCommon = { stroke: axisStroke, font: axisFont, ticks: { show: false } };

    const axes = [
      timeBased
        ? { ...axisCommon, size: 26, gap: 6, splits: timeSplits, values: timeValues, grid: { show: false } }
        : { ...axisCommon, size: 26, gap: 6, grid: { show: false } },
      { ...axisCommon, size: Y_AXIS_SIZE, scale: "y", grid: { stroke: gridStroke },
        label: spec.yLabel || null, labelSize: spec.yLabel ? 14 : 0 },
    ];
    if (spec.axes2) {
      axes.push({ ...axisCommon, size: Y_AXIS_SIZE - 8, scale: "y2", side: 1, grid: { show: false } });
    }

    const scales = { x: { time: false } };
    scales.y = spec.yRange
      ? { range: () => spec.yRange }
      : { range: (u, lo, hi) => [lo, hi] };
    if (spec.axes2) scales.y2 = { range: (u, lo, hi) => [lo, hi] };

    const opts = {
      width: spec.width || 600,
      height: spec.height || 200,
      // reserve a band BELOW the plot (under the x-axis labels) for the
      // availability lanes (drawn in drawExtra), so band + series share one
      // canvas and one x-scale
      padding: [6, 12, spec.bottomBand || 4, 4],
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
        ready: [(u) => { wireInteractions(u); u.over.appendChild(tooltip); }],
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
