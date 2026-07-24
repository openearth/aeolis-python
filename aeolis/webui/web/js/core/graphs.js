/* Graphs panel: ONE time-series chart + a data-availability strip that
 * share the same time axis.
 *
 * Data sources (boundary-condition series, raw downloads, output probe
 * series) register here; the user picks which to plot via the series
 * picker. Series with the same unit can be combined in one plot (chips
 * act as the legend); a series with a different unit cannot be added
 * until the conflicting ones are removed. Input series shorter than the
 * simulation are shown with their cyclic repetition (the way AeoLiS
 * wraps them), shaded as "repeated".
 *
 * The availability strip can be hidden with a header button; the chart
 * takes the freed space. The playbar clock draws a shared cursor. */
"use strict";

const Graphs = (() => {

  // key -> {group, label, unit, points?, data? [t,v], load? ()=>Promise,
  //         range? [t0,t1], repeatFrom? epoch}
  const sources = new Map();
  let selection = [];        // ordered source keys currently plotted
  let ts = null;             // TSPlot instance
  let card = null;
  let chartHost = null;
  let availLast = null;       // last {conditions, output, domain} payloads
  let timeWindow = null;     // [t0,t1] epoch or null (= full range)
  let pickerEl = null;
  let building = false;
  let rebuildQueued = false;

  function _wrapEl() { return document.getElementById("graphs-wrap"); }
  function _container() { return document.getElementById("graphs"); }

  function _selValid() { return selection.filter((k) => sources.has(k)); }

  function _unitOf(key) {
    const s = sources.get(key);
    return s ? (s.unit || "") : "";
  }

  /* The direction source that shares this magnitude's origin (same key
   * prefix, e.g. cond-wind-0 ↔ cond-wind-1, rawcond-<id>-0 ↔ …-1). */
  function _companionDir(magKey) {
    const base = magKey.replace(/-\d+$/, "");
    for (const [k, s] of sources) {
      if (k !== magKey && k.replace(/-\d+$/, "") === base && /°|deg/.test(s.unit || "")) {
        return k;
      }
    }
    return null;
  }

  function _currentUnit() {
    const sel = _selValid();
    return sel.length ? _unitOf(sel[0]) : null;
  }

  function _persistSelection() {
    App.state.ui.graphSeries = [...selection];
    App.touchUi();
  }

  /* ================= source registry ================= */

  // registration happens in bursts (project open, tab refresh):
  // coalesce the expensive chart rebuilds
  const scheduleRebuild = U.debounce(() => _rebuild(), 80);

  /* def: {group, label, unit, points?, data?, load?, range?, repeatFrom?} */
  function registerSource(key, def, opts = {}) {
    const existing = sources.get(key);
    sources.set(key, { ...existing, ...def });
    const saved = App.state.ui.graphSeries || [];
    if (opts.select) {
      select(key);
      return;
    }
    if (!selection.includes(key) && saved.includes(key)) {
      // restore a persisted selection when its source (re)appears
      const unit = _currentUnit();
      if (unit === null || unit === (def.unit || "")) selection.push(key);
    }
    scheduleRebuild();
  }

  function unregisterSource(key) {
    sources.delete(key);
    if (selection.includes(key)) {
      selection = selection.filter((k) => k !== key);
    }
    scheduleRebuild();
  }

  function clearAll() {
    sources.clear();
    selection = [];
    availLast = null;
    timeWindow = null;
    scheduleRebuild();
  }

  /* Add a series; on unit conflict the selection is replaced. */
  function select(key) {
    if (!sources.has(key)) return;
    const unit = _currentUnit();
    if (unit !== null && unit !== _unitOf(key)) {
      selection = [key];
      U.toast(`Different unit [${_unitOf(key) || "–"}] — replaced the plotted series`, "");
    } else if (!selection.includes(key)) {
      selection.push(key);
    }
    _persistSelection();
    _rebuild();
  }

  function toggle(key) {
    if (selection.includes(key)) {
      selection = selection.filter((k) => k !== key);
      _persistSelection();
      _rebuild();
    } else {
      select(key);
    }
  }

  /* ================= shared time window ================= */

  function _fullWindow() {
    let t0 = Infinity, t1 = -Infinity;
    for (const key of _selValid()) {
      const s = sources.get(key);
      const r = s.range || (s.data ? [s.data[0][0], s.data[0][s.data[0].length - 1]] : null);
      if (r) { t0 = Math.min(t0, r[0]); t1 = Math.max(t1, r[1]); }
    }
    if (availLast) {
      const r = _availFullRange(availLast);
      if (r) { t0 = Math.min(t0, r[0]); t1 = Math.max(t1, r[1]); }
    }
    if (!Number.isFinite(t0) || t1 <= t0) return null;
    return [t0, t1];
  }

  function setTimeWindow(win) {
    timeWindow = win;
    if (ts) ts.setWindow(win);   // redraws the canvas (band + series) together
    _reloadWindow(win);
    if (typeof Playbar !== "undefined") Playbar.setViewWindow(win);
  }

  /* Zoom-aware refetch: pull denser data for the visible window from the
   * backend (which decimates within [t0,t1]) and swap it in, so zooming
   * reveals the true cadence instead of the coarse full-range decimation. */
  const _reloadWindow = U.debounce(async (win) => {
    if (!ts) return;
    const plotted = _selValid().filter((k) => {
      const s = sources.get(k);
      return s && s.data && typeof s.loadWindow === "function";
    });
    if (!plotted.length) return;
    await Promise.all(plotted.map(async (k) => {
      const s = sources.get(k);
      if (!win) { s.windowData = null; return; }
      try {
        const d = await s.loadWindow(win[0], win[1]);
        s.windowData = (d && d[0] && d[0].length) ? d : null;
      } catch { s.windowData = null; }
    }));
    _rebuild();
  }, 260);

  function _zoomBy(factor) {
    const full = _fullWindow();
    if (!full) return;
    const win = timeWindow || full;
    const mid = (win[0] + win[1]) / 2;
    const half = Math.max(30, (win[1] - win[0]) / 2 * factor);
    let lo = mid - half, hi = mid + half;
    if (hi - lo >= full[1] - full[0]) { setTimeWindow(null); return; }
    if (lo < full[0]) { hi += full[0] - lo; lo = full[0]; }
    if (hi > full[1]) { lo -= hi - full[1]; hi = full[1]; }
    setTimeWindow([lo, hi]);
  }

  /* ================= chart build ================= */

  function _plotSize() {
    const wrap = _wrapEl();
    const width = Math.max(220, _container().clientWidth - 36);
    // the availability band lives INSIDE the canvas (uPlot top padding), so
    // the plot keeps the full height — no separate strip to subtract
    const height = Math.max(90, (wrap ? wrap.clientHeight : 240) - 40);
    return { width, height };
  }

  function _availVisible() { return App.state.ui.graphAvail !== false; }
  // height (CSS px) reserved at the BOTTOM of the canvas (below the x-axis
  // labels) for the availability band when visible; a small pad otherwise
  const BAND_PX = 70;
  function _bandPx() { return _availVisible() && availLast ? BAND_PX : 4; }

  async function _rebuild() {
    if (building) { rebuildQueued = true; return; }
    building = true;
    try {
      await _rebuildNow();
    } finally {
      building = false;
      if (rebuildQueued) { rebuildQueued = false; _rebuild(); }
    }
  }

  async function _rebuildNow() {
    const container = _container();
    if (!container) return;

    const sel = _selValid();
    // lazy-load series data
    for (const key of sel) {
      const s = sources.get(key);
      if (!s.data && s.load) {
        try {
          s.data = await s.load();
          if (s.data && !s.range) {
            const t = s.data[0];
            s.range = t.length ? [t[0], t[t.length - 1]] : null;
          }
        } catch (err) {
          console.warn(`series ${key} failed: ${err.message}`);
          s.error = err.message;
        }
      }
    }
    const plotted = sel.filter((k) => sources.get(k).data);

    if (ts) { ts.destroy(); ts = null; }
    if (card) { card.remove(); card = null; }

    const note = document.getElementById("graphs-empty");
    if (!sources.size && !availLast) {
      if (note) note.style.display = "block";
      return;
    }
    if (note) note.style.display = "none";

    card = U.el("div", { class: "graph-card graph-single" });
    container.append(card);
    chartHost = U.el("div", { class: "graph-chart" });
    card.append(chartHost);

    /* ---- header buttons: variable/source picker on the LEFT; only a
     *      "home" (reset zoom) on the RIGHT. Zoom/step buttons removed —
     *      scroll-to-zoom + drag-to-pan replace them. The availability
     *      show/hide lives on the strip itself (logically grouped). ---- */
    const addBtn = U.el("button", { class: "gbtn gs-add gs-tools-left",
      title: "Choose variable & sources to plot" },
      U.icon("chart", 13), " series ▾");
    addBtn.addEventListener("click", (ev) => { ev.stopPropagation(); _togglePicker(addBtn); });

    const homeBtn = U.el("button", { class: "gbtn",
      title: "Reset to the full time range (or double-click the chart)",
      onclick: () => setTimeWindow(null) }, U.icon("home", 13));
    // the availability show/hide eye lives on the band itself now (next to the
    // lane labels), added as an overlay after the chart is built
    const rightActions = [homeBtn];
    const leftActions = [addBtn];
    // (the windrose lives in the topbar now — one entry point with a source
    // + period chooser — so no per-graph windrose button here.)

    if (!plotted.length) {
      chartHost.append(U.el("div", { class: "graph-pick-note" },
        U.el("span", { class: "muted" }, "No series plotted — "),
        U.el("button", { class: "gbtn", onclick: () => _togglePicker(addBtn) }, "＋ choose series"),
      ));
      // still show the header row so the buttons are reachable
      chartHost.prepend(U.el("div", { class: "ts-head" },
        ...leftActions,
        U.el("span", { class: "ts-title" }, "Time series"),
        U.el("span", { class: "grow" }), ...rightActions));
      return;
    }

    /* ---- assemble uPlot data ---- */
    const specs = plotted.map((key, i) => {
      const s = sources.get(key);
      return {
        key,
        label: s.label,
        unit: s.unit || "",
        color: _seriesColor(i),
        points: Boolean(s.points),
        spanGaps: plotted.length > 1 && !s.points,
      };
    });
    // use the zoom-windowed data when present (denser detail on zoom-in)
    const dataOf = (k) => sources.get(k).windowData || sources.get(k).data;
    let data;
    if (plotted.length === 1) {
      data = TSPlot.withGaps(dataOf(plotted[0]));
    } else {
      // align different time bases on one x array (uPlot.join)
      data = uPlot.join(plotted.map(dataOf));
    }

    const unit = specs[0].unit;
    const isDir = specs.every((s) => /°|deg/.test(s.unit)) &&
      plotted.every((k) => sources.get(k).points);
    const title = plotted.length === 1
      ? specs[0].label + (unit ? ` [${unit}]` : "")
      : (unit ? `[${unit}]` : "Time series");

    const { width, height } = _plotSize();
    ts = TSPlot.create(chartHost, {
      title,
      data,
      series: specs,
      yRange: isDir ? [0, 360] : null,
      yZeroFloor: !isDir,
      fullRange: _fullWindow(),
      width, height,
      bottomBand: _bandPx(),
      actionsLeft: leftActions,
      actions: rightActions,
      onWindow: (win) => {
        timeWindow = win;
        _reloadWindow(win);
        if (typeof Playbar !== "undefined") Playbar.setViewWindow(win);
      },
      // everything below the series lives in one canvas: full-height x-grid,
      // the availability band, sim start/finish lines and the playbar cursor
      drawExtra: (u) => {
        _drawGridLines(u);
        _drawAvailBand(u);
        _drawRepeats(u, plotted);
        _drawBaseline(u);
        _drawSimLines(u);
        _drawTimeCursor(u);
      },
    });
    if (timeWindow) ts.setWindow(timeWindow);

    // availability show/hide eye, overlaid at the left of the chart just
    // ABOVE the band's first ("Simulation") lane label
    chartHost.style.position = "relative";
    const availBtn = U.el("span", {
      class: `eye avail-eye ${_availVisible() ? "" : "off"}`,
      title: _availVisible() ? "Hide the data-availability band" : "Show the data-availability band",
      onclick: () => { App.state.ui.graphAvail = !_availVisible(); App.touchUi(); _rebuild(); },
    }, "👁");
    // sit just above the band top when the band is shown, else near the base
    availBtn.style.bottom = `${(_availVisible() && availLast) ? BAND_PX + 6 : 6}px`;
    chartHost.append(availBtn);

    // the first build happens before the flex layout has settled, so the
    // plot can render at the wrong height; re-measure on the next frame
    requestAnimationFrame(() => { if (ts) resizeAll(); });
  }

  /* ================= series picker ================= */

  let pickerAnchor = null;

  function _closePicker() {
    if (pickerEl) { pickerEl.remove(); pickerEl = null; }
    pickerAnchor = null;
    document.removeEventListener("mousedown", _pickerOutside, true);
  }

  function _pickerOutside(ev) {
    // clicks on the anchor button are handled by its own click handler;
    // closing here too would make the button reopen instead of toggle
    if (pickerAnchor && pickerAnchor.contains(ev.target)) return;
    if (pickerEl && !pickerEl.contains(ev.target)) _closePicker();
  }

  function _togglePicker(anchor) {
    if (pickerEl) { _closePicker(); return; }
    pickerAnchor = anchor;
    // fixed to the viewport (the graphs panel clips its overflow); the
    // panel is short, so open upward from the button, over the map
    pickerEl = U.el("div", { class: "gs-pop" });
    pickerEl.style.position = "fixed";
    pickerEl.style.zIndex = "360";
    _renderPicker();
    document.body.append(pickerEl);
    const r = anchor.getBoundingClientRect();
    const ph = pickerEl.offsetHeight;
    const pw = pickerEl.offsetWidth;
    let top = r.top - 6 - ph;
    if (top < 8) top = Math.min(r.bottom + 6, window.innerHeight - ph - 8);
    const left = Math.min(r.left, window.innerWidth - pw - 8);
    pickerEl.style.top = `${Math.max(8, top)}px`;
    pickerEl.style.left = `${Math.max(8, left)}px`;
    document.addEventListener("mousedown", _pickerOutside, true);
  }

  function _renderPicker() {
    if (!pickerEl) return;
    U.clear(pickerEl);
    const unit = _currentUnit();
    // group by physical VARIABLE (wind speed, direction, water level, …);
    // the sources for that variable (wind.txt, raw waterinfo, ERA5, …) are
    // listed under it and can be multi-selected
    const groups = new Map();
    for (const [key, s] of sources) {
      const group = s.variable || s.group || "Other";
      if (!groups.has(group)) groups.set(group, []);
      groups.get(group).push([key, s]);
    }
    if (!groups.size) {
      pickerEl.append(U.el("div", { class: "muted", style: "padding:8px 10px" },
        "No series available yet."));
      return;
    }
    pickerEl.append(U.el("div", { class: "gs-hint" }, "Pick a variable, then its sources"));
    // keep "Output" groups last, everything else in insertion order
    const ordered = [...groups.entries()].sort((a, b) =>
      (a[0] === "Output" ? 1 : 0) - (b[0] === "Output" ? 1 : 0));
    for (const [groupName, items] of ordered) _pickerGroup(groupName, items, unit);

    if (selection.length) {
      const clear = U.el("button", { class: "ghost gs-clear" }, "clear selection");
      clear.addEventListener("click", () => {
        selection = [];
        _persistSelection();
        _rebuild();
        _renderPicker();
      });
      pickerEl.append(clear);
    }
  }

  function _pickerGroup(name, items, unit) {
    pickerEl.append(U.el("div", { class: "gs-group" }, name));
    for (const [key, s] of items) {
      const on = selection.includes(key);
      const conflict = !on && unit !== null && (s.unit || "") !== unit;
      const row = U.el("div", {
        class: `gs-row ${on ? "on" : ""} ${conflict ? "conflict" : ""}`,
        title: conflict
          ? `Unit [${s.unit || "–"}] differs from the plotted [${unit || "–"}] — click to plot this instead`
          : (on ? "Click to remove from the plot" : "Click to add to the plot"),
      },
        U.el("span", { class: "gs-dot", style: on
          ? `background:${TSPlot.COLORS[_selValid().indexOf(key) % TSPlot.COLORS.length]}`
          : "" }),
        U.el("span", { class: "gs-label" }, s.source || s.label),
        U.el("span", { class: "gs-unit" }, s.unit || ""));
      row.addEventListener("click", () => {
        toggle(key);
        _renderPicker();
      });
      pickerEl.append(row);
    }
  }

  /* ================= repeated-input shading ================= */

  function _drawRepeats(u, plotted) {
    const { min, max } = u.scales.x;
    const ctx = u.ctx;
    let firstRepeat = Infinity;
    for (const key of plotted) {
      const s = sources.get(key);
      if (Number.isFinite(s.repeatFrom)) firstRepeat = Math.min(firstRepeat, s.repeatFrom);
    }
    if (!Number.isFinite(firstRepeat) || firstRepeat >= max) return;

    // faint boundary marker at where AeoLiS starts repeating the file
    const x0 = u.valToPos(Math.max(firstRepeat, min), "x", true);
    ctx.save();
    if (firstRepeat >= min) {
      ctx.strokeStyle = "rgba(160, 130, 60, .5)";
      ctx.lineWidth = 1;
      ctx.setLineDash([2, 3]);
      ctx.beginPath();
      ctx.moveTo(x0, u.bbox.top);
      ctx.lineTo(x0, u.bbox.top + u.bbox.height);
      ctx.stroke();
    }
    ctx.restore();

    // the repeated series itself: the file's polyline tiled forward,
    // drawn thinner and dotted so it reads as "AeoLiS repeats this"
    plotted.forEach((key, i) => {
      const s = sources.get(key);
      if (!Number.isFinite(s.repeatFrom)) return;
      _drawRepeatLine(u, s, TSPlot.COLORS[i % TSPlot.COLORS.length]);
    });
  }

  /* Tile a source's own (original) polyline forward across the repeated
   * region [repeatFrom, xmax], thin + dotted. */
  function _drawRepeatLine(u, s, color) {
    const d = s.data;
    if (!d || !d[0] || d[0].length < 2) return;
    const t = d[0], v = d[1];
    const t0 = t[0], tEnd = t[t.length - 1];
    const period = tEnd - t0;
    if (!(period > 0)) return;
    const { min, max } = u.scales.x;
    const from = Math.max(s.repeatFrom, min);
    const ctx = u.ctx;
    const dpr = devicePixelRatio;
    ctx.save();
    ctx.beginPath();
    ctx.rect(u.bbox.left, u.bbox.top, u.bbox.width, u.bbox.height);
    ctx.clip();
    ctx.strokeStyle = color;
    ctx.globalAlpha = 0.55;
    ctx.lineWidth = 1 * dpr;
    ctx.setLineDash([1.5 * dpr, 2.5 * dpr]);
    for (let k = 1; ; k += 1) {
      const shift = k * period;
      if (t0 + shift > max) break;
      if (tEnd + shift < from) continue;
      ctx.beginPath();
      let started = false;
      for (let idx = 0; idx < t.length; idx += 1) {
        const tt = t[idx] + shift;
        const vv = v[idx];
        if (tt < from || tt > max) { started = false; continue; }
        if (vv === null || vv === undefined || Number.isNaN(vv)) { started = false; continue; }
        const px = u.valToPos(tt, "x", true);
        const py = u.valToPos(vv, s.scale || "y", true);
        if (!started) { ctx.moveTo(px, py); started = true; } else ctx.lineTo(px, py);
      }
      ctx.stroke();
    }
    ctx.restore();
  }

  /* ================= simulation start/end boundaries ================= */

  // Fun themes re-skin the series lines too. When the active theme has a
  // palette here, series cycle through it; otherwise the default TSPlot set.
  const SERIES_PALETTES = {
    darkside:  ["#ff2d2d", "#ff9500", "#ffe000", "#38d430", "#2f6bff", "#8a3ffb"],
    discovery: ["#e6eaf0", "#9aa2af", "#2ae0cf", "#c4ccd6", "#6fe07a", "#7d8593"],
    nevermind: ["#05314e", "#0a4a72", "#106a9c", "#2f9e63", "#4f9fd0", "#8bbf3c"],
    rumours:   ["#9a7b56", "#3b3327", "#a8584a", "#6e8a5a", "#cfa961", "#836542"],
    velvet:    ["#141414", "#f7c600", "#d0322f", "#3f9e57", "#c9a400", "#7a7565"],
  };
  function _seriesColor(i) {
    const p = SERIES_PALETTES[document.documentElement.dataset.theme];
    return p ? p[i % p.length] : TSPlot.COLORS[i % TSPlot.COLORS.length];
  }

  // read the theme accent once; invalidated on theme change (see init)
  let _accentCache = null;
  let _dangerCache = null;
  function _accent() {
    if (_accentCache === null) {
      _accentCache = (getComputedStyle(document.documentElement)
        .getPropertyValue("--accent").trim()) || "#0f766e";
    }
    return _accentCache;
  }
  function _danger() {
    if (_dangerCache === null) {
      _dangerCache = (getComputedStyle(document.documentElement)
        .getPropertyValue("--danger").trim()) || "#b4423b";
    }
    return _dangerCache;
  }
  let _mutedCache = null;
  function _muted() {
    if (_mutedCache === null) {
      _mutedCache = (getComputedStyle(document.documentElement)
        .getPropertyValue("--muted").trim()) || "#6b7686";
    }
    return _mutedCache;
  }

  /* Device-pixel span of the availability band, at the BOTTOM of the canvas
   * (below the x-axis labels). When hidden, the band collapses and the
   * grid/sim lines simply reach the plot baseline. */
  function _bandGeom(u) {
    const dpr = window.devicePixelRatio || 1;
    const H = u.ctx.canvas.height;
    const visible = _availVisible() && Boolean(availLast);
    if (!visible) {
      const y = u.bbox.top + u.bbox.height;
      return { top: y, bottom: y, visible: false, dpr };
    }
    return { top: H - BAND_PX * dpr + 2 * dpr, bottom: H - 3 * dpr, visible: true, dpr };
  }

  /* Vertical grid lines at the x-axis ticks. Drawn through the plot AND the
   * availability band (skipping the label strip between them) so the raster
   * reads as spanning the whole panel. */
  function _drawGridLines(u) {
    const scale = u.scales.x;
    if (!Number.isFinite(scale.min) || !Number.isFinite(scale.max)) return;
    const splits = TSPlot.timeSplits(u, 0, scale.min, scale.max) || [];
    const ctx = u.ctx;
    const plotTop = u.bbox.top;
    const plotBottom = u.bbox.top + u.bbox.height;
    const band = _bandGeom(u);
    ctx.save();
    ctx.strokeStyle = "rgba(128,128,128,.18)";
    ctx.lineWidth = 1;
    for (const t of splits) {
      if (t < scale.min || t > scale.max) continue;
      const x = Math.round(u.valToPos(t, "x", true)) + 0.5;
      ctx.beginPath();
      ctx.moveTo(x, plotTop);
      ctx.lineTo(x, plotBottom);
      ctx.stroke();
      if (band.visible) {
        ctx.beginPath();
        ctx.moveTo(x, band.top);
        ctx.lineTo(x, band.bottom);
        ctx.stroke();
      }
    }
    ctx.restore();
  }

  /* A slightly thicker baseline at the bottom of the plot so y=0 reads as
   * the floor of the graph (there are no protruding ticks below it). */
  function _drawBaseline(u) {
    const ctx = u.ctx;
    const y = u.bbox.top + u.bbox.height + 0.5;
    ctx.save();
    ctx.strokeStyle = _muted();
    ctx.globalAlpha = 0.55;
    ctx.lineWidth = 1.5 * (window.devicePixelRatio || 1);
    ctx.beginPath();
    ctx.moveTo(u.bbox.left, y);
    ctx.lineTo(u.bbox.left + u.bbox.width, y);
    ctx.stroke();
    ctx.restore();
  }

  function _roundRect(ctx, x, y, w, h, r) {
    const rr = Math.min(r, w / 2, h / 2);
    ctx.beginPath();
    ctx.moveTo(x + rr, y);
    ctx.arcTo(x + w, y, x + w, y + h, rr);
    ctx.arcTo(x + w, y + h, x, y + h, rr);
    ctx.arcTo(x, y + h, x, y, rr);
    ctx.arcTo(x, y, x + w, y, rr);
    ctx.closePath();
  }

  /* Data-availability lanes, drawn on the canvas in the reserved top band so
   * they share the chart's exact x-scale (perfect alignment, live pan/zoom).
   * Every lane uses the same accent blue; survey events are diamonds. */
  function _drawAvailBand(u) {
    if (!availLast || !_availVisible()) return;
    const { rows } = _availRows(availLast);
    if (!rows.length) return;
    const ctx = u.ctx;
    const dpr = window.devicePixelRatio || 1;
    const scale = u.scales.x;
    const left = u.bbox.left;
    const right = u.bbox.left + u.bbox.width;
    const band = _bandGeom(u);
    const top = band.top;
    const bottom = band.bottom;
    const laneH = (bottom - top) / rows.length;
    const accent = _accent();
    const xOf = (t) => u.valToPos(t, "x", true);
    ctx.save();
    ctx.font = `${10.5 * dpr}px ${(getComputedStyle(document.documentElement)
      .getPropertyValue("--font-ui") || "system-ui").trim()}`;
    ctx.textBaseline = "middle";
    rows.forEach((row, i) => {
      const cy = top + i * laneH;
      const mid = cy + laneH / 2;
      // lane label in the left gutter (where the y-axis sits below)
      ctx.fillStyle = _muted();
      ctx.textAlign = "left";
      ctx.fillText(row.label, 3 * dpr, mid, left - 7 * dpr);
      const barH = Math.max(3 * dpr, laneH * 0.56);
      const barTop = mid - barH / 2;
      for (const [a, b, cls] of row.spans || []) {
        if (b < scale.min || a > scale.max) continue;
        const xa = Math.max(xOf(Math.max(a, scale.min)), left);
        const xb = Math.min(xOf(Math.min(b, scale.max)), right);
        if (xb <= xa + 0.5) continue;
        ctx.fillStyle = accent;
        ctx.globalAlpha = cls === "repeat" ? 0.32 : (cls === "sim" ? 0.3 : 0.8);
        _roundRect(ctx, xa, barTop, xb - xa, barH, 3 * dpr);
        ctx.fill();
      }
      ctx.globalAlpha = 1;
      for (const m of row.marks || []) {
        if (m < scale.min || m > scale.max) continue;
        const xm = xOf(m);
        const r = Math.min(laneH * 0.34, 5 * dpr);
        ctx.fillStyle = accent;
        ctx.beginPath();
        ctx.moveTo(xm, mid - r);
        ctx.lineTo(xm + r, mid);
        ctx.lineTo(xm, mid + r);
        ctx.lineTo(xm - r, mid);
        ctx.closePath();
        ctx.fill();
      }
    });
    ctx.restore();
  }

  function _simBounds() {
    if (!availLast) return null;
    const { simT0, simT1 } = _availRows(availLast);
    if (!Number.isFinite(simT0) || !Number.isFinite(simT1)) return null;
    return [simT0, simT1];
  }

  /* Solid vertical lines at the simulation start & end, so the chart and
   * the availability strip below are tied together by the same markers. */
  function _drawSimLines(u) {
    const bounds = _simBounds();
    if (!bounds) return;
    const scale = u.scales.x;
    const ctx = u.ctx;
    const band = _bandGeom(u);
    const top = u.bbox.top;
    const bottom = band.visible ? band.bottom : (u.bbox.top + u.bbox.height);
    ctx.save();
    ctx.globalAlpha = 0.55;
    ctx.lineWidth = 1.5;
    bounds.forEach((t, i) => {
      if (!Number.isFinite(t) || t < scale.min || t > scale.max) return;
      // start = accent, end = danger; one continuous line spanning the plot
      // and the availability band below (ties the two together, no gap)
      const color = i === 0 ? _accent() : _danger();
      ctx.strokeStyle = color;
      ctx.globalAlpha = 0.6;
      const x = u.valToPos(t, "x", true);
      ctx.beginPath();
      ctx.moveTo(x, top);
      ctx.lineTo(x, bottom);
      ctx.stroke();
      // a solid "pin" symbol at the very top so start/finish pop out
      ctx.globalAlpha = 1;
      ctx.fillStyle = color;
      const r = 5 * (window.devicePixelRatio || 1);
      ctx.beginPath();
      ctx.moveTo(x - r, top);
      ctx.lineTo(x + r, top);
      ctx.lineTo(x, top + r * 1.5);
      ctx.closePath();
      ctx.fill();
    });
    ctx.restore();
  }

  /* ================= playbar time cursor ================= */

  function _drawTimeCursor(u) {
    const t = App.state.clock.t;
    if (!Number.isFinite(t)) return;
    const scale = u.scales.x;
    if (t < scale.min || t > scale.max) return;
    const x = u.valToPos(t, "x", true);
    const ctx = u.ctx;
    const band = _bandGeom(u);
    ctx.save();
    ctx.strokeStyle = _accent();
    ctx.lineWidth = 1.5;
    ctx.setLineDash([5, 4]);
    ctx.beginPath();
    ctx.moveTo(x, u.bbox.top);
    ctx.lineTo(x, band.visible ? band.bottom : (u.bbox.top + u.bbox.height));
    ctx.stroke();
    ctx.restore();
  }

  function redrawCursors() {
    if (ts) ts.redraw();   // drawExtra repaints the band + cursor on the canvas
  }

  function resizeAll() {
    if (!ts) { if (card) _rebuild(); return; }
    const { width, height } = _plotSize();
    if (width > 50) ts.setSize({ width, height });
  }

  /* =================================================================
   * Data availability (drawn inside the chart canvas, same time axis)
   * ================================================================= */

  const refreshAvailability = U.debounce(async () => {
    if (!App.state.project) return;
    let conditions = null, output = null, domain = null;
    try { conditions = await Api.get("/api/conditions"); } catch { /* no config */ }
    try { output = await Api.get("/api/output/meta"); } catch { /* none */ }
    try { domain = await Api.get("/api/domain"); } catch { /* none */ }
    if (!conditions) return;
    const hadData = Boolean(availLast);
    availLast = { conditions, output, domain };
    // widen the chart range FIRST so the band renders against the same span
    const fw = _fullWindow();
    if (ts) ts.setFullRange(fw);
    // keep the playbar slider usable off the graph's own time span (the
    // simulation range at minimum), so it scrubs/animates the time cursor
    // even before a run produces output
    if (typeof Playbar !== "undefined") {
      if (fw) Playbar.setSource("graph", fw[0], fw[1]);
      else Playbar.removeSource("graph");
    }
    // when availability first appears the reserved band height changes, so
    // the plot must be rebuilt; afterwards a redraw repaints the band
    if (ts && hadData) ts.redraw();
    else _rebuild();
  }, 400);

  function _availRows({ conditions, output, domain }) {
    const rows = [];
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

    for (const [kind, title] of [["wind", "Wind"], ["tide", "Water lvl"], ["wave", "Waves"]]) {
      const info = conditions.kinds[kind];
      if (!info || !info.series) continue;
      const s = info.series;
      // split the coverage bar at real data gaps (large dt / missing values)
      const spans = _dataSegments(s).map(([a, b]) => [a, b, "data"]);
      if (s.t1_epoch < simT1 - 1) spans.push([s.t1_epoch, simT1, "repeat"]);
      rows.push({ label: title, spans });
    }

    if (output && output.exists && output.times_epoch && output.times_epoch.length) {
      rows.push({
        label: "Output",
        spans: [[output.times_epoch[0], output.times_epoch[output.times_epoch.length - 1], "output"]],
      });
    }

    if (domain && domain.entries && domain.entries.length) {
      const marks = [];
      for (const entry of domain.entries) {
        if (entry.date) marks.push(Date.parse(entry.date) / 1000);
        else if (entry.year) marks.push(Date.UTC(entry.year, 6, 1) / 1000);
      }
      if (marks.length) rows.push({ label: "Surveys", marks });
    }
    return { rows, simT0, simT1 };
  }

  /* Contiguous data-coverage segments of a (decimated) series: break on
   * a time step much larger than the median cadence or on missing values,
   * so genuine gaps show as gaps instead of one solid bar. */
  function _dataSegments(s) {
    const t = s.t_epoch;
    if (!t || t.length < 2) return [[s.t0_epoch, s.t1_epoch]];
    const col = (s.columns && s.columns[0]) || null;
    const dts = [];
    for (let i = 1; i < t.length; i += 1) dts.push(t[i] - t[i - 1]);
    const sorted = [...dts].sort((a, b) => a - b);
    const median = sorted[Math.floor(sorted.length / 2)] || 1;
    const cut = Math.max(median * 6, 1);
    const present = (i) => (col ? (col[i] !== null && col[i] !== undefined) : true);
    const segs = [];
    let start = null, prev = null;
    for (let i = 0; i < t.length; i += 1) {
      if (!present(i)) { if (start !== null) { segs.push([start, prev]); start = null; } continue; }
      if (start === null) start = t[i];
      else if (t[i] - prev > cut) { segs.push([start, prev]); start = t[i]; }
      prev = t[i];
    }
    if (start !== null) segs.push([start, prev]);
    return segs.length ? segs : [[s.t0_epoch, s.t1_epoch]];
  }

  function _availFullRange(payloads) {
    const { rows, simT0, simT1 } = _availRows(payloads);
    let t0 = simT0, t1 = simT1;
    for (const row of rows) {
      for (const [a, b] of row.spans || []) { t0 = Math.min(t0, a); t1 = Math.max(t1, b); }
      for (const m of row.marks || []) { t0 = Math.min(t0, m); t1 = Math.max(t1, m); }
    }
    if (!Number.isFinite(t0) || t1 <= t0) return null;
    return [t0, t1];
  }

  /* ================= init ================= */

  function init() {
    App.on("clock-tick", redrawCursors);
    App.on("theme", () => {
      _accentCache = null; _dangerCache = null; _mutedCache = null;
      _rebuild(); redrawCursors();
    });
    window.addEventListener("resize", U.debounce(resizeAll, 150));
    App.on("project", refreshAvailability);
    App.on("config-changed", refreshAvailability);
    App.on("run-finished", refreshAvailability);
    window.addEventListener("keydown", (ev) => {
      if (ev.key === "Escape") _closePicker();
    });
  }

  return { init, registerSource, unregisterSource, clearAll, select, toggle,
    redrawCursors, resizeAll, refreshAvailability, setTimeWindow,
    timeWindow: () => timeWindow,
    timeAxis: () => TSPlot.timeAxis() };
})();
