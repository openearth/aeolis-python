/* Windrose popup for any magnitude + direction pair (wind, waves — raw
 * downloads, wind.txt, etc.).
 *
 * Directions follow the meteorological convention (bearing the flow comes
 * FROM: 0 = N, 90 = E). Petals point toward that bearing and stack speed
 * bins from the centre outward. All controls (sectors, bin width, colour
 * scheme, calm threshold) recompute and redraw live.
 */
"use strict";

const Windrose = (() => {

  const SCHEMES = {
    viridis: ["#440154", "#414487", "#2a788e", "#22a884", "#7ad151", "#fde725"],
    blues:   ["#eff3ff", "#c6dbef", "#9ecae1", "#6baed6", "#3182bd", "#08519c"],
    warm:    ["#ffffb2", "#fed976", "#feb24c", "#fd8d3c", "#f03b20", "#bd0026"],
    cool:    ["#f7fcf0", "#bae4bc", "#7bccc4", "#43a2ca", "#2b8cbe", "#0868ac"],
  };

  function _hex(c) {
    return [parseInt(c.slice(1, 3), 16), parseInt(c.slice(3, 5), 16), parseInt(c.slice(5, 7), 16)];
  }
  /* sample a ramp at t in [0,1] */
  function _ramp(stops, t) {
    t = Math.max(0, Math.min(1, t));
    const seg = t * (stops.length - 1);
    const i = Math.min(stops.length - 2, Math.floor(seg));
    const f = seg - i;
    const a = _hex(stops[i]), b = _hex(stops[i + 1]);
    const mix = a.map((v, k) => Math.round(v + (b[k] - v) * f));
    return `rgb(${mix[0]},${mix[1]},${mix[2]})`;
  }

  /* bin the samples into [sector][speedbin] frequencies (% of all valid) */
  function _compute(mag, dir, opts) {
    const { sectors, binWidth, calm } = opts;
    const sw = 360 / sectors;
    let total = 0, calmCount = 0, maxMag = 0;
    const pairs = [];
    for (let i = 0; i < mag.length; i += 1) {
      const m = mag[i], d = dir[i];
      if (m === null || d === null || !Number.isFinite(m) || !Number.isFinite(d)) continue;
      total += 1;
      if (m < calm) { calmCount += 1; continue; }
      if (m > maxMag) maxMag = m;
      pairs.push([m, ((d % 360) + 360) % 360]);
    }
    const nBins = Math.max(1, Math.min(7, Math.ceil((maxMag || 1) / binWidth)));
    const counts = Array.from({ length: sectors }, () => new Array(nBins).fill(0));
    for (const [m, d] of pairs) {
      const s = Math.round(d / sw) % sectors;
      const bin = Math.min(nBins - 1, Math.floor(m / binWidth));
      counts[s][bin] += 1;
    }
    // to percentages of total
    const freq = counts.map((row) => row.map((c) => total ? 100 * c / total : 0));
    const sectorTotals = freq.map((row) => row.reduce((a, b) => a + b, 0));
    return {
      freq, nBins, sectors, sw,
      calmPct: total ? 100 * calmCount / total : 0,
      total,
      maxSector: Math.max(0.001, ...sectorTotals),
      binWidth,
    };
  }

  const SIZE = 320, CX = SIZE / 2, CY = SIZE / 2, R0 = 14, RMAX = 130;
  const DIRLABELS = [["N", 0], ["E", 90], ["S", 180], ["W", 270]];

  function _pt(r, bearingDeg) {
    const a = bearingDeg * Math.PI / 180;
    return [CX + r * Math.sin(a), CY - r * Math.cos(a)];
  }
  function _wedge(r0, r1, a0, a1) {
    const [x0, y0] = _pt(r0, a0), [x1, y1] = _pt(r1, a0);
    const [x2, y2] = _pt(r1, a1), [x3, y3] = _pt(r0, a1);
    return `M${x0} ${y0} L${x1} ${y1} A${r1} ${r1} 0 0 1 ${x2} ${y2} `
      + `L${x3} ${y3} A${r0} ${r0} 0 0 0 ${x0} ${y0} Z`;
  }
  const SVGNS = "http://www.w3.org/2000/svg";
  function _svgEl(tag, attrs) {
    const n = document.createElementNS(SVGNS, tag);
    for (const [k, v] of Object.entries(attrs)) n.setAttribute(k, v);
    return n;
  }

  function _draw(host, data, scheme) {
    U.clear(host);
    const stops = SCHEMES[scheme] || SCHEMES.viridis;
    const svg = _svgEl("svg", { viewBox: `0 0 ${SIZE} ${SIZE}`, width: "100%", height: "100%" });
    // radial grid: 4 rings scaled to the busiest sector, "nice" %
    const ringMax = data.maxSector;
    const step = _niceStep(ringMax / 4);
    const muted = getComputedStyle(document.documentElement).getPropertyValue("--border") || "#ddd";
    const textC = getComputedStyle(document.documentElement).getPropertyValue("--muted") || "#888";
    for (let v = step; v <= ringMax * 1.001; v += step) {
      const r = R0 + (RMAX - R0) * (v / ringMax);
      svg.append(_svgEl("circle", { cx: CX, cy: CY, r, fill: "none", stroke: muted, "stroke-width": 1 }));
      const lbl = _svgEl("text", { x: CX + 3, y: CY - r + 3, "font-size": 9, fill: textC });
      lbl.textContent = `${v.toFixed(step < 1 ? 1 : 0)}%`;
      svg.append(lbl);
    }
    // spokes + compass labels
    for (const [name, bearing] of DIRLABELS) {
      const [lx, ly] = _pt(RMAX + 14, bearing);
      const t = _svgEl("text", {
        x: lx, y: ly, "font-size": 12, "font-weight": 700, fill: textC,
        "text-anchor": "middle", "dominant-baseline": "middle",
      });
      t.textContent = name;
      svg.append(t);
    }
    // petals: stacked speed bins per sector
    const half = data.sw / 2 * 0.86;   // small gap between petals
    for (let s = 0; s < data.sectors; s += 1) {
      const center = s * data.sw;
      let acc = 0;
      for (let bin = 0; bin < data.nBins; bin += 1) {
        const f = data.freq[s][bin];
        if (f <= 0) { continue; }
        const r0 = R0 + (RMAX - R0) * (acc / ringMax);
        const r1 = R0 + (RMAX - R0) * ((acc + f) / ringMax);
        acc += f;
        const color = _ramp(stops, data.nBins === 1 ? 0.5 : bin / (data.nBins - 1));
        svg.append(_svgEl("path", {
          d: _wedge(r0, r1, center - half, center + half),
          fill: color, stroke: "rgba(0,0,0,.12)", "stroke-width": 0.5,
        }));
      }
    }
    // calm centre
    svg.append(_svgEl("circle", { cx: CX, cy: CY, r: R0, fill: "var(--panel-2)", stroke: muted }));
    const calm = _svgEl("text", {
      x: CX, y: CY, "font-size": 8, fill: textC,
      "text-anchor": "middle", "dominant-baseline": "middle",
    });
    calm.textContent = `${data.calmPct.toFixed(0)}%`;
    svg.append(calm);
    host.append(svg);
  }

  function _niceStep(x) {
    if (x <= 0) return 1;
    const pow = Math.pow(10, Math.floor(Math.log10(x)));
    const n = x / pow;
    const nice = n < 1.5 ? 1 : n < 3 ? 2 : n < 7 ? 5 : 10;
    return nice * pow;
  }

  function _legend(host, data, scheme, unit) {
    U.clear(host);
    const stops = SCHEMES[scheme] || SCHEMES.viridis;
    host.append(U.el("div", { class: "wr-legend-title" }, `Speed [${unit || "–"}]`));
    for (let bin = 0; bin < data.nBins; bin += 1) {
      const lo = (bin * data.binWidth);
      const hi = (bin + 1) * data.binWidth;
      const label = bin === data.nBins - 1 ? `≥ ${U.fmtNum(lo, 3)}` : `${U.fmtNum(lo, 3)}–${U.fmtNum(hi, 3)}`;
      const color = _ramp(stops, data.nBins === 1 ? 0.5 : bin / (data.nBins - 1));
      host.append(U.el("div", { class: "wr-legend-row" },
        U.el("span", { class: "wr-swatch", style: `background:${color}` }),
        U.el("span", {}, label)));
    }
    host.append(U.el("div", { class: "wr-legend-calm" },
      `calm < ${U.fmtNum(data.calm ?? 0, 3)} ${unit || ""} · n = ${data.total}`));
  }

  /* a sensible default speed-bin width from the data spread */
  function _defaultBinWidth(mag) {
    let mx = 0;
    for (const m of mag) { if (Number.isFinite(m) && m > mx) mx = m; }
    if (mx <= 0) return 1;
    return _niceStep(mx / 5);
  }

  /* epoch seconds → a datetime-local value (UTC, matched by the "Z" parse
   * below so the period filter lines up with the UTC data epochs). */
  function _isoLocal(epoch) {
    return Number.isFinite(epoch) ? new Date(epoch * 1000).toISOString().slice(0, 16) : "";
  }
  function _parseLocal(v) {
    return v ? Date.parse(`${v}:00Z`) / 1000 : null;
  }

  /* The one windrose window: it shows the rose immediately for the first
   * source over its full period, and every control (source, period, sectors,
   * bin, calm, colours) recomputes it live.
   *
   * opts.sources = [{label, magName, magUnit, t0_epoch?, t1_epoch?,
   *                  load: async () => ({magnitude, direction, t_epoch})}]
   * or the legacy direct form { title, magnitude, direction, magName, magUnit }. */
  async function open(opts = {}) {
    let sources = opts.sources;
    if (!sources) {
      // legacy direct-arrays form → wrap as a single source
      sources = [{
        label: opts.title || opts.magName || "series",
        magName: opts.magName, magUnit: opts.magUnit,
        load: async () => ({
          magnitude: opts.magnitude || [], direction: opts.direction || [],
          t_epoch: opts.t_epoch || [],
        }),
      }];
    }

    const popup = Popup.open({ title: "Windrose", width: 840 });
    if (!sources.length) {
      popup.body.append(U.el("div", { class: "muted" },
        "No magnitude + direction data available yet. Download or fill a wind or "
        + "wave series in the Conditions tab first."));
      return popup;
    }

    const state = { sectors: 16, binWidth: 1, calm: 0.5, scheme: "viridis" };
    let cur = null;   // { mag, dir, t, unit, magName } of the loaded source

    const roseHost = U.el("div", { class: "wr-rose" });
    const legendHost = U.el("div", { class: "wr-legend" });
    const infoLine = U.el("div", { class: "muted", style: "font-size:11.5px;margin-top:6px" });

    // --- source + period ---
    const srcSel = U.el("select", {}, ...sources.map((s, i) =>
      U.el("option", { value: String(i) }, s.label)));
    const fromIn = U.el("input", { type: "datetime-local" });
    const toIn = U.el("input", { type: "datetime-local" });

    // --- rose controls ---
    const sectorsSel = U.el("select", {},
      ...[8, 16, 32].map((n) => U.el("option", { value: n, selected: n === state.sectors ? "" : null }, `${n} sectors`)));
    sectorsSel.value = String(state.sectors);
    sectorsSel.addEventListener("change", () => { state.sectors = Number(sectorsSel.value); redraw(); });

    const binInput = U.numField(state.binWidth, (v) => { if (v > 0) { state.binWidth = v; redraw(); } },
      { style: "width:70px" });
    const calmInput = U.numField(state.calm, (v) => { if (v >= 0) { state.calm = v; redraw(); } },
      { style: "width:70px" });
    const binLabel = U.el("label", {}, "Speed bin");
    const calmLabel = U.el("label", {}, "Calm below");

    const schemeSel = U.el("select", {},
      ...Object.keys(SCHEMES).map((s) => U.el("option", { value: s }, s)));
    schemeSel.value = state.scheme;
    schemeSel.addEventListener("change", () => { state.scheme = schemeSel.value; redraw(); });

    fromIn.addEventListener("change", redraw);
    toIn.addEventListener("change", redraw);

    function redraw() {
      if (!cur) return;
      const t0 = _parseLocal(fromIn.value);
      const t1 = _parseLocal(toIn.value);
      let mag = cur.mag, dir = cur.dir;
      if ((t0 !== null || t1 !== null) && cur.t.length === mag.length && cur.t.length) {
        const m2 = [], d2 = [];
        for (let i = 0; i < cur.t.length; i += 1) {
          if (t0 !== null && cur.t[i] < t0) continue;
          if (t1 !== null && cur.t[i] > t1) continue;
          m2.push(mag[i]); d2.push(dir[i]);
        }
        mag = m2; dir = d2;
      }
      const data = _compute(mag, dir, state);
      data.calm = state.calm;
      _draw(roseHost, data, state.scheme);
      _legend(legendHost, data, state.scheme, cur.unit);
      const u = cur.unit ? ` ${cur.unit}` : "";
      binLabel.textContent = `Speed bin [${cur.unit || "–"}]`;
      calmLabel.textContent = `Calm below [${cur.unit || "–"}]`;
      infoLine.textContent = `${cur.magName || "magnitude"}${u} · ${data.total} samples in view`;
    }

    async function loadSource(i, fillPeriod) {
      const src = sources[i] || sources[0];
      U.clear(roseHost);
      roseHost.append(U.el("div", { class: "muted", style: "padding:20px" }, "loading…"));
      let d;
      try { d = await src.load(); }
      catch (err) { U.toast(err.message, "error"); return; }
      const t = d.t_epoch || [];
      cur = { mag: d.magnitude || [], dir: d.direction || [], t,
        unit: src.magUnit || "", magName: src.magName };
      state.binWidth = _defaultBinWidth(cur.mag);
      binInput.value = String(state.binWidth);
      if (fillPeriod) {
        const lo = Number.isFinite(src.t0_epoch) ? src.t0_epoch : (t.length ? t[0] : null);
        const hi = Number.isFinite(src.t1_epoch) ? src.t1_epoch : (t.length ? t[t.length - 1] : null);
        fromIn.value = _isoLocal(lo);
        toIn.value = _isoLocal(hi);
      }
      redraw();
    }
    srcSel.addEventListener("change", () => loadSource(Number(srcSel.value), true));

    const controls = U.el("div", { class: "wr-controls" },
      U.el("div", { class: "form-row" }, U.el("label", {}, "Source"), srcSel),
      U.el("div", { class: "form-row" }, U.el("label", {}, "From"), fromIn),
      U.el("div", { class: "form-row" }, U.el("label", {}, "To"), toIn),
      U.el("div", { class: "form-row" }, U.el("label", {}, "Sectors"), sectorsSel),
      U.el("div", { class: "form-row" }, binLabel, binInput),
      U.el("div", { class: "form-row" }, calmLabel, calmInput),
      U.el("div", { class: "form-row" }, U.el("label", {}, "Colours"), schemeSel),
      infoLine);

    // three aligned columns: legend (left) · rose (center) · settings (right)
    popup.body.append(U.el("div", { class: "wr-wrap" },
      legendHost,
      roseHost,
      controls));
    await loadSource(0, true);   // draw immediately for the first source
    return popup;
  }

  // back-compat alias: the topbar calls openChooser({sources})
  const openChooser = open;

  return { open, openChooser };
})();
