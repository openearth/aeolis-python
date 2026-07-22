/* CondHud: a small draggable window floating over the map showing the
 * boundary conditions at the current playbar time - wind speed +
 * direction arrow, wave height/period, and water level - so the user
 * can quickly verify the input conditions are interpreted correctly.
 *
 * Toggled by a small map button; visibility persists per project. */
"use strict";

const CondHud = (() => {

  let panel = null;
  let toggleBtn = null;
  let series = null;    // {wind: {t, cols}, tide: ..., wave: ...}
  let els = {};

  function init() {
    App.on("project", _reload);
    App.on("config-changed", U.debounce(_maybeReload, 600));
    App.on("clock-tick", _update);
    _buildToggle();
  }

  function _maybeReload(key) {
    if (!key || /wind_file|tide_file|wave_file|refdate/.test(String(key))) _reload();
  }

  async function _reload() {
    if (!App.state.project) return;
    try {
      const overview = await Api.get("/api/conditions");
      series = {};
      for (const [kind, info] of Object.entries(overview.kinds)) {
        if (info.series) {
          series[kind] = {
            t: info.series.t_epoch,
            cols: info.series.columns,
            labels: info.labels,
          };
        }
      }
    } catch { series = null; }
    _sync();
    _update();
  }

  /* ---- UI ---- */

  function _buildToggle() {
    toggleBtn = U.el("button", {
      id: "condhud-toggle",
      title: "Show/hide input conditions at the current time",
    }, _windIcon());
    toggleBtn.addEventListener("click", () => {
      App.state.ui.condHud = !App.state.ui.condHud;
      App.touchUi();
      _sync();
      _update();
    });
    document.getElementById("map-wrap").append(toggleBtn);
  }

  function _windIcon() {
    const svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");
    svg.setAttribute("viewBox", "0 0 24 24");
    svg.setAttribute("width", "16");
    svg.setAttribute("height", "16");
    const path = document.createElementNS("http://www.w3.org/2000/svg", "path");
    path.setAttribute("fill", "currentColor");
    path.setAttribute("d", "M4 6h9a3 3 0 1 0-2.8-4l1.9.7A1 1 0 1 1 13 4H4v2zm0 5h13a3 3 0 1 1-2.8 4l1.9-.7A1 1 0 1 0 17 13H4v-2zm0 7h7a2.5 2.5 0 1 1-2.4 3.2l1.9-.6a.6.6 0 1 0 .5-.6H4v-2z");
    svg.append(path);
    return svg;
  }

  function _build() {
    if (panel) return;
    els = {};
    const grid = U.el("div", { class: "condhud-grid" });
    for (const [key, label] of [["wind", "Wind"], ["wave", "Waves"], ["tide", "Level"]]) {
      const value = U.el("span", { class: "ch-value" }, "–");
      const sub = U.el("span", { class: "ch-sub" }, "");
      els[key] = { value, sub };
      grid.append(
        U.el("span", { class: "ch-label" }, label),
        U.el("span", { class: "ch-cell" }, value, sub));
    }

    // big wind compass: N mark + rotating arrow, spans all rows so the
    // wind direction is readable against the map at a glance
    const arrow = document.createElementNS("http://www.w3.org/2000/svg", "svg");
    arrow.setAttribute("viewBox", "0 0 24 24");
    arrow.setAttribute("class", "ch-arrow");
    const path = document.createElementNS("http://www.w3.org/2000/svg", "path");
    path.setAttribute("fill", "currentColor");
    path.setAttribute("d", "M12 2.5 17 15l-5-2.6L7 15 12 2.5zM12 17.6a1.6 1.6 0 1 1 0 3.2 1.6 1.6 0 0 1 0-3.2z");
    arrow.append(path);
    const compass = U.el("div", { class: "ch-compass", title: "Direction the wind blows towards (map-north up)" },
      U.el("span", { class: "ch-n" }, "N"), arrow);
    els.arrow = arrow;
    els.compass = compass;
    grid.append(compass);

    els.time = U.el("span", { class: "condhud-time" }, "–");
    const closeBtn = U.el("button", { class: "ghost popup-close", title: "Hide" }, "✕");
    closeBtn.addEventListener("click", () => {
      App.state.ui.condHud = false;
      App.touchUi();
      _sync();
    });
    const head = U.el("div", { class: "condhud-head" },
      U.el("span", { class: "condhud-title" }, "Conditions"), els.time, closeBtn);
    panel = U.el("div", { id: "condhud" }, head, grid);
    document.getElementById("map-wrap").append(panel);
    _makeDraggable(head);
  }

  function _makeDraggable(handle) {
    handle.addEventListener("mousedown", (ev) => {
      if (ev.target.tagName === "BUTTON") return;
      ev.preventDefault();
      const wrap = document.getElementById("map-wrap").getBoundingClientRect();
      const rect = panel.getBoundingClientRect();
      const dx = ev.clientX - rect.left, dy = ev.clientY - rect.top;
      const move = (mv) => {
        const x = U.clamp(mv.clientX - wrap.left - dx, 0, wrap.width - rect.width);
        const y = U.clamp(mv.clientY - wrap.top - dy, 0, wrap.height - rect.height);
        panel.style.left = `${x}px`;
        panel.style.top = `${y}px`;
        // the panel is CSS-anchored bottom-right; releasing BOTH the right
        // and bottom anchors is required or the fixed bottom edge stretches
        // the panel while the top follows the cursor
        panel.style.right = "auto";
        panel.style.bottom = "auto";
      };
      const up = () => {
        window.removeEventListener("mousemove", move);
        window.removeEventListener("mouseup", up);
      };
      window.addEventListener("mousemove", move);
      window.addEventListener("mouseup", up);
    });
  }

  function _sync() {
    const on = Boolean(App.state.ui.condHud);
    if (on) _build();
    if (panel) panel.style.display = on ? "" : "none";
    if (toggleBtn) toggleBtn.classList.toggle("active", on);
  }

  /* ---- values at the current clock time ---- */

  function _at(kind, col) {
    const s = series && series[kind];
    let t = App.state.clock.t;
    if (!s || !Number.isFinite(t) || !s.t.length) return null;
    const ts = s.t;
    // outside the series: wrap cyclically, exactly like AeoLiS repeats
    // a short boundary-condition file (and like the graphs display it)
    const span = ts[ts.length - 1] - ts[0];
    if ((t < ts[0] || t > ts[ts.length - 1]) && span > 0) {
      t = ts[0] + (((t - ts[0]) % span) + span) % span;
    }
    if (t <= ts[0]) return s.cols[col] ? s.cols[col][0] : null;
    if (t >= ts[ts.length - 1]) {
      return s.cols[col] ? s.cols[col][ts.length - 1] : null;
    }
    let lo = 0, hi = ts.length - 1;
    while (hi - lo > 1) {
      const mid = (lo + hi) >> 1;
      if (ts[mid] <= t) lo = mid; else hi = mid;
    }
    const c = s.cols[col];
    if (!c) return null;
    const f = (t - ts[lo]) / Math.max(ts[hi] - ts[lo], 1e-9);
    const a = c[lo], b = c[hi];
    if (a === null || b === null) return a ?? b;
    return a + (b - a) * f;
  }

  function _update() {
    if (!panel || !App.state.ui.condHud) return;
    const t = App.state.clock.t;
    els.time.textContent = Number.isFinite(t) ? U.fmtDate(t) : "–";

    // fixed decimals per quantity so the panel never jumps around
    // wind: speed + direction (nautical: coming FROM, CW from north;
    // the compass arrow points where the wind blows TOWARDS)
    const speed = _at("wind", 0);
    const dir = _at("wind", 1);
    if (speed !== null) {
      els.wind.value.textContent = `${speed.toFixed(1)} m/s`;
      els.wind.sub.textContent = dir !== null ? `from ${String(Math.round(dir)).padStart(3, "0")}°` : "";
    } else {
      els.wind.value.textContent = "–";
      els.wind.sub.textContent = "";
    }
    if (dir !== null) {
      const convention = (App.state.config || {}).wind_convention || "nautical";
      // arrow glyph points up (north); rotate to the blow-towards
      // compass angle (map-north up)
      const towards = convention === "cartesian" ? 90 - dir : dir + 180;
      els.arrow.style.transform = `rotate(${Math.round(towards)}deg)`;
      els.compass.classList.remove("idle");
    } else {
      els.compass.classList.add("idle");
    }

    // waves: Hs + Tp (wave.txt has no direction column)
    const hs = _at("wave", 0);
    const tp = _at("wave", 1);
    if (hs !== null) {
      els.wave.value.textContent = `${hs.toFixed(2)} m`;
      els.wave.sub.textContent = tp !== null ? `Tp ${tp.toFixed(1)} s` : "";
    } else {
      els.wave.value.textContent = "–";
      els.wave.sub.textContent = "";
    }

    // water level
    const eta = _at("tide", 0);
    els.tide.value.textContent = eta !== null ? `${eta.toFixed(2)} m` : "–";
    els.tide.sub.textContent = "";
  }

  return { init };
})();
