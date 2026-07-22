/* Playbar: the shared time slider + play/step controls at the bottom.
 *
 * The clock is absolute model time in seconds since the config refdate
 * (epoch seconds). Sources (boundary-condition series, netCDF output)
 * register their time ranges; the slider spans the union. Every change
 * emits "clock-tick" which the map field layers and graphs listen to.
 */
"use strict";

const Playbar = (() => {

  const sources = new Map();   // id -> {t0, t1}
  let raf = null;
  let lastWall = null;

  const els = {};

  function _clock() { return App.state.clock; }

  /* ---- source registration ---- */

  function setSource(id, t0, t1) {
    sources.set(id, { t0, t1 });
    _recomputeRange();
  }

  function removeSource(id) {
    sources.delete(id);
    _recomputeRange();
  }

  function clearSources() {
    sources.clear();
    _recomputeRange();
  }

  function _recomputeRange() {
    const clock = _clock();
    if (!sources.size) {
      clock.t0 = clock.t1 = clock.t = null;
      _updateUI();
      return;
    }
    clock.t0 = Math.min(...Array.from(sources.values(), (s) => s.t0));
    clock.t1 = Math.max(...Array.from(sources.values(), (s) => s.t1));
    if (!Number.isFinite(clock.t) || clock.t < clock.t0 || clock.t > clock.t1) {
      clock.t = clock.t0;
    }
    _updateUI();
    App.emit("clock-tick", clock.t);
  }

  /* ---- time control ---- */

  function setTime(t, fromSlider = false) {
    const clock = _clock();
    if (!Number.isFinite(clock.t0)) return;
    clock.t = U.clamp(t, clock.t0, clock.t1);
    _updateUI(fromSlider);
    App.emit("clock-tick", clock.t);
  }

  function step(direction) {
    const clock = _clock();
    if (!Number.isFinite(clock.t)) return;
    const stepSize = App.state.outputTimes || 3600;
    setTime(clock.t + direction * stepSize);
  }

  function togglePlay(force = null) {
    const clock = _clock();
    clock.playing = force !== null ? force : !clock.playing;
    if (clock.playing && !Number.isFinite(clock.t)) clock.playing = false;
    els.play.textContent = clock.playing ? "⏸" : "▶";
    if (clock.playing) {
      lastWall = performance.now();
      raf = requestAnimationFrame(_tick);
    } else if (raf) {
      cancelAnimationFrame(raf);
      raf = null;
    }
  }

  function _tick(now) {
    const clock = _clock();
    if (!clock.playing) return;
    const dtWall = (now - lastWall) / 1000;
    lastWall = now;
    let t = clock.t + clock.speed * dtWall;
    if (t >= clock.t1) {
      t = clock.t0;   // loop
    }
    setTime(t);
    raf = requestAnimationFrame(_tick);
  }

  /* speed slider is logarithmic: 10^v model-seconds per wall-second */
  function _applySpeed() {
    const clock = _clock();
    clock.speed = Math.pow(10, parseFloat(els.speed.value));
  }

  /* ---- UI ---- */

  function _updateUI(fromSlider = false) {
    const clock = _clock();
    const hasRange = Number.isFinite(clock.t0) && clock.t1 > clock.t0;
    els.timeline.disabled = !hasRange;
    els.play.disabled = !hasRange;
    if (!hasRange) {
      els.label.textContent = "–";
      return;
    }
    if (!fromSlider) {
      const frac = (clock.t - clock.t0) / (clock.t1 - clock.t0);
      els.timeline.value = Math.round(frac * 1000);
    }
    els.label.textContent = U.fmtDate(clock.t);
  }

  /* Highlight the graphs' current zoom window on the timeline. */
  function setViewWindow(window) {
    const clock = _clock();
    let band = document.getElementById("timeline-window");
    if (!band) {
      const wrap = U.el("span", { id: "timeline-wrap" });
      els.timeline.replaceWith(wrap);
      wrap.append(U.el("span", { id: "timeline-window" }), els.timeline);
      band = document.getElementById("timeline-window");
    }
    if (!window || !Number.isFinite(clock.t0) || clock.t1 <= clock.t0) {
      band.style.display = "none";
      return;
    }
    const frac = (t) => U.clamp((t - clock.t0) / (clock.t1 - clock.t0), 0, 1);
    band.style.display = "";
    band.style.left = `${(frac(window[0]) * 100).toFixed(2)}%`;
    band.style.width = `${((frac(window[1]) - frac(window[0])) * 100).toFixed(2)}%`;
  }

  function init() {
    els.play = document.getElementById("btn-play");
    els.timeline = document.getElementById("timeline");
    els.label = document.getElementById("time-label");
    els.speed = document.getElementById("rng-speed");

    els.play.addEventListener("click", () => togglePlay());
    els.timeline.addEventListener("input", () => {
      const clock = _clock();
      if (!Number.isFinite(clock.t0)) return;
      const frac = parseInt(els.timeline.value, 10) / 1000;
      setTime(clock.t0 + frac * (clock.t1 - clock.t0), true);
    });
    els.speed.addEventListener("input", _applySpeed);
    _applySpeed();

    window.addEventListener("keydown", (ev) => {
      if (ev.code === "Space" && !["INPUT", "TEXTAREA", "SELECT"].includes(ev.target.tagName)) {
        ev.preventDefault();
        togglePlay();
      }
    });

    _updateUI();
  }

  return { init, setSource, removeSource, clearSources, setTime, togglePlay, setViewWindow };
})();
