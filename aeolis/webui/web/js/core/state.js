/* Global application state + a minimal pub/sub event bus.
 *
 * Modules read/write App.state directly and emit named events so
 * dependent views can refresh. Persistable UI state is mirrored to the
 * backend (gui/project.json) with a debounce.
 */
"use strict";

const App = (() => {

  const state = {
    project: null,          // {configfile, root, name} or null
    schema: null,           // config schema (sections/params) from backend
    config: null,           // current aeolis.txt values {key: value}
    configDirty: false,

    crs: { mode: "projected", epsg: 28992 },  // mode: projected | local

    clock: {                // shared time model (seconds since refdate epoch)
      playing: false,
      t: null,              // current absolute time [s since epoch]
      t0: null, t1: null,   // range
      speed: 3600,          // model seconds per wall second
    },

    layers: [],             // registered viewer layers (Phase 7)
    objects: [],            // polygons/transects (objects store mirror)

    ui: {                   // persisted per-project UI state
      basemap: "gray",
      theme: "light",       // album-cover theme id (see theme.js)
      // responsive defaults; clamped to the actual screen on restore
      graphsHeight: Math.round(Math.min(280, Math.max(150, window.innerHeight * 0.24))),
      sidebarWidth: Math.round(Math.min(400, Math.max(280, window.innerWidth * 0.24))),
      collapsed: {},
    },
  };

  const listeners = new Map();

  function on(event, fn) {
    if (!listeners.has(event)) listeners.set(event, new Set());
    listeners.get(event).add(fn);
    return () => listeners.get(event).delete(fn);
  }

  function emit(event, payload = null) {
    for (const fn of listeners.get(event) || []) {
      try { fn(payload); } catch (err) { console.error(`listener for ${event} failed`, err); }
    }
    for (const fn of listeners.get("*") || []) {
      try { fn(event, payload); } catch (err) { console.error("wildcard listener failed", err); }
    }
  }

  /* ---- persisted UI state (gui/project.json) ---- */

  const saveUiState = U.debounce(async () => {
    if (!state.project) return;
    try {
      const st = await Api.get("/api/project/state") || {};
      st.ui = state.ui;
      st.crs = state.crs;
      await Api.post("/api/project/state", st);
    } catch (err) {
      console.warn("could not persist UI state", err);
    }
  }, 800);

  function touchUi() { saveUiState(); }

  return { state, on, emit, touchUi };
})();
