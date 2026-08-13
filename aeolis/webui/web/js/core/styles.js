/* Styles: user-defined colormap/style presets, shared across projects.
 *
 * A preset = {id, name, cmap ("name" or "name!r" = reversed), mode
 * ("cells"|"dots"), dotSize, opacity, min, max} where min/max null
 * means "auto from the layer's data range". Layers reference presets
 * by id (per-project, in ui.layerStyles); editing a preset restyles
 * every layer that uses it. Persisted app-globally via /api/styles.
 */
"use strict";

const Styles = (() => {

  let list = [];

  async function load() {
    try {
      list = (await Api.get("/api/styles")).styles || [];
    } catch (err) {
      console.warn("styles load failed", err.message);
      list = [];
    }
    App.emit("styles");
  }

  function all() { return list; }

  function get(id) { return list.find((s) => s.id === id) || null; }

  async function save(preset) {
    const idx = list.findIndex((s) => s.id === preset.id);
    if (idx >= 0) list[idx] = preset; else list.push(preset);
    await _persist();
  }

  async function remove(id) {
    list = list.filter((s) => s.id !== id);
    await _persist();
  }

  async function _persist() {
    try {
      await Api.post("/api/styles", { styles: list });
    } catch (err) {
      U.toast(`Saving colormaps failed: ${err.message}`, "error");
    }
    App.emit("styles");
  }

  function newId() {
    return `st_${Date.now().toString(36)}_${Math.floor(Math.random() * 1e4)}`;
  }

  return { load, all, get, save, remove, newId };
})();
