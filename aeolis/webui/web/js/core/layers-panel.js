/* Central layer registry.
 *
 * Tabs and the viewer register displayable layers here; the Viewer tab
 * renders the manageable tree (visibility + order + styling). The
 * registry keeps insertion order; the Viewer may reorder entries, which
 * drives the maplibre draw order.
 */
"use strict";

const Layers = (() => {

  function register(layer) {
    // layer: {id, group, title, subtitle?, visible?, entry?, fieldId?}
    const existing = App.state.layers.findIndex((l) => l.id === layer.id);
    if (existing >= 0) App.state.layers[existing] = { ...App.state.layers[existing], ...layer };
    else App.state.layers.push({ visible: true, ...layer });
    App.emit("layers", layer.id);
  }

  function unregister(id) {
    const idx = App.state.layers.findIndex((l) => l.id === id);
    if (idx >= 0) {
      App.state.layers.splice(idx, 1);
      App.emit("layers", id);
    }
  }

  function get(id) { return App.state.layers.find((l) => l.id === id) || null; }

  function byGroup(group) { return App.state.layers.filter((l) => l.group === group); }

  function move(id, direction) {
    const list = App.state.layers;
    const idx = list.findIndex((l) => l.id === id);
    if (idx < 0) return;
    const target = idx + direction;
    if (target < 0 || target >= list.length) return;
    // only swap within the same group so the tree stays organized
    if (list[target].group !== list[idx].group) return;
    [list[idx], list[target]] = [list[target], list[idx]];
    App.emit("layers", id);
    App.emit("layer-order");
  }

  return { register, unregister, get, byGroup, move };
})();
