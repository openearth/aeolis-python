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
    // update in place so in-flight references (e.g. a layer that is
    // still loading) keep seeing the current state. Only emit when
    // something actually changed: register() is also called from render
    // passes that themselves re-render on the "layers" event, and an
    // unconditional emit would recurse.
    const existing = App.state.layers.find((l) => l.id === layer.id);
    if (existing) {
      let changed = false;
      for (const [key, val] of Object.entries(layer)) {
        if (existing[key] !== val) { existing[key] = val; changed = true; }
      }
      if (changed) App.emit("layers", layer.id);
      return;
    }
    App.state.layers.push({ visible: true, ...layer });
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

  /* Replace the order of one group's layers (drag & drop reorder). */
  function reorderGroup(group, orderedIds) {
    const list = App.state.layers;
    const entries = list.filter((l) => l.group === group);
    const byId = new Map(entries.map((l) => [l.id, l]));
    const reordered = orderedIds.map((id) => byId.get(id)).filter(Boolean);
    for (const entry of entries) {
      if (!reordered.includes(entry)) reordered.push(entry);
    }
    let k = 0;
    for (let i = 0; i < list.length; i += 1) {
      if (list[i].group === group) list[i] = reordered[k++];
    }
    App.emit("layers", group);
    App.emit("layer-order");
  }

  return { register, unregister, get, byGroup, move, reorderGroup };
})();
