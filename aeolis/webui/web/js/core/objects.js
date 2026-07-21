/* Objects store: the single source of truth for user-drawn shapes
 * (polygons, boxes, transects). Other modules (grid, domain, viewer)
 * reference objects from this store by id; the Layers panel only
 * toggles visibility. Persisted to gui/polygons.json via the backend.
 *
 * Object shape:
 *   { id, name, kind: 'polygon'|'box'|'transect', color, visible,
 *     coords: [[x,y], ...]   // model CRS coordinates
 *     props: {...} }
 */
"use strict";

const Objects = (() => {

  let seq = 1;

  function all() { return App.state.objects; }

  function get(id) { return App.state.objects.find((o) => o.id === id) || null; }

  function byKind(kind) { return App.state.objects.filter((o) => o.kind === kind); }

  function add(obj) {
    obj.id = obj.id || `obj-${Date.now().toString(36)}-${seq++}`;
    obj.name = obj.name || _defaultName(obj.kind);
    obj.color = obj.color || _nextColor();
    obj.visible = obj.visible !== false;
    obj.props = obj.props || {};
    App.state.objects.push(obj);
    App.emit("objects", { type: "add", id: obj.id });
    persist();
    return obj;
  }

  function update(id, patch) {
    const obj = get(id);
    if (!obj) return null;
    Object.assign(obj, patch);
    App.emit("objects", { type: "update", id });
    persist();
    return obj;
  }

  function remove(id) {
    const idx = App.state.objects.findIndex((o) => o.id === id);
    if (idx >= 0) {
      App.state.objects.splice(idx, 1);
      App.emit("objects", { type: "remove", id });
      persist();
    }
  }

  function _defaultName(kind) {
    const count = byKind(kind).length + 1;
    return `${kind} ${count}`;
  }

  const PALETTE = ["#e6552f", "#2f7fe6", "#27a355", "#a034c6", "#e0a020", "#12a5b5", "#d1387f"];
  function _nextColor() { return PALETTE[App.state.objects.length % PALETTE.length]; }

  /* ---- persistence ---- */

  const persist = U.debounce(async () => {
    if (!App.state.project) return;
    try {
      await Api.post("/api/objects/save", { objects: App.state.objects });
    } catch (err) {
      console.warn("objects not persisted (endpoint pending)", err.message);
    }
  }, 600);

  async function load() {
    if (!App.state.project) return;
    try {
      const res = await Api.get("/api/objects");
      App.state.objects = res.objects || [];
      App.emit("objects", { type: "load" });
    } catch (err) {
      App.state.objects = [];
      console.warn("objects not loaded (endpoint pending)", err.message);
    }
  }

  /* GeoJSON (in model CRS -> map lnglat) for rendering */
  function toGeojson(filter = null) {
    const features = [];
    for (const obj of App.state.objects) {
      if (!obj.visible) continue;
      if (filter && !filter(obj)) continue;
      const ring = obj.coords.map((xy) => CRS.toLngLat(xy));
      if (obj.kind === "transect") {
        features.push({ type: "Feature", properties: { id: obj.id, color: obj.color, kind: obj.kind },
          geometry: { type: "LineString", coordinates: ring } });
      } else {
        features.push({ type: "Feature", properties: { id: obj.id, color: obj.color, kind: obj.kind },
          geometry: { type: "Polygon", coordinates: [[...ring, ring[0]]] } });
      }
    }
    return { type: "FeatureCollection", features };
  }

  return { all, get, byKind, add, update, remove, load, toGeojson };
})();
