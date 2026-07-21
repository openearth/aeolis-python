/* Drawing tools (Terra Draw wrapper) + rendering of the Objects store.
 *
 * Any tab can request a drawn shape:
 *   const obj = await Draw.polygon({name: "veg area"});
 *   const obj = await Draw.transect();
 * The finished shape is stored in the Objects store (model CRS
 * coordinates) and rendered by this module on every objects change.
 */
"use strict";

const Draw = (() => {

  let td = null;              // TerraDraw instance
  let active = null;          // {resolve, reject, kind, opts}

  const SRC_OBJECTS = "objects-fill";
  const SRC_OBJECT_LINES = "objects-lines";

  function _ensure() {
    if (td) return td;
    const map = MapView.instance();
    const adapter = new terraDrawMaplibreGlAdapter.TerraDrawMapLibreGLAdapter({ map });
    td = new terraDraw.TerraDraw({
      adapter,
      modes: [
        new terraDraw.TerraDrawPolygonMode({
          styles: {
            fillColor: "#0f766e", fillOpacity: 0.15,
            outlineColor: "#0f766e", outlineWidth: 2,
            closingPointColor: "#ffffff", closingPointOutlineColor: "#0f766e",
          },
        }),
        new terraDraw.TerraDrawLineStringMode({
          styles: { lineStringColor: "#b4423b", lineStringWidth: 2.5 },
        }),
      ],
    });
    td.start();
    td.setMode("static");

    td.on("finish", (id, context) => {
      if (!active || (context && context.action && context.action !== "draw")) return;
      const snapshot = td.getSnapshot().find((f) => f.id === id);
      if (!snapshot) return;
      const { resolve, kind, opts } = active;
      active = null;
      td.setMode("static");
      td.removeFeatures([id]);
      MapView.instance().getCanvas().style.cursor = "";

      let coords;
      if (snapshot.geometry.type === "Polygon") {
        coords = snapshot.geometry.coordinates[0].slice(0, -1).map(CRS.fromLngLat);
      } else {
        coords = snapshot.geometry.coordinates.map(CRS.fromLngLat);
      }
      const obj = Objects.add({ kind, coords, ...(opts || {}) });
      resolve(obj);
    });
    return td;
  }

  function _begin(mode, kind, opts) {
    _ensure();
    cancel();
    MapView.instance().getCanvas().style.cursor = "crosshair";
    td.setMode(mode);
    return new Promise((resolve, reject) => {
      active = { resolve, reject, kind, opts };
    });
  }

  function polygon(opts = {}) { return _begin("polygon", "polygon", opts); }
  function transect(opts = {}) { return _begin("linestring", "transect", opts); }

  function cancel() {
    if (!td) return;
    if (active) {
      active.reject(new Error("draw cancelled"));
      active = null;
    }
    td.setMode("static");
    td.clear();
    const map = MapView.instance();
    if (map) map.getCanvas().style.cursor = "";
  }

  /* ---- render the Objects store ---- */

  function renderObjects() {
    const map = MapView.instance();
    if (!map || !map.isStyleLoaded()) { map && map.once("idle", renderObjects); return; }

    const polys = Objects.toGeojson((o) => o.kind !== "transect");
    const lines = Objects.toGeojson((o) => o.kind === "transect");

    MapView.upsertGeojson(SRC_OBJECTS, polys);
    MapView.ensureLayer({
      id: SRC_OBJECTS, type: "fill", source: SRC_OBJECTS,
      paint: { "fill-color": ["get", "color"], "fill-opacity": 0.18 },
    });
    MapView.ensureLayer({
      id: SRC_OBJECTS + "-outline", type: "line", source: SRC_OBJECTS,
      paint: { "line-color": ["get", "color"], "line-width": 2 },
    });

    MapView.upsertGeojson(SRC_OBJECT_LINES, lines);
    MapView.ensureLayer({
      id: SRC_OBJECT_LINES, type: "line", source: SRC_OBJECT_LINES,
      paint: { "line-color": ["get", "color"], "line-width": 2.5, "line-dasharray": [3, 2] },
    });
  }

  function init() {
    App.on("objects", renderObjects);
    App.on("crs", renderObjects);
    window.addEventListener("keydown", (ev) => {
      if (ev.key === "Escape" && active) cancel();
    });
  }

  return { init, polygon, transect, cancel, renderObjects };
})();
