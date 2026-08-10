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

    // keep the in-progress drawing above the WebGL field layers
    td.on("change", () => { if (active) _raiseDrawLayers(); });

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
      _hideBanner();
      resolve(obj);
    });
    return td;
  }

  /* on-map instruction banner shown while drawing (Terra Draw already
   * finishes a line on Enter/double-click and cancels on Escape - the
   * banner just makes that discoverable, with clickable Finish/Cancel). */
  let banner = null;

  function _showBanner(kind) {
    _hideBanner();
    const isLine = kind === "transect";
    const text = isLine
      ? "Click to add points along the transect — press Enter or double-click to finish"
      : "Click to add corners — click the first point (or double-click) to close the area";
    const finishBtn = U.el("button", { class: "primary" }, "Finish");
    finishBtn.addEventListener("click", _finishActive);
    const cancelBtn = U.el("button", { class: "ghost" }, "Cancel");
    cancelBtn.addEventListener("click", cancel);
    banner = U.el("div", { class: "draw-banner" },
      U.el("span", {}, text),
      U.el("span", { class: "draw-banner-keys" }, "Esc to cancel"),
      finishBtn, cancelBtn);
    document.body.append(banner);
  }

  function _hideBanner() {
    if (banner) { banner.remove(); banner = null; }
  }

  // move Terra Draw's own layers above ours (field/point/object layers)
  // so the shape being drawn is always visible on top
  function _raiseDrawLayers() {
    const map = MapView.instance();
    if (!map) return;
    let style;
    try { style = map.getStyle(); } catch (e) { return; }
    const ours = /^(field-|pts-|objects-|grid-|basemap|bg$)/;
    for (const l of (style.layers || [])) {
      if (!ours.test(l.id)) { try { map.moveLayer(l.id); } catch (e) { /* ignore */ } }
    }
  }

  // finalize the in-progress feature directly from the snapshot (reliable
  // even if the synthetic Enter key does not reach Terra Draw)
  function _finishActive() {
    if (!td || !active) return;
    const feat = td.getSnapshot().find((f) => f.geometry &&
      (f.geometry.type === "LineString" || f.geometry.type === "Polygon"));
    if (!feat) return;
    let coords = feat.geometry.type === "Polygon"
      ? feat.geometry.coordinates[0].slice(0, -1).map(CRS.fromLngLat)
      : feat.geometry.coordinates.map(CRS.fromLngLat);
    // drop consecutive duplicate points (incl. a trailing cursor point)
    coords = coords.filter((c, i) => i === 0 || c[0] !== coords[i - 1][0] || c[1] !== coords[i - 1][1]);
    const need = feat.geometry.type === "Polygon" ? 3 : 2;
    if (coords.length < need) { U.toast(`Add at least ${need} points`, ""); return; }
    const { resolve, kind, opts } = active;
    active = null;
    td.setMode("static");
    td.clear();
    MapView.instance().getCanvas().style.cursor = "";
    _hideBanner();
    resolve(Objects.add({ kind, coords, ...(opts || {}) }));
  }

  function _begin(mode, kind, opts) {
    _ensure();
    cancel();
    MapView.instance().getCanvas().style.cursor = "crosshair";
    td.setMode(mode);
    _showBanner(kind);
    return new Promise((resolve, reject) => {
      active = { resolve, reject, kind, opts };
    });
  }

  function polygon(opts = {}) { return _begin("polygon", "polygon", opts); }
  function transect(opts = {}) { return _begin("linestring", "transect", opts); }

  function cancel() {
    _hideBanner();
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
    // transects always sit on top of the data layers
    if (map.getLayer(SRC_OBJECT_LINES)) map.moveLayer(SRC_OBJECT_LINES);
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
