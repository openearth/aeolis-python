/* Map view: MapLibre GL wrapper.
 *
 * - Basemaps: grey (CARTO light) / satellite (Esri World Imagery) / none.
 *   In local (conceptual) CRS mode the basemap is forced to "none".
 * - All labels are HTML markers (no glyph server needed -> fully offline
 *   capable except for basemap tiles).
 * - Model <-> map conversion goes through the CRS module.
 */
"use strict";

const MapView = (() => {

  let map = null;
  const markers = new Map();   // id -> maplibregl.Marker

  const BASEMAPS = {
    gray: {
      tiles: [
        "https://a.basemaps.cartocdn.com/light_all/{z}/{x}/{y}.png",
        "https://b.basemaps.cartocdn.com/light_all/{z}/{x}/{y}.png",
        "https://c.basemaps.cartocdn.com/light_all/{z}/{x}/{y}.png",
      ],
      attribution: "© OpenStreetMap contributors © CARTO",
    },
    sat: {
      tiles: [
        "https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}",
      ],
      attribution: "Esri, Maxar, Earthstar Geographics",
    },
  };

  function _style() {
    return {
      version: 8,
      sources: {
        "basemap-gray": { type: "raster", tiles: BASEMAPS.gray.tiles, tileSize: 256, attribution: BASEMAPS.gray.attribution },
        "basemap-sat": { type: "raster", tiles: BASEMAPS.sat.tiles, tileSize: 256, attribution: BASEMAPS.sat.attribution },
      },
      layers: [
        { id: "bg", type: "background", paint: { "background-color": "#e8eaed" } },
        { id: "basemap-gray", type: "raster", source: "basemap-gray", layout: { visibility: "visible" } },
        { id: "basemap-sat", type: "raster", source: "basemap-sat", layout: { visibility: "none" } },
      ],
    };
  }

  function init(containerId) {
    map = new maplibregl.Map({
      container: containerId,
      style: _style(),
      center: [5.2, 52.7],
      zoom: 6.2,
      attributionControl: { compact: true },
      dragRotate: false,
      pitchWithRotate: false,
      touchPitch: false,
    });
    map.addControl(new maplibregl.NavigationControl({ showCompass: false }), "top-left");
    // no scale bar (bottom-left holds the coordinate/CRS readout,
    // bottom-right the conditions toggle + attribution)
    map.touchZoomRotate.disableRotation();

    _buildBasemapSwitch();

    map.on("mousemove", (ev) => {
      const xy = CRS.fromLngLat(ev.lngLat);
      const out = document.getElementById("coord-readout");
      if (out) {
        const digits = App.state.crs.epsg === 4326 ? 5 : 1;
        out.textContent = `x ${xy[0].toFixed(digits)}  y ${xy[1].toFixed(digits)}`;
      }
    });

    App.on("crs", () => {
      _applyCrsMode();
      const crsOut = document.getElementById("crs-readout");
      if (crsOut) crsOut.textContent = CRS.label();
    });

    const crsOut = document.getElementById("crs-readout");
    if (crsOut) {
      crsOut.style.cursor = "pointer";
      crsOut.title = "Click to view/change the coordinate reference system";
      crsOut.addEventListener("click", () => CRS.openSelector());
    }

    return new Promise((resolve) => map.on("load", () => resolve(map)));
  }

  function instance() { return map; }

  /* ---- basemap ---- */

  function setBasemap(name) {
    App.state.ui.basemap = name;
    _refreshBasemapLayers();
    for (const btn of document.querySelectorAll("#basemap-switch button")) {
      btn.classList.toggle("active", btn.dataset.base === name);
    }
    App.touchUi();
    App.emit("basemap", name);
  }

  function _refreshBasemapLayers() {
    if (!map || !map.isStyleLoaded()) {
      if (map) map.once("styledata", _refreshBasemapLayers);
      return;
    }
    const name = CRS.isLocal() ? "none" : App.state.ui.basemap;
    map.setLayoutProperty("basemap-gray", "visibility", name === "gray" ? "visible" : "none");
    map.setLayoutProperty("basemap-sat", "visibility", name === "sat" ? "visible" : "none");
  }

  function _buildBasemapSwitch() {
    const wrap = U.el("div", { id: "basemap-switch" });
    for (const [key, label] of [["gray", "Grey"], ["sat", "Satellite"], ["none", "None"]]) {
      wrap.append(U.el("button", {
        dataset: { base: key },
        class: key === App.state.ui.basemap ? "active" : "",
        onclick: () => setBasemap(key),
      }, label));
    }
    document.getElementById("map-wrap").append(wrap);
  }

  function _applyCrsMode() {
    _refreshBasemapLayers();
    const switchEl = document.getElementById("basemap-switch");
    if (switchEl) switchEl.style.display = CRS.isLocal() ? "none" : "flex";
  }

  /* ---- view helpers ---- */

  function fitModelBounds(minX, minY, maxX, maxY, padding = 60) {
    const a = CRS.toLngLat([minX, minY]);
    const b = CRS.toLngLat([maxX, maxY]);
    map.fitBounds([
      [Math.min(a[0], b[0]), Math.min(a[1], b[1])],
      [Math.max(a[0], b[0]), Math.max(a[1], b[1])],
    ], { padding, duration: 400, maxZoom: 17 });
  }

  /* ---- HTML label markers (no glyphs needed) ---- */

  function setLabel(id, modelXY, text, className = "") {
    removeLabel(id);
    // a misconfigured CRS can map model coords outside the valid
    // lat/lng range; skip the label instead of throwing (an uncaught
    // throw here would break whole render passes)
    const lngLat = CRS.toLngLat(modelXY);
    if (!Number.isFinite(lngLat[0]) || !Number.isFinite(lngLat[1])
        || Math.abs(lngLat[1]) > 89.9) return null;
    const node = U.el("div", { class: `map-label ${className}` }, text);
    const marker = new maplibregl.Marker({ element: node, anchor: "center" })
      .setLngLat(lngLat)
      .addTo(map);
    markers.set(id, marker);
    return marker;
  }

  function removeLabel(id) {
    const marker = markers.get(id);
    if (marker) { marker.remove(); markers.delete(id); }
  }

  function removeLabels(prefix) {
    for (const [id, marker] of Array.from(markers.entries())) {
      if (id.startsWith(prefix)) { marker.remove(); markers.delete(id); }
    }
  }

  /* ---- geojson source/layer helpers ---- */

  function upsertGeojson(sourceId, data) {
    const source = map.getSource(sourceId);
    if (source) source.setData(data);
    else map.addSource(sourceId, { type: "geojson", data });
  }

  function ensureLayer(def) {
    if (!map.getLayer(def.id)) map.addLayer(def);
  }

  function removeLayerAndSource(id) {
    if (map.getLayer(id)) map.removeLayer(id);
    if (map.getSource(id)) map.removeSource(id);
  }

  return {
    init, instance, setBasemap, fitModelBounds,
    setLabel, removeLabel, removeLabels,
    upsertGeojson, ensureLayer, removeLayerAndSource,
  };
})();
