/* Coordinate reference system handling.
 *
 * Model data lives in a projected CRS (default Dutch RD New, EPSG:28992)
 * or in a "local" (fictional/conceptual) frame in plain meters.
 *
 * The map itself is standard Web-Mercator. Chain:
 *   model (projected CRS) <-> WGS84 <-> map
 * In local mode, model meters are interpreted directly as Web-Mercator
 * meters (true-to-scale at the equator), and no basemap is shown.
 */
"use strict";

const CRS = (() => {

  const DEFS = {
    28992: {
      name: "Amersfoort / RD New",
      proj: "+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 " +
        "+k=0.9999079 +x_0=155000 +y_0=463000 +ellps=bessel " +
        "+towgs84=565.417,50.3319,465.552,-0.398957,0.343988,-1.8774,4.0725 " +
        "+units=m +no_defs",
    },
    4326: { name: "WGS 84 (lon/lat)", proj: "+proj=longlat +datum=WGS84 +no_defs" },
    3857: { name: "Web Mercator", proj: "+proj=merc +a=6378137 +b=6378137 +lat_ts=0 +lon_0=0 " +
      "+x_0=0 +y_0=0 +k=1 +units=m +nadgrids=@null +no_defs" },
    25831: { name: "ETRS89 / UTM 31N", proj: "+proj=utm +zone=31 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" },
    25832: { name: "ETRS89 / UTM 32N", proj: "+proj=utm +zone=32 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" },
    32631: { name: "WGS 84 / UTM 31N", proj: "+proj=utm +zone=31 +datum=WGS84 +units=m +no_defs" },
    32632: { name: "WGS 84 / UTM 32N", proj: "+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs" },
  };

  for (const [code, def] of Object.entries(DEFS)) {
    proj4.defs(`EPSG:${code}`, def.proj);
  }

  function list() {
    return Object.entries(DEFS).map(([code, def]) => ({ epsg: Number(code), name: def.name }));
  }

  function isLocal() { return App.state.crs.mode === "local"; }

  function set(mode, epsg = null) {
    App.state.crs = { mode, epsg: mode === "projected" ? (epsg || 28992) : null };
    App.emit("crs", App.state.crs);
    App.touchUi();
  }

  /* model [x,y] -> [lng,lat] for the map */
  function toLngLat(xy) {
    if (isLocal()) return proj4("EPSG:3857", "EPSG:4326", [xy[0], xy[1]]);
    const epsg = App.state.crs.epsg;
    if (epsg === 4326) return [xy[0], xy[1]];
    return proj4(`EPSG:${epsg}`, "EPSG:4326", [xy[0], xy[1]]);
  }

  /* map [lng,lat] -> model [x,y] */
  function fromLngLat(lngLat) {
    const ll = Array.isArray(lngLat) ? lngLat : [lngLat.lng, lngLat.lat];
    if (isLocal()) return proj4("EPSG:4326", "EPSG:3857", ll);
    const epsg = App.state.crs.epsg;
    if (epsg === 4326) return ll;
    return proj4("EPSG:4326", `EPSG:${epsg}`, ll);
  }

  /* Heuristic CRS detection from model coordinate ranges. */
  function detect(xs, ys) {
    const x = median(xs), y = median(ys);
    if (Math.abs(x) <= 180 && Math.abs(y) <= 90 && Math.abs(x) > 2 && Math.abs(y) > 2) {
      return { mode: "projected", epsg: 4326 };
    }
    if (x > -7000 && x < 300000 && y > 289000 && y < 629000) {
      return { mode: "projected", epsg: 28992 };  // RD New validity window
    }
    if (x > 100000 && x < 900000 && y > 3000000 && y < 9000000) {
      return { mode: "projected", epsg: 25831 };  // some UTM zone; 31N default
    }
    return { mode: "local", epsg: null };
  }

  function median(values) {
    const sorted = Array.from(values).filter(Number.isFinite).sort((a, b) => a - b);
    return sorted.length ? sorted[Math.floor(sorted.length / 2)] : NaN;
  }

  function label() {
    if (isLocal()) return "local (conceptual)";
    const epsg = App.state.crs.epsg;
    const def = DEFS[epsg];
    return def ? `EPSG:${epsg} — ${def.name}` : `EPSG:${epsg}`;
  }

  /* Detect from the current grid (if any) and apply. */
  function autodetectFromGrid() {
    const p = (typeof GridTab !== "undefined") ? GridTab.params() : null;
    if (!p) return null;
    const detected = detect([p.x0, p.x0 + p.nx * p.dx], [p.y0, p.y0 + p.ny * p.dx]);
    set(detected.mode, detected.epsg);
    return detected;
  }

  /* Small selector popup, opened from the CRS chip on the map. */
  function openSelector() {
    const popup = Popup.open({ title: "Coordinate reference system", width: 440 });
    const current = App.state.crs;

    const modeProjected = U.el("input", { type: "radio", name: "crs-mode", id: "crs-projected" });
    const modeLocal = U.el("input", { type: "radio", name: "crs-mode", id: "crs-local" });
    (current.mode === "local" ? modeLocal : modeProjected).checked = true;

    const epsgSelect = U.el("select", {});
    for (const item of list()) {
      epsgSelect.append(U.el("option", {
        value: item.epsg,
        selected: item.epsg === current.epsg ? "" : null,
      }, `EPSG:${item.epsg} — ${item.name}`));
    }

    const detectBtn = U.el("button", { class: "ghost" }, "Auto-detect from grid");
    const detectNote = U.el("span", { class: "muted", style: "font-size:12px" });
    detectBtn.addEventListener("click", () => {
      const detected = autodetectFromGrid();
      if (!detected) { detectNote.textContent = "no grid available"; return; }
      if (detected.mode === "local") { modeLocal.checked = true; }
      else {
        modeProjected.checked = true;
        epsgSelect.value = String(detected.epsg);
      }
      detectNote.textContent = detected.mode === "local"
        ? "detected: local/conceptual" : `detected: EPSG:${detected.epsg}`;
    });

    const applyBtn = U.el("button", { class: "primary" }, "Apply");
    applyBtn.addEventListener("click", () => {
      if (modeLocal.checked) set("local");
      else set("projected", Number(epsgSelect.value));
      popup.close();
      U.toast(`CRS: ${label()}`, "ok");
    });

    popup.body.append(
      U.el("div", { class: "form-row" },
        modeProjected, U.el("label", { for: "crs-projected" }, "Projected (real-world)"), epsgSelect),
      U.el("div", { class: "form-row" },
        modeLocal, U.el("label", { for: "crs-local" },
          "Local / conceptual (model meters, no basemap)")),
      U.el("div", { class: "btn-row" }, detectBtn, detectNote),
      U.el("div", { class: "muted", style: "font-size:11.5px;margin-top:6px" },
        "The CRS tells the GUI how grid coordinates map onto the world. ",
        "Changing it re-projects all layers; it does not alter any files."),
      U.el("div", { class: "btn-row" }, applyBtn),
    );
  }

  return { list, set, isLocal, toLngLat, fromLngLat, detect, label,
    autodetectFromGrid, openSelector };
})();
