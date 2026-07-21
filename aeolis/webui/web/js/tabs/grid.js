/* Grid tab: generate the rectangular, equidistant, square-celled model
 * grid by drawing a (rotatable) box on the map or editing a parameter
 * table. Writes x.grd / y.grd through the backend and updates the
 * configuration (xgrid_file, ygrid_file, nx, ny; alfa stays 0).
 *
 * Grid convention (from aeolis.inout.visualize_grid):
 *   column 0 = offshore boundary, last column = onshore,
 *   first/last rows = lateral boundaries; (x0, y0) is the (0,0) corner.
 * Rotation is the CCW angle of the cross-shore axis w.r.t. east.
 */
"use strict";

const GridTab = (() => {

  // draft grid parameters (model CRS)
  let draft = null;            // {x0, y0, dx, nx, ny, rotation}
  let committed = null;        // params of the grid on disk
  let drawArmed = false;
  let editing = false;         // handles visible
  let draggingHandle = false;  // a handle drag is in progress
  let handles = [];            // maplibre markers
  let shearVisible = false;
  let shearUdir = 270;

  const SRC_FILL = "grid-fill";
  const SRC_LINES = "grid-lines";
  const SRC_OUTLINE = "grid-outline";
  const SRC_SHEAR = "grid-shear";

  /* ================= geometry (model CRS) ================= */

  function axes(rotation) {
    const t = rotation * Math.PI / 180;
    return {
      ex: [Math.cos(t), Math.sin(t)],       // cross-shore (i)
      ey: [-Math.sin(t), Math.cos(t)],      // alongshore (j)
    };
  }

  function corner(p, i, j) {
    const { ex, ey } = axes(p.rotation);
    return [
      p.x0 + ex[0] * i * p.dx + ey[0] * j * p.dx,
      p.y0 + ex[1] * i * p.dx + ey[1] * j * p.dx,
    ];
  }

  function outlineRing(p) {
    return [
      corner(p, 0, 0), corner(p, p.nx, 0),
      corner(p, p.nx, p.ny), corner(p, 0, p.ny),
      corner(p, 0, 0),
    ];
  }

  function gridLines(p, maxLines = 50) {
    const lines = [];
    const si = Math.max(1, Math.round(p.nx / maxLines));
    const sj = Math.max(1, Math.round(p.ny / maxLines));
    for (let i = 0; i <= p.nx; i += si) {
      lines.push([corner(p, i, 0), corner(p, i, p.ny)]);
    }
    if (p.nx % si) lines.push([corner(p, p.nx, 0), corner(p, p.nx, p.ny)]);
    for (let j = 0; j <= p.ny; j += sj) {
      lines.push([corner(p, 0, j), corner(p, p.nx, j)]);
    }
    if (p.ny % sj) lines.push([corner(p, 0, p.ny), corner(p, p.nx, p.ny)]);
    return lines;
  }

  function center(p) {
    return corner(p, p.nx / 2, p.ny / 2);
  }

  function isSaved() {
    if (!draft || !committed) return false;
    const close = (a, b, tol) => Math.abs(a - b) <= tol;
    return draft.nx === committed.nx && draft.ny === committed.ny
      && close(draft.dx, committed.dx, 1e-6)
      && close(draft.x0, committed.x0, 1e-4)
      && close(draft.y0, committed.y0, 1e-4)
      && close(draft.rotation, committed.rotation, 1e-6);
  }

  /* ================= rendering ================= */

  function _lineColor() {
    return App.state.ui.basemap === "sat" && !CRS.isLocal() ? "#ffffff" : "#1c2430";
  }

  function render() {
    const map = MapView.instance();
    if (!map || !map.isStyleLoaded()) { map && map.once("idle", render); return; }

    if (!draft) {
      _clearAll();
      return;
    }
    _renderGeometry();
    if (!draggingHandle) _renderHandles();
    _renderShear();
  }

  function _renderGeometry() {
    const map = MapView.instance();
    const ring = outlineRing(draft).map(CRS.toLngLat);
    MapView.upsertGeojson(SRC_FILL, {
      type: "Feature", properties: {},
      geometry: { type: "Polygon", coordinates: [ring] },
    });
    MapView.upsertGeojson(SRC_OUTLINE, {
      type: "Feature", properties: {},
      geometry: { type: "LineString", coordinates: ring },
    });
    MapView.upsertGeojson(SRC_LINES, {
      type: "FeatureCollection",
      features: gridLines(draft).map((line) => ({
        type: "Feature", properties: {},
        geometry: { type: "LineString", coordinates: line.map(CRS.toLngLat) },
      })),
    });

    const color = _lineColor();
    const saved = isSaved();
    MapView.ensureLayer({
      id: SRC_FILL, type: "fill", source: SRC_FILL,
      paint: { "fill-color": color, "fill-opacity": 0.06 },
    });
    MapView.ensureLayer({
      id: SRC_LINES, type: "line", source: SRC_LINES,
      paint: { "line-color": color, "line-width": 0.7, "line-opacity": 0.4 },
    });
    MapView.ensureLayer({
      id: SRC_OUTLINE, type: "line", source: SRC_OUTLINE,
      paint: { "line-color": color, "line-width": 2.2 },
    });
    map.setPaintProperty(SRC_FILL, "fill-color", color);
    map.setPaintProperty(SRC_LINES, "line-color", color);
    map.setPaintProperty(SRC_OUTLINE, "line-color", color);
    // unsaved drafts render dashed
    map.setPaintProperty(SRC_OUTLINE, "line-dasharray", saved ? [1, 0] : [2.5, 1.8]);

    const layerVisible = _gridLayerVisible();
    for (const id of [SRC_FILL, SRC_LINES, SRC_OUTLINE]) {
      map.setLayoutProperty(id, "visibility", layerVisible ? "visible" : "none");
    }
    _renderBoundaryLabels(layerVisible);
  }

  function _gridLayerVisible() {
    const layer = Layers.get("grid-main");
    return !layer || layer.visible !== false;
  }

  function _renderBoundaryLabels(visible) {
    MapView.removeLabels("grid-");
    if (!draft || !visible) return;
    const cfg = App.state.config || {};
    const mid = (a, b) => [(a[0] + b[0]) / 2, (a[1] + b[1]) / 2];
    const p = draft;
    const labels = [
      ["grid-b-offshore", mid(corner(p, 0, 0), corner(p, 0, p.ny)),
        `Offshore: ${cfg.boundary_offshore || "?"}`, "boundary-offshore"],
      ["grid-b-onshore", mid(corner(p, p.nx, 0), corner(p, p.nx, p.ny)),
        `Onshore: ${cfg.boundary_onshore || "?"}`, "boundary-onshore"],
      ["grid-b-lat-a", mid(corner(p, 0, 0), corner(p, p.nx, 0)),
        `Lateral: ${cfg.boundary_lateral || "?"}`, "boundary-lateral"],
      ["grid-b-lat-b", mid(corner(p, 0, p.ny), corner(p, p.nx, p.ny)),
        `Lateral: ${cfg.boundary_lateral || "?"}`, "boundary-lateral"],
      ["grid-c-00", corner(p, 0, 0), "(0,0)", ""],
      ["grid-c-n0", corner(p, p.nx, 0), `(0,${p.nx})`, ""],
      ["grid-c-0m", corner(p, 0, p.ny), `(${p.ny},0)`, ""],
      ["grid-c-nm", corner(p, p.nx, p.ny), `(${p.ny},${p.nx})`, ""],
    ];
    for (const [id, xy, text, cls] of labels) MapView.setLabel(id, xy, text, cls);
  }

  /* ---- edit handles ---- */

  function _clearHandles() {
    for (const marker of handles) marker.remove();
    handles = [];
  }

  function _handlePositions() {
    const { ex } = axes(draft.rotation);
    const midOn = [(corner(draft, draft.nx, 0)[0] + corner(draft, draft.nx, draft.ny)[0]) / 2,
      (corner(draft, draft.nx, 0)[1] + corner(draft, draft.nx, draft.ny)[1]) / 2];
    const offset = Math.max(draft.dx * 2, draft.nx * draft.dx * 0.12);
    return {
      corner: corner(draft, draft.nx, draft.ny),
      origin: corner(draft, 0, 0),
      rotate: [midOn[0] + ex[0] * offset, midOn[1] + ex[1] * offset],
    };
  }

  function _renderHandles() {
    _clearHandles();
    if (!draft || !editing) return;
    const map = MapView.instance();
    const pos = _handlePositions();

    handles.push(_handle("corner", pos.corner, "grid-handle corner", (xy) => {
      const { ex, ey } = axes(draft.rotation);
      const vx = xy[0] - draft.x0, vy = xy[1] - draft.y0;
      draft.nx = Math.max(1, Math.round((vx * ex[0] + vy * ex[1]) / draft.dx));
      draft.ny = Math.max(1, Math.round((vx * ey[0] + vy * ey[1]) / draft.dx));
    }));

    handles.push(_handle("origin", pos.origin, "grid-handle origin", (xy) => {
      draft.x0 = xy[0];
      draft.y0 = xy[1];
    }));

    handles.push(_handle("rotate", pos.rotate, "grid-handle rotate", (xy) => {
      const c = center(draft);
      const angle = Math.atan2(xy[1] - c[1], xy[0] - c[0]) * 180 / Math.PI;
      _rotateAbout(c, angle);
    }, "↻"));

    for (const marker of handles) marker.addTo(map);
  }

  function _rotateAbout(c, newRotation) {
    const p = draft;
    const { ex, ey } = axes(newRotation);
    const hx = p.nx * p.dx / 2, hy = p.ny * p.dx / 2;
    p.rotation = newRotation;
    p.x0 = c[0] - ex[0] * hx - ey[0] * hy;
    p.y0 = c[1] - ex[1] * hx - ey[1] * hy;
  }

  function _handle(kind, modelXY, cls, applyDrag, text = "") {
    const node = U.el("div", { class: cls }, text);
    const marker = new maplibregl.Marker({ element: node, draggable: true, anchor: "center" })
      .setLngLat(CRS.toLngLat(modelXY));
    marker._kind = kind;
    marker.on("dragstart", () => { draggingHandle = true; });
    marker.on("drag", () => {
      applyDrag(CRS.fromLngLat(marker.getLngLat()));
      _syncTable();
      _renderGeometryThrottled();
      _repositionHandles(marker);
    });
    marker.on("dragend", () => {
      draggingHandle = false;
      render();
    });
    return marker;
  }

  /* While dragging one handle, keep the others in place without
   * rebuilding them (rebuilding would kill the active drag). */
  function _repositionHandles(activeMarker) {
    const pos = _handlePositions();
    for (const marker of handles) {
      if (marker === activeMarker) continue;
      marker.setLngLat(CRS.toLngLat(pos[marker._kind]));
    }
  }

  const _renderGeometryThrottled = (() => {
    let pending = false;
    return () => {
      if (pending) return;
      pending = true;
      requestAnimationFrame(() => { pending = false; if (draft) _renderGeometry(); });
    };
  })();

  const renderThrottled = (() => {
    let pending = false;
    return () => {
      if (pending) return;
      pending = true;
      requestAnimationFrame(() => { pending = false; render(); });
    };
  })();

  function _clearAll() {
    _clearHandles();
    MapView.removeLabels("grid-");
    for (const id of [SRC_FILL, SRC_LINES, SRC_OUTLINE, SRC_SHEAR, SRC_SHEAR + "-inner"]) {
      MapView.removeLayerAndSource(id);
    }
  }

  /* ================= shear (2nd computational) grid ================= */

  async function _renderShear() {
    if (!shearVisible || !committed || !isSaved()) {
      // an unsaved draft no longer matches the grid on disk -> the
      // computational grid preview would be misleading
      MapView.removeLayerAndSource(SRC_SHEAR);
      MapView.removeLayerAndSource(SRC_SHEAR + "-inner");
      MapView.removeLabel("grid-shear-info");
      _syncShearInfo(null);
      return;
    }
    try {
      const cfg = App.state.config || {};
      const q = new URLSearchParams({
        udir: shearUdir,
        dx: cfg.dx ?? 1, dy: cfg.dy ?? 1,
        buffer_width: cfg.buffer_width ?? 10,
      });
      const shear = await Api.get(`/api/grid/shear?${q}`);
      const ring = [...shear.ring, shear.ring[0]].map(CRS.toLngLat);
      MapView.upsertGeojson(SRC_SHEAR, {
        type: "Feature", properties: {},
        geometry: { type: "Polygon", coordinates: [ring] },
      });
      MapView.ensureLayer({
        id: SRC_SHEAR, type: "line", source: SRC_SHEAR,
        paint: { "line-color": "#b4423b", "line-width": 1.6, "line-dasharray": [4, 3] },
      });
      if (shear.inner.length) {
        const innerRing = [...shear.inner, shear.inner[0]].map(CRS.toLngLat);
        MapView.upsertGeojson(SRC_SHEAR + "-inner", {
          type: "Feature", properties: {},
          geometry: { type: "Polygon", coordinates: [innerRing] },
        });
        MapView.ensureLayer({
          id: SRC_SHEAR + "-inner", type: "line", source: SRC_SHEAR + "-inner",
          paint: { "line-color": "#b4423b", "line-width": 1, "line-dasharray": [2, 3], "line-opacity": 0.6 },
        });
      }
      MapView.setLabel("grid-shear-info", shear.ring[1],
        `shear grid ${shear.n_cells[0]}×${shear.n_cells[1]} @ udir ${Math.round(shear.udir)}°`, "boundary-lateral");
      _syncShearInfo(shear);
    } catch (err) {
      console.warn("shear preview failed", err.message);
    }
  }

  /* ================= draw interaction ================= */

  function armDraw(force = null) {
    const map = MapView.instance();
    const arm = force !== null ? force : !drawArmed;
    if (!arm) { _disarmDraw(); _syncButtons(); return; }

    drawArmed = true;
    editing = false;
    map.getCanvas().style.cursor = "crosshair";
    map.dragPan.disable();
    _syncButtons();
    U.toast("Click and drag to draw the grid area (Esc to cancel)");

    const onDown = (ev) => {
      if (!drawArmed) return;
      const anchor = CRS.fromLngLat(ev.lngLat);
      const dx = draft ? draft.dx : _defaultDx();
      draft = { x0: anchor[0], y0: anchor[1], dx, nx: 1, ny: 1, rotation: draft ? draft.rotation : 0 };

      const onMove = (mv) => {
        const xy = CRS.fromLngLat(mv.lngLat);
        const { ex, ey } = axes(draft.rotation);
        const vx = xy[0] - anchor[0], vy = xy[1] - anchor[1];
        const w = vx * ex[0] + vy * ex[1];
        const h = vx * ey[0] + vy * ey[1];
        const w0 = Math.min(w, 0), h0 = Math.min(h, 0);
        draft.x0 = anchor[0] + ex[0] * w0 + ey[0] * h0;
        draft.y0 = anchor[1] + ex[1] * w0 + ey[1] * h0;
        draft.nx = Math.max(1, Math.round(Math.abs(w) / draft.dx));
        draft.ny = Math.max(1, Math.round(Math.abs(h) / draft.dx));
        _syncTable();
        _renderGeometryThrottled();
      };
      const onUp = () => {
        map.off("mousemove", onMove);
        map.off("mouseup", onUp);
        _disarmDraw();
        editing = true;
        _syncButtons();
        _syncTable();
        render();
      };
      map.on("mousemove", onMove);
      map.on("mouseup", onUp);
    };
    map.once("mousedown", onDown);
  }

  function _disarmDraw() {
    const map = MapView.instance();
    drawArmed = false;
    map.getCanvas().style.cursor = "";
    map.dragPan.enable();
  }

  function _defaultDx() {
    const map = MapView.instance();
    const bounds = map.getBounds();
    const a = CRS.fromLngLat([bounds.getWest(), bounds.getSouth()]);
    const b = CRS.fromLngLat([bounds.getEast(), bounds.getNorth()]);
    const span = Math.abs(b[0] - a[0]);
    const raw = span / 200;
    const mag = Math.pow(10, Math.floor(Math.log10(raw)));
    return Math.max(mag, Math.round(raw / mag) * mag) || 1;
  }

  /* ================= panel UI ================= */

  let tableInputs = {};
  let buttons = {};
  let shearInfoEl = null;

  function _buildPanel() {
    const panel = document.getElementById("grid-panel");
    U.clear(panel);

    const drawBtn = U.el("button", { class: "toggle" }, "Draw grid on map");
    drawBtn.addEventListener("click", () => armDraw());
    const editBtn = U.el("button", { class: "toggle" }, "Edit handles");
    editBtn.addEventListener("click", () => {
      editing = !editing;
      if (editing && drawArmed) _disarmDraw();
      _syncButtons();
      render();
    });
    buttons = { draw: drawBtn, edit: editBtn };
    panel.append(U.el("div", { class: "btn-row" }, drawBtn, editBtn));

    // parameters table
    const rows = [
      ["x0", "Origin x (0,0 corner)", "m"],
      ["y0", "Origin y", "m"],
      ["dx", "Cell size dx = dy", "m"],
      ["nx", "Cells cross-shore (nx)", "-"],
      ["ny", "Cells alongshore (ny)", "-"],
      ["rotation", "Rotation (CCW from east)", "°"],
    ];
    const table = U.el("table", { class: "data" });
    table.append(U.el("tr", {}, U.el("th", {}, "Parameter"), U.el("th", {}, "Value")));
    tableInputs = {};
    for (const [key, label, unit] of rows) {
      const input = U.el("input", { type: "text" });
      input.addEventListener("keydown", (ev) => {
        if (ev.key === "Enter") input.blur();
      });
      input.addEventListener("blur", () => _commitTable(key, input));
      tableInputs[key] = input;
      table.append(U.el("tr", {},
        U.el("th", { title: `[${unit}]` }, label),
        U.el("td", {}, input)));
    }
    panel.append(U.el("div", { class: "form-group" },
      U.el("span", { class: "fg-label" }, "Grid parameters"), table));

    const derived = U.el("div", { class: "muted", id: "grid-derived" });
    panel.append(derived);

    const saveBtn = U.el("button", { class: "primary" }, "Generate & save x.grd / y.grd");
    saveBtn.addEventListener("click", _saveGrid);
    panel.append(U.el("div", { class: "btn-row" }, saveBtn));

    panel.append(U.el("div", { class: "muted", style: "font-size:12px" },
      "Boundary convention: column 0 = offshore, last column = onshore, first/last row = lateral. ",
      "Boundary types are set in Settings → Boundary conditions."));

    // shear subgrid section (visible when process_shear is on)
    const shearWrap = U.el("div", { class: "form-group", id: "grid-shear-section" });
    shearWrap.append(U.el("span", { class: "fg-label" }, "2nd computational grid (wind shear)"));
    const shearToggle = U.el("input", { type: "checkbox", id: "shear-toggle" });
    shearToggle.checked = shearVisible;
    shearToggle.addEventListener("change", () => { shearVisible = shearToggle.checked; _renderShear(); });
    shearWrap.append(U.el("div", { class: "form-row" },
      U.el("label", { for: "shear-toggle" }, "Show on map"), shearToggle));

    shearInfoEl = U.el("div", { class: "shear-info" });
    shearWrap.append(shearInfoEl);

    const settingsLink = U.el("button", { class: "ghost", style: "font-size:12px" },
      "edit dx / dy / buffer_width in Settings →");
    settingsLink.addEventListener("click", () => {
      Tabs.activate("settings");
      document.getElementById("cfg-search").value = "buffer_width";
      SchemaForm.filter("buffer");
    });
    shearWrap.append(U.el("div", { class: "btn-row" }, settingsLink));

    const udirSlider = U.el("input", { type: "range", min: 0, max: 360, step: 5, value: shearUdir });
    const udirLabel = U.el("label", {}, `example wind dir ${shearUdir}°`);
    udirSlider.addEventListener("input", () => {
      shearUdir = Number(udirSlider.value);
      udirLabel.textContent = `example wind dir ${shearUdir}°`;
      _renderShear();
    });
    shearWrap.append(U.el("div", { class: "form-row" }, udirLabel, udirSlider));
    panel.append(shearWrap);

    _syncShearSection();
    _syncShearInfo(null);
    _syncButtons();
    _syncTable();
  }

  function _syncButtons() {
    if (buttons.draw) buttons.draw.classList.toggle("active", drawArmed);
    if (buttons.edit) buttons.edit.classList.toggle("active", editing);
  }

  function _syncShearSection() {
    const section = document.getElementById("grid-shear-section");
    if (!section) return;
    const cfg = App.state.config || {};
    section.style.display = cfg.process_shear ? "" : "none";
  }

  function _syncShearInfo(shear) {
    if (!shearInfoEl) return;
    const cfg = App.state.config || {};
    const parts = [
      `dx = ${cfg.dx ?? "?"} m, dy = ${cfg.dy ?? "?"} m, buffer = ${cfg.buffer_width ?? "?"} m`,
    ];
    if (shear) {
      parts.push(`≈ ${shear.n_cells[0]} × ${shear.n_cells[1]} cells ` +
        `(${U.fmtNum(shear.length, 4)} × ${U.fmtNum(shear.width, 4)} m at udir ${Math.round(shear.udir)}°)`);
    } else if (!isSaved()) {
      parts.push("preview hidden: grid not saved yet");
    }
    U.clear(shearInfoEl);
    for (const text of parts) {
      shearInfoEl.append(U.el("div", { class: "muted", style: "font-size:12px" }, text));
    }
  }

  function _commitTable(key, input) {
    if (!draft) draft = { x0: 0, y0: 0, dx: 10, nx: 50, ny: 50, rotation: 0 };
    const v = Number(input.value);
    if (!Number.isFinite(v)) { _syncTable(); return; }
    if (key === "nx" || key === "ny") draft[key] = Math.max(1, Math.round(v));
    else if (key === "dx") draft.dx = Math.max(1e-6, v);
    else draft[key] = v;
    editing = true;
    _syncButtons();
    _syncTable();
    render();
  }

  function _syncTable() {
    for (const [key, input] of Object.entries(tableInputs)) {
      if (document.activeElement === input) continue;
      input.value = draft ? U.fmtNum(draft[key], 6) : "";
    }
    const derived = document.getElementById("grid-derived");
    if (derived) {
      if (draft) {
        const lx = draft.nx * draft.dx, ly = draft.ny * draft.dx;
        const saved = isSaved();
        derived.textContent =
          `extent ${U.fmtNum(lx, 4)} × ${U.fmtNum(ly, 4)} m — ${draft.nx * draft.ny} cells` +
          (saved ? "" : "  ⚠ not saved");
        derived.style.color = saved ? "" : "var(--danger)";
      } else {
        derived.textContent = "No grid yet — draw a box or enter parameters.";
      }
    }
  }

  function _zoomToGrid() {
    if (!draft) return;
    const ring = outlineRing(draft);
    const xs = ring.map((c) => c[0]), ys = ring.map((c) => c[1]);
    MapView.fitModelBounds(Math.min(...xs), Math.min(...ys), Math.max(...xs), Math.max(...ys));
  }

  async function _saveGrid() {
    if (!draft) { U.toast("Draw or define a grid first", "error"); return; }
    try {
      const res = await Api.post("/api/grid/save", draft);
      committed = res.params;
      draft = { x0: res.params.x0, y0: res.params.y0, dx: res.params.dx,
        nx: res.params.nx, ny: res.params.ny, rotation: res.params.rotation };
      editing = false;
      const cfg = await Api.get("/api/config");
      App.state.config = cfg.values;
      for (const key of ["xgrid_file", "ygrid_file", "nx", "ny"]) App.emit("config-changed", key);
      Layers.register({ id: "grid-main", group: "grid", title: "Model grid",
        subtitle: `${draft.nx}×${draft.ny}` });
      _syncButtons();
      _syncTable();
      render();
      U.toast(`Saved ${res.files.xgrid_file} / ${res.files.ygrid_file}`, "ok");
    } catch (err) {
      U.toast(`Grid save failed: ${err.message}`, "error");
    }
  }

  /* ================= lifecycle ================= */

  async function _loadExisting() {
    try {
      const res = await Api.get("/api/grid");
      if (res.exists && res.params) {
        committed = res.params;
        draft = { x0: res.params.x0, y0: res.params.y0, dx: res.params.dx,
          nx: res.params.nx, ny: res.params.ny, rotation: res.params.rotation };
        Layers.register({ id: "grid-main", group: "grid", title: "Model grid",
          subtitle: `${draft.nx}×${draft.ny}` });
        if (!res.params.uniform) {
          U.toast("Existing grid is not uniform; table shows an approximation", "error");
        }
        if (!App.state.crsFromState) {
          const detected = CRS.autodetectFromGrid();
          if (detected) {
            App.state.crsFromState = true;
            U.toast(`Coordinate system detected: ${CRS.label()} (click the CRS chip to change)`, "");
          }
        }
        _syncTable();
        render();
      }
    } catch (err) {
      console.warn("grid load failed", err.message);
    }
  }

  function init() {
    Tabs.register("grid", {
      enter: () => { _syncShearSection(); render(); },
      leave: () => { if (drawArmed) { _disarmDraw(); _syncButtons(); } },
    });
    _buildPanel();
    App.on("project", async () => {
      draft = null; committed = null; editing = false;
      _buildPanel();
      await _loadExisting();
      if (draft) _zoomToGrid();
    });
    App.on("config-changed", (key) => {
      if (["boundary_offshore", "boundary_onshore", "boundary_lateral"].includes(key)) render();
      if (key === "process_shear") _syncShearSection();
      if (["dx", "dy", "buffer_width"].includes(key)) { _syncShearInfo(null); _renderShear(); }
    });
    App.on("layer-visibility", (layer) => {
      if (layer.id === "grid-main") render();
    });
    App.on("basemap", () => renderThrottled());
    App.on("crs", () => renderThrottled());
    window.addEventListener("keydown", (ev) => {
      if (ev.key === "Escape" && drawArmed) { _disarmDraw(); _syncButtons(); }
    });
  }

  return { init, params: () => committed };
})();
