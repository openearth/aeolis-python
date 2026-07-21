/* Domain tab: fill the model grid with bathymetry, vegetation and
 * non-erodible layer data.
 *
 * Workflow: pick an area -> check availability -> download raw data
 * (LiDAR / JarKus / Vaklodingen / custom xyz) -> modify with polygons
 * or cell indices -> interpolate onto the grid (writes the .grd file
 * and updates the config). Raw data remains available as layers.
 */
"use strict";

const DomainTab = (() => {

  let overview = null;         // /api/domain payload
  let bounds = null;           // [minx,miny,maxx,maxy] model CRS
  let availability = null;     // check results per source
  let els = {};

  function init() {
    Tabs.register("domain", { enter: _refresh });
    App.on("project", () => { availability = null; bounds = null; _refresh(); });
  }

  async function _refresh() {
    if (!App.state.project) return;
    try {
      overview = await Api.get("/api/domain");
    } catch (err) {
      U.toast(err.message, "error");
      return;
    }
    _build();
    _registerRawLayers();
  }

  function _registerRawLayers() {
    for (const entry of overview.entries) {
      const existing = Layers.get(`raw-${entry.id}`);
      Layers.register({
        id: `raw-${entry.id}`, group: "rawdata",
        title: entry.label || entry.path,
        subtitle: entry.source, entry,
        visible: existing ? existing.visible : false,
      });
    }
    for (const [name, info] of Object.entries(overview.targets)) {
      if (!info.exists) continue;
      const existing = Layers.get(`domain-${name}`);
      Layers.register({
        id: `domain-${name}`, group: "domain",
        title: `${name} (${info.file})`,
        visible: existing ? existing.visible : false,
      });
    }
  }

  /* ================= UI ================= */

  function _build() {
    const panel = document.getElementById("domain-panel");
    U.clear(panel);
    els = {};

    panel.append(_sectionRawData());
    panel.append(_sectionModify());
    panel.append(_sectionInterpolate());
    panel.append(_sectionQuickFlows());
    panel.append(_sectionHistory());
  }

  function _section(title, ...children) {
    const body = U.el("div", { class: "section-body" }, ...children);
    const head = U.el("header", {}, U.el("span", { class: "caret" }, "▾"), title);
    const wrap = U.el("div", { class: "section" }, head, body);
    head.addEventListener("click", () => wrap.classList.toggle("collapsed"));
    return wrap;
  }

  /* ---- 1. raw data ---- */

  function _sectionRawData() {
    const boundsLabel = U.el("div", { class: "muted" }, _boundsText());
    els.boundsLabel = boundsLabel;

    const useGrid = U.el("button", { class: "ghost" }, "Use grid extent");
    useGrid.addEventListener("click", () => {
      const p = GridTab.params();
      if (!p) { U.toast("No grid yet - create one in the Grid tab", "error"); return; }
      const margin = (p.nx * p.dx + p.ny * p.dx) / 2 * 0.1;
      // grid corners (rotation-aware) via its outline
      const ring = _gridRing(p);
      const xs = ring.map((c) => c[0]), ys = ring.map((c) => c[1]);
      bounds = [Math.min(...xs) - margin, Math.min(...ys) - margin,
        Math.max(...xs) + margin, Math.max(...ys) + margin];
      boundsLabel.textContent = _boundsText();
    });

    const drawArea = U.el("button", { class: "ghost" }, "Draw area");
    drawArea.addEventListener("click", async () => {
      try {
        const obj = await Draw.polygon({ name: "data area" });
        const xs = obj.coords.map((c) => c[0]), ys = obj.coords.map((c) => c[1]);
        bounds = [Math.min(...xs), Math.min(...ys), Math.max(...xs), Math.max(...ys)];
        boundsLabel.textContent = _boundsText();
      } catch { /* cancelled */ }
    });

    const checkBtn = U.el("button", { class: "primary" }, "Check availability");
    const progress = _progressBar();
    checkBtn.addEventListener("click", async () => {
      if (!bounds) { U.toast("Set an area first", "error"); return; }
      try {
        checkBtn.disabled = true;
        const res = await Api.post("/api/domain/check", { bounds });
        availability = await Api.waitJob(res.job, progress.update);
        progress.done();
        _renderAvailability();
      } catch (err) {
        U.toast(err.message, "error");
      } finally {
        checkBtn.disabled = false;
      }
    });

    const importXyz = U.el("button", { class: "ghost" }, "Import *.xyz…");
    importXyz.addEventListener("click", async () => {
      const path = await Api.pickFile({
        title: "Import sample file",
        patterns: [["Sample files", "*.xyz;*.txt;*.csv"], ["All files", "*.*"]],
      }).catch(() => null);
      if (!path) return;
      try {
        await Api.post("/api/domain/import_xyz", { path });
        U.toast("Samples imported", "ok");
        _refresh();
      } catch (err) {
        U.toast(err.message, "error");
      }
    });

    els.availability = U.el("div");
    els.rawList = U.el("div");
    _renderRawList();

    return _section("Raw data",
      U.el("div", { class: "btn-row" }, useGrid, drawArea),
      boundsLabel,
      U.el("div", { class: "btn-row" }, checkBtn, importXyz),
      progress.el,
      els.availability,
      U.el("span", { class: "fg-label", style: "margin-top:8px" }, "Downloaded layers"),
      els.rawList,
    );
  }

  function _gridRing(p) {
    const t = p.rotation * Math.PI / 180;
    const ex = [Math.cos(t), Math.sin(t)], ey = [-Math.sin(t), Math.cos(t)];
    const c = (i, j) => [p.x0 + ex[0] * i * p.dx + ey[0] * j * p.dx,
      p.y0 + ex[1] * i * p.dx + ey[1] * j * p.dx];
    return [c(0, 0), c(p.nx, 0), c(p.nx, p.ny), c(0, p.ny)];
  }

  function _boundsText() {
    if (!bounds) return "No area selected yet.";
    return `area: x ${U.fmtNum(bounds[0], 6)} … ${U.fmtNum(bounds[2], 6)}, ` +
      `y ${U.fmtNum(bounds[1], 6)} … ${U.fmtNum(bounds[3], 6)}`;
  }

  function _renderAvailability() {
    const box = els.availability;
    U.clear(box);
    if (!availability) return;
    for (const source of overview.sources) {
      if (source.id === "xyz") continue;
      const res = availability[source.id];
      if (!res) continue;
      const title = U.el("div", { class: "fg-label" }, source.title);
      box.append(title);
      if (res.error) {
        box.append(U.el("div", { class: "muted" }, `⚠ ${res.error}`));
        continue;
      }
      if (!res.available) {
        box.append(U.el("div", { class: "muted" }, res.notes || "not available here"));
        continue;
      }
      const selected = new Set();
      const chips = U.el("div", { class: "year-chips" });
      let estNode;
      for (const info of res.years) {
        const chip = U.el("button", { class: "year-chip", title: U.fmtBytes(info.est_bytes) },
          String(info.year));
        chip.addEventListener("click", () => {
          if (selected.has(info.year)) { selected.delete(info.year); chip.classList.remove("on"); }
          else { selected.add(info.year); chip.classList.add("on"); }
          const est = res.years.filter((y) => selected.has(y.year))
            .reduce((sum, y) => sum + (y.est_bytes || 0), 0);
          estNode.textContent = selected.size
            ? `${selected.size} year(s), ~${U.fmtBytes(est)}` : "";
        });
        chips.append(chip);
      }
      estNode = U.el("span", { class: "muted", style: "margin-left:8px" });
      const dlBtn = U.el("button", { class: "primary" }, "Download");
      const progress = _progressBar();
      dlBtn.addEventListener("click", async () => {
        if (!selected.size) { U.toast("Select years first", "error"); return; }
        try {
          dlBtn.disabled = true;
          const res2 = await Api.post("/api/domain/download", {
            source: source.id, bounds, years: [...selected],
          });
          const out = await Api.waitJob(res2.job, progress.update);
          progress.done();
          U.toast(`Downloaded ${out.entries.length} dataset(s)`, "ok");
          await _refreshEntries();
        } catch (err) {
          U.toast(err.message, "error");
        } finally {
          dlBtn.disabled = false;
        }
      });
      box.append(chips, U.el("div", { class: "btn-row" }, dlBtn, estNode), progress.el);
      if (res.notes) box.append(U.el("div", { class: "muted", style: "font-size:11.5px" }, res.notes));
    }
  }

  async function _refreshEntries() {
    overview = await Api.get("/api/domain");
    _renderRawList();
    _registerRawLayers();
    _renderInterpolateLayers();
  }

  function _renderRawList() {
    const list = els.rawList;
    if (!list) return;
    U.clear(list);
    if (!overview.entries.length) {
      list.append(U.el("div", { class: "muted" }, "Nothing downloaded yet."));
      return;
    }
    for (const entry of overview.entries) {
      const del = U.el("span", { class: "lp-mini", style: "cursor:pointer", title: "Remove" }, "✕");
      del.addEventListener("click", async () => {
        if (!window.confirm(`Remove ${entry.label}? (file stays on disk)`)) return;
        await Api.post("/api/domain/forget", { id: entry.id });
        Layers.unregister(`raw-${entry.id}`);
        _refreshEntries();
      });
      list.append(U.el("div", { class: "lp-row" },
        U.el("span", { class: "lp-name" }, entry.label || entry.path),
        U.el("span", { class: "lp-mini" }, entry.source),
        del));
    }
  }

  /* ---- 2. modify ---- */

  function _sectionModify() {
    const target = _targetSelect();
    const op = U.el("select", {},
      ...["set", "add", "subtract", "multiply", "min", "max"].map((o) =>
        U.el("option", { value: o }, o)));
    const value = U.el("input", { type: "text", placeholder: "value" });
    const init = U.el("input", { type: "text", placeholder: "0 (if file missing)" });

    const polySelect = U.el("select", {});
    const refreshPolys = () => {
      U.clear(polySelect);
      polySelect.append(U.el("option", { value: "" }, "— whole grid —"));
      for (const obj of Objects.byKind("polygon")) {
        polySelect.append(U.el("option", { value: obj.id }, obj.name));
      }
    };
    refreshPolys();
    App.on("objects", refreshPolys);

    const drawBtn = U.el("button", { class: "ghost" }, "Draw new polygon");
    drawBtn.addEventListener("click", async () => {
      try {
        const obj = await Draw.polygon({ name: "modification area" });
        refreshPolys();
        polySelect.value = obj.id;
      } catch { /* cancelled */ }
    });

    const idx = ["j0", "j1", "i0", "i1"].map((ph) =>
      U.el("input", { type: "text", placeholder: ph, style: "width:52px" }));

    const apply = U.el("button", { class: "primary" }, "Apply modification");
    apply.addEventListener("click", async () => {
      const body = {
        target: target.value, op: op.value, value: value.value,
      };
      if (init.value.trim() !== "") body.init = Number(init.value);
      if (polySelect.value) body.polygon = polySelect.value;
      const idxVals = idx.map((n) => n.value.trim());
      if (idxVals.every((v) => v !== "")) body.indices = idxVals.map(Number);
      try {
        const res = await Api.post("/api/domain/modify", body);
        U.toast(`Modified ${res.cells} cells in ${res.file} ` +
          `(range ${U.fmtNum(res.min)} … ${U.fmtNum(res.max)})`, "ok");
        _refresh();
      } catch (err) {
        U.toast(err.message, "error");
      }
    });

    return _section("Modify",
      U.el("div", { class: "form-row" }, U.el("label", {}, "Target"), target),
      U.el("div", { class: "form-row" }, U.el("label", {}, "Operation"), op, value),
      U.el("div", { class: "form-row" }, U.el("label", {}, "Polygon"), polySelect),
      U.el("div", { class: "btn-row" }, drawBtn),
      U.el("div", { class: "form-row" }, U.el("label", {}, "…or cell indices"), ...idx),
      U.el("div", { class: "form-row" }, U.el("label", {}, "Init value"), init),
      U.el("div", { class: "btn-row" }, apply),
    );
  }

  function _targetSelect() {
    const select = U.el("select", {});
    for (const [name, info] of Object.entries(overview.targets)) {
      select.append(U.el("option", { value: name },
        `${name} (${info.file}${info.exists ? "" : " — missing"})`));
    }
    return select;
  }

  /* ---- 3. interpolate ---- */

  function _sectionInterpolate() {
    const target = _targetSelect();
    els.interpTarget = target;
    els.interpLayers = U.el("div");
    const fill = U.el("input", { type: "text", placeholder: "e.g. -20 (optional)" });
    const progress = _progressBar();

    const run = U.el("button", { class: "primary" }, "Interpolate → save .grd");
    run.addEventListener("click", async () => {
      const layers = [...els.interpLayers.querySelectorAll("input:checked")]
        .map((cb) => cb.dataset.entry);
      if (!layers.length) { U.toast("Select at least one source layer", "error"); return; }
      const body = { target: target.value, layers };
      if (fill.value.trim() !== "") body.fill = Number(fill.value);
      try {
        run.disabled = true;
        const res = await Api.post("/api/domain/interpolate", body);
        const out = await Api.waitJob(res.job, progress.update);
        progress.done();
        U.toast(`Wrote ${out.file} (${U.fmtNum(out.min)} … ${U.fmtNum(out.max)} m)`, "ok");
        const cfg = await Api.get("/api/config");
        App.state.config = cfg.values;
        App.emit("config-changed", overview.targets[target.value].config_key);
        _refresh();
      } catch (err) {
        U.toast(err.message, "error");
      } finally {
        run.disabled = false;
      }
    });

    const section = _section("Interpolate to grid",
      U.el("div", { class: "form-row" }, U.el("label", {}, "Target"), target),
      U.el("span", { class: "fg-label" }, "Source layers (priority order = list order)"),
      els.interpLayers,
      U.el("div", { class: "form-row" }, U.el("label", {}, "Fill remaining"), fill),
      U.el("div", { class: "btn-row" }, run),
      progress.el,
    );
    _renderInterpolateLayers();
    return section;
  }

  function _renderInterpolateLayers() {
    const box = els.interpLayers;
    if (!box) return;
    U.clear(box);
    if (!overview.entries.length) {
      box.append(U.el("div", { class: "muted" }, "Download raw data first."));
      return;
    }
    for (const entry of overview.entries) {
      const cb = U.el("input", { type: "checkbox", dataset: { entry: entry.id } });
      box.append(U.el("div", { class: "lp-row" }, cb,
        U.el("span", { class: "lp-name" }, entry.label || entry.path),
        U.el("span", { class: "lp-mini" }, entry.kind)));
    }
  }

  /* ---- 4. quick flows ---- */

  function _sectionQuickFlows() {
    // duplicate bed -> ne with offset
    const offset = U.el("input", { type: "text", value: "-0.5", style: "width:70px" });
    const dupBtn = U.el("button", { class: "ghost" }, "bed → ne-layer");
    dupBtn.addEventListener("click", async () => {
      try {
        const res = await Api.post("/api/domain/duplicate", {
          from: "bed", to: "ne", offset: Number(offset.value) || 0,
        });
        U.toast(`Wrote ${res.file} (bed ${res.offset >= 0 ? "+" : ""}${res.offset} m)`, "ok");
        _refresh();
      } catch (err) {
        U.toast(err.message, "error");
      }
    });

    // vegetation polygon fill
    const vegTarget = U.el("select", {},
      ...Object.keys(overview.targets).filter((t) => ["veg", "hveg", "Nt"].includes(t))
        .map((t) => U.el("option", { value: t }, t)));
    const vegValue = U.el("input", { type: "text", placeholder: "density / height", style: "width:90px" });
    const vegBtn = U.el("button", { class: "ghost" }, "Draw & fill polygon");
    vegBtn.addEventListener("click", async () => {
      const value = Number(vegValue.value);
      if (!Number.isFinite(value)) { U.toast("Enter a fill value first", "error"); return; }
      try {
        const obj = await Draw.polygon({ name: "vegetation" });
        const res = await Api.post("/api/domain/modify", {
          target: vegTarget.value, op: "set", value, polygon: obj.id, init: 0,
        });
        U.toast(`Vegetation set on ${res.cells} cells of ${res.file}`, "ok");
        _refresh();
      } catch (err) {
        if (err.message !== "draw cancelled") U.toast(err.message, "error");
      }
    });

    return _section("Quick flows",
      U.el("div", { class: "form-row" },
        U.el("label", {}, "Duplicate with offset [m]"), offset, dupBtn),
      U.el("div", { class: "form-row" },
        U.el("label", {}, "Vegetation fill"), vegTarget, vegValue),
      U.el("div", { class: "btn-row" }, vegBtn),
    );
  }

  /* ---- 5. history ---- */

  function _sectionHistory() {
    const rows = (overview.history || []).slice(-12).reverse().map((h) =>
      U.el("div", { class: "lp-row" },
        U.el("span", { class: "lp-name" },
          h.action === "interpolate" ? `interpolate → ${h.file}` :
            h.action === "duplicate" ? `${h.from} → ${h.to} (${h.offset} m)` :
              `${h.op} ${h.value} on ${h.target} (${h.cells} cells)`),
        U.el("span", { class: "lp-mini" }, h.time || "")));
    const section = _section("History",
      rows.length ? U.el("div", {}, ...rows)
        : U.el("div", { class: "muted" }, "No operations yet."));
    section.classList.add("collapsed");
    return section;
  }

  /* ---- helpers ---- */

  function _progressBar() {
    const bar = U.el("div");
    const wrap = U.el("div", { class: "progress", style: "display:none" }, bar);
    const msg = U.el("div", { class: "muted", style: "font-size:11.5px" });
    const el = U.el("div", {}, wrap, msg);
    return {
      el,
      update: (job) => {
        wrap.style.display = "";
        bar.style.width = `${Math.round(Math.max(0, job.progress) * 100)}%`;
        msg.textContent = job.message || "";
      },
      done: () => {
        wrap.style.display = "none";
        msg.textContent = "";
      },
    };
  }

  return { init };
})();
