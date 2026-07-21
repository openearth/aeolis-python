/* Domain tab: bathymetry / vegetation / ne-layer data.
 *
 * Two sections:
 *  - Sample data: downloaded/imported datasets with visibility, rename,
 *    reorder, duplicate, remove and a Modify… popup (scope: all /
 *    polygon / indices; op: set/add/subtract/multiply/min/max; save or
 *    save-as-new-layer).
 *  - Interpolated data: the model .grd files with visibility, an
 *    Interpolate… popup (sample priority + fill value) and a
 *    convert-to-samples action. Staleness vs the current grid is
 *    badged.
 * Downloading happens in a popup wizard: pick area (grid extent +
 * buffer, or draw), check availability, multi-select years across
 * sources (click / Shift+click / drag), download all at once.
 */
"use strict";

const DomainTab = (() => {

  let overview = null;

  function init() {
    Tabs.register("domain", { enter: _refresh });
    App.on("project", _refresh);
  }

  async function _refresh() {
    if (!App.state.project) return;
    try {
      overview = await Api.get("/api/domain");
    } catch (err) {
      U.toast(err.message, "error");
      return;
    }
    _registerLayers();
    _build();
  }

  function _registerLayers() {
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

  /* ================= panel ================= */

  function _build() {
    const panel = document.getElementById("domain-panel");
    U.clear(panel);

    // --- sample data ---
    const samples = U.section("Sample data", { count: overview.entries.length });
    samples.body.append(
      U.el("div", { class: "tbtn-row" },
        U.tbtn("download", "Download", { primary: true, title: "Download data from Dutch coastal sources", onclick: _downloadWizard }),
        U.tbtn("upload", "Import", { title: "Import a *.xyz sample file", onclick: _importXyz })),
      _sampleList(),
    );
    panel.append(samples.wrap);

    // --- interpolated data ---
    const targetCount = Object.values(overview.targets).filter((t) => t.exists).length;
    const targets = U.section("Interpolated data (.grd)", { count: targetCount });
    targets.body.append(_targetList());
    if (!overview.grid_available) {
      targets.body.append(U.el("div", { class: "muted", style: "margin-top:8px" },
        "⚠ No model grid yet — create one in the Grid tab before interpolating."));
    }
    panel.append(targets.wrap);
  }

  /* ---- sample list ---- */

  function _sampleList() {
    const list = U.el("div", { class: "layer-tree" });
    if (!overview.entries.length) {
      list.append(U.el("div", { class: "muted" }, "Nothing downloaded or imported yet."));
      return list;
    }
    overview.entries.forEach((entry, idx) => {
      const layerId = `raw-${entry.id}`;
      const layer = Layers.get(layerId) || { visible: false };

      const eye = U.el("span", { class: `eye ${layer.visible ? "" : "off"}`, title: "Show/hide" }, "👁");
      eye.addEventListener("click", () => {
        layer.visible = !layer.visible;
        App.emit("layer-visibility", layer);
        _build();
      });

      const name = U.el("span", { class: "lp-name", title: "Double-click to rename" },
        entry.label || entry.path);
      name.addEventListener("dblclick", () => _renameSample(entry, name));

      const row = U.el("div", { class: "lp-row" }, eye, name,
        U.el("span", { class: "lp-mini" }, entry.source));

      const up = U.miniBtn("up", "Raise (drawn on top)", () => _reorderSample(idx, idx - 1));
      const down = U.miniBtn("down", "Lower", () => _reorderSample(idx, idx + 1));
      up.disabled = idx === 0;
      down.disabled = idx === overview.entries.length - 1;
      row.append(up, down,
        U.miniBtn("copy", "Duplicate", async () => {
          await Api.post("/api/domain/sample_duplicate", { id: entry.id });
          _refresh();
        }),
        U.miniBtn("modify", "Modify…", () => _modifyWizard(entry)),
        _deleteSampleBtn(entry, layerId));
      list.append(row);
    });
    return list;
  }

  function _deleteSampleBtn(entry, layerId) {
    const btn = U.miniBtn("trash", "Remove…", () => {
      const popup = Popup.open({ title: `Remove ${entry.label || entry.path}`, width: 420 });
      const delFile = U.el("input", { type: "checkbox", id: "del-file", checked: "" });
      const okBtn = U.el("button", { class: "danger" }, "Remove");
      okBtn.addEventListener("click", async () => {
        try {
          await Api.post("/api/domain/forget", { id: entry.id, delete_file: delFile.checked });
          Layers.unregister(layerId);
          popup.close();
          _refresh();
        } catch (err) {
          U.toast(err.message, "error");
        }
      });
      const cancelBtn = U.el("button", { class: "ghost" }, "Cancel");
      cancelBtn.addEventListener("click", popup.close);
      popup.body.append(
        U.el("div", { style: "font-size:13px" }, "Remove this sample layer from the project?"),
        U.el("div", { class: "choice-row", style: "margin-top:8px" },
          delFile, U.el("label", { for: "del-file" }, `also delete the file from disk (${entry.path})`)),
        U.el("div", { class: "btn-row", style: "justify-content:flex-end" }, cancelBtn, okBtn),
      );
    });
    btn.classList.add("danger-hover");
    return btn;
  }

  async function _reorderSample(from, to) {
    const ids = overview.entries.map((e) => e.id);
    [ids[from], ids[to]] = [ids[to], ids[from]];
    await Api.post("/api/domain/sample_order", { ids });
    // mirror the order in the layer registry
    await _refresh();
    App.emit("layer-order");
  }

  function _renameSample(entry, nameNode) {
    const input = U.el("input", { type: "text", value: entry.label || "", style: "flex:1;font-size:12px" });
    nameNode.replaceWith(input);
    input.focus(); input.select();
    const commit = async () => {
      const name = input.value.trim();
      if (name && name !== entry.label) {
        await Api.post("/api/domain/sample_rename", { id: entry.id, name });
      }
      _refresh();
    };
    input.addEventListener("blur", commit);
    input.addEventListener("keydown", (ev) => {
      if (ev.key === "Enter") input.blur();
      if (ev.key === "Escape") { input.value = entry.label; input.blur(); }
    });
  }

  async function _importXyz() {
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
  }

  /* ---- interpolated targets ---- */

  function _targetList() {
    const list = U.el("div", { class: "layer-tree" });
    for (const [name, info] of Object.entries(overview.targets)) {
      const layerId = `domain-${name}`;
      const layer = Layers.get(layerId) || { visible: false, id: layerId };

      const eye = U.el("span", {
        class: `eye ${layer.visible ? "" : "off"}`,
        title: info.exists ? "Show/hide" : "File does not exist yet",
      }, "👁");
      if (info.exists) {
        eye.addEventListener("click", () => {
          layer.visible = !layer.visible;
          App.emit("layer-visibility", layer);
          _build();
        });
      } else {
        eye.style.opacity = "0.25";
      }

      const row = U.el("div", { class: "lp-row" }, eye,
        U.el("span", { class: "lp-name" }, `${name} — ${info.file}`));

      if (!info.exists) {
        row.append(U.el("span", { class: "lp-mini" }, "missing"));
      } else if (info.stale) {
        row.append(U.el("span", {
          class: "lp-mini stale-badge",
          title: info.shape_ok
            ? "The grid changed after this file was interpolated"
            : "Shape does not match the current grid - re-interpolate",
        }, "⚠ grid changed"));
      }

      row.append(U.miniBtn("interp", "Interpolate…", () => _interpolateWizard(name, info)));

      if (info.exists && !name.endsWith("_mask")) {
        row.append(U.miniBtn("copy", "Convert to sample data (for modification)", async () => {
          try {
            await Api.post("/api/domain/to_sample", { target: name });
            U.toast(`${name} converted to a sample layer`, "ok");
            _refresh();
          } catch (err) {
            U.toast(err.message, "error");
          }
        }));
      }
      list.append(row);
    }
    return list;
  }

  /* ================= download wizard ================= */

  function _downloadWizard() {
    const popup = Popup.open({ title: "Download data", width: 620 });
    let bounds = null;
    let availability = null;
    const selections = new Map();   // source -> Set(years)

    // --- area choice ---
    const useGrid = U.el("input", { type: "radio", name: "dl-area", id: "dl-grid", checked: "" });
    const useDraw = U.el("input", { type: "radio", name: "dl-area", id: "dl-draw" });
    const buffer = U.el("input", { type: "text", value: "500", style: "width:70px" });
    const drawBtn = U.el("button", { class: "ghost" }, "Draw area on map");
    const areaNote = U.el("div", { class: "muted", style: "font-size:12px" });

    const computeBounds = () => {
      if (useGrid.checked) {
        const p = GridTab.params();
        if (!p) { areaNote.textContent = "⚠ no grid yet — draw an area instead"; return null; }
        const t = p.rotation * Math.PI / 180;
        const ex = [Math.cos(t), Math.sin(t)], ey = [-Math.sin(t), Math.cos(t)];
        const c = (i, j) => [p.x0 + ex[0] * i * p.dx + ey[0] * j * p.dx,
          p.y0 + ex[1] * i * p.dx + ey[1] * j * p.dx];
        const ring = [c(0, 0), c(p.nx, 0), c(p.nx, p.ny), c(0, p.ny)];
        const xs = ring.map((q) => q[0]), ys = ring.map((q) => q[1]);
        const b = Number(buffer.value) || 0;
        return [Math.min(...xs) - b, Math.min(...ys) - b, Math.max(...xs) + b, Math.max(...ys) + b];
      }
      return bounds;
    };

    drawBtn.addEventListener("click", async () => {
      useDraw.checked = true;
      popup.hide();
      try {
        const obj = await Draw.polygon({ name: "download area" });
        const xs = obj.coords.map((c) => c[0]), ys = obj.coords.map((c) => c[1]);
        bounds = [Math.min(...xs), Math.min(...ys), Math.max(...xs), Math.max(...ys)];
        areaNote.textContent = `drawn area: ${U.fmtNum(bounds[2] - bounds[0], 5)} × ${U.fmtNum(bounds[3] - bounds[1], 5)} m`;
      } catch { /* cancelled */ }
      popup.show();
    });

    // --- availability ---
    const checkBtn = U.el("button", { class: "primary" }, "Check availability");
    const progress = U.el("div", { class: "muted", style: "font-size:12px" });
    const results = U.el("div");
    const dlAllBtn = U.el("button", { class: "primary", disabled: "" }, "Download selected");
    const estNote = U.el("span", { class: "muted", style: "margin-left:8px" });

    checkBtn.addEventListener("click", async () => {
      const area = computeBounds();
      if (!area) { U.toast("Define an area first", "error"); return; }
      checkBtn.disabled = true;
      try {
        const res = await Api.post("/api/domain/check", { bounds: area });
        availability = await Api.waitJob(res.job, (j) => { progress.textContent = j.message || ""; });
        progress.textContent = "";
        _renderAvailability(area);
      } catch (err) {
        U.toast(err.message, "error");
      } finally {
        checkBtn.disabled = false;
      }
    });

    function _updateEstimate() {
      let bytes = 0, count = 0;
      for (const [source, years] of selections) {
        const res = availability[source];
        for (const info of (res && res.years) || []) {
          if (years.has(info.year)) { bytes += info.est_bytes || 0; count += 1; }
        }
      }
      dlAllBtn.disabled = count === 0;
      estNote.textContent = count ? `${count} dataset(s), ~${U.fmtBytes(bytes)}` : "";
    }

    function _renderAvailability(area) {
      U.clear(results);
      selections.clear();
      for (const source of overview.sources) {
        if (source.id === "xyz") continue;
        const res = availability[source.id];
        if (!res) continue;
        results.append(U.el("div", { class: "fg-label", style: "margin-top:8px" }, source.title));
        if (res.error) {
          results.append(U.el("div", { class: "muted" }, `⚠ ${res.error}`));
          continue;
        }
        if (!res.available) {
          results.append(U.el("div", { class: "muted" }, res.notes || "not available here"));
          continue;
        }
        const selected = new Set();
        selections.set(source.id, selected);
        results.append(_yearChips(res.years, selected, _updateEstimate));
        if (res.notes) {
          results.append(U.el("div", { class: "muted", style: "font-size:11px" }, res.notes));
        }
      }
      results.append(U.el("div", { class: "muted", style: "font-size:11.5px;margin-top:6px" },
        "Select years by clicking, Shift+click for a range, or drag across chips."));
      dlAllBtn.onclick = () => _downloadAll(area);
    }

    async function _downloadAll(area) {
      dlAllBtn.disabled = true;
      try {
        let total = 0;
        for (const [source, years] of selections) {
          if (!years.size) continue;
          progress.textContent = `downloading ${source}…`;
          const res = await Api.post("/api/domain/download", {
            source, bounds: area, years: [...years],
          });
          const out = await Api.waitJob(res.job, (j) => {
            progress.textContent = `${source}: ${j.message || ""}`;
          });
          total += out.entries.length;
        }
        progress.textContent = "";
        U.toast(`Downloaded ${total} dataset(s)`, "ok");
        popup.close();
        _refresh();
      } catch (err) {
        U.toast(err.message, "error");
        dlAllBtn.disabled = false;
      }
    }

    popup.body.append(
      U.el("span", { class: "fg-label" }, "Area"),
      U.el("div", { class: "choice-row" }, useGrid,
        U.el("label", { for: "dl-grid" }, "Grid extent + buffer of"), buffer,
        U.el("span", { class: "muted" }, "m")),
      U.el("div", { class: "choice-row" }, useDraw,
        U.el("label", { for: "dl-draw" }, "Drawn area"),
        U.el("span", { class: "grow" }), drawBtn),
      areaNote,
      U.el("div", { style: "border-top:1px solid var(--border);margin:12px 0 10px" }),
      U.el("div", { class: "btn-row" }, checkBtn),
      progress,
      results,
      U.el("div", { class: "btn-row", style: "margin-top:10px;justify-content:flex-end;align-items:center" },
        estNote, dlAllBtn),
    );
  }

  /* year chips with click / shift+click / drag selection */
  function _yearChips(years, selected, onChange) {
    const wrap = U.el("div", { class: "year-chips" });
    const chips = [];
    let lastIndex = null;
    let dragging = false;
    let dragMode = true;   // select or deselect during drag

    const sync = () => {
      chips.forEach((chip, i) => chip.classList.toggle("on", selected.has(years[i].year)));
      onChange();
    };
    const setSel = (i, on) => {
      if (on) selected.add(years[i].year);
      else selected.delete(years[i].year);
    };

    years.forEach((info, i) => {
      const chip = U.el("button", {
        class: "year-chip",
        title: info.est_bytes ? U.fmtBytes(info.est_bytes) : "",
      }, String(info.year));
      chip.addEventListener("click", (ev) => {
        if (ev.shiftKey && lastIndex !== null) {
          const [a, b] = [Math.min(lastIndex, i), Math.max(lastIndex, i)];
          for (let k = a; k <= b; k += 1) setSel(k, true);
        } else {
          setSel(i, !selected.has(info.year));
        }
        lastIndex = i;
        sync();
      });
      chip.addEventListener("mousedown", (ev) => {
        if (ev.shiftKey) return;
        dragging = true;
        dragMode = !selected.has(info.year);
      });
      chip.addEventListener("mouseenter", () => {
        if (!dragging) return;
        setSel(i, dragMode);
        sync();
      });
      chips.push(chip);
      wrap.append(chip);
    });
    window.addEventListener("mouseup", () => { dragging = false; });
    return wrap;
  }

  /* ================= modify wizard (samples) ================= */

  function _modifyWizard(entry) {
    const popup = Popup.open({ title: `Modify — ${entry.label}`, width: 520 });

    const scopeAll = U.el("input", { type: "radio", name: "mod-scope", id: "ms-all", checked: "" });
    const scopePoly = U.el("input", { type: "radio", name: "mod-scope", id: "ms-poly" });
    const scopeIdx = U.el("input", { type: "radio", name: "mod-scope", id: "ms-idx" });

    const polySelect = U.el("select", {});
    const refreshPolys = () => {
      U.clear(polySelect);
      for (const obj of Objects.byKind("polygon")) {
        polySelect.append(U.el("option", { value: obj.id }, obj.name));
      }
      if (!polySelect.children.length) {
        polySelect.append(U.el("option", { value: "" }, "— none drawn yet —"));
      }
    };
    refreshPolys();
    const drawBtn = U.el("button", { class: "ghost" }, "Draw new");
    drawBtn.addEventListener("click", async () => {
      scopePoly.checked = true;
      popup.hide();
      try {
        const obj = await Draw.polygon({ name: "modify area" });
        refreshPolys();
        polySelect.value = obj.id;
      } catch { /* cancelled */ }
      popup.show();
    });

    const idxInputs = ["j0", "j1", "i0", "i1"].map((ph) =>
      U.el("input", { type: "text", placeholder: ph, style: "width:52px" }));

    const op = U.el("select", {},
      ...["set", "add", "subtract", "multiply", "min", "max"].map((o) =>
        U.el("option", { value: o }, o)));
    const value = U.el("input", { type: "text", placeholder: "value", style: "width:90px" });

    const saveOver = U.el("input", { type: "radio", name: "mod-save", id: "msv-over", checked: "" });
    const saveNew = U.el("input", { type: "radio", name: "mod-save", id: "msv-new" });
    const newName = U.el("input", { type: "text", placeholder: "new layer name" });
    newName.addEventListener("input", () => { if (newName.value) saveNew.checked = true; });

    const applyBtn = U.el("button", { class: "primary" }, "Apply");
    applyBtn.addEventListener("click", async () => {
      const body = { id: entry.id, op: op.value, value: value.value };
      if (scopePoly.checked) {
        if (!polySelect.value) { U.toast("Select or draw a polygon", "error"); return; }
        body.scope = { type: "polygon", polygon: polySelect.value };
      } else if (scopeIdx.checked) {
        body.scope = { type: "indices", indices: idxInputs.map((n) => Number(n.value || 0)) };
      } else {
        body.scope = { type: "all" };
      }
      if (saveNew.checked) {
        if (!newName.value.trim()) { U.toast("Enter a name for the new layer", "error"); return; }
        body.save_as = newName.value.trim();
      }
      try {
        applyBtn.disabled = true;
        const res = await Api.post("/api/domain/sample_modify", body);
        U.toast(`Modified ${res.cells} samples (${U.fmtNum(res.min)} … ${U.fmtNum(res.max)})`, "ok");
        popup.close();
        _refresh();
      } catch (err) {
        U.toast(err.message, "error");
        applyBtn.disabled = false;
      }
    });

    popup.body.append(
      U.el("span", { class: "fg-label" }, "Which samples"),
      U.el("div", { class: "choice-row" }, scopeAll, U.el("label", { for: "ms-all" }, "All")),
      U.el("div", { class: "choice-row" }, scopePoly, U.el("label", { for: "ms-poly" }, "Inside polygon"),
        polySelect, drawBtn),
      U.el("div", { class: "choice-row" }, scopeIdx, U.el("label", { for: "ms-idx" }, "Index range"),
        ...idxInputs),
      U.el("span", { class: "fg-label", style: "margin-top:8px" }, "Operation"),
      U.el("div", { class: "form-row" }, op, value),
      U.el("span", { class: "fg-label", style: "margin-top:8px" }, "Save"),
      U.el("div", { class: "choice-row" }, saveOver, U.el("label", { for: "msv-over" }, "Overwrite this layer")),
      U.el("div", { class: "choice-row" }, saveNew, U.el("label", { for: "msv-new" }, "Save as new layer"),
        U.el("span", { class: "grow" }), newName),
      U.el("div", { class: "btn-row", style: "justify-content:flex-end" }, applyBtn),
    );
  }

  /* ================= interpolate wizard ================= */

  function _interpolateWizard(target, info) {
    if (!overview.grid_available) {
      U.toast("Create a model grid first (Grid tab)", "error");
      return;
    }
    const popup = Popup.open({ title: `Interpolate → ${target} (${info.file})`, width: 540 });

    let order = overview.entries.map((e) => e.id);
    const checked = new Set();
    const listEl = U.el("div", { class: "layer-tree" });

    const renderList = () => {
      U.clear(listEl);
      if (!order.length) {
        listEl.append(U.el("div", { class: "muted" }, "No sample data yet — download or import first."));
        return;
      }
      order.forEach((id, idx) => {
        const entry = overview.entries.find((e) => e.id === id);
        if (!entry) return;
        const cb = U.el("input", { type: "checkbox" });
        cb.checked = checked.has(id);
        cb.addEventListener("change", () => {
          if (cb.checked) checked.add(id); else checked.delete(id);
        });
        const row = U.el("div", { class: "lp-row" }, cb,
          U.el("span", { class: "lp-name" }, entry.label || entry.path),
          U.el("span", { class: "lp-mini" }, entry.kind));
        if (idx > 0) {
          const up = U.el("span", { class: "lp-mini lp-btn" }, "↑");
          up.addEventListener("click", () => {
            [order[idx - 1], order[idx]] = [order[idx], order[idx - 1]];
            renderList();
          });
          row.append(up);
        }
        if (idx < order.length - 1) {
          const down = U.el("span", { class: "lp-mini lp-btn" }, "↓");
          down.addEventListener("click", () => {
            [order[idx + 1], order[idx]] = [order[idx], order[idx + 1]];
            renderList();
          });
          row.append(down);
        }
        listEl.append(row);
      });
    };
    renderList();

    const fill = U.el("input", { type: "text", placeholder: "e.g. -20 (optional)", style: "width:120px" });
    const progress = U.el("div", { class: "muted", style: "font-size:12px" });
    const runBtn = U.el("button", { class: "primary" }, "Interpolate & save");
    runBtn.addEventListener("click", async () => {
      const layers = order.filter((id) => checked.has(id));
      if (!layers.length) { U.toast("Select at least one sample layer", "error"); return; }
      const body = { target, layers };
      if (fill.value.trim() !== "") body.fill = Number(fill.value);
      try {
        runBtn.disabled = true;
        const res = await Api.post("/api/domain/interpolate", body);
        const out = await Api.waitJob(res.job, (j) => { progress.textContent = j.message || ""; });
        U.toast(`Wrote ${out.file} (${U.fmtNum(out.min)} … ${U.fmtNum(out.max)})`, "ok");
        const cfg = await Api.get("/api/config");
        App.state.config = cfg.values;
        App.emit("config-changed", overview.targets[target].config_key);
        popup.close();
        _refresh();
      } catch (err) {
        U.toast(err.message, "error");
        runBtn.disabled = false;
      }
    });

    popup.body.append(
      U.el("span", { class: "fg-label" }, "Sample layers (top = highest priority)"),
      listEl,
      U.el("div", { class: "form-row", style: "margin-top:8px" },
        U.el("label", {}, "Fill remaining cells"), fill),
      U.el("div", { class: "btn-row" }, runBtn),
      progress,
    );
  }

  return { init };
})();
