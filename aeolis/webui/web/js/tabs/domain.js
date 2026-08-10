/* Domain tab: bathymetry / vegetation / ne-layer data.
 *
 * Two sections, both rendered as one card per dataset:
 *  - Sample data: downloaded/imported datasets with visibility, rename,
 *    drag-to-reorder, duplicate, remove and a Modify… popup.
 *  - Interpolated data: the model .grd files with visibility, an
 *    Interpolate… popup and a duplicate-to-samples action. Staleness vs
 *    the current grid shows as an orange warning icon (hover for the
 *    reason); a stale file cannot be displayed.
 * Downloading happens in a popup wizard with a real progress bar.
 */
"use strict";

const DomainTab = (() => {

  let overview = null;
  const loading = new Set();      // layer ids currently loading on the map
  const selectedIds = new Set();  // multi-selected sample cards
  let lastClickedId = null;       // shift+click range anchor
  const revealedTargets = new Set(); // optional domain files the user opted to add
  const hiddenTargets = new Set();   // optional files the user chose to hide from the list

  function init() {
    Tabs.register("domain", { enter: _refresh });
    App.on("project", _refresh);
    App.on("layer-loading", ({ id, busy }) => {
      if (busy) loading.add(id); else loading.delete(id);
      if (overview) _build();
    });
    App.on("layer-visibility", () => { if (overview) _build(); });
    // clicking anywhere outside a card / toolbar clears the selection
    U.deselectOnOutside(() => selectedIds.size,
      () => { selectedIds.clear(); if (overview) _build(); });
  }

  async function _refresh() {
    if (!App.state.project) return;
    try {
      overview = await Api.get("/api/domain");
    } catch (err) {
      U.toast(err.message, "error");
      return;
    }
    // prune selection of removed entries
    const known = new Set(overview.entries.map((e) => e.id));
    for (const id of [...selectedIds]) {
      if (!known.has(id)) selectedIds.delete(id);
    }
    _registerLayers();
    _build();
  }

  /* Re-load a layer's map display after its data changed on disk. */
  function _reloadLayer(layerId) {
    const layer = Layers.get(layerId);
    if (layer && layer.visible) App.emit("layer-visibility", layer);
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
        stale: info.stale, shape_ok: info.shape_ok,
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
    const removeBtn = U.tbtn("trash", "Remove", {
      title: selectedIds.size
        ? `Remove ${selectedIds.size} selected dataset(s)…`
        : "Select cards first (click, Ctrl+click, Shift+click for a range)",
      onclick: _removeSelected,
    });
    removeBtn.disabled = !selectedIds.size;
    samples.body.append(
      U.el("div", { class: "tbtn-row" },
        U.tbtn("download", "Download", { primary: true, title: "Download data from Dutch coastal sources", onclick: _downloadWizard }),
        U.tbtn("upload", "Import", { title: "Import a *.tif topography raster or *.xyz sample file", onclick: _importFile }),
        U.tbtn("link", "Link", { title: "Reference raw data from another project (no copy / re-download)", onclick: _linkExternal }),
        removeBtn),
      _sampleList(),
    );
    if (overview.entries.length > 1) {
      samples.body.append(U.el("div", { class: "muted", style: "font-size:11px;margin-top:4px" },
        "Click to select, Ctrl+click to add, Shift+click for a range — selected cards move and are removed together."));
    }
    panel.append(samples.wrap);

    // --- interpolated data ---
    const targetCount = Object.values(overview.targets).filter((t) => t.exists).length;
    const targets = U.section("Interpolated data (.grd)", { count: targetCount });
    targets.body.append(_targetList());
    const addOptional = _addOptionalControl();
    if (addOptional) targets.body.append(addOptional);
    if (!overview.grid_available) {
      targets.body.append(U.el("div", { class: "muted", style: "margin-top:8px" },
        "⚠ No model grid yet — create one in the Grid tab before interpolating."));
    }
    panel.append(targets.wrap);
  }

  /* ---- generic card helpers ---- */

  function _eyeOrSpinner(layerId, layer, { disabled = false, disabledTitle = "" } = {}) {
    if (loading.has(layerId)) {
      return U.el("span", { class: "spin", title: "Loading…" });
    }
    // a stale layer cannot be turned ON, but a visible one can
    // still be hidden
    if (layer.visible) disabled = false;
    const eye = U.el("span", {
      class: `eye ${layer.visible ? "" : "off"} ${disabled ? "disabled" : ""}`,
      title: disabled ? disabledTitle : "Show/hide on the map",
    }, "👁");
    if (!disabled) {
      eye.addEventListener("click", () => {
        layer.visible = !layer.visible;
        App.emit("layer-visibility", layer);
        _build();
      });
    }
    return eye;
  }

  function _warnIcon(title) {
    return U.el("span", { class: "warn-icon", title }, "⚠");
  }

  /* Make the cards of a list draggable; onDrop(fromIdx, toIdx).
   * Delegates to the shared, de-lagged sortable (drag starts from the
   * grip so dblclick/selection on the card body keeps working). */
  function _wireCardDrag(list, onDrop) {
    U.wireSortable(list, onDrop);
  }

  /* ---- sample list ---- */

  function _sampleList() {
    const list = U.el("div", { class: "obj-list" });
    if (!overview.entries.length) {
      list.append(U.el("div", { class: "muted" }, "Nothing downloaded or imported yet."));
      return list;
    }
    overview.entries.forEach((entry, idx) => {
      const layerId = `raw-${entry.id}`;
      const layer = Layers.get(layerId) || { visible: false };

      // always-editable name field; Save commits the label AND renames the
      // underlying .npz file on disk (enabled only when the text changed)
      const stored = entry.label || entry.path;
      const nameInput = U.el("input", {
        class: "lp-name-edit", type: "text", value: stored,
        title: "Edit the sample name — Save also renames the file in gui/rawdata",
      });
      const saveName = U.miniBtn("save", "Save name (renames the file on disk)", async () => {
        const nm = nameInput.value.trim();
        if (!nm || nm === stored) return;
        try {
          await Api.post("/api/domain/sample_rename", { id: entry.id, name: nm, rename_file: true });
          U.toast("Renamed", "ok");
          _refresh();
        } catch (err) { U.toast(err.message, "error"); }
      });
      const syncSave = () => {
        const nm = nameInput.value.trim();
        saveName.disabled = !nm || nm === stored;
      };
      nameInput.addEventListener("input", syncSave);
      nameInput.addEventListener("click", (ev) => ev.stopPropagation());
      nameInput.addEventListener("keydown", (ev) => {
        if (ev.key === "Enter") { ev.preventDefault(); saveName.click(); }
        if (ev.key === "Escape") { nameInput.value = stored; syncSave(); }
      });
      syncSave();

      const card = U.el("div", {
        class: `obj-card ${selectedIds.has(entry.id) ? "selected" : ""}`,
        dataset: { idx, eid: entry.id },
      },
        U.el("span", { class: "drag-grip", draggable: "true", title: "Drag to reorder (top = highest priority); selected cards move together" }, "⠿"),
        _eyeOrSpinner(layerId, layer),
        nameInput,
        saveName,
        U.el("span", { class: "lp-mini" }, entry.source),
        _bandSelect(entry),
        U.miniBtn("copy", "Duplicate", async () => {
          await Api.post("/api/domain/sample_duplicate", { id: entry.id });
          _refresh();
        }),
        U.miniBtn("modify", "Modify…", () => _modifyWizard(entry)));
      // click = select (Ctrl toggles, Shift selects a range)
      card.addEventListener("click", (ev) => {
        if (ev.target.closest("button, .eye, .spin, input, .drag-grip")) return;
        const ids = overview.entries.map((e) => e.id);
        if (ev.shiftKey && lastClickedId && ids.includes(lastClickedId)) {
          const a = ids.indexOf(lastClickedId), b = ids.indexOf(entry.id);
          for (let k = Math.min(a, b); k <= Math.max(a, b); k += 1) selectedIds.add(ids[k]);
        } else if (ev.ctrlKey || ev.metaKey) {
          if (selectedIds.has(entry.id)) selectedIds.delete(entry.id);
          else selectedIds.add(entry.id);
        } else if (selectedIds.size === 1 && selectedIds.has(entry.id)) {
          selectedIds.clear();     // click the only selected card = deselect
        } else {
          selectedIds.clear();
          selectedIds.add(entry.id);
        }
        lastClickedId = entry.id;
        _build();
      });
      list.append(card);
    });
    _wireCardDrag(list, async (from, to) => {
      const ids = overview.entries.map((e) => e.id);
      const dragId = ids[from];
      // selected cards move as one block when a selected card is dragged
      const group = (selectedIds.has(dragId) && selectedIds.size > 1)
        ? ids.filter((id) => selectedIds.has(id))
        : [dragId];
      // anchor = the element that ends up at the drop position when
      // only the dragged card is removed (single-move semantics)
      const single = ids.filter((id) => id !== dragId);
      let anchor = to < single.length ? single[to] : null;
      const rest = ids.filter((id) => !group.includes(id));
      // if the anchor is part of the moving group, insert after the
      // nearest following card that stays put
      while (anchor && group.includes(anchor)) {
        const next = single.indexOf(anchor) + 1;
        anchor = next < single.length ? single[next] : null;
      }
      const at = anchor ? rest.indexOf(anchor) : rest.length;
      rest.splice(at, 0, ...group);
      await Api.post("/api/domain/sample_order", { ids: rest });
      await _refresh();
      App.emit("layer-order");
    });
    return list;
  }

  /* Reference raw data that already lives in another project instead of
   * re-downloading it: pick that project (or its gui/rawdata), choose
   * datasets, and add manifest entries that point at the files in place. */
  async function _linkExternal() {
    const path = await Api.pickFolder({
      title: "Select another AeoLiS project (or its gui/rawdata folder)",
    }).catch(() => null);
    if (!path) return;
    let scan;
    try {
      scan = await Api.post("/api/domain/scan_external", { path });
    } catch (err) { U.toast(err.message, "error"); return; }

    const linkable = scan.entries.filter((e) => e.exists && !e.is_local);
    if (!linkable.length) {
      U.toast(scan.entries.length ? "Nothing new to link here" : "No datasets found there", "error");
      return;
    }

    const popup = Popup.open({ title: "Link external raw data", width: 540 });
    const checks = new Map();
    const list = U.el("div", { class: "obj-list", style: "max-height:320px;overflow-y:auto" });
    for (const e of linkable) {
      const cb = U.el("input", { type: "checkbox", checked: "" });
      checks.set(e.src_id, { cb, e });
      const meta = [e.source, e.date || (e.year != null ? String(e.year) : null), U.fmtBytes(e.size)]
        .filter((v) => v != null && v !== "").join(" · ");
      list.append(U.el("label", { class: "obj-card", style: "cursor:pointer" },
        cb,
        U.el("span", { class: "lp-name" }, e.label),
        U.el("span", { class: "lp-mini" }, meta)));
    }
    const okBtn = U.el("button", { class: "primary" }, "Link selected");
    okBtn.addEventListener("click", async () => {
      const ids = [...checks.values()].filter(({ cb }) => cb.checked).map(({ e }) => e.src_id);
      if (!ids.length) { U.toast("Select at least one dataset", "error"); return; }
      okBtn.disabled = true;
      try {
        const res = await Api.post("/api/domain/link_external", { path, ids });
        U.toast(`Linked ${res.linked} dataset(s)`, "ok");
        popup.close();
        _refresh();
      } catch (err) { U.toast(err.message, "error"); okBtn.disabled = false; }
    });
    const cancelBtn = U.el("button", { class: "ghost" }, "Cancel");
    cancelBtn.addEventListener("click", popup.close);
    popup.body.append(
      U.el("div", { class: "muted", style: "font-size:12px;margin-bottom:6px" },
        `From ${scan.root} — referenced in place (no copy, no re-download).`),
      list,
      U.el("div", { class: "btn-row", style: "justify-content:flex-end;margin-top:8px" }, cancelBtn, okBtn),
    );
  }

  /* band picker for multi-channel rasters (RGB/CIR tiffs): which band is
   * displayed and used in interpolation. Returns null for single-band. */
  function _bandSelect(entry) {
    if (!entry.bands || entry.bands <= 1) return null;
    const sel = U.el("select", { class: "cmap-mini", title: "Which band (channel) to display & use" });
    for (let b = 1; b <= entry.bands; b += 1) {
      sel.append(U.el("option", { value: b, selected: b === (entry.band || 1) ? "" : null }, `b${b}`));
    }
    sel.addEventListener("click", (ev) => ev.stopPropagation());
    sel.addEventListener("change", async () => {
      try {
        await Api.post("/api/domain/sample_band", { id: entry.id, band: Number(sel.value) });
        await _refresh();
        _reloadLayer(`raw-${entry.id}`);
      } catch (err) { U.toast(err.message, "error"); }
    });
    return sel;
  }

  function _removeSelected() {
    const entries = overview.entries.filter((e) => selectedIds.has(e.id));
    if (!entries.length) return;
    const popup = Popup.open({
      title: `Remove ${entries.length} sample dataset(s)`, width: 460 });
    const delFile = U.el("input", { type: "checkbox", id: "del-file", checked: "" });
    const okBtn = U.el("button", { class: "danger" }, `Remove ${entries.length}`);
    okBtn.addEventListener("click", async () => {
      try {
        okBtn.disabled = true;
        for (const entry of entries) {
          await Api.post("/api/domain/forget", { id: entry.id, delete_file: delFile.checked });
          FieldLayer.remove(`field-raw-${entry.id}`);
          MapView.removeLayerAndSource(`pts-raw-${entry.id}`);
          Layers.unregister(`raw-${entry.id}`);
          selectedIds.delete(entry.id);
        }
        popup.close();
        _refresh();
      } catch (err) {
        U.toast(err.message, "error");
        okBtn.disabled = false;
      }
    });
    const cancelBtn = U.el("button", { class: "ghost" }, "Cancel");
    cancelBtn.addEventListener("click", popup.close);
    popup.body.append(
      U.el("div", { style: "font-size:13px" }, "Remove these sample layers from the project?"),
      U.el("div", { class: "obj-list", style: "margin:8px 0;max-height:180px;overflow-y:auto" },
        ...entries.map((entry) => U.el("div", { class: "obj-card" },
          U.el("span", { class: "lp-name" }, entry.label || entry.path),
          U.el("span", { class: "lp-mini" }, entry.source)))),
      U.el("div", { class: "choice-row", style: "margin-top:8px" },
        delFile, U.el("label", { for: "del-file" }, "also delete the files from disk")),
      U.el("div", { class: "btn-row", style: "justify-content:flex-end" }, cancelBtn, okBtn),
    );
  }

  async function _importFile() {
    const path = await Api.pickFile({
      title: "Import topography or samples",
      patterns: [
        ["Topography & samples", "*.tif;*.tiff;*.xyz;*.txt;*.csv"],
        ["GeoTIFF topography", "*.tif;*.tiff"],
        ["Sample files", "*.xyz;*.txt;*.csv"],
        ["All files", "*.*"],
      ],
    }).catch(() => null);
    if (!path) return;
    const isTiff = /\.tiff?$/i.test(path);
    try {
      await Api.post(isTiff ? "/api/domain/import_tiff" : "/api/domain/import_xyz", { path });
      U.toast(isTiff ? "Topography imported" : "Samples imported", "ok");
      _refresh();
    } catch (err) {
      U.toast(err.message, "error");
    }
  }

  /* ---- interpolated targets ---- */

  function _targetOrder() {
    const names = Object.keys(overview.targets);
    const saved = App.state.ui.targetOrder || [];
    return [...saved.filter((n) => names.includes(n)),
      ...names.filter((n) => !saved.includes(n))];
  }

  function _targetList() {
    const list = U.el("div", { class: "obj-list" });
    // hide opt-in files the user hasn't added yet (offered via "Add…")
    const order = _targetOrder().filter((name) => {
      const info = overview.targets[name];
      return info && (!info.hidden || revealedTargets.has(name)) && !hiddenTargets.has(name);
    });
    order.forEach((name, idx) => {
      const info = overview.targets[name];
      const layerId = `domain-${name}`;
      const layer = Layers.get(layerId) || { visible: false, id: layerId };

      const staleTitle = info.shape_ok
        ? "The grid changed after this file was interpolated — re-interpolate before using it"
        : "Shape does not match the current grid — re-interpolate";

      // renderable when a saved file exists OR there is an unsaved draft
      const renderable = info.exists || info.has_draft;
      const eye = renderable
        ? _eyeOrSpinner(layerId, layer, {
          disabled: info.stale,
          disabledTitle: staleTitle,
        })
        : U.el("span", { class: "eye off disabled", title: "Nothing to show yet — interpolate first" }, "👁");

      // not required by the current config (e.g. veg under the grass
      // method, or an opt-in mask) → de-emphasise, keep it usable
      const card = U.el("div", {
        class: `obj-card ${info.needed ? "" : "not-needed"}`,
        dataset: { idx },
      },
        U.el("span", { class: "drag-grip", draggable: "true", title: "Drag to reorder" }, "⠿"),
        eye,
        U.el("span", { class: "lp-name" }, `${name} — ${info.file}`));

      if (!info.needed) {
        card.append(U.el("span", {
          class: "opt-badge",
          title: info.note || "not required by the current configuration",
        }, info.optional ? "optional" : "not needed"));
      }
      // clear state indicator: unsaved draft > missing file > not created
      if (info.has_draft) {
        card.append(U.el("span", { class: "draft-badge", title: "Interpolated but not saved yet — click Save to write the file" }, "unsaved draft"));
      } else if (!info.exists) {
        card.append(info.configured
          ? _warnIcon(`missing: ${info.file} — the configured file was not found`)
          : U.el("span", { class: "lp-mini" }, "not created yet"));
      } else if (info.stale) {
        card.append(_warnIcon(staleTitle));
      }

      const actions = U.el("span", { class: "obj-actions" });
      actions.append(U.miniBtn("interp", "Interpolate…", () => _interpolateWizard(name, info)));
      actions.append(U.miniBtn("wand", "New grid from a constant value (e.g. a mask of 1s)", () => _targetConstant(name)));
      if (info.exists || info.has_draft) {
        actions.append(U.miniBtn("modify", "Edit values (set/add/… over a polygon or index box)", () => _targetModifyWizard(name)));
      }
      if (info.has_draft) {
        actions.append(U.miniBtn("save", "Save the draft to the configured file", () => _saveTargetDraft(name)));
      }
      actions.append(U.miniBtn("open", "Load an existing .grd on top (as a draft)", () => _loadTarget(name)));
      if (info.exists || info.has_draft) {
        actions.append(U.miniBtn("saveas", "Save to a chosen location…", () => _saveTargetAs(name, info)));
      }
      if (info.exists && !name.endsWith("_mask")) {
        actions.append(U.miniBtn("copy", "Duplicate to sample data (for modification)", async () => {
          try {
            await Api.post("/api/domain/to_sample", { target: name });
            U.toast(`${name} duplicated to a sample layer`, "ok");
            _refresh();
          } catch (err) {
            U.toast(err.message, "error");
          }
        }));
      }
      // an added-but-empty optional file can be dismissed again
      if (info.hidden && revealedTargets.has(name) && !info.exists) {
        actions.append(U.miniBtn("trash", "Remove from the list (no file is written)", () => {
          revealedTargets.delete(name);
          _build();
        }));
      } else if (info.optional) {
        // optional masks shown by default (configured/exist) can be hidden
        actions.append(U.miniBtn("trash", "Hide from this list (does not delete the file)", () => {
          hiddenTargets.add(name);
          _build();
        }));
      }
      card.append(actions);
      list.append(card);
    });
    _wireCardDrag(list, (from, to) => {
      const order2 = _targetOrder();
      const [moved] = order2.splice(from, 1);
      order2.splice(to, 0, moved);
      App.state.ui.targetOrder = order2;
      App.touchUi();
      _build();
    });
    return list;
  }

  /* Load an existing .grd on top of a target as an unsaved draft. */
  async function _loadTarget(name) {
    const path = await Api.pickFile({
      title: `Load a .grd file for ${name}`,
      patterns: [["Grid files", "*.grd"], ["All files", "*.*"]],
      initial: App.state.project ? App.state.project.root : "",
    }).catch((err) => { U.toast(err.message, "error"); return null; });
    if (!path) return;
    try {
      const res = await Api.post("/api/domain/target_load", { target: name, path });
      U.toast(`Loaded as a draft${res.shape_ok ? "" : " (shape differs from grid)"} — Save to commit`,
        res.shape_ok ? "ok" : "error");
      await _refresh();
      _reloadLayer(`domain-${name}`);
    } catch (err) { U.toast(err.message, "error"); }
  }

  /* Commit the current draft to the target's configured file. */
  async function _saveTargetDraft(name) {
    try {
      const res = await Api.post("/api/domain/target_save", { target: name });
      U.toast(`Saved ${res.file}`, "ok");
      const cfg = await Api.get("/api/config");
      App.state.config = cfg.values;
      App.emit("config-changed", overview.targets[name].config_key);
      await _refresh();
    } catch (err) { U.toast(err.message, "error"); }
  }

  /* Save the draft (or the saved file) to a chosen location + repoint config. */
  async function _saveTargetAs(name, info) {
    const path = await Api.pickFile({
      title: `Save ${name} as…`,
      save: true,
      patterns: [["Grid files", "*.grd"], ["All files", "*.*"]],
      initial: App.state.project ? App.state.project.root : "",
      filename: info.file || `${name}.grd`,
    }).catch((err) => { U.toast(err.message, "error"); return null; });
    if (!path) return;
    try {
      const res = await Api.post("/api/domain/target_save_as", { target: name, path });
      U.toast(`Saved ${res.file}`, "ok");
      const cfg = await Api.get("/api/config");
      App.state.config = cfg.values;
      App.emit("config-changed", overview.targets[name].config_key);
      await _refresh();
    } catch (err) { U.toast(err.message, "error"); }
  }

  /* "Add optional file…" — reveal an opt-in domain file (mask etc.) so it
   * can be interpolated; most domain files need no input and stay hidden. */
  function _addOptionalControl() {
    const hidden = Object.entries(overview.targets)
      .filter(([name, info]) => info.optional
        && (info.hidden || hiddenTargets.has(name))
        && !(revealedTargets.has(name) && !hiddenTargets.has(name)));
    if (!hidden.length) return null;
    const sel = U.el("select", { style: "flex:1;min-width:0" },
      U.el("option", { value: "" }, "— add an optional domain file —"),
      ...hidden.map(([name, info]) =>
        U.el("option", { value: name }, `${name} (${info.file})`)));
    sel.addEventListener("change", () => {
      if (!sel.value) return;
      hiddenTargets.delete(sel.value);
      revealedTargets.add(sel.value);
      _build();
    });
    return U.el("div", { class: "form-row", style: "margin-top:8px" }, sel);
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
    const progress = U.progressBar();
    const results = U.el("div");
    const dlAllBtn = U.el("button", { class: "primary", disabled: "" }, "Download selected");
    const estNote = U.el("span", { class: "muted", style: "margin-left:8px" });

    checkBtn.addEventListener("click", async () => {
      const area = computeBounds();
      if (!area) { U.toast("Define an area first", "error"); return; }
      checkBtn.disabled = true;
      progress.start("checking availability…");
      try {
        const res = await Api.post("/api/domain/check", { bounds: area });
        availability = await Api.waitJob(res.job, (j) => progress.update(j));
        progress.done();
        _renderAvailability(area);
      } catch (err) {
        progress.done();
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
        if (source.id === "xyz" || source.id === "tiff") continue;
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
      const jobs = [];
      for (const [source, years] of selections) {
        if (years.size) jobs.push([source, [...years]]);
      }
      try {
        let total = 0;
        for (let k = 0; k < jobs.length; k += 1) {
          const [source, years] = jobs[k];
          progress.start(`downloading ${source}…`);
          const res = await Api.post("/api/domain/download", {
            source, bounds: area, years,
          });
          const out = await Api.waitJob(res.job, (j) => {
            // overall bar: finished sources + progress within this one
            const frac = j.progress >= 0 ? j.progress : null;
            progress.update({
              progress: frac === null ? -1 : (k + frac) / jobs.length,
              message: `${source}: ${j.message || ""}`,
            });
          });
          total += out.entries.length;
        }
        progress.done();
        U.toast(`Downloaded ${total} dataset(s)`, "ok");
        popup.close();
        _refresh();
      } catch (err) {
        progress.done();
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
      results,
      U.el("div", { class: "btn-row", style: "margin-top:10px;justify-content:flex-end;align-items:center" },
        estNote, dlAllBtn),
      progress.el,
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
    // self-cleaning: drop the listener once this chip row is gone
    const onUp = () => {
      dragging = false;
      if (!wrap.isConnected) window.removeEventListener("mouseup", onUp);
    };
    window.addEventListener("mouseup", onUp);
    return wrap;
  }

  /* ================= modify wizard (samples) ================= */

  const MODIFY_OPS = [
    ["set", "set to value"],
    ["add", "add value"],
    ["subtract", "subtract value"],
    ["multiply", "multiply by value"],
    ["clip_max", "cap at max (values above → value)"],
    ["clip_min", "cap at min (values below → value)"],
  ];

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
      ...MODIFY_OPS.map(([value, label]) => U.el("option", { value }, label)));
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
        await _refresh();
        // the layer data changed on disk: refresh the map display too
        if (!body.save_as) _reloadLayer(`raw-${entry.id}`);
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

    // ---- band math & classify (rasters): NDVI-style expressions ----
    if (entry.kind === "raster" || entry.kind === "raster_nc") {
      _appendBandMath(popup, entry);
    }
  }

  /* Band math + threshold classification producing a NEW dataset, e.g.
   * NDVI = (b4-b1)/(b4+b1) then "NDVI > 0.3 → 1 else 0" for rho_veg. */
  function _appendBandMath(popup, entry) {
    const bands = entry.bands || 1;
    const expr = U.el("input", {
      type: "text", style: "flex:1",
      placeholder: bands > 1 ? "e.g. (b4-b1)/(b4+b1)" : "e.g. b1",
      value: bands >= 4 ? "(b4-b1)/(b4+b1)" : "b1",
    });
    const thrOn = U.el("input", { type: "checkbox", id: "bm-thr", checked: "" });
    const thrOp = U.el("select", {},
      U.el("option", { value: ">" }, ">"), U.el("option", { value: "<" }, "<"));
    const thrX = U.el("input", { type: "text", value: "0.3", style: "width:64px" });
    const thrThen = U.el("input", { type: "text", value: "1", style: "width:52px" });
    const thrElse = U.el("input", { type: "text", value: "0", style: "width:52px" });
    const outName = U.el("input", { type: "text", placeholder: "new dataset name", style: "flex:1" });

    const runBtn = U.el("button", { class: "primary" }, "Create dataset");
    runBtn.addEventListener("click", async () => {
      const name = outName.value.trim();
      if (!name) { U.toast("Name the new dataset", "error"); return; }
      const body = { id: entry.id, expr: expr.value, save_as: name };
      if (thrOn.checked) {
        body.threshold = { op: thrOp.value, x: Number(thrX.value),
          then: Number(thrThen.value), else: Number(thrElse.value) };
      }
      try {
        runBtn.disabled = true;
        const res = await Api.post("/api/domain/raster_derive", body);
        U.toast(`"${name}" created (${U.fmtNum(res.min)} … ${U.fmtNum(res.max)})`, "ok");
        popup.close();
        _refresh();
      } catch (err) { U.toast(err.message, "error"); runBtn.disabled = false; }
    });

    popup.body.append(
      U.el("hr", { style: "margin:12px 0;border:none;border-top:1px solid var(--border)" }),
      U.el("span", { class: "fg-label" }, "Band math & classify → new dataset"),
      U.el("div", { class: "muted", style: "font-size:11.5px" },
        `Available bands: ${Array.from({ length: bands }, (_, i) => `b${i + 1}`).join(", ")}`
        + " — e.g. NDVI = (NIR−Red)/(NIR+Red)."),
      U.el("div", { class: "form-row" }, U.el("label", {}, "Expression"), expr),
      U.el("div", { class: "choice-row" }, thrOn,
        U.el("label", { for: "bm-thr" }, "classify: if result"), thrOp, thrX,
        U.el("span", {}, "→"), thrThen, U.el("span", {}, "else"), thrElse),
      U.el("div", { class: "form-row" }, U.el("label", {}, "Save as"), outName, runBtn),
    );
  }

  /* ================= interpolated-target editing ================= */

  /* New draft filled with one constant value — e.g. a mask of all 1s
   * that polygon edits then carve 0s into. */
  function _targetConstant(target) {
    const popup = Popup.open({ title: `New ${target} grid — constant value`, width: 380 });
    const value = U.el("input", { type: "text", value: "1", style: "width:90px" });
    const okBtn = U.el("button", { class: "primary" }, "Create draft");
    okBtn.addEventListener("click", async () => {
      const v = Number(value.value);
      if (!Number.isFinite(v)) { U.toast("Enter a numeric value", "error"); return; }
      try {
        okBtn.disabled = true;
        await Api.post("/api/domain/target_constant", { target, value: v });
        U.toast(`${target} draft created (all ${v}) — edit it, then Save`, "ok");
        popup.close();
        await _refresh();
        _reloadLayer(`domain-${target}`);
      } catch (err) { U.toast(err.message, "error"); okBtn.disabled = false; }
    });
    const cancelBtn = U.el("button", { class: "ghost" }, "Cancel");
    cancelBtn.addEventListener("click", popup.close);
    popup.body.append(
      U.el("div", { class: "form-row" }, U.el("label", {}, "Fill value"), value),
      U.el("div", { class: "muted", style: "font-size:11.5px" },
        "Creates an unsaved draft on the model grid (e.g. 1 everywhere for a mask). "
        + "Use Edit values to set areas, then Save to write the .grd."),
      U.el("div", { class: "btn-row", style: "justify-content:flex-end;margin-top:8px" }, cancelBtn, okBtn),
    );
  }

  /* Edit the target draft (or saved .grd) — same scope options as the
   * sample modify wizard, but writing to the target draft. */
  function _targetModifyWizard(target) {
    const popup = Popup.open({ title: `Edit values — ${target}`, width: 520 });

    const scopeAll = U.el("input", { type: "radio", name: "tm-scope", id: "tms-all", checked: "" });
    const scopePoly = U.el("input", { type: "radio", name: "tm-scope", id: "tms-poly" });
    const scopeIdx = U.el("input", { type: "radio", name: "tm-scope", id: "tms-idx" });

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
        const obj = await Draw.polygon({ name: `${target} area` });
        refreshPolys();
        polySelect.value = obj.id;
      } catch { /* cancelled */ }
      popup.show();
    });

    const idxInputs = ["j0", "j1", "i0", "i1"].map((ph) =>
      U.el("input", { type: "text", placeholder: ph, style: "width:52px" }));

    const op = U.el("select", {},
      ...MODIFY_OPS.map(([value, label]) => U.el("option", { value }, label)));
    const value = U.el("input", { type: "text", placeholder: "value", style: "width:90px" });

    const applyBtn = U.el("button", { class: "primary" }, "Apply to draft");
    applyBtn.addEventListener("click", async () => {
      const body = { target, op: op.value, value: value.value };
      if (scopePoly.checked) {
        if (!polySelect.value) { U.toast("Select or draw a polygon", "error"); return; }
        body.scope = { type: "polygon", polygon: polySelect.value };
      } else if (scopeIdx.checked) {
        body.scope = { type: "indices", indices: idxInputs.map((n) => Number(n.value || 0)) };
      } else {
        body.scope = { type: "all" };
      }
      try {
        applyBtn.disabled = true;
        const res = await Api.post("/api/domain/target_modify", body);
        U.toast(`Edited ${res.cells} cells (${U.fmtNum(res.min)} … ${U.fmtNum(res.max)}) — draft updated`, "ok");
        popup.close();
        await _refresh();
        _reloadLayer(`domain-${target}`);
      } catch (err) {
        U.toast(err.message, "error");
        applyBtn.disabled = false;
      }
    });

    popup.body.append(
      U.el("span", { class: "fg-label" }, "Which cells"),
      U.el("div", { class: "choice-row" }, scopeAll, U.el("label", { for: "tms-all" }, "All")),
      U.el("div", { class: "choice-row" }, scopePoly, U.el("label", { for: "tms-poly" }, "Inside polygon"),
        polySelect, drawBtn),
      U.el("div", { class: "choice-row" }, scopeIdx, U.el("label", { for: "tms-idx" }, "Index range"),
        ...idxInputs),
      U.el("span", { class: "fg-label", style: "margin-top:8px" }, "Operation"),
      U.el("div", { class: "form-row" }, op, value),
      U.el("div", { class: "muted", style: "font-size:11.5px" },
        "Edits the unsaved draft (shown on the map). Save the card afterwards to write the .grd."),
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
    const listEl = U.el("div", { class: "obj-list" });

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
        listEl.append(U.el("div", {
          class: "obj-card", dataset: { idx },
        },
          U.el("span", { class: "drag-grip", draggable: "true", title: "Drag to reorder (top = highest priority)" }, "⠿"),
          cb,
          U.el("span", { class: "lp-name" }, entry.label || entry.path),
          U.el("span", { class: "lp-mini" }, entry.kind)));
      });
    };
    renderList();
    _wireCardDrag(listEl, (from, to) => {
      const [moved] = order.splice(from, 1);
      order.splice(to, 0, moved);
      renderList();
    });

    const extrap = U.el("input", { type: "checkbox", id: "interp-extrap" });
    const fill = U.el("input", { type: "text", placeholder: "e.g. -20 (optional)", style: "width:120px" });
    const progress = U.progressBar();
    const runBtn = U.el("button", { class: "primary" }, "Interpolate");
    runBtn.addEventListener("click", async () => {
      const layers = order.filter((id) => checked.has(id));
      if (!layers.length) { U.toast("Select at least one sample layer", "error"); return; }
      const body = { target, layers };
      if (extrap.checked) body.extrapolate = true;
      if (fill.value.trim() !== "") body.fill = Number(fill.value);
      try {
        runBtn.disabled = true;
        progress.start("interpolating…");
        const res = await Api.post("/api/domain/interpolate", body);
        const out = await Api.waitJob(res.job, (j) => progress.update(j));
        progress.done();
        U.toast(`Draft ready (${U.fmtNum(out.min)} … ${U.fmtNum(out.max)}) — Save to write the file`, "ok");
        popup.close();
        await _refresh();
        _reloadLayer(`domain-${target}`);   // previews the draft via gridfield
      } catch (err) {
        progress.done();
        U.toast(err.message, "error");
        runBtn.disabled = false;
      }
    });

    popup.body.append(
      U.el("span", { class: "fg-label" }, "Sample layers (top = highest priority)"),
      listEl,
      U.el("span", { class: "fg-label", style: "margin-top:10px" }, "Cells without data"),
      U.el("div", { class: "choice-row" }, extrap,
        U.el("label", { for: "interp-extrap" }, "Extrapolate to nearest neighbour"),
        U.el("span", { class: "muted", style: "font-size:11px" }, "(fill gaps from the closest sample)")),
      U.el("div", { class: "form-row" },
        U.el("label", {}, "…or fill with a constant"), fill),
      U.el("div", { class: "btn-row" }, runBtn),
      progress.el,
    );
  }

  return { init };
})();
