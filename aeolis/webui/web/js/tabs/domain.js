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
        U.tbtn("upload", "Import", { title: "Import a *.xyz sample file", onclick: _importXyz }),
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

  /* Make the cards of a list draggable; onDrop(fromIdx, toIdx). */
  function _wireCardDrag(list, onDrop) {
    let fromIdx = null;
    list.addEventListener("dragstart", (ev) => {
      const card = ev.target.closest(".obj-card");
      if (!card) return;
      fromIdx = Number(card.dataset.idx);
      card.classList.add("dragging");
      ev.dataTransfer.effectAllowed = "move";
      ev.dataTransfer.setData("text/plain", "");   // Firefox needs data to drag
    });
    list.addEventListener("dragend", () => {
      fromIdx = null;
      for (const c of list.querySelectorAll(".obj-card")) {
        c.classList.remove("dragging", "drop-above", "drop-below");
      }
    });
    list.addEventListener("dragover", (ev) => {
      if (fromIdx === null) return;
      ev.preventDefault();
      const card = ev.target.closest(".obj-card");
      for (const c of list.querySelectorAll(".obj-card")) {
        c.classList.remove("drop-above", "drop-below");
      }
      if (!card) return;
      const rect = card.getBoundingClientRect();
      const below = ev.clientY > rect.top + rect.height / 2;
      card.classList.add(below ? "drop-below" : "drop-above");
    });
    list.addEventListener("drop", (ev) => {
      if (fromIdx === null) return;
      ev.preventDefault();
      const card = ev.target.closest(".obj-card");
      if (!card) return;
      const rect = card.getBoundingClientRect();
      const below = ev.clientY > rect.top + rect.height / 2;
      let toIdx = Number(card.dataset.idx) + (below ? 1 : 0);
      if (toIdx > fromIdx) toIdx -= 1;
      if (toIdx !== fromIdx) onDrop(fromIdx, toIdx);
      fromIdx = null;
    });
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

      const name = U.el("span", { class: "lp-name", title: "Double-click to rename" },
        entry.label || entry.path);
      name.addEventListener("dblclick", () => _renameSample(entry, name));

      const card = U.el("div", {
        class: `obj-card ${selectedIds.has(entry.id) ? "selected" : ""}`,
        draggable: "true", dataset: { idx, eid: entry.id },
      },
        U.el("span", { class: "drag-grip", title: "Drag to reorder (top = highest priority); selected cards move together" }, "⠿"),
        _eyeOrSpinner(layerId, layer),
        name,
        U.el("span", { class: "lp-mini" }, entry.source),
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
      return info && (!info.hidden || revealedTargets.has(name));
    });
    order.forEach((name, idx) => {
      const info = overview.targets[name];
      const layerId = `domain-${name}`;
      const layer = Layers.get(layerId) || { visible: false, id: layerId };

      const staleTitle = info.shape_ok
        ? "The grid changed after this file was interpolated — re-interpolate before using it"
        : "Shape does not match the current grid — re-interpolate";

      const eye = info.exists
        ? _eyeOrSpinner(layerId, layer, {
          disabled: info.stale,
          disabledTitle: staleTitle,
        })
        : U.el("span", { class: "eye off disabled", title: "File does not exist yet" }, "👁");

      // not required by the current config (e.g. veg under the grass
      // method, or an opt-in mask) → de-emphasise, keep it usable
      const card = U.el("div", {
        class: `obj-card ${info.needed ? "" : "not-needed"}`,
        draggable: "true", dataset: { idx },
      },
        U.el("span", { class: "drag-grip", title: "Drag to reorder" }, "⠿"),
        eye,
        U.el("span", { class: "lp-name" }, `${name} — ${info.file}`));

      if (!info.needed) {
        card.append(U.el("span", {
          class: "opt-badge",
          title: info.note || "not required by the current configuration",
        }, info.optional ? "optional" : "not needed"));
      }
      if (!info.exists) {
        card.append(U.el("span", { class: "lp-mini" }, "no file"));
      } else if (info.stale) {
        card.append(_warnIcon(staleTitle));
      }

      const actions = U.el("span", { class: "obj-actions" });
      actions.append(U.miniBtn("interp", "Interpolate…", () => _interpolateWizard(name, info)));
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

  /* "Add optional file…" — reveal an opt-in domain file (mask etc.) so it
   * can be interpolated; most domain files need no input and stay hidden. */
  function _addOptionalControl() {
    const hidden = Object.entries(overview.targets)
      .filter(([name, info]) => info.hidden && !revealedTargets.has(name));
    if (!hidden.length) return null;
    const sel = U.el("select", { style: "flex:1;min-width:0" },
      U.el("option", { value: "" }, "— add an optional domain file —"),
      ...hidden.map(([name, info]) =>
        U.el("option", { value: name }, `${name} (${info.file})`)));
    sel.addEventListener("change", () => {
      if (!sel.value) return;
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
          class: "obj-card", draggable: "true", dataset: { idx },
        },
          U.el("span", { class: "drag-grip", title: "Drag to reorder (top = highest priority)" }, "⠿"),
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
    const runBtn = U.el("button", { class: "primary" }, "Interpolate & save");
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
        U.toast(`Wrote ${out.file} (${U.fmtNum(out.min)} … ${U.fmtNum(out.max)})`, "ok");
        const cfg = await Api.get("/api/config");
        App.state.config = cfg.values;
        App.emit("config-changed", overview.targets[target].config_key);
        popup.close();
        await _refresh();
        _reloadLayer(`domain-${target}`);
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
