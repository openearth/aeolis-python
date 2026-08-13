/* Conditions tab: wind, water level and wave boundary conditions.
 *
 * Two ways to get data, per quantity:
 *  - Generate… (synthetic): segment tables (constant / linear /
 *    harmonic pieces with durations) with a live preview. When all
 *    segments have durations the written file covers just that time and
 *    AeoLiS repeats it cyclically - the preview shows that repetition.
 *  - Download… (from source): waterinfo.rws.nl or ERA5. Downloads land
 *    as separate RAW SERIES objects (not written to wind.txt directly):
 *    the user can inspect/clean them (drop NaN, fill, crop, resample,
 *    …) and then apply one to the AeoLiS input file. Downloads keep
 *    running in the background when the wizard is closed.
 */
"use strict";

const ConditionsTab = (() => {

  const KIND_TITLES = { wind: "Wind", tide: "Water levels", wave: "Waves" };
  const COLORS = TSPlot.COLORS;

  // The six physical variables the raw-data section works with. Each maps
  // to an AeoLiS input "kind" (file) and to the primary column of that
  // kind's series (used to seed pickers / windrose).
  const VARIABLES = [
    { id: "wind_speed",   label: "Wind speed",     kind: "wind", col: 0 },
    { id: "wind_dir",     label: "Wind direction", kind: "wind", col: 1 },
    { id: "water_level",  label: "Water level",    kind: "tide", col: 0 },
    { id: "wave_height",  label: "Wave height",    kind: "wave", col: 0 },
    { id: "wave_period",  label: "Wave period",    kind: "wave", col: 1 },
    { id: "wave_dir",     label: "Wave direction", kind: "wave", col: 1 },
  ];
  function _variableById(id) { return VARIABLES.find((v) => v.id === id) || VARIABLES[0]; }

  let overview = null;
  let rawEntries = [];
  let stationMarkers = [];
  const activeFetch = new Map();     // job id -> {kind, progressBar}
  const rawSeriesCache = new Map();  // entry id -> Promise(payload)
  const selectedRaw = new Set();     // multi-selected raw-series ids
  let lastRawClicked = null;         // shift+click range anchor

  function init() {
    Tabs.register("conditions", { enter: _refresh, leave: _clearStations });
    App.on("project", _refresh);
    // clicking anywhere outside a card / toolbar clears the raw selection
    U.deselectOnOutside(() => selectedRaw.size,
      () => { selectedRaw.clear(); if (overview) _build(); });
  }

  async function _refresh() {
    if (!App.state.project) return;
    try {
      overview = await Api.get("/api/conditions");
      rawEntries = (await Api.get("/api/conditions/raw")).entries || [];
    } catch (err) {
      U.toast(err.message, "error");
      return;
    }
    // prune selection of removed entries
    const knownRaw = new Set(rawEntries.map((e) => e.id));
    for (const id of [...selectedRaw]) if (!knownRaw.has(id)) selectedRaw.delete(id);
    _build();
    _registerGraphSources();
  }

  /* ================= graph-panel sources ================= */

  function _unitOfLabel(label) {
    return (String(label).match(/\[(.*)\]/) || [])[1] || "";
  }

  function _nameOfLabel(label) {
    return String(label).replace(/\s*\[.*\]$/, "");
  }

  /* Canonical physical variable for a column, so the graph picker can
   * group every source of the same quantity (wind.txt, raw waterinfo,
   * ERA5, …) together. */
  function _variableOf(kind, label) {
    const n = _nameOfLabel(label).toLowerCase();
    if (/dir|rtg/.test(n)) return `${KIND_TITLES[kind]} direction`;
    if (kind === "wind" || /speed|snelheid|u10|wind/.test(n)) return "Wind speed";
    if (/hs|hm0|height|hoogte|golfh/.test(n)) return "Wave height";
    if (/tp|tm|period|periode/.test(n)) return "Wave period";
    if (kind === "tide" || /level|wathte|water|tide|getij/.test(n)) return "Water level";
    return `${KIND_TITLES[kind]}: ${_nameOfLabel(label)}`;
  }

  /* Column indices of a magnitude+direction pair (for the windrose), or
   * null when the series has no direction column. */
  function _windroseCols(labels) {
    labels = labels || [];
    let mag = -1, dir = -1;
    labels.forEach((lab, i) => {
      const n = _nameOfLabel(lab).toLowerCase();
      if (/dir|rtg/.test(n)) { if (dir < 0) dir = i; }
      else if (/speed|snelheid|hs|hm0|height|hoogte|u10|golfh/.test(n)) { if (mag < 0) mag = i; }
    });
    return (mag >= 0 && dir >= 0)
      ? { mag, dir, magName: _nameOfLabel(labels[mag]), magUnit: _unitOfLabel(labels[mag]) }
      : null;
  }

  /* Every magnitude+direction pair available for a windrose, drawn from the
   * input files (wind/wave) and the raw series. Fetches fresh so the topbar
   * windrose works even before the Conditions tab has been opened. Each entry
   * is {label, magName, magUnit, load()} where load() yields the arrays. */
  async function windroseSources() {
    if (!App.state.project) return [];
    const ov = await Api.get("/api/conditions");
    const raws = (await Api.get("/api/conditions/raw")).entries || [];
    const out = [];
    for (const kind of ["wind", "wave"]) {
      const info = ov.kinds[kind];
      const cols = info && _windroseCols(info.labels);
      if (info && info.series && cols) {
        out.push({
          label: `${info.file || KIND_TITLES[kind]} (input file)`,
          magName: cols.magName, magUnit: cols.magUnit,
          t0_epoch: info.series.t0_epoch, t1_epoch: info.series.t1_epoch,
          load: async () => ({
            magnitude: info.series.columns[cols.mag],
            direction: info.series.columns[cols.dir],
            t_epoch: info.series.t_epoch,
          }),
        });
      }
    }
    for (const entry of raws) {
      const cols = _windroseCols(entry.labels);
      if (!cols) continue;
      out.push({
        label: `${entry.label} (raw)`,
        magName: cols.magName, magUnit: cols.magUnit,
        t0_epoch: entry.t0_epoch, t1_epoch: entry.t1_epoch,
        load: async () => {
          const [t, magArr] = await _loadRawColumn(entry, cols.mag);
          const [, dirArr] = await _loadRawColumn(entry, cols.dir);
          return { magnitude: magArr, direction: dirArr, t_epoch: t };
        },
      });
    }
    return out;
  }

  /* Repeat a series' sample points cyclically until *until* - exactly
   * what AeoLiS does with a boundary-condition file shorter than the
   * simulation (interp_circular / interp_circular_nearest). */
  function _tile(t, v, until) {
    const n = t.length;
    if (n < 2) return [t, v];
    const period = t[n - 1] - t[0];
    if (period <= 0 || t[n - 1] >= until - 1) return [t, v];
    const T = t.slice(), V = v.slice();
    for (let k = 1; t[0] + k * period < until && k < 4000; k += 1) {
      for (let i = 1; i < n; i += 1) {
        const tt = t[i] + k * period;
        if (tt > until) break;
        T.push(tt);
        V.push(v[i]);
      }
    }
    return [T, V];
  }

  function _simRange() {
    const t0 = overview.refdate_epoch + (overview.tstart || 0);
    const t1 = overview.refdate_epoch + (overview.tstop || 0);
    return [t0, t1];
  }

  function _registerGraphSources() {
    const [, simT1] = _simRange();

    for (const [kind, info] of Object.entries(overview.kinds)) {
      if (!info.series) {
        // input file removed/unreadable: drop its stale sources
        for (let i = 0; i < 3; i += 1) Graphs.unregisterSource(`cond-${kind}-${i}`);
        Playbar.removeSource(`cond-${kind}`);
        continue;
      }
      const s = info.series;
      const repeats = s.t1_epoch < simT1 - 1;
      info.labels.forEach((label, i) => {
        if (!s.columns[i]) return;
        const [t, v] = repeats
          ? _tile(s.t_epoch, s.columns[i], simT1)
          : [s.t_epoch, s.columns[i]];
        const isDir = /direction/i.test(label);
        Graphs.registerSource(`cond-${kind}-${i}`, {
          group: "Input",
          variable: _variableOf(kind, label),
          source: `${info.file || kind + ".txt"} (input)`,
          label: `${KIND_TITLES[kind]}: ${_nameOfLabel(label)}`,
          unit: _unitOfLabel(label),
          points: isDir,
          data: [t, v],
          // zoom-aware refetch for denser detail within the file's own range —
          // even when the series repeats (a long measured file that ends before
          // tstop is still marked "repeats", but its samples are far denser than
          // the full-span decimation, so windowing must stay enabled)
          loadWindow: (t0, t1) => _loadInputWindow(kind, i, t0, t1),
          range: [s.t0_epoch, Math.max(s.t1_epoch, repeats ? simT1 : s.t1_epoch)],
          repeatFrom: repeats ? s.t1_epoch : null,
        });
      });
      Playbar.setSource(`cond-${kind}`, s.t0_epoch, s.t1_epoch);
    }

    for (const entry of rawEntries) {
      (entry.labels || []).forEach((label, i) => {
        Graphs.registerSource(_rawSeriesKey(entry, i), {
          group: "Raw data",
          variable: _variableOf(entry.kind, label),
          source: `${entry.label} (raw)`,
          label: `${entry.label}: ${_nameOfLabel(label)}`,
          unit: _unitOfLabel(label),
          points: /direction/i.test(label),
          data: null,
          load: () => _loadRawColumn(entry, i),
          loadWindow: (t0, t1) => _loadRawColumn(entry, i, [t0, t1]),
          range: (entry.t0_epoch !== null && entry.t1_epoch !== null)
            ? [entry.t0_epoch, entry.t1_epoch] : null,
        });
      });
    }
  }

  function _rawSeriesKey(entry, col) { return `rawcond-${entry.id}-${col}`; }

  async function _loadRawColumn(entry, col, win) {
    // windowed fetch (zoom-aware): denser detail within [t0,t1], not cached
    if (win) {
      const res = await Api.get(
        `/api/conditions/raw_series?id=${entry.id}&tmin=${win[0]}&tmax=${win[1]}`);
      return [res.series.t_epoch, res.series.columns[col] || []];
    }
    if (!rawSeriesCache.has(entry.id)) {
      // don't cache failures: a transient error would stick forever
      const promise = Api.get(`/api/conditions/raw_series?id=${entry.id}`)
        .catch((err) => { rawSeriesCache.delete(entry.id); throw err; });
      rawSeriesCache.set(entry.id, promise);
    }
    const res = await rawSeriesCache.get(entry.id);
    return [res.series.t_epoch, res.series.columns[col] || []];
  }

  /* Higher-resolution slice of an input file within a zoom window. */
  async function _loadInputWindow(kind, col, t0, t1) {
    const res = await Api.get(`/api/conditions/series?kind=${kind}&tmin=${t0}&tmax=${t1}`);
    return [res.series.t_epoch, res.series.columns[col] || []];
  }

  /* ================= panel ================= */

  function _build() {
    const panel = document.getElementById("conditions-panel");
    U.keepScroll(panel);
    U.clear(panel);

    panel.append(U.el("div", { class: "muted", style: "font-size:12px" },
      `refdate ${overview.refdate} — simulation ${U.fmtDuration(overview.tstop - overview.tstart)}`));

    // ---- Raw data: download / generate / clean any variable ----
    const rawSec = U.section("Raw data", { count: rawEntries.length || null });
    rawSec.body.append(U.el("div", { class: "muted", style: "font-size:12px" },
      "Measured, reanalysis or synthetic series for any variable. Clean them here, "
      + "then build the input timeseries below."));
    rawSec.body.append(U.el("div", { class: "tbtn-row" },
      U.tbtn("download", "Download", {
        title: "Download measured/reanalysis data (pick a variable and source)",
        onclick: () => _sourceWizard(),
      }),
      U.tbtn("wand", "Generate", {
        title: "Generate a synthetic series (pick a variable)",
        onclick: () => _synthWizard(),
      }),
      (() => {
        const sel = rawEntries.filter((e) => selectedRaw.has(e.id));
        const btn = U.tbtn("trash", "Remove", {
          title: sel.length ? `Remove ${sel.length} selected series…`
            : "Select series first (click, Ctrl+click, Shift+click for a range)",
          onclick: () => _removeSelectedRaw(),
        });
        btn.disabled = !sel.length;
        return btn;
      })()));
    for (const fetching of activeFetch.values()) rawSec.body.append(fetching.row || fetching.progressBar.el);
    if (rawEntries.length) rawSec.body.append(_rawList(rawEntries));
    else rawSec.body.append(U.el("div", { class: "muted" }, "Nothing downloaded or generated yet."));
    panel.append(rawSec.wrap);

    // ---- Input timeseries: the wind/tide/wave files AeoLiS reads ----
    const inSec = U.section("Input timeseries", { count: null });
    inSec.body.append(U.el("div", { class: "muted", style: "font-size:12px" },
      "The wind / water-level / wave files AeoLiS reads — build (fill) each from the raw series above."));
    const fileList = U.el("div", { class: "obj-list" });
    for (const kind of ["wind", "tide", "wave"]) {
      fileList.append(_inputFileCard(kind));
    }
    inSec.body.append(fileList);
    panel.append(inSec.wrap);
  }

  function _inputFileCard(kind) {
    const info = overview.kinds[kind];
    // clear state: present & readable / present but broken / configured but
    // missing (e.g. a link that broke on duplicate) / never created
    let status, cls;
    if (info.exists && info.series) {
      status = `✔ ${info.file} (${info.series.n} rows)`; cls = "";
    } else if (info.exists) {
      status = `⚠ ${info.file} — ${info.error ? "unreadable" : "empty"}`; cls = "warn";
    } else if (info.file) {
      status = `⚠ missing: ${info.file}`; cls = "warn";
    } else {
      status = `${KIND_TITLES[kind]} file — not created yet`; cls = "muted";
    }

    // fixed button set (same order as the Domain cards); inapplicable
    // buttons are greyed out, never hidden
    const actions = U.el("span", { class: "obj-actions" });
    actions.append(U.miniBtn("interp", "Fill / create from raw series…", () => _fillWizard(kind)));
    actions.append(U.miniBtn("saveas", "Save as… (choose location; the config is repointed)",
      () => _saveFileAs(kind),
      { disabled: !info.series, disabledTitle: "no file to save yet" }));
    actions.append(U.miniBtn("open", "Load an existing file…", () => _loadFile(kind)));

    return U.el("div", { class: "obj-card" },
      U.el("span", { class: "lp-name", title: info.file || "" }, KIND_TITLES[kind]),
      U.el("span", { class: `lp-mini ${cls}`, title: info.error || info.file || "" }, status),
      actions);
  }

  function _rawList(entries) {
    const wrap = U.el("div", {});
    const list = U.el("div", { class: "obj-list" });
    const ids = entries.map((e) => e.id);
    for (const entry of entries) {
      const name = U.el("span", { class: "lp-name", title: "Double-click to rename" },
        entry.label);
      name.addEventListener("dblclick", () => _renameRaw(entry, name));

      const meta = `${KIND_TITLES[entry.kind] || entry.kind} · ${entry.rows || 0} rows` +
        (entry.nan ? ` · ${entry.nan} NaN` : "");

      const actions = U.el("span", { class: "obj-actions" });
      actions.append(U.miniBtn("copy", "Duplicate this series", () => _duplicateRaw(entry)));
      actions.append(U.miniBtn("modify", "Modify… (clean NaN, crop, arithmetic, …)",
        () => _rawModifyWizard(entry)));
      actions.append(U.miniBtn("search", "Clean up… (detect sentinel values & stuck-sensor stretches)",
        () => _rawCleanWizard(entry)));

      const card = U.el("div", {
        class: `obj-card ${selectedRaw.has(entry.id) ? "selected" : ""}`,
      },
        U.el("span", { class: "drag-grip", draggable: "true", title: "Drag to reorder" }, "⠿"),
        name,
        entry.nan
          ? U.el("span", { class: "warn-icon", title: `${entry.nan} NaN value(s) — clean before use` }, "⚠")
          : null,
        U.el("span", { class: "lp-mini" }, meta),
        actions);
      // click = select (Ctrl toggles, Shift selects a range)
      card.addEventListener("click", (ev) => {
        if (ev.target.closest("button, .eye, input, .drag-grip")) return;
        if (ev.shiftKey && lastRawClicked && ids.includes(lastRawClicked)) {
          const a = ids.indexOf(lastRawClicked), b = ids.indexOf(entry.id);
          for (let k = Math.min(a, b); k <= Math.max(a, b); k += 1) selectedRaw.add(ids[k]);
        } else if (ev.ctrlKey || ev.metaKey) {
          if (selectedRaw.has(entry.id)) selectedRaw.delete(entry.id);
          else selectedRaw.add(entry.id);
        } else if (selectedRaw.size === 1 && selectedRaw.has(entry.id)) {
          selectedRaw.clear();
        } else {
          selectedRaw.clear();
          selectedRaw.add(entry.id);
        }
        lastRawClicked = entry.id;
        _build();
      });
      list.append(card);
    }
    // reorder is cosmetic (raw order), but keep it consistent with Domain
    U.wireSortable(list, (from, to) => {
      const arr = rawEntries.slice();
      const [moved] = arr.splice(from, 1);
      arr.splice(to, 0, moved);
      rawEntries = arr;
      _build();
    });
    wrap.append(list);
    if (entries.length > 1) {
      wrap.append(U.el("div", { class: "muted", style: "font-size:11px;margin-top:4px" },
        "Click to select, Ctrl+click to add, Shift+click for a range — selected series are removed together."));
    }
    return wrap;
  }

  async function _duplicateRaw(entry) {
    try {
      await Api.post("/api/conditions/raw_duplicate", { id: entry.id });
      U.toast("Series duplicated", "ok");
      _refresh();
    } catch (err) {
      U.toast(err.message, "error");
    }
  }

  function _removeSelectedRaw() {
    const entries = rawEntries.filter((e) => selectedRaw.has(e.id));
    if (!entries.length) return;
    const popup = Popup.open({ title: `Remove ${entries.length} raw series`, width: 440 });
    const okBtn = U.el("button", { class: "danger" }, `Remove ${entries.length}`);
    okBtn.addEventListener("click", async () => {
      try {
        okBtn.disabled = true;
        for (const entry of entries) {
          await Api.post("/api/conditions/raw_delete", { id: entry.id });
          (entry.labels || []).forEach((_, i) => Graphs.unregisterSource(_rawSeriesKey(entry, i)));
          rawSeriesCache.delete(entry.id);
          selectedRaw.delete(entry.id);
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
      U.el("div", { style: "font-size:13px" },
        "Remove these raw series from the project? Files already applied to the input file are not affected."),
      U.el("div", { class: "obj-list", style: "margin:8px 0;max-height:160px;overflow-y:auto" },
        ...entries.map((entry) => U.el("div", { class: "obj-card" },
          U.el("span", { class: "lp-name" }, entry.label),
          U.el("span", { class: "lp-mini" }, entry.source || "")))),
      U.el("div", { class: "btn-row", style: "justify-content:flex-end" }, cancelBtn, okBtn));
  }

  function _renameRaw(entry, nameNode) {
    const input = U.el("input", { type: "text", value: entry.label || "", style: "flex:1;font-size:12px" });
    nameNode.replaceWith(input);
    input.focus(); input.select();
    const commit = async () => {
      const name = input.value.trim();
      if (name && name !== entry.label) {
        await Api.post("/api/conditions/raw_rename", { id: entry.id, name });
      }
      _refresh();
    };
    input.addEventListener("blur", commit);
    input.addEventListener("keydown", (ev) => {
      if (ev.key === "Enter") input.blur();
      if (ev.key === "Escape") { input.value = entry.label; input.blur(); }
    });
  }

  async function _applyRaw(entry) {
    try {
      const res = await Api.post("/api/conditions/raw_apply", { id: entry.id });
      U.toast(`Wrote ${res.file} (${res.rows} rows)`, "ok");
      await _reloadConfig(entry.kind);
      _refresh();
    } catch (err) {
      U.toast(err.message, "error");
    }
  }

  function _deleteRaw(entry) {
    const popup = Popup.open({ title: `Delete ${entry.label}`, width: 400 });
    const okBtn = U.el("button", { class: "danger" }, "Delete");
    okBtn.addEventListener("click", async () => {
      try {
        await Api.post("/api/conditions/raw_delete", { id: entry.id });
        (entry.labels || []).forEach((_, i) =>
          Graphs.unregisterSource(_rawSeriesKey(entry, i)));
        rawSeriesCache.delete(entry.id);
        popup.close();
        _refresh();
      } catch (err) {
        U.toast(err.message, "error");
      }
    });
    const cancelBtn = U.el("button", { class: "ghost" }, "Cancel");
    cancelBtn.addEventListener("click", popup.close);
    popup.body.append(
      U.el("div", { style: "font-size:13px" },
        "Delete this raw series from the project? Files already applied to ",
        "wind.txt/tide.txt/waves.txt are not affected."),
      U.el("div", { class: "btn-row", style: "justify-content:flex-end" }, cancelBtn, okBtn));
  }

  /* ================= raw modify wizard ================= */

  const RAW_OPS = [
    ["drop_nan", "remove rows with NaN"],
    ["fill_nan", "fill NaN with value"],
    ["fill_from", "fill NaN from another series"],
    ["crop", "crop to date range"],
    ["add", "add value"],
    ["subtract", "subtract value"],
    ["multiply", "multiply by value"],
    ["clip_max", "cap at max (values above → value)"],
    ["clip_min", "cap at min (values below → value)"],
    ["set", "set to value"],
  ];
  const VALUE_OPS = new Set(["fill_nan", "add", "subtract", "multiply",
    "clip_max", "clip_min", "set"]);
  const COLUMN_OPS = new Set(["fill_nan", "add", "subtract", "multiply",
    "clip_max", "clip_min", "set"]);

  function _rawModifyWizard(entry) {
    const popup = Popup.open({ title: `Modify — ${entry.label}`, width: 520 });

    const iso = (epoch) => Number.isFinite(epoch)
      ? new Date(epoch * 1000).toISOString().slice(0, 10) : "";

    const opSel = U.el("select", {},
      ...RAW_OPS.map(([v, l]) => U.el("option", { value: v }, l)));
    const value = U.el("input", { type: "text", placeholder: "value", style: "width:90px" });
    const colSel = U.el("select", {},
      U.el("option", { value: "all" }, "all columns"),
      ...(entry.labels || []).map((label, i) => U.el("option", { value: i }, label)));
    const otherSel = U.el("select", {});
    for (const other of rawEntries.filter((e) => e.kind === entry.kind && e.id !== entry.id)) {
      otherSel.append(U.el("option", { value: other.id }, other.label));
    }
    const date0 = U.el("input", { type: "text", value: iso(entry.t0_epoch), title: "YYYY-MM-DD" });
    const date1 = U.el("input", { type: "text", value: iso(entry.t1_epoch), title: "YYYY-MM-DD" });
    const rsVal = U.el("input", { type: "text", value: "1", style: "width:60px" });
    const rsUnit = U.el("select", {},
      U.el("option", { value: 60 }, "minutes"),
      U.el("option", { value: 3600, selected: "" }, "hours"),
      U.el("option", { value: 86400 }, "days"));

    const rows = {
      value: U.el("div", { class: "form-row" }, U.el("label", {}, "value"), value),
      column: U.el("div", { class: "form-row" }, U.el("label", {}, "apply to"), colSel),
      other: U.el("div", { class: "form-row" }, U.el("label", {}, "fill from"), otherSel),
      crop: U.el("div", { class: "form-row" }, U.el("label", {}, "from / to"), date0, date1),
      resample: U.el("div", { class: "form-row" }, U.el("label", {}, "interval"), rsVal, rsUnit),
    };
    const syncRows = () => {
      const op = opSel.value;
      rows.value.style.display = VALUE_OPS.has(op) ? "" : "none";
      rows.column.style.display = COLUMN_OPS.has(op) ? "" : "none";
      rows.other.style.display = op === "fill_from" ? "" : "none";
      rows.crop.style.display = op === "crop" ? "" : "none";
      rows.resample.style.display = op === "resample" ? "" : "none";
    };
    opSel.addEventListener("change", syncRows);

    const saveNew = U.el("input", { type: "checkbox", id: "rm-new" });
    const newName = U.el("input", { type: "text", placeholder: "new series name" });
    newName.addEventListener("input", () => { if (newName.value) saveNew.checked = true; });

    const applyBtn = U.el("button", { class: "primary" }, "Apply");
    applyBtn.addEventListener("click", async () => {
      const body = { id: entry.id, op: opSel.value };
      if (VALUE_OPS.has(opSel.value)) body.value = Number(value.value);
      if (COLUMN_OPS.has(opSel.value) && colSel.value !== "all") body.column = Number(colSel.value);
      if (opSel.value === "fill_from") body.other = otherSel.value;
      if (opSel.value === "crop") { body.date0 = date0.value.trim(); body.date1 = date1.value.trim(); }
      if (opSel.value === "resample") {
        body.value = Number(rsVal.value) * Number(rsUnit.value);
      }
      if (saveNew.checked) {
        if (!newName.value.trim()) { U.toast("Enter a name for the new series", "error"); return; }
        body.save_as = newName.value.trim();
      }
      try {
        applyBtn.disabled = true;
        const res = await Api.post("/api/conditions/raw_modify", body);
        rawSeriesCache.delete(entry.id);
        rawSeriesCache.delete(res.entry.id);
        U.toast(`Saved ${res.entry.label} (${res.entry.rows} rows` +
          (res.entry.nan ? `, ${res.entry.nan} NaN left)` : ")"), "ok");
        popup.close();
        await _refresh();
        _reregisterRaw(res.entry);
      } catch (err) {
        U.toast(err.message, "error");
        applyBtn.disabled = false;
      }
    });

    popup.body.append(
      U.el("div", { class: "muted", style: "font-size:12px" },
        `${entry.rows || 0} rows, ${entry.nan || 0} NaN — ` +
        `${U.fmtDate(entry.t0_epoch)} — ${U.fmtDate(entry.t1_epoch)}`),
      U.el("div", { class: "form-row", style: "margin-top:8px" },
        U.el("label", {}, "operation"), opSel),
      ...Object.values(rows),
      U.el("div", { class: "choice-row", style: "margin-top:8px" },
        saveNew, U.el("label", { for: "rm-new" }, "Save as a new series"),
        U.el("span", { class: "grow" }), newName),
      U.el("div", { class: "btn-row", style: "justify-content:flex-end" }, applyBtn),
    );
    syncRows();
  }

  /* Clean-up wizard: sensors report flaws as sentinel values (999 m
   * waves) or long constant stretches (0 m/s wind, 72 deg direction for
   * days). Detect them, preview the counts, then replace with NaN,
   * interpolate over, or remove the samples. */
  function _rawCleanWizard(entry) {
    const popup = Popup.open({ title: `Clean up — ${entry.label}`, width: 540 });

    const colSel = U.el("select", {},
      U.el("option", { value: "all" }, "all columns"),
      ...(entry.labels || []).map((label, i) => U.el("option", { value: i }, label)));

    const sentOn = U.el("input", { type: "checkbox", id: "cl-sent", checked: "" });
    const sentVals = U.el("input", { type: "text", value: "999, -999, 9999", style: "flex:1",
      title: "Comma-separated exact values treated as 'no data'" });

    const runOn = U.el("input", { type: "checkbox", id: "cl-run", checked: "" });
    const runMin = U.el("input", { type: "number", value: "24", style: "width:64px",
      title: "Flag stretches where the value does not change for at least this many time steps" });
    const runVal = U.el("input", { type: "text", placeholder: "any value", style: "width:90px",
      title: "Only flag constant stretches of this value (empty = any repeated value)" });

    const actionSel = U.el("select", {},
      U.el("option", { value: "nan" }, "replace with NaN"),
      U.el("option", { value: "interp" }, "replace + interpolate over"),
      U.el("option", { value: "remove" }, "remove the samples"));

    const saveNew = U.el("input", { type: "checkbox", id: "cl-new" });
    const newName = U.el("input", { type: "text", placeholder: "new series name" });
    newName.addEventListener("input", () => { if (newName.value) saveNew.checked = true; });

    const result = U.el("div", { class: "muted", style: "font-size:12px;margin-top:6px" });

    const buildBody = (preview) => {
      const rules = {};
      if (sentOn.checked) {
        rules.sentinels = sentVals.value.split(",").map((s) => Number(s.trim()))
          .filter((v) => Number.isFinite(v));
      }
      if (runOn.checked) {
        rules.run_min = Number(runMin.value) || 24;
        if (runVal.value.trim() !== "") rules.run_value = Number(runVal.value);
      }
      const body = { id: entry.id, rules, action: actionSel.value, preview };
      if (colSel.value !== "all") body.column = Number(colSel.value);
      if (!preview && saveNew.checked) body.save_as = newName.value.trim();
      return body;
    };

    const previewBtn = U.el("button", { class: "ghost" }, "Preview");
    previewBtn.addEventListener("click", async () => {
      try {
        previewBtn.disabled = true;
        const res = await Api.post("/api/conditions/raw_clean", buildBody(true));
        U.clear(result);
        result.append(U.el("div", {},
          `${res.total} of ${res.rows * res.report.length} samples flagged:`));
        for (const r of res.report) {
          const parts = Object.entries(r.rules).map(([k, n]) => `${k}: ${n}`).join(" · ");
          result.append(U.el("div", {}, `— ${r.column}: ${r.flagged} (${parts || "none"})`));
        }
      } catch (err) { U.toast(err.message, "error"); }
      previewBtn.disabled = false;
    });

    const applyBtn = U.el("button", { class: "primary" }, "Clean");
    applyBtn.addEventListener("click", async () => {
      if (saveNew.checked && !newName.value.trim()) {
        U.toast("Enter a name for the new series", "error"); return;
      }
      try {
        applyBtn.disabled = true;
        const res = await Api.post("/api/conditions/raw_clean", buildBody(false));
        rawSeriesCache.delete(entry.id);
        rawSeriesCache.delete(res.entry.id);
        U.toast(`Cleaned ${res.total} samples → ${res.entry.label} ` +
          `(${res.entry.rows} rows${res.entry.nan ? `, ${res.entry.nan} NaN` : ""})`, "ok");
        popup.close();
        await _refresh();
        _reregisterRaw(res.entry);
      } catch (err) {
        U.toast(err.message, "error");
        applyBtn.disabled = false;
      }
    });

    popup.body.append(
      U.el("div", { class: "muted", style: "font-size:12px" },
        `${entry.rows || 0} rows, ${entry.nan || 0} NaN — ` +
        `${U.fmtDate(entry.t0_epoch)} — ${U.fmtDate(entry.t1_epoch)}`),
      U.el("div", { class: "form-row", style: "margin-top:8px" },
        U.el("label", {}, "check"), colSel),
      U.el("span", { class: "fg-label", style: "margin-top:8px" }, "Detect"),
      U.el("div", { class: "choice-row" }, sentOn,
        U.el("label", { for: "cl-sent" }, "sentinel values"), sentVals),
      U.el("div", { class: "choice-row" }, runOn,
        U.el("label", { for: "cl-run" }, "constant for ≥"), runMin,
        U.el("span", {}, "steps, value"), runVal),
      U.el("span", { class: "fg-label", style: "margin-top:8px" }, "Then"),
      U.el("div", { class: "form-row" }, U.el("label", {}, "action"), actionSel),
      U.el("div", { class: "muted", style: "font-size:11.5px" },
        "Tip: interpolating across a direction column can cut the 360° wrap — "
        + "prefer NaN + the Fill wizard for wind direction."),
      U.el("div", { class: "choice-row", style: "margin-top:8px" },
        saveNew, U.el("label", { for: "cl-new" }, "Save as a new series"),
        U.el("span", { class: "grow" }), newName),
      result,
      U.el("div", { class: "btn-row", style: "justify-content:flex-end" }, previewBtn, applyBtn),
    );
  }

  /* re-register a modified raw entry's series with fresh data */
  function _reregisterRaw(entry) {
    (entry.labels || []).forEach((label, i) => {
      Graphs.registerSource(_rawSeriesKey(entry, i), {
        group: "Raw data",
        variable: _variableOf(entry.kind, label),
        source: `${entry.label} (raw)`,
        label: `${entry.label}: ${_nameOfLabel(label)}`,
        unit: _unitOfLabel(label),
        points: /direction/i.test(label),
        data: null,
        load: () => _loadRawColumn(entry, i),
        range: (entry.t0_epoch !== null && entry.t1_epoch !== null)
          ? [entry.t0_epoch, entry.t1_epoch] : null,
      });
    });
  }

  /* ================= synthetic wizard ================= */

  const SEG_TYPES = ["constant", "linear", "harmonic"];
  const SEG_PARAMS = {
    constant: [["value", "value"]],
    linear: [["start", "start (empty = continue)"], ["end", "end"]],
    harmonic: [["mean", "mean"], ["amplitude", "ampl."], ["period", "period [h]"], ["phase", "phase [°]"]],
  };
  const HOUR_KEYS = new Set(["period", "duration"]);

  /* One editable segment table for a quantity (e.g. wind speed).
   * Returns {el, spec} - spec() gives {type:"segments", segments:[...]}. */
  function _segmentTable(label, defaults, onInput) {
    const rows = [];   // each: {typeSel, inputs: {key: input}, durInput, el}
    const body = U.el("div", { class: "segtable-body" });

    const addRow = (preset = {}) => {
      const typeSel = U.el("select", { class: "seg-type" },
        ...SEG_TYPES.map((t) => U.el("option", { value: t, selected: t === (preset.type || "constant") ? "" : null }, t)));
      const paramBox = U.el("span", { class: "seg-params" });
      const inputs = {};
      const renderParams = () => {
        U.clear(paramBox);
        for (const key of Object.keys(inputs)) delete inputs[key];
        for (const [key, ph] of SEG_PARAMS[typeSel.value]) {
          const input = U.el("input", {
            type: "text", class: "seg-in", placeholder: ph,
            value: preset[key] ?? "",
            title: ph,
          });
          input.addEventListener("input", onInput);
          inputs[key] = input;
          paramBox.append(input);
        }
        onInput();
      };
      typeSel.addEventListener("change", renderParams);

      const durInput = U.el("input", {
        type: "text", class: "seg-in seg-dur", placeholder: "rest",
        value: preset.duration ?? "",
        title: "Segment duration in hours (empty = to the end of the simulation)",
      });
      durInput.addEventListener("input", onInput);

      const delBtn = U.el("button", { class: "mini-btn danger-hover", title: "Remove segment" }, "✕");
      const rowEl = U.el("div", { class: "segtable-row" },
        typeSel, paramBox, durInput, delBtn);
      const row = { typeSel, inputs, durInput, el: rowEl };
      delBtn.addEventListener("click", () => {
        if (rows.length <= 1) return;
        rows.splice(rows.indexOf(row), 1);
        rowEl.remove();
        onInput();
      });
      rows.push(row);
      body.append(rowEl);
      renderParams();
    };

    addRow(defaults);

    const addBtn = U.el("button", { class: "ghost seg-add" }, "＋ add segment");
    addBtn.addEventListener("click", () => addRow({}));

    const el = U.el("div", { class: "form-group segtable" },
      U.el("span", { class: "fg-label" }, label),
      U.el("div", { class: "segtable-head" },
        U.el("span", {}, "type"), U.el("span", {}, "parameters"),
        U.el("span", { title: "duration in hours (empty = to the end)" }, "dur. [h]"),
        U.el("span", {})),
      body, addBtn);

    const spec = () => {
      const segments = [];
      for (const row of rows) {
        const seg = { type: row.typeSel.value };
        for (const [key] of SEG_PARAMS[row.typeSel.value]) {
          const raw = row.inputs[key].value.trim();
          if (raw === "") continue;
          let v = Number(raw);
          if (!Number.isFinite(v)) continue;
          if (HOUR_KEYS.has(key)) v *= 3600;
          seg[key] = v;
        }
        const dur = Number(row.durInput.value);
        if (Number.isFinite(dur) && dur > 0) seg.duration = dur * 3600;
        segments.push(seg);
      }
      return { type: "segments", segments };
    };

    return { el, spec };
  }

  // Segment-table preset + label for each variable (one column per variable).
  const SYNTH_FORMS = {
    wind_speed:  ["Wind speed [m/s]",     { value: 10 }],
    wind_dir:    ["Wind direction [deg]", { value: 270 }],
    water_level: ["Water level [m]",      { type: "harmonic", mean: 0, amplitude: 1, period: 12.42 }],
    wave_height: ["Wave height Hs [m]",   { value: 1 }],
    wave_period: ["Wave period Tp [s]",   { value: 6 }],
    wave_dir:    ["Wave direction [deg]", { value: 300 }],
  };

  function _synthWizard(initialVar) {
    const popup = Popup.open({ title: "Generate synthetic series", width: 700 });
    let variable = _variableById(initialVar || "wind_speed").id;
    let form = null;
    const previewPlots = [];
    const previewEl = U.el("div", { class: "wizard-preview" });
    const aliasNote = U.el("div", { class: "muted", style: "font-size:11.5px;color:#b45309" });
    const formHost = U.el("div", {});
    // guards out-of-order async previews (rapid variable switching used to
    // let a stale preview resolve into a rebuilt/closed wizard and crash)
    let previewToken = 0;

    const varSel = U.el("select", {},
      ...VARIABLES.map((v) => U.el("option", {
        value: v.id, selected: v.id === variable ? "" : null,
      }, v.label)));

    const renderPreview = (res, token) => {
      if (token !== previewToken || !previewEl.isConnected) return;
      for (const p of previewPlots.splice(0)) p.destroy();
      U.clear(previewEl);
      const s = res.series;
      const repeatFrom = res.repeat_from_epoch;
      res.labels.forEach((label, i) => {
        const unit = _unitOfLabel(label);
        previewPlots.push(TSPlot.create(previewEl, {
          title: label,
          data: [s.t_epoch, s.columns[i]],
          series: [{ label: _nameOfLabel(label), unit,
            color: COLORS[i % COLORS.length],
            points: /direction/i.test(label) }],
          yRange: /direction/i.test(label) ? [0, 360] : null,
          width: Math.min(620, window.innerWidth * 0.55),
          height: 170,
          drawExtra: Number.isFinite(repeatFrom)
            ? (u) => _drawRepeatShade(u, repeatFrom) : null,
        }));
      });
      if (Number.isFinite(repeatFrom)) {
        previewEl.append(U.el("div", { class: "muted", style: "font-size:11.5px" },
          "The shaded part is not written to the file — AeoLiS repeats the ",
          "series automatically when it is shorter than the simulation."));
      }
    };

    const collectBody = () => ({
      variable,
      dt: Number(dtInput.value) * 3600 || 3600,
      segments: form.spec(),
    });

    const checkAliasing = (body) => {
      const dt = body.dt;
      let worst = null;
      for (const seg of (body.segments.segments || [])) {
        if (seg.type === "harmonic" && seg.period && seg.period < 4 * dt) worst = seg.period;
      }
      aliasNote.textContent = worst !== null
        ? `⚠ the output step (${U.fmtNum(dt / 3600, 3)} h) is coarse for a ` +
          `${U.fmtNum(worst / 3600, 3)} h period — the sine will look jagged; use a smaller output step`
        : "";
    };

    const updatePreview = U.debounce(async () => {
      const body = collectBody();
      checkAliasing(body);
      const token = ++previewToken;
      try {
        renderPreview(await Api.post("/api/conditions/preview", body), token);
      } catch (err) {
        console.warn("preview failed", err.message);
      }
    }, 350);

    const buildForm = () => {
      U.clear(formHost);
      const [label, preset] = SYNTH_FORMS[variable] || SYNTH_FORMS.wind_speed;
      form = _segmentTable(label, preset, updatePreview);
      formHost.append(form.el);
      updatePreview();
    };
    varSel.addEventListener("change", () => {
      variable = varSel.value;
      buildForm();
    });

    const dtInput = U.el("input", { type: "text", value: "1", style: "width:70px" });
    dtInput.addEventListener("input", updatePreview);

    const saveBtn = U.el("button", { class: "primary" }, "Generate raw series");
    saveBtn.addEventListener("click", async () => {
      try {
        saveBtn.disabled = true;
        const res = await Api.post("/api/conditions/synthetic_raw", collectBody());
        U.toast(`Generated ${res.entry.label} (${res.entry.rows} rows)`, "ok");
        popup.close();
        _refresh();
      } catch (err) {
        U.toast(err.message, "error");
        saveBtn.disabled = false;
      }
    });

    popup.body.append(
      U.el("div", { class: "form-row" }, U.el("label", {}, "variable"), varSel),
      formHost,
      U.el("div", { class: "form-row" }, U.el("label", {}, "output step [h]"), dtInput),
      aliasNote,
      U.el("span", { class: "fg-label" }, "Preview (incl. repetition over the simulation)"),
      previewEl,
      U.el("div", { class: "btn-row" }, saveBtn),
    );
    buildForm();
  }

  /* ================= fill / save-as / load input files ================= */

  /* One output column's source picker: a priority-ordered, drag-sortable,
   * checkbox list of raw-series columns (top = highest priority; lower ones
   * only fill samples still NaN), plus how to fill whatever gaps remain. */
  function _fillColumnBlock(colLabel, opts, defaultKey, note = null) {
    const keyOf = (o) => `${o.id}:${o.column}`;
    const byKey = new Map(opts.map((o) => [keyOf(o), o]));
    let order = opts.map(keyOf);
    if (defaultKey) order = [defaultKey, ...order.filter((k) => k !== defaultKey)];
    const checked = new Set(defaultKey ? [defaultKey] : []);

    const listEl = U.el("div", { class: "obj-list" });
    const render = () => {
      U.clear(listEl);
      order.forEach((k) => {
        const o = byKey.get(k);
        if (!o) return;
        const cb = U.el("input", { type: "checkbox" });
        cb.checked = checked.has(k);
        cb.addEventListener("change", () => { if (cb.checked) checked.add(k); else checked.delete(k); });
        listEl.append(U.el("div", { class: "obj-card" },
          U.el("span", { class: "drag-grip", draggable: "true", title: "Drag to reorder (top = highest priority)" }, "⠿"),
          cb,
          U.el("span", { class: "lp-name" }, o.text)));
      });
    };
    render();
    U.wireSortable(listEl, (from, to) => {
      const [m] = order.splice(from, 1);
      order.splice(to, 0, m);
      render();
    });

    const methodSel = U.el("select", {},
      U.el("option", { value: "linear" }, "linear interpolation"),
      U.el("option", { value: "nearest" }, "nearest value"),
      U.el("option", { value: "value" }, "fill value"),
      U.el("option", { value: "series" }, "another series"));
    const valInput = U.el("input", { type: "text", value: "0", style: "width:80px" });
    const otherSel = U.el("select", { style: "flex:1;min-width:0" },
      ...opts.map((o) => U.el("option", { value: keyOf(o) }, o.text)));
    const valRow = U.el("div", { class: "form-row" }, U.el("label", {}, "value"), valInput);
    const otherRow = U.el("div", { class: "form-row" }, U.el("label", {}, "from series"), otherSel);
    const syncMethod = () => {
      valRow.style.display = methodSel.value === "value" ? "" : "none";
      otherRow.style.display = methodSel.value === "series" ? "" : "none";
    };
    methodSel.addEventListener("change", syncMethod);
    syncMethod();

    const el = U.el("div", { class: "form-group" },
      U.el("span", { class: "fg-label" }, `${colLabel} — sources (top = highest priority)`),
      note ? U.el("div", { class: "muted", style: "font-size:11px" }, `⚠ ${note}`) : "",
      listEl,
      U.el("div", { class: "form-row" }, U.el("label", {}, "fill remaining gaps"), methodSel),
      valRow, otherRow);

    const collect = () => {
      const sources = order.filter((k) => checked.has(k))
        .map((k) => { const o = byKey.get(k); return { id: o.id, column: o.column }; });
      const fill = { method: methodSel.value };
      if (methodSel.value === "value") fill.value = Number(valInput.value);
      if (methodSel.value === "series") {
        const o = byKey.get(otherSel.value);
        if (o) fill.other = { id: o.id, column: o.column };
      }
      return { sources, fill };
    };
    return { el, collect };
  }

  function _fillWizard(kind) {
    const info = overview.kinds[kind];
    const labels = info.labels || [];
    const popup = Popup.open({ title: `Fill ${info.file || KIND_TITLES[kind]} from raw series`, width: 600 });
    if (!rawEntries.length) {
      popup.body.append(U.el("div", { class: "muted" },
        "No raw series yet — download or generate some in the Raw data section first."));
      return;
    }
    // flat list of every (raw series, column) as a pickable option, tagged
    // with its canonical physical variable ("Wave height", "Wind speed", …)
    const opts = [];
    for (const e of rawEntries) {
      (e.labels || []).forEach((lab, ci) => {
        opts.push({ id: e.id, column: ci, entry: e,
          variable: _variableOf(e.kind, lab),
          text: `${e.label}: ${_nameOfLabel(lab)}` });
      });
    }
    // one prioritized-source block per output column, offering ONLY raw
    // columns of the SAME physical variable (an Hs column must not list Tp
    // or wind series); default-select the best match (same kind + column)
    const blocks = labels.map((lab, i) => {
      const want = _variableOf(kind, lab);
      let colOpts = opts.filter((o) => o.variable === want);
      let note = null;
      if (!colOpts.length) {   // nothing of this quantity → show all, flagged
        colOpts = opts;
        note = `no raw series with ${want.toLowerCase()} found — showing everything`;
      }
      let match = colOpts.find((o) => o.entry.kind === kind && o.column === i);
      if (!match) match = colOpts[0];
      const defaultKey = match ? `${match.id}:${match.column}` : null;
      return { label: lab, block: _fillColumnBlock(lab, colOpts, defaultKey, note) };
    });

    const rsCb = U.el("input", { type: "checkbox", id: "fill-rs" });
    const rsVal = U.el("input", { type: "text", value: "1", style: "width:60px", disabled: "" });
    const rsUnit = U.el("select", { disabled: "" },
      U.el("option", { value: 60 }, "minutes"),
      U.el("option", { value: 3600, selected: "" }, "hours"),
      U.el("option", { value: 86400 }, "days"));
    rsCb.addEventListener("change", () => {
      rsVal.disabled = !rsCb.checked;
      rsUnit.disabled = !rsCb.checked;
    });
    const fname = U.el("input", { type: "text", value: info.file || `${kind}.txt` });

    const fillBtn = U.el("button", { class: "primary" }, "Fill & save");
    fillBtn.addEventListener("click", async () => {
      const columns = blocks.map((b) => b.block.collect());
      if (columns.some((c) => !c.sources.length)) {
        U.toast("Select at least one source for each column", "error");
        return;
      }
      const body = { kind, columns };
      if (rsCb.checked) {
        const w = Number(rsVal.value) * Number(rsUnit.value);
        if (Number.isFinite(w) && w > 0) body.resample = w;
      }
      if (fname.value.trim()) body.filename = fname.value.trim();
      try {
        fillBtn.disabled = true;
        const res = await Api.post("/api/conditions/fill", body);
        U.toast(`Wrote ${res.file} (${res.rows} rows)`, "ok");
        await _reloadConfig(kind);
        popup.close();
        _refresh();
      } catch (err) {
        U.toast(err.message, "error");
        fillBtn.disabled = false;
      }
    });

    popup.body.append(
      U.el("div", { class: "muted", style: "font-size:12px" },
        `Build each column of ${KIND_TITLES[kind].toLowerCase()} from one or more raw series. `
        + "Higher in the list = higher priority; where the top series has gaps "
        + "(missing samples or NaN), the next series takes over, and the chosen "
        + "method fills whatever remains."),
      ...blocks.map((b) => b.block.el),
      U.el("div", { class: "form-row" },
        U.el("label", { for: "fill-rs" }, "resample to interval means"), rsCb, rsVal, rsUnit),
      U.el("div", { class: "form-row" }, U.el("label", {}, "save as"), fname),
      U.el("div", { class: "btn-row", style: "justify-content:flex-end" }, fillBtn),
    );
  }

  /* The folder the kind's current file lives in (pickers start there). */
  function _condFileDir(kind) {
    const root = App.state.project ? App.state.project.root : "";
    const file = ((overview.kinds || {})[kind] || {}).file || "";
    if (/^[A-Za-z]:[\\/]/.test(file) || file.startsWith("/")) {
      return file.replace(/[\\/][^\\/]*$/, "") || root;
    }
    const dir = file.includes("/") || file.includes("\\")
      ? file.replace(/[\\/][^\\/]*$/, "") : "";
    return dir ? `${root}\\${dir}` : root;
  }

  async function _saveFileAs(kind) {
    const info = overview.kinds[kind];
    const path = await Api.pickFile({
      title: `Save ${KIND_TITLES[kind]} file as…`,
      save: true,
      patterns: [["Text files", "*.txt"], ["All files", "*.*"]],
      initial: _condFileDir(kind),
      filename: info.file || `${kind}.txt`,
    }).catch((err) => { U.toast(err.message, "error"); return null; });
    if (!path) return;
    try {
      const res = await Api.post("/api/conditions/save_file_as", { kind, path });
      U.toast(`Saved ${res.file}`, "ok");
      await _reloadConfig(kind);
      _refresh();
    } catch (err) { U.toast(err.message, "error"); }
  }

  async function _loadFile(kind) {
    const path = await Api.pickFile({
      title: `Load a ${KIND_TITLES[kind].toLowerCase()} file`,
      patterns: [["Text files", "*.txt"], ["All files", "*.*"]],
      initial: _condFileDir(kind),
    }).catch((err) => { U.toast(err.message, "error"); return null; });
    if (!path) return;
    try {
      const res = await Api.post("/api/conditions/load_file", { kind, path });
      U.toast(`Loaded ${res.file} (${res.rows} rows)`, "ok");
      await _reloadConfig(kind);
      _refresh();
    } catch (err) { U.toast(err.message, "error"); }
  }

  function _drawRepeatShade(u, repeatFrom) {
    const { min, max } = u.scales.x;
    if (repeatFrom >= max) return;
    const ctx = u.ctx;
    const x0 = u.valToPos(Math.max(repeatFrom, min), "x", true);
    const x1 = u.bbox.left + u.bbox.width;
    ctx.save();
    ctx.fillStyle = "rgba(217, 169, 78, .12)";
    ctx.fillRect(x0, u.bbox.top, x1 - x0, u.bbox.height);
    if (repeatFrom >= min) {
      ctx.strokeStyle = "rgba(217, 169, 78, .8)";
      ctx.setLineDash([4, 4]);
      ctx.beginPath();
      ctx.moveTo(x0, u.bbox.top);
      ctx.lineTo(x0, u.bbox.top + u.bbox.height);
      ctx.stroke();
    }
    ctx.restore();
  }

  /* ================= from-source wizard ================= */

  function _simDates() {
    const [t0, t1] = _simRange();
    const iso = (t) => new Date(t * 1000).toISOString().slice(0, 10);
    return [iso(t0), iso(t1)];
  }

  function _sourceWizard(initialVar) {
    // a background fetch finishing must not re-close this popup after
    // the user already closed it (Popup.close re-runs onClose, which
    // would clear a LATER wizard's station markers / pick mode)
    let wizardOpen = true;
    const popup = Popup.open({ title: "Download raw series",
      width: 700, onClose: () => { wizardOpen = false; _clearStations(); } });
    const box = popup.body;

    // pick a source, then tick one or more quantities to download from a
    // single station in one action (waterinfo has wind/water level/waves;
    // ERA5 provides wind only). Stations are discovered for the FIRST ticked
    // quantity; the others are pulled from that same station.
    const KIND_LIST = [["wind", "Wind (speed + direction)"],
      ["tide", "Water level"], ["wave", "Waves (Hs + Tp)"]];
    const selectedKinds = new Set([_variableById(initialVar || "wind_speed").kind]);
    const primaryKind = () => (selectedKinds.size ? [...selectedKinds][0] : "wind");
    let kind = primaryKind();   // station discovery + period probe use this

    const sourceSel = U.el("select", {});
    const rebuildSources = () => {
      U.clear(sourceSel);
      // ERA5 only makes sense while the selection is wind-only
      const onlyWind = selectedKinds.size === 1 && selectedKinds.has("wind");
      const sources = onlyWind ? ["waterinfo", "era5"] : ["waterinfo"];
      for (const s of sources) {
        sourceSel.append(U.el("option", { value: s },
          s === "era5" ? "ERA5 reanalysis (CDS)" : "waterinfo.rws.nl (measurements)"));
      }
    };
    rebuildSources();

    const kindBox = U.el("div", { class: "choice-col" });
    const rebuildKindBox = () => {
      U.clear(kindBox);
      const eraOnly = sourceSel.value === "era5";
      for (const [k, label] of KIND_LIST) {
        const cb = U.el("input", { type: "checkbox", id: `dl-k-${k}` });
        cb.checked = selectedKinds.has(k);
        cb.disabled = eraOnly && k !== "wind";
        cb.addEventListener("change", () => {
          if (cb.checked) selectedKinds.add(k); else selectedKinds.delete(k);
          kind = primaryKind();
          rebuildSources();
          _clearStations();
          _renderStations([], null);
          _syncPickBtn();
          syncCds();
        });
        kindBox.append(U.el("div", { class: "choice-row" }, cb,
          U.el("label", { for: `dl-k-${k}` }, label)));
      }
    };

    const cdsBox = U.el("div", { class: "muted", style: "font-size:12px" });
    const syncCds = async () => {
      U.clear(cdsBox);
      if (sourceSel.value !== "era5") return;
      cdsBox.append(U.el("div", { style: "margin-bottom:4px" },
        "ⓘ ERA5 requests wait in the shared Copernicus (CDS) server queue — "
        + "typically minutes, but it can take hours when the service is busy. "
        + "That wait is on Copernicus' side and cannot be sped up from here; "
        + "the download runs in the background and finished years are cached, "
        + "so a retry resumes where it left off."));
      const status = await Api.get("/api/conditions/cds").catch(() => null);
      if (status && status.configured) {
        const changeBtn = U.el("button", { class: "ghost", style: "font-size:11.5px;padding:1px 8px" },
          "Change…");
        changeBtn.addEventListener("click", () => _cdsKeyDialog(syncCds));
        const removeBtn = U.el("button", { class: "ghost danger-hover", style: "font-size:11.5px;padding:1px 8px" },
          "Remove");
        removeBtn.addEventListener("click", async () => {
          await Api.post("/api/conditions/cds_clear");
          U.toast("CDS key removed (~/.cdsapirc deleted)", "ok");
          syncCds();
        });
        cdsBox.append(
          U.el("span", {}, "✔ CDS API key configured "), changeBtn, " ", removeBtn);
      } else {
        const setupBtn = U.el("button", { class: "ghost" },
          U.icon("key", 13), " Set up CDS key…");
        setupBtn.classList.add("btn-ict");
        setupBtn.addEventListener("click", () => _cdsKeyDialog(syncCds));
        cdsBox.append(
          U.el("div", {}, `⚠ ${status ? status.reason : "CDS status unknown"}`),
          setupBtn);
      }
    };
    sourceSel.addEventListener("change", () => {
      rebuildKindBox();   // ERA5 disables non-wind quantities
      _clearStations();
      _renderStations([], null);
      _syncPickBtn();
      syncCds();
      _updateEra5Cache();
    });
    rebuildKindBox();
    syncCds();

    let selectedStation = null;
    const stationList = U.el("div", { class: "layer-tree", style: "max-height:180px;overflow-y:auto" });
    const selectedLine = U.el("div", { class: "muted", style: "font-size:12px" }, "no station selected");
    const periodLine = U.el("div", { class: "muted", style: "font-size:12px" });

    const findBtn = U.el("button", { class: "ghost btn-ict" },
      U.icon("search", 14), U.el("span", {}, "Find stations near grid"));
    findBtn.addEventListener("click", async () => {
      findBtn.disabled = true;
      try {
        const res = await Api.get(`/api/conditions/stations?source=${sourceSel.value}&kind=${kind}`);
        const payload = res.job ? await Api.waitJob(res.job) : res;
        _renderStations(payload.stations || [], select);
      } catch (err) {
        U.toast(err.message, "error");
      } finally {
        findBtn.disabled = false;
      }
    });

    const pickBtn = U.el("button", { class: "ghost btn-ict", style: "display:none" },
      U.icon("target", 14), U.el("span", {}, "Pick on map"));
    const _syncPickBtn = () => {
      pickBtn.style.display = stationMarkers.length ? "" : "none";
    };
    pickBtn.addEventListener("click", () => {
      if (!stationMarkers.length) return;
      popup.hide();
      U.toast("Click a station marker on the map (Esc to cancel)");
      _mapPickMode = (station) => {
        _mapPickMode = null;
        popup.show();
        if (station) select(station);
      };
      const cancel = (ev) => {
        if (ev.key === "Escape" && _mapPickMode) {
          window.removeEventListener("keydown", cancel);
          _mapPickMode(null);
        }
      };
      window.addEventListener("keydown", cancel);
    });

    const select = (station) => {
      selectedStation = station;
      selectedLine.textContent = `selected: ${station.name || station.id}` +
        (station.dist_km !== undefined ? ` (${station.dist_km.toFixed(0)} km)` : "");
      periodLine.textContent = "";
      for (const row of stationList.querySelectorAll(".lp-row")) {
        row.classList.toggle("selected", row.dataset.sid === station.id);
      }
      _highlightMarker(station.id);
      _updateEra5Cache();
    };

    let lastEraCells = [];   // era5 cells currently listed (for cache marks)

    function _renderStations(stations, onSelect) {
      U.clear(stationList);
      _clearStations();
      if (!stations.length) {
        lastEraCells = [];
        stationList.append(U.el("div", { class: "muted" }, "—"));
        _syncPickBtn();
        return;
      }
      for (const st of stations.slice(0, 20)) {
        const row = U.el("div", { class: "lp-row", style: "cursor:pointer", dataset: { sid: st.id } },
          U.el("span", { class: "lp-name" }, st.name || st.id),
          U.el("span", { class: "lp-mini" },
            st.dist_km !== undefined ? `${st.dist_km.toFixed(0)} km` : ""));
        row.addEventListener("click", () => onSelect && onSelect(st));
        row.addEventListener("mouseenter", () => _highlightMarker(st.id, true));
        row.addEventListener("mouseleave", () => _highlightMarker(selectedStation ? selectedStation.id : null));
        stationList.append(row);

        if (st.lon !== undefined && st.lat !== undefined && st.lon !== null) {
          _addStationMarker(st, () => {
            if (_mapPickMode) _mapPickMode(st);
            else if (onSelect) onSelect(st);
          });
        }
      }
      _syncPickBtn();
      if (sourceSel.value === "era5") {
        lastEraCells = stations.slice(0, 20);
        _markEraCells(onSelect);
      } else {
        lastEraCells = [];
      }
    }

    /* The ERA5 cache is PER CELL: mark each listed cell with its cached
     * years, and auto-select the best-cached cell so the download reuses
     * what is already on disk instead of silently starting a fresh cell. */
    async function _markEraCells(onSelect) {
      const d0 = date0.value.trim(), d1 = date1.value.trim();
      if (!lastEraCells.length || !d0 || !d1) return;
      let res;
      try {
        res = await Api.post("/api/conditions/era5_cells_cached", {
          cells: lastEraCells.map((s) => ({ id: s.id, lon: s.lon, lat: s.lat })),
          date0: d0, date1: d1,
        });
      } catch (e) { return; }
      let best = null;
      for (const c of res.cells || []) {
        if (!best || c.cached > best.cached) best = c;
        if (!c.cached) continue;
        const row = stationList.querySelector(`.lp-row[data-sid="${c.id}"]`);
        const mini = row && row.querySelector(".lp-mini");
        if (mini) {
          mini.textContent = `${c.cached}/${c.total} yrs cached`;
          mini.style.color = "var(--accent)";
          mini.style.fontWeight = "600";
        }
      }
      if (best && best.cached > 0 && !selectedStation && onSelect) {
        const st = lastEraCells.find((s) => s.id === best.id);
        if (st) {
          onSelect(st);
          selectedLine.textContent += " — auto-selected: this cell has cached data";
        }
      }
    }

    const periodBtn = U.el("button", { class: "ghost btn-ict" },
      U.icon("clock", 14), U.el("span", {}, "Check available period"));
    periodBtn.addEventListener("click", async () => {
      if (!selectedStation) { U.toast("Select a station first", "error"); return; }
      if (sourceSel.value === "era5") {
        periodLine.textContent = "ERA5: 1940 — present (global reanalysis)";
        return;
      }
      periodBtn.disabled = true;
      periodLine.textContent = "probing…";
      try {
        const res = await Api.post("/api/conditions/station_period",
          { station: selectedStation.id, kind });
        const period = await Api.waitJob(res.job, (j) => {
          periodLine.textContent = `probing… ${j.message || ""}`;
        });
        periodLine.textContent = period
          ? `data available ≈ ${period.from} — ${period.to || "?"} (${period.note})`
          : "no data found for this station";
      } catch (err) {
        periodLine.textContent = "";
        U.toast(err.message, "error");
      } finally {
        periodBtn.disabled = false;
      }
    });

    const [simFrom, simTo] = _simDates();
    const date0 = U.el("input", { type: "text", value: simFrom, title: "YYYY-MM-DD" });
    const date1 = U.el("input", { type: "text", value: simTo, title: "YYYY-MM-DD" });

    // ERA5: report which years are already cached on disk for the picked
    // cell + period, BEFORE the user starts the download
    const cacheLine = U.el("div", { class: "muted", style: "font-size:12px" });
    let cacheReq = 0;
    async function _updateEra5Cache() {
      cacheLine.textContent = "";
      if (sourceSel.value !== "era5" || !selectedStation) return;
      const d0 = date0.value.trim(), d1 = date1.value.trim();
      if (!d0 || !d1) return;
      const req = ++cacheReq;
      try {
        const res = await Api.post("/api/conditions/era5_cached",
          { lon: selectedStation.lon, lat: selectedStation.lat, date0: d0, date1: d1 });
        if (req !== cacheReq || !res.total) return;
        if (!res.missing.length) {
          cacheLine.textContent = `✔ all ${res.total} year(s) already on disk for this cell — `
            + "no CDS request needed, the series is assembled from the cache";
        } else if (res.cached.length) {
          cacheLine.textContent = `✔ ${res.cached.length} of ${res.total} years already on disk `
            + `for this cell — only ${res.missing.join(", ")} will be requested from CDS`;
        } else {
          cacheLine.textContent = `no cached ERA5 data for this cell yet — `
            + `${res.total} year(s) will be requested from CDS`;
        }
        // the cache is per cell: point at a neighbouring cell with more
        if (res.missing.length && lastEraCells.length) {
          const alt = await Api.post("/api/conditions/era5_cells_cached", {
            cells: lastEraCells.map((s) => ({ id: s.id, lon: s.lon, lat: s.lat })),
            date0: d0, date1: d1,
          });
          if (req !== cacheReq) return;
          const better = (alt.cells || [])
            .filter((c) => c.id !== selectedStation.id && c.cached > res.cached.length)
            .sort((a, b) => b.cached - a.cached)[0];
          if (better) {
            cacheLine.append(U.el("div", { style: "color:var(--accent)" },
              `tip: cell (${better.lat.toFixed(2)}N, ${better.lon.toFixed(2)}E) has `
              + `${better.cached} cached year(s) for this period — select that cell to reuse them`));
          }
        }
      } catch (e) { /* project closed / bad dates - just show nothing */ }
    }
    const _onDates = () => { _updateEra5Cache(); _markEraCells(null); };
    date0.addEventListener("change", _onDates);
    date1.addEventListener("change", _onDates);

    const progress = U.progressBar();
    const fetchBtn = U.el("button", { class: "primary btn-ict" },
      U.icon("download", 14), U.el("span", {}, "Download raw series"));
    const stopBtn = U.el("button", { class: "danger", style: "display:none" }, "■ Stop");
    stopBtn.title = "Stop the download — everything fetched so far is kept";
    const bgNote = U.el("div", { class: "muted", style: "font-size:11.5px" });
    fetchBtn.addEventListener("click", async () => {
      if (!selectedStation) { U.toast("Select a station/cell first", "error"); return; }
      if (!selectedKinds.size) { U.toast("Tick at least one quantity to download", "error"); return; }
      const body = { source: sourceSel.value, kinds: [...selectedKinds] };
      body.station = sourceSel.value === "waterinfo" ? selectedStation.id : selectedStation;
      body.station_name = selectedStation.name || selectedStation.id;
      if (date0.value.trim()) body.date0 = date0.value.trim();
      if (date1.value.trim()) body.date1 = date1.value.trim();
      try {
        fetchBtn.disabled = true;
        progress.start("starting download…");
        bgNote.textContent = "The download continues in the background — " +
          "you can close this window (progress + Stop stay visible in the Conditions tab).";
        const res = await Api.post("/api/conditions/fetch", body);
        stopBtn.style.display = "";
        stopBtn.disabled = false;
        stopBtn.onclick = () => {
          stopBtn.disabled = true;
          Api.post(`/api/job/cancel/${res.job}`).catch(() => {});
        };
        _watchFetch(primaryKind(), res.job, progress, (err) => {
          if (!wizardOpen) return;   // user closed the wizard meanwhile
          fetchBtn.disabled = false;
          stopBtn.style.display = "none";
          progress.done();
          _updateEra5Cache();        // partial years may now be cached
          // close on success, keep open on failure
          if (!err) popup.close();
        });
      } catch (err) {
        progress.done();
        U.toast(err.message, "error");
        fetchBtn.disabled = false;
      }
    });

    box.append(
      U.el("div", { class: "form-row" }, U.el("label", {}, "source"), sourceSel),
      U.el("span", { class: "fg-label" }, "Quantities to download"),
      kindBox,
      cdsBox,
      U.el("div", { class: "btn-row" }, findBtn, pickBtn),
      stationList,
      selectedLine,
      U.el("div", { class: "btn-row" }, periodBtn),
      periodLine,
      U.el("div", { class: "form-row" }, U.el("label", {}, "from"), date0),
      U.el("div", { class: "form-row" }, U.el("label", {}, "to"), date1),
      cacheLine,
      U.el("div", { class: "muted", style: "font-size:11.5px" },
        "The download is stored as a raw series — inspect and clean it, ",
        "then build the input file with Fill (resampling happens there)."),
      U.el("div", { class: "btn-row" }, fetchBtn, stopBtn),
      progress.el,
      bgNote,
    );
  }

  /* Track a fetch job; the section progress bar keeps updating even
   * when the wizard popup is closed. */
  function _watchFetch(kind, jobId, popupProgress, onDone) {
    const bar = U.progressBar();
    bar.start("downloading…");
    // stop button next to the tab's progress row, so a long download can
    // be cancelled even after the wizard popup was closed
    const stop = U.el("button", {
      class: "ghost danger-hover",
      style: "font-size:11.5px;padding:1px 8px;flex:none;align-self:center",
      title: "Stop this download — everything fetched so far is kept",
    }, "■ Stop");
    stop.addEventListener("click", () => {
      stop.disabled = true;
      Api.post(`/api/job/cancel/${jobId}`).catch(() => {});
    });
    const row = U.el("div", { style: "display:flex;gap:8px;align-items:center" },
      U.el("div", { style: "flex:1;min-width:0" }, bar.el), stop);
    activeFetch.set(jobId, { kind, progressBar: bar, row });
    _build();
    Api.waitJob(jobId, (j) => {
      bar.update(j);
      if (popupProgress) popupProgress.update(j);
    }).then(async (out) => {
      activeFetch.delete(jobId);
      // one job can now yield several raw series (multiple quantities)
      const entries = out.entries || (out.entry ? [out.entry] : []);
      if (entries.length === 1) {
        U.toast(`Downloaded ${entries[0].label} (${entries[0].rows} rows)`, "ok");
      } else if (entries.length) {
        U.toast(`Downloaded ${entries.length} raw series`, "ok");
      }
      if (out.errors && out.errors.length) {
        U.toast(`Some quantities failed: ${out.errors.join("; ")}`, "error");
      }
      if (onDone) onDone();
      await _refresh();
    }).catch((err) => {
      activeFetch.delete(jobId);
      U.toast(err.message, "error");
      if (onDone) onDone(err);
      _build();
    });
  }

  function _cdsKeyDialog(onDone) {
    const popup = Popup.open({ title: "Set up your CDS API key (ERA5)", width: 520 });
    const key = U.el("input", { type: "text", placeholder: "paste your Personal Access Token" });
    const saveBtn = U.el("button", { class: "primary" }, "Save key");
    const note = U.el("div", { class: "muted", style: "font-size:12px" });
    saveBtn.addEventListener("click", async () => {
      try {
        const res = await Api.post("/api/conditions/cds_key", { key: key.value });
        if (res.ok) {
          U.toast("CDS key stored locally (~/.cdsapirc)", "ok");
          popup.close();
          if (onDone) onDone();
        } else {
          note.textContent = res.reason || "key stored, but cdsapi still not ready";
        }
      } catch (err) {
        note.textContent = err.message;
      }
    });
    popup.body.append(
      U.el("div", { style: "font-size:13px;line-height:1.5" },
        "ERA5 downloads need a free personal API key from the Copernicus ",
        "Climate Data Store:", U.el("br"),
        "1. Create an account at ",
        U.el("a", { href: "https://cds.climate.copernicus.eu", target: "_blank" },
          "cds.climate.copernicus.eu"), U.el("br"),
        "2. Accept the ERA5 licence (once, on the dataset page)", U.el("br"),
        "3. Copy the API token from your profile page and paste it below.", U.el("br"),
        "The key is stored only on this computer (~/.cdsapirc), never in the project."),
      U.el("div", { class: "form-row", style: "margin-top:10px" }, key),
      U.el("div", { class: "btn-row" }, saveBtn),
      note,
    );
  }

  /* ---- station markers (dots, name on hover) ---- */

  let _mapPickMode = null;

  function _addStationMarker(station, onClick) {
    const dot = U.el("div", { class: "station-dot", dataset: { sid: station.id } },
      U.el("span", { class: "station-tip" }, station.name || station.id));
    dot.addEventListener("click", (ev) => { ev.stopPropagation(); onClick(); });
    const marker = new maplibregl.Marker({ element: dot, anchor: "center" })
      .setLngLat([station.lon, station.lat])
      .addTo(MapView.instance());
    stationMarkers.push(marker);
  }

  function _highlightMarker(stationId, hover = false) {
    for (const marker of stationMarkers) {
      const el = marker.getElement();
      el.classList.toggle("selected", el.dataset.sid === stationId);
    }
  }

  function _clearStations() {
    for (const marker of stationMarkers) marker.remove();
    stationMarkers = [];
    _mapPickMode = null;
  }

  async function _reloadConfig(kind) {
    const cfg = await Api.get("/api/config");
    App.state.config = cfg.values;
    App.emit("config-changed", overview.kinds[kind].config_key);
  }

  return { init, windroseSources };
})();
