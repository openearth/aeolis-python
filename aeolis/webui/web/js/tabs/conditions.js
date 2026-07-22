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

  async function _openWindroseRaw(entry) {
    const cols = _windroseCols(entry.labels);
    if (!cols) return;
    try {
      const [, magArr] = await _loadRawColumn(entry, cols.mag);
      const [, dirArr] = await _loadRawColumn(entry, cols.dir);
      Windrose.open({ title: entry.label, magnitude: magArr, direction: dirArr,
        magName: cols.magName, magUnit: cols.magUnit });
    } catch (err) { U.toast(err.message, "error"); }
  }

  function _openWindroseInput(kind) {
    const info = overview.kinds[kind];
    const cols = _windroseCols(info && info.labels);
    if (!info || !info.series || !cols) return;
    Windrose.open({
      title: info.file || kind,
      magnitude: info.series.columns[cols.mag],
      direction: info.series.columns[cols.dir],
      magName: cols.magName, magUnit: cols.magUnit,
    });
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
          range: (entry.t0_epoch !== null && entry.t1_epoch !== null)
            ? [entry.t0_epoch, entry.t1_epoch] : null,
        });
      });
    }
  }

  function _rawSeriesKey(entry, col) { return `rawcond-${entry.id}-${col}`; }

  async function _loadRawColumn(entry, col) {
    if (!rawSeriesCache.has(entry.id)) {
      // don't cache failures: a transient error would stick forever
      const promise = Api.get(`/api/conditions/raw_series?id=${entry.id}`)
        .catch((err) => { rawSeriesCache.delete(entry.id); throw err; });
      rawSeriesCache.set(entry.id, promise);
    }
    const res = await rawSeriesCache.get(entry.id);
    return [res.series.t_epoch, res.series.columns[col] || []];
  }

  /* ================= panel ================= */

  function _build() {
    const panel = document.getElementById("conditions-panel");
    U.clear(panel);

    panel.append(U.el("div", { class: "muted", style: "font-size:12px" },
      `refdate ${overview.refdate} — simulation ${U.fmtDuration(overview.tstop - overview.tstart)}`));

    for (const kind of ["wind", "tide", "wave"]) {
      const info = overview.kinds[kind];
      const kindRaw = rawEntries.filter((e) => e.kind === kind);
      const section = U.section(KIND_TITLES[kind],
        { count: kindRaw.length || null });

      const status = info.exists
        ? U.el("div", { class: "muted" },
          `✔ ${info.file} (${info.series ? info.series.n + " rows" : "unreadable"})`)
        : U.el("div", { class: "muted" }, `${info.file || "input file"} not written yet`);

      const tools = U.el("div", { class: "tbtn-row" },
        U.tbtn("wand", "Generate", {
          primary: !info.exists && !kindRaw.length,
          title: "Generate a synthetic series (written directly to the input file)",
          onclick: () => _synthWizard(kind),
        }),
        U.tbtn("download", "Download", {
          title: "Download measured/reanalysis data as a raw series",
          onclick: () => _sourceWizard(kind),
        }));
      // the input file has both magnitude + direction → offer a windrose
      if (info.series && _windroseCols(info.labels)) {
        tools.append(U.tbtn("compass", "Windrose", {
          title: `Windrose of ${info.file || kind}`,
          onclick: () => _openWindroseInput(kind),
        }));
      }
      section.body.append(status, tools);

      // active background downloads (possibly several per kind)
      for (const fetching of activeFetch.values()) {
        if (fetching.kind === kind) section.body.append(fetching.progressBar.el);
      }

      // raw series cards
      if (kindRaw.length) {
        section.body.append(U.el("span", { class: "fg-label" }, "Raw series"));
        section.body.append(_rawList(kind, kindRaw));
      }
      panel.append(section.wrap);
    }
  }

  function _rawList(kind, entries) {
    const wrap = U.el("div", {});

    // top toolbar: Remove acts on the multi-selection (mirrors Domain tab)
    const selInKind = entries.filter((e) => selectedRaw.has(e.id));
    const removeBtn = U.tbtn("trash", "Remove", {
      title: selInKind.length
        ? `Remove ${selInKind.length} selected series…`
        : "Select series first (click, Ctrl+click, Shift+click for a range)",
      onclick: () => _removeSelectedRaw(kind),
    });
    removeBtn.disabled = !selInKind.length;
    wrap.append(U.el("div", { class: "tbtn-row", style: "margin:2px 0 6px" }, removeBtn));

    const list = U.el("div", { class: "obj-list" });
    const ids = entries.map((e) => e.id);
    for (const entry of entries) {
      const name = U.el("span", { class: "lp-name", title: "Double-click to rename" },
        entry.label);
      name.addEventListener("dblclick", () => _renameRaw(entry, name));

      const meta = `${entry.rows || 0} rows` +
        (entry.nan ? ` · ${entry.nan} NaN` : "");

      const actions = U.el("span", { class: "obj-actions" });
      actions.append(U.miniBtn("chart", "Show in the graph panel",
        () => Graphs.select(_rawSeriesKey(entry, 0))));
      if (_windroseCols(entry.labels)) {
        actions.append(U.miniBtn("compass", "Windrose (magnitude + direction)",
          () => _openWindroseRaw(entry)));
      }
      actions.append(U.miniBtn("copy", "Duplicate this series", () => _duplicateRaw(entry)));
      actions.append(U.miniBtn("modify", "Modify… (clean NaN, crop, resample, …)",
        () => _rawModifyWizard(entry)));
      actions.append(U.miniBtn("check",
        `Use as ${overview.kinds[kind].file || "the input file"} (converts to seconds since refdate)`,
        () => _applyRaw(entry)));

      const card = U.el("div", {
        class: `obj-card ${selectedRaw.has(entry.id) ? "selected" : ""}`,
      },
        name,
        entry.nan
          ? U.el("span", { class: "warn-icon", title: `${entry.nan} NaN value(s) — clean before applying` }, "⚠")
          : null,
        U.el("span", { class: "lp-mini" }, meta),
        actions);
      // click = select (Ctrl toggles, Shift selects a range)
      card.addEventListener("click", (ev) => {
        if (ev.target.closest("button, .eye, input")) return;
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

  function _removeSelectedRaw(kind) {
    const entries = rawEntries.filter((e) => selectedRaw.has(e.id) && e.kind === kind);
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
    ["resample", "resample to interval means"],
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

  function _synthWizard(kind) {
    const popup = Popup.open({ title: `Generate ${KIND_TITLES[kind].toLowerCase()}`, width: 700 });
    const forms = {};
    const previewPlots = [];
    const previewEl = U.el("div", { class: "wizard-preview" });
    const aliasNote = U.el("div", { class: "muted", style: "font-size:11.5px;color:#b45309" });

    const renderPreview = (res) => {
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
          height: res.labels.length > 1 ? 130 : 170,
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

    const collectBody = () => {
      const body = { kind, dt: Number(dtInput.value) * 3600 || 3600 };
      for (const [name, form] of Object.entries(forms)) body[name] = form.spec();
      return body;
    };

    const checkAliasing = (body) => {
      const dt = body.dt;
      let worst = null;
      for (const key of Object.keys(forms)) {
        for (const seg of (body[key].segments || [])) {
          if (seg.type === "harmonic" && seg.period && seg.period < 4 * dt) {
            worst = seg.period;
          }
        }
      }
      aliasNote.textContent = worst !== null
        ? `⚠ the output step (${U.fmtNum(dt / 3600, 3)} h) is coarse for a ` +
          `${U.fmtNum(worst / 3600, 3)} h period — the sine will look jagged; use a smaller output step`
        : "";
    };

    const updatePreview = U.debounce(async () => {
      const body = collectBody();
      checkAliasing(body);
      try {
        renderPreview(await Api.post("/api/conditions/preview", body));
      } catch (err) {
        console.warn("preview failed", err.message);
      }
    }, 350);

    if (kind === "wind") {
      forms.speed = _segmentTable("Wind speed [m/s]", { value: 10 }, updatePreview);
      forms.direction = _segmentTable("Wind direction [deg]", { value: 270 }, updatePreview);
    } else if (kind === "tide") {
      forms.level = _segmentTable("Water level [m]",
        { type: "harmonic", mean: 0, amplitude: 1, period: 12.42 }, updatePreview);
    } else {
      forms.hs = _segmentTable("Wave height Hs [m]", { value: 1 }, updatePreview);
      forms.tp = _segmentTable("Wave period Tp [s]", { value: 6 }, updatePreview);
    }

    const dtInput = U.el("input", { type: "text", value: "1", style: "width:70px" });
    dtInput.addEventListener("input", updatePreview);

    const saveBtn = U.el("button", { class: "primary" }, "Generate & save");
    saveBtn.addEventListener("click", async () => {
      try {
        saveBtn.disabled = true;
        const res = await Api.post("/api/conditions/synthetic", collectBody());
        U.toast(`Wrote ${res.file} (${res.rows} rows)`, "ok");
        await _reloadConfig(kind);
        popup.close();
        _refresh();
      } catch (err) {
        U.toast(err.message, "error");
        saveBtn.disabled = false;
      }
    });

    for (const form of Object.values(forms)) popup.body.append(form.el);
    popup.body.append(
      U.el("div", { class: "form-row" }, U.el("label", {}, "output step [h]"), dtInput),
      aliasNote,
      U.el("span", { class: "fg-label" }, "Preview (incl. repetition over the simulation)"),
      previewEl,
      U.el("div", { class: "btn-row" }, saveBtn),
    );
    updatePreview();
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

  function _sourceWizard(kind) {
    // a background fetch finishing must not re-close this popup after
    // the user already closed it (Popup.close re-runs onClose, which
    // would clear a LATER wizard's station markers / pick mode)
    let wizardOpen = true;
    const popup = Popup.open({ title: `Download ${KIND_TITLES[kind].toLowerCase()} (raw series)`,
      width: 700, onClose: () => { wizardOpen = false; _clearStations(); } });
    const box = popup.body;

    const sources = kind === "wind" ? ["waterinfo", "era5"] : ["waterinfo"];
    const sourceSel = U.el("select", {},
      ...sources.map((s) => U.el("option", { value: s },
        s === "era5" ? "ERA5 reanalysis (CDS)" : "waterinfo.rws.nl (measurements)")));

    const cdsBox = U.el("div", { class: "muted", style: "font-size:12px" });
    const syncCds = async () => {
      U.clear(cdsBox);
      if (sourceSel.value !== "era5") return;
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
      _clearStations();
      _renderStations([], null);
      _syncPickBtn();
      syncCds();
    });
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
      U.icon("eye", 14), U.el("span", {}, "Pick on map"));
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
    };

    function _renderStations(stations, onSelect) {
      U.clear(stationList);
      _clearStations();
      if (!stations.length) {
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

    // resample: off by default; on -> custom interval value + unit
    const resampleCb = U.el("input", { type: "checkbox", id: "rs-on" });
    const resampleVal = U.el("input", { type: "text", value: "1", style: "width:60px", disabled: "" });
    const resampleUnit = U.el("select", { disabled: "" },
      U.el("option", { value: 60 }, "minutes"),
      U.el("option", { value: 3600, selected: "" }, "hours"),
      U.el("option", { value: 86400 }, "days"));
    resampleCb.addEventListener("change", () => {
      resampleVal.disabled = !resampleCb.checked;
      resampleUnit.disabled = !resampleCb.checked;
    });

    const progress = U.progressBar();
    const fetchBtn = U.el("button", { class: "primary btn-ict" },
      U.icon("download", 14), U.el("span", {}, "Download raw series"));
    const bgNote = U.el("div", { class: "muted", style: "font-size:11.5px" });
    fetchBtn.addEventListener("click", async () => {
      if (!selectedStation) { U.toast("Select a station/cell first", "error"); return; }
      const body = { source: sourceSel.value, kind };
      body.station = sourceSel.value === "waterinfo" ? selectedStation.id : selectedStation;
      body.station_name = selectedStation.name || selectedStation.id;
      if (date0.value.trim()) body.date0 = date0.value.trim();
      if (date1.value.trim()) body.date1 = date1.value.trim();
      if (resampleCb.checked) {
        const width = Number(resampleVal.value) * Number(resampleUnit.value);
        if (Number.isFinite(width) && width > 0) body.resample = width;
      }
      try {
        fetchBtn.disabled = true;
        progress.start("starting download…");
        bgNote.textContent = "The download continues in the background — " +
          "you can close this window (progress stays visible in the Conditions tab).";
        const res = await Api.post("/api/conditions/fetch", body);
        _watchFetch(kind, res.job, progress, (err) => {
          if (!wizardOpen) return;   // user closed the wizard meanwhile
          fetchBtn.disabled = false;
          progress.done();
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
      cdsBox,
      U.el("div", { class: "btn-row" }, findBtn, pickBtn),
      stationList,
      selectedLine,
      U.el("div", { class: "btn-row" }, periodBtn),
      periodLine,
      U.el("div", { class: "form-row" }, U.el("label", {}, "from"), date0),
      U.el("div", { class: "form-row" }, U.el("label", {}, "to"), date1),
      U.el("div", { class: "form-row" },
        U.el("label", { for: "rs-on" }, "resample to interval means"),
        resampleCb, resampleVal, resampleUnit),
      U.el("div", { class: "muted", style: "font-size:11.5px" },
        "The download is stored as a raw series first — inspect and clean it, ",
        "then apply it to the input file with the ✓ button."),
      U.el("div", { class: "btn-row" }, fetchBtn),
      progress.el,
      bgNote,
    );
  }

  /* Track a fetch job; the section progress bar keeps updating even
   * when the wizard popup is closed. */
  function _watchFetch(kind, jobId, popupProgress, onDone) {
    const bar = U.progressBar();
    bar.start("downloading…");
    activeFetch.set(jobId, { kind, progressBar: bar });
    _build();
    Api.waitJob(jobId, (j) => {
      bar.update(j);
      if (popupProgress) popupProgress.update(j);
    }).then(async (out) => {
      activeFetch.delete(jobId);
      U.toast(`Downloaded ${out.entry.label} (${out.entry.rows} rows)`, "ok");
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

  return { init };
})();
