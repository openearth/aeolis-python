/* Conditions tab: wind, water level and wave boundary conditions.
 *
 * Each series has one "Generate…" button opening a wizard popup:
 *  - Synthetic: profile forms with a live preview chart (computed
 *    server-side, nothing written until Save).
 *  - Measured / modelled: waterinfo.rws.nl stations or ERA5 wind
 *    (personal CDS key, set up in-app). Stations appear as dot markers
 *    on the map (name on hover, hover-synced with the list); the popup
 *    steps aside while picking on the map. Long periods are fetched in
 *    chunks and can be resampled to hourly/daily means.
 */
"use strict";

const ConditionsTab = (() => {

  const KIND_TITLES = { wind: "Wind", tide: "Water levels", wave: "Waves" };
  const SERIES_COLORS = ["#0f766e", "#b4423b", "#2f7fe6"];

  let overview = null;
  let stationMarkers = [];

  function init() {
    Tabs.register("conditions", { enter: _refresh, leave: _clearStations });
    App.on("project", _refresh);
  }

  async function _refresh() {
    if (!App.state.project) return;
    try {
      overview = await Api.get("/api/conditions");
    } catch (err) {
      U.toast(err.message, "error");
      return;
    }
    _build();
    _plotExisting();
  }

  /* ================= plotting (graphs panel) ================= */

  function _plotExisting() {
    for (const [kind, info] of Object.entries(overview.kinds)) {
      if (info.series) _plotSeries(kind, info);
    }
  }

  function _plotSeries(kind, info) {
    const s = info.series;
    if (kind === "wind") {
      // wind speed and direction as two separate graphs
      Graphs.add("cond-wind-speed", {
        title: `Wind speed — ${info.file}`,
        data: [s.t_epoch, s.columns[0]],
        series: [{}, { label: "speed [m/s]", stroke: SERIES_COLORS[0], width: 1.5 }],
        timeBased: true,
      });
      Graphs.add("cond-wind-dir", {
        title: `Wind direction — ${info.file}`,
        data: [s.t_epoch, s.columns[1]],
        series: [{}, { label: "direction [°]", stroke: SERIES_COLORS[1], width: 1.2 }],
        timeBased: true,
        scales: { x: { time: false }, y: { range: [0, 360] } },
      });
    } else {
      const series = [{}];
      info.labels.forEach((label, i) => {
        series.push({ label, stroke: SERIES_COLORS[i % SERIES_COLORS.length], width: 1.5,
          scale: i === 0 ? "y" : "y2" });
      });
      Graphs.add(`cond-${kind}`, {
        title: `${KIND_TITLES[kind]} — ${info.file}`,
        data: [s.t_epoch, ...s.columns],
        series,
        timeBased: true,
        axes: [
          {}, { size: 52, scale: "y" },
          { size: 52, scale: "y2", side: 1, grid: { show: false } },
        ],
        scales: { x: { time: false }, y: {}, y2: {} },
      });
    }
    Playbar.setSource(`cond-${kind}`, s.t0_epoch, s.t1_epoch);
  }

  /* ================= panel ================= */

  function _build() {
    const panel = document.getElementById("conditions-panel");
    U.clear(panel);

    panel.append(U.el("div", { class: "muted", style: "font-size:12px" },
      `refdate ${overview.refdate} — simulation ${U.fmtDuration(overview.tstop - overview.tstart)}`));

    for (const kind of ["wind", "tide", "wave"]) {
      const info = overview.kinds[kind];
      const status = info.exists
        ? U.el("div", { class: "muted" },
          `✔ ${info.file} (${info.series ? info.series.n + " rows" : "unreadable"})`)
        : U.el("div", { class: "muted" }, "not configured yet");
      const genBtn = U.el("button", { class: "primary" }, "Generate…");
      genBtn.addEventListener("click", () => _wizard(kind));
      panel.append(U.el("div", { class: "form-group" },
        U.el("span", { class: "fg-label" }, KIND_TITLES[kind]),
        status,
        U.el("div", { class: "btn-row" }, genBtn)));
    }
  }

  /* ================= wizard ================= */

  const PROFILE_PARAMS = {
    constant: [["value", "value"]],
    blocks: [["values", "values (comma sep.)"], ["block_duration", "block dur. [h]"]],
    harmonic: [["mean", "mean"], ["amplitude", "amplitude"], ["period", "period [h]"], ["phase", "phase [deg]"]],
    linear: [["start", "start"], ["end", "end"]],
    rotational: [["start", "start [deg]"], ["rate", "rate [deg/h]"]],
  };
  const HOUR_PARAMS = new Set(["block_duration", "period"]);

  function _profileForm(label, types, defaults = {}, onInput = null) {
    const typeSel = U.el("select", {}, ...types.map((t) => U.el("option", { value: t }, t)));
    const paramBox = U.el("div");
    const inputs = {};

    const renderParams = () => {
      U.clear(paramBox);
      for (const [key, title] of PROFILE_PARAMS[typeSel.value]) {
        const input = U.el("input", { type: "text", value: defaults[key] ?? "" });
        if (onInput) input.addEventListener("input", onInput);
        inputs[key] = input;
        paramBox.append(U.el("div", { class: "form-row" },
          U.el("label", {}, title), input));
      }
      if (onInput) onInput();
    };
    typeSel.addEventListener("change", renderParams);
    renderParams();

    return {
      el: U.el("div", { class: "form-group" },
        U.el("span", { class: "fg-label" }, label),
        U.el("div", { class: "form-row" }, U.el("label", {}, "type"), typeSel),
        paramBox),
      spec: () => {
        const spec = { type: typeSel.value };
        for (const [key] of PROFILE_PARAMS[typeSel.value]) {
          const raw = inputs[key].value.trim();
          if (raw === "") continue;
          if (key === "values") {
            spec[key] = raw.split(/[\s,]+/).map(Number).filter(Number.isFinite);
          } else {
            let v = Number(raw);
            if (HOUR_PARAMS.has(key)) v *= 3600;
            spec[key] = v;
          }
        }
        return spec;
      },
    };
  }

  function _wizard(kind) {
    const popup = Popup.open({ title: `Generate ${KIND_TITLES[kind].toLowerCase()}`, width: 640,
      onClose: _clearStations });

    const typeSyn = U.el("input", { type: "radio", name: "cw-type", id: "cw-syn", checked: "" });
    const typeMeas = U.el("input", { type: "radio", name: "cw-type", id: "cw-meas" });
    const synBox = U.el("div");
    const measBox = U.el("div", { style: "display:none" });
    const syncType = () => {
      synBox.style.display = typeSyn.checked ? "" : "none";
      measBox.style.display = typeMeas.checked ? "" : "none";
    };
    typeSyn.addEventListener("change", syncType);
    typeMeas.addEventListener("change", syncType);

    popup.body.append(
      U.el("div", { class: "form-row" },
        typeSyn, U.el("label", { for: "cw-syn" }, "Synthetic"),
        typeMeas, U.el("label", { for: "cw-meas" }, "Measured / modelled")),
      synBox, measBox,
    );

    _buildSynthetic(kind, synBox, popup);
    _buildMeasured(kind, measBox, popup);
  }

  /* ---- synthetic branch with live preview ---- */

  function _buildSynthetic(kind, box, popup) {
    const forms = {};
    let previewPlot = null;

    const updatePreview = U.debounce(async () => {
      const body = { kind, dt: Number(dtInput.value) * 3600 || 3600 };
      for (const [name, form] of Object.entries(forms)) body[name] = form.spec();
      try {
        const res = await Api.post("/api/conditions/preview", body);
        const s = res.series;
        const data = [s.t_epoch, ...s.columns];
        if (previewPlot) previewPlot.destroy();
        previewPlot = new uPlot({
          width: 560, height: 150,
          series: [{}, ...res.labels.map((label, i) => ({
            label, stroke: SERIES_COLORS[i % SERIES_COLORS.length], width: 1.4,
            scale: i === 0 ? "y" : "y2",
          }))],
          axes: [{ values: (u, t) => t.map((v) => U.fmtDate(v).slice(5, 16)) },
            { size: 50, scale: "y" }, { size: 50, scale: "y2", side: 1, grid: { show: false } }],
          scales: { x: { time: false }, y: {}, y2: {} },
          legend: { show: true },
        }, data, previewEl);
      } catch (err) {
        console.warn("preview failed", err.message);
      }
    }, 350);

    if (kind === "wind") {
      forms.speed = _profileForm("Wind speed [m/s]",
        ["constant", "blocks", "harmonic", "linear"], { value: 10 }, updatePreview);
      forms.direction = _profileForm("Wind direction [deg]",
        ["constant", "blocks", "harmonic", "rotational"], { value: 270 }, updatePreview);
    } else if (kind === "tide") {
      forms.level = _profileForm("Water level [m]",
        ["constant", "harmonic", "blocks", "linear"],
        { mean: 0, amplitude: 1, period: 12.42, value: 0 }, updatePreview);
    } else {
      forms.hs = _profileForm("Wave height Hs [m]",
        ["constant", "blocks", "harmonic", "linear"], { value: 1 }, updatePreview);
      forms.tp = _profileForm("Wave period Tp [s]",
        ["constant", "blocks", "harmonic", "linear"], { value: 6 }, updatePreview);
    }

    const dtInput = U.el("input", { type: "text", value: "1", style: "width:70px" });
    dtInput.addEventListener("input", updatePreview);

    const previewEl = U.el("div", { class: "wizard-preview" });

    const saveBtn = U.el("button", { class: "primary" }, "Generate & save");
    saveBtn.addEventListener("click", async () => {
      const body = { kind, dt: Number(dtInput.value) * 3600 || 3600 };
      for (const [name, form] of Object.entries(forms)) body[name] = form.spec();
      try {
        saveBtn.disabled = true;
        const res = await Api.post("/api/conditions/synthetic", body);
        U.toast(`Wrote ${res.file} (${res.rows} rows)`, "ok");
        await _reloadConfig(kind);
        popup.close();
        _refresh();
      } catch (err) {
        U.toast(err.message, "error");
        saveBtn.disabled = false;
      }
    });

    for (const form of Object.values(forms)) box.append(form.el);
    box.append(
      U.el("div", { class: "form-row" }, U.el("label", {}, "output step [h]"), dtInput),
      U.el("span", { class: "fg-label" }, "Preview"),
      previewEl,
      U.el("div", { class: "btn-row" }, saveBtn),
    );
    updatePreview();
  }

  /* ---- measured / modelled branch ---- */

  function _buildMeasured(kind, box, popup) {
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
        cdsBox.append("✔ CDS API key configured");
      } else {
        const setupBtn = U.el("button", { class: "ghost" }, "Set up CDS key…");
        setupBtn.addEventListener("click", () => _cdsKeyDialog(syncCds));
        cdsBox.append(
          U.el("div", {}, `⚠ ${status ? status.reason : "CDS status unknown"}`),
          setupBtn);
      }
    };
    sourceSel.addEventListener("change", () => { _clearStations(); _renderStations([], null); syncCds(); });
    syncCds();

    let selectedStation = null;
    const stationList = U.el("div", { class: "layer-tree", style: "max-height:180px;overflow-y:auto" });
    const selectedLine = U.el("div", { class: "muted", style: "font-size:12px" }, "no station selected");
    const periodLine = U.el("div", { class: "muted", style: "font-size:12px" });

    const findBtn = U.el("button", { class: "ghost" }, "Find stations near grid");
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

    const pickBtn = U.el("button", { class: "ghost" }, "Pick on map");
    pickBtn.addEventListener("click", () => {
      if (!stationMarkers.length) { U.toast("Find stations first", "error"); return; }
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
    }

    const periodBtn = U.el("button", { class: "ghost" }, "Check available period");
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

    const date0 = U.el("input", { type: "text", placeholder: "YYYY-MM-DD (default: sim start)" });
    const date1 = U.el("input", { type: "text", placeholder: "YYYY-MM-DD (default: sim end)" });
    const resample = U.el("select", {},
      U.el("option", { value: "" }, "raw (as measured)"),
      U.el("option", { value: "hour" }, "hourly means"),
      U.el("option", { value: "day" }, "daily means"));

    const progress = U.el("div", { class: "muted", style: "font-size:12px" });
    const fetchBtn = U.el("button", { class: "primary" }, "Download & convert");
    fetchBtn.addEventListener("click", async () => {
      if (!selectedStation) { U.toast("Select a station/cell first", "error"); return; }
      const body = { source: sourceSel.value, kind };
      body.station = sourceSel.value === "waterinfo" ? selectedStation.id : selectedStation;
      if (date0.value.trim()) body.date0 = date0.value.trim();
      if (date1.value.trim()) body.date1 = date1.value.trim();
      if (resample.value) body.resample = resample.value;
      try {
        fetchBtn.disabled = true;
        const res = await Api.post("/api/conditions/fetch", body);
        const out = await Api.waitJob(res.job, (j) => { progress.textContent = j.message || ""; });
        progress.textContent = "";
        U.toast(`Wrote ${out.file} (${out.rows} rows)`, "ok");
        await _reloadConfig(kind);
        popup.close();
        _refresh();
      } catch (err) {
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
      U.el("div", { class: "form-row" }, U.el("label", {}, "resample"), resample),
      U.el("div", { class: "btn-row" }, fetchBtn),
      progress,
    );
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
