/* Conditions tab: wind, water level and wave boundary conditions.
 *
 * Three ways to fill each series: synthetic generators (constant /
 * blocks / harmonic / linear / rotational), measured data from
 * waterinfo.rws.nl stations, or ERA5 reanalysis wind (CDS API). All
 * write the AeoLiS text files and update the config; series preview in
 * the graphs panel below the map.
 */
"use strict";

const ConditionsTab = (() => {

  const KIND_TITLES = { wind: "Wind", tide: "Water levels", wave: "Waves" };
  const SERIES_COLORS = ["#0f766e", "#b4423b", "#2f7fe6"];

  let overview = null;
  let stationMarkers = [];
  let selectedStation = null;

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

  /* ================= plotting ================= */

  function _plotExisting() {
    for (const [kind, info] of Object.entries(overview.kinds)) {
      if (info.series) _plotSeries(kind, info);
    }
  }

  function _plotSeries(kind, info) {
    const s = info.series;
    const data = [s.t_epoch, ...s.columns];
    const series = [{}];
    info.labels.forEach((label, i) => {
      series.push({ label, stroke: SERIES_COLORS[i % SERIES_COLORS.length], width: 1.5,
        scale: i === 0 ? "y" : "y2" });
    });
    Graphs.add(`cond-${kind}`, {
      title: `${KIND_TITLES[kind]} — ${info.file}`,
      height: 150,
      data,
      series,
      timeBased: true,
      axes: [
        { values: (u, ticks) => ticks.map((t) => U.fmtDate(t).slice(5, 16)) },
        { size: 52, scale: "y" },
        { size: 52, scale: "y2", side: 1, grid: { show: false } },
      ],
      scales: { x: { time: false }, y: {}, y2: {} },
    });
    Playbar.setSource(`cond-${kind}`, s.t0_epoch, s.t1_epoch);
  }

  /* ================= UI ================= */

  function _build() {
    const panel = document.getElementById("conditions-panel");
    U.clear(panel);

    panel.append(U.el("div", { class: "muted", style: "font-size:12px" },
      `refdate ${overview.refdate} — simulation ${U.fmtDuration(overview.tstop - overview.tstart)} `,
      `(t = ${overview.tstart} … ${overview.tstop} s)`));

    for (const kind of ["wind", "tide", "wave"]) {
      panel.append(_kindSection(kind, overview.kinds[kind]));
    }
  }

  function _kindSection(kind, info) {
    const status = info.exists
      ? U.el("div", { class: "muted" }, `current: ${info.file} (${info.series ? info.series.n + " rows" : "unreadable"})`)
      : U.el("div", { class: "muted" }, "no file configured yet");

    const body = U.el("div", { class: "section-body" },
      status,
      _syntheticForm(kind),
      _measuredForm(kind),
    );
    const head = U.el("header", {}, U.el("span", { class: "caret" }, "▾"),
      KIND_TITLES[kind]);
    const wrap = U.el("div", { class: `section ${kind === "wind" ? "" : "collapsed"}` }, head, body);
    head.addEventListener("click", () => wrap.classList.toggle("collapsed"));
    return wrap;
  }

  /* ---- synthetic ---- */

  const PROFILE_PARAMS = {
    constant: [["value", "value"]],
    blocks: [["values", "values (comma sep.)"], ["block_duration", "block dur. [h]"]],
    harmonic: [["mean", "mean"], ["amplitude", "amplitude"], ["period", "period [h]"], ["phase", "phase [deg]"]],
    linear: [["start", "start"], ["end", "end"]],
    rotational: [["start", "start [deg]"], ["rate", "rate [deg/h]"]],
  };
  const HOUR_PARAMS = new Set(["block_duration", "period"]);

  function _profileForm(label, types, defaults = {}) {
    const typeSel = U.el("select", {}, ...types.map((t) => U.el("option", { value: t }, t)));
    const paramBox = U.el("div");
    const inputs = {};

    const renderParams = () => {
      U.clear(paramBox);
      for (const [key, title] of PROFILE_PARAMS[typeSel.value]) {
        const input = U.el("input", { type: "text", value: defaults[key] ?? "" });
        inputs[key] = input;
        paramBox.append(U.el("div", { class: "form-row" },
          U.el("label", {}, title), input));
      }
    };
    typeSel.addEventListener("change", renderParams);
    renderParams();

    const el = U.el("div", { class: "form-group" },
      U.el("span", { class: "fg-label" }, label),
      U.el("div", { class: "form-row" }, U.el("label", {}, "type"), typeSel),
      paramBox);

    return {
      el,
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

  function _syntheticForm(kind) {
    const forms = {};
    const parts = [];
    if (kind === "wind") {
      forms.speed = _profileForm("Wind speed [m/s]",
        ["constant", "blocks", "harmonic", "linear"], { value: 10 });
      forms.direction = _profileForm("Wind direction [deg]",
        ["constant", "blocks", "harmonic", "rotational"], { value: 270 });
      parts.push(forms.speed.el, forms.direction.el);
    } else if (kind === "tide") {
      forms.level = _profileForm("Water level [m]",
        ["constant", "harmonic", "blocks", "linear"], { mean: 0, amplitude: 1, period: 12.42, value: 0 });
      parts.push(forms.level.el);
    } else {
      forms.hs = _profileForm("Wave height Hs [m]",
        ["constant", "blocks", "harmonic", "linear"], { value: 1 });
      forms.tp = _profileForm("Wave period Tp [s]",
        ["constant", "blocks", "harmonic", "linear"], { value: 6 });
      parts.push(forms.hs.el, forms.tp.el);
    }

    const dtInput = U.el("input", { type: "text", value: "1", style: "width:70px" });
    const genBtn = U.el("button", { class: "primary" }, "Generate & save");
    genBtn.addEventListener("click", async () => {
      const body = { kind, dt: Number(dtInput.value) * 3600 };
      for (const [name, form] of Object.entries(forms)) body[name] = form.spec();
      try {
        genBtn.disabled = true;
        const res = await Api.post("/api/conditions/synthetic", body);
        U.toast(`Wrote ${res.file} (${res.rows} rows)`, "ok");
        await _reloadConfig(kind);
        overview = await Api.get("/api/conditions");
        _plotSeries(kind, overview.kinds[kind]);
      } catch (err) {
        U.toast(err.message, "error");
      } finally {
        genBtn.disabled = false;
      }
    });

    const details = U.el("div", { class: "form-group" },
      U.el("span", { class: "fg-label" }, "Synthetic"),
      ...parts,
      U.el("div", { class: "form-row" }, U.el("label", {}, "output step [h]"), dtInput),
      U.el("div", { class: "btn-row" }, genBtn));
    return details;
  }

  /* ---- measured / reanalysis ---- */

  function _measuredForm(kind) {
    const sources = kind === "wind" ? ["waterinfo", "era5"] : ["waterinfo"];
    const sourceSel = U.el("select", {},
      ...sources.map((s) => U.el("option", { value: s },
        s === "era5" ? "ERA5 reanalysis (CDS)" : "waterinfo.rws.nl")));

    const stationBox = U.el("div", { class: "muted" }, "—");
    const findBtn = U.el("button", { class: "ghost" }, "Find stations near grid");
    findBtn.addEventListener("click", () => _findStations(sourceSel.value, kind, stationBox));

    const date0 = U.el("input", { type: "text", placeholder: "YYYY-MM-DD (default: sim start)" });
    const date1 = U.el("input", { type: "text", placeholder: "YYYY-MM-DD (default: sim end)" });

    const fetchBtn = U.el("button", { class: "primary" }, "Download & convert");
    const progress = U.el("div", { class: "muted", style: "font-size:11.5px" });
    fetchBtn.addEventListener("click", async () => {
      if (!selectedStation) { U.toast("Select a station/cell first", "error"); return; }
      const body = { source: sourceSel.value, kind, station: selectedStation };
      if (sourceSel.value === "waterinfo") body.station = selectedStation.id;
      if (date0.value.trim()) body.date0 = date0.value.trim();
      if (date1.value.trim()) body.date1 = date1.value.trim();
      try {
        fetchBtn.disabled = true;
        const res = await Api.post("/api/conditions/fetch", body);
        const out = await Api.waitJob(res.job, (j) => { progress.textContent = j.message || ""; });
        progress.textContent = "";
        U.toast(`Wrote ${out.file} (${out.rows} rows)`, "ok");
        await _reloadConfig(kind);
        overview = await Api.get("/api/conditions");
        _plotSeries(kind, overview.kinds[kind]);
      } catch (err) {
        U.toast(err.message, "error");
      } finally {
        fetchBtn.disabled = false;
      }
    });

    return U.el("div", { class: "form-group" },
      U.el("span", { class: "fg-label" }, "Measured / reanalysis"),
      U.el("div", { class: "form-row" }, U.el("label", {}, "source"), sourceSel),
      U.el("div", { class: "btn-row" }, findBtn),
      stationBox,
      U.el("div", { class: "form-row" }, U.el("label", {}, "from"), date0),
      U.el("div", { class: "form-row" }, U.el("label", {}, "to"), date1),
      U.el("div", { class: "btn-row" }, fetchBtn),
      progress);
  }

  async function _findStations(source, kind, box) {
    _clearStations();
    box.textContent = "searching…";
    try {
      const res = await Api.get(`/api/conditions/stations?source=${source}&kind=${kind}`);
      let payload = res;
      if (res.job) payload = await Api.waitJob(res.job);
      if (source === "era5" && res.configured === false) {
        box.textContent = "";
        box.append(U.el("div", {}, `⚠ ${res.reason}`));
      }
      const stations = payload.stations || [];
      if (!stations.length) {
        box.textContent = "no stations found nearby";
        return;
      }
      U.clear(box);
      for (const st of stations.slice(0, 15)) {
        const row = U.el("div", { class: "lp-row", style: "cursor:pointer" },
          U.el("span", { class: "lp-name" }, st.name || st.id),
          U.el("span", { class: "lp-mini" },
            st.dist_km !== undefined ? `${st.dist_km.toFixed(0)} km` : ""));
        row.addEventListener("click", () => {
          selectedStation = st;
          for (const r of box.querySelectorAll(".lp-row")) r.classList.remove("selected");
          row.classList.add("selected");
        });
        box.append(row);

        if (st.lon !== undefined || st.lat !== undefined) {
          const node = U.el("div", { class: "map-label" }, st.name || st.id);
          node.style.pointerEvents = "auto";
          node.style.cursor = "pointer";
          node.addEventListener("click", () => row.click());
          const marker = new maplibregl.Marker({ element: node, anchor: "bottom" })
            .setLngLat([st.lon, st.lat]).addTo(MapView.instance());
          stationMarkers.push(marker);
        }
      }
    } catch (err) {
      box.textContent = "";
      U.toast(err.message, "error");
    }
  }

  function _clearStations() {
    for (const marker of stationMarkers) marker.remove();
    stationMarkers = [];
    selectedStation = null;
  }

  async function _reloadConfig(kind) {
    const cfg = await Api.get("/api/config");
    App.state.config = cfg.values;
    App.emit("config-changed", overview.kinds[kind].config_key);
  }

  return { init };
})();
