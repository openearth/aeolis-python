/* Run tab: execute the simulation with live progress and log.
 *
 * A runner-backend selector (local now, Deltares HPC planned) with a
 * pre-run checklist, play/stop, progress bar, elapsed/remaining time
 * and an incrementally-tailed log panel (1 s polling).
 */
"use strict";

const RunTab = (() => {

  let pollTimer = null;
  let logOffset = 0;
  let lastState = "idle";
  // where the backend reads output from, as last seen in the status
  // poll (undefined = not seen yet). When it changes server-side (e.g.
  // the automatic switch after an HPC submission) the viewer and the
  // availability strip are refreshed via the output-source event.
  let lastOutputDir;
  let els = {};
  let activeBackend = "local";
  // HPC form state; getters read live values, base keeps the fields the
  // form does not expose (modules, conda_setup, extra_sbatch). The
  // password and the connected/jobs state live here so they survive
  // panel rebuilds within the session (never persisted).
  const hpc = { fields: {}, els: {}, base: {}, manual: false, timer: null,
    password: "", connected: false, jobs: null, selectedJobId: null };
  const ACTIVE_STATES = new Set(["submitting", "queued", "running"]);

  function init() {
    Tabs.register("run", { enter: _enter, leave: _leave });
    App.on("project", () => {
      lastOutputDir = undefined;      // new project — source not seen yet
      if (Tabs.current() === "run") _enter();
    });
  }

  async function _enter() {
    if (!App.state.project) return;
    await _build();
    _startPolling();
  }

  function _leave() {
    // keep polling if a run is active (incl. HPC queued) so we notice completion
    if (!ACTIVE_STATES.has(lastState)) _stopPolling();
  }

  async function _build() {
    const panel = document.getElementById("run-panel");
    const restoreScroll = U.keepScroll(panel);
    U.clear(panel);
    els = {};

    // backend selector
    let backends = { backends: [], active: "local" };
    try { backends = await Api.get("/api/run/backends"); } catch { /* ignore */ }
    const backendSel = U.el("select", {});
    for (const b of backends.backends) {
      backendSel.append(U.el("option", {
        value: b.id, selected: b.id === backends.active ? "" : null,
        disabled: b.available ? null : "",
      }, b.title + (b.available ? "" : " — " + (b.note || "unavailable"))));
    }
    backendSel.addEventListener("change", async () => {
      try {
        await Api.post("/api/run/backend", { id: backendSel.value });
        await _build();          // swap in/out the HPC settings panel
      } catch (err) {
        U.toast(err.message, "error");
        backendSel.value = "local";
      }
    });
    activeBackend = backends.active;
    panel.append(U.el("div", { class: "form-row" },
      U.el("label", {}, "Run on"), backendSel));

    // HPC settings (only for the Deltares HYDRAX backend)
    if (activeBackend === "hpc") {
      await _hpcSettings(panel);
    } else {
      await _hpcReminder(panel, backends, backendSel);
    }

    // checklist
    els.checklist = U.el("div", { class: "form-group" });
    panel.append(els.checklist);
    await _refreshChecklist();

    // controls
    els.startBtn = U.el("button", { class: "primary", style: "min-width:90px" }, "▶ Start");
    els.stopBtn = U.el("button", { class: "danger", disabled: "" }, "■ Stop");
    els.startBtn.addEventListener("click", _start);
    els.stopBtn.addEventListener("click", _stop);
    panel.append(U.el("div", { class: "btn-row" }, els.startBtn, els.stopBtn));

    // progress
    els.progressBar = U.el("div");
    els.progressText = U.el("div", { class: "muted", style: "font-size:12px" }, "idle");
    panel.append(
      U.el("div", { class: "progress" }, els.progressBar),
      els.progressText,
    );

    // log
    els.log = U.el("pre", { id: "run-log" });
    panel.append(els.log);
    logOffset = 0;
    restoreScroll();   // content landed after awaits — restore again
  }

  const _CHECK_ICON = { ok: "✔", warn: "⚠", error: "✖" };

  function _checkPill(level, n) {
    return U.el("span", { class: `check-pill ${level}` }, `${_CHECK_ICON[level]} ${n}`);
  }

  function _checkGroup(title, list, level, collapsed) {
    if (!list.length) return null;
    const body = U.el("div", { class: "check-group-body" },
      ...list.map((c) => U.el("div", { class: `check-row ${level}` },
        U.el("span", { class: "check-ic" }, _CHECK_ICON[level]),
        U.el("span", { class: "check-txt", title: c.text }, c.text))));
    const head = U.el("header", {},
      U.el("span", { class: "caret" }, "▾"), title,
      U.el("span", { class: "count" }, String(list.length)));
    const wrap = U.el("div", { class: `check-group ${level} ${collapsed ? "collapsed" : ""}` },
      head, body);
    head.addEventListener("click", () => wrap.classList.toggle("collapsed"));
    return wrap;
  }

  async function _refreshChecklist() {
    try {
      const res = await Api.get("/api/run/checklist");
      U.clear(els.checklist);
      const checks = res.checks || [];
      const by = { error: [], warn: [], ok: [] };
      for (const c of checks) (by[c.level] || by.ok).push(c);

      els.checklist.append(U.el("div", { class: "checks-summary" },
        U.el("span", { class: "fg-label" }, "Pre-run checks"),
        U.el("span", { class: "grow" }),
        by.error.length ? _checkPill("error", by.error.length) : null,
        by.warn.length ? _checkPill("warn", by.warn.length) : null,
        _checkPill("ok", by.ok.length)));

      // _checkGroup returns null for empty groups; append(null) would
      // render the literal text "null"
      for (const grp of [_checkGroup("Must fix", by.error, "error", false),
                         _checkGroup("Warnings", by.warn, "warn", false),
                         _checkGroup("Passed", by.ok, "ok", true)]) {
        if (grp) els.checklist.append(grp);
      }
      els.ready = res.ready;
    } catch (err) {
      U.clear(els.checklist);
      els.checklist.append(U.el("div", { class: "muted" }, err.message));
    }
  }

  async function _start() {
    await _refreshChecklist();
    if (els.ready === false &&
        !window.confirm("Some required inputs look missing. Start anyway?")) return;
    if (activeBackend === "hpc") { await _startHpc(); return; }
    try {
      els.log.textContent = "";
      logOffset = 0;
      await Api.post("/api/run/start", {});
      U.toast("Simulation started", "ok");
      _startPolling();
    } catch (err) {
      U.toast(err.message, "error");
    }
  }

  /* ---------------- Deltares HYDRAX (HPC) ---------------- */

  /* An HPC job keeps running when the GUI closes. When this project has
   * recorded submissions and the backend is (back to) local — i.e. the
   * GUI was reopened — point the user at the reattach flow. */
  async function _hpcReminder(panel, backends, backendSel) {
    const hpcBackend = (backends.backends || []).find((b) => b.id === "hpc");
    if (!hpcBackend || !hpcBackend.available) return;
    let cfg;
    try { cfg = await Api.get("/api/run/hpc"); } catch { return; }
    const last = (cfg.jobs || [])[0];
    if (!last) return;
    const when = last.submitted
      ? ` on ${new Date(last.submitted * 1000).toLocaleString()}` : "";
    const openBtn = U.el("button", { class: "ghost", style: "font-size:12px" },
      "Open HPC settings");
    openBtn.addEventListener("click", async () => {
      try {
        await Api.post("/api/run/backend", { id: "hpc" });
        await _build();
      } catch (err) { U.toast(err.message, "error"); }
    });
    panel.append(U.el("div", { class: "hpc-remind" },
      U.el("span", {},
        `Job ${last.job_id} was submitted to HYDRAX from this project${when} `
        + "— it may still be running. Reconnect via the HPC settings."),
      openBtn));
  }

  const _winToLinux = (p) => {
    const m = /^([A-Za-z]):[\\/](.*)$/.exec(p || "");
    return m ? `/${m[1].toLowerCase()}/${m[2].replace(/\\/g, "/")}` : (p || "").replace(/\\/g, "/");
  };

  /* Shared field builders, scoped to one section box (the login box and
   * the job-settings box both use them). */
  function _formHelpers(section, p) {
    const helpIcon = (text) => U.el("span", { class: "help-i", title: text }, U.icon("info", 14));
    const rowNode = (label, node, help) => section.body.append(
      U.el("div", { class: "form-row" },
        U.el("label", {}, label, help ? helpIcon(help) : null), node));
    const textField = (label, key, opts = {}) => {
      const input = U.el("input", {
        type: opts.type || "text", value: p[key] != null ? p[key] : "",
        placeholder: opts.placeholder || "",
      });
      hpc.fields[key] = () => input.value;
      input.addEventListener("input", _scheduleHpcPreview);
      rowNode(label, input, opts.help);
      return input;
    };
    // a path field with a Browse button that converts a Windows/mounted
    // path to its /p Linux form
    const pathField = (label, key, help, placeholder) => {
      const input = U.el("input", { type: "text", value: p[key] != null ? p[key] : "", placeholder, style: "flex:1" });
      hpc.fields[key] = () => input.value;
      input.addEventListener("input", _scheduleHpcPreview);
      const browse = U.el("button", { class: "ghost", title: "Browse — converts the picked folder to a /p path" }, "Browse…");
      browse.addEventListener("click", async () => {
        const win = await Api.pickFolder({ title: `Select ${label}` }).catch(() => null);
        if (!win) return;
        input.value = _winToLinux(win);
        _scheduleHpcPreview();
      });
      rowNode(label, U.el("div", { style: "display:flex;gap:6px;flex:1" }, input, browse), help);
      return input;
    };

    const subhead = (text) => section.body.append(
      U.el("div", { class: "hpc-subhead" }, text));

    return { rowNode, textField, pathField, subhead };
  }

  /* The HPC panel, top to bottom: (1) login box, (2) once connected the
   * simulations on the cluster — selectable for monitoring, with the
   * output-source decision — and (3) the job settings for a new run. */
  async function _hpcSettings(panel) {
    let cfg;
    try { cfg = await Api.get("/api/run/hpc"); }
    catch (err) { panel.append(U.el("div", { class: "muted" }, err.message)); return; }
    hpc.fields = {};
    hpc.base = cfg.profile;
    hpc.manual = false;
    _hpcLoginBox(panel, cfg);
    if (hpc.connected) _hpcJobsBox(panel, cfg);
    _hpcJobSettings(panel, cfg);
  }

  /* --- login: connect first, then your cluster jobs appear below --- */

  function _hpcLoginBox(panel, cfg) {
    const section = U.section("Deltares HYDRAX — login", {});
    const H = _formHelpers(section, cfg.profile);
    H.textField("Username", "user", { placeholder: "samAccountName, e.g. weste_bt" });
    hpc.els.password = U.el("input", {
      type: "password", placeholder: "AD password (not stored)",
      value: hpc.password || "",
    });
    hpc.els.password.addEventListener("input", () => {
      hpc.password = hpc.els.password.value;
    });
    H.rowNode("Password", hpc.els.password,
      "Your Deltares (AD) password — used only to open the SSH login. Kept in memory for this session, never saved.");
    H.textField("Login node", "host");
    const connectBtn = U.el("button", { class: "primary", style: "min-width:110px" },
      hpc.connected ? "⟳ Reconnect" : "Connect");
    const status = U.el("span", { class: "muted", style: "font-size:12px" },
      hpc.connected
        ? `Connected — ${(hpc.jobs || []).length} job(s) listed below.`
        : "Connect to list your simulations on the cluster.");
    connectBtn.addEventListener("click", () => _hpcConnect(connectBtn, status));
    hpc.els.password.addEventListener("keydown", (ev) => {
      if (ev.key === "Enter") connectBtn.click();
    });
    H.rowNode("", U.el("div",
      { style: "display:flex;gap:10px;align-items:center;flex-wrap:wrap" },
      connectBtn, status));
    panel.append(section.wrap);
  }

  async function _hpcConnect(btn, status) {
    if (!hpc.password) { U.toast("Enter your password first", "error"); return; }
    btn.disabled = true;
    if (status) status.textContent = "Connecting (SSH) …";
    try {
      // persist the typed settings first, so the rebuild keeps them
      await Api.post("/api/run/hpc", { profile: _hpcProfile() }).catch(() => {});
      const res = await Api.post("/api/run/hpc/jobs",
        { profile: _hpcProfile(), password: hpc.password });
      hpc.connected = true;
      hpc.jobs = res.jobs || [];
      await _build();
    } catch (err) {
      btn.disabled = false;
      if (status) status.textContent = `✖ ${err.message}`;
      U.toast(err.message, "error");
    }
  }

  /* --- the user's simulations on the cluster (live squeue merged with
   * this project's recorded submissions). Selecting one offers the
   * monitor plus, for a run living outside the open project, the
   * decision to read THAT run folder's output in the Viewer instead of
   * the project's own output file. --- */

  function _hpcJobsBox(panel, cfg) {
    const jobs = hpc.jobs || [];
    const section = U.section("Simulations on the cluster",
      { count: jobs.length || null });
    const body = section.body;
    if (!jobs.length) {
      body.append(U.el("div", { class: "muted" },
        "No jobs are queued or running for this account, and no earlier "
        + "submissions are recorded for this project."));
      panel.append(section.wrap);
      return;
    }
    let selected = jobs.find((j) => String(j.job_id) === hpc.selectedJobId)
      || jobs.find((j) => _JOB_LIVE.has(j.state)) || jobs[0];
    hpc.selectedJobId = String(selected.job_id);
    const actions = U.el("div", { style: "margin-top:8px" });
    const norm = (s) => String(s || "").replace(/\/+$/, "").toLowerCase();
    const renderActions = () => {
      U.clear(actions);
      const j = selected;
      const workdir = (j.workdir || "").trim();
      const external = workdir && norm(workdir) !== norm(cfg.project_dir_linux);
      let followChk = null;
      if (external) {
        followChk = U.el("input", { type: "checkbox", checked: "" });
        actions.append(U.el("label", { class: "choice-row", style: "font-size:12px" },
          followChk, U.el("span", {},
            `Show this run's output in the Viewer — read it from ${workdir} `
            + "(the mounted P-drive) instead of the project's own output file")));
      } else if (workdir) {
        actions.append(U.el("div", { class: "muted", style: "font-size:11.5px" },
          "This job runs in the open project folder — the Viewer reads its output already."));
      }
      const goBtn = U.el("button", { class: "primary" }, "Monitor this job");
      goBtn.addEventListener("click", () =>
        _hpcMonitor(j, Boolean(followChk && followChk.checked)));
      actions.append(U.el("div", { class: "btn-row", style: "margin-top:8px" }, goBtn));
    };
    for (const j of jobs) {
      const radio = U.el("input", { type: "radio", name: "hpc-job",
        checked: j === selected ? "" : null });
      radio.addEventListener("change", () => {
        selected = j;
        hpc.selectedJobId = String(j.job_id);
        renderActions();
      });
      const bits = [`job ${j.job_id}`, j.name, j.state || "state unknown",
        j.time, j.partition].filter(Boolean).join(" · ");
      body.append(U.el("label", { class: "choice-row hpc-job-row" }, radio,
        U.el("span", { style: "min-width:0" },
          U.el("div", {}, bits + (j.recorded ? "  — submitted from this project" : "")),
          U.el("div", { class: "muted", style: "font-size:11px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap" },
            (j.workdir || "") + (j.submitted
              ? ` · submitted ${new Date(j.submitted * 1000).toLocaleString()}` : "")))));
    }
    renderActions();
    body.append(actions);
    panel.append(section.wrap);
  }

  async function _hpcMonitor(job, followOutput) {
    if (!hpc.password) { U.toast("Enter your HPC password first", "error"); return; }
    // 1) the output-source decision — purely local (the mounted P:), so
    //    it is applied first and works even when the SSH monitor fails
    try {
      const res = await Api.post("/api/output/source", followOutput && job.workdir
        ? { dir: job.workdir, config: job.config || null }
        : { dir: null });
      lastOutputDir = res.dir || null;   // keep _poll from double-toasting
      if (followOutput && job.workdir) {
        U.toast(res.exists
          ? "Viewer output: following the run folder"
          : "Viewer will follow the run folder (no output file written yet)", "ok");
      }
      App.emit("output-source");
    } catch (err) {
      U.toast(`Output source: ${err.message}`, "error");
    }
    // 2) (re)attach the monitor — a no-op when this job is already
    //    monitored, a clean switch when another one was
    try {
      await Api.post("/api/run/hpc/attach",
        { profile: _hpcProfile(), password: hpc.password, job_id: job.job_id });
    } catch (err) { U.toast(err.message, "error"); return; }
    els.log.textContent = "";
    logOffset = 0;
    U.toast(_JOB_LIVE.has(job.state)
      ? `Monitoring job ${job.job_id}`
      : `Job ${job.job_id} already ended — fetching its final log`, "ok");
    _startPolling();
  }

  /* --- job settings for submitting a new run --- */

  function _hpcJobSettings(panel, cfg) {
    const p = cfg.profile;
    const section = U.section("Deltares HYDRAX — job settings", {});
    const { rowNode, textField, pathField, subhead } = _formHelpers(section, p);

    // --- resources ---
    subhead("Job resources");
    textField("Job name", "job_name");
    const partSel = U.el("select", {});
    for (const part of cfg.partitions) {
      partSel.append(U.el("option", { value: part, selected: part === p.partition ? "" : null }, part));
    }
    partSel.addEventListener("change", _scheduleHpcPreview);
    hpc.fields.partition = () => partSel.value;
    rowNode("Partition", partSel,
      "HYDRAX node type: 1vcpu = 1 core/8 GB … 60vcpu = 60 cores/480 GB; 'test' is a 30-min quick check.");
    textField("Tasks (--ntasks)", "ntasks", { type: "number",
      help: "Number of parallel tasks (MPI processes). AeoLiS normally uses 1." });
    textField("CPUs / task", "cpus_per_task", { type: "number",
      help: "CPU cores (threads) per task. Usually 1 for AeoLiS." });
    textField("Walltime limit", "walltime", {
      help: "Maximum run time, format days-hours:minutes:seconds (e.g. 5-00:00:00). The job is KILLED after this, so set it generously. Max 32 days on 1/4vcpu, 24 days on larger partitions." });

    textField("E-mail (optional)", "mail_user", { placeholder: "you@deltares.nl",
      help: "Slurm mails job begin/end/fail notifications to this address." });

    /* --- run location: one clear choice ---
     * Either the job runs directly in the project folder (only possible
     * when that folder lives on the shared P-drive the cluster mounts),
     * or the project is copied to a chosen /p folder first. */
    subhead("Run location");
    const projLinux = cfg.project_dir_linux || "";
    const onP = /^\/p\//i.test(projLinux);
    const wantCopy = (p.run_mode === "copy") || !onP;
    const inplaceRadio = U.el("input", { type: "radio", name: "hpc-runloc",
      checked: wantCopy ? null : "", disabled: onP ? null : "" });
    const copyRadio = U.el("input", { type: "radio", name: "hpc-runloc",
      checked: wantCopy ? "" : null });
    const copyDir = U.el("input", { type: "text",
      value: (p.run_dir && p.run_dir !== projLinux) ? p.run_dir : "",
      placeholder: "/p/<project>/03_simulations/<run>", style: "flex:1" });
    copyDir.addEventListener("input", _scheduleHpcPreview);
    const copyBrowse = U.el("button", { class: "ghost",
      title: "Browse — converts the picked folder to its /p path" }, "Browse…");
    copyBrowse.addEventListener("click", async () => {
      const win = await Api.pickFolder({ title: "Select run folder on the P-drive" }).catch(() => null);
      if (!win) return;
      copyDir.value = _winToLinux(win);
      _scheduleHpcPreview();
    });
    hpc.fields.run_mode = () => (copyRadio.checked ? "copy" : "inplace");
    hpc.fields.run_dir = () => (copyRadio.checked ? copyDir.value : projLinux);
    const copyRow = U.el("div", { style: "display:flex;gap:6px;margin:2px 0 0 24px" },
      copyDir, copyBrowse);
    const locNote = U.el("div", { class: "muted", style: "font-size:11.5px;margin:4px 0 2px" });
    const syncLoc = () => {
      copyRow.style.display = copyRadio.checked ? "" : "none";
      locNote.textContent = copyRadio.checked
        ? "The project inputs are copied to this /p folder before submitting; the job runs and writes its output there (GUI cache and old outputs are skipped)."
        : `The job runs directly in the project folder (${projLinux}) — nothing is copied.`;
    };
    for (const r of [inplaceRadio, copyRadio]) {
      r.addEventListener("change", () => { syncLoc(); _scheduleHpcPreview(); });
    }
    section.body.append(
      U.el("label", { class: "choice-row", style: "font-size:12.5px" }, inplaceRadio,
        U.el("span", {}, onP
          ? "Run in the project folder — it is on the P-drive"
          : "Run in the project folder — unavailable: this project is not on the P-drive")),
      U.el("label", { class: "choice-row", style: "font-size:12.5px" }, copyRadio,
        U.el("span", {}, "Copy the project to a run folder on /p")),
      copyRow, locNote);
    syncLoc();
    textField("Config file", "config", {
      help: "The aeolis configuration file the job runs, relative to the run folder." });

    // --- environment ---
    subhead("Environment");
    const kindSel = U.el("select", {});
    kindSel.append(
      U.el("option", { value: "conda", selected: p.env_kind === "venv" ? null : "" },
        "conda / miniforge"),
      U.el("option", { value: "venv", selected: p.env_kind === "venv" ? "" : null },
        "Python venv (uv)"));
    kindSel.addEventListener("change", _scheduleHpcPreview);
    hpc.fields.env_kind = () => kindSel.value;
    // a venv needs no module/conda bootstrap lines
    hpc.fields.modules = () => (kindSel.value === "venv" ? [] : (hpc.base.modules || []));
    rowNode("Environment kind", kindSel,
      "conda: module load + 'conda activate <path>'. venv: 'source <path>/bin/activate' — pick this for the uv environments in 00_environments.");
    pathField("Environment", "env_path",
      "Path to the environment on /p (conda prefix or venv root). Browse to pick it as a folder.",
      "/p/<project>/00_environments/<env>");
    // the model does not create missing output directories, so the job
    // script must (e.g. output_file = output/aeolis.nc -> mkdir -p output)
    hpc.fields.pre_lines = () => (cfg.output_subdir
      ? [`mkdir -p ${cfg.output_subdir}`] : (hpc.base.pre_lines || []));

    // editable script preview
    subhead("Job script");
    hpc.els.script = U.el("textarea", { class: "hpc-script", spellcheck: "false", readonly: "" });
    hpc.els.script.value = cfg.script;
    hpc.els.script.addEventListener("input", () => { hpc.manual = true; });
    const manualToggle = U.el("input", { type: "checkbox" });
    manualToggle.addEventListener("change", () => {
      hpc.manual = manualToggle.checked;
      if (hpc.manual) hpc.els.script.removeAttribute("readonly");
      else { hpc.els.script.setAttribute("readonly", ""); _scheduleHpcPreview(true); }
    });
    const saveScriptBtn = U.el("button", { class: "ghost", style: "font-size:12px" }, "Save .sh…");
    saveScriptBtn.addEventListener("click", _saveHpcScript);
    section.body.append(
      U.el("div", { class: "form-row", style: "margin-top:6px" },
        U.el("label", {}, "Job script"),
        U.el("div", { style: "display:flex;gap:10px;align-items:center;flex:1" },
          U.el("label", { class: "choice-row", style: "font-size:12px" },
            manualToggle, U.el("span", {}, "edit manually")),
          U.el("span", { class: "grow" }),
          saveScriptBtn)),
      hpc.els.script,
      U.el("div", { class: "muted", style: "font-size:11.5px" },
        cfg.available
          ? "Only the job script is written over SSH; the model reads its inputs from the run dir on /p."
          : `⚠ ${cfg.note}`));

    panel.append(section.wrap);
  }

  function _hpcProfile() {
    const p = { ...hpc.base };
    for (const [key, get] of Object.entries(hpc.fields)) p[key] = get();
    p.ntasks = Number(p.ntasks) || 1;
    p.cpus_per_task = Number(p.cpus_per_task) || 1;
    return p;
  }

  // regenerate + persist the script preview from the form (debounced),
  // unless the user has taken manual control of the textarea
  function _scheduleHpcPreview(force) {
    if (hpc.manual && force !== true) return;
    clearTimeout(hpc.timer);
    hpc.timer = setTimeout(async () => {
      try {
        const res = await Api.post("/api/run/hpc", { profile: _hpcProfile() });
        if (!hpc.manual && hpc.els.script) hpc.els.script.value = res.script;
      } catch { /* ignore preview errors */ }
    }, 350);
  }

  async function _saveHpcScript() {
    const script = hpc.manual ? hpc.els.script.value
      : (await Api.post("/api/run/hpc", { profile: _hpcProfile() })).script;
    const jobName = (hpc.fields.job_name && hpc.fields.job_name()) || "aeolis";
    const path = await Api.pickFile({ save: true, title: "Save job script",
      filename: `${jobName}.sh`, patterns: [["Shell script", "*.sh"], ["All files", "*.*"]] }).catch(() => null);
    if (!path) return;
    try {
      await Api.post("/api/run/hpc/save_script", { path, script });
      U.toast("Script saved", "ok");
    } catch (err) { U.toast(err.message, "error"); }
  }

  const _JOB_LIVE = new Set(["PENDING", "RUNNING", "COMPLETING", "CONFIGURING"]);

  async function _startHpc() {
    const password = hpc.password || "";
    if (!password) { U.toast("Enter your HPC password", "error"); return; }
    try {
      els.log.textContent = "";
      logOffset = 0;
      await Api.post("/api/run/hpc/start", {
        profile: _hpcProfile(),
        password,
        script: hpc.manual ? hpc.els.script.value : null,
      });
      U.toast("Submitting to HYDRAX — follow the log below", "ok");
      _startPolling();
    } catch (err) {
      U.toast(err.message, "error");
    }
  }

  async function _stop() {
    try {
      await Api.post("/api/run/stop", {});
      U.toast("Stop requested", "");
    } catch (err) {
      U.toast(err.message, "error");
    }
  }

  function _startPolling() {
    if (pollTimer) return;
    pollTimer = setInterval(_poll, 1000);
    _poll();
  }

  function _stopPolling() {
    clearInterval(pollTimer);
    pollTimer = null;
  }

  async function _poll() {
    let status;
    try {
      status = await Api.get("/api/run/status");
    } catch {
      return;
    }

    const running = status.running;
    if (els.startBtn) {
      els.startBtn.disabled = running;
      els.stopBtn.disabled = !running;
    }

    // backend-side output-source switches (e.g. automatically after an
    // HPC submission): refresh the viewer + availability strip once
    const o = status.output;
    if (o) {
      const dir = o.dir || null;
      if (lastOutputDir !== undefined && dir !== lastOutputDir) {
        U.toast(o.override
          ? "Viewer output now follows the HPC run folder on P:"
          : "Viewer output: back to the project's own file", "ok");
        App.emit("output-source");
      }
      lastOutputDir = dir;
    }

    const p = status.progress || {};
    if (els.progressBar) {
      els.progressBar.style.width = `${p.percent || 0}%`;
      const parts = [];
      // HPC queue/job monitor line
      const h = status.hpc;
      if (h && h.job_id) {
        parts.push(`job ${h.job_id}`);
        if (h.state) parts.push(h.reason && h.reason !== "None" ? `${h.state} (${h.reason})` : h.state);
        if (h.partition) parts.push(h.partition);
        if (h.time) parts.push(h.time);
      }
      if (running || status.state === "finished") {
        parts.push(`${(p.percent || 0).toFixed(1)}%`);
        if (p.elapsed) parts.push(`elapsed ${p.elapsed}`);
        if (p.remaining && running) parts.push(`remaining ${p.remaining}`);
        if (p.avg_dt) parts.push(`avg dt ${p.avg_dt}s`);
      }
      els.progressText.textContent = parts.length
        ? `${status.state} — ${parts.join(" · ")}` : status.state;
    }

    // log tail
    if (els.log) {
      try {
        const tail = await Api.get(`/api/run/log?offset=${logOffset}`);
        if (tail.lines.length) {
          const atBottom = els.log.scrollTop + els.log.clientHeight >= els.log.scrollHeight - 30;
          els.log.textContent += tail.lines.join("\n") + "\n";
          logOffset = tail.offset;
          const maxChars = 400000;
          if (els.log.textContent.length > maxChars) {
            els.log.textContent = els.log.textContent.slice(-maxChars);
          }
          if (atBottom) els.log.scrollTop = els.log.scrollHeight;
        } else {
          logOffset = tail.offset;
        }
      } catch { /* ignore */ }
    }

    // completion transitions (local: running→finished; HPC also queued→…)
    if (ACTIVE_STATES.has(lastState) && !ACTIVE_STATES.has(status.state)) {
      if (status.state === "finished") {
        U.toast("Simulation finished", "ok");
        App.emit("run-finished");
      } else if (status.state === "error") {
        const ex = status.exit_code != null ? ` (exit ${status.exit_code})` : "";
        U.toast(`Simulation failed${ex}`, "error");
      } else if (status.state === "stopped") {
        U.toast("Simulation stopped", "");
      }
      if (Tabs.current() !== "run") _stopPolling();
    }
    lastState = status.state;
  }

  return { init };
})();
