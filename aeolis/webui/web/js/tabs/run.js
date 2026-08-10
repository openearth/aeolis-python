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
  let els = {};
  let activeBackend = "local";
  // HPC form state; getters read live values, base keeps the fields the
  // form does not expose (modules, conda_setup, extra_sbatch)
  const hpc = { fields: {}, els: {}, base: {}, manual: false, timer: null };
  const ACTIVE_STATES = new Set(["submitting", "queued", "running"]);

  function init() {
    Tabs.register("run", { enter: _enter, leave: _leave });
    App.on("project", () => { if (Tabs.current() === "run") _enter(); });
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

      els.checklist.append(_checkGroup("Must fix", by.error, "error", false));
      els.checklist.append(_checkGroup("Warnings", by.warn, "warn", false));
      els.checklist.append(_checkGroup("Passed", by.ok, "ok", true));
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

  const _winToLinux = (p) => {
    const m = /^([A-Za-z]):[\\/](.*)$/.exec(p || "");
    return m ? `/${m[1].toLowerCase()}/${m[2].replace(/\\/g, "/")}` : (p || "").replace(/\\/g, "/");
  };

  async function _hpcSettings(panel) {
    let cfg;
    try { cfg = await Api.get("/api/run/hpc"); }
    catch (err) { panel.append(U.el("div", { class: "muted" }, err.message)); return; }

    hpc.fields = {};
    hpc.base = cfg.profile;
    hpc.manual = false;
    const p = cfg.profile;
    const section = U.section("Deltares HYDRAX — job settings", {});

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

    // --- connection ---
    textField("Username", "user", { placeholder: "samAccountName, e.g. weste_bt" });
    hpc.els.password = U.el("input", { type: "password", placeholder: "AD password (not stored)" });
    rowNode("Password", hpc.els.password,
      "Your Deltares (AD) password — used only to open the SSH login. Kept in memory for this session, never saved.");
    textField("Login node", "host");

    // --- resources ---
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

    // --- environment & run location ---
    pathField("Conda env", "env_path",
      "Path to the conda environment to activate (on /p). Browse to pick it as a folder.",
      "/p/<project>/00_environments/<env>");

    const modeSel = U.el("select", {});
    for (const [v, l] of [["inplace", "Run in the project folder (already on /p)"],
                          ["copy", "Copy the project to a /p run folder"]]) {
      modeSel.append(U.el("option", { value: v, selected: v === (p.run_mode || "inplace") ? "" : null }, l));
    }
    hpc.fields.run_mode = () => modeSel.value;
    rowNode("Run directory", modeSel,
      "Run in place when the GUI already runs on the P-drive (WCF); copy when the GUI runs locally and the project must be moved to /p first.");

    const runDir = pathField("Run dir", "run_dir",
      "Folder on /p where the job runs (and where the .sh + outputs go).",
      "/p/<project>/03_simulations/<run>");
    const useProj = U.el("button", { class: "ghost" }, "Use current project folder");
    useProj.addEventListener("click", () => { runDir.value = cfg.project_dir_linux || ""; _scheduleHpcPreview(); });
    const modeNote = U.el("div", { class: "muted", style: "font-size:11.5px" });
    const syncMode = () => {
      modeNote.textContent = modeSel.value === "copy"
        ? "The project will be copied to the run dir on /p before submitting (GUI cache and old outputs are skipped)."
        : "The run dir should be the project's own folder on /p — nothing is copied.";
      if (modeSel.value === "inplace" && !runDir.value) runDir.value = cfg.project_dir_linux || "";
    };
    modeSel.addEventListener("change", () => { syncMode(); _scheduleHpcPreview(); });
    section.body.append(U.el("div", { class: "form-row" }, U.el("label", {}, ""), useProj), modeNote);
    syncMode();

    textField("Config file", "config");
    textField("E-mail (optional)", "mail_user", { placeholder: "you@deltares.nl" });

    // editable script preview
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

  async function _startHpc() {
    const password = hpc.els.password ? hpc.els.password.value : "";
    if (!password) { U.toast("Enter your HPC password", "error"); return; }
    try {
      els.log.textContent = "";
      logOffset = 0;
      const res = await Api.post("/api/run/hpc/start", {
        profile: _hpcProfile(),
        password,
        script: hpc.manual ? hpc.els.script.value : null,
      });
      U.toast(`Submitted to HYDRAX (job ${res.job_id})`, "ok");
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
