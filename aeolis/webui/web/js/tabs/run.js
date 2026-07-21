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
    // keep polling if a run is active so we notice completion
    if (lastState !== "running") _stopPolling();
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
      } catch (err) {
        U.toast(err.message, "error");
        backendSel.value = "local";
      }
    });
    panel.append(U.el("div", { class: "form-row" },
      U.el("label", {}, "Run on"), backendSel));

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

  async function _refreshChecklist() {
    try {
      const res = await Api.get("/api/run/checklist");
      U.clear(els.checklist);
      els.checklist.append(U.el("span", { class: "fg-label" }, "Pre-run checks"));
      for (const check of res.checks) {
        const level = check.level || (check.ok ? "ok" : "error");
        const icon = level === "ok" ? "✔" : level === "warn" ? "⚠" : "✖";
        const style = level === "error" ? "color:var(--danger)"
          : level === "warn" ? "color:#a06a00" : "";
        els.checklist.append(U.el("div", { class: "lp-row" },
          U.el("span", { class: "eye", style }, icon),
          U.el("span", { class: "lp-name", style, title: check.text }, check.text)));
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
      if (running || status.state === "finished") {
        const parts = [];
        parts.push(`${(p.percent || 0).toFixed(1)}%`);
        if (p.elapsed) parts.push(`elapsed ${p.elapsed}`);
        if (p.remaining && running) parts.push(`remaining ${p.remaining}`);
        if (p.avg_dt) parts.push(`avg dt ${p.avg_dt}s`);
        els.progressText.textContent = `${status.state} — ${parts.join(" · ")}`;
      } else {
        els.progressText.textContent = status.state;
      }
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

    // completion transitions
    if (lastState === "running" && status.state !== "running") {
      if (status.state === "finished") {
        U.toast("Simulation finished", "ok");
        App.emit("run-finished");
      } else if (status.state === "error") {
        U.toast(`Simulation failed (exit ${status.exit_code})`, "error");
      } else if (status.state === "stopped") {
        U.toast("Simulation stopped", "");
      }
      if (Tabs.current() !== "run") _stopPolling();
    }
    lastState = status.state;
  }

  return { init };
})();
