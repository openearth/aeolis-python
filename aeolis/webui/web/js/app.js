/* Application bootstrap: wires everything together and handles the
 * project open/create flow. */
"use strict";

(async function main() {

  Theme.init();
  Tabs.init();
  Playbar.init();
  Graphs.init();
  CondHud.init();

  await MapView.init("map");
  await Styles.load();   // colormap presets (global, used by the Viewer)

  Draw.init();

  SettingsTab.init();
  GridTab.init();
  DomainTab.init();
  ConditionsTab.init();
  RunTab.init();
  ViewerTab.init();
  Transect.init();

  _wireChrome();

  /* ---- project bootstrap ---- */
  try {
    const info = await Api.get("/api/project");
    if (info.open) {
      await _projectOpened(info);
    } else {
      _showProjectModal(info.recent || []);
    }
  } catch (err) {
    U.toast(`Backend error: ${err.message}`, "error");
  }

  /* =================================================================
   * chrome: panel resizing, collapsibles, project modal
   * ================================================================= */

  function _wireChrome() {
    U.initTruncationTips();

    // sidebar resize (the playbar starts right of the sidebar, so keep
    // its inset in sync via a CSS variable)
    const sidebar = document.getElementById("sidebar");
    const syncSidebarWidth = (width) => {
      sidebar.style.width = `${width}px`;
      document.documentElement.style.setProperty("--sidebar-live-w", `${width}px`);
    };
    _dragResize(document.getElementById("sidebar-grip"), (ev) => {
      const width = U.clamp(ev.clientX, 260, 640);
      syncSidebarWidth(width);
      App.state.ui.sidebarWidth = width;
    }, () => App.touchUi());
    syncSidebarWidth(sidebar.getBoundingClientRect().width || 340);

    // graphs resize
    const graphs = document.getElementById("graphs-wrap");
    _dragResize(document.getElementById("split-grip"), (ev) => {
      const total = document.getElementById("viewport").getBoundingClientRect();
      const height = U.clamp(total.bottom - ev.clientY, 60, total.height - 120);
      graphs.style.height = `${height}px`;
      App.state.ui.graphsHeight = height;
    }, () => {
      App.touchUi();
      Graphs.resizeAll();
      const map = MapView.instance();
      if (map) map.resize();
    });

    // project button
    document.getElementById("btn-project").addEventListener("click", async () => {
      const info = await Api.get("/api/project");
      _showProjectModal(info.recent || []);
    });

    // documentation (one entry point; pages are browsable inside the popup)
    document.getElementById("btn-docs").addEventListener("click", () => {
      DocsPopup.open("https://aeolis.readthedocs.io/en/update_documentation/", "AeoLiS documentation");
    });

    // windrose — the single entry point (source + optional period chooser)
    document.getElementById("btn-windrose").addEventListener("click", async () => {
      if (!App.state.project) { U.toast("Open a project first", "error"); return; }
      let sources = [];
      try {
        sources = await ConditionsTab.windroseSources();
      } catch (err) { U.toast(err.message, "error"); return; }
      Windrose.openChooser({ sources });
    });

    // config file viewer + reveal-in-explorer
    document.getElementById("btn-cfg-view").addEventListener("click", async () => {
      try {
        const res = await Api.get("/api/config/raw");
        const popup = Popup.open({ title: res.path, width: 720 });
        const pre = U.el("pre", { class: "cfg-view" }, res.text);
        const revealBtn = U.el("button", { class: "ghost" },
          U.icon("open", 14), " Show in Explorer");
        revealBtn.addEventListener("click", () => Api.post("/api/project/reveal"));
        const copyBtn = U.el("button", { class: "ghost" }, U.icon("copy", 14), " Copy");
        copyBtn.addEventListener("click", () => {
          navigator.clipboard.writeText(res.text)
            .then(() => U.toast("Copied to clipboard", "ok"))
            .catch(() => U.toast("Copy failed", "error"));
        });
        popup.body.append(pre,
          U.el("div", { class: "btn-row", style: "justify-content:flex-end" },
            copyBtn, revealBtn));
      } catch (err) {
        U.toast(err.message, "error");
      }
    });
    document.getElementById("btn-cfg-reveal").addEventListener("click", () => {
      Api.post("/api/project/reveal").catch((err) => U.toast(err.message, "error"));
    });
    document.getElementById("btn-duplicate").addEventListener("click", _duplicateModel);
    document.getElementById("btn-backup").addEventListener("click", _backupModel);
    App.on("project", () => {
      document.getElementById("btn-cfg-view").disabled = false;
      document.getElementById("btn-cfg-reveal").disabled = false;
      document.getElementById("btn-duplicate").disabled = false;
      document.getElementById("btn-backup").disabled = false;
    });
  }

  /* Clone the whole model (config + input files + GUI state) into a new
   * folder the user picks, then open it — a safe starting point for a
   * variant that never touches the original files. */
  async function _duplicateModel() {
    if (!App.state.project) return;
    const parent = await Api.pickFolder({ title: "Choose where to create the model copy" })
      .catch((err) => { U.toast(err.message, "error"); return null; });
    if (!parent) return;

    const popup = Popup.open({ title: "Duplicate model into a new folder", width: 420 });
    const inp = U.el("input", { type: "text", value: `${App.state.project.name || "model"}_copy` });
    inp.addEventListener("keydown", (ev) => { if (ev.key === "Enter") go(); });
    const outCb = U.el("input", { type: "checkbox", id: "dup-outputs", checked: "" });
    const gatherRadio = U.el("input", { type: "radio", name: "dup-inputs", id: "dup-gather", value: "gather", checked: "" });
    const keepRadio = U.el("input", { type: "radio", name: "dup-inputs", id: "dup-keep", value: "keep" });
    const btn = U.el("button", { class: "primary" }, "Create copy");
    btn.addEventListener("click", go);
    popup.body.append(
      U.el("div", { class: "muted", style: "font-size:12px;margin-bottom:6px" },
        `A new folder is created inside ${parent}. All model files (config, grids, `
        + "timeseries, GUI state) are copied; the gui cache is never copied."),
      U.el("div", { class: "form-row" }, U.el("label", {}, "New folder name"), inp),
      U.el("div", { class: "muted", style: "font-size:11.5px;margin:8px 0 3px" }, "Input files referenced from outside this folder:"),
      U.el("div", { class: "choice-row" }, gatherRadio,
        U.el("label", { for: "dup-gather" }, "Gather all into the new folder (self-contained)")),
      U.el("div", { class: "choice-row" }, keepRadio,
        U.el("label", { for: "dup-keep" }, "Keep originals in place (absolute links)")),
      U.el("div", { class: "choice-row", style: "margin-top:6px" }, outCb,
        U.el("label", { for: "dup-outputs" }, "Include run outputs (aeolis.nc, logs)")),
      U.el("div", { class: "btn-row", style: "justify-content:flex-end" }, btn));

    async function go() {
      const name = inp.value.trim();
      if (!name) { U.toast("Enter a folder name", "error"); return; }
      btn.disabled = true;
      let info;
      try {
        info = await Api.post("/api/project/duplicate",
          { parent, name, include_outputs: outCb.checked,
            input_mode: keepRadio.checked ? "keep" : "gather" });
      } catch (err) {
        btn.disabled = false;
        U.toast(err.message, "error");
        return;
      }
      popup.close();
      try {
        await _projectOpened(info);
        U.toast(`Copied to ${info.name} — now editing the copy`, "ok");
        if (info.skipped && info.skipped.length) {
          U.toast(`Could not find ${info.skipped.length} referenced file(s): ${info.skipped.join(", ")}`, "error");
        }
      } catch (err) {
        U.toast(`Copied, but loading state failed: ${err.message}`, "error");
      }
    }
  }

  /* Snapshot the entire model setup (config + inputs + GUI state; raw
   * data optional, run outputs optional) into a timestamped zip in
   * <project>/backups/. */
  async function _backupModel() {
    if (!App.state.project) return;
    const popup = Popup.open({ title: "Backup model setup", width: 440 });
    const rawCb = U.el("input", { type: "checkbox", id: "bk-raw", checked: "" });
    const outCb = U.el("input", { type: "checkbox", id: "bk-out" });
    const progress = U.progressBar();
    const btn = U.el("button", { class: "primary" }, "Create backup");
    btn.addEventListener("click", go);
    popup.body.append(
      U.el("div", { class: "muted", style: "font-size:12px;margin-bottom:6px" },
        "Saves a timestamped zip in this project's backups folder with the "
        + "configuration, all model input files and the GUI state. Input files "
        + "referenced from outside the project folder are gathered into the zip "
        + "as well, so the backup is self-contained."),
      U.el("div", { class: "choice-row" }, rawCb,
        U.el("label", { for: "bk-raw" }, "Include raw data (downloads / imports)")),
      U.el("div", { class: "choice-row" }, outCb,
        U.el("label", { for: "bk-out" }, "Include run output file(s) (aeolis.nc, logs)")),
      progress.el,
      U.el("div", { class: "btn-row", style: "justify-content:flex-end" }, btn));

    async function go() {
      btn.disabled = true;
      progress.start("archiving…");
      try {
        const res = await Api.post("/api/project/backup",
          { include_rawdata: rawCb.checked, include_outputs: outCb.checked });
        const out = await Api.waitJob(res.job, (j) => progress.update(j));
        popup.close();
        U.toast(`Backup saved: ${out.file} (${U.fmtBytes(out.bytes)})`, "ok");
      } catch (err) {
        btn.disabled = false;
        progress.done();
        U.toast(err.message, "error");
      }
    }
  }

  function _dragResize(grip, onMove, onDone) {
    grip.addEventListener("mousedown", (ev) => {
      ev.preventDefault();
      const move = (mv) => onMove(mv);
      const up = () => {
        window.removeEventListener("mousemove", move);
        window.removeEventListener("mouseup", up);
        document.body.style.cursor = "";
        if (onDone) onDone();
      };
      window.addEventListener("mousemove", move);
      window.addEventListener("mouseup", up);
    });
  }

  /* =================================================================
   * project open/create
   * ================================================================= */

  function _showProjectModal(recent) {
    const backdrop = document.getElementById("modal-backdrop");
    const modal = document.getElementById("modal-project");
    backdrop.hidden = false;
    modal.hidden = false;

    // never trap the user: backdrop click / Escape close the modal
    backdrop.onclick = (ev) => { if (ev.target === backdrop) _hideProjectModal(); };
    document.addEventListener("keydown", _escClose);

    const list = document.getElementById("recent-list");
    U.clear(list);
    recent = recent.filter((item) => item.exists).slice(0, 8);
    if (recent.length) {
      list.append(U.el("div", { class: "muted", style: "margin-top:10px" }, "Recent projects"));
      for (const item of recent) {
        const btn = U.el("button", { class: `recent-item ${item.exists ? "" : "missing"}` },
          U.el("span", { class: "ri-name" }, item.name),
          U.el("span", { class: "ri-path" }, item.configfile),
        );
        btn.addEventListener("click", () => item.exists && _openPath(item.configfile));
        list.append(btn);
      }
    }

    document.getElementById("btn-open-existing").onclick = async () => {
      const path = await Api.pickFile({
        title: "Open AeoLiS configuration",
        patterns: [["AeoLiS config", "*.txt"], ["All files", "*.*"]],
      }).catch((err) => { U.toast(err.message, "error"); return null; });
      if (path) _openPath(path);
    };

    document.getElementById("btn-new-project").onclick = async () => {
      const folder = await Api.pickFolder({ title: "Choose a folder for the new project" })
        .catch((err) => { U.toast(err.message, "error"); return null; });
      if (!folder) return;
      let info;
      try {
        info = await Api.post("/api/project/new", { folder });
      } catch (err) {
        U.toast(err.message, "error");
        return;
      }
      _hideProjectModal();
      try {
        await _projectOpened(info);
      } catch (err) {
        U.toast(`Project created, but loading state failed: ${err.message}`, "error");
      }
    };

    document.getElementById("btn-open-manual").onclick = () => {
      const path = document.getElementById("manual-path").value.trim();
      if (path) _openPath(path);
    };
  }

  function _escClose(ev) {
    if (ev.key === "Escape") _hideProjectModal();
  }

  function _hideProjectModal() {
    document.getElementById("modal-backdrop").hidden = true;
    document.getElementById("modal-project").hidden = true;
    document.removeEventListener("keydown", _escClose);
  }

  async function _openPath(path) {
    let info;
    try {
      info = await Api.post("/api/project/open", { path });
    } catch (err) {
      U.toast(err.message, "error");
      return;
    }
    // hide the modal as soon as the project is open on the backend, so
    // a hiccup while loading state can never leave the UI blocked
    _hideProjectModal();
    try {
      await _projectOpened(info);
    } catch (err) {
      U.toast(`Project opened, but loading state failed: ${err.message}`, "error");
    }
  }

  async function _projectOpened(info) {
    // a project switch must not leak the previous project's layers,
    // charts or playbar sources into the new one
    for (const layer of App.state.layers) {
      FieldLayer.remove(`field-${layer.id}`);
      MapView.removeLayerAndSource(`pts-${layer.id}`);
    }
    App.state.layers.length = 0;
    Graphs.clearAll();
    Playbar.clearSources();

    App.state.project = info;
    const pathEl = document.getElementById("topbar-path");
    U.clear(pathEl);
    // <bdi> keeps the rtl-ellipsis trick from mirroring the text
    pathEl.append(U.el("bdi", {}, info.configfile));
    pathEl.title = info.configfile;

    // restore persisted UI state
    App.state.crsFromState = false;
    try {
      const st = await Api.get("/api/project/state");
      if (st && st.ui) Object.assign(App.state.ui, st.ui);
      if (st && st.crs) { App.state.crs = st.crs; App.state.crsFromState = true; }
    } catch { /* fresh project */ }

    // clamp persisted sizes so a layout saved on a large monitor never
    // squeezes the map on a smaller screen
    const sbw = U.clamp(App.state.ui.sidebarWidth || 340, 260, window.innerWidth * 0.32);
    const grh = U.clamp(App.state.ui.graphsHeight || 220, 120, window.innerHeight * 0.38);
    App.state.ui.sidebarWidth = sbw;
    App.state.ui.graphsHeight = grh;
    document.getElementById("sidebar").style.width = `${sbw}px`;
    document.documentElement.style.setProperty("--sidebar-live-w", `${sbw}px`);
    document.getElementById("graphs-wrap").style.height = `${grh}px`;
    MapView.setBasemap(App.state.ui.basemap);
    Theme.restoreFromProject();
    App.emit("crs", App.state.crs);

    await Objects.load();
    App.emit("project", info);
    U.toast(`Project ${info.name} opened`, "ok");
  }

})();
