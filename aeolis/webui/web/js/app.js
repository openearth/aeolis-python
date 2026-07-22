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

  Draw.init();

  SettingsTab.init();
  GridTab.init();
  DomainTab.init();
  ConditionsTab.init();
  RunTab.init();
  ViewerTab.init();

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
    App.on("project", () => {
      document.getElementById("btn-cfg-view").disabled = false;
      document.getElementById("btn-cfg-reveal").disabled = false;
    });
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
