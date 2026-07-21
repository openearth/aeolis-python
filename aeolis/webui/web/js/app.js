/* Application bootstrap: wires everything together and handles the
 * project open/create flow. */
"use strict";

(async function main() {

  Tabs.init();
  Playbar.init();
  Graphs.init();

  await MapView.init("map");

  Draw.init();
  LayersPanel.init();
  ObjectsPanel.init();

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
    // sidebar resize
    const sidebar = document.getElementById("sidebar");
    _dragResize(document.getElementById("sidebar-grip"), (ev) => {
      const width = U.clamp(ev.clientX, 260, 640);
      sidebar.style.width = `${width}px`;
      App.state.ui.sidebarWidth = width;
    }, () => App.touchUi());

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

    // floating panel collapse
    for (const header of document.querySelectorAll(".float-panel > header")) {
      header.addEventListener("click", () => {
        header.parentElement.classList.toggle("collapsed");
      });
    }

    // project button
    document.getElementById("btn-project").addEventListener("click", async () => {
      const info = await Api.get("/api/project");
      _showProjectModal(info.recent || []);
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

    const list = document.getElementById("recent-list");
    U.clear(list);
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
      try {
        const info = await Api.post("/api/project/new", { folder });
        await _projectOpened(info);
        _hideProjectModal();
      } catch (err) {
        U.toast(err.message, "error");
      }
    };

    document.getElementById("btn-open-manual").onclick = () => {
      const path = document.getElementById("manual-path").value.trim();
      if (path) _openPath(path);
    };
  }

  function _hideProjectModal() {
    document.getElementById("modal-backdrop").hidden = true;
    document.getElementById("modal-project").hidden = true;
  }

  async function _openPath(path) {
    try {
      const info = await Api.post("/api/project/open", { path });
      await _projectOpened(info);
      _hideProjectModal();
    } catch (err) {
      U.toast(err.message, "error");
    }
  }

  async function _projectOpened(info) {
    App.state.project = info;
    const chip = document.getElementById("project-chip");
    chip.textContent = info.name;
    chip.title = info.configfile;

    // restore persisted UI state
    try {
      const st = await Api.get("/api/project/state");
      if (st && st.ui) Object.assign(App.state.ui, st.ui);
      if (st && st.crs) App.state.crs = st.crs;
    } catch { /* fresh project */ }

    document.getElementById("sidebar").style.width = `${App.state.ui.sidebarWidth}px`;
    document.getElementById("graphs-wrap").style.height = `${App.state.ui.graphsHeight}px`;
    MapView.setBasemap(App.state.ui.basemap);
    App.emit("crs", App.state.crs);

    await Objects.load();
    App.emit("project", info);
    U.toast(`Project ${info.name} opened`, "ok");
  }

})();
