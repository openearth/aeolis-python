/* Settings tab: view and modify the full aeolis.txt configuration.
 *
 * The form is generated from /api/schema (parsed from constants.py) so
 * it always matches the installed AeoLiS version. Load / Save / Save-as
 * go through aeolis.inout on the backend.
 */
"use strict";

const SettingsTab = (() => {

  let loaded = false;

  function init() {
    Tabs.register("settings", { enter: _enter });

    document.getElementById("btn-cfg-save").addEventListener("click", () => _save());
    document.getElementById("btn-cfg-saveas").addEventListener("click", _saveAs);
    document.getElementById("cfg-search").addEventListener("input",
      U.debounce((ev) => SchemaForm.filter(ev.target.value), 120));

    App.on("project", () => { loaded = false; _load(); });
    App.on("config-changed", (key) => SchemaForm.refreshParam(key));

    const panel = document.getElementById("settings-form");
    panel.append(U.el("div", { class: "muted" }, "Open a project to edit its configuration."));
  }

  function _enter() {
    if (!loaded && App.state.project) _load();
  }

  async function _load() {
    if (!App.state.project) return;
    const panel = document.getElementById("settings-form");
    const restoreScroll = U.keepScroll(panel);
    try {
      if (!App.state.schema) App.state.schema = await Api.get("/api/schema");
      const cfg = await Api.get("/api/config");
      App.state.config = cfg.values;
      App.state.configDirty = false;
      loaded = true;

      SchemaForm.build(panel, App.state.schema, App.state.config, _onEdit);
      const first = panel.querySelector(".section");
      if (first) first.classList.remove("collapsed");
      _injectShearToggle(panel);
      panel.prepend(_prefsSection().wrap);
      _updateButtons();
      restoreScroll();   // form landed after awaits — restore again
    } catch (err) {
      U.clear(panel);
      panel.append(U.el("div", { class: "muted" }, `Could not load configuration: ${err.message}`));
    }
  }

  function _onEdit() {
    App.state.configDirty = true;
    _updateButtons();
  }

  /* GUI preferences (not aeolis parameters): how file references are
   * written into aeolis.txt when the GUI saves/repoints a file. */
  function _prefsSection() {
    const section = U.section("GUI preferences", { collapsed: true });
    const sel = U.el("select", {});
    for (const [v, label] of [["relative", "relative (portable — recommended)"],
                              ["absolute", "absolute (full Windows path)"]]) {
      sel.append(U.el("option", {
        value: v, selected: v === (App.state.ui.pathFormat || "relative") ? "" : null,
      }, label));
    }
    sel.addEventListener("change", () => {
      App.state.ui.pathFormat = sel.value;
      App.touchUi();
    });
    section.body.append(
      U.el("div", { class: "form-row" }, U.el("label", {}, "File references"), sel),
      U.el("div", { class: "muted", style: "font-size:11.5px" },
        "How the GUI writes file paths into aeolis.txt when saving files. "
        + "Relative references use forward slashes, so the project also runs "
        + "unchanged on the cluster (/p); absolute Windows paths are converted "
        + "to their /p form automatically when the project is copied for an HPC run."),
    );
    return section;
  }

  /* Mirror the grid tab's shear-grid show/hide eye into the shear section
   * header here, so the computational grid can be toggled from Settings too. */
  function _injectShearToggle(panel) {
    if (typeof GridTab === "undefined" || !GridTab.setShearVisible) return;
    const head = panel.querySelector('[data-section="Topographic steering (shear)"] > header');
    if (!head || head.querySelector(".shear-eye-settings")) return;
    const vis = GridTab.shearVisible();
    const eye = U.el("span", {
      class: `eye group-eye shear-eye-settings ${vis ? "" : "off"}`,
      title: "Show/hide the computational (shear) grid on the map",
    }, "👁");
    eye.addEventListener("click", (ev) => {
      ev.stopPropagation();   // don't collapse the section
      const next = !GridTab.shearVisible();
      GridTab.setShearVisible(next);
      eye.classList.toggle("off", !next);
    });
    const count = head.querySelector(".count");
    head.insertBefore(eye, count);
  }

  function _updateButtons() {
    const saveBtn = document.getElementById("btn-cfg-save");
    saveBtn.disabled = !loaded || !App.state.configDirty;
    // unsaved-changes dot on the save button
    saveBtn.classList.toggle("dirty", loaded && Boolean(App.state.configDirty));
    saveBtn.title = App.state.configDirty
      ? "Save configuration (unsaved changes)" : "Save configuration";
    document.getElementById("btn-cfg-saveas").disabled = !loaded;
  }

  async function _save(path = null) {
    try {
      const res = await Api.post("/api/config/save", {
        values: App.state.config,
        ...(path ? { path } : {}),
      });
      App.state.configDirty = false;
      _updateButtons();
      U.toast(`Saved ${res.path}`, "ok");
      if (path) {
        const info = await Api.get("/api/project");
        App.state.project = info;
        const pathEl = document.getElementById("topbar-path");
        U.clear(pathEl);
        pathEl.append(U.el("bdi", {}, info.configfile));
        pathEl.title = info.configfile;
      }
    } catch (err) {
      U.toast(`Save failed: ${err.message}`, "error");
    }
  }

  async function _saveAs() {
    const path = await Api.pickFile({
      title: "Save configuration as",
      save: true,
      patterns: [["AeoLiS config", "*.txt"]],
      initial: App.state.project ? App.state.project.root : "",
    }).catch((err) => { U.toast(err.message, "error"); return null; });
    if (path) _save(path);
  }

  /* Other tabs call this after writing file parameters (e.g. grid
   * generation sets xgrid_file/ygrid_file) so the form and the file on
   * disk stay in sync. */
  async function setConfigValues(patch, save = true) {
    Object.assign(App.state.config, patch);
    for (const key of Object.keys(patch)) App.emit("config-changed", key);
    if (save) await _save();
    else { App.state.configDirty = true; _updateButtons(); }
  }

  return { init, setConfigValues };
})();
