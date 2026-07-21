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
    try {
      if (!App.state.schema) App.state.schema = await Api.get("/api/schema");
      const cfg = await Api.get("/api/config");
      App.state.config = cfg.values;
      App.state.configDirty = false;
      loaded = true;

      SchemaForm.build(panel, App.state.schema, App.state.config, _onEdit);
      const first = panel.querySelector(".section");
      if (first) first.classList.remove("collapsed");
      _updateButtons();
    } catch (err) {
      U.clear(panel);
      panel.append(U.el("div", { class: "muted" }, `Could not load configuration: ${err.message}`));
    }
  }

  function _onEdit() {
    App.state.configDirty = true;
    _updateButtons();
  }

  function _updateButtons() {
    document.getElementById("btn-cfg-save").disabled = !loaded || !App.state.configDirty;
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
        const chip = document.getElementById("project-chip");
        chip.textContent = info.name;
        chip.title = info.configfile;
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
