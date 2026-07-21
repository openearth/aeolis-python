/* Settings tab: aeolis.txt editor (implemented in Phase 2). */
"use strict";

const SettingsTab = (() => {
  function init() {
    Tabs.register("settings", {});
    const panel = document.getElementById("settings-form");
    panel.append(U.el("div", { class: "muted" }, "Open a project to edit its configuration."));
  }
  return { init };
})();
