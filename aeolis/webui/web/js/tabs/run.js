/* Run tab: model execution (implemented in Phase 6). */
"use strict";

const RunTab = (() => {
  function init() {
    Tabs.register("run", {});
    const panel = document.getElementById("run-panel");
    panel.append(U.el("div", { class: "muted" }, "Run controls arrive in Phase 6."));
  }
  return { init };
})();
