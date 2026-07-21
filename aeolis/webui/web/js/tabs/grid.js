/* Grid tab: model grid generation (implemented in Phase 3). */
"use strict";

const GridTab = (() => {
  function init() {
    Tabs.register("grid", {});
    const panel = document.getElementById("grid-panel");
    panel.append(U.el("div", { class: "muted" }, "Grid generation tools arrive in Phase 3."));
  }
  return { init };
})();
