/* Conditions tab: wind / waves / water levels (implemented in Phase 5). */
"use strict";

const ConditionsTab = (() => {
  function init() {
    Tabs.register("conditions", {});
    const panel = document.getElementById("conditions-panel");
    panel.append(U.el("div", { class: "muted" }, "Boundary condition tools arrive in Phase 5."));
  }
  return { init };
})();
