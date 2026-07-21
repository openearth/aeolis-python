/* Domain tab: bathymetry / vegetation / ne-layer (implemented in Phase 4). */
"use strict";

const DomainTab = (() => {
  function init() {
    Tabs.register("domain", {});
    const panel = document.getElementById("domain-panel");
    panel.append(U.el("div", { class: "muted" }, "Domain data tools arrive in Phase 4."));
  }
  return { init };
})();
