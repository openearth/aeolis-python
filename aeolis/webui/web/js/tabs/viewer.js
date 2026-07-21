/* Viewer tab: layers, styling, netCDF output (implemented in Phase 7). */
"use strict";

const ViewerTab = (() => {
  function init() {
    Tabs.register("viewer", {});
    const panel = document.getElementById("viewer-panel");
    panel.append(U.el("div", { class: "muted" }, "The output viewer arrives in Phase 7."));
  }
  return { init };
})();
