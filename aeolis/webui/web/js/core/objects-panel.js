/* Objects floating panel: manage user-drawn shapes (rename, recolor,
 * delete, zoom-to). Editing vertices happens via the draw tools (Phase 3+). */
"use strict";

const ObjectsPanel = (() => {

  let selected = null;

  function refresh() {
    const body = document.getElementById("objects-body");
    if (!body) return;
    U.clear(body);

    if (!App.state.objects.length) {
      body.append(U.el("div", { class: "muted" }, "No objects yet. Draw polygons or transects from the tab panels."));
      return;
    }

    for (const obj of App.state.objects) {
      const row = U.el("div", { class: `lp-row ${selected === obj.id ? "selected" : ""}` });
      const swatch = U.el("span", {
        class: "eye", title: "Change colour",
        style: `color:${obj.color}`,
      }, "■");
      swatch.addEventListener("click", () => _cycleColor(obj));

      const name = U.el("span", { class: "lp-name", title: "Double-click to rename" }, obj.name);
      name.addEventListener("dblclick", () => _rename(obj, name));

      const zoom = U.el("span", { class: "lp-mini", style: "cursor:pointer", title: "Zoom to" }, "⌖");
      zoom.addEventListener("click", () => _zoomTo(obj));

      const del = U.el("span", { class: "lp-mini", style: "cursor:pointer", title: "Delete" }, "✕");
      del.addEventListener("click", () => {
        if (window.confirm(`Delete ${obj.name}?`)) Objects.remove(obj.id);
      });

      row.addEventListener("click", () => {
        selected = obj.id;
        App.emit("object-select", obj.id);
        refresh();
      });

      row.append(swatch, name, zoom, del);
      body.append(row);
    }
  }

  function _rename(obj, nameNode) {
    const input = U.el("input", { type: "text", value: obj.name, style: "flex:1;font-size:12px" });
    nameNode.replaceWith(input);
    input.focus();
    input.select();
    const commit = () => Objects.update(obj.id, { name: input.value.trim() || obj.name });
    input.addEventListener("blur", commit);
    input.addEventListener("keydown", (ev) => {
      if (ev.key === "Enter") input.blur();
      if (ev.key === "Escape") { input.value = obj.name; input.blur(); }
    });
  }

  const PALETTE = ["#e6552f", "#2f7fe6", "#27a355", "#a034c6", "#e0a020", "#12a5b5", "#d1387f"];
  function _cycleColor(obj) {
    const idx = PALETTE.indexOf(obj.color);
    Objects.update(obj.id, { color: PALETTE[(idx + 1) % PALETTE.length] });
  }

  function _zoomTo(obj) {
    const xs = obj.coords.map((c) => c[0]);
    const ys = obj.coords.map((c) => c[1]);
    MapView.fitModelBounds(Math.min(...xs), Math.min(...ys), Math.max(...xs), Math.max(...ys));
  }

  function init() {
    App.on("objects", refresh);
    refresh();
  }

  return { init, refresh, selectedId: () => selected };
})();
