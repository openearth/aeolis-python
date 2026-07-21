/* Layers floating panel: read-only view over registered layers grouped
 * by category, with visibility (eye) toggles. Layers themselves are
 * registered by the tabs/viewer via Layers.register(). */
"use strict";

const LayersPanel = (() => {

  const GROUP_ORDER = ["output", "domain", "rawdata", "grid", "objects", "background"];
  const GROUP_TITLES = {
    output: "Model output",
    domain: "Domain (.grd)",
    rawdata: "Raw data",
    grid: "Grid",
    objects: "Objects",
    background: "Background",
  };

  function refresh() {
    const body = document.getElementById("layers-body");
    if (!body) return;
    U.clear(body);

    const groups = new Map();
    for (const layer of App.state.layers) {
      if (!groups.has(layer.group)) groups.set(layer.group, []);
      groups.get(layer.group).push(layer);
    }

    for (const group of GROUP_ORDER) {
      const layers = groups.get(group);
      if (group === "objects") {
        _objectsGroup(body);
        continue;
      }
      if (group === "background") {
        _backgroundGroup(body);
        continue;
      }
      if (!layers || !layers.length) continue;
      const wrap = U.el("div", { class: "lp-group" });
      wrap.append(U.el("div", { class: "lp-group-head" }, GROUP_TITLES[group] || group));
      for (const layer of layers) wrap.append(_layerRow(layer));
      body.append(wrap);
    }

    if (!body.children.length) {
      body.append(U.el("div", { class: "muted" }, "No layers yet."));
    }
  }

  function _layerRow(layer) {
    const eye = U.el("span", { class: `eye ${layer.visible ? "" : "off"}`, title: "Show/hide" }, "👁");
    eye.addEventListener("click", () => {
      layer.visible = !layer.visible;
      App.emit("layer-visibility", layer);
      refresh();
    });
    const row = U.el("div", { class: "lp-row" },
      eye,
      U.el("span", { class: "lp-name", title: layer.title }, layer.title),
      layer.subtitle ? U.el("span", { class: "lp-mini" }, layer.subtitle) : null,
    );
    row.addEventListener("click", (ev) => {
      if (ev.target === eye) return;
      App.emit("layer-select", layer);
    });
    return row;
  }

  function _objectsGroup(body) {
    if (!App.state.objects.length) return;
    const wrap = U.el("div", { class: "lp-group" });
    wrap.append(U.el("div", { class: "lp-group-head" }, GROUP_TITLES.objects));
    for (const obj of App.state.objects) {
      const eye = U.el("span", { class: `eye ${obj.visible ? "" : "off"}` }, "👁");
      eye.addEventListener("click", () => Objects.update(obj.id, { visible: !obj.visible }));
      wrap.append(U.el("div", { class: "lp-row" },
        eye,
        U.el("span", { class: "lp-name", style: `color:${obj.color}` }, obj.name),
        U.el("span", { class: "lp-mini" }, obj.kind),
      ));
    }
    body.append(wrap);
  }

  function _backgroundGroup(body) {
    if (CRS.isLocal()) return;
    const wrap = U.el("div", { class: "lp-group" });
    wrap.append(U.el("div", { class: "lp-group-head" }, GROUP_TITLES.background));
    for (const [key, label] of [["gray", "Grey map"], ["sat", "Satellite"], ["none", "None"]]) {
      const active = App.state.ui.basemap === key;
      const row = U.el("div", { class: `lp-row ${active ? "selected" : ""}` },
        U.el("span", { class: "eye" }, active ? "●" : "○"),
        U.el("span", { class: "lp-name" }, label));
      row.addEventListener("click", () => MapView.setBasemap(key));
      wrap.append(row);
    }
    body.append(wrap);
  }

  function init() {
    App.on("objects", refresh);
    App.on("layers", refresh);
    App.on("basemap", refresh);
    App.on("crs", refresh);
    refresh();
  }

  return { init, refresh };
})();

/* Central layer registry used by tabs and the viewer. */
const Layers = (() => {

  function register(layer) {
    // layer: {id, group, title, subtitle?, visible, render hooks added later}
    const existing = App.state.layers.findIndex((l) => l.id === layer.id);
    if (existing >= 0) App.state.layers[existing] = { ...App.state.layers[existing], ...layer };
    else App.state.layers.push({ visible: true, ...layer });
    App.emit("layers", layer.id);
  }

  function unregister(id) {
    const idx = App.state.layers.findIndex((l) => l.id === id);
    if (idx >= 0) {
      App.state.layers.splice(idx, 1);
      App.emit("layers", id);
    }
  }

  function get(id) { return App.state.layers.find((l) => l.id === id) || null; }

  return { register, unregister, get };
})();
