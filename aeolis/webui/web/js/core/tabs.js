/* Tab controller: top tab bar <-> left panels, with per-tab lifecycle hooks. */
"use strict";

const Tabs = (() => {

  const hooks = new Map();   // name -> {enter, leave}
  let active = "settings";

  function register(name, tabHooks) { hooks.set(name, tabHooks || {}); }

  function activate(name) {
    if (name === active) return;
    const prev = hooks.get(active);
    if (prev && prev.leave) prev.leave();

    active = name;
    document.body.className = `tab-${name}`;

    for (const btn of document.querySelectorAll("#tabbar .tab-btn")) {
      btn.classList.toggle("active", btn.dataset.tab === name);
    }
    for (const panel of document.querySelectorAll("#sidebar .tabpanel")) {
      panel.hidden = panel.dataset.tab !== name;
    }

    const next = hooks.get(name);
    if (next && next.enter) next.enter();
    App.emit("tab", name);
  }

  function current() { return active; }

  function init() {
    for (const btn of document.querySelectorAll("#tabbar .tab-btn")) {
      btn.addEventListener("click", () => activate(btn.dataset.tab));
    }
  }

  return { register, activate, current, init };
})();
