/* Album-cover theming.
 *
 * Every surface routes through the CSS design tokens in :root; a theme is
 * just a token override block keyed on html[data-theme]. Palettes are
 * inspired by iconic album covers (no copyrighted artwork embedded).
 *
 * The chosen theme is a *global* preference: mirrored to localStorage so
 * the app opens themed before any project loads, and also into the
 * per-project ui state so it travels with the project.
 */
"use strict";

const Theme = (() => {

  const LS_KEY = "aeolis.theme";

  // Plain themes (light/dark) show a minimal colour chip in the picker;
  // album themes show their actual cover art (served from web/covers/).
  const THEMES = [
    { id: "light",     name: "Light",                     sub: "AeoLiS default", chip: "#f4f5f7" },
    { id: "dark",      name: "Dark",                      sub: "Night mode",     chip: "#12181f" },
    { id: "discovery", name: "Discovery",                 sub: "Daft Punk",      cover: "covers/discovery.jpg" },
    { id: "rumours",   name: "Rumours",                   sub: "Fleetwood Mac",  cover: "covers/rumours.jpg" },
    { id: "nevermind", name: "Nevermind",                 sub: "Nirvana",        cover: "covers/nevermind.jpg" },
    { id: "velvet",    name: "The Velvet Underground",    sub: "& Nico",         cover: "covers/velvet.jpg" },
    { id: "darkside",  name: "The Dark Side of the Moon", sub: "Pink Floyd",     cover: "covers/darkside.jpg" },
  ];

  let current = "light";
  let menu = null;

  function byId(id) { return THEMES.find((t) => t.id === id) || THEMES[0]; }

  /* Apply a theme id to the document (light = no attribute). */
  function apply(id, persist = true) {
    const theme = byId(id);
    current = theme.id;
    if (theme.id === "light") delete document.documentElement.dataset.theme;
    else document.documentElement.dataset.theme = theme.id;
    document.getElementById("btn-theme")
      ?.setAttribute("title", `Theme: ${theme.name} — click to change`);
    if (persist) {
      try { localStorage.setItem(LS_KEY, theme.id); } catch { /* ignore */ }
      if (App.state.ui) { App.state.ui.theme = theme.id; App.touchUi(); }
    }
    App.emit("theme", theme.id);
    if (menu) renderMenu();
  }

  function init() {
    // apply the last global choice immediately (before any project opens)
    let saved = "light";
    try { saved = localStorage.getItem(LS_KEY) || "light"; } catch { /* ignore */ }
    apply(saved, false);

    const btn = document.getElementById("btn-theme");
    if (btn) btn.addEventListener("click", (ev) => { ev.stopPropagation(); toggleMenu(btn); });
  }

  /* When a project is opened its stored theme wins over the global one. */
  function restoreFromProject() {
    const id = App.state.ui && App.state.ui.theme;
    if (id && id !== current) apply(id, true);
  }

  function toggleMenu(anchor) {
    if (menu) { closeMenu(); return; }
    anchor.classList.add("on");
    menu = U.el("div", { class: "theme-menu" });
    renderMenu();
    document.body.append(menu);
    const r = anchor.getBoundingClientRect();
    // right-align the menu under the button, kept inside the viewport
    const width = menu.offsetWidth || 260;
    menu.style.top = `${r.bottom + 6}px`;
    menu.style.left = `${Math.max(6, Math.min(r.right - width, window.innerWidth - width - 6))}px`;
    // close on any outside interaction (defer so this very click doesn't fire it)
    setTimeout(() => {
      document.addEventListener("mousedown", onOutside, true);
      document.addEventListener("keydown", onEsc, true);
    }, 0);
  }

  function renderMenu() {
    if (!menu) return;
    U.clear(menu);
    menu.append(U.el("div", { class: "theme-menu-head" }, "Theme"));
    for (const t of THEMES) {
      // album themes show their cover; light/dark show a plain colour chip
      const art = t.cover
        ? U.el("img", { class: "theme-cover", src: t.cover, alt: t.name, loading: "lazy" })
        : U.el("span", { class: "theme-chip", style: `background:${t.chip}` });
      const item = U.el("button", {
        class: `theme-item ${t.cover ? "album" : "plain"} ${t.id === current ? "on" : ""}`,
      },
        art,
        U.el("span", { class: "theme-label" },
          U.el("span", { class: "theme-name" }, t.name),
          U.el("span", { class: "theme-sub" }, t.sub)),
        t.id === current ? U.icon("check", 15) : null);
      item.addEventListener("click", () => { apply(t.id); closeMenu(); });
      menu.append(item);
    }
  }

  function onOutside(ev) {
    if (menu && !menu.contains(ev.target) &&
        !ev.target.closest("#btn-theme")) closeMenu();
  }
  function onEsc(ev) { if (ev.key === "Escape") closeMenu(); }

  function closeMenu() {
    if (menu) { menu.remove(); menu = null; }
    document.getElementById("btn-theme")?.classList.remove("on");
    document.removeEventListener("mousedown", onOutside, true);
    document.removeEventListener("keydown", onEsc, true);
  }

  return { init, apply, restoreFromProject, THEMES };
})();
