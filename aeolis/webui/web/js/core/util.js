/* Generic helpers: DOM, formatting, toasts, tooltip. */
"use strict";

const U = (() => {

  function el(tag, attrs = {}, ...children) {
    const node = document.createElement(tag);
    for (const [k, v] of Object.entries(attrs)) {
      if (k === "class") node.className = v;
      else if (k === "dataset") Object.assign(node.dataset, v);
      else if (k.startsWith("on") && typeof v === "function") {
        node.addEventListener(k.slice(2), v);
      } else if (v !== null && v !== undefined) node.setAttribute(k, v);
    }
    for (const child of children.flat()) {
      if (child === null || child === undefined) continue;
      node.append(child.nodeType ? child : document.createTextNode(child));
    }
    return node;
  }

  function clear(node) { while (node.firstChild) node.removeChild(node.firstChild); }

  function toast(message, kind = "") {
    const box = document.getElementById("toasts");
    const t = el("div", { class: `toast ${kind}` }, message);
    box.append(t);
    setTimeout(() => { t.style.opacity = "0"; t.style.transition = "opacity .3s"; }, 3200);
    setTimeout(() => t.remove(), 3600);
  }

  /* ---- shared tooltip ---- */
  let tipNode = null;
  function showTip(html, x, y) {
    if (!tipNode) {
      tipNode = el("div", { id: "tooltip" });
      document.body.append(tipNode);
    }
    tipNode.innerHTML = html;
    tipNode.style.display = "block";
    const rect = tipNode.getBoundingClientRect();
    const left = Math.min(x + 14, window.innerWidth - rect.width - 8);
    const top = Math.min(y + 14, window.innerHeight - rect.height - 8);
    tipNode.style.left = `${Math.max(4, left)}px`;
    tipNode.style.top = `${Math.max(4, top)}px`;
  }
  function hideTip() { if (tipNode) tipNode.style.display = "none"; }

  /* ---- formatting ---- */
  function fmtNum(v, digits = 3) {
    if (v === null || v === undefined || Number.isNaN(v)) return "–";
    if (v === 0) return "0";
    const a = Math.abs(v);
    if (a >= 1e6 || a < 1e-3) return v.toExponential(digits - 1);
    return Number(v.toFixed(digits)).toString();
  }

  function fmtBytes(n) {
    if (!Number.isFinite(n)) return "–";
    const units = ["B", "kB", "MB", "GB", "TB"];
    let i = 0;
    while (n >= 1024 && i < units.length - 1) { n /= 1024; i += 1; }
    return `${n.toFixed(n >= 10 || i === 0 ? 0 : 1)} ${units[i]}`;
  }

  function fmtDuration(seconds) {
    if (!Number.isFinite(seconds)) return "–";
    seconds = Math.max(0, Math.round(seconds));
    const h = Math.floor(seconds / 3600);
    const m = Math.floor((seconds % 3600) / 60);
    const s = seconds % 60;
    if (h > 48) return `${(h / 24).toFixed(1)} d`;
    return `${String(h).padStart(2, "0")}:${String(m).padStart(2, "0")}:${String(s).padStart(2, "0")}`;
  }

  function fmtDate(epochSeconds) {
    if (!Number.isFinite(epochSeconds)) return "–";
    const d = new Date(epochSeconds * 1000);
    const pad = (v) => String(v).padStart(2, "0");
    return `${d.getUTCFullYear()}-${pad(d.getUTCMonth() + 1)}-${pad(d.getUTCDate())} ` +
      `${pad(d.getUTCHours())}:${pad(d.getUTCMinutes())}`;
  }

  function debounce(fn, ms) {
    let handle = null;
    return (...args) => {
      clearTimeout(handle);
      handle = setTimeout(() => fn(...args), ms);
    };
  }

  function clamp(v, lo, hi) { return Math.min(hi, Math.max(lo, v)); }

  /* ---- shared inline SVG icons (17x17 monochrome) ---- */
  const ICONS = {
    draw: "M3 17.2 13.9 6.3l3.8 3.8L6.8 21H3v-3.8zM19.7 8 16 4.2l1.6-1.6a1 1 0 0 1 1.4 0l2.4 2.4a1 1 0 0 1 0 1.4L19.7 8z",
    edit: "M12 2l3.2 3.2h-2.2v4.6h4.6V7.6L20.8 12l-3.2 3.2v-2.2h-4.6v4.6h2.2L12 20.8l-3.2-3.2h2.2v-4.6H6.4v2.2L3.2 12l3.2-3.2v2.2H11V5.2H8.8L12 2z",
    save: "M5 3h11l5 5v12a1 1 0 0 1-1 1H5a1 1 0 0 1-1-1V4a1 1 0 0 1 1-1zm2 2v5h9V5H7zm10 15v-7H7v7h10z",
    download: "M12 3v10.2l3.6-3.6 1.4 1.4-6 6-6-6 1.4-1.4L10 13.2V3h2zM4 19h16v2H4v-2z",
    upload: "M12 21V10.8l-3.6 3.6L7 13l6-6 6 6-1.4 1.4-3.6-3.6V21h-2zM4 3h16v2H4V3z",
    check: "M9.2 16.6 4.8 12.2l-1.6 1.6 6 6L21 8l-1.6-1.6-10.2 10.2z",
    trash: "M9 3h6l1 2h5v2H3V5h5l1-2zM5 9h14l-1 12a1 1 0 0 1-1 .9H7a1 1 0 0 1-1-.9L5 9zm4 3v7h2v-7H9zm4 0v7h2v-7h-2z",
    book: "M5 3h13a2 2 0 0 1 2 2v14a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2V5a2 2 0 0 1 2-2zm0 2v12.5A3.5 3.5 0 0 1 6.5 17H18V5H5zm1 14a1 1 0 0 0 0 2h12v-2H6.5H6z",
    gear: "M12 8.5A3.5 3.5 0 1 0 12 15.5 3.5 3.5 0 0 0 12 8.5zM20 12c0-.5 0-1-.1-1.4l2-1.6-2-3.4-2.4 1a8 8 0 0 0-2.4-1.4L14.7 2h-4l-.4 2.6a8 8 0 0 0-2.4 1.4l-2.4-1-2 3.4 2 1.6a8 8 0 0 0 0 2.8l-2 1.6 2 3.4 2.4-1a8 8 0 0 0 2.4 1.4l.4 2.6h4l.4-2.6a8 8 0 0 0 2.4-1.4l2.4 1 2-3.4-2-1.6c.1-.4.1-.9.1-1.4z",
    layers: "M12 3 2 9l10 6 10-6-10-6zm-6.5 9.9L2 15l10 6 10-6-3.5-2.1L12 17l-6.5-4.1z",
    wand: "M6 3l1 2.4L9.5 6 7 7l-1 2.4L5 7 2.5 6 5 5.4 6 3zm12.7 2.3a1 1 0 0 1 0 1.4L8.4 17 7 15.6 17.3 5.3a1 1 0 0 1 1.4 0zM19 12l.8 1.8 1.7.7-1.7.8L19 17l-.8-1.7-1.7-.8 1.7-.7L19 12z",
    modify: "M14.1 5.9 18.1 9.9 8 20H4v-4L14.1 5.9zm2.8-2.8a1 1 0 0 1 1.4 0l2.6 2.6a1 1 0 0 1 0 1.4l-1.8 1.8-4-4 1.8-1.8z",
    copy: "M8 2h11a1 1 0 0 1 1 1v13h-2V4H8V2zM4 6h11a1 1 0 0 1 1 1v14a1 1 0 0 1-1 1H4a1 1 0 0 1-1-1V7a1 1 0 0 1 1-1zm1 2v12h9V8H5z",
    up: "M12 5l7 8h-4.5v6h-5v-6H5l7-8z",
    down: "M12 19l-7-8h4.5V5h5v6H19l-7 8z",
    eye: "M12 5c5 0 9 4.4 10 7-1 2.6-5 7-10 7S3 14.6 2 12c1-2.6 5-7 10-7zm0 3.2A3.8 3.8 0 1 0 12 15.8 3.8 3.8 0 0 0 12 8.2z",
    interp: "M3 3h4.5v4.5H3V3zm13.5 0H21v4.5h-4.5V3zM3 16.5h4.5V21H3v-4.5zm13.5 0H21V21h-4.5v-4.5zM9.8 9.8h4.4v4.4H9.8V9.8z",
    clock: "M12 2a10 10 0 1 0 0 20 10 10 0 0 0 0-20zm0 2a8 8 0 1 1 0 16 8 8 0 0 1 0-16zm1 3h-2v6l4.5 2.7 1-1.6-3.5-2.1V7z",
    open: "M4 4h6l2 2h8a1 1 0 0 1 1 1v2H3V5a1 1 0 0 1 1-1zm-1 6h19l-2.2 9.2a1 1 0 0 1-1 .8H5.2a1 1 0 0 1-1-.8L3 10z",
  };

  function icon(name, size = 17) {
    const svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");
    svg.setAttribute("viewBox", "0 0 24 24");
    svg.setAttribute("width", size);
    svg.setAttribute("height", size);
    const path = document.createElementNS("http://www.w3.org/2000/svg", "path");
    path.setAttribute("fill", "currentColor");
    path.setAttribute("d", ICONS[name] || ICONS.gear);
    svg.append(path);
    return svg;
  }

  /* Toolbar button: icon with a small text label below.
   * opts: {primary, toggle, title, onclick} */
  function tbtn(iconName, label, opts = {}) {
    const btn = el("button", {
      class: `tbtn ${opts.primary ? "primary" : ""} ${opts.toggle ? "toggle" : ""}`,
      title: opts.title || label,
    }, icon(iconName), el("span", { class: "tbtn-label" }, label));
    if (opts.onclick) btn.addEventListener("click", opts.onclick);
    return btn;
  }

  /* Small square icon button for list rows. */
  function miniBtn(iconName, title, onclick) {
    const btn = el("button", { class: "mini-btn", title }, icon(iconName, 13));
    if (onclick) btn.addEventListener("click", onclick);
    return btn;
  }

  /* Collapsible section (same look as the Settings sections).
   * Returns {wrap, body, head}. */
  function section(title, { collapsed = false, count = null } = {}) {
    const body = el("div", { class: "section-body" });
    const head = el("header", {},
      el("span", { class: "caret" }, "▾"), title,
      count !== null ? el("span", { class: "count" }, String(count)) : null);
    const wrap = el("div", { class: `section ${collapsed ? "collapsed" : ""}` }, head, body);
    head.addEventListener("click", () => wrap.classList.toggle("collapsed"));
    return { wrap, body, head };
  }

  /* Buffered numeric input: commits on Enter/blur only. */
  function numField(value, onCommit, attrs = {}) {
    const input = el("input", { type: "text", value: value ?? "", ...attrs });
    const commit = () => {
      const parsed = parseFloat(input.value);
      if (Number.isFinite(parsed)) onCommit(parsed);
      else input.value = value ?? "";
    };
    input.addEventListener("blur", commit);
    input.addEventListener("keydown", (ev) => {
      if (ev.key === "Enter") { commit(); input.blur(); }
      if (ev.key === "Escape") { input.value = value ?? ""; input.blur(); }
    });
    return input;
  }

  return { el, clear, toast, showTip, hideTip, fmtNum, fmtBytes, fmtDuration, fmtDate,
    debounce, clamp, numField, icon, tbtn, miniBtn, section };
})();
