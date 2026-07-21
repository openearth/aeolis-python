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

  return { el, clear, toast, showTip, hideTip, fmtNum, fmtBytes, fmtDuration, fmtDate, debounce, clamp, numField };
})();
