/* Popup infrastructure.
 *
 * Popup.open() — generic modal popup used by the Domain/Conditions
 * wizards; supports hide()/show() so a popup can step aside while the
 * user draws on the map.
 * DocsPopup.open(url) — floating scrollable documentation panel
 * (readthedocs pages embed fine: no X-Frame-Options).
 * TimeTool.open() — date/duration → seconds helper for time settings.
 */
"use strict";

const Popup = (() => {

  function open({ title = "", width = 560, onClose = null, closable = true } = {}) {
    const body = U.el("div", { class: "popup-body" });
    const closeBtn = closable
      ? U.el("button", { class: "ghost popup-close", title: "Close" }, "✕")
      : null;
    const box = U.el("div", { class: "popup-box", style: `width:${width}px` },
      U.el("header", { class: "popup-head" }, U.el("span", {}, title), closeBtn),
      body);
    const backdrop = U.el("div", { class: "popup-backdrop" }, box);

    const api = {
      body,
      close() {
        backdrop.remove();
        document.removeEventListener("keydown", onKey);
        if (onClose) onClose();
      },
      hide() { backdrop.style.display = "none"; },
      show() { backdrop.style.display = "flex"; },
      setTitle(text) { box.querySelector(".popup-head span").textContent = text; },
    };

    const onKey = (ev) => {
      if (ev.key === "Escape" && closable) { ev.stopPropagation(); api.close(); }
    };
    if (closeBtn) closeBtn.addEventListener("click", api.close);
    backdrop.addEventListener("mousedown", (ev) => {
      if (ev.target === backdrop && closable) api.close();
    });
    document.addEventListener("keydown", onKey);
    document.body.append(backdrop);
    return api;
  }

  return { open };
})();


const DocsPopup = (() => {

  let panel = null;

  function open(url, title = "AeoLiS documentation") {
    close();
    // proxied same-origin so readthedocs rate-limit pages (which set
    // X-Frame-Options) can never break the embed
    const frame = U.el("iframe", {
      src: `/api/docs?url=${encodeURIComponent(url)}`, class: "docs-frame",
    });
    const external = U.el("a", {
      href: url, target: "_blank", class: "muted",
      style: "font-size:11.5px;margin-right:8px",
    }, "open in browser ⧉");
    const closeBtn = U.el("button", { class: "ghost popup-close" }, "✕");
    closeBtn.addEventListener("click", close);
    panel = U.el("div", { class: "docs-popup" },
      U.el("header", { class: "popup-head" },
        U.el("span", {}, title),
        U.el("span", {}, external, closeBtn)),
      frame);
    document.body.append(panel);
  }

  function close() {
    if (panel) { panel.remove(); panel = null; }
  }

  return { open, close };
})();


const TimeTool = (() => {

  /* Helper to compute seconds for time parameters.
   * refdateEpoch: epoch seconds of the config refdate.
   * onApply(seconds) writes the value back into the form. */
  function open(paramKey, refdateEpoch, onApply) {
    const popup = Popup.open({ title: `Compute ${paramKey} [s]`, width: 480 });

    const result = U.el("input", { type: "text", readonly: "", style: "font-weight:700" });
    const setResult = (seconds) => {
      if (Number.isFinite(seconds)) result.value = String(Math.round(seconds));
    };

    // --- from a calendar date (relative to refdate) ---
    const dateInput = U.el("input", { type: "datetime-local", step: 60 });
    const dateBtn = U.el("button", { class: "ghost" }, "→ seconds since refdate");
    dateBtn.addEventListener("click", () => {
      if (!dateInput.value) return;
      const epoch = Date.parse(dateInput.value + "Z") / 1000;
      setResult(epoch - refdateEpoch);
    });

    // --- from a duration ---
    const amount = U.el("input", { type: "text", value: "1", style: "width:80px" });
    const unit = U.el("select", {},
      ...[["hours", 3600], ["days", 86400], ["weeks", 7 * 86400],
        ["months (30 d)", 30 * 86400], ["years (365 d)", 365 * 86400]]
        .map(([label, s]) => U.el("option", { value: s }, label)));
    const durBtn = U.el("button", { class: "ghost" }, "→ seconds");
    durBtn.addEventListener("click", () => {
      const v = Number(amount.value);
      if (Number.isFinite(v)) setResult(v * Number(unit.value));
    });

    const applyBtn = U.el("button", { class: "primary" }, `Apply to ${paramKey}`);
    applyBtn.addEventListener("click", () => {
      const v = Number(result.value);
      if (Number.isFinite(v)) { onApply(v); popup.close(); }
    });

    popup.body.append(
      U.el("div", { class: "muted", style: "margin-bottom:8px" },
        `refdate: ${U.fmtDate(refdateEpoch)} (UTC) — times are seconds since refdate`),
      U.el("div", { class: "form-group" },
        U.el("span", { class: "fg-label" }, "From a date"),
        U.el("div", { class: "form-row" }, dateInput, dateBtn)),
      U.el("div", { class: "form-group" },
        U.el("span", { class: "fg-label" }, "From a duration"),
        U.el("div", { class: "form-row" }, amount, unit, durBtn)),
      U.el("div", { class: "form-group" },
        U.el("span", { class: "fg-label" }, "Result [s]"),
        U.el("div", { class: "form-row" }, result, applyBtn)),
    );
  }

  return { open };
})();
