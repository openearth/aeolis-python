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
      // while the popup steps aside (map draw / station pick) Escape
      // belongs to that interaction, not to the hidden popup
      if (backdrop.style.display === "none") return;
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

  // parameters that are points in time (a date makes sense) vs plain
  // durations (dt, output interval, restart interval)
  const POINT_PARAMS = new Set(["tstart", "tstop"]);

  /* Helper to compute seconds for time parameters. The result updates
   * live while typing; one Apply button at the bottom. */
  function open(paramKey, refdateEpoch, onApply) {
    const popup = Popup.open({ title: `Set ${paramKey}`, width: 440 });
    const isPoint = POINT_PARAMS.has(paramKey);
    let seconds = null;

    const resultLine = U.el("div", {
      style: "font-size:15px;font-weight:700;margin:2px 0 0",
    }, "–");
    const resultNote = U.el("div", { class: "muted", style: "font-size:12px" });

    const setResult = (value) => {
      seconds = Number.isFinite(value) ? Math.round(value) : null;
      applyBtn.disabled = seconds === null;
      if (seconds === null) { resultLine.textContent = "–"; resultNote.textContent = ""; return; }
      resultLine.textContent = `${paramKey} = ${seconds.toLocaleString("en-US").replace(/,/g, " ")} s`;
      const human = U.fmtDuration(Math.abs(seconds));
      resultNote.textContent = isPoint
        ? `= ${human} after refdate → ${U.fmtDate(refdateEpoch + seconds)} UTC`
        : `= ${human}`;
    };

    // --- as a calendar date/time (points in time only) ---
    const dateInput = U.el("input", { type: "datetime-local", step: 60, class: "grow" });
    dateInput.addEventListener("input", () => {
      if (!dateInput.value) return;
      dateRadio.checked = true;
      setResult(Date.parse(dateInput.value + "Z") / 1000 - refdateEpoch);
    });

    // --- as a duration ---
    const amount = U.el("input", { type: "text", value: "", placeholder: "e.g. 30", style: "width:80px" });
    const unit = U.el("select", {},
      ...[["seconds", 1], ["minutes", 60], ["hours", 3600], ["days", 86400],
        ["weeks", 7 * 86400], ["months (30 d)", 30 * 86400], ["years (365 d)", 365 * 86400]]
        .map(([label, s]) => U.el("option", { value: s, selected: s === 86400 ? "" : null }, label)));
    const durChanged = () => {
      const v = Number(amount.value);
      if (!Number.isFinite(v) || amount.value.trim() === "") return;
      durRadio.checked = true;
      setResult(v * Number(unit.value));
    };
    amount.addEventListener("input", durChanged);
    unit.addEventListener("change", durChanged);

    const dateRadio = U.el("input", { type: "radio", name: "tt-mode", id: "tt-date" });
    const durRadio = U.el("input", { type: "radio", name: "tt-mode", id: "tt-dur" });
    (isPoint ? dateRadio : durRadio).checked = true;

    const applyBtn = U.el("button", { class: "primary", disabled: "" }, "Apply");
    applyBtn.addEventListener("click", () => {
      if (seconds !== null) { onApply(seconds); popup.close(); }
    });
    const cancelBtn = U.el("button", { class: "ghost" }, "Cancel");
    cancelBtn.addEventListener("click", popup.close);

    popup.body.append(
      U.el("div", { class: "muted", style: "font-size:12px;margin-bottom:10px" },
        isPoint
          ? `Times are in seconds since the refdate (${U.fmtDate(refdateEpoch)} UTC).`
          : `${paramKey} is a duration in seconds.`),
      isPoint ? U.el("div", { class: "choice-row" },
        dateRadio, U.el("label", { for: "tt-date" }, "At date / time"), dateInput) : null,
      U.el("div", { class: "choice-row" },
        durRadio, U.el("label", { for: "tt-dur" }, isPoint ? "After refdate" : "Duration"),
        amount, unit),
      U.el("div", { style: "border-top:1px solid var(--border);margin:12px 0 8px" }),
      resultLine, resultNote,
      U.el("div", { class: "btn-row", style: "justify-content:flex-end;margin-top:12px" },
        cancelBtn, applyBtn),
    );
  }

  return { open };
})();


/* Output variables picker: choose which spatial variables go into the
 * netCDF output and, per variable, which statistics (instantaneous
 * snapshot and/or avg/sum/var/min/max — written as var_stat). */
const OutputVarsPicker = (() => {

  let catalog = null;   // {variables: [{name, dims, desc}], stats: [...]}

  async function open(onDone) {
    if (!catalog) {
      try {
        catalog = await Api.get("/api/schema/output_vars");
      } catch (err) {
        U.toast(err.message, "error");
        return;
      }
    }
    const popup = Popup.open({ title: "Output variables", width: 620 });

    // current selection -> {base: Set("" | stat)}
    const selection = new Map();
    for (const raw of (App.state.config.output_vars || [])) {
      const m = String(raw).match(/^(.*?)(?:[._](avg|sum|var|min|max))?$/);
      const base = m[1], stat = m[2] || "";
      if (!selection.has(base)) selection.set(base, new Set());
      selection.get(base).add(stat);
    }

    const search = U.el("input", { type: "search", placeholder: "Filter variables…", style: "width:100%;margin-bottom:8px" });
    const list = U.el("div", {
      class: "layer-tree", style: "max-height:46vh;overflow-y:auto",
    });

    const renderRow = (v) => {
      const chosen = selection.get(v.name);
      const cb = U.el("input", { type: "checkbox" });
      cb.checked = Boolean(chosen);
      const statBox = U.el("span", { class: "ov-stats" });
      const renderStats = () => {
        U.clear(statBox);
        if (!selection.has(v.name)) return;
        const sel = selection.get(v.name);
        for (const stat of ["", ...catalog.stats]) {
          const chip = U.el("button", {
            class: `year-chip ${sel.has(stat) ? "on" : ""}`,
            style: "padding:1px 7px;font-size:11px",
            title: stat === "" ? "instantaneous snapshot" : `${v.name}_${stat}`,
          }, stat === "" ? "inst" : stat);
          chip.addEventListener("click", () => {
            if (sel.has(stat)) { if (sel.size > 1) sel.delete(stat); }
            else sel.add(stat);
            renderStats();
          });
          statBox.append(chip);
        }
      };
      cb.addEventListener("change", () => {
        if (cb.checked) selection.set(v.name, new Set([""]));
        else selection.delete(v.name);
        renderStats();
      });
      renderStats();
      // fixed grid columns keep names, units/descriptions and the
      // statistic chips vertically aligned across all rows
      const unit = (v.desc || "").match(/^\[[^\]]*\]/);
      const desc = unit ? (v.desc || "").slice(unit[0].length).trim() : (v.desc || "");
      const row = U.el("div", { class: "ov-row", title: v.desc || v.name },
        cb,
        U.el("span", { class: "ov-name" }, v.name),
        U.el("span", { class: "ov-unit" }, unit ? unit[0] : ""),
        U.el("span", { class: "ov-desc" }, desc),
        statBox);
      return row;
    };

    const renderList = () => {
      U.clear(list);
      const term = search.value.trim().toLowerCase();
      const vars = catalog.variables.filter((v) =>
        !term || v.name.toLowerCase().includes(term) || (v.desc || "").toLowerCase().includes(term));
      // selected first, then the rest
      for (const v of vars.filter((x) => selection.has(x.name))) list.append(renderRow(v));
      for (const v of vars.filter((x) => !selection.has(x.name))) list.append(renderRow(v));
      if (!vars.length) list.append(U.el("div", { class: "muted" }, "no matches"));
    };
    search.addEventListener("input", U.debounce(renderList, 150));
    renderList();

    const applyBtn = U.el("button", { class: "primary" }, "Apply");
    applyBtn.addEventListener("click", () => {
      const out = [];
      for (const v of catalog.variables) {
        const sel = selection.get(v.name);
        if (!sel) continue;
        for (const stat of ["", ...catalog.stats]) {
          if (sel.has(stat)) out.push(stat ? `${v.name}_${stat}` : v.name);
        }
      }
      App.state.config.output_vars = out;
      popup.close();
      if (onDone) onDone();
    });

    popup.body.append(
      U.el("div", { class: "muted", style: "font-size:12px;margin-bottom:6px" },
        "Tick a variable, then choose per variable: instantaneous value (inst) ",
        "and/or a statistic over each output interval (avg, sum, var, min, max)."),
      search, list,
      U.el("div", { class: "btn-row", style: "justify-content:flex-end" }, applyBtn),
    );
  }

  return { open };
})();
