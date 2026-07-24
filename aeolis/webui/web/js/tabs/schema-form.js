/* Schema-driven configuration form.
 *
 * Renders the sections/parameters served by /api/schema as collapsible
 * groups with typed inputs, permanent unit suffixes, info tooltips,
 * changed-from-default highlighting, cross-tab links for file
 * parameters, docs links per section, a date/duration helper for time
 * parameters, and live conditional visibility (e.g. transport constants
 * follow method_transport; a section whose parameters are all hidden is
 * hidden entirely).
 */
"use strict";

const SchemaForm = (() => {

  const TAB_LABELS = {
    grid: "Grid tab",
    domain: "Domain tab",
    conditions: "Conditions tab",
  };

  let container = null;
  let onChange = null;
  let schema = null;
  let searchTerm = "";
  let divider = null;          // "Disabled" divider element
  let lastDisabledKey = "";    // change detection for DOM reordering

  function build(containerEl, schemaObj, values, onChangeCb) {
    container = containerEl;
    schema = schemaObj;
    onChange = onChangeCb;
    lastDisabledKey = "";
    U.clear(container);

    for (const section of schema.sections) {
      container.append(_section(section, values));
    }
    divider = U.el("div", { class: "disabled-divider", hidden: "" }, "Disabled");
    container.append(divider);
    applyRules();
  }

  function _section(section, values) {
    const body = U.el("div", { class: "section-body" });
    for (const param of section.params) {
      body.append(_paramRow(param, values, section));
    }
    // sections whose files come from a dedicated tab link there from
    // the header (the rows themselves just show the configured path)
    let tabLink = null;
    if (section.tab) {
      tabLink = U.el("button", { class: "section-tablink", title: `Open the ${TAB_LABELS[section.tab]}` },
        `${TAB_LABELS[section.tab]} →`);
      tabLink.addEventListener("click", (ev) => {
        ev.stopPropagation();
        Tabs.activate(section.tab);
      });
    }
    const head = U.el("header", {},
      U.el("span", { class: "caret" }, "▾"),
      section.name,
      tabLink,
      U.el("span", { class: "count" }, String(section.params.length)),
    );
    const wrap = U.el("div", {
      class: "section collapsed",
      dataset: { section: section.name },
    }, head, body);
    head.addEventListener("click", () => {
      wrap.classList.toggle("collapsed");
    });
    return wrap;
  }

  function _paramRow(param, values, section = null) {
    const value = values[param.key];
    const row = U.el("div", { class: "param-row", dataset: { param: param.key } });

    const info = U.el("span", { class: "info-dot" }, "i");
    info.addEventListener("mouseenter", (ev) => {
      const unit = param.unit ? `<span class="tt-unit">[${param.unit}]</span> ` : "";
      const desc = param.desc || "(no description)";
      const def = `<div class="tt-default">default: ${_displayValue(param.default)}</div>`;
      U.showTip(`<b>${param.key}</b><br>${unit}${desc}${def}`, ev.clientX, ev.clientY);
    });
    info.addEventListener("mouseleave", U.hideTip);

    const label = U.el("label", { for: `param-${param.key}` }, info, param.key);
    row.append(label);

    if (param.link && section && section.tab) {
      // the section header links to the tab; the row just shows the
      // configured path, read-only
      const pathBox = U.el("input", {
        type: "text", class: "readonly filepath", disabled: "",
        value: _displayValue(value) || "",
        placeholder: "not set",
        title: `Set in the ${TAB_LABELS[param.link]}`,
      });
      row.append(pathBox);
    } else if (param.link) {
      const link = U.el("button", { class: "filelink", title: `Configured in the ${TAB_LABELS[param.link]}` },
        `${_displayValue(value) || "not set"} → ${TAB_LABELS[param.link]}`);
      link.addEventListener("click", () => Tabs.activate(param.link));
      row.append(link);
    } else if (param.picker === "output_vars") {
      const btn = U.el("button", { class: "filelink", title: "Select output variables and statistics" });
      const syncLabel = () => {
        const v = App.state.config ? App.state.config[param.key] : value;
        const n = Array.isArray(v) ? v.length : 0;
        btn.textContent = n ? `${n} variable${n === 1 ? "" : "s"} — click to edit` : "select variables…";
      };
      syncLabel();
      btn.addEventListener("click", () => OutputVarsPicker.open(() => {
        syncLabel();
        _commit(param, App.state.config[param.key]);
      }));
      row.append(btn);
    } else {
      row.append(_inputWrap(param, value));
    }

    _markChanged(row, param, value);
    return row;
  }

  /* Human word for the entries of a list parameter. */
  function _listNoun(key) {
    if (/^grain_/.test(key)) return "fraction";
    if (/species|veg|tiller|Hveg|Nt/i.test(key)) return "species";
    return "value";
  }

  /* Small gray badge showing how many entries a list value holds. */
  function _listBadge(param) {
    const badge = U.el("span", {
      class: "list-badge", dataset: { listBadge: param.key },
      title: "This parameter accepts multiple space-separated values",
    });
    _syncListBadge(badge, param);
    return badge;
  }

  function _syncListBadge(badge, param) {
    const v = App.state.config ? App.state.config[param.key] : param.default;
    const n = Array.isArray(v) ? v.length : (v === null || v === undefined ? 0 : 1);
    const noun = _listNoun(param.key);
    badge.textContent = n > 1 ? `${n} ${noun}s` : (n === 1 ? `1 ${noun}` : "");
    badge.style.display = badge.textContent ? "" : "none";
  }

  function _inputWrap(param, value) {
    const input = _input(param, value);
    const showUnit = param.unit && param.unit !== "-" && param.type !== "bool";
    const extras = [];
    if (param.type === "list") extras.push(_listBadge(param));
    if (showUnit) extras.push(U.el("span", { class: "unit-suffix" }, param.unit));
    if (param.time_tool) {
      const clock = U.el("button", { class: "ghost time-tool-btn", title: "Compute from date/duration" }, "🕒");
      clock.addEventListener("click", () => {
        TimeTool.open(param.key, _refdateEpoch(), (seconds) => {
          _commit(param, seconds);
          const inp = container.querySelector(`[data-param="${param.key}"] input`);
          if (inp) inp.value = _displayValue(seconds);
        });
      });
      extras.push(clock);
    }
    if (!extras.length) return input;
    return U.el("span", { class: "input-wrap" }, input, ...extras);
  }

  function _refdateEpoch() {
    const raw = (App.state.config && App.state.config.refdate) || "2020-01-01 00:00";
    const parsed = Date.parse(raw.replace(" ", "T") + (raw.length <= 16 ? ":00Z" : "Z"));
    return Number.isFinite(parsed) ? parsed / 1000 : 0;
  }

  function _input(param, value) {
    const id = `param-${param.key}`;

    if (param.type === "bool") {
      const input = U.el("input", { type: "checkbox", id });
      input.checked = Boolean(value);
      input.addEventListener("change", () => _commit(param, input.checked));
      return input;
    }

    if (param.type === "datetime") {
      const input = U.el("input", { type: "datetime-local", id });
      input.value = _toDatetimeLocal(value);
      input.addEventListener("change", () => _commit(param, _fromDatetimeLocal(input.value)));
      return input;
    }

    if (param.options) {
      const select = U.el("select", { id });
      const opts = [...param.options];
      if (value !== null && value !== undefined && !opts.includes(value)) opts.unshift(value);
      for (const opt of opts) {
        select.append(U.el("option", { value: opt, selected: opt === value ? "" : null }, opt));
      }
      select.addEventListener("change", () => _commit(param, select.value));
      return select;
    }

    const input = U.el("input", { type: "text", id, value: _displayValue(value) });
    if (param.readonly) {
      input.disabled = true;
      input.classList.add("readonly");
      input.title = "Derived from the grid files (set in the Grid tab)";
      return input;
    }
    const commit = () => {
      const parsed = _parseValue(param, input.value);
      if (parsed.ok) {
        _commit(param, parsed.value);
        input.value = _displayValue(parsed.value);
        input.style.borderColor = "";
      } else {
        input.style.borderColor = "var(--danger)";
      }
    };
    input.addEventListener("blur", commit);
    input.addEventListener("keydown", (ev) => {
      if (ev.key === "Enter") input.blur();
      if (ev.key === "Escape") { input.value = _displayValue(value); input.blur(); }
    });
    return input;
  }

  function _commit(param, value) {
    App.state.config[param.key] = value;
    const row = container.querySelector(`[data-param="${param.key}"]`);
    if (row) {
      _markChanged(row, param, value);
      const badge = row.querySelector(`[data-list-badge="${CSS.escape(param.key)}"]`);
      if (badge) _syncListBadge(badge, param);
    }
    applyRules();
    if (onChange) onChange(param.key, value);
    App.emit("config-changed", param.key);
  }

  function _markChanged(row, param, value) {
    const changed = !_equals(value, param.default);
    row.classList.toggle("changed", changed);
    // explain the teal dot rendered by .param-row.changed::before
    if (changed) {
      row.title = `differs from the default (${_displayValue(param.default) || "none"})`;
    } else {
      row.removeAttribute("title");
    }
  }

  function _equals(a, b) {
    if (Array.isArray(a) && Array.isArray(b)) {
      return a.length === b.length && a.every((v, i) => _equals(v, b[i]));
    }
    return a === b || (a === null && b === null);
  }

  function _displayValue(value) {
    if (value === null || value === undefined) return "";
    if (Array.isArray(value)) return value.join(" ");
    return String(value);
  }

  /* AeoLiS stores dates as "YYYY-MM-DD HH:MM"; <input datetime-local> wants a "T". */
  function _toDatetimeLocal(value) {
    if (!value) return "";
    const s = String(value).trim().replace(" ", "T");
    return s.slice(0, 16); // drop any seconds
  }
  function _fromDatetimeLocal(value) {
    if (!value) return null;
    return value.replace("T", " ").slice(0, 16);
  }

  function _parseValue(param, text) {
    text = text.trim();
    if (text === "" || text.toLowerCase() === "none") return { ok: true, value: null };

    if (param.type === "int") {
      const v = Number(text);
      return Number.isInteger(v) ? { ok: true, value: v } : { ok: false };
    }
    if (param.type === "float") {
      const v = Number(text);
      return Number.isFinite(v) ? { ok: true, value: v } : { ok: false };
    }
    if (param.type === "list") {
      const parts = text.split(/[\s,]+/).filter(Boolean);
      const nums = parts.map(Number);
      if (nums.every(Number.isFinite)) return { ok: true, value: nums };
      return { ok: true, value: parts };
    }
    if (param.type === "any") {
      const v = Number(text);
      if (Number.isFinite(v) && /^[\d.eE+-]+$/.test(text)) return { ok: true, value: v };
    }
    return { ok: true, value: text };
  }

  /* ================= conditional visibility ================= */

  function _ruleSatisfied(rule) {
    if (!rule) return true;
    if (rule.any) return rule.any.some(_ruleSatisfied);
    if (rule.all) return rule.all.every(_ruleSatisfied);
    const current = App.state.config ? App.state.config[rule.key] : undefined;
    return (rule.in || []).some((v) => v === current);
  }

  function applyRules() {
    if (!container || !schema) return;
    const disabledNames = [];
    for (const section of schema.sections) {
      const sectionEl = container.querySelector(`[data-section="${CSS.escape(section.name)}"]`);
      if (!sectionEl) continue;
      const sectionVisible = _ruleSatisfied(section.visible_if);
      const enabled = _ruleSatisfied(section.enabled_if);
      // the parameter controlling a section's visibility must never be
      // hidden by it
      const controller = section.visible_if ? section.visible_if.key : null;
      let anyVisible = false;
      let visibleCount = 0;
      for (const param of section.params) {
        const row = sectionEl.querySelector(`[data-param="${param.key}"]`);
        if (!row) continue;
        const ruleOk = (sectionVisible || param.key === controller)
          && _ruleSatisfied(param.visible_if);
        const searchOk = !searchTerm || param.key.toLowerCase().includes(searchTerm);
        const show = ruleOk && searchOk;
        row.style.display = show ? "" : "none";
        if (show) { anyVisible = true; visibleCount += 1; }
      }
      // the header count follows the currently applicable parameters
      // (e.g. switching vegetation method changes it)
      const countEl = sectionEl.querySelector(":scope > header .count");
      if (countEl) countEl.textContent = String(visibleCount);
      sectionEl.style.display = anyVisible ? "" : "none";
      sectionEl.classList.toggle("disabled-group", !enabled);
      if (!enabled) disabledNames.push(section.name);
      if (searchTerm && anyVisible) sectionEl.classList.remove("collapsed");
    }
    _parkDisabled(disabledNames);
  }

  /* Move sections whose process is switched off below the "Disabled"
   * divider (keeping schema order in both zones). */
  function _parkDisabled(disabledNames) {
    if (!divider) return;
    const key = disabledNames.join("|");
    if (key === lastDisabledKey) return;
    lastDisabledKey = key;
    const disabled = new Set(disabledNames);
    for (const section of schema.sections) {
      const el = container.querySelector(`[data-section="${CSS.escape(section.name)}"]`);
      if (el && !disabled.has(section.name)) container.insertBefore(el, divider);
    }
    divider.hidden = disabledNames.length === 0;
    for (const section of schema.sections) {
      const el = container.querySelector(`[data-section="${CSS.escape(section.name)}"]`);
      if (el && disabled.has(section.name)) {
        el.classList.add("collapsed");
        container.append(el);
      }
    }
  }

  /* ---- search filter ---- */

  function filter(term) {
    searchTerm = term.trim().toLowerCase();
    applyRules();
    if (!searchTerm) {
      for (const sectionEl of container.querySelectorAll(".section")) {
        sectionEl.classList.add("collapsed");
      }
      const first = container.querySelector(".section");
      if (first) first.classList.remove("collapsed");
    }
  }

  /* Refresh the displayed value of a single parameter. */
  function refreshParam(key) {
    if (!container || !schema) return;
    const row = container.querySelector(`[data-param="${key}"]`);
    if (!row) return;
    const param = _findParam(key);
    const value = App.state.config[key];
    if (param && param.link) {
      const link = row.querySelector(".filelink");
      if (link) link.textContent = `${_displayValue(value) || "not set"} → ${TAB_LABELS[param.link]}`;
      const pathBox = row.querySelector("input.filepath");
      if (pathBox) pathBox.value = _displayValue(value) || "";
    } else {
      const input = row.querySelector("input, select");
      if (input && input.type !== "checkbox") input.value = _displayValue(value);
      else if (input) input.checked = Boolean(value);
      const badge = row.querySelector(`[data-list-badge="${CSS.escape(key)}"]`);
      if (badge && param) _syncListBadge(badge, param);
    }
    if (param) _markChanged(row, param, value);
    applyRules();
  }

  function _findParam(key) {
    for (const section of schema.sections) {
      for (const param of section.params) {
        if (param.key === key) return param;
      }
    }
    return null;
  }

  return { build, filter, refreshParam, applyRules };
})();
