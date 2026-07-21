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
      body.append(_paramRow(param, values));
    }
    const head = U.el("header", {},
      U.el("span", { class: "caret" }, "▾"),
      section.name,
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

  function _paramRow(param, values) {
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

    if (param.link) {
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

  function _inputWrap(param, value) {
    const input = _input(param, value);
    const showUnit = param.unit && param.unit !== "-" && param.type !== "bool";
    const extras = [];
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
    if (row) _markChanged(row, param, value);
    applyRules();
    if (onChange) onChange(param.key, value);
    App.emit("config-changed", param.key);
  }

  function _markChanged(row, param, value) {
    row.classList.toggle("changed", !_equals(value, param.default));
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
      for (const param of section.params) {
        const row = sectionEl.querySelector(`[data-param="${param.key}"]`);
        if (!row) continue;
        const ruleOk = (sectionVisible || param.key === controller)
          && _ruleSatisfied(param.visible_if);
        const searchOk = !searchTerm || param.key.toLowerCase().includes(searchTerm);
        const show = ruleOk && searchOk;
        row.style.display = show ? "" : "none";
        if (show) anyVisible = true;
      }
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
    } else {
      const input = row.querySelector("input, select");
      if (input && input.type !== "checkbox") input.value = _displayValue(value);
      else if (input) input.checked = Boolean(value);
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
