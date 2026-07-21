/* Schema-driven configuration form.
 *
 * Renders the sections/parameters served by /api/schema as collapsible
 * groups, with typed inputs, an info tooltip per parameter (description
 * + unit + default from constants.py), changed-from-default
 * highlighting, and cross-tab links for file parameters produced by the
 * Grid / Domain / Conditions tabs.
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

  function build(containerEl, schema, values, onChangeCb) {
    container = containerEl;
    onChange = onChangeCb;
    U.clear(container);

    for (const section of schema.sections) {
      container.append(_section(section, values));
    }
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
    const wrap = U.el("div", { class: "section collapsed", dataset: { section: section.name } }, head, body);
    head.addEventListener("click", () => wrap.classList.toggle("collapsed"));
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
    } else {
      row.append(_input(param, value));
    }

    _markChanged(row, param, value);
    return row;
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

    // numbers, lists, strings, generic -> buffered text input
    const input = U.el("input", { type: "text", id, value: _displayValue(value) });
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
    if (onChange) onChange(param.key, value);
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
      return { ok: true, value: parts };   // list of strings (e.g. output_vars)
    }
    // str / file / any: numbers pass through as numbers for 'any'
    if (param.type === "any") {
      const v = Number(text);
      if (Number.isFinite(v) && /^[\d.eE+-]+$/.test(text)) return { ok: true, value: v };
    }
    return { ok: true, value: text };
  }

  /* ---- search filter ---- */

  function filter(term) {
    term = term.trim().toLowerCase();
    for (const section of container.querySelectorAll(".section")) {
      let any = false;
      for (const row of section.querySelectorAll(".param-row")) {
        const key = row.dataset.param.toLowerCase();
        const hit = !term || key.includes(term);
        row.style.display = hit ? "" : "none";
        if (hit) any = true;
      }
      section.style.display = any ? "" : "none";
      if (term && any) section.classList.remove("collapsed");
      if (!term) section.classList.add("collapsed");
    }
    if (!term) {
      const first = container.querySelector(".section");
      if (first) first.classList.remove("collapsed");
    }
  }

  /* Refresh the displayed value of a single parameter (e.g. after
   * another tab wrote a file parameter). */
  function refreshParam(key) {
    if (!container || !App.state.schema) return;
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
  }

  function _findParam(key) {
    for (const section of App.state.schema.sections) {
      for (const param of section.params) {
        if (param.key === key) return param;
      }
    }
    return null;
  }

  return { build, filter, refreshParam };
})();
