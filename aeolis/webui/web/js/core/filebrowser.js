/* In-app file/folder picker modal.
 *
 * Replaces native OS dialogs, which can open behind the pywebview
 * window and leave the app disabled. Browses through /api/browse.
 *
 *   const path = await FileBrowser.pick({
 *     title, patterns: [["AeoLiS config", "*.txt"]],
 *     mode: "open" | "save" | "folder", initial, filename,
 *   });   // -> path string or null on cancel
 */
"use strict";

const FileBrowser = (() => {

  let modal = null;
  let resolver = null;
  let state = { path: "", mode: "open", patterns: ["*"], filename: "" };

  function pick(options = {}) {
    const patterns = (options.patterns || [["All files", "*.*"]])
      .flatMap(([, glob]) => String(glob).split(";"))
      .map((g) => g.trim())
      .filter(Boolean);
    state = {
      mode: options.save ? "save" : (options.mode || "open"),
      patterns,
      filename: options.filename || "aeolis.txt",
      title: options.title ||
        (options.mode === "folder" ? "Select folder"
          : options.save ? "Save file" : "Open file"),
    };
    _build();
    _navigate(options.initial || localStorage.getItem("aeolis-browse-path") || "");
    return new Promise((resolve) => { resolver = resolve; });
  }

  function _finish(result) {
    if (modal) { modal.remove(); modal = null; }
    document.removeEventListener("keydown", _onKey);
    const resolve = resolver;
    resolver = null;
    if (resolve) resolve(result);
  }

  function _onKey(ev) {
    if (ev.key === "Escape") { ev.stopPropagation(); _finish(null); }
  }

  /* ---- UI ---- */

  let els = {};

  function _build() {
    if (modal) modal.remove();
    els = {};

    els.pathInput = U.el("input", { type: "text", class: "fb-path" });
    els.pathInput.addEventListener("keydown", (ev) => {
      if (ev.key === "Enter") _navigate(els.pathInput.value.trim());
    });
    const upBtn = U.el("button", { class: "ghost", title: "Parent folder" }, "↑");
    upBtn.addEventListener("click", () => { if (state.parent) _navigate(state.parent); });

    const newFolderBtn = U.el("button", { class: "ghost", title: "Create a new folder here" }, "New folder");
    newFolderBtn.addEventListener("click", _promptNewFolder);

    els.drives = U.el("select", { class: "fb-drives", title: "Drive" });
    els.drives.addEventListener("change", () => _navigate(els.drives.value));

    els.list = U.el("div", { class: "fb-list" });

    els.nameInput = U.el("input", { type: "text", class: "fb-name", value: state.filename });

    const okLabel = state.mode === "folder" ? "Select this folder"
      : state.mode === "save" ? "Save here" : "Open";
    const okBtn = U.el("button", { class: "primary" }, okLabel);
    okBtn.addEventListener("click", () => {
      if (state.mode === "folder") _finish(state.path);
      else if (state.mode === "save") {
        const name = els.nameInput.value.trim();
        if (name) _finish(state.path + state.sep + name);
      } else if (state.selected) {
        _finish(state.path + state.sep + state.selected);
      } else {
        U.toast("Select a file first", "error");
      }
    });
    const cancelBtn = U.el("button", { class: "ghost" }, "Cancel");
    cancelBtn.addEventListener("click", () => _finish(null));

    const nameRow = state.mode === "save"
      ? U.el("div", { class: "fb-name-row" }, U.el("label", {}, "File name"), els.nameInput)
      : null;

    const box = U.el("div", { class: "modal fb-modal" },
      U.el("h2", {}, state.title),
      U.el("div", { class: "fb-toolbar" }, upBtn, newFolderBtn, els.drives, els.pathInput),
      els.list,
      nameRow,
      U.el("div", { class: "modal-actions" }, okBtn, cancelBtn),
    );
    modal = U.el("div", { class: "fb-backdrop" }, box);
    modal.addEventListener("mousedown", (ev) => {
      if (ev.target === modal) _finish(null);
    });
    document.body.append(modal);
    document.addEventListener("keydown", _onKey);
  }

  /* ---- navigation ---- */

  async function _navigate(path) {
    let res;
    try {
      res = await Api.post("/api/browse", { path, patterns: state.patterns });
    } catch (err) {
      U.toast(err.message, "error");
      return;
    }
    state.path = res.path;
    state.parent = res.parent;
    state.sep = res.sep;
    state.selected = null;
    localStorage.setItem("aeolis-browse-path", res.path);

    els.pathInput.value = res.path;
    U.clear(els.drives);
    for (const drive of res.drives) {
      els.drives.append(U.el("option", {
        value: drive,
        selected: res.path.toLowerCase().startsWith(drive.toLowerCase()) ? "" : null,
      }, drive));
    }

    U.clear(els.list);
    for (const dir of res.dirs) {
      const row = U.el("div", { class: "fb-row" },
        U.el("span", { class: "fb-icon" }, "📁"),
        U.el("span", { class: "fb-rowname" }, dir));
      row.addEventListener("dblclick", () => _navigate(res.path + res.sep + dir));
      row.addEventListener("click", () => _highlight(row));
      els.list.append(row);
    }
    for (const file of res.files) {
      const row = U.el("div", { class: "fb-row" },
        U.el("span", { class: "fb-icon" }, "📄"),
        U.el("span", { class: "fb-rowname" }, file.name),
        U.el("span", { class: "fb-meta" }, U.fmtBytes(file.size)));
      row.addEventListener("click", () => {
        _highlight(row);
        state.selected = file.name;
        if (state.mode === "save") els.nameInput.value = file.name;
      });
      row.addEventListener("dblclick", () => {
        if (state.mode !== "folder") _finish(res.path + res.sep + file.name);
      });
      els.list.append(row);
    }
    if (!res.dirs.length && !res.files.length) {
      els.list.append(U.el("div", { class: "muted", style: "padding:10px" }, "Empty folder"));
    }
  }

  /* Inline "new folder" prompt: an editable row at the top of the list.
   * Enter creates the folder (via /api/mkdir) and navigates into it so
   * the user can immediately save there; Escape cancels just the prompt. */
  function _promptNewFolder() {
    if (!state.path) return;
    const existing = els.list.querySelector(".fb-newfolder input");
    if (existing) { existing.focus(); existing.select(); return; }

    const input = U.el("input", { type: "text", class: "fb-newfolder-input", value: "New folder" });
    const row = U.el("div", { class: "fb-row fb-newfolder" },
      U.el("span", { class: "fb-icon" }, "📁"), input);
    els.list.prepend(row);
    input.focus();
    input.select();

    let done = false;
    const commit = async () => {
      if (done) return;
      const name = input.value.trim();
      if (!name) { row.remove(); return; }
      done = true;
      try {
        const res = await Api.post("/api/mkdir", { path: state.path, name });
        _navigate(res.path);      // enter the freshly-made folder
      } catch (err) {
        done = false;
        U.toast(err.message, "error");
        input.focus();
        input.select();
      }
    };
    input.addEventListener("keydown", (ev) => {
      if (ev.key === "Enter") { ev.preventDefault(); commit(); }
      else if (ev.key === "Escape") { ev.stopPropagation(); done = true; row.remove(); }
    });
    input.addEventListener("blur", () => { if (!done) row.remove(); });
  }

  function _highlight(row) {
    for (const r of els.list.querySelectorAll(".fb-row")) r.classList.remove("selected");
    row.classList.add("selected");
    if (!row.querySelector(".fb-meta")) state.selected = null;   // a directory
  }

  return { pick };
})();
