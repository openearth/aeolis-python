"""Project open/create/state API routes."""

from aeolis.webui.backend import project
from aeolis.webui.backend.httpd import route
from aeolis.webui.backend.util import send_error_json, send_json


@route("GET", "/api/project")
def _info(handler, query, tail):
    current = project.current()
    if current is None:
        send_json(handler, {"open": False, "recent": project.recent_projects()})
    else:
        info = current.info()
        info["open"] = True
        info["recent"] = project.recent_projects()
        send_json(handler, info)


@route("POST", "/api/project/open")
def _open(handler, body, tail):
    path = body.get("path")
    if not path:
        send_error_json(handler, "missing 'path'")
        return
    opened = project.open_project(path)
    if not opened.configfile.is_file():
        send_error_json(handler, f"no aeolis.txt found at {opened.configfile}", 404)
        return
    send_json(handler, opened.info())


@route("POST", "/api/project/new")
def _new(handler, body, tail):
    folder = body.get("folder")
    if not folder:
        send_error_json(handler, "missing 'folder'")
        return
    created = project.new_project(folder)
    send_json(handler, created.info())


def _referenced_inputs(values, root):
    """Yield (key, raw_value, abs_path, is_internal, exists) for every config
    key that holds an input-file path. ``is_internal`` is True when the file
    lives under the project *root* (so a plain copytree already brings it
    along); external refs (``../…`` or absolute paths elsewhere) are what break
    when a project is duplicated."""
    import os
    from pathlib import Path

    from aeolis.webui.backend.schema_api import LINKS

    root = Path(os.path.normpath(str(root)))
    for key in LINKS:
        raw = values.get(key)
        if not isinstance(raw, str) or not raw.strip():
            continue
        raw = raw.strip()
        p = Path(raw)
        abs_path = p if p.is_absolute() else root / p
        abs_path = Path(os.path.normpath(str(abs_path)))
        try:
            abs_path.relative_to(root)
            is_internal = True
        except ValueError:
            is_internal = False
        yield key, raw, abs_path, is_internal, abs_path.is_file()


@route("POST", "/api/project/duplicate")
def _duplicate(handler, body, tail):
    """Clone the current model into a fresh folder the user picks, so it can
    be a starting point for a new run without touching existing files. Model
    inputs (config, .grd, .txt, GUI state, raw data) are always copied; run
    outputs (aeolis.nc, *.log) come along only when ``include_outputs`` is
    set. The gui/cache scratch dir is never copied.

    Input files referenced from *outside* the project root (e.g. a shared
    ``../input/z.grd``) are handled per ``input_mode``:

    - ``gather`` (default): copy each external file into the new folder and
      rewrite aeolis.txt to the local filename, so the copy is self-contained.
    - ``keep``: leave the originals in place and rewrite aeolis.txt to their
      absolute paths so the links still resolve from the new location."""
    import os
    import shutil
    from pathlib import Path

    from aeolis.webui.backend.config_api import load_config

    current = project.require()
    parent = body.get("parent")
    name = (body.get("name") or "").strip()
    include_outputs = bool(body.get("include_outputs"))
    input_mode = body.get("input_mode") or "gather"
    if input_mode not in ("gather", "keep"):
        input_mode = "gather"
    if not parent or not name:
        send_error_json(handler, "missing destination folder or name")
        return
    dest = (Path(parent).expanduser() / name).resolve()
    if dest.exists():
        send_error_json(handler, f"'{dest}' already exists - choose another name", 409)
        return

    # Snapshot the source config (file-path values as strings) before copying.
    src_values = load_config(current.configfile)
    src_root = current.root
    cache_dir = current.cache_dir.resolve()

    def _ignore(dirpath, names):
        skip = set()
        here = Path(dirpath).resolve()
        for n in names:
            low = n.lower()
            if not include_outputs and low.endswith((".nc", ".log")):
                skip.add(n)
            elif (here / n).resolve() == cache_dir:
                skip.add(n)
        return skip

    try:
        shutil.copytree(current.root, dest, ignore=_ignore)
    except Exception as exc:  # noqa: BLE001 - report copy failures to the UI
        send_error_json(handler, f"copy failed: {exc}")
        return

    # Relink input files. Internal refs stored as absolute paths get
    # relativized so the copy is portable; external refs are gathered or
    # kept per input_mode. Missing sources are reported, not fatal.
    patch = {}
    skipped = []
    used_names = {}   # abs source path -> local filename already assigned
    dest = dest.resolve()
    for key, raw, abs_path, is_internal, exists in _referenced_inputs(src_values, src_root):
        if is_internal:
            if os.path.isabs(raw):
                rel = os.path.relpath(str(abs_path), str(src_root)).replace("\\", "/")
                patch[key] = rel
            continue
        if not exists:
            skipped.append(raw)
            continue
        if input_mode == "keep":
            patch[key] = str(abs_path)
            continue
        # gather: copy the external file into the new folder (unique basename)
        src_key = str(abs_path)
        if src_key in used_names:
            patch[key] = used_names[src_key]
            continue
        base = abs_path.name
        target = dest / base
        stem, suffix = os.path.splitext(base)
        n = 2
        while target.exists() and target.resolve() != abs_path:
            base = f"{stem}_{n}{suffix}"
            target = dest / base
            n += 1
        try:
            shutil.copy2(abs_path, target)
        except Exception as exc:  # noqa: BLE001
            skipped.append(f"{raw} ({exc})")
            continue
        used_names[src_key] = base
        patch[key] = base

    opened = project.open_project(dest / "aeolis.txt")
    if patch:
        # patch_config targets the now-current (new) project's aeolis.txt
        from aeolis.webui.backend.grid_api import patch_config
        patch_config(patch)
    info = opened.info()
    info["open"] = True
    info["recent"] = project.recent_projects()
    info["skipped"] = skipped
    send_json(handler, info)


@route("POST", "/api/project/reveal")
def _reveal(handler, body, tail):
    """Open the OS file explorer at the config file location."""
    import subprocess
    import sys

    current = project.require()
    target = current.configfile if current.configfile.is_file() else current.root
    if sys.platform == "win32":
        # explorer /select highlights the file inside its folder
        subprocess.Popen(["explorer", "/select,", str(target)])
    elif sys.platform == "darwin":
        subprocess.Popen(["open", "-R", str(target)])
    else:
        subprocess.Popen(["xdg-open", str(target.parent)])
    send_json(handler, {"ok": True})


@route("GET", "/api/project/state")
def _load_state(handler, query, tail):
    send_json(handler, project.require().load_state())


@route("POST", "/api/project/state")
def _save_state(handler, body, tail):
    project.require().save_state(body)
    send_json(handler, {"ok": True})
