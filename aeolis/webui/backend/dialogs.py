"""Native OS file/folder dialogs, served to the frontend.

When the GUI runs in its pywebview window, pywebview's own
create_file_dialog is used — it is parented to the app window and safe
to call from the server thread. In browser mode (no pywebview window)
a tkinter dialog is the fallback. If neither works the frontend lets
the user type a path.
"""

import threading

from aeolis.webui.backend.httpd import route
from aeolis.webui.backend.util import send_error_json, send_json

_TK_LOCK = threading.Lock()


def _webview_window():
    try:
        import webview
        return webview.windows[0] if webview.windows else None
    except ImportError:
        return None


def _via_webview(save, folder, title, patterns, initial, filename="aeolis.txt"):
    import webview

    window = _webview_window()
    if window is None:
        return None, False

    if folder:
        result = window.create_file_dialog(webview.FOLDER_DIALOG, directory=initial or "")
    elif save:
        result = window.create_file_dialog(
            webview.SAVE_DIALOG, directory=initial or "", save_filename=filename)
    else:
        file_types = tuple(
            f"{name} ({ext})" for name, ext in patterns
        ) or ("All files (*.*)",)
        result = window.create_file_dialog(
            webview.OPEN_DIALOG, directory=initial or "", file_types=file_types)

    if not result:
        return None, True
    path = result[0] if isinstance(result, (list, tuple)) else result
    return str(path), True


def _via_tkinter(save, folder, title, patterns, initial):
    import tkinter as tk
    from tkinter import filedialog

    with _TK_LOCK:
        root = tk.Tk()
        root.withdraw()
        root.attributes("-topmost", True)
        try:
            kwargs = {"title": title}
            if initial:
                kwargs["initialdir"] = initial
            if folder:
                path = filedialog.askdirectory(**kwargs)
            else:
                kwargs["filetypes"] = [tuple(p) for p in patterns]
                if save:
                    path = filedialog.asksaveasfilename(**kwargs)
                else:
                    path = filedialog.askopenfilename(**kwargs)
        finally:
            root.destroy()
    return path or None


def _pick(handler, body, folder=False):
    title = body.get("title") or ("Select folder" if folder else "Select file")
    patterns = body.get("patterns") or [["All files", "*.*"]]
    save = bool(body.get("save"))
    initial = body.get("initial") or ""

    try:
        path, handled = _via_webview(save, folder, title, patterns, initial,
                                     filename=body.get("filename") or "aeolis.txt")
        if handled:
            send_json(handler, {"path": path})
            return
    except Exception:  # noqa: BLE001 - fall through to tkinter
        pass

    try:
        path = _via_tkinter(save, folder, title, patterns, initial)
    except Exception as exc:  # noqa: BLE001 - tkinter may be unavailable
        send_error_json(handler, f"native dialog unavailable: {exc}", 501)
        return
    send_json(handler, {"path": path})


@route("POST", "/api/pickfile")
def _pickfile(handler, body, tail):
    _pick(handler, body, folder=False)


@route("POST", "/api/pickfolder")
def _pickfolder(handler, body, tail):
    _pick(handler, body, folder=True)


# ---------------------------------------------------------------------
# In-app file browser (the default picker).
#
# Native modal dialogs opened from the server thread can appear BEHIND
# the pywebview window on Windows, leaving the app greyed-out and
# unclickable. The frontend therefore browses the filesystem through
# this endpoint and renders its own picker modal.
# ---------------------------------------------------------------------

import fnmatch
import os
from pathlib import Path


def _list_drives():
    if os.name != "nt":
        return ["/"]
    drives = []
    for letter in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
        if Path(f"{letter}:\\").exists():
            drives.append(f"{letter}:\\")
    return drives


@route("POST", "/api/browse")
def _browse(handler, body, tail):
    raw = body.get("path") or str(Path.home())
    patterns = body.get("patterns") or ["*"]
    path = Path(raw).expanduser()
    if path.is_file():
        path = path.parent
    if not path.is_dir():
        path = Path.home()
    path = path.resolve()

    dirs = []
    files = []
    try:
        for entry in sorted(path.iterdir(), key=lambda e: e.name.lower()):
            name = entry.name
            if name.startswith((".", "$")) or name.lower() in ("system volume information",):
                continue
            try:
                if entry.is_dir():
                    dirs.append(name)
                elif any(fnmatch.fnmatch(name.lower(), p.lower()) for p in patterns):
                    stat = entry.stat()
                    files.append({
                        "name": name,
                        "size": stat.st_size,
                        "mtime": int(stat.st_mtime),
                    })
            except OSError:
                continue
    except OSError as exc:
        send_error_json(handler, f"cannot list {path}: {exc}", 403)
        return

    parent = str(path.parent) if path.parent != path else None
    send_json(handler, {
        "path": str(path),
        "parent": parent,
        "sep": os.sep,
        "dirs": dirs,
        "files": files,
        "drives": _list_drives(),
    })


@route("POST", "/api/mkdir")
def _mkdir(handler, body, tail):
    """Create a single subfolder inside an existing directory, so the
    save-as picker can make a new folder to save into."""
    parent = body.get("path")
    name = (body.get("name") or "").strip()
    if not parent or not name:
        send_error_json(handler, "missing 'path' or 'name'")
        return
    # single path segment only - never traverse
    if name in (".", "..") or any(sep in name for sep in ("/", "\\")):
        send_error_json(handler, "invalid folder name")
        return
    base = Path(parent).expanduser()
    if not base.is_dir():
        send_error_json(handler, f"not a folder: {parent}", 404)
        return
    target = (base / name).resolve()
    try:
        target.mkdir(parents=False, exist_ok=False)
    except FileExistsError:
        send_error_json(handler, f"'{name}' already exists")
        return
    except OSError as exc:
        send_error_json(handler, f"cannot create folder: {exc}", 403)
        return
    send_json(handler, {"path": str(target)})
