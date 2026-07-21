"""Native OS file/folder dialogs, served to the frontend.

Uses tkinter (stdlib) on the server side, so the pick happens in a real
OS dialog on the machine running the GUI — same pattern as the
SedTRAILS GUI. Falls back with an error the frontend can handle by
letting the user type a path.
"""

import threading

from aeolis.webui.backend.httpd import route
from aeolis.webui.backend.util import send_error_json, send_json

_TK_LOCK = threading.Lock()


def _with_tk(fn):
    """Run a tkinter dialog function with a temporary hidden root."""
    import tkinter as tk

    with _TK_LOCK:
        root = tk.Tk()
        root.withdraw()
        root.attributes("-topmost", True)
        try:
            return fn()
        finally:
            root.destroy()


@route("POST", "/api/pickfile")
def _pickfile(handler, body, tail):
    title = body.get("title") or "Select file"
    patterns = body.get("patterns") or [["All files", "*.*"]]
    save = bool(body.get("save"))
    initial = body.get("initial") or ""
    try:
        from tkinter import filedialog

        def _ask():
            kwargs = {
                "title": title,
                "filetypes": [tuple(p) for p in patterns],
            }
            if initial:
                kwargs["initialdir"] = initial
            if save:
                return filedialog.asksaveasfilename(**kwargs)
            return filedialog.askopenfilename(**kwargs)

        path = _with_tk(_ask)
    except Exception as exc:  # noqa: BLE001 - tkinter may be unavailable
        send_error_json(handler, f"native dialog unavailable: {exc}", 501)
        return
    send_json(handler, {"path": path or None})


@route("POST", "/api/pickfolder")
def _pickfolder(handler, body, tail):
    title = body.get("title") or "Select folder"
    initial = body.get("initial") or ""
    try:
        from tkinter import filedialog

        def _ask():
            kwargs = {"title": title}
            if initial:
                kwargs["initialdir"] = initial
            return filedialog.askdirectory(**kwargs)

        path = _with_tk(_ask)
    except Exception as exc:  # noqa: BLE001 - tkinter may be unavailable
        send_error_json(handler, f"native dialog unavailable: {exc}", 501)
        return
    send_json(handler, {"path": path or None})
