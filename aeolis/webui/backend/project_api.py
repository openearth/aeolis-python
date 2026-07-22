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
