"""Run tab API routes."""

from aeolis.webui.backend import project, run_manager
from aeolis.webui.backend.config_api import load_config
from aeolis.webui.backend.httpd import route
from aeolis.webui.backend.util import send_error_json, send_json


@route("GET", "/api/run/backends")
def _backends(handler, query, tail):
    send_json(handler, {
        "active": run_manager.active.id,
        "backends": [
            {"id": b.id, "title": b.title, "available": b.available, "note": b.note}
            for b in run_manager.BACKENDS.values()
        ],
    })


@route("POST", "/api/run/backend")
def _select_backend(handler, body, tail):
    try:
        backend = run_manager.select(body.get("id"))
    except (ValueError, RuntimeError) as exc:
        send_error_json(handler, exc)
        return
    send_json(handler, {"ok": True, "active": backend.id})


@route("GET", "/api/run/checklist")
def _checklist(handler, query, tail):
    current = project.require()
    values = load_config(current.configfile)
    checks = []

    def _file_check(key, required):
        filename = values.get(key)
        if not filename:
            checks.append({"key": key, "ok": not required,
                           "text": f"{key} not set" + (" (required)" if required else " (optional)")})
        else:
            exists = (current.root / str(filename)).is_file()
            checks.append({"key": key, "ok": exists,
                           "text": f"{key} = {filename}" + ("" if exists else " (file missing!)")})

    for key in ("xgrid_file", "ygrid_file", "bed_file", "wind_file"):
        _file_check(key, required=True)
    for key in ("ne_file", "veg_file", "tide_file", "wave_file", "meteo_file"):
        _file_check(key, required=False)

    duration = (values.get("tstop") or 0) - (values.get("tstart") or 0)
    checks.append({"key": "time", "ok": duration > 0,
                   "text": f"simulation duration {duration:g} s"})
    send_json(handler, {"checks": checks, "ready": all(
        c["ok"] for c in checks if "(optional)" not in c["text"])})


@route("POST", "/api/run/start")
def _start(handler, body, tail):
    current = project.require()
    try:
        run_manager.active.start(current)
    except RuntimeError as exc:
        send_error_json(handler, exc)
        return
    send_json(handler, {"ok": True})


@route("POST", "/api/run/stop")
def _stop(handler, body, tail):
    stopped = run_manager.active.stop()
    send_json(handler, {"ok": True, "stopped": stopped})


@route("GET", "/api/run/status")
def _status(handler, query, tail):
    send_json(handler, run_manager.active.status())


@route("GET", "/api/run/log")
def _log(handler, query, tail):
    offset = int(query.get("offset", 0))
    send_json(handler, run_manager.active.log_tail(offset))
