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
    from aeolis.webui.backend.validation import run_checks
    checks = run_checks(project.require())
    send_json(handler, {
        "checks": checks,
        "ready": all(c["level"] != "error" for c in checks),
    })


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


# ---------------------------------------------------------------------
# Deltares HYDRAX (HPC) settings + submission
# ---------------------------------------------------------------------

def _hpc_backend():
    return run_manager.BACKENDS["hpc"]


def _load_hpc_profile(current):
    profile = dict(run_manager.DEFAULT_HPC_PROFILE)
    profile.update(current.load_state().get("hpc_profile") or {})
    # sensible per-project defaults when unset
    if not profile.get("job_name"):
        profile["job_name"] = current.root.name[:24] or "aeolis"
    if not profile.get("config"):
        profile["config"] = current.configfile.name
    return profile


def _save_hpc_profile(current, profile):
    state = current.load_state()
    state["hpc_profile"] = profile          # holds no secrets (password is never persisted)
    current.save_state(state)


@route("GET", "/api/run/hpc")
def _hpc_get(handler, query, tail):
    current = project.require()
    profile = _load_hpc_profile(current)
    backend = _hpc_backend()
    send_json(handler, {
        "profile": profile,
        "partitions": run_manager.HPC_PARTITIONS,
        "script": run_manager.build_job_script(profile),
        "available": backend.available,
        "note": backend.note,
        "project_dir": str(current.root),
        "project_dir_linux": run_manager.local_to_linux(str(current.root)),
    })


@route("POST", "/api/run/hpc")
def _hpc_save(handler, body, tail):
    current = project.require()
    profile = {**run_manager.DEFAULT_HPC_PROFILE, **(body.get("profile") or {})}
    _save_hpc_profile(current, profile)
    _hpc_backend().configure(profile=profile)
    send_json(handler, {"ok": True, "script": run_manager.build_job_script(profile)})


@route("POST", "/api/run/hpc/save_script")
def _hpc_save_script(handler, body, tail):
    """Write the job script to a file the user picked (like saving
    aeolis.txt), so it can be inspected or version-controlled."""
    from pathlib import Path
    path = body.get("path")
    script = body.get("script")
    if not path or script is None:
        send_error_json(handler, "missing 'path' or 'script'")
        return
    try:
        target = Path(path).expanduser()
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(script, newline="\n")
    except OSError as exc:
        send_error_json(handler, f"cannot write script: {exc}", 403)
        return
    send_json(handler, {"ok": True, "path": str(target)})


@route("POST", "/api/run/hpc/start")
def _hpc_start(handler, body, tail):
    current = project.require()
    profile = {**run_manager.DEFAULT_HPC_PROFILE, **(body.get("profile") or {})}
    password = body.get("password")
    script = body.get("script")     # optional raw-script override from the editor
    if not password:
        send_error_json(handler, "enter your HPC password")
        return
    _save_hpc_profile(current, profile)
    backend = _hpc_backend()
    backend.configure(profile=profile, password=password, script=script)
    try:
        run_manager.select("hpc")
        backend.start(current)
    except (ValueError, RuntimeError) as exc:
        send_error_json(handler, exc)
        return
    send_json(handler, {"ok": True, "job_id": backend._job_id})
