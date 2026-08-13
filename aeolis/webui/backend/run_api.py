"""Run tab API routes."""

import posixpath

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
    # a local run writes the project's own output file — stop following a
    # remote run folder so the Viewer shows what is being computed here
    from aeolis.webui.backend import output_api
    try:
        output_api.set_source(current, None)
    except Exception:  # noqa: BLE001 - the run itself already started
        pass
    send_json(handler, {"ok": True})


@route("POST", "/api/run/stop")
def _stop(handler, body, tail):
    stopped = run_manager.active.stop()
    send_json(handler, {"ok": True, "stopped": stopped})


@route("GET", "/api/run/status")
def _status(handler, query, tail):
    status = run_manager.active.status()
    # where the Viewer reads output from, so the frontend can react to
    # backend-side switches (e.g. the automatic one after an HPC submit)
    from aeolis.webui.backend import output_api
    try:
        status["output"] = output_api.source_summary()
    except Exception:  # noqa: BLE001 - no open project
        pass
    send_json(handler, status)


@route("GET", "/api/run/log")
def _log(handler, query, tail):
    offset = int(query.get("offset", 0))
    send_json(handler, run_manager.active.log_tail(offset))


# ---------------------------------------------------------------------
# Deltares HYDRAX (HPC) settings + submission
# ---------------------------------------------------------------------

def _hpc_backend():
    return run_manager.BACKENDS["hpc"]


def _output_subdir(current):
    """Directory part of the config's output_file (POSIX form), e.g.
    'output' for output/aeolis.nc. The model does not create missing
    directories, so the job script must mkdir it before the run."""
    try:
        out = str(load_config(current.configfile).get("output_file") or "")
    except Exception:
        return ""
    out = out.replace("\\", "/")
    return posixpath.dirname(out)


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
        "output_subdir": _output_subdir(current),
        # jobs previously submitted from this project (for reattaching
        # after the GUI was closed and reopened)
        "jobs": run_manager.recorded_hpc_jobs(current),
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
    # submission continues in the background - stages stream into the log
    send_json(handler, {"ok": True})


@route("POST", "/api/run/hpc/jobs")
def _hpc_jobs(handler, body, tail):
    """The user's jobs on the cluster (live via squeue, merged with the
    jobs recorded for this project) — POST because it carries the
    password needed to open the SSH connection."""
    current = project.require()
    password = body.get("password")
    if not password:
        send_error_json(handler, "enter your HPC password")
        return
    profile = {**run_manager.DEFAULT_HPC_PROFILE, **(body.get("profile") or {})}
    backend = _hpc_backend()
    backend.configure(profile=profile, password=password)
    try:
        jobs = backend.list_jobs(run_manager.recorded_hpc_jobs(current))
    except Exception as exc:  # noqa: BLE001 - ssh/squeue errors go to the user
        send_error_json(handler, exc)
        return
    send_json(handler, {"jobs": jobs})


@route("POST", "/api/run/hpc/attach")
def _hpc_attach(handler, body, tail):
    """Re-attach the monitor to an already-submitted SLURM job, e.g.
    after the GUI was closed and reopened while the job kept running."""
    current = project.require()
    job_id = str(body.get("job_id") or "").strip()
    password = body.get("password")
    if not job_id:
        send_error_json(handler, "missing 'job_id'")
        return
    if not password:
        send_error_json(handler, "enter your HPC password")
        return
    profile = {**run_manager.DEFAULT_HPC_PROFILE, **(body.get("profile") or {})}
    _save_hpc_profile(current, profile)
    backend = _hpc_backend()
    backend.configure(profile=profile, password=password)
    # a job recorded from this project knows its exact stdout path;
    # other jobs are resolved via scontrol on the attach thread
    out_file = None
    for entry in run_manager.recorded_hpc_jobs(current):
        if str(entry.get("job_id")) == job_id and entry.get("run_dir"):
            name = entry.get("job_name") or "aeolis"
            out_file = posixpath.join(entry["run_dir"], f"{name}.o{job_id}")
            break
    try:
        run_manager.select("hpc")
        backend.attach(job_id, out_file=out_file)
    except (ValueError, RuntimeError) as exc:
        send_error_json(handler, exc)
        return
    send_json(handler, {"ok": True})
