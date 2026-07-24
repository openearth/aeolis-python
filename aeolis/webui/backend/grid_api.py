"""Grid generation API (Grid tab).

The frontend computes live previews itself while drawing; this API is
the single writer of x.grd / y.grd and the reader of existing grids.
"""

import os
from pathlib import Path

import aeolis.inout
from aeolis.webui.backend import grd_io, project
from aeolis.webui.backend.config_api import load_config
from aeolis.webui.backend.httpd import route
from aeolis.webui.backend.util import send_error_json, send_json


def patch_config(patch):
    """Merge *patch* into the project's aeolis.txt (via aeolis.inout)."""
    current = project.require()
    values = aeolis.inout.read_configfile(
        str(current.configfile), parse_files=False, load_defaults=True
    )
    values.update(patch)
    aeolis.inout.write_configfile(str(current.configfile), values)
    return values


def resolve_target(current, raw, default):
    """Resolve a save target that may be an absolute path (from the file
    browser) or a project-relative name. Returns (write_path, config_ref)
    where config_ref is relative when the target lives inside the project
    root, else absolute — so aeolis.txt stays portable when possible.

    Both the target and the root are normalized first (collapsing any ``..``
    / ``.`` segments) before the relative check. The file browser can hand
    back paths like ``…/input/../input/wind.txt`` when the user navigates up
    and back down; without normalization ``relative_to`` keeps the ``..`` in
    the remainder and the config ends up with ``..\\input\\wind.txt``."""
    raw = (raw or default).strip() if isinstance(raw, str) else default
    root = Path(os.path.normpath(str(current.root)))
    path = Path(raw)
    if not path.is_absolute():
        path = root / path
    path = Path(os.path.normpath(str(path)))   # collapse .. and .
    try:
        return path, str(path.relative_to(root))   # clean relative, no ..
    except ValueError:
        return path, str(path)                      # truly outside → absolute


def _load_current_grid():
    """Load X/Y from the config's grid files, or None."""
    current = project.require()
    values = load_config(current.configfile)
    xf, yf = values.get("xgrid_file"), values.get("ygrid_file")
    if not xf or not yf:
        return None, values
    xpath = current.root / xf
    ypath = current.root / yf
    if not xpath.is_file() or not ypath.is_file():
        return None, values
    X = grd_io.read_grd(xpath)
    Y = grd_io.read_grd(ypath)
    return (X, Y), values


@route("GET", "/api/grid")
def _get(handler, query, tail):
    grids, values = _load_current_grid()
    if grids is None:
        send_json(handler, {"exists": False})
        return
    X, Y = grids
    payload = {
        "exists": True,
        "params": grd_io.derive_params(X, Y),
        "geometry": grd_io.geometry(X, Y),
        "boundary_types": {
            "offshore": values.get("boundary_offshore"),
            "onshore": values.get("boundary_onshore"),
            "lateral": values.get("boundary_lateral"),
        },
        "files": {
            "xgrid_file": values.get("xgrid_file"),
            "ygrid_file": values.get("ygrid_file"),
        },
    }
    send_json(handler, payload)


@route("POST", "/api/grid/save")
def _save(handler, body, tail):
    current = project.require()
    try:
        params = {k: float(body[k]) for k in ("x0", "y0", "dx", "rotation")}
        params["nx"] = int(body["nx"])
        params["ny"] = int(body["ny"])
    except (KeyError, TypeError, ValueError) as exc:
        send_error_json(handler, f"invalid grid parameters: {exc}")
        return

    if params["nx"] * params["ny"] > 4_000_000:
        send_error_json(handler, "grid too large (> 4M cells)")
        return

    X, Y = grd_io.generate(
        params["x0"], params["y0"], params["dx"],
        params["nx"], params["ny"], params["rotation"],
    )

    xpath, xfile = resolve_target(current, body.get("xgrid_file"), "x.grd")
    ypath, yfile = resolve_target(current, body.get("ygrid_file"), "y.grd")
    xpath.parent.mkdir(parents=True, exist_ok=True)
    ypath.parent.mkdir(parents=True, exist_ok=True)
    grd_io.write_grd(xpath, X)
    grd_io.write_grd(ypath, Y)

    patch_config({
        "xgrid_file": xfile,
        "ygrid_file": yfile,
        "nx": params["nx"],
        "ny": params["ny"],
        "alfa": 0,
    })

    # persist the accepted params as the grid draft in project state
    state = current.load_state()
    state["grid"] = params
    current.save_state(state)

    send_json(handler, {
        "ok": True,
        "params": grd_io.derive_params(X, Y),
        "geometry": grd_io.geometry(X, Y),
        "files": {"xgrid_file": xfile, "ygrid_file": yfile},
    })


@route("POST", "/api/grid/load")
def _load(handler, body, tail):
    """Load an existing x/y .grd pair from disk, copy it into the project
    root under its basenames and point the config at it."""
    current = project.require()
    xpath = body.get("xgrid_file")
    ypath = body.get("ygrid_file")
    if not xpath or not ypath:
        send_error_json(handler, "need both an x-grid and a y-grid file")
        return
    try:
        X = grd_io.read_grd(Path(xpath))
        Y = grd_io.read_grd(Path(ypath))
    except Exception as exc:  # noqa: BLE001 - report any read failure to the UI
        send_error_json(handler, f"could not read grid: {exc}")
        return
    if X.shape != Y.shape:
        send_error_json(handler, "x-grid and y-grid have different shapes")
        return
    params = grd_io.derive_params(X, Y)
    if params is None:
        send_error_json(handler, "grid is too small or not a valid rectangular grid")
        return

    xfile = os.path.basename(xpath)
    yfile = os.path.basename(ypath)
    grd_io.write_grd(current.root / xfile, X)
    grd_io.write_grd(current.root / yfile, Y)

    patch_config({
        "xgrid_file": xfile,
        "ygrid_file": yfile,
        "nx": int(params["nx"]),
        "ny": int(params["ny"]),
        "alfa": 0,
    })
    state = current.load_state()
    state["grid"] = {k: params[k] for k in ("x0", "y0", "dx", "nx", "ny", "rotation")}
    current.save_state(state)

    send_json(handler, {
        "ok": True,
        "params": params,
        "geometry": grd_io.geometry(X, Y),
        "files": {"xgrid_file": xfile, "ygrid_file": yfile},
    })


@route("GET", "/api/grid/shear")
def _shear(handler, query, tail):
    grids, values = _load_current_grid()
    if grids is None:
        send_error_json(handler, "no grid available", 404)
        return
    X, Y = grids
    udir = float(query.get("udir", 270.0))
    dx_c = float(query.get("dx", values.get("dx") or 1.0))
    dy_c = float(query.get("dy", values.get("dy") or 1.0))
    buffer_width = float(query.get("buffer_width", values.get("buffer_width") or 10.0))
    send_json(handler, grd_io.shear_preview(X, Y, dx_c, dy_c, buffer_width, udir))
