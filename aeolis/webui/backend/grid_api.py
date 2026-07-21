"""Grid generation API (Grid tab).

The frontend computes live previews itself while drawing; this API is
the single writer of x.grd / y.grd and the reader of existing grids.
"""

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

    xfile = body.get("xgrid_file") or "x.grd"
    yfile = body.get("ygrid_file") or "y.grd"
    grd_io.write_grd(current.root / xfile, X)
    grd_io.write_grd(current.root / yfile, Y)

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
