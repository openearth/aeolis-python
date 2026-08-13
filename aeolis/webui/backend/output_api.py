"""Field-data streaming API for the map viewer.

Serves three kinds of fields, all rendered by the same WebGL layer:

- model output (netCDF, per-timestep float32 slabs with an LRU cache so
  scrubbing the time slider is instantaneous),
- domain .grd files (static),
- downloaded raw rasters (decimated to a display resolution).

Binary layout of a "gridfield" payload (application/octet-stream)::

    uint32 n, uint32 s                  # rows, cols
    float32 x[n*s], float32 y[n*s]      # model CRS coordinates
    float32 v[n*s]                      # values (NaN -> -1e30 sentinel)

Per-timestep "field" payloads are only ``float32 v[n*s]`` with the value
range in the ``X-Data-Range`` response header.
"""

import re
import struct
from collections import OrderedDict
from pathlib import Path

import numpy as np

from aeolis.webui.backend import grd_io, project
from aeolis.webui.backend.config_api import load_config
from aeolis.webui.backend.conditions_api import parse_refdate
from aeolis.webui.backend.domain_api import (
    TARGETS, get_entry, get_target_draft, load_raw)
from aeolis.webui.backend.grid_api import _load_current_grid
from aeolis.webui.backend.httpd import route
from aeolis.webui.backend.run_manager import linux_to_local, recorded_hpc_jobs
from aeolis.webui.backend.util import NC_LOCK, send_bytes, send_error_json, send_json

SENTINEL = -1.0e30
_slab_cache = OrderedDict()   # (path, var, t, extra) -> (bytes, vmin, vmax)
SLAB_CACHE_MAX = 80
_meta_cache = {}              # path -> (mtime, meta)

# Output-source override: the Viewer normally reads the netCDF named by
# the project's config, from the project folder. When a simulation of
# this project runs elsewhere (e.g. the /p copy submitted to the HPC),
# the Viewer can point at that run folder instead. The override is
# persisted in the project's state file, so reopening the GUI keeps
# following the run — until the user switches back or starts a local
# run (which clears it). It is set from three places: the Viewer's
# output-source switch, the Run tab's job monitor, and automatically by
# the HPC backend right after a copy-mode submission.
_source_cache = None          # {"root": str, "dir": Path|None, "config": str|None}


def _load_source():
    """The persisted override of the open project (dir None = read the
    project's own output file). Cached per project root."""
    global _source_cache
    try:
        current = project.require()
    except RuntimeError:
        return None
    root = str(current.root)
    if _source_cache is None or _source_cache.get("root") != root:
        saved = current.load_state().get("output_source") or {}
        raw = saved.get("dir")
        _source_cache = {"root": root,
                         "dir": Path(raw) if raw else None,
                         "config": saved.get("config")}
    return _source_cache


def set_source(current, dir_str, config=None):
    """Point *current*'s Viewer at another run folder (``dir_str`` may
    be in the cluster's /p form), or back at the project itself
    (``dir_str=None``). Persisted in the project state; raises
    FileNotFoundError when the folder is not accessible."""
    global _source_cache
    folder = None
    if dir_str:
        raw = str(dir_str).strip()
        if re.match(r"^/[A-Za-z](/|$)", raw):
            raw = linux_to_local(raw)
        folder = Path(raw).expanduser()
        if not folder.is_dir():
            raise FileNotFoundError(
                f"'{folder}' is not an accessible folder (is the drive mounted?)")
    config = str(config or "").strip() or None
    state = current.load_state()
    if folder is None:
        state.pop("output_source", None)
    else:
        state["output_source"] = {"dir": str(folder), "config": config}
    current.save_state(state)
    _source_cache = {"root": str(current.root), "dir": folder,
                     "config": config if folder else None}


def _source_active():
    src = _load_source()
    return bool(src and src["dir"])


def source_summary():
    """Cheap {override, dir} for change detection (no config parsing) —
    embedded in /api/run/status so the frontend notices backend-side
    switches (e.g. the automatic one after an HPC submission)."""
    src = _load_source()
    active = bool(src and src["dir"])
    return {"override": active, "dir": str(src["dir"]) if active else None}


def _output_path():
    current = project.require()
    src = _load_source()
    if src and src["dir"]:
        cfg = src["dir"] / (src.get("config") or current.configfile.name)
        if not cfg.is_file():
            cfg = current.configfile
        values = load_config(cfg)
        name = str(values.get("output_file") or (cfg.stem + ".nc"))
        if re.match(r"^/[A-Za-z](/|$)", name):
            # a linuxified copied config may hold the /p form of the path
            name = linux_to_local(name)
        return src["dir"] / name, values
    values = load_config(current.configfile)
    name = values.get("output_file") or (current.configfile.stem + ".nc")
    return current.root / str(name), values


def _known_remote():
    """The most recent HPC submission of this project whose run folder
    is not the project folder itself — offered as the "HPC run" side of
    the Viewer's output-source switch."""
    current = project.require()
    for entry in recorded_hpc_jobs(current):
        raw = entry.get("run_dir")
        if not raw:
            continue
        local = linux_to_local(raw) if re.match(r"^/[A-Za-z](/|$)", raw) else raw
        if Path(local) != current.root:
            return {"dir": raw, "config": entry.get("config"),
                    "job_id": entry.get("job_id")}
    return None


def _source_info():
    """Where the output is effectively read from, the active override
    and the known remote run folder the user could switch to."""
    path, _values = _output_path()
    src = _load_source()
    active = bool(src and src["dir"])
    return {
        "override": active,
        "dir": str(src["dir"]) if active else None,
        "path": str(path),
        "exists": path.is_file(),
        "known": _known_remote(),
    }


@route("GET", "/api/output/source")
def _source_get(handler, query, tail):
    send_json(handler, _source_info())


@route("POST", "/api/output/source")
def _source_set(handler, body, tail):
    """Point the Viewer at another run folder of this project (body
    ``{"dir": ..., "config": ...}``; the dir may be given in the
    cluster's /p form) or back at the project itself (``{"dir": null}``)."""
    current = project.require()
    try:
        set_source(current, body.get("dir"), body.get("config"))
    except FileNotFoundError as exc:
        send_error_json(handler, str(exc), 404)
        return
    send_json(handler, _source_info())


def _open_nc(path):
    import netCDF4
    return netCDF4.Dataset(str(path), mode="r")


@route("GET", "/api/output/meta")
def _meta(handler, query, tail):
    path, values = _output_path()
    if not path.is_file():
        send_json(handler, {"exists": False, "file": str(path.name),
                            "path": str(path), "override": _source_active()})
        return

    mtime = path.stat().st_mtime
    cached = _meta_cache.get(str(path))
    if cached and cached[0] == mtime:
        send_json(handler, {**cached[1], "override": _source_active()})
        return

    refdate = parse_refdate(values)
    with NC_LOCK:
        ds = _open_nc(path)
        try:
            n = ds.dimensions["n"].size
            s = ds.dimensions["s"].size
            tvar = ds.variables["time"]
            times = np.asarray(tvar[:], dtype="float64")
            variables = []
            for name, var in ds.variables.items():
                dims = var.dimensions
                if "time" not in dims or name == "time":
                    continue
                if "n" not in dims or "s" not in dims:
                    continue
                extra = [d for d in dims if d not in ("time", "n", "s")]
                variables.append({
                    "name": name,
                    "units": getattr(var, "units", ""),
                    "long_name": getattr(var, "long_name", name),
                    "extra_dims": [
                        {"name": d, "size": ds.dimensions[d].size} for d in extra
                    ],
                })
        finally:
            ds.close()

    epoch0 = refdate.timestamp()
    meta = {
        "exists": True,
        "file": path.name,
        "path": str(path),
        "shape": [int(n), int(s)],
        "times": times.tolist(),
        "times_epoch": (epoch0 + times).tolist(),
        "variables": variables,
    }
    _meta_cache[str(path)] = (mtime, meta)
    send_json(handler, {**meta, "override": _source_active()})


@route("GET", "/api/output/mesh")
def _mesh(handler, query, tail):
    path, values = _output_path()
    if not path.is_file():
        send_error_json(handler, "no output file", 404)
        return
    with NC_LOCK:
        ds = _open_nc(path)
        try:
            x = np.asarray(ds.variables["x"][:], dtype="float32")
            y = np.asarray(ds.variables["y"][:], dtype="float32")
        finally:
            ds.close()
    n, s = x.shape
    payload = struct.pack("<II", n, s) + x.tobytes() + y.tobytes()
    send_bytes(handler, payload)


@route("GET", "/api/output/field")
def _field(handler, query, tail):
    path, values = _output_path()
    if not path.is_file():
        send_error_json(handler, "no output file", 404)
        return
    var = query.get("var")
    t = int(query.get("t", 0))
    extra = query.get("k", "0")   # comma-separated indices for extra dims
    if not var:
        send_error_json(handler, "missing 'var'")
        return

    key = (str(path), path.stat().st_mtime, var, t, extra)
    cached = _slab_cache.get(key)
    if cached is None:
        with NC_LOCK:
            ds = _open_nc(path)
            try:
                if var not in ds.variables:
                    send_error_json(handler, f"unknown variable '{var}'", 404)
                    return
                ncvar = ds.variables[var]
                dims = ncvar.dimensions
                index = []
                extra_idx = [int(v) for v in str(extra).split(",") if v != ""]
                pos = 0
                for d in dims:
                    if d == "time":
                        index.append(min(t, ncvar.shape[0] - 1))
                    elif d in ("n", "s"):
                        index.append(slice(None))
                    else:
                        k = extra_idx[pos] if pos < len(extra_idx) else 0
                        index.append(min(k, ds.dimensions[d].size - 1))
                        pos += 1
                data = ncvar[tuple(index)]
            finally:
                ds.close()
        arr = np.ma.filled(np.ma.masked_invalid(data), np.nan).astype("float32")
        finite = np.isfinite(arr)
        vmin = float(np.nanmin(arr)) if finite.any() else 0.0
        vmax = float(np.nanmax(arr)) if finite.any() else 1.0
        arr[~finite] = SENTINEL
        cached = (arr.tobytes(), vmin, vmax)
        _slab_cache[key] = cached
        while len(_slab_cache) > SLAB_CACHE_MAX:
            _slab_cache.popitem(last=False)
    else:
        _slab_cache.move_to_end(key)

    payload, vmin, vmax = cached
    send_bytes(handler, payload, extra_headers={"X-Data-Range": f"{vmin},{vmax}"})


@route("GET", "/api/output/series")
def _series(handler, query, tail):
    """Time series of one variable at a single grid cell (probe)."""
    path, values = _output_path()
    if not path.is_file():
        send_error_json(handler, "no output file", 404)
        return
    var = query.get("var")
    j = int(query.get("j", 0))
    i = int(query.get("i", 0))
    extra = [int(v) for v in str(query.get("k", "0")).split(",") if v != ""]
    refdate = parse_refdate(values)

    with NC_LOCK:
        ds = _open_nc(path)
        try:
            if var not in ds.variables:
                send_error_json(handler, f"unknown variable '{var}'", 404)
                return
            ncvar = ds.variables[var]
            index = []
            pos = 0
            for d in ncvar.dimensions:
                if d == "time":
                    index.append(slice(None))
                elif d == "n":
                    index.append(min(j, ds.dimensions["n"].size - 1))
                elif d == "s":
                    index.append(min(i, ds.dimensions["s"].size - 1))
                else:
                    k = extra[pos] if pos < len(extra) else 0
                    index.append(min(k, ds.dimensions[d].size - 1))
                    pos += 1
            data = np.ma.filled(np.ma.masked_invalid(ncvar[tuple(index)]), np.nan)
            times = np.asarray(ds.variables["time"][:], dtype="float64")
        finally:
            ds.close()

    epoch0 = refdate.timestamp()
    send_json(handler, {
        "var": var, "j": j, "i": i,
        "t_epoch": (epoch0 + times).tolist(),
        "values": [None if not np.isfinite(v) else float(v) for v in np.atleast_1d(data)],
    })


# ---------------------------------------------------------------------
# static gridfields: domain .grd files and raw rasters
# ---------------------------------------------------------------------

def _pack_gridfield(X, Y, V):
    V = np.asarray(V, dtype="float32").copy()
    finite = np.isfinite(V)
    vmin = float(V[finite].min()) if finite.any() else 0.0
    vmax = float(V[finite].max()) if finite.any() else 1.0
    V[~finite] = SENTINEL
    n, s = X.shape
    payload = (struct.pack("<II", n, s)
               + np.asarray(X, dtype="float32").tobytes()
               + np.asarray(Y, dtype="float32").tobytes()
               + V.tobytes())
    return payload, vmin, vmax


@route("GET", "/api/domain/gridfield")
def _gridfield(handler, query, tail):
    target = query.get("target")
    if target not in TARGETS:
        send_error_json(handler, f"unknown target '{target}'")
        return
    grids, values = _load_current_grid()
    if grids is None:
        send_error_json(handler, "no model grid", 404)
        return
    X, Y = grids
    key, default_name = TARGETS[target]
    filename = values.get(key) or default_name
    # an unsaved interpolation draft takes precedence, so it previews on the
    # map before the user commits it to a .grd
    draft = get_target_draft(target)
    if draft is not None:
        V = np.asarray(draft, dtype="float64")
    else:
        path = project.require().root / str(filename)
        if not path.is_file():
            send_error_json(handler, f"{filename} does not exist", 404)
            return
        V = grd_io.read_grd(path)
    n_species = 1
    if V.shape != X.shape:
        # species-stacked file (hveg/Nt: flat ny*nx*nspecies, cf.
        # grass.py reshape) -> slice out the requested species
        if target in ("hveg", "Nt") and V.size % X.size == 0 and V.size // X.size >= 1:
            n_species = V.size // X.size
            k = max(0, min(n_species - 1, int(query.get("k", 0) or 0)))
            V = V.reshape(X.shape[0], X.shape[1], n_species)[:, :, k]
        else:
            send_error_json(
                handler,
                f"{filename} was made for a different grid "
                f"({V.shape[0]}x{V.shape[1] if V.ndim > 1 else 1} vs "
                f"{X.shape[0]}x{X.shape[1]}) - re-interpolate it in the Domain tab",
                409)
            return
    payload, vmin, vmax = _pack_gridfield(X, Y, V)
    send_bytes(handler, payload, extra_headers={
        "X-Data-Range": f"{vmin},{vmax}",
        "X-Species": str(n_species),
    })


@route("GET", "/api/domain/rawfield")
def _rawfield(handler, query, tail):
    entry = get_entry(query.get("id", ""))
    if entry is None:
        send_error_json(handler, "unknown raw layer", 404)
        return
    kind, x, y, z = load_raw(entry)
    if kind == "points":
        # decimated point cloud as JSON (rendered as circles). Points
        # that originate from a 2D grid (e.g. a .grd duplicated to
        # samples) are decimated per row/column so the display keeps a
        # regular pattern - a flat stride on the raveled array would
        # produce a staggered checkerboard. The per-axis steps are chosen
        # so the *physical* spacing is roughly equal in both directions,
        # so a grid with dx != dy (or a rotated grid) still shows an
        # even, square-looking dot lattice.
        shape = entry.get("shape")
        if shape and int(shape[0]) * int(shape[1]) == x.size:
            ny, nx = int(shape[0]), int(shape[1])
            Xg = x.reshape(ny, nx)
            Yg = y.reshape(ny, nx)
            drow = float(np.nanmean(np.hypot(np.diff(Xg, axis=0),
                                             np.diff(Yg, axis=0)))) or 1.0
            dcol = float(np.nanmean(np.hypot(np.diff(Xg, axis=1),
                                             np.diff(Yg, axis=1)))) or 1.0
            drow, dcol = abs(drow) or 1.0, abs(dcol) or 1.0
            target = 60000.0
            # equal-physical-spacing S with total points ~ target
            s_phys = np.sqrt(max(ny * nx * drow * dcol / target, 1e-9))
            step_r = max(1, int(round(s_phys / drow)))
            step_c = max(1, int(round(s_phys / dcol)))
            idx = (np.arange(0, ny, step_r)[:, None] * nx
                   + np.arange(0, nx, step_c)[None, :]).ravel()
            x, y, z = x[idx], y[idx], z[idx]
        else:
            stride = max(1, x.size // 20000)
            x, y, z = x[::stride], y[::stride], z[::stride]
        send_json(handler, {
            "kind": "points",
            "x": x.astype(float).tolist(),
            "y": y.astype(float).tolist(),
            # bare NaN is invalid JSON for the browser's JSON.parse
            "z": [float(v) if np.isfinite(v) else None for v in z],
        })
        return
    # raster: decimate to display resolution and pack as gridfield
    max_cells = 700
    sy = max(1, z.shape[0] // max_cells)
    sx = max(1, z.shape[1] // max_cells)
    z2 = z[::sy, ::sx]
    x2 = x[::sx]
    y2 = y[::sy]
    X, Y = np.meshgrid(x2, y2)
    payload, vmin, vmax = _pack_gridfield(X, Y, z2)
    send_bytes(handler, payload, extra_headers={"X-Data-Range": f"{vmin},{vmax}"})
