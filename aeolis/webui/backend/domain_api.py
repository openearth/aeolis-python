"""Domain data API: raw data acquisition and gridded domain files.

Workflow: select area -> check availability -> download to gui/rawdata
(manifest.json records provenance) -> optionally modify -> interpolate
onto the model grid -> write .grd + update the config. Raw data stays
available as separate layers after interpolation.
"""

import hashlib
from datetime import datetime, timezone

import numpy as np

from aeolis.constants import DEFAULT_CONFIG
from aeolis.webui.backend import datasources, grd_io, jobs, project, settings
from aeolis.webui.backend.config_api import load_config
from aeolis.webui.backend.grid_api import _load_current_grid, patch_config
from aeolis.webui.backend.httpd import route
from aeolis.webui.backend.util import load_json, save_json, send_error_json, send_json

# GUI target name -> (config key, default filename)
TARGETS = {
    "bed": ("bed_file", "zb.grd"),
    "ne": ("ne_file", "ne.grd"),
    "veg": ("veg_file", "veg.grd"),
    "hveg": ("hveg_file", "hveg.grd"),
    "Nt": ("Nt_file", "Nt.grd"),
}
TARGETS = {k: v for k, v in TARGETS.items() if v[0] in DEFAULT_CONFIG}


# ---------------------------------------------------------------------
# manifest helpers
# ---------------------------------------------------------------------

def _manifest_path():
    return project.require().rawdata_dir / settings.MANIFEST_FILE


def load_manifest():
    return load_json(_manifest_path(), default={"entries": []})


def save_manifest(manifest):
    save_json(_manifest_path(), manifest)


def add_entries(new_entries):
    manifest = load_manifest()
    known = {e["id"] for e in manifest["entries"] if "id" in e}
    for entry in new_entries:
        entry["id"] = hashlib.sha1(entry["path"].encode()).hexdigest()[:10]
        if entry["id"] not in known:
            manifest["entries"].append(entry)
            known.add(entry["id"])
    save_manifest(manifest)
    return manifest


def get_entry(entry_id):
    for entry in load_manifest()["entries"]:
        if entry.get("id") == entry_id:
            return entry
    return None


# ---------------------------------------------------------------------
# routes: overview / availability / download / import
# ---------------------------------------------------------------------

@route("GET", "/api/domain")
def _overview(handler, query, tail):
    current = project.require()
    values = load_config(current.configfile)
    targets = {}
    for name, (key, default_name) in TARGETS.items():
        filename = values.get(key) or default_name
        targets[name] = {
            "config_key": key,
            "file": filename,
            "exists": (current.root / filename).is_file(),
            "configured": bool(values.get(key)),
        }
    send_json(handler, {
        "sources": datasources.INFO,
        "entries": load_manifest()["entries"],
        "targets": targets,
        "history": current.load_state().get("domain_history", []),
    })


@route("POST", "/api/domain/check")
def _check(handler, body, tail):
    bounds = body.get("bounds")
    source_ids = body.get("sources") or [s["id"] for s in datasources.INFO if s["id"] != "xyz"]
    if not bounds or len(bounds) != 4:
        send_error_json(handler, "missing 'bounds' [minx,miny,maxx,maxy]")
        return

    def _run(job):
        results = {}
        for source_id in source_ids:
            job.update(message=f"checking {source_id}")
            try:
                results[source_id] = datasources.get(source_id).check(tuple(bounds), job)
            except Exception as exc:  # noqa: BLE001 - reported per source
                results[source_id] = {"available": False, "years": [], "error": str(exc)}
        return results

    send_json(handler, {"job": jobs.start("check availability", _run)})


@route("POST", "/api/domain/download")
def _download(handler, body, tail):
    current = project.require()
    source_id = body.get("source")
    bounds = body.get("bounds")
    years = body.get("years") or []
    if not source_id or not bounds:
        send_error_json(handler, "missing 'source' or 'bounds'")
        return

    def _run(job):
        module = datasources.get(source_id)
        entries = module.download(tuple(bounds), years, current.rawdata_dir, job)
        add_entries(entries)
        return {"entries": entries}

    send_json(handler, {"job": jobs.start(f"download {source_id}", _run)})


@route("POST", "/api/domain/import_xyz")
def _import_xyz(handler, body, tail):
    current = project.require()
    path = body.get("path")
    if not path:
        send_error_json(handler, "missing 'path'")
        return
    from aeolis.webui.backend.datasources import xyz_import
    entry = xyz_import.import_file(path, current.rawdata_dir, crs=body.get("crs"))
    add_entries([entry])
    send_json(handler, {"ok": True, "entry": entry})


@route("POST", "/api/domain/forget")
def _forget(handler, body, tail):
    entry_id = body.get("id")
    manifest = load_manifest()
    entry = next((e for e in manifest["entries"] if e.get("id") == entry_id), None)
    if entry is None:
        send_error_json(handler, f"unknown entry {entry_id}", 404)
        return
    manifest["entries"] = [e for e in manifest["entries"] if e.get("id") != entry_id]
    save_manifest(manifest)
    if body.get("delete_file"):
        target = project.require().root / entry["path"]
        if target.is_file():
            target.unlink()
    send_json(handler, {"ok": True})


# ---------------------------------------------------------------------
# raw-data loading (shared with interpolation and the viewer)
# ---------------------------------------------------------------------

def load_raw(entry):
    """Return ('raster', x1d, y1d, Z) or ('points', x, y, z)."""
    path = project.require().root / entry["path"]
    kind = entry.get("kind")
    if kind == "raster":
        import rasterio
        with rasterio.open(path) as ds:
            z = ds.read(1).astype("float32")
            nodata = ds.nodata
            if nodata is not None:
                z[z == np.float32(nodata)] = np.nan
            tr = ds.transform
            x = tr.c + tr.a * (np.arange(ds.width) + 0.5)
            y = tr.f + tr.e * (np.arange(ds.height) + 0.5)
        return "raster", x, y, z
    if kind == "raster_nc":
        data = np.load(path)
        return "raster", np.asarray(data["x"]), np.asarray(data["y"]), np.asarray(data["z"])
    if kind == "points":
        data = np.load(path)
        return "points", np.asarray(data["x"]), np.asarray(data["y"]), np.asarray(data["z"])
    raise ValueError(f"unknown raw data kind '{kind}'")


def _sample_raster(x, y, Z, XI, YI):
    """Bilinear sample of a raster (x ascending; y any order) at XI/YI."""
    from scipy.interpolate import RegularGridInterpolator
    if y[0] > y[-1]:
        y = y[::-1]
        Z = Z[::-1, :]
    interp = RegularGridInterpolator(
        (y, x), Z, method="linear", bounds_error=False, fill_value=np.nan
    )
    return interp(np.column_stack([YI.ravel(), XI.ravel()])).reshape(XI.shape)


# ---------------------------------------------------------------------
# interpolation onto the model grid
# ---------------------------------------------------------------------

@route("POST", "/api/domain/interpolate")
def _interpolate(handler, body, tail):
    current = project.require()
    target = body.get("target")
    layer_ids = body.get("layers") or []
    fill = body.get("fill")
    if target not in TARGETS:
        send_error_json(handler, f"unknown target '{target}'")
        return
    if not layer_ids:
        send_error_json(handler, "no source layers selected")
        return

    grids, values = _load_current_grid()
    if grids is None:
        send_error_json(handler, "no model grid - generate one in the Grid tab first")
        return
    X, Y = grids

    def _run(job):
        result = np.full(X.shape, np.nan, dtype="float64")
        for n, entry_id in enumerate(layer_ids):
            entry = get_entry(entry_id)
            if entry is None:
                continue
            job.update(progress=n / len(layer_ids), message=f"interpolating {entry.get('label', entry_id)}")
            kind, x, y, z = load_raw(entry)
            hole = ~np.isfinite(result)
            if not hole.any():
                break
            if kind == "raster":
                sampled = _sample_raster(x, y, z, X, Y)
                result[hole] = sampled[hole]
            else:
                from scipy.interpolate import griddata
                pts = np.column_stack([x, y])
                sampled = griddata(pts, z, (X[hole], Y[hole]), method="linear")
                result[hole] = sampled

        missing = int(np.sum(~np.isfinite(result)))
        if missing and fill is not None:
            result[~np.isfinite(result)] = float(fill)
            missing = 0

        if missing:
            raise RuntimeError(
                f"{missing} of {result.size} grid cells have no data - "
                "add more sources or set a fill value"
            )

        key, default_name = TARGETS[target]
        filename = values.get(key) or default_name
        grd_io.write_grd(current.root / filename, result)
        patch_config({key: filename})
        _log_history({
            "action": "interpolate", "target": target, "file": filename,
            "layers": layer_ids, "fill": fill,
        })
        return {"target": target, "file": filename,
                "min": float(np.nanmin(result)), "max": float(np.nanmax(result))}

    send_json(handler, {"job": jobs.start(f"interpolate {target}", _run)})


# ---------------------------------------------------------------------
# modifications (polygon / index-range) and duplication
# ---------------------------------------------------------------------

def _load_target(values, target, X, create_value=None):
    current = project.require()
    key, default_name = TARGETS[target]
    filename = values.get(key) or default_name
    path = current.root / filename
    if path.is_file():
        Z = grd_io.read_grd(path)
        if Z.shape != X.shape:
            raise RuntimeError(
                f"{filename} shape {Z.shape} does not match the grid {X.shape}"
            )
    elif create_value is not None:
        Z = np.full(X.shape, float(create_value))
    else:
        raise RuntimeError(f"{filename} does not exist yet - interpolate or initialize first")
    return key, filename, path, Z


@route("POST", "/api/domain/modify")
def _modify(handler, body, tail):
    target = body.get("target")
    op = body.get("op")
    value = body.get("value")
    if target not in TARGETS:
        send_error_json(handler, f"unknown target '{target}'")
        return
    if op not in ("set", "add", "subtract", "multiply", "min", "max"):
        send_error_json(handler, f"unknown op '{op}'")
        return
    try:
        value = float(value)
    except (TypeError, ValueError):
        send_error_json(handler, "missing numeric 'value'")
        return

    grids, values = _load_current_grid()
    if grids is None:
        send_error_json(handler, "no model grid yet")
        return
    X, Y = grids

    try:
        key, filename, path, Z = _load_target(values, target, X, create_value=body.get("init"))
    except RuntimeError as exc:
        send_error_json(handler, exc)
        return

    # build the selection mask
    mask = np.ones(X.shape, dtype=bool)
    polygon_id = body.get("polygon")
    indices = body.get("indices")
    if polygon_id:
        polygons = project.require().load_polygons()
        obj = next((o for o in polygons.get("objects", []) if o.get("id") == polygon_id), None)
        if obj is None:
            send_error_json(handler, f"unknown polygon {polygon_id}", 404)
            return
        from matplotlib.path import Path as MplPath
        poly = MplPath(np.asarray(obj["coords"], dtype=float))
        mask = poly.contains_points(
            np.column_stack([X.ravel(), Y.ravel()])
        ).reshape(X.shape)
    elif indices:
        j0, j1, i0, i1 = [int(v) for v in indices]
        mask = np.zeros(X.shape, dtype=bool)
        mask[j0:j1 + 1, i0:i1 + 1] = True

    if not mask.any():
        send_error_json(handler, "selection covers no grid cells")
        return

    ops = {
        "set": lambda a: np.full_like(a, value),
        "add": lambda a: a + value,
        "subtract": lambda a: a - value,
        "multiply": lambda a: a * value,
        "min": lambda a: np.minimum(a, value),
        "max": lambda a: np.maximum(a, value),
    }
    Z[mask] = ops[op](Z[mask])
    grd_io.write_grd(path, Z)
    if not values.get(key):
        patch_config({key: filename})

    _log_history({
        "action": "modify", "target": target, "op": op, "value": value,
        "polygon": polygon_id, "indices": indices, "cells": int(mask.sum()),
    })
    send_json(handler, {
        "ok": True, "cells": int(mask.sum()), "file": filename,
        "min": float(np.nanmin(Z)), "max": float(np.nanmax(Z)),
    })


@route("POST", "/api/domain/duplicate")
def _duplicate(handler, body, tail):
    source = body.get("from", "bed")
    target = body.get("to")
    offset = float(body.get("offset", 0.0))
    if source not in TARGETS or target not in TARGETS:
        send_error_json(handler, "unknown 'from'/'to' target")
        return

    grids, values = _load_current_grid()
    if grids is None:
        send_error_json(handler, "no model grid yet")
        return
    X, _ = grids
    try:
        _, src_file, _, Z = _load_target(values, source, X)
    except RuntimeError as exc:
        send_error_json(handler, exc)
        return

    key, default_name = TARGETS[target]
    filename = values.get(key) or default_name
    grd_io.write_grd(project.require().root / filename, Z + offset)
    patch_config({key: filename})
    _log_history({
        "action": "duplicate", "from": source, "to": target,
        "offset": offset, "file": filename,
    })
    send_json(handler, {"ok": True, "file": filename, "offset": offset, "from": src_file})


def _log_history(record):
    current = project.require()
    state = current.load_state()
    history = state.setdefault("domain_history", [])
    record["time"] = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M")
    history.append(record)
    state["domain_history"] = history[-100:]
    current.save_state(state)
