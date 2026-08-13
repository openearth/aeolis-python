"""Domain data API: raw data acquisition and gridded domain files.

Workflow: select area -> check availability -> download to gui/rawdata
(manifest.json records provenance) -> optionally modify -> interpolate
onto the model grid -> write .grd + update the config. Raw data stays
available as separate layers after interpolation.
"""

import hashlib
import os
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from aeolis.constants import DEFAULT_CONFIG
from aeolis.webui.backend import datasources, grd_io, jobs, project, settings
from aeolis.webui.backend.config_api import load_config
from aeolis.webui.backend.grid_api import _load_current_grid, patch_config, resolve_target
from aeolis.webui.backend.httpd import route
from aeolis.webui.backend.util import load_json, save_json, send_error_json, send_json

# GUI target name -> (config key, default filename)
TARGETS = {
    "bed": ("bed_file", "zb.grd"),
    "ne": ("ne_file", "ne.grd"),
    "veg": ("veg_file", "veg.grd"),
    "hveg": ("hveg_file", "hveg.grd"),
    "Nt": ("Nt_file", "Nt.grd"),
    # 2D spatial files / masks (masks may hold complex values;
    # bedcomp_file is 4D and deliberately not offered here)
    "threshold": ("threshold_file", "threshold.grd"),
    "fence": ("fence_file", "fence.grd"),
    "supply": ("supply_file", "supply.grd"),
    "wave_mask": ("wave_mask", "wave_mask.grd"),
    "tide_mask": ("tide_mask", "tide_mask.grd"),
    "runup_mask": ("runup_mask", "runup_mask.grd"),
    "threshold_mask": ("threshold_mask", "threshold_mask.grd"),
    "gw_mask": ("gw_mask", "gw_mask.grd"),
    "vver_mask": ("vver_mask", "vver_mask.grd"),
}
TARGETS = {k: v for k, v in TARGETS.items() if v[0] in DEFAULT_CONFIG}

# targets shown by default (masks appear when configured or on request)
PRIMARY_TARGETS = ["bed", "ne", "veg", "hveg", "Nt"]


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

def grid_signature(values=None):
    """Cheap fingerprint of the current model grid (for staleness)."""
    grids, _ = _load_current_grid()
    if grids is None:
        return None
    X, Y = grids
    return f"{X.shape[0]}x{X.shape[1]}:{float(X[0, 0]):.3f}:{float(Y[0, 0]):.3f}:{float(X[-1, -1]):.3f}"


@route("GET", "/api/domain")
def _overview(handler, query, tail):
    current = project.require()
    values = load_config(current.configfile)
    state = current.load_state()
    signatures = state.get("interp_signatures", {})
    current_sig = grid_signature()

    # which domain files are actually needed given the current config
    method = str(values.get("method_vegetation", "duran") or "duran")
    proc_veg = bool(values.get("process_vegetation", False))
    VEG = {"veg", "hveg", "Nt"}
    ALWAYS = {"bed", "ne"}

    targets = {}
    for name, (key, default_name) in TARGETS.items():
        filename = values.get(key) or default_name
        exists = (current.root / str(filename)).is_file()
        configured = bool(values.get(key))
        # optional = not a base file and not a vegetation file
        optional = name not in ALWAYS and name not in VEG
        needed, note = True, None
        if name in VEG:
            if not proc_veg:
                needed, note = False, "vegetation process is off"
            elif name == "veg" and method != "duran":
                needed, note = False, f"method_vegetation = {method} (veg not used)"
            elif name in ("hveg", "Nt") and method != "grass":
                needed, note = False, f"method_vegetation = {method} (grass files not used)"
        elif optional:
            needed = False  # masks / threshold / fence / supply are opt-in
        has_draft = has_target_draft(name)
        # optional files with no data yet stay hidden until the user adds
        # them; an unsaved draft counts as data (else the card would
        # disappear from the list on a restart)
        hidden = optional and not configured and not exists and not has_draft
        # vegetation files the current config will never use: hide them too
        if name in VEG and not needed:
            hidden = True
        entry = {
            "config_key": key,
            "file": filename,
            "exists": exists,
            "configured": configured,
            "optional": optional,
            "needed": needed,
            "note": note,
            "hidden": hidden,
            "stale": False,
            "shape_ok": True,
            "has_draft": has_draft,
        }
        if exists and current_sig:
            stored = signatures.get(name)
            entry["stale"] = bool(stored) and stored != current_sig
            grids, _ = _load_current_grid()
            if grids is not None:
                try:
                    Z = grd_io.read_grd(current.root / str(filename))
                    # species-stacked files (hveg/Nt) are fine as long as
                    # the total size is a multiple of the grid size
                    species_ok = (name in ("hveg", "Nt")
                                  and Z.size % grids[0].size == 0
                                  and Z.size // grids[0].size >= 1)
                    entry["shape_ok"] = Z.shape == grids[0].shape or species_ok
                except (ValueError, OSError):
                    entry["shape_ok"] = False
            entry["stale"] = entry["stale"] or not entry["shape_ok"]
        targets[name] = entry

    send_json(handler, {
        "sources": datasources.INFO,
        "entries": load_manifest()["entries"],
        "targets": targets,
        "grid_available": current_sig is not None,
    })


@route("POST", "/api/domain/check")
def _check(handler, body, tail):
    bounds = body.get("bounds")
    source_ids = body.get("sources") or [s["id"] for s in datasources.INFO if s["id"] in datasources.SOURCES]
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


@route("POST", "/api/domain/import_tiff")
def _import_tiff(handler, body, tail):
    current = project.require()
    path = body.get("path")
    if not path:
        send_error_json(handler, "missing 'path'")
        return
    from aeolis.webui.backend.datasources import tiff_import
    entry = tiff_import.import_file(path, current.rawdata_dir, crs=body.get("crs"))
    add_entries([entry])
    send_json(handler, {"ok": True, "entry": entry})


# ---------------------------------------------------------------------
# linking raw data that already lives in another project (no re-download)
# ---------------------------------------------------------------------

def _find_external_manifest(raw):
    """Locate a rawdata manifest from a path the user picked - it may be
    the manifest.json itself, a project root, or a gui/rawdata folder."""
    base = Path(raw).expanduser()
    candidates = [base]
    if base.name != settings.MANIFEST_FILE:
        candidates += [
            base / settings.MANIFEST_FILE,
            base / "gui" / "rawdata" / settings.MANIFEST_FILE,
        ]
    for c in candidates:
        if c.is_file() and c.name == settings.MANIFEST_FILE:
            return c
    return None


def _external_root(manifest_path):
    """Folder that the manifest's relative entry paths resolve against.
    A project manifest lives at <root>/gui/rawdata/manifest.json and its
    entries carry gui/rawdata/... paths relative to <root>; a bare data
    folder (e.g. a shared 01_data/lidar) keeps its manifest next to the
    files with plain filenames as paths."""
    parent = manifest_path.parent
    if parent.name == settings.RAWDATA_DIRNAME and parent.parent.name == settings.GUI_DIRNAME:
        return parent.parent.parent
    return parent


@route("POST", "/api/domain/scan_external")
def _scan_external(handler, body, tail):
    """List the datasets of another project so the user can pick which to
    link. Resolves each entry's file to an absolute path and reports
    whether it still exists on disk."""
    manifest_path = _find_external_manifest(body.get("path") or "")
    if manifest_path is None:
        send_error_json(handler, "no rawdata manifest found there "
                                 "(pick a project folder or its gui/rawdata)", 404)
        return
    src_root = _external_root(manifest_path)
    current = project.require()
    data = load_json(manifest_path, default={"entries": []})
    items = []
    for e in data.get("entries", []):
        abspath = (src_root / e.get("path", "")).resolve()
        exists = abspath.is_file()
        # already part of this project? (file lives under our root)
        is_local = abspath.is_relative_to(current.root.resolve())
        items.append({
            "src_id": e.get("id"),
            "label": e.get("label") or e.get("path"),
            "source": e.get("source"),
            "kind": e.get("kind"),
            "year": e.get("year"),
            "date": e.get("date"),
            "abspath": str(abspath),
            "size": abspath.stat().st_size if exists else 0,
            "exists": exists,
            "is_local": is_local,
        })
    send_json(handler, {"manifest": str(manifest_path), "root": str(src_root), "entries": items})


@route("POST", "/api/domain/link_external")
def _link_external(handler, body, tail):
    """Add manifest entries that POINT AT another project's files
    (absolute path, ``linked`` flag) instead of copying/re-downloading."""
    manifest_path = _find_external_manifest(body.get("path") or "")
    if manifest_path is None:
        send_error_json(handler, "no rawdata manifest found there", 404)
        return
    wanted = set(body.get("ids") or [])
    src_root = _external_root(manifest_path)
    data = load_json(manifest_path, default={"entries": []})

    new_entries = []
    for e in data.get("entries", []):
        if wanted and e.get("id") not in wanted:
            continue
        abspath = (src_root / e.get("path", "")).resolve()
        if not abspath.is_file():
            continue
        entry = dict(e)
        entry["path"] = str(abspath)         # absolute -> resolves from any project
        entry["linked"] = True
        entry["origin"] = str(src_root)
        entry["downloaded"] = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M")
        entry.pop("id", None)                # add_entries re-keys by path
        new_entries.append(entry)

    if not new_entries:
        send_error_json(handler, "no linkable files selected")
        return
    add_entries(new_entries)
    send_json(handler, {"ok": True, "linked": len(new_entries)})


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
        # never delete data that physically lives outside this project
        # (a linked/external reference) - only unlink our own rawdata copy
        if not entry.get("linked") and _under_rawdata(target):
            if target.is_file():
                target.unlink()
    send_json(handler, {"ok": True})


def _under_rawdata(path):
    """True if ``path`` resolves to somewhere inside this project's
    gui/rawdata folder (so it is safe for us to delete/rename)."""
    try:
        rd = project.require().rawdata_dir.resolve()
        p = Path(path).resolve()
        return p == rd or rd in p.parents
    except Exception:  # noqa: BLE001
        return False


# ---------------------------------------------------------------------
# raw-data loading (shared with interpolation and the viewer)
# ---------------------------------------------------------------------

def load_raw(entry):
    """Return ('raster', x1d, y1d, Z) or ('points', x, y, z)."""
    path = project.require().root / entry["path"]
    kind = entry.get("kind")
    if kind == "raster":
        import rasterio
        band = int(entry.get("band", 1))
        with rasterio.open(path) as ds:
            band = min(max(1, band), ds.count)
            z = ds.read(band).astype("float32")
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


def _fill_edge_holes(result, iterations=2):
    """Fill lone NaN cells from their neighbours' mean.

    Linear interpolation (qhull) can leave isolated NaN cells exactly on
    the convex-hull edge of the samples - e.g. when a grid is converted
    to samples and interpolated back onto itself. Two dilation passes
    fill holes at most two cells from valid data; real data gaps stay
    NaN and are still reported. Returns the number of filled cells."""
    import warnings
    filled = 0
    for _ in range(iterations):
        hole = ~np.isfinite(result)
        if not hole.any():
            break
        padded = np.pad(result, 1, mode="constant", constant_values=np.nan)
        shifts = [padded[1 + dj:padded.shape[0] - 1 + dj,
                         1 + di:padded.shape[1] - 1 + di]
                  for dj in (-1, 0, 1) for di in (-1, 0, 1) if (dj, di) != (0, 0)]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            neighbour_mean = np.nanmean(np.stack(shifts), axis=0)
        fillable = hole & np.isfinite(neighbour_mean)
        if not fillable.any():
            break
        result[fillable] = neighbour_mean[fillable]
        filled += int(fillable.sum())
    return filled


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

# ---------------------------------------------------------------------
# interpolated-target drafts
#
# Interpolating (or loading a .grd on top) produces an in-memory "draft"
# grid stored in the gui cache — a first-class entity the user can preview
# on the map before choosing to Save it to the configured path, Save-as a
# new file, or Load a different file on top. Nothing touches aeolis.txt or
# writes a .grd until the user explicitly Saves, so the workflow never
# errors on a target file that does not exist yet.
# ---------------------------------------------------------------------

def _draft_dir():
    d = project.require().cache_dir / "target_drafts"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _draft_path(target):
    return _draft_dir() / f"{target}.npy"


def has_target_draft(target):
    return _draft_path(target).is_file()


def get_target_draft(target):
    """The unsaved interpolation array for *target*, or None."""
    p = _draft_path(target)
    if not p.is_file():
        return None
    try:
        return np.load(p)
    except Exception:  # noqa: BLE001 - a corrupt draft just means "no draft"
        return None


def _set_target_draft(target, Z):
    np.save(_draft_path(target), np.asarray(Z, dtype="float64"))


def _clear_target_draft(target):
    p = _draft_path(target)
    if p.is_file():
        try:
            p.unlink()
        except OSError:
            pass


def _commit_draft(current, target, Z, write_path, config_ref):
    """Write *Z* to disk at *write_path*, repoint the config to *config_ref*,
    record the interpolation signature, and drop the draft."""
    write_path.parent.mkdir(parents=True, exist_ok=True)
    grd_io.write_grd(write_path, Z)
    key, _ = TARGETS[target]
    patch_config({key: config_ref})
    state = current.load_state()
    state.setdefault("interp_signatures", {})[target] = grid_signature()
    state.get("draft_signatures", {}).pop(target, None)
    current.save_state(state)
    _clear_target_draft(target)


@route("POST", "/api/domain/interpolate")
def _interpolate(handler, body, tail):
    current = project.require()
    target = body.get("target")
    layer_ids = body.get("layers") or []
    fill = body.get("fill")
    extrapolate = bool(body.get("extrapolate"))
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
            # "tgt:<name>" = another interpolated grid on the same mesh:
            # copy its values straight into the remaining holes
            if isinstance(entry_id, str) and entry_id.startswith("tgt:"):
                src_t = entry_id[4:]
                if src_t not in TARGETS or src_t == target:
                    continue
                job.update(progress=n / len(layer_ids), message=f"copying {src_t}")
                Z = get_target_draft(src_t)
                if Z is None:
                    try:
                        _, _, _, Z = _load_target(values, src_t, X)
                    except RuntimeError:
                        continue
                Z = np.asarray(Z, dtype="float64")
                if Z.shape != X.shape:
                    continue
                hole = ~np.isfinite(result)
                result[hole] = Z[hole]
                continue
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
        if missing:
            edge_filled = _fill_edge_holes(result)
            if edge_filled:
                job.update(message=f"filled {edge_filled} edge cells from neighbours")
            missing = int(np.sum(~np.isfinite(result)))
        # extrapolate: fill every remaining hole with its nearest valid cell
        if missing and extrapolate:
            valid = np.isfinite(result)
            if valid.any():
                from scipy.interpolate import griddata
                job.update(message=f"extrapolating {missing} cells to nearest neighbour")
                nn = griddata(
                    np.column_stack([X[valid], Y[valid]]), result[valid],
                    (X[~valid], Y[~valid]), method="nearest")
                result[~valid] = nn
                missing = int(np.sum(~np.isfinite(result)))
        if missing and fill is not None:
            result[~np.isfinite(result)] = float(fill)
            missing = 0

        if missing:
            raise RuntimeError(
                f"{missing} of {result.size} grid cells have no data - "
                "add more sources or set a fill value"
            )

        # produce a draft (not written to disk) — the user Saves it explicitly
        _set_target_draft(target, result)
        state = current.load_state()
        state.setdefault("draft_signatures", {})[target] = grid_signature()
        current.save_state(state)
        _log_history({
            "action": "interpolate_draft", "target": target,
            "layers": layer_ids, "fill": fill, "extrapolate": extrapolate,
        })
        return {"target": target, "draft": True,
                "min": float(np.nanmin(result)), "max": float(np.nanmax(result))}

    send_json(handler, {"job": jobs.start(f"interpolate {target}", _run)})


# ---------------------------------------------------------------------
# draft commit / load: save to the configured path, save-as, load on top
# ---------------------------------------------------------------------

@route("POST", "/api/domain/target_save")
def _target_save(handler, body, tail):
    """Write the current draft to the target's *configured* path."""
    current = project.require()
    target = body.get("target")
    if target not in TARGETS:
        send_error_json(handler, f"unknown target '{target}'")
        return
    Z = get_target_draft(target)
    if Z is None:
        send_error_json(handler, "no interpolation draft to save - interpolate first")
        return
    values = load_config(current.configfile)
    key, default_name = TARGETS[target]
    filename = str(values.get(key) or default_name)
    write_path, config_ref = resolve_target(current, filename, str(default_name))
    _commit_draft(current, target, Z, write_path, config_ref)
    _log_history({"action": "save_draft", "target": target, "file": config_ref})
    send_json(handler, {"ok": True, "target": target, "file": config_ref})


@route("POST", "/api/domain/target_save_as")
def _target_save_as(handler, body, tail):
    """Write the current draft (or the saved file, if there is no draft) to a
    chosen path and repoint the config to it."""
    current = project.require()
    target = body.get("target")
    raw = (body.get("path") or body.get("filename") or "").strip()
    if target not in TARGETS:
        send_error_json(handler, f"unknown target '{target}'")
        return
    if not raw:
        send_error_json(handler, "missing path")
        return
    if not raw.lower().endswith(".grd"):
        raw += ".grd"
    key, default_name = TARGETS[target]
    Z = get_target_draft(target)
    if Z is None:
        src_name = load_config(current.configfile).get(key) or default_name
        src = current.root / str(src_name)
        if not src.is_file():
            send_error_json(handler, "nothing to save - interpolate first")
            return
        Z = grd_io.read_grd(src)
    write_path, config_ref = resolve_target(current, raw, str(default_name))
    _commit_draft(current, target, Z, write_path, config_ref)
    _log_history({"action": "save_as", "target": target, "file": config_ref})
    send_json(handler, {"ok": True, "target": target, "file": config_ref})


@route("POST", "/api/domain/target_load")
def _target_load(handler, body, tail):
    """Load an existing .grd from disk into the target *draft* (replace on
    top). Nothing is written to the config until the user Saves."""
    current = project.require()
    target = body.get("target")
    path = body.get("path")
    if target not in TARGETS:
        send_error_json(handler, f"unknown target '{target}'")
        return
    if not path:
        send_error_json(handler, "missing path")
        return
    try:
        Z = grd_io.read_grd(Path(path))
    except Exception as exc:  # noqa: BLE001 - report any read failure to the UI
        send_error_json(handler, f"could not read {path}: {exc}")
        return
    _set_target_draft(target, Z)
    state = current.load_state()
    state.setdefault("draft_signatures", {})[target] = grid_signature()
    current.save_state(state)
    grids, _ = _load_current_grid()
    shape_ok = grids is None or Z.shape == grids[0].shape
    _log_history({"action": "load_draft", "target": target, "file": os.path.basename(path)})
    send_json(handler, {"ok": True, "target": target, "draft": True,
                        "shape_ok": bool(shape_ok),
                        "min": float(np.nanmin(Z)), "max": float(np.nanmax(Z))})


# ---------------------------------------------------------------------
# load an existing .grd into a target / save a target under a new name
# ---------------------------------------------------------------------

@route("POST", "/api/domain/save_target_as")
def _save_target_as(handler, body, tail):
    current = project.require()
    target = body.get("target")
    # accept a full path from the file browser, or a legacy bare filename
    raw = (body.get("path") or body.get("filename") or "").strip()
    if target not in TARGETS:
        send_error_json(handler, f"unknown target '{target}'")
        return
    if not raw:
        send_error_json(handler, "missing path")
        return
    if not raw.lower().endswith(".grd"):
        raw += ".grd"
    key, default_name = TARGETS[target]
    src_name = load_config(current.configfile).get(key) or default_name
    src = current.root / str(src_name)
    if not src.is_file():
        send_error_json(handler, f"{src_name} does not exist yet - interpolate or initialize it first")
        return
    write_path, config_ref = resolve_target(current, raw, str(default_name))
    write_path.parent.mkdir(parents=True, exist_ok=True)
    Z = grd_io.read_grd(src)
    grd_io.write_grd(write_path, Z)
    patch_config({key: config_ref})
    _log_history({"action": "save_as", "target": target, "file": config_ref})
    send_json(handler, {"ok": True, "target": target, "file": config_ref})


@route("POST", "/api/domain/load_target")
def _load_target_file(handler, body, tail):
    current = project.require()
    target = body.get("target")
    path = body.get("path")
    if target not in TARGETS:
        send_error_json(handler, f"unknown target '{target}'")
        return
    if not path:
        send_error_json(handler, "missing path")
        return
    try:
        Z = grd_io.read_grd(Path(path))
    except Exception as exc:  # noqa: BLE001 - report any read failure to the UI
        send_error_json(handler, f"could not read {path}: {exc}")
        return
    filename = os.path.basename(path)
    if not filename.endswith(".grd"):
        filename += ".grd"
    grd_io.write_grd(current.root / filename, Z)
    key, _ = TARGETS[target]
    patch_config({key: filename})
    # this file now defines the target for the current grid
    state = current.load_state()
    state.setdefault("interp_signatures", {})[target] = grid_signature()
    current.save_state(state)
    _log_history({"action": "load", "target": target, "file": filename})
    grids, _ = _load_current_grid()
    shape_ok = grids is None or Z.shape == grids[0].shape
    send_json(handler, {"ok": True, "target": target, "file": filename,
                        "shape_ok": bool(shape_ok)})


# ---------------------------------------------------------------------
# draft-based creation / editing of interpolated targets
#   - make a constant grid (e.g. all 1s for a mask)
#   - apply set/add/multiply/clip over the whole grid, a polygon, or an
#     index range; result stays a DRAFT until the user Saves it
# ---------------------------------------------------------------------

@route("POST", "/api/domain/target_constant")
def _target_constant(handler, body, tail):
    current = project.require()
    target = body.get("target")
    if target not in TARGETS:
        send_error_json(handler, f"unknown target '{target}'")
        return
    try:
        value = float(body.get("value", 0))
    except (TypeError, ValueError):
        send_error_json(handler, "missing numeric value")
        return
    grids, _ = _load_current_grid()
    if grids is None:
        send_error_json(handler, "no model grid yet - create one in the Grid tab")
        return
    X, _ = grids
    Z = np.full(X.shape, value, dtype="float64")
    _set_target_draft(target, Z)
    state = current.load_state()
    state.setdefault("draft_signatures", {})[target] = grid_signature()
    current.save_state(state)
    _log_history({"action": "constant_draft", "target": target, "value": value})
    send_json(handler, {"ok": True, "target": target, "draft": True, "min": value, "max": value})


@route("POST", "/api/domain/target_modify")
def _target_modify(handler, body, tail):
    """Edit the target DRAFT (or the saved .grd, loaded into a draft) with a
    set/add/multiply/clip op over the whole grid, a polygon or an index box."""
    current = project.require()
    target = body.get("target")
    op = body.get("op")
    if target not in TARGETS:
        send_error_json(handler, f"unknown target '{target}'")
        return
    if op not in _OPS:
        send_error_json(handler, f"unknown op '{op}'")
        return
    try:
        value = float(body.get("value"))
    except (TypeError, ValueError):
        send_error_json(handler, "missing numeric value")
        return
    grids, values = _load_current_grid()
    if grids is None:
        send_error_json(handler, "no model grid yet")
        return
    X, Y = grids

    Z = get_target_draft(target)
    if Z is None:
        try:
            _, _, _, Z = _load_target(values, target, X)   # fall back to the saved file
        except RuntimeError as exc:
            send_error_json(handler, exc)
            return
    Z = np.array(Z, dtype="float64", copy=True)
    if Z.shape != X.shape:
        send_error_json(handler, f"grid data shape {Z.shape} does not match the grid {X.shape}")
        return

    scope = body.get("scope") or {"type": "all"}
    mask = np.ones(X.shape, dtype=bool)
    if scope.get("type") == "polygon":
        try:
            mask = _polygon_mask(X, Y, scope.get("polygon"))
        except RuntimeError as exc:
            send_error_json(handler, exc, 404)
            return
    elif scope.get("type") == "indices":
        idx = [int(v) for v in (scope.get("indices") or [])]
        mask = np.zeros(X.shape, dtype=bool)
        if len(idx) >= 4:
            mask[idx[0]:idx[1] + 1, idx[2]:idx[3] + 1] = True
    if not mask.any():
        send_error_json(handler, "selection covers no grid cells")
        return

    Z[mask] = _OPS[op](Z[mask], value)
    _set_target_draft(target, Z)
    state = current.load_state()
    state.setdefault("draft_signatures", {})[target] = grid_signature()
    current.save_state(state)
    _log_history({"action": "modify_draft", "target": target, "op": op,
                  "value": value, "cells": int(mask.sum())})
    send_json(handler, {"ok": True, "target": target, "draft": True, "cells": int(mask.sum()),
                        "min": float(np.nanmin(Z)), "max": float(np.nanmax(Z))})


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
    if op not in _OPS:
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

    Z[mask] = _OPS[op](Z[mask], value)
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


# ---------------------------------------------------------------------
# sample-level operations (modify / duplicate / rename / reorder /
# convert an interpolated .grd back into a sample layer)
# ---------------------------------------------------------------------

_OPS = {
    "set": lambda a, v: np.full_like(a, v),
    "add": lambda a, v: a + v,
    "subtract": lambda a, v: a - v,
    "multiply": lambda a, v: a * v,
    # clip_max caps values ABOVE v down to v; clip_min raises values
    # BELOW v up to v ("min"/"max" kept as legacy aliases)
    "clip_max": lambda a, v: np.minimum(a, v),
    "clip_min": lambda a, v: np.maximum(a, v),
    "min": lambda a, v: np.minimum(a, v),
    "max": lambda a, v: np.maximum(a, v),
}


def _polygon_mask(xs, ys, polygon_id):
    polygons = project.require().load_polygons()
    obj = next((o for o in polygons.get("objects", []) if o.get("id") == polygon_id), None)
    if obj is None:
        raise RuntimeError(f"unknown polygon {polygon_id}")
    from matplotlib.path import Path as MplPath
    poly = MplPath(np.asarray(obj["coords"], dtype=float))
    return poly.contains_points(np.column_stack([np.ravel(xs), np.ravel(ys)])).reshape(np.shape(xs))


@route("POST", "/api/domain/sample_modify")
def _sample_modify(handler, body, tail):
    current = project.require()
    entry = get_entry(body.get("id", ""))
    if entry is None:
        send_error_json(handler, "unknown sample layer", 404)
        return
    op = body.get("op")
    if op not in _OPS:
        send_error_json(handler, f"unknown op '{op}'")
        return
    try:
        value = float(body.get("value"))
    except (TypeError, ValueError):
        send_error_json(handler, "missing numeric 'value'")
        return
    scope = body.get("scope") or {"type": "all"}
    save_as = (body.get("save_as") or "").strip()

    kind, x, y, z = load_raw(entry)
    z = np.array(z, dtype="float64", copy=True)

    if kind == "raster":
        X, Y = np.meshgrid(x, y)
    else:
        X, Y = x, y

    if scope.get("type") == "polygon":
        try:
            mask = _polygon_mask(X, Y, scope.get("polygon"))
        except RuntimeError as exc:
            send_error_json(handler, exc, 404)
            return
    elif scope.get("type") == "indices":
        idx = [int(v) for v in (scope.get("indices") or [])]
        mask = np.zeros(z.shape, dtype=bool)
        if z.ndim == 2 and len(idx) >= 4:
            mask[idx[0]:idx[1] + 1, idx[2]:idx[3] + 1] = True
        elif z.ndim == 1 and len(idx) >= 2:
            mask[idx[0]:idx[1] + 1] = True
        else:
            send_error_json(handler, "invalid indices for this sample layer")
            return
    else:
        mask = np.ones(z.shape, dtype=bool)

    if not mask.any():
        send_error_json(handler, "selection matches no samples")
        return
    valid = mask & np.isfinite(z)
    z[valid] = _OPS[op](z[valid], value)

    # write result: overwrite or save as a new sample layer
    src_path = current.root / entry["path"]
    if save_as:
        out_name = f"mod_{save_as.replace(' ', '_')}.npz"
        out_path = current.rawdata_dir / out_name
    else:
        out_path = src_path
        if src_path.suffix == ".tif":
            # rewrite tif in place needs rasterio; convert to npz instead
            out_path = src_path.with_suffix(".npz")

    if kind == "raster":
        np.savez_compressed(out_path, x=np.asarray(x), y=np.asarray(y), z=z.astype("float32"))
        new_kind = "raster_nc"
    else:
        np.savez_compressed(out_path, x=np.asarray(x), y=np.asarray(y), z=z.astype("float32"))
        new_kind = "points"

    rel = f"gui/rawdata/{out_path.name}"
    if save_as:
        new_entry = {**entry, "kind": new_kind, "path": rel,
                     "label": save_as, "source": entry.get("source", "modified")}
        new_entry.pop("id", None)
        add_entries([new_entry])
    else:
        manifest = load_manifest()
        for item in manifest["entries"]:
            if item.get("id") == entry["id"]:
                item["kind"] = new_kind
                item["path"] = rel
        save_manifest(manifest)

    send_json(handler, {"ok": True, "cells": int(valid.sum()),
                        "min": float(np.nanmin(z)), "max": float(np.nanmax(z))})


# ---------------------------------------------------------------------
# band math on multi-channel rasters (e.g. NDVI) + threshold classify.
# Produces a new single-band raster dataset, e.g. for vegetation masks.
# ---------------------------------------------------------------------

def _read_all_bands(entry):
    """Every band of a raster entry as float32 arrays, plus x/y coords."""
    path = project.require().root / entry["path"]
    if entry.get("kind") == "raster_nc":
        data = np.load(path)
        return (np.asarray(data["x"]), np.asarray(data["y"]),
                [np.asarray(data["z"], dtype="float32")])
    import rasterio
    with rasterio.open(path) as ds:
        bands = [ds.read(i + 1).astype("float32") for i in range(ds.count)]
        nodata = ds.nodata
        if nodata is not None:
            for b in bands:
                b[b == np.float32(nodata)] = np.nan
        tr = ds.transform
        x = tr.c + tr.a * (np.arange(ds.width) + 0.5)
        y = tr.f + tr.e * (np.arange(ds.height) + 0.5)
    return x, y, bands


_EXPR_OK = re.compile(r"^[\s0-9.eE+\-*/()bB]+$")


def _eval_band_expr(expr, bands):
    """Evaluate an arithmetic expression over band arrays b1..bN. Only
    numbers and + - * / ( ) are allowed (regex-guarded, no builtins)."""
    expr = (expr or "").strip() or "b1"
    if not _EXPR_OK.match(expr):
        raise RuntimeError("expression may use only b1..bN, numbers and + - * / ( )")
    env = {f"b{i + 1}": bands[i] for i in range(len(bands))}
    with np.errstate(divide="ignore", invalid="ignore"):
        try:
            out = eval(expr, {"__builtins__": {}}, env)   # noqa: S307 - regex-restricted
        except Exception as exc:  # noqa: BLE001
            raise RuntimeError(f"bad expression: {exc}")
    return np.asarray(out, dtype="float32")


@route("POST", "/api/domain/sample_band")
def _sample_band(handler, body, tail):
    """Select which band of a multi-channel raster is displayed/used."""
    entry = get_entry(body.get("id", ""))
    if entry is None:
        send_error_json(handler, "unknown dataset", 404)
        return
    try:
        band = int(body.get("band"))
    except (TypeError, ValueError):
        send_error_json(handler, "missing integer 'band'")
        return
    band = max(1, min(band, int(entry.get("bands", 1))))
    manifest = load_manifest()
    for item in manifest["entries"]:
        if item.get("id") == entry["id"]:
            item["band"] = band
    save_manifest(manifest)
    send_json(handler, {"ok": True, "band": band})


@route("POST", "/api/domain/raster_derive")
def _raster_derive(handler, body, tail):
    """New single-band raster from a band expression (e.g. NDVI
    ``(b4-b1)/(b4+b1)``) with an optional threshold that classifies it
    (value cmp x -> then, else other). Saved as a new raw dataset."""
    current = project.require()
    entry = get_entry(body.get("id", ""))
    if entry is None:
        send_error_json(handler, "unknown dataset", 404)
        return
    if entry.get("kind") not in ("raster", "raster_nc"):
        send_error_json(handler, "band math needs a raster dataset")
        return
    save_as = (body.get("save_as") or "").strip() or "derived"
    try:
        x, y, bands = _read_all_bands(entry)
        z = _eval_band_expr(body.get("expr"), bands)
    except RuntimeError as exc:
        send_error_json(handler, exc)
        return

    thr = body.get("threshold")
    if thr:
        try:
            xcut = float(thr.get("x"))
            then_v = float(thr.get("then", 1))
            else_v = float(thr.get("else", 0))
        except (TypeError, ValueError):
            send_error_json(handler, "threshold needs numeric x / then / else")
            return
        cmp = z >= xcut if str(thr.get("op", ">")) in (">", ">=") else z <= xcut
        z = np.where(np.isfinite(z) & cmp, then_v, else_v).astype("float32")

    slug = re.sub(r"[^A-Za-z0-9._-]+", "_", save_as).strip("_") or "raster"
    out_name = f"derived_{slug}.npz"
    out_path = current.rawdata_dir / out_name
    np.savez_compressed(out_path, x=np.asarray(x), y=np.asarray(y), z=z)
    new_entry = {
        "source": "derived", "kind": "raster_nc",
        "path": f"gui/rawdata/{out_name}", "crs": entry.get("crs"),
        "res": entry.get("res"), "bounds": entry.get("bounds"),
        "label": save_as,
        "downloaded": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M"),
    }
    add_entries([new_entry])
    finite = z[np.isfinite(z)]
    send_json(handler, {"ok": True, "entry": new_entry,
                        "min": float(finite.min()) if finite.size else 0.0,
                        "max": float(finite.max()) if finite.size else 0.0})


@route("POST", "/api/domain/sample_duplicate")
def _sample_duplicate(handler, body, tail):
    import shutil
    current = project.require()
    entry = get_entry(body.get("id", ""))
    if entry is None:
        send_error_json(handler, "unknown sample layer", 404)
        return
    src = current.root / entry["path"]
    stem = src.stem
    n = 2
    while (current.rawdata_dir / f"{stem}_copy{n}{src.suffix}").exists():
        n += 1
    dst = current.rawdata_dir / f"{stem}_copy{n}{src.suffix}"
    shutil.copy2(src, dst)
    new_entry = {**entry, "path": f"gui/rawdata/{dst.name}",
                 "label": f"{entry.get('label', stem)} (copy)"}
    new_entry.pop("id", None)
    add_entries([new_entry])
    send_json(handler, {"ok": True})


@route("POST", "/api/domain/sample_rename")
def _sample_rename(handler, body, tail):
    """Rename a sample. Always updates the display ``label``; when
    ``rename_file`` is set, the underlying .npz in gui/rawdata is renamed to a
    slug of the new name (collisions suffixed), and the entry's ``path``/``id``
    are updated to match so the file and card stay in sync."""
    import re

    name = (body.get("name") or "").strip()
    if not name:
        send_error_json(handler, "missing 'name'")
        return
    rename_file = bool(body.get("rename_file"))
    current = project.require()
    manifest = load_manifest()
    updated = None
    for item in manifest["entries"]:
        if item.get("id") != body.get("id"):
            continue
        item["label"] = name
        # a linked/external entry only gets a new label - never move the
        # file it points at (it belongs to another project)
        if rename_file and not item.get("linked"):
            old_path = current.root / item.get("path", "")
            suffix = old_path.suffix or ".npz"
            slug = re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("_.") or "sample"
            base = f"{slug}{suffix}"
            target = current.rawdata_dir / base
            n = 2
            while target.exists() and target.resolve() != old_path.resolve():
                base = f"{slug}_{n}{suffix}"
                target = current.rawdata_dir / base
                n += 1
            if target.resolve() != old_path.resolve():
                try:
                    if old_path.is_file():
                        old_path.rename(target)
                except Exception as exc:  # noqa: BLE001
                    send_error_json(handler, f"rename failed: {exc}")
                    return
                item["path"] = f"gui/rawdata/{base}"
                item["id"] = hashlib.sha1(item["path"].encode()).hexdigest()[:10]
        updated = item
        break
    save_manifest(manifest)
    send_json(handler, {"ok": True, "entry": updated})


@route("POST", "/api/domain/sample_order")
def _sample_order(handler, body, tail):
    order = body.get("ids") or []
    manifest = load_manifest()
    by_id = {e.get("id"): e for e in manifest["entries"]}
    reordered = [by_id[i] for i in order if i in by_id]
    reordered += [e for e in manifest["entries"] if e.get("id") not in order]
    manifest["entries"] = reordered
    save_manifest(manifest)
    send_json(handler, {"ok": True})


@route("POST", "/api/domain/to_sample")
def _to_sample(handler, body, tail):
    """Convert an interpolated .grd into a (point) sample layer so it
    can be modified like any other sample data."""
    current = project.require()
    target = body.get("target")
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
    path = current.root / str(filename)
    if not path.is_file():
        send_error_json(handler, f"{filename} does not exist", 404)
        return
    Z = grd_io.read_grd(path)
    if np.iscomplexobj(Z):
        send_error_json(handler, "complex-valued masks cannot be converted to samples")
        return
    if Z.shape != X.shape:
        send_error_json(handler, f"{filename} does not match the current grid", 409)
        return
    # never overwrite an earlier conversion: bump a version suffix so the
    # user gets a fresh "(2)", "(3)", … layer each time (both the file on
    # disk and the manifest entry are keyed by the same unique path)
    known_paths = {e.get("path") for e in load_manifest()["entries"]}
    version = 1
    out_name = f"from_{target}.npz"
    while (f"gui/rawdata/{out_name}" in known_paths
           or (current.rawdata_dir / out_name).exists()):
        version += 1
        out_name = f"from_{target}_{version}.npz"
    label = f"{target} ({filename}) → samples"
    if version > 1:
        label += f" ({version})"
    np.savez_compressed(current.rawdata_dir / out_name,
                        x=X.ravel(), y=Y.ravel(), z=Z.ravel().astype("float32"))
    add_entries([{
        "source": "converted",
        "kind": "points",
        "path": f"gui/rawdata/{out_name}",
        "crs": None,
        "res": None,
        "bounds": [float(X.min()), float(Y.min()), float(X.max()), float(Y.max())],
        # remember the grid shape so the viewer can decimate the point
        # display per row/column (keeps a regular dot pattern)
        "shape": [int(X.shape[0]), int(X.shape[1])],
        "label": label,
        "downloaded": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M"),
    }])
    send_json(handler, {"ok": True, "label": label})


def _log_history(record):
    current = project.require()
    state = current.load_state()
    history = state.setdefault("domain_history", [])
    record["time"] = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M")
    history.append(record)
    state["domain_history"] = history[-100:]
    current.save_state(state)
