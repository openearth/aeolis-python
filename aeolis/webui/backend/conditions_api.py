"""Boundary conditions API (wind / water levels / waves).

Synthetic generation, measured data (waterinfo.rws.nl), and reanalysis
wind (ERA5) are all converted to the AeoLiS text formats:

    wind_file  columns [t, speed, direction]
    tide_file  columns [t, water level]
    wave_file  columns [t, Hs, Tp]

with t in seconds since the configuration ``refdate``. Raw downloads
stay cached in gui/rawdata.
"""

import hashlib
import os
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np

from aeolis.webui.backend import jobs, project
from aeolis.webui.backend.config_api import load_config
from aeolis.webui.backend.datasources import era5, synthetic, waterinfo
from aeolis.webui.backend.grid_api import _load_current_grid, patch_config, resolve_target
from aeolis.webui.backend.httpd import route
from aeolis.webui.backend.util import load_json, save_json, send_error_json, send_json

KINDS = {
    "wind": {"key": "wind_file", "default": "wind.txt", "cols": ["speed [m/s]", "direction [deg]"]},
    "tide": {"key": "tide_file", "default": "tide.txt", "cols": ["water level [m]"]},
    "wave": {"key": "wave_file", "default": "waves.txt", "cols": ["Hs [m]", "Tp [s]"]},
}

KIND_TITLES = {"wind": "Wind", "tide": "Water levels", "wave": "Waves"}

MAX_PREVIEW = 3000

# The six physical variables the Generate wizard works with. Each is a
# single column of its kind's file, generated on its own so a raw series can
# hold just that one quantity (they recombine later in Fill).
VARIABLE_META = {
    "wind_speed":  {"kind": "wind", "label": "speed [m/s]",     "clip0": True,  "wrap360": False},
    "wind_dir":    {"kind": "wind", "label": "direction [deg]", "clip0": False, "wrap360": True},
    "water_level": {"kind": "tide", "label": "water level [m]", "clip0": False, "wrap360": False},
    "wave_height": {"kind": "wave", "label": "Hs [m]",          "clip0": True,  "wrap360": False},
    "wave_period": {"kind": "wave", "label": "Tp [s]",          "clip0": True,  "wrap360": False},
    "wave_dir":    {"kind": "wave", "label": "direction [deg]", "clip0": False, "wrap360": True},
}


def parse_refdate(values):
    raw = str(values.get("refdate") or "2020-01-01 00:00")
    for fmt in ("%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M", "%Y-%m-%d"):
        try:
            return datetime.strptime(raw, fmt).replace(tzinfo=timezone.utc)
        except ValueError:
            continue
    raise ValueError(f"cannot parse refdate '{raw}'")


def _series_payload(data, refdate, tmin=None, tmax=None):
    """Decimated series + epoch times for the graphs. When a [tmin, tmax]
    epoch window is given, the series is restricted to it BEFORE decimating,
    so zooming in re-fetches denser detail (n/t0/t1 stay the full extent so
    the availability band still spans the whole record)."""
    data = np.atleast_2d(np.asarray(data, dtype=float))
    epoch0 = refdate.timestamp()
    t_all = epoch0 + data[:, 0]
    full_n = int(data.shape[0])
    full_t0 = float(t_all[0]) if full_n else float(epoch0)
    full_t1 = float(t_all[-1]) if full_n else float(epoch0)
    if tmin is not None or tmax is not None:
        lo = -np.inf if tmin is None else float(tmin)
        hi = np.inf if tmax is None else float(tmax)
        # keep one sample of padding each side so lines reach the edges
        idx = np.where((t_all >= lo) & (t_all <= hi))[0]
        if idx.size:
            a = max(0, idx[0] - 1)
            b = min(full_n, idx[-1] + 2)
            data = data[a:b]
    stride = max(1, data.shape[0] // MAX_PREVIEW)
    d = data[::stride]
    return {
        "t_epoch": (epoch0 + d[:, 0]).tolist(),
        # NaN (e.g. in a hand-edited wind.txt) is invalid JSON -> null
        "columns": [[float(v) if np.isfinite(v) else None for v in d[:, i]]
                    for i in range(1, d.shape[1])],
        "n": full_n,
        "t0_epoch": full_t0,
        "t1_epoch": full_t1,
    }


def _query_window(query):
    """Parse optional tmin/tmax epoch-second query params."""
    def _num(key):
        v = query.get(key)
        try:
            return float(v) if v not in (None, "") else None
        except (TypeError, ValueError):
            return None
    return _num("tmin"), _num("tmax")


@route("GET", "/api/conditions")
def _overview(handler, query, tail):
    current = project.require()
    values = load_config(current.configfile)
    refdate = parse_refdate(values)
    out = {
        "refdate": refdate.strftime("%Y-%m-%d %H:%M"),
        "refdate_epoch": refdate.timestamp(),
        "tstart": values.get("tstart"),
        "tstop": values.get("tstop"),
        "wind_convention": values.get("wind_convention"),
        "kinds": {},
    }
    for kind, info in KINDS.items():
        filename = values.get(info["key"])
        entry = {"config_key": info["key"], "file": filename, "exists": False,
                 "labels": info["cols"]}
        if filename:
            path = current.root / filename
            if path.is_file():
                entry["exists"] = True
                try:
                    data = np.atleast_2d(np.loadtxt(path))
                    entry["series"] = _series_payload(data, refdate)
                except (ValueError, OSError) as exc:
                    entry["error"] = str(exc)
        out["kinds"][kind] = entry
    send_json(handler, out)


def _synthetic_specs(kind, body):
    if kind == "wind":
        return [body.get("speed", {}), body.get("direction", {})]
    if kind == "tide":
        return [body.get("level", {})]
    return [body.get("hs", {}), body.get("tp", {})]


def _generate_synthetic(kind, body, values):
    tstart = float(body.get("tstart", values.get("tstart") or 0.0))
    tstop = float(body.get("tstop", values.get("tstop") or 3600.0))
    dt = float(body.get("dt", 3600.0))

    # when every segment of every quantity has an explicit duration the
    # series is generated over that time only - AeoLiS repeats a shorter
    # boundary-condition file cyclically over the simulation
    specs = _synthetic_specs(kind, body)
    totals = [synthetic.segments_duration(s) for s in specs]
    if totals and all(t is not None for t in totals):
        total = max(totals)
        if total > 0:
            tstop = min(tstop, tstart + total)

    if kind == "wind":
        return synthetic.wind(tstart, tstop, dt, body.get("speed", {}), body.get("direction", {}))
    if kind == "tide":
        return synthetic.tide(tstart, tstop, dt, body.get("level", {}))
    return synthetic.waves(tstart, tstop, dt, body.get("hs", {}), body.get("tp", {}))


def _generate_variable(variable, spec, values, body):
    """Generate a single-column [t, value] series for one physical variable
    from its segment spec. Mirrors the repeat-when-fully-timed behaviour of
    the kind-based generator."""
    meta = VARIABLE_META[variable]
    tstart = float(body.get("tstart", values.get("tstart") or 0.0))
    tstop = float(body.get("tstop", values.get("tstop") or 3600.0))
    dt = float(body.get("dt", 3600.0))
    total = synthetic.segments_duration(spec)
    if total is not None and total > 0:
        tstop = min(tstop, tstart + total)
    t = synthetic.time_axis(tstart, tstop, dt)
    col = synthetic.profile(t, spec)
    if meta["clip0"]:
        col = np.clip(col, 0.0, None)
    if meta["wrap360"]:
        col = np.mod(col, 360.0)
    return np.column_stack([t, col]), meta["label"]


def _tile_series(data, t_end):
    """Repeat a [t, cols...] series cyclically until *t_end*, the way
    AeoLiS wraps boundary conditions (interp_circular*). Returns the
    tiled array, or None when the series already covers t_end."""
    t = data[:, 0]
    period = float(t[-1] - t[0])
    if period <= 0 or t[-1] >= t_end - 1e-9:
        return None
    parts = [data]
    k = 1
    while t[0] + k * period < t_end and k < 4000:
        shifted = data[1:].copy()      # drop the duplicate first sample
        shifted[:, 0] += k * period
        parts.append(shifted)
        k += 1
    tiled = np.vstack(parts)
    return tiled[tiled[:, 0] <= t_end + 1e-9]


@route("POST", "/api/conditions/preview")
def _preview(handler, body, tail):
    """Synthetic series preview - computed only, nothing written. Accepts a
    single `variable` (one column) or the legacy `kind` (all columns)."""
    current = project.require()
    variable = body.get("variable")
    values = load_config(current.configfile)
    refdate = parse_refdate(values)
    try:
        if variable is not None:
            if variable not in VARIABLE_META:
                send_error_json(handler, f"unknown variable '{variable}'")
                return
            data, label = _generate_variable(variable, body.get("segments", {}), values, body)
            labels = [label]
        else:
            kind = body.get("kind")
            if kind not in KINDS:
                send_error_json(handler, f"unknown kind '{kind}'")
                return
            data = _generate_synthetic(kind, body, values)
            labels = KINDS[kind]["cols"]
    except ValueError as exc:
        send_error_json(handler, exc)
        return
    out = {"series": _series_payload(data, refdate), "labels": labels}
    # preview the cyclic repetition AeoLiS applies when the series is
    # shorter than the simulation
    tstop = float(values.get("tstop") or 0.0)
    tiled = _tile_series(data, tstop)
    if tiled is not None:
        out["series"] = _series_payload(tiled, refdate)
        out["repeat_from_epoch"] = float(refdate.timestamp() + data[-1, 0])
    send_json(handler, out)


@route("POST", "/api/conditions/synthetic")
def _synthetic(handler, body, tail):
    current = project.require()
    kind = body.get("kind")
    if kind not in KINDS:
        send_error_json(handler, f"unknown kind '{kind}'")
        return
    values = load_config(current.configfile)
    refdate = parse_refdate(values)

    try:
        data = _generate_synthetic(kind, body, values)
    except ValueError as exc:
        send_error_json(handler, exc)
        return

    info = KINDS[kind]
    filename = body.get("filename") or values.get(info["key"]) or info["default"]
    synthetic.write_series(current.root / filename, data)
    patch_config({info["key"]: filename})
    send_json(handler, {
        "ok": True, "file": filename, "rows": int(data.shape[0]),
        "series": _series_payload(data, refdate),
    })


# ---------------------------------------------------------------------
# measured / reanalysis sources
# ---------------------------------------------------------------------

def _grid_center_lonlat():
    grids, values = _load_current_grid()
    if grids is None:
        raise RuntimeError("no model grid - generate one in the Grid tab first")
    X, Y = grids
    x, y = float(np.mean(X)), float(np.mean(Y))
    state = project.require().load_state()
    crs = (state.get("crs") or {})
    if crs.get("mode") == "local":
        raise RuntimeError("project uses a local/conceptual coordinate system - "
                           "measured data sources need a georeferenced grid")
    epsg = int(crs.get("epsg") or 28992)
    if epsg == 4326:
        return x, y
    try:
        import pyproj
        tf = pyproj.Transformer.from_crs(epsg, 4326, always_xy=True)
        return tf.transform(x, y)
    except ImportError:
        if epsg == 28992:
            return _rd_to_wgs84(x, y)
        raise RuntimeError(
            f"pyproj is required to georeference EPSG:{epsg} - "
            "install with: pip install aeolis[webui-data]"
        )


def _rd_to_wgs84(x, y):
    """RD New -> WGS84 (polynomial approximation, ~0.3 m)."""
    dx = (x - 155000.0) * 1e-5
    dy = (y - 463000.0) * 1e-5
    lat = 52.15517440 + (
        3235.65389 * dy - 32.58297 * dx ** 2 - 0.24750 * dy ** 2
        - 0.84978 * dx ** 2 * dy - 0.06550 * dy ** 3 - 0.01709 * dx ** 2 * dy ** 2
        - 0.00738 * dx + 0.00530 * dx ** 4 - 0.00039 * dx ** 2 * dy ** 3
        + 0.00033 * dx ** 4 * dy - 0.00012 * dx * dy) / 3600.0
    lon = 5.38720621 + (
        5260.52916 * dx + 105.94684 * dx * dy + 2.45656 * dx * dy ** 2
        - 0.81885 * dx ** 3 + 0.05594 * dx * dy ** 3 - 0.05607 * dx ** 3 * dy
        + 0.01199 * dy - 0.00256 * dx ** 3 * dy ** 2 + 0.00128 * dx * dy ** 4
        + 0.00022 * dy ** 2 - 0.00022 * dx ** 2 + 0.00026 * dx ** 5) / 3600.0
    return lon, lat


@route("GET", "/api/conditions/stations")
def _stations(handler, query, tail):
    source = query.get("source", "waterinfo")
    kind = query.get("kind", "wind")
    try:
        lon, lat = _grid_center_lonlat()
    except RuntimeError as exc:
        send_error_json(handler, exc)
        return

    if source == "era5":
        ok, reason = era5.configured()
        send_json(handler, {
            "center": [lon, lat], "configured": ok, "reason": reason,
            "stations": era5.nearby_cells(lon, lat) if ok else era5.nearby_cells(lon, lat),
        })
        return

    def _run(job):
        job.update(message="loading waterinfo catalog")
        stations = waterinfo.stations_for(kind, lon=lon, lat=lat)
        return {"center": [lon, lat], "stations": stations}

    send_json(handler, {"job": jobs.start("waterinfo stations", _run)})


def _resample(seconds, cols, kind, interval):
    """Bin-average a series to a fixed interval. ``interval`` is a bin
    width in seconds (legacy "hour"/"day" strings still work). Wind
    direction is averaged circularly via unit-vector components.

    Bins are stamped at their LEFT edge, so the first output sample keeps
    the exact start time of the input (stamping at bin centres shifted
    the whole series by half an interval, and a series meant to start at
    tstart began at tstart + width/2 instead)."""
    width = {"hour": 3600.0, "day": 86400.0}.get(interval)
    if width is None:
        try:
            width = float(interval)
        except (TypeError, ValueError):
            width = 0.0
    if width <= 0 or seconds.size == 0:
        return seconds, cols
    bins = np.floor((seconds - seconds[0]) / width).astype(int)
    unique_bins, inverse = np.unique(bins, return_inverse=True)
    counts = np.bincount(inverse)
    t_out = seconds[0] + unique_bins * width

    if kind == "wind" and len(cols) >= 2:
        speed, direction = cols[0], cols[1]
        rad = np.deg2rad(direction)
        u = np.bincount(inverse, weights=speed * np.sin(rad)) / counts
        v = np.bincount(inverse, weights=speed * np.cos(rad)) / counts
        mean_speed = np.bincount(inverse, weights=speed) / counts
        mean_dir = np.mod(np.rad2deg(np.arctan2(u, v)), 360.0)
        out_cols = [mean_speed, mean_dir] + [
            np.bincount(inverse, weights=c) / counts for c in cols[2:]
        ]
    else:
        out_cols = [np.bincount(inverse, weights=c) / counts for c in cols]
    return t_out, out_cols


@route("POST", "/api/conditions/station_period")
def _station_period(handler, body, tail):
    station = body.get("station")
    kind = body.get("kind")
    if not station or kind not in KINDS:
        send_error_json(handler, "missing 'station'/'kind'")
        return

    def _run(job):
        return waterinfo.probe_period(station, kind, job)

    send_json(handler, {"job": jobs.start("probe data period", _run)})


# ---------------------------------------------------------------------
# CDS (ERA5) API key management - the key is personal and stored
# locally in ~/.cdsapirc, never in the project or repository.
# ---------------------------------------------------------------------

@route("GET", "/api/conditions/cds")
def _cds_status(handler, query, tail):
    ok, reason = era5.configured()
    send_json(handler, {"configured": ok, "reason": reason})


@route("POST", "/api/conditions/cds_key")
def _cds_key(handler, body, tail):
    from pathlib import Path
    key = (body.get("key") or "").strip()
    if not key or len(key) < 10:
        send_error_json(handler, "paste the full API token from your CDS profile page")
        return
    url = (body.get("url") or "https://cds.climate.copernicus.eu/api").strip()
    rc = Path.home() / ".cdsapirc"
    rc.write_text(f"url: {url}\nkey: {key}\n", encoding="utf-8")
    ok, reason = era5.configured()
    send_json(handler, {"ok": ok, "reason": reason, "path": str(rc)})


@route("POST", "/api/conditions/era5_cached")
def _era5_cached(handler, body, tail):
    """Which years of the requested period are already cached on disk for
    an ERA5 cell — shown in the download wizard before starting."""
    current = project.require()
    try:
        lon = float(body.get("lon"))
        lat = float(body.get("lat"))
        date0 = datetime.fromisoformat(str(body.get("date0"))).replace(tzinfo=timezone.utc)
        date1 = datetime.fromisoformat(str(body.get("date1"))).replace(tzinfo=timezone.utc)
    except (TypeError, ValueError) as exc:
        send_error_json(handler, f"missing/invalid lon, lat or dates: {exc}")
        return
    send_json(handler, era5.cached_summary(lon, lat, date0, date1, current.rawdata_dir))


@route("POST", "/api/conditions/era5_cells_cached")
def _era5_cells_cached(handler, body, tail):
    """Cached-year counts for a LIST of ERA5 cells — lets the wizard mark
    (and auto-select) the cell that was downloaded before."""
    current = project.require()
    try:
        date0 = datetime.fromisoformat(str(body.get("date0"))).replace(tzinfo=timezone.utc)
        date1 = datetime.fromisoformat(str(body.get("date1"))).replace(tzinfo=timezone.utc)
    except (TypeError, ValueError) as exc:
        send_error_json(handler, f"invalid dates: {exc}")
        return
    out = []
    for cell in (body.get("cells") or []):
        try:
            lon, lat = float(cell.get("lon")), float(cell.get("lat"))
        except (TypeError, ValueError):
            continue
        s = era5.cached_summary(lon, lat, date0, date1, current.rawdata_dir)
        out.append({"id": cell.get("id"), "lon": lon, "lat": lat,
                    "cached": len(s["cached"]), "total": s["total"]})
    send_json(handler, {"cells": out})


@route("POST", "/api/conditions/cds_clear")
def _cds_clear(handler, body, tail):
    from pathlib import Path
    rc = Path.home() / ".cdsapirc"
    if rc.is_file():
        rc.unlink()
    send_json(handler, {"ok": True})


# ---------------------------------------------------------------------
# raw time-series store: downloads land here first (gui/rawdata/*.npz +
# conditions_manifest.json), can be inspected/modified, and are only
# written to wind.txt/tide.txt/waves.txt when the user applies them.
# ---------------------------------------------------------------------

def _cond_manifest_path():
    return project.require().rawdata_dir / "conditions_manifest.json"


def _load_cond_manifest():
    return load_json(_cond_manifest_path(), default={"entries": []})


def _save_cond_manifest(manifest):
    save_json(_cond_manifest_path(), manifest)


def _get_raw_entry(entry_id):
    for entry in _load_cond_manifest()["entries"]:
        if entry.get("id") == entry_id:
            return entry
    return None


def _load_raw_series(entry):
    """Returns (t_epoch float64 [n], cols float64 [n, m])."""
    data = np.load(project.require().root / entry["path"])
    t = np.asarray(data["t"], dtype="float64")
    cols = np.atleast_2d(np.asarray(data["cols"], dtype="float64"))
    if cols.shape[0] != t.size:
        cols = cols.T
    return t, cols


def _save_raw_entry(kind, source, label, t, cols, entry=None):
    """Write/overwrite an npz + manifest entry; returns the entry."""
    current = project.require()
    cols = np.column_stack(cols) if isinstance(cols, (list, tuple)) else np.atleast_2d(cols)
    if cols.shape[0] != t.size:
        cols = cols.T
    manifest = _load_cond_manifest()
    if entry is None:
        stamp = hashlib.sha1(f"{kind}|{source}|{label}|{t[0] if t.size else 0}|{len(manifest['entries'])}"
                             .encode()).hexdigest()[:10]
        path = current.rawdata_dir / f"cond_{kind}_{stamp}.npz"
        entry = {
            "id": stamp, "kind": kind, "source": source, "label": label,
            "path": f"gui/rawdata/{path.name}",
            "labels": KINDS[kind]["cols"][:cols.shape[1]],
            "created": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M"),
        }
        manifest["entries"].append(entry)
    else:
        path = current.root / entry["path"]
    np.savez_compressed(path, t=t.astype("float64"), cols=cols.astype("float64"))
    entry["rows"] = int(t.size)
    entry["nan"] = int(np.sum(~np.isfinite(cols)))
    entry["t0_epoch"] = float(t[0]) if t.size else None
    entry["t1_epoch"] = float(t[-1]) if t.size else None
    for item in manifest["entries"]:
        if item.get("id") == entry["id"]:
            item.update(entry)
    _save_cond_manifest(manifest)
    return entry


@route("GET", "/api/conditions/raw")
def _raw_list(handler, query, tail):
    project.require()
    entries = _load_cond_manifest()["entries"]
    kind = query.get("kind")
    if kind:
        entries = [e for e in entries if e.get("kind") == kind]
    send_json(handler, {"entries": entries})


@route("GET", "/api/conditions/series")
def _input_series(handler, query, tail):
    """A single input file's series, optionally windowed (for zoom-aware
    higher-resolution graph data). Returns the same shape as overview series."""
    current = project.require()
    kind = query.get("kind")
    if kind not in KINDS:
        send_error_json(handler, f"unknown kind '{kind}'")
        return
    values = load_config(current.configfile)
    refdate = parse_refdate(values)
    filename = values.get(KINDS[kind]["key"])
    path = current.root / str(filename) if filename else None
    if not path or not path.is_file():
        send_error_json(handler, f"{kind} file not written yet", 404)
        return
    tmin, tmax = _query_window(query)
    data = np.atleast_2d(np.loadtxt(path))
    send_json(handler, {"series": _series_payload(data, refdate, tmin, tmax),
                        "labels": KINDS[kind]["cols"]})


@route("GET", "/api/conditions/raw_series")
def _raw_series(handler, query, tail):
    entry = _get_raw_entry(query.get("id", ""))
    if entry is None:
        send_error_json(handler, "unknown raw series", 404)
        return
    t, cols = _load_raw_series(entry)
    full_n = int(t.size)
    full_t0 = float(t[0]) if t.size else None
    full_t1 = float(t[-1]) if t.size else None
    # optional zoom window: decimate within it for denser detail on zoom-in
    tmin, tmax = _query_window(query)
    tw, cw = t, cols
    if (tmin is not None or tmax is not None) and t.size:
        lo = -np.inf if tmin is None else tmin
        hi = np.inf if tmax is None else tmax
        idx = np.where((t >= lo) & (t <= hi))[0]
        if idx.size:
            a = max(0, idx[0] - 1)
            b = min(t.size, idx[-1] + 2)
            tw, cw = t[a:b], cols[a:b]
    stride = max(1, tw.size // MAX_PREVIEW)
    payload = {
        "t_epoch": tw[::stride].tolist(),
        "columns": [[float(v) if np.isfinite(v) else None for v in cw[::stride, i]]
                    for i in range(cw.shape[1])],
        "labels": entry.get("labels") or KINDS[entry["kind"]]["cols"][:cols.shape[1]],
        "n": full_n,
        "t0_epoch": full_t0,
        "t1_epoch": full_t1,
        "nan": int(np.sum(~np.isfinite(cols))),
    }
    send_json(handler, {"entry": entry, "series": payload})


_RAW_OPS = {"add", "multiply", "clip_max", "clip_min", "fill_nan"}


@route("POST", "/api/conditions/raw_modify")
def _raw_modify(handler, body, tail):
    entry = _get_raw_entry(body.get("id", ""))
    if entry is None:
        send_error_json(handler, "unknown raw series", 404)
        return
    op = body.get("op")
    t, cols = _load_raw_series(entry)
    column = body.get("column")     # None = all columns
    col_idx = range(cols.shape[1]) if column in (None, "", "all") else [int(column)]

    try:
        if op == "drop_nan":
            keep = np.all(np.isfinite(cols), axis=1)
            t, cols = t[keep], cols[keep]
        elif op == "crop":
            d0 = datetime.fromisoformat(body["date0"]).replace(tzinfo=timezone.utc).timestamp()
            d1 = datetime.fromisoformat(body["date1"]).replace(tzinfo=timezone.utc).timestamp()
            keep = (t >= d0) & (t <= d1 + 86399)
            t, cols = t[keep], cols[keep]
        elif op == "resample":
            width = float(body.get("value") or 3600)
            t, out = _resample(t, [cols[:, i] for i in range(cols.shape[1])],
                               entry["kind"], width)
            cols = np.column_stack(out)
        elif op == "fill_from":
            other = _get_raw_entry(body.get("other", ""))
            if other is None:
                raise ValueError("select another raw series to fill from")
            to, ocols = _load_raw_series(other)
            m = min(cols.shape[1], ocols.shape[1])
            for i in range(m):
                good = np.isfinite(ocols[:, i])
                if not good.any():
                    continue
                bad = ~np.isfinite(cols[:, i])
                cols[bad, i] = np.interp(t[bad], to[good], ocols[good, i])
        elif op == "fill_nan":
            value = float(body.get("value"))
            for i in col_idx:
                col = cols[:, i]
                col[~np.isfinite(col)] = value
        elif op in _RAW_OPS or op in ("set", "subtract"):
            value = float(body.get("value"))
            fns = {
                "set": lambda a, v: np.full_like(a, v),
                "add": lambda a, v: a + v,
                "subtract": lambda a, v: a - v,
                "multiply": lambda a, v: a * v,
                "clip_max": lambda a, v: np.minimum(a, v),
                "clip_min": lambda a, v: np.maximum(a, v),
            }
            for i in col_idx:
                good = np.isfinite(cols[:, i])
                cols[good, i] = fns[op](cols[good, i], value)
        else:
            send_error_json(handler, f"unknown op '{op}'")
            return
    except (KeyError, TypeError, ValueError) as exc:
        send_error_json(handler, exc)
        return

    if not t.size:
        send_error_json(handler, "the operation removed every sample - nothing saved")
        return

    save_as = (body.get("save_as") or "").strip()
    if save_as:
        entry = _save_raw_entry(entry["kind"], "modified", save_as, t, cols)
    else:
        entry = _save_raw_entry(entry["kind"], entry.get("source"), entry.get("label"),
                                t, cols, entry=entry)
    send_json(handler, {"ok": True, "entry": entry})


# ---------------------------------------------------------------------
# data-flaw detection & clean-up: sensors report flaws as sentinel
# values (999 m waves, -999) or as long constant stretches (0 m/s wind,
# 72 deg direction for days). Detect them and replace / interpolate /
# remove.
# ---------------------------------------------------------------------

def _flag_flaws(col, rules):
    """Boolean mask of flawed samples in one column + per-rule counts."""
    flags = np.zeros(col.size, dtype=bool)
    counts = {}
    for s in (rules.get("sentinels") or []):
        m = np.isclose(col, float(s), rtol=0.0, atol=1e-9)
        counts[f"= {s}"] = int(m.sum())
        flags |= m
    run_min = rules.get("run_min")
    if run_min and col.size >= 2:
        finite = np.isfinite(col)
        change = np.ones(col.size, dtype=bool)
        change[1:] = ~((col[1:] == col[:-1]) & finite[1:] & finite[:-1])
        run_id = np.cumsum(change) - 1
        lengths = np.bincount(run_id)[run_id]
        m = (lengths >= int(run_min)) & finite
        rv = rules.get("run_value")
        if rv not in (None, ""):
            m &= np.isclose(col, float(rv), rtol=0.0, atol=1e-9)
        counts[f"constant ≥ {int(run_min)} steps"] = int(m.sum())
        flags |= m
    return flags, counts


@route("POST", "/api/conditions/raw_clean")
def _raw_clean(handler, body, tail):
    """Detect flawed stretches in a raw series and clean them up.

    body: {id, column: "all"|idx, rules: {sentinels: [999, ...],
    run_min: N, run_value: v|null}, action: "nan"|"interp"|"remove",
    preview: bool, save_as: str}
    """
    entry = _get_raw_entry(body.get("id", ""))
    if entry is None:
        send_error_json(handler, "unknown raw series", 404)
        return
    rules = body.get("rules") or {}
    action = body.get("action") or "nan"
    t, cols = _load_raw_series(entry)
    labels = entry.get("labels") or KINDS[entry["kind"]]["cols"][:cols.shape[1]]

    column = body.get("column")
    col_idx = list(range(cols.shape[1])) if column in (None, "", "all") else [int(column)]

    report = []
    masks = np.zeros(cols.shape, dtype=bool)
    for i in col_idx:
        flags, counts = _flag_flaws(cols[:, i], rules)
        masks[:, i] = flags
        report.append({"column": labels[i] if i < len(labels) else f"col {i}",
                       "flagged": int(flags.sum()), "rules": counts})
    total = int(masks.sum())

    if body.get("preview"):
        send_json(handler, {"ok": True, "preview": True, "total": total,
                            "rows": int(t.size), "report": report})
        return
    if total == 0:
        send_error_json(handler, "nothing detected with these rules - nothing to clean")
        return

    if action == "remove":
        keep = ~masks.any(axis=1)
        t, cols = t[keep], cols[keep]
        if not t.size:
            send_error_json(handler, "cleaning removed every sample - nothing saved")
            return
    else:
        cols = np.array(cols, copy=True)
        cols[masks] = np.nan
        if action == "interp":
            for i in col_idx:
                col = cols[:, i]
                good = np.isfinite(col)
                if good.sum() >= 2:
                    bad = ~good
                    col[bad] = np.interp(t[bad], t[good], col[good])

    save_as = (body.get("save_as") or "").strip()
    if save_as:
        entry = _save_raw_entry(entry["kind"], "cleaned", save_as, t, cols)
    else:
        entry = _save_raw_entry(entry["kind"], entry.get("source"), entry.get("label"),
                                t, cols, entry=entry)
    send_json(handler, {"ok": True, "entry": entry, "total": total, "report": report})


@route("POST", "/api/conditions/raw_rename")
def _raw_rename(handler, body, tail):
    name = (body.get("name") or "").strip()
    if not name:
        send_error_json(handler, "missing 'name'")
        return
    manifest = _load_cond_manifest()
    for item in manifest["entries"]:
        if item.get("id") == body.get("id"):
            item["label"] = name
    _save_cond_manifest(manifest)
    send_json(handler, {"ok": True})


@route("POST", "/api/conditions/raw_duplicate")
def _raw_duplicate(handler, body, tail):
    """Copy a raw series to a new independent entry (same workflow as the
    Domain sample duplicate) - handy when several series share most of
    their processing."""
    entry = _get_raw_entry(body.get("id", ""))
    if entry is None:
        send_error_json(handler, "unknown raw series", 404)
        return
    t, cols = _load_raw_series(entry)
    label = f"{entry.get('label', 'series')} (copy)"
    new_entry = _save_raw_entry(
        entry["kind"], entry.get("source", "copy"), label,
        t, [cols[:, i] for i in range(cols.shape[1])])
    # carry the column labels of the source (may differ from the defaults)
    if entry.get("labels"):
        manifest = _load_cond_manifest()
        for item in manifest["entries"]:
            if item.get("id") == new_entry["id"]:
                item["labels"] = entry["labels"]
        _save_cond_manifest(manifest)
        new_entry["labels"] = entry["labels"]
    send_json(handler, {"entry": new_entry})


@route("POST", "/api/conditions/raw_delete")
def _raw_delete(handler, body, tail):
    manifest = _load_cond_manifest()
    entry = next((e for e in manifest["entries"] if e.get("id") == body.get("id")), None)
    if entry is None:
        send_error_json(handler, "unknown raw series", 404)
        return
    manifest["entries"] = [e for e in manifest["entries"] if e.get("id") != entry["id"]]
    _save_cond_manifest(manifest)
    target = project.require().root / entry["path"]
    if target.is_file():
        target.unlink()
    send_json(handler, {"ok": True})


@route("POST", "/api/conditions/raw_apply")
def _raw_apply(handler, body, tail):
    """Convert a raw series to the AeoLiS text format (seconds since
    refdate) and write it as the configured wind/tide/wave file."""
    current = project.require()
    entry = _get_raw_entry(body.get("id", ""))
    if entry is None:
        send_error_json(handler, "unknown raw series", 404)
        return
    kind = entry["kind"]
    values = load_config(current.configfile)
    refdate = parse_refdate(values)
    t, cols = _load_raw_series(entry)
    n_bad = int(np.sum(~np.isfinite(cols)))
    if n_bad:
        send_error_json(
            handler,
            f"series contains {n_bad} NaN value(s) - clean them first "
            "(Modify: drop NaN rows, fill with a value, or fill from another series)")
        return
    info = KINDS[kind]
    filename = body.get("filename") or values.get(info["key"]) or info["default"]
    seconds = t - refdate.timestamp()
    data = np.column_stack([seconds, cols])
    synthetic.write_series(current.root / filename, data)
    patch_config({info["key"]: filename})
    send_json(handler, {
        "ok": True, "file": filename, "rows": int(data.shape[0]),
        "series": _series_payload(data, refdate),
    })


@route("POST", "/api/conditions/synthetic_raw")
def _synthetic_raw(handler, body, tail):
    """Generate a synthetic series and store it as a RAW series (so it joins
    downloads in the raw-data section). Accepts a single `variable` (one
    column) or the legacy `kind` (all columns)."""
    current = project.require()
    variable = body.get("variable")
    values = load_config(current.configfile)
    refdate = parse_refdate(values)
    try:
        if variable is not None:
            if variable not in VARIABLE_META:
                send_error_json(handler, f"unknown variable '{variable}'")
                return
            data, col_label = _generate_variable(variable, body.get("segments", {}), values, body)
            kind = VARIABLE_META[variable]["kind"]
            labels = [col_label]
        else:
            kind = body.get("kind")
            if kind not in KINDS:
                send_error_json(handler, f"unknown kind '{kind}'")
                return
            data = _generate_synthetic(kind, body, values)
            labels = None
    except ValueError as exc:
        send_error_json(handler, exc)
        return
    t = refdate.timestamp() + data[:, 0]
    cols = [data[:, i] for i in range(1, data.shape[1])]
    default_label = f"synthetic {variable}" if variable else f"synthetic {kind}"
    label = (body.get("label") or "").strip() or default_label
    entry = _save_raw_entry(kind, "synthetic", label,
                            np.asarray(t, dtype="float64"), cols)
    # a single-variable series carries its own column label (not the kind's
    # default first column, which _save_raw_entry would otherwise assign)
    if labels is not None:
        manifest = _load_cond_manifest()
        for item in manifest["entries"]:
            if item.get("id") == entry["id"]:
                item["labels"] = labels
        _save_cond_manifest(manifest)
        entry["labels"] = labels
    send_json(handler, {"ok": True, "entry": entry})


def _nearest_fill(master, mask_good, values):
    """Nearest-neighbour hold for the NaN samples of *values* on time base
    *master*, using only the samples flagged good in *mask_good*."""
    out = values.copy()
    good_idx = np.where(mask_good)[0]
    bad_idx = np.where(~mask_good)[0]
    if good_idx.size == 0 or bad_idx.size == 0:
        return out
    mt = master[good_idx]
    pos = np.searchsorted(mt, master[bad_idx])
    pos = np.clip(pos, 1, mt.size - 1)
    left = pos - 1
    choose_left = (master[bad_idx] - mt[left]) <= (mt[pos] - master[bad_idx])
    src = np.where(choose_left, good_idx[left], good_idx[pos])
    if mt.size == 1:
        src = np.full(bad_idx.shape, good_idx[0])
    out[bad_idx] = values[src]
    return out


def _paint_from_source(master, out, t, v):
    """Fill still-NaN samples of *out* (on time base *master*) from source
    (t, v) — but ONLY where the source genuinely has data: a master sample
    is claimed when its bracketing source samples are at most ~3 sampling
    steps apart. Bigger holes (data gaps, sensor outage) stay NaN so a
    lower-priority source or the gap-fill method can take them; the old
    blanket np.interp silently drew straight lines across such gaps and
    lower sources never got a turn."""
    good = np.isfinite(v)
    if not good.any():
        return
    tg, vg = t[good], v[good]
    need = (~np.isfinite(out)) & (master >= tg[0]) & (master <= tg[-1])
    if not need.any():
        return
    tm = master[need]
    if tg.size < 2:
        ok = tm == tg[0]
    else:
        gap_limit = 3.0 * float(np.median(np.diff(tg)))
        idx = np.clip(np.searchsorted(tg, tm, side="right"), 1, tg.size - 1)
        ok = (tg[idx] - tg[idx - 1]) <= gap_limit
        ok |= (tm == tg[idx]) | (tm == tg[idx - 1])   # exact hits always count
    if ok.any():
        sel = np.flatnonzero(need)[ok]
        out[sel] = np.interp(tm[ok], tg, vg)


def _column_from_sources(master, series_list, fill, cache_loader):
    """Build one output column on time base *master* by layering *series_list*
    (priority order, top first — lower ones only fill samples still NaN), then
    apply the remaining-NaN *fill* policy (value / nearest / linear / series)."""
    out = np.full(master.shape, np.nan, dtype=float)
    for t, v in series_list:
        _paint_from_source(master, out, t, v)

    bad = ~np.isfinite(out)
    if bad.any():
        method = (fill or {}).get("method")
        if method == "value":
            out[bad] = float((fill or {}).get("value", 0.0))
        elif method == "nearest":
            out = _nearest_fill(master, ~bad, out)
        elif method == "linear":
            good = ~bad
            if good.any():
                out[bad] = np.interp(master[bad], master[good], out[good])
        elif method == "series":
            other = (fill or {}).get("other") or {}
            ot, ocols = cache_loader(other.get("id", ""))
            oc = int(other.get("column") or 0)
            if oc >= ocols.shape[1]:
                oc = 0
            _paint_from_source(master, out, ot, ocols[:, oc])
    return out


@route("POST", "/api/conditions/fill")
def _fill(handler, body, tail):
    """Assemble a wind/tide/wave input file. Each output column is built from a
    PRIORITY-ORDERED list of raw-series columns (top wins; lower ones only fill
    samples still NaN), then a chosen remaining-NaN method fills the rest. The
    legacy one-source-per-column `sources` body is still accepted."""
    current = project.require()
    kind = body.get("kind")
    if kind not in KINDS:
        send_error_json(handler, f"unknown kind '{kind}'")
        return
    ncol = len(KINDS[kind]["cols"])

    columns_spec = body.get("columns")
    if columns_spec is None:
        # legacy: one source per column, no NaN-fill policy
        sources = body.get("sources") or []
        columns_spec = [{"sources": [s], "fill": {}} for s in sources]
    if len(columns_spec) != ncol:
        send_error_json(handler, f"{kind} needs {ncol} column(s), got {len(columns_spec)}")
        return

    values = load_config(current.configfile)
    refdate = parse_refdate(values)

    _cache = {}

    def _load_cached(entry_id):
        if entry_id not in _cache:
            entry = _get_raw_entry(entry_id)
            if entry is None:
                raise ValueError("unknown raw series in selection")
            _cache[entry_id] = _load_raw_series(entry)
        return _cache[entry_id]

    def _column_series(spec):
        out = []
        for s in (spec.get("sources") or []):
            t, cols = _load_cached(s.get("id", ""))
            ci = int(s.get("column") or 0)
            if ci >= cols.shape[1]:
                ci = 0
            if t.size:
                out.append((t, cols[:, ci]))
        return out

    try:
        # Master time base: the UNION of every selected source's timestamps,
        # clipped to the span of the top-priority source of the first column.
        # Using only the primary's own timestamps (the old behaviour) meant a
        # data gap in the primary — a stretch with no samples at all — could
        # never be filled by a lower-priority station: those moments simply
        # did not exist in the output.
        first_series = _column_series(columns_spec[0])
        if not first_series:
            send_error_json(handler, "select at least one source for each column")
            return
        primary = first_series[0][0]
        if primary.size < 2:
            send_error_json(handler, "the primary series is too short")
            return
        all_t = [s[0] for spec in columns_spec for s in _column_series(spec)]
        master = np.unique(np.concatenate(all_t))
        master = master[(master >= primary[0]) & (master <= primary[-1])]

        out_cols = []
        for spec in columns_spec:
            series_list = _column_series(spec)
            if not series_list:
                send_error_json(handler, "select at least one source for each column")
                return
            out_cols.append(_column_from_sources(master, series_list, spec.get("fill"), _load_cached))
    except ValueError as exc:
        send_error_json(handler, exc, 404)
        return

    resample = body.get("resample")
    if resample:
        master, out_cols = _resample(master, out_cols, kind, resample)

    n_bad = int(np.sum(~np.isfinite(np.column_stack(out_cols))))
    if n_bad:
        send_error_json(handler, f"result still has {n_bad} NaN value(s) - add a lower-priority "
                        "source or choose a fill method for the gaps")
        return

    seconds = master - refdate.timestamp()
    data = np.column_stack([seconds] + list(out_cols))
    info = KINDS[kind]
    filename = body.get("filename") or values.get(info["key"]) or info["default"]
    synthetic.write_series(current.root / filename, data)
    patch_config({info["key"]: filename})
    send_json(handler, {
        "ok": True, "file": filename, "rows": int(data.shape[0]),
        "series": _series_payload(data, refdate),
    })


@route("POST", "/api/conditions/save_file_as")
def _save_file_as(handler, body, tail):
    """Save the current input file under a new name and repoint the config."""
    current = project.require()
    kind = body.get("kind")
    # accept a full path from the file browser, or a legacy bare filename
    raw = (body.get("path") or body.get("filename") or "").strip()
    if kind not in KINDS:
        send_error_json(handler, f"unknown kind '{kind}'")
        return
    if not raw:
        send_error_json(handler, "missing path")
        return
    info = KINDS[kind]
    values = load_config(current.configfile)
    src_name = values.get(info["key"]) or info["default"]
    src = current.root / str(src_name)
    if not src.is_file():
        send_error_json(handler, f"{src_name} does not exist yet - fill or generate it first")
        return
    write_path, config_ref = resolve_target(current, raw, str(info["default"]))
    write_path.parent.mkdir(parents=True, exist_ok=True)
    data = np.atleast_2d(np.loadtxt(src))
    synthetic.write_series(write_path, data)
    patch_config({info["key"]: config_ref})
    send_json(handler, {"ok": True, "file": config_ref})


@route("POST", "/api/conditions/load_file")
def _load_file(handler, body, tail):
    """Import an existing condition file from disk (copy + repoint config)."""
    current = project.require()
    kind = body.get("kind")
    path = body.get("path")
    if kind not in KINDS:
        send_error_json(handler, f"unknown kind '{kind}'")
        return
    if not path:
        send_error_json(handler, "missing path")
        return
    try:
        data = np.atleast_2d(np.loadtxt(Path(path)))
    except Exception as exc:  # noqa: BLE001 - report any read failure to the UI
        send_error_json(handler, f"could not read {path}: {exc}")
        return
    filename = os.path.basename(path)
    synthetic.write_series(current.root / filename, data)
    info = KINDS[kind]
    patch_config({info["key"]: filename})
    send_json(handler, {"ok": True, "file": filename, "rows": int(data.shape[0])})


@route("POST", "/api/conditions/fetch")
def _fetch(handler, body, tail):
    current = project.require()
    source = body.get("source")
    # accept several quantities (kinds) to download from one station in one
    # action; fall back to the legacy single `kind`
    kinds = body.get("kinds") or ([body.get("kind")] if body.get("kind") else [])
    kinds = [k for k in kinds if k in KINDS]
    if not kinds:
        send_error_json(handler, "select at least one quantity to download")
        return
    values = load_config(current.configfile)
    refdate = parse_refdate(values)
    tstart = float(values.get("tstart") or 0.0)
    tstop = float(values.get("tstop") or 3600.0)
    date0 = refdate + timedelta(seconds=tstart)
    date1 = refdate + timedelta(seconds=tstop)
    if body.get("date0"):
        date0 = datetime.fromisoformat(body["date0"]).replace(tzinfo=timezone.utc)
    if body.get("date1"):
        date1 = datetime.fromisoformat(body["date1"]).replace(tzinfo=timezone.utc)

    warnings = []   # non-fatal notes (partial downloads etc.) shown in the UI

    def _fetch_one(kind, job):
        if source == "era5":
            if kind != "wind":
                raise RuntimeError("ERA5 source currently provides wind only")
            cell = body.get("station") or {}
            t, u, v, cache = era5.download_wind(
                float(cell.get("lon")), float(cell.get("lat")),
                date0, date1, current.rawdata_dir, job,
            )
            speed, direction = era5.to_speed_direction(
                u, v, values.get("wind_convention") or "nautical")
            epoch = t.astype("float64")
            cols = [speed, direction]
            station_name = (cell.get("name") or "ERA5")
        elif source == "waterinfo":
            station = body.get("station")
            series, warns = waterinfo.fetch_series(station, kind, date0, date1, job)
            for w in warns:
                warnings.append(f"{KIND_TITLES.get(kind, kind)}: {w}")
            epoch, base = series[0]
            cols = [base]
            for epoch_i, values_i in series[1:]:
                if len(values_i) == 0:
                    # secondary quantity absent/failed: keep the column as NaN
                    cols.append(np.full(epoch.shape, np.nan))
                    continue
                # nearest-neighbour align extra quantities on the first one
                idx = np.searchsorted(epoch_i, epoch).clip(0, len(values_i) - 1)
                cols.append(values_i[idx])
            station_name = body.get("station_name") or str(station)
        else:
            raise RuntimeError(f"unknown source '{source}'")

        resample = body.get("resample")
        if resample:
            job.update(message="resampling to interval means")
            epoch, cols = _resample(np.asarray(epoch, dtype="float64"), cols, kind, resample)

        label = body.get("label") if len(kinds) == 1 and body.get("label") else \
            f"{station_name} {KIND_TITLES.get(kind, kind)} {date0:%Y-%m-%d} — {date1:%Y-%m-%d}"
        return _save_raw_entry(kind, source, label,
                               np.asarray(epoch, dtype="float64"), cols)

    def _run(job):
        entries, errors = [], []
        for i, kind in enumerate(kinds):
            job.update(message=f"downloading {KIND_TITLES.get(kind, kind)} "
                               f"({i + 1}/{len(kinds)})")
            try:
                entries.append(_fetch_one(kind, job))
            except Exception as exc:  # noqa: BLE001 - one quantity may be absent
                errors.append(f"{KIND_TITLES.get(kind, kind)}: {exc}")
        if not entries:
            raise RuntimeError("; ".join(errors) or "nothing downloaded")
        # partial-download warnings surface next to real errors in the UI
        return {"entries": entries, "errors": errors + warnings}

    send_json(handler, {"job": jobs.start(f"fetch {source} {'+'.join(kinds)}", _run)})
