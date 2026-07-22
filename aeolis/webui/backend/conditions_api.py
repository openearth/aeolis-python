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
from datetime import datetime, timedelta, timezone

import numpy as np

from aeolis.webui.backend import jobs, project
from aeolis.webui.backend.config_api import load_config
from aeolis.webui.backend.datasources import era5, synthetic, waterinfo
from aeolis.webui.backend.grid_api import _load_current_grid, patch_config
from aeolis.webui.backend.httpd import route
from aeolis.webui.backend.util import load_json, save_json, send_error_json, send_json

KINDS = {
    "wind": {"key": "wind_file", "default": "wind.txt", "cols": ["speed [m/s]", "direction [deg]"]},
    "tide": {"key": "tide_file", "default": "tide.txt", "cols": ["water level [m]"]},
    "wave": {"key": "wave_file", "default": "waves.txt", "cols": ["Hs [m]", "Tp [s]"]},
}

MAX_PREVIEW = 3000


def parse_refdate(values):
    raw = str(values.get("refdate") or "2020-01-01 00:00")
    for fmt in ("%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M", "%Y-%m-%d"):
        try:
            return datetime.strptime(raw, fmt).replace(tzinfo=timezone.utc)
        except ValueError:
            continue
    raise ValueError(f"cannot parse refdate '{raw}'")


def _series_payload(data, refdate):
    """Decimated series + epoch times for the graphs."""
    data = np.atleast_2d(np.asarray(data, dtype=float))
    stride = max(1, data.shape[0] // MAX_PREVIEW)
    d = data[::stride]
    epoch0 = refdate.timestamp()
    return {
        "t_epoch": (epoch0 + d[:, 0]).tolist(),
        # NaN (e.g. in a hand-edited wind.txt) is invalid JSON -> null
        "columns": [[float(v) if np.isfinite(v) else None for v in d[:, i]]
                    for i in range(1, d.shape[1])],
        "n": int(data.shape[0]),
        "t0_epoch": float(epoch0 + data[0, 0]),
        "t1_epoch": float(epoch0 + data[-1, 0]),
    }


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
    """Synthetic series preview - computed only, nothing written."""
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
    out = {"series": _series_payload(data, refdate), "labels": KINDS[kind]["cols"]}
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
    direction is averaged circularly via unit-vector components."""
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
    t_out = seconds[0] + (unique_bins + 0.5) * width

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


@route("GET", "/api/conditions/raw_series")
def _raw_series(handler, query, tail):
    entry = _get_raw_entry(query.get("id", ""))
    if entry is None:
        send_error_json(handler, "unknown raw series", 404)
        return
    t, cols = _load_raw_series(entry)
    stride = max(1, t.size // MAX_PREVIEW)
    payload = {
        "t_epoch": t[::stride].tolist(),
        "columns": [[float(v) if np.isfinite(v) else None for v in cols[::stride, i]]
                    for i in range(cols.shape[1])],
        "labels": entry.get("labels") or KINDS[entry["kind"]]["cols"][:cols.shape[1]],
        "n": int(t.size),
        "t0_epoch": float(t[0]) if t.size else None,
        "t1_epoch": float(t[-1]) if t.size else None,
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


@route("POST", "/api/conditions/fetch")
def _fetch(handler, body, tail):
    current = project.require()
    source = body.get("source")
    kind = body.get("kind")
    if kind not in KINDS:
        send_error_json(handler, f"unknown kind '{kind}'")
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

    def _run(job):
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
            series = waterinfo.fetch_series(station, kind, date0, date1, job)
            epoch, base = series[0]
            cols = [base]
            for epoch_i, values_i in series[1:]:
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

        # store as a raw series object; the user inspects/cleans it and
        # then applies it to wind.txt/tide.txt/waves.txt explicitly
        label = body.get("label") or \
            f"{station_name} {date0:%Y-%m-%d} — {date1:%Y-%m-%d}"
        job.update(message="storing raw series")
        entry = _save_raw_entry(kind, source, label,
                                np.asarray(epoch, dtype="float64"), cols)
        return {"entry": entry}

    send_json(handler, {"job": jobs.start(f"fetch {source} {kind}", _run)})
