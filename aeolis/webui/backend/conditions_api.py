"""Boundary conditions API (wind / water levels / waves).

Synthetic generation, measured data (waterinfo.rws.nl), and reanalysis
wind (ERA5) are all converted to the AeoLiS text formats:

    wind_file  columns [t, speed, direction]
    tide_file  columns [t, water level]
    wave_file  columns [t, Hs, Tp]

with t in seconds since the configuration ``refdate``. Raw downloads
stay cached in gui/rawdata.
"""

from datetime import datetime, timedelta, timezone

import numpy as np

from aeolis.webui.backend import jobs, project
from aeolis.webui.backend.config_api import load_config
from aeolis.webui.backend.datasources import era5, synthetic, waterinfo
from aeolis.webui.backend.grid_api import _load_current_grid, patch_config
from aeolis.webui.backend.httpd import route
from aeolis.webui.backend.util import send_error_json, send_json

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
        "columns": [d[:, i].tolist() for i in range(1, d.shape[1])],
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


def _generate_synthetic(kind, body, values):
    tstart = float(body.get("tstart", values.get("tstart") or 0.0))
    tstop = float(body.get("tstop", values.get("tstop") or 3600.0))
    dt = float(body.get("dt", 3600.0))
    if kind == "wind":
        return synthetic.wind(tstart, tstop, dt, body.get("speed", {}), body.get("direction", {}))
    if kind == "tide":
        return synthetic.tide(tstart, tstop, dt, body.get("level", {}))
    return synthetic.waves(tstart, tstop, dt, body.get("hs", {}), body.get("tp", {}))


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
    send_json(handler, {"series": _series_payload(data, refdate),
                        "labels": KINDS[kind]["cols"]})


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
    """Bin-average a series to hourly/daily values. Wind direction is
    averaged circularly via unit-vector components."""
    if interval not in ("hour", "day") or seconds.size == 0:
        return seconds, cols
    width = 3600.0 if interval == "hour" else 86400.0
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

    info = KINDS[kind]
    filename = body.get("filename") or values.get(info["key"]) or info["default"]

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
            seconds = t.astype("float64") - refdate.timestamp()
            data = np.column_stack([seconds, speed, direction])
        elif source == "waterinfo":
            station = body.get("station")
            series = waterinfo.fetch_series(station, kind, date0, date1, job)
            epoch0, base = series[0]
            cols = [base]
            for epoch_i, values_i in series[1:]:
                # nearest-neighbour align extra quantities on the first one
                idx = np.searchsorted(epoch_i, epoch0).clip(0, len(values_i) - 1)
                cols.append(values_i[idx])
            resample = body.get("resample")
            if resample:
                job.update(message=f"resampling to {resample}ly means")
                epoch0, cols = _resample(epoch0, cols, kind, resample)
            seconds = epoch0 - refdate.timestamp()
            data = np.column_stack([seconds] + list(cols))
        else:
            raise RuntimeError(f"unknown source '{source}'")

        job.update(message=f"writing {filename}")
        synthetic.write_series(current.root / filename, data)
        patch_config({info["key"]: filename})
        return {
            "file": filename, "rows": int(data.shape[0]),
            "series": _series_payload(data, refdate),
        }

    send_json(handler, {"job": jobs.start(f"fetch {source} {kind}", _run)})
