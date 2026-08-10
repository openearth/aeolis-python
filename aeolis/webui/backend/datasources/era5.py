"""ERA5 reanalysis wind source (Copernicus Climate Data Store).

Requires the ``cdsapi`` package and a configured CDS API key
(~/.cdsapirc, see https://cds.climate.copernicus.eu/how-to-api).
ERA5 lives on a regular 0.25 deg lat/lon grid; the GUI shows nearby
cell centres for the user to pick, then downloads 10 m u/v wind for the
simulation window and converts to an AeoLiS wind.txt (speed +
nautical direction).
"""

from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np


def configured():
    try:
        import cdsapi  # noqa: F401
    except ImportError:
        return False, "cdsapi not installed - pip install aeolis[webui-data]"
    if not (Path.home() / ".cdsapirc").is_file():
        return False, ("no CDS API key found (~/.cdsapirc) - create an account at "
                       "cds.climate.copernicus.eu and store your key")
    return True, None


def nearby_cells(lon, lat, n=3):
    """The n*n ERA5 0.25 deg cell centres around (lon, lat)."""
    lon0 = round(lon * 4) / 4
    lat0 = round(lat * 4) / 4
    cells = []
    half = n // 2
    for dj in range(-half, half + 1):
        for di in range(-half, half + 1):
            cells.append({
                "id": f"era5_{lat0 + dj * 0.25:.2f}_{lon0 + di * 0.25:.2f}",
                "lon": lon0 + di * 0.25,
                "lat": lat0 + dj * 0.25,
                "name": f"ERA5 cell ({lat0 + dj * 0.25:.2f}N, {lon0 + di * 0.25:.2f}E)",
            })
    return cells


def check(bounds_lonlat, job=None):
    ok, reason = configured()
    return {
        "available": ok,
        "years": [{"year": y} for y in range(1940, datetime.now(timezone.utc).year + 1)] if ok else [],
        "notes": reason or "ERA5 hourly 10 m wind, 0.25 deg grid, 1940-present.",
    }


def download_wind(lon, lat, date0, date1, dest_dir, job=None):
    """Download hourly 10 m u/v wind at one cell; returns (times, u, v)
    with times as datetimes (UTC). Result cached as npz."""
    ok, reason = configured()
    if not ok:
        raise RuntimeError(reason)
    import cdsapi

    stamp = f"{lon:.2f}_{lat:.2f}_{date0:%Y%m%d}_{date1:%Y%m%d}"
    cache = dest_dir / f"era5_wind_{stamp}.npz"
    if cache.exists():
        data = np.load(cache, allow_pickle=False)
        return data["t"], data["u"], data["v"], cache

    client = cdsapi.Client(quiet=True)

    # One CDS request per calendar year: a single request for a long
    # period enumerates the year x month x day cross-product and blows
    # past the CDS cost limit ("Your request is too large"). Yearly
    # chunks (max ~17.5k fields) always fit, and finished years stay
    # cached on disk so a retry only fetches what is missing.
    year0, year1 = date0.year, date1.year
    n_chunks = year1 - year0 + 1
    parts = []
    for k, year in enumerate(range(year0, year1 + 1)):
        d0 = max(date0, datetime(year, 1, 1, tzinfo=timezone.utc))
        d1 = min(date1, datetime(year, 12, 31, tzinfo=timezone.utc))
        part = dest_dir / f"era5_wind_{lon:.2f}_{lat:.2f}_{d0:%Y%m%d}_{d1:%Y%m%d}.nc"
        parts.append(part)
        if part.exists():
            continue
        if job:
            job.update(progress=k / n_chunks,
                       message=f"ERA5 {year} ({k + 1}/{n_chunks}) — waiting in the "
                               "Copernicus (CDS) server queue; this is outside our "
                               "control and can take minutes to hours when busy. "
                               "Finished years are cached, so retrying later resumes.")
        days = (d1 - d0).days + 1
        dates = [d0 + timedelta(days=i) for i in range(days)]
        months = sorted({f"{d.month:02d}" for d in dates})
        day_list = sorted({f"{d.day:02d}" for d in dates})
        client.retrieve(
            "reanalysis-era5-single-levels",
            {
                "product_type": "reanalysis",
                "variable": ["10m_u_component_of_wind", "10m_v_component_of_wind"],
                "year": [f"{year}"], "month": months, "day": day_list,
                "time": [f"{h:02d}:00" for h in range(24)],
                "area": [lat + 0.01, lon - 0.01, lat - 0.01, lon + 0.01],  # N W S E
                "format": "netcdf",
            },
            str(part),
        )

    if job:
        job.update(progress=0.95, message="reading ERA5 data")

    import netCDF4
    from aeolis.webui.backend.util import NC_LOCK
    all_t, all_u, all_v = [], [], []
    with NC_LOCK:
        for part in parts:
            ds = netCDF4.Dataset(part)
            try:
                tvar = ds.variables.get("time") or ds.variables.get("valid_time")
                times = netCDF4.num2date(tvar[:], tvar.units)
                u = np.asarray(ds.variables["u10"][:]).reshape(len(times), -1)[:, 0]
                v = np.asarray(ds.variables["v10"][:]).reshape(len(times), -1)[:, 0]
            finally:
                ds.close()
            all_t.append(np.array([np.datetime64(str(x)) for x in times]))
            all_u.append(u)
            all_v.append(v)

    t = np.concatenate(all_t)
    u = np.concatenate(all_u)
    v = np.concatenate(all_v)
    order = np.argsort(t)
    t, u, v = t[order], u[order], v[order]
    np.savez_compressed(cache, t=t.astype("datetime64[s]").astype("int64"),
                        u=u.astype("float32"), v=v.astype("float32"))
    data = np.load(cache)
    return data["t"], data["u"], data["v"], cache


def to_speed_direction(u, v, convention="nautical"):
    """u/v components -> (speed, direction). Nautical = direction the
    wind comes FROM, clockwise from north."""
    speed = np.hypot(u, v)
    if convention == "cartesian":
        direction = np.mod(np.degrees(np.arctan2(v, u)), 360.0)
    else:
        direction = np.mod(270.0 - np.degrees(np.arctan2(v, u)), 360.0)
    return speed, direction
