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

    if job:
        job.update(progress=-1, message="requesting ERA5 subset (CDS queue, may take minutes)")

    grib = dest_dir / f"era5_wind_{stamp}.nc"
    client = cdsapi.Client(quiet=True)
    days = (date1 - date0).days + 1
    dates = [date0 + timedelta(days=k) for k in range(days)]
    years = sorted({f"{d.year}" for d in dates})
    months = sorted({f"{d.month:02d}" for d in dates})
    day_list = sorted({f"{d.day:02d}" for d in dates})
    client.retrieve(
        "reanalysis-era5-single-levels",
        {
            "product_type": "reanalysis",
            "variable": ["10m_u_component_of_wind", "10m_v_component_of_wind"],
            "year": years, "month": months, "day": day_list,
            "time": [f"{h:02d}:00" for h in range(24)],
            "area": [lat + 0.01, lon - 0.01, lat - 0.01, lon + 0.01],  # N W S E
            "format": "netcdf",
        },
        str(grib),
    )

    import netCDF4
    from aeolis.webui.backend.util import NC_LOCK
    with NC_LOCK:
        ds = netCDF4.Dataset(grib)
        try:
            tvar = ds.variables.get("time") or ds.variables.get("valid_time")
            times = netCDF4.num2date(tvar[:], tvar.units)
            u = np.asarray(ds.variables["u10"][:]).reshape(len(times), -1)[:, 0]
            v = np.asarray(ds.variables["v10"][:]).reshape(len(times), -1)[:, 0]
        finally:
            ds.close()

    t = np.array([np.datetime64(str(x)) for x in times])
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
