"""JarKus transect source (Deltares OPeNDAP).

Annual cross-shore profiles of the entire Dutch coast since 1965. The
combined dataset lives in a single netCDF on the Deltares OPeNDAP
server; we subset transects by their RD coordinates and store the
selected years as x/y/z point sets (.npz) in the project rawdata.
"""

from datetime import datetime, timezone

import numpy as np

from aeolis.webui.backend.util import NC_LOCK

URL = ("https://opendap.deltares.nl/thredds/dodsC/opendap/"
       "rijkswaterstaat/jarkus/profiles/transect_r20180914.nc")
URL_FALLBACK = ("https://opendap.deltares.nl/thredds/dodsC/opendap/"
                "rijkswaterstaat/jarkus/profiles/transect.nc")


def _open():
    import netCDF4
    for url in (URL_FALLBACK, URL):
        try:
            return netCDF4.Dataset(url)
        except OSError:
            continue
    raise RuntimeError(
        "could not reach the Deltares OPeNDAP server (JarKus); "
        "check your network connection"
    )


def _mask_area(ds, bounds):
    """Boolean mask over (alongshore, cross_shore) points within bounds."""
    minx, miny, maxx, maxy = bounds
    x = ds.variables["x"][:]     # (alongshore, cross_shore) RD coordinates
    y = ds.variables["y"][:]
    inside = (x >= minx) & (x <= maxx) & (y >= miny) & (y <= maxy)
    return np.asarray(inside), np.asarray(x), np.asarray(y)


def check(bounds, job=None):
    with NC_LOCK:
        ds = _open()
        try:
            if job:
                job.update(progress=0.2, message="locating transects")
            inside, _, _ = _mask_area(ds, bounds)
            rows = np.where(inside.any(axis=1))[0]
            if rows.size == 0:
                return {"available": False, "years": [], "notes": "no JarKus transects in this area"}

            times = ds.variables["time"]
            import netCDF4
            years_all = [d.year for d in netCDF4.num2date(times[:], times.units)]
            npoints = int(inside.sum())
            years = []
            for i, year in enumerate(years_all):
                if job and i % 10 == 0:
                    job.update(progress=0.2 + 0.8 * i / len(years_all), message=f"scanning {year}")
                years.append({"year": int(year), "est_bytes": int(npoints * 12)})
            return {
                "available": True,
                "years": years,
                "notes": f"{rows.size} transects intersect the area; per-year data "
                         "gaps are dropped on download.",
            }
        finally:
            ds.close()


def download(bounds, years, dest_dir, job=None):
    import netCDF4
    entries = []
    with NC_LOCK:
        ds = _open()
        try:
            inside, x, y = _mask_area(ds, bounds)
            rows = np.where(inside.any(axis=1))[0]
            if rows.size == 0:
                return []
            times = ds.variables["time"]
            dates_all = netCDF4.num2date(times[:], times.units)
            years_all = np.array([d.year for d in dates_all])
            zvar = ds.variables["altitude"]

            for n, year in enumerate(years):
                if job:
                    job.update(progress=n / max(1, len(years)), message=f"JarKus {year}")
                    if job.cancel_requested:
                        break
                t_idx = np.where(years_all == year)[0]
                if t_idx.size == 0:
                    continue
                out_name = f"jarkus_{year}.npz"
                out_path = dest_dir / out_name
                if not out_path.exists():
                    z = zvar[t_idx[0], rows, :]
                    z = np.ma.filled(z, np.nan).astype("float32")
                    xs = x[rows, :].astype("float64")
                    ys = y[rows, :].astype("float64")
                    valid = np.isfinite(z) & inside[rows, :]
                    if not valid.any():
                        continue
                    np.savez_compressed(
                        out_path,
                        x=xs[valid], y=ys[valid], z=z[valid],
                    )
                entries.append({
                    "source": "jarkus",
                    "year": int(year),
                    "date": str(dates_all[t_idx[0]])[:10],   # actual survey date
                    "kind": "points",
                    "path": f"gui/rawdata/{out_name}",
                    "crs": 28992,
                    "res": None,
                    "bounds": list(bounds),
                    "label": f"JarKus {year}",
                    "downloaded": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M"),
                })
        finally:
            ds.close()
    return entries
