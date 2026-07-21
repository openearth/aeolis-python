"""Vaklodingen source (Deltares OPeNDAP).

Bathymetric surveys of the Dutch coastal waters on a 20 m grid, stored
as one netCDF per "kaartblad" tile on the Deltares OPeNDAP server. The
tile extents are not encoded in the filenames, so on first use a tile
index (name -> x/y extent) is built by probing each tile once and cached
globally in ~/.aeolis_webui/vaklodingen_index.json.
"""

import re
import urllib.request
from datetime import datetime, timezone

import numpy as np

from aeolis.webui.backend import settings
from aeolis.webui.backend.util import NC_LOCK, load_json, save_json

CATALOG = ("https://opendap.deltares.nl/thredds/catalog/opendap/"
           "rijkswaterstaat/vaklodingen/catalog.xml")
DODS = ("https://opendap.deltares.nl/thredds/dodsC/opendap/"
        "rijkswaterstaat/vaklodingen/")
INDEX_FILE = settings.APP_DIR / "vaklodingen_index.json"


def _list_tiles():
    with urllib.request.urlopen(CATALOG, timeout=20) as response:
        xml = response.read().decode("utf-8", "replace")
    names = sorted(set(re.findall(r'name="(vaklodingen[^"]*\.nc)"', xml)))
    if not names:
        raise RuntimeError("no vaklodingen tiles found in the THREDDS catalog")
    return names


def _tile_index(job=None):
    """name -> {x0, x1, y0, y1} for every tile (cached globally)."""
    index = load_json(INDEX_FILE, default=None)
    names = _list_tiles()
    if index and set(index) >= set(names):
        return index

    import netCDF4
    index = index or {}
    todo = [n for n in names if n not in index]
    for i, name in enumerate(todo):
        if job:
            job.update(progress=i / len(todo),
                       message=f"indexing tiles {i + 1}/{len(todo)} (one-time)")
            if job.cancel_requested:
                break
        try:
            with NC_LOCK:
                ds = netCDF4.Dataset(DODS + name)
                try:
                    x = ds.variables["x"][:]
                    y = ds.variables["y"][:]
                    index[name] = {
                        "x0": float(np.min(x)), "x1": float(np.max(x)),
                        "y0": float(np.min(y)), "y1": float(np.max(y)),
                    }
                finally:
                    ds.close()
        except OSError:
            continue
        if i % 10 == 9:
            save_json(INDEX_FILE, index)
    save_json(INDEX_FILE, index)
    return index


def _tiles_for(bounds, job=None):
    minx, miny, maxx, maxy = bounds
    index = _tile_index(job)
    hits = []
    for name, ext in index.items():
        if ext["x1"] >= minx and ext["x0"] <= maxx and ext["y1"] >= miny and ext["y0"] <= maxy:
            hits.append(name)
    return hits


def check(bounds, job=None):
    import netCDF4
    tiles = _tiles_for(bounds, job)
    if not tiles:
        return {"available": False, "years": [],
                "notes": "no vaklodingen tiles intersect this area"}

    years = {}
    for i, name in enumerate(tiles):
        if job:
            job.update(progress=i / len(tiles), message=f"scanning {name}")
            if job.cancel_requested:
                break
        try:
            with NC_LOCK:
                ds = netCDF4.Dataset(DODS + name)
                try:
                    times = ds.variables["time"]
                    dates = netCDF4.num2date(times[:], times.units)
                    for d in dates:
                        key = int(d.year)
                        years.setdefault(key, {"year": key, "surveys": 0, "est_bytes": 0})
                        years[key]["surveys"] += 1
                        years[key]["est_bytes"] += 500 * 625 * 4
                finally:
                    ds.close()
        except OSError:
            continue
    return {
        "available": bool(years),
        "years": [years[k] for k in sorted(years)],
        "notes": f"{len(tiles)} tiles intersect the area; survey dates vary per tile.",
    }


def download(bounds, years, dest_dir, job=None):
    import netCDF4
    minx, miny, maxx, maxy = bounds
    tiles = _tiles_for(bounds)
    entries = []
    years = set(int(y) for y in years)

    for n, name in enumerate(tiles):
        if job:
            job.update(progress=n / max(1, len(tiles)), message=f"tile {name}")
            if job.cancel_requested:
                break
        with NC_LOCK:
            try:
                ds = netCDF4.Dataset(DODS + name)
            except OSError:
                continue
            try:
                x = np.asarray(ds.variables["x"][:])
                y = np.asarray(ds.variables["y"][:])
                ix = np.where((x >= minx) & (x <= maxx))[0]
                iy = np.where((y >= miny) & (y <= maxy))[0]
                if ix.size == 0 or iy.size == 0:
                    continue
                times = ds.variables["time"]
                dates = netCDF4.num2date(times[:], times.units)
                zvar = ds.variables["z"]
                for t, date in enumerate(dates):
                    if int(date.year) not in years:
                        continue
                    stamp = f"{date.year:04d}{date.month:02d}"
                    tile_id = name.replace("vaklodingen", "").replace(".nc", "")
                    out_name = f"vaklodingen{tile_id}_{stamp}.npz"
                    out_path = dest_dir / out_name
                    if not out_path.exists():
                        z = zvar[t, iy.min():iy.max() + 1, ix.min():ix.max() + 1]
                        z = np.ma.filled(z, np.nan).astype("float32")
                        if not np.isfinite(z).any():
                            continue
                        np.savez_compressed(
                            out_path,
                            x=x[ix.min():ix.max() + 1],
                            y=y[iy.min():iy.max() + 1],
                            z=z,
                            date=str(date)[:10],
                        )
                    entries.append({
                        "source": "vaklodingen",
                        "year": int(date.year),
                        "date": str(date)[:10],
                        "kind": "raster_nc",
                        "path": f"gui/rawdata/{out_name}",
                        "crs": 28992,
                        "res": 20,
                        "bounds": list(bounds),
                        "label": f"Vaklodingen {tile_id} {str(date)[:10]}",
                        "downloaded": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M"),
                    })
            finally:
                ds.close()
    return entries
