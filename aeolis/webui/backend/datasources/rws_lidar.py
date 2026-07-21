"""RWS coastal LiDAR (kusthoogte) source.

Annual airborne LiDAR DTMs of the Dutch coast, served as large striped
(Big)GeoTIFFs on https://downloads.rijkswaterstaatdata.nl. Windowed
reads via GDAL /vsicurl (rasterio) fetch only the blocks intersecting
the requested area.

The per-year file table follows the download scripts in the boschplaat
project; filename conventions changed over the years. Availability is
checked with cheap HTTP HEAD requests; whether a year actually covers
the requested area is only known after the (windowed) download, which
reports "no coverage" when the window is empty.
"""

import hashlib
import urllib.request
from datetime import datetime, timezone

import numpy as np


def bounds_tag(bounds):
    """Short area fingerprint used in cached filenames, so data
    downloaded for one grid location is never reused for another."""
    key = "_".join(f"{v:.0f}" for v in bounds)
    return hashlib.sha1(key.encode()).hexdigest()[:6]

BASE = "https://downloads.rijkswaterstaatdata.nl"
NODATA = -9999.0

# year -> (resolution m, [relative paths]; first is primary, rest fill gaps)
FILES = {
    1997: (5, ["hoogte_1997/kust_1997.tif"]),
    1998: (5, ["hoogte_1998/kust_1998.tif"]),
    1999: (5, ["hoogte_1999/kust_1999.tif"]),
    2000: (5, ["hoogte_2000/kust_2000.tif"]),
    2001: (5, ["hoogte_2001/kust_2001.tif"]),
    2002: (5, ["hoogte_2002/kust_2002.tif"]),
    2003: (5, ["hoogte_2003/kust_2003.tif"]),
    2004: (5, ["hoogte_2004/kust_2004.tif", "hoogte_2004/wadden_2004.tif"]),
    2005: (5, ["hoogte_2005/kust_2005.tif"]),
    2006: (5, ["hoogte_2006/kust_2006.tif", "hoogte_2006/kust_2006_2.tif"]),
    2007: (5, ["hoogte_2007/kust_2007.tif", "hoogte_2007/wadden_2007.tif"]),
    2008: (5, ["hoogte_2008/kust_2008.tif", "hoogte_2008/wadden_2008.tif"]),
    2009: (5, ["hoogte_2009/kust_2009.tif", "hoogte_2009/wadden_2009.tif"]),
    2010: (5, ["hoogte_2010/kust_2010.tif", "hoogte_2010/vliestroom_2010.tif"]),
    2011: (5, ["hoogte_2011/kust_2011.tif", "hoogte_2011/waddenzee_mid_2011.tif"]),
    2012: (5, ["hoogte_2012/kust_2012.tif"]),
    2013: (2, ["hoogte_2013/kust_2013.tif"]),
    2014: (2, ["hoogte_2014/kust_2014.tif"]),
    2015: (2, ["hoogte_2015/kust_2015.tif"]),
    2016: (2, ["hoogte_2016/kust_2016_dtm_2.tif", "hoogte_2016/waddenzee_2016_dtm_2.tif"]),
    2017: (2, ["hoogte_2017/kust_2017_dtm_2.tif", "hoogte_2017/kust_2017_dtm_2_2.tif",
               "hoogte_2017/waddenzee_2017_dtm_2.tif"]),
    2018: (2, ["hoogte_2018/kust_2018_dtm_2.tif"]),
    2019: (2, ["hoogte_2019/kust_2019_dtm_2m.tif"]),
    2020: (2, ["hoogte_2020/hoogte_kust_2020_dtm_2m.tif"]),
    2021: (2, ["hoogte_2021/hoogte_kust_2021_dtm_2m.tif"]),
    2022: (2, ["hoogte_2022/hoogte_kust_2022_dtm_2m.tif",
               "hoogte_2022/hoogte_wadden_vliestroom_2022_dtm_2m.tif"]),
    2023: (2, ["hoogte_2023/hoogte_kust_2023_dtm_2m.tif",
               "hoogte_2023/hoogte_wadden_eieramelanderzeegat_2023_dtm_2m.tif"]),
    2024: (2, ["hoogte_2024/hoogte_kust_2024_dtm_2m.tif"]),
    2025: (2, ["hoogte_2025/hoogte_kust_2025_dtm_2m.tif"]),
    2026: (2, ["hoogte_2026/hoogte_kust_2026_dtm_2m.tif"]),
}

# GDAL/CURL tuning: required for these large striped BigTIFFs, otherwise
# range reads can silently return all-nodata.
GDAL_OPTS = dict(
    GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
    GDAL_HTTP_MULTIRANGE="YES",
    GDAL_HTTP_MERGE_CONSECUTIVE_RANGES="YES",
    CPL_VSIL_CURL_USE_HEAD="NO",
    GDAL_HTTP_VERSION="2",
    VSI_CACHE="TRUE",
    CPL_VSIL_CURL_CHUNK_SIZE="10485760",
)


def _rasterio():
    try:
        import rasterio  # noqa: F401
        return rasterio
    except ImportError as exc:
        raise RuntimeError(
            "rasterio is required for LiDAR downloads - "
            "install with: pip install aeolis[webui-data]"
        ) from exc


def _head_exists(url, timeout=6):
    request = urllib.request.Request(url, method="HEAD")
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            return 200 <= response.status < 300
    except Exception:
        return False


def check(bounds, job=None):
    minx, miny, maxx, maxy = bounds
    area = max(0.0, maxx - minx) * max(0.0, maxy - miny)
    years = []
    total = len(FILES)
    for i, (year, (res, paths)) in enumerate(sorted(FILES.items())):
        if job:
            job.update(progress=i / total, message=f"checking {year}")
            if job.cancel_requested:
                break
        if _head_exists(f"{BASE}/{paths[0]}"):
            years.append({
                "year": year,
                "res": res,
                "est_bytes": int(area / res / res * 4),
            })
    return {
        "available": bool(years),
        "years": years,
        "notes": "Yearly flights; actual coverage of your area is only known "
                 "after download (window may be empty for some years).",
    }


def download(bounds, years, dest_dir, job=None):
    rasterio = _rasterio()
    from rasterio.env import Env
    from rasterio.transform import Affine
    from rasterio.windows import Window, from_bounds

    minx, miny, maxx, maxy = bounds
    entries = []

    with Env(**GDAL_OPTS):
        for n, year in enumerate(years):
            if year not in FILES:
                continue
            if job:
                job.update(progress=n / max(1, len(years)), message=f"downloading {year}")
                if job.cancel_requested:
                    break
            res, paths = FILES[year]
            # the filename carries an area tag: re-downloading for a new
            # grid location must never reuse a cached file of the old one
            out_name = f"rws_lidar_{year}_{res}m_{bounds_tag(bounds)}.tif"
            out_path = dest_dir / out_name
            if out_path.exists():
                entries.append(_entry(year, res, out_name, bounds))
                continue

            canvas = None
            canvas_tr = None
            for k, rel in enumerate(paths):
                try:
                    with rasterio.open(f"/vsicurl/{BASE}/{rel}") as ds:
                        win = from_bounds(minx, miny, maxx, maxy, transform=ds.transform)
                        win = win.intersection(Window(0, 0, ds.width, ds.height))
                        if win.width <= 0 or win.height <= 0:
                            continue
                        win = win.round_offsets().round_lengths()
                        data = ds.read(1, window=win).astype("float32")
                        if ds.nodata is not None:
                            data[data == np.float32(ds.nodata)] = NODATA
                        data[~np.isfinite(data)] = NODATA
                        tr = ds.window_transform(win)
                except Exception as exc:  # noqa: BLE001 - per-file, keep going
                    if k == 0:
                        raise RuntimeError(f"LiDAR {year}: {exc}") from exc
                    continue

                if canvas is None:
                    canvas, canvas_tr = data, tr
                else:
                    # paste secondary files into gaps only
                    col = round((tr.c - canvas_tr.c) / canvas_tr.a)
                    row = round((tr.f - canvas_tr.f) / canvas_tr.e)
                    r0, c0 = max(row, 0), max(col, 0)
                    r1 = min(row + data.shape[0], canvas.shape[0])
                    c1 = min(col + data.shape[1], canvas.shape[1])
                    if r1 > r0 and c1 > c0:
                        sub = canvas[r0:r1, c0:c1]
                        src = data[r0 - row:r1 - row, c0 - col:c1 - col]
                        gaps = sub == NODATA
                        sub[gaps] = src[gaps]

            if canvas is None or np.all(canvas == NODATA):
                continue  # no coverage this year

            profile = {
                "driver": "GTiff", "height": canvas.shape[0], "width": canvas.shape[1],
                "count": 1, "dtype": "float32", "crs": "EPSG:28992",
                "transform": canvas_tr, "nodata": NODATA,
                "compress": "lzw", "tiled": True,
            }
            with rasterio.open(out_path, "w", **profile) as out:
                out.write(canvas, 1)
            entries.append(_entry(year, res, out_name, bounds))

    return entries


def _entry(year, res, filename, bounds):
    return {
        "source": "rws_lidar",
        "year": year,
        "kind": "raster",
        "path": f"gui/rawdata/{filename}",
        "crs": 28992,
        "res": res,
        "bounds": list(bounds),
        "label": f"LiDAR {year} ({res} m)",
        "downloaded": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M"),
    }
