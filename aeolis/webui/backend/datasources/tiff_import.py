"""Custom GeoTIFF import (topography raster).

Copies a local GeoTIFF into the project rawdata folder and records a
``kind: "raster"`` manifest entry — exactly the shape the RWS LiDAR /
Vaklodingen downloads produce, so interpolation and the map viewer
(``domain_api.load_raw``) handle it with no extra code.

Some RD-New GeoTIFFs (e.g. the Witteveen+Bos design variants) carry a
``LOCAL_CS["Netherlands-RD"]`` tag instead of a proper EPSG code, so
``crs.to_epsg()`` returns ``None``; we then fall back to EPSG:28992,
which is what those coordinates actually are.
"""

import re
import shutil
from datetime import datetime, timezone
from pathlib import Path


def _safe_stem(stem):
    """Filesystem/URL-safe version of the source filename stem."""
    return re.sub(r"[^0-9A-Za-z._-]+", "_", stem).strip("_") or "raster"


def import_file(path, dest_dir, crs=None):
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"file not found: {path}")

    import rasterio

    with rasterio.open(path) as ds:
        if ds.count < 1:
            raise ValueError("GeoTIFF has no raster bands")
        bands = int(ds.count)
        b = ds.bounds
        res = float(abs(ds.transform.a))
        detected = None
        try:
            detected = ds.crs.to_epsg() if ds.crs is not None else None
        except Exception:  # noqa: BLE001 - malformed/LOCAL_CS tags
            detected = None

    # LOCAL_CS RD-New tiffs report no EPSG; those Dutch coordinates are 28992.
    epsg = int(crs or detected or 28992)

    out_name = f"tiff_{_safe_stem(path.stem)}.tif"
    shutil.copyfile(path, dest_dir / out_name)

    return {
        "source": "tiff",
        "kind": "raster",
        "path": f"gui/rawdata/{out_name}",
        "origin": str(path),
        "crs": epsg,
        "res": res,
        "bounds": [float(b.left), float(b.bottom), float(b.right), float(b.top)],
        "bands": bands,        # number of channels (RGB/CIR/… kept for band-math)
        "band": 1,             # which band is displayed / used by default
        "label": f"TIFF {path.name}",
        "downloaded": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M"),
    }
