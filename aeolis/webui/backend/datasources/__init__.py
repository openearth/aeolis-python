"""Raw-data sources for the Domain and Conditions tabs.

Every source implements the same small interface so the GUI can offer a
uniform "select area -> check availability -> download" workflow:

    check(bounds, job)            -> {available, years, est_bytes, notes}
    download(bounds, years, dest_dir, job) -> [manifest entries]

``bounds`` is (minx, miny, maxx, maxy) in the source CRS (EPSG:28992 for
the Dutch bathymetry/elevation sources). Downloads are cached in the
project's ``gui/rawdata`` folder; a manifest entry looks like::

    {"id": "...", "source": "rws_lidar", "year": 2023, "kind": "raster",
     "path": "rawdata/rws_lidar_2023.tif", "crs": 28992,
     "bounds": [..], "res": 2.0, "downloaded": "..."}

Sources with heavy dependencies (rasterio) or network requirements
raise a clear error message the frontend shows verbatim.
"""

from aeolis.webui.backend.datasources import (
    jarkus, rws_lidar, tiff_import, vaklodingen, xyz_import,
)

SOURCES = {
    "rws_lidar": rws_lidar,
    "jarkus": jarkus,
    "vaklodingen": vaklodingen,
}

INFO = [
    {
        "id": "rws_lidar",
        "title": "RWS coastal LiDAR (kusthoogte)",
        "description": "Annual airborne LiDAR DTM of the Dutch coast "
                       "(downloads.rijkswaterstaatdata.nl); 2-5 m grid. Open data.",
        "kind": "raster",
        "crs": 28992,
    },
    {
        "id": "jarkus",
        "title": "JarKus transects",
        "description": "Annual cross-shore profile measurements of the entire "
                       "Dutch coast since 1965 (Deltares OPeNDAP). Open data.",
        "kind": "points",
        "crs": 28992,
    },
    {
        "id": "vaklodingen",
        "title": "Vaklodingen",
        "description": "Bathymetric surveys of the Dutch coastal waters and "
                       "estuaries, 20 m grid (Deltares OPeNDAP). Open data.",
        "kind": "raster",
        "crs": 28992,
    },
    {
        "id": "xyz",
        "title": "Custom sample file (*.xyz)",
        "description": "Import your own x/y/z samples (whitespace or comma "
                       "separated text).",
        "kind": "points",
        "crs": None,
    },
    {
        "id": "tiff",
        "title": "Custom topography (*.tif)",
        "description": "Import your own GeoTIFF elevation raster (e.g. a "
                       "survey DTM or design variant). RD-New assumed when "
                       "the file carries no EPSG code.",
        "kind": "raster",
        "crs": 28992,
    },
]


def get(source_id):
    if source_id not in SOURCES:
        raise ValueError(f"unknown data source '{source_id}'")
    return SOURCES[source_id]
