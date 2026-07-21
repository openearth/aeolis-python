"""Rijkswaterstaat waterinfo source (WaterWebservices / DDL 2.0).

Measured wind, water levels and waves from the RWS monitoring network
(waterinfo.rws.nl) via the 2025 "DD-API 2.0" JSON services. No API key
is required yet; a dummy X-API-KEY header is sent for forward
compatibility. Station coordinates are ETRS89 lat/lon.

Quantities:
    wind  -> WINDSHD (speed, m/s) + WINDRTG (direction, deg)
    tide  -> WATHTE (water level, cm -> m)
    wave  -> Hm0 (significant wave height, cm -> m) + Tm02 (s)
"""

import json
import urllib.error
import urllib.request
from datetime import datetime, timezone

import numpy as np

BASE = "https://ddapi20-waterwebservices.rijkswaterstaat.nl"
CATALOG_URL = f"{BASE}/METADATASERVICES/OphalenCatalogus"
DATA_URL = f"{BASE}/ONLINEWAARNEMINGENSERVICES/OphalenWaarnemingen"

QUANTITIES = {
    "wind": {"codes": ["WINDSHD", "WINDRTG"], "scale": [1.0, 1.0]},
    "tide": {"codes": ["WATHTE"], "scale": [0.01]},
    "wave": {"codes": ["Hm0", "Tm02"], "scale": [0.01, 1.0]},
}

_catalog_cache = None


def _post(url, payload, timeout=60):
    request = urllib.request.Request(
        url,
        data=json.dumps(payload).encode("utf-8"),
        headers={
            "Content-Type": "application/json",
            "X-API-KEY": "aeolis-webui",
        },
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            return json.loads(response.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        raise RuntimeError(f"waterinfo: HTTP {exc.code} from {url}") from exc
    except OSError as exc:
        raise RuntimeError(f"waterinfo unreachable: {exc}") from exc


def catalog():
    """Full station/quantity catalog (cached per process)."""
    global _catalog_cache
    if _catalog_cache is None:
        _catalog_cache = _post(CATALOG_URL, {"CatalogusFilter": {"Grootheden": True}})
        if not _catalog_cache.get("Succesvol", False):
            _catalog_cache = None
            raise RuntimeError("waterinfo: catalog request unsuccessful")
    return _catalog_cache


def stations_for(kind, lon=None, lat=None, max_km=75.0):
    """Stations measuring *kind*, sorted/filtered by distance to
    (lon, lat) when given."""
    if kind not in QUANTITIES:
        raise ValueError(f"unknown quantity kind '{kind}'")
    code = QUANTITIES[kind]["codes"][0]
    cat = catalog()

    meta_ids = {
        m["AquoMetadata_MessageID"]
        for m in cat.get("AquoMetadataLijst", [])
        if m.get("Grootheid", {}).get("Code") == code
    }
    loc_ids = {
        link["Locatie_MessageID"]
        for link in cat.get("AquoMetadataLocatieLijst", [])
        if link.get("AquoMetaData_MessageID") in meta_ids
    }

    stations = []
    for loc in cat.get("LocatieLijst", []):
        if loc.get("Locatie_MessageID") not in loc_ids:
            continue
        station = {
            "id": loc.get("Code"),
            "name": loc.get("Naam") or loc.get("Code"),
            "lon": loc.get("Lon"),
            "lat": loc.get("Lat"),
        }
        if lon is not None and lat is not None and station["lon"] is not None:
            station["dist_km"] = _haversine(lon, lat, station["lon"], station["lat"])
        stations.append(station)

    if lon is not None and lat is not None:
        stations = [s for s in stations if s.get("dist_km", 1e9) <= max_km]
        stations.sort(key=lambda s: s.get("dist_km", 1e9))
    return stations


def _haversine(lon1, lat1, lon2, lat2):
    r = 6371.0
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dl = np.radians(lon2 - lon1)
    a = np.sin((p2 - p1) / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return float(2 * r * np.arcsin(np.sqrt(a)))


def fetch_series(station, kind, date0, date1, job=None):
    """Measured series per quantity code -> [(epoch_seconds, values), ...]."""
    if kind not in QUANTITIES:
        raise ValueError(f"unknown quantity kind '{kind}'")
    spec = QUANTITIES[kind]
    results = []
    for code, scale in zip(spec["codes"], spec["scale"]):
        if job:
            job.update(message=f"waterinfo {station} {code}")
        payload = {
            "Locatie": {"Code": station},
            "AquoPlusWaarnemingMetadata": {
                "AquoMetadata": {"Grootheid": {"Code": code}},
            },
            "Periode": {
                "Begindatumtijd": date0.strftime("%Y-%m-%dT%H:%M:%S.000+00:00"),
                "Einddatumtijd": date1.strftime("%Y-%m-%dT%H:%M:%S.000+00:00"),
            },
        }
        data = _post(DATA_URL, payload, timeout=180)
        if not data.get("Succesvol", False):
            raise RuntimeError(
                f"waterinfo: {data.get('Foutmelding', 'request failed')} "
                f"(station {station}, {code})"
            )
        times, values = [], []
        for block in data.get("WaarnemingenLijst", []):
            for m in block.get("MetingenLijst", []):
                meetwaarde = m.get("Meetwaarde") or {}
                val = meetwaarde.get("Waarde_Numeriek")
                if val is None or val > 9e8:
                    continue
                stamp = m.get("Tijdstip")
                if not stamp:
                    continue
                dt = datetime.fromisoformat(stamp)
                if dt.tzinfo is None:
                    dt = dt.replace(tzinfo=timezone.utc)
                times.append(dt.astimezone(timezone.utc))
                values.append(float(val) * scale)
        if not times:
            raise RuntimeError(f"waterinfo: no data for {station}/{code} in this period")
        order = np.argsort(np.array([t.timestamp() for t in times]))
        epoch = np.array([t.timestamp() for t in times])[order]
        results.append((epoch, np.asarray(values)[order]))
    return results
