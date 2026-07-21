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


class LimitExceeded(RuntimeError):
    """The API's max-observations-per-request cap was hit (263 088)."""


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
            body = response.read().decode("utf-8", "replace")
    except urllib.error.HTTPError as exc:
        body = exc.read().decode("utf-8", "replace")
        try:
            data = json.loads(body)
        except json.JSONDecodeError:
            raise RuntimeError(
                f"waterinfo: HTTP {exc.code} from {url}: {body[:160]}"
            ) from exc
        message = data.get("Foutmelding", "")
        if "overschreden" in message:
            raise LimitExceeded(message)
        raise RuntimeError(f"waterinfo: HTTP {exc.code}: {message or body[:160]}") from exc
    except OSError as exc:
        raise RuntimeError(f"waterinfo unreachable: {exc}") from exc
    try:
        return json.loads(body)
    except json.JSONDecodeError as exc:
        raise RuntimeError(
            f"waterinfo: unexpected non-JSON response ({body[:120]}...)"
        ) from exc


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


def _fetch_window(station, code, scale, date0, date1):
    """One raw request -> (epoch array, values array); may raise
    LimitExceeded when the window holds too many observations."""
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
        message = data.get("Foutmelding", "request failed")
        if "overschreden" in message:
            raise LimitExceeded(message)
        # "no data" style failures return unsuccessful too
        return np.array([]), np.array([])
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
            times.append(dt.astimezone(timezone.utc).timestamp())
            values.append(float(val) * scale)
    if not times:
        return np.array([]), np.array([])
    order = np.argsort(times)
    return np.asarray(times)[order], np.asarray(values)[order]


def fetch_series(station, kind, date0, date1, job=None):
    """Measured series per quantity code -> [(epoch_seconds, values), ...].

    The DD-API caps a single request at 263 088 observations, so long
    periods are fetched in adaptive chunks (start ~6 months, halve on a
    limit error) and concatenated.
    """
    from datetime import timedelta

    if kind not in QUANTITIES:
        raise ValueError(f"unknown quantity kind '{kind}'")
    spec = QUANTITIES[kind]
    results = []
    total_seconds = max(1.0, (date1 - date0).total_seconds())

    for code, scale in zip(spec["codes"], spec["scale"]):
        chunk = timedelta(days=180)
        epochs, values = [], []
        cursor = date0
        while cursor < date1:
            if job:
                frac = (cursor - date0).total_seconds() / total_seconds
                job.update(progress=frac,
                           message=f"waterinfo {station} {code} "
                                   f"({cursor:%Y-%m}, chunk {chunk.days} d)")
                if job.cancel_requested:
                    break
            end = min(cursor + chunk, date1)
            try:
                e, v = _fetch_window(station, code, scale, cursor, end)
            except LimitExceeded:
                if chunk.days <= 7:
                    raise RuntimeError(
                        "waterinfo: observation limit hit even for a 7-day "
                        f"window ({station}/{code})"
                    )
                chunk = timedelta(days=max(7, chunk.days // 2))
                continue
            epochs.append(e)
            values.append(v)
            cursor = end
        epoch = np.concatenate(epochs) if epochs else np.array([])
        vals = np.concatenate(values) if values else np.array([])
        if epoch.size == 0:
            raise RuntimeError(f"waterinfo: no data for {station}/{code} in this period")
        results.append((epoch, vals))
    return results


def probe_period(station, kind, job=None):
    """Best-effort estimate of the available data period for a station:
    coarse 5-year probes with 30-day windows, refined to the year."""
    from datetime import timedelta

    code = QUANTITIES[kind]["codes"][0]
    scale = QUANTITIES[kind]["scale"][0]
    now = datetime.now(timezone.utc)

    def has_data(year):
        d0 = datetime(year, 6, 1, tzinfo=timezone.utc)
        if d0 > now:
            d0 = now - timedelta(days=30)
        try:
            e, _ = _fetch_window(station, code, scale, d0, d0 + timedelta(days=30))
        except (LimitExceeded, RuntimeError):
            return False
        return e.size > 0

    # find the earliest 5-year block with data
    earliest = None
    for year in range(1950, now.year + 1, 5):
        if job:
            job.update(message=f"probing {year}")
            if job.cancel_requested:
                break
        if has_data(year):
            earliest = year
            break
    if earliest is None:
        return None
    # refine within the block
    for year in range(max(1950, earliest - 4), earliest + 1):
        if has_data(year):
            earliest = year
            break
    has_recent = has_data(now.year) or has_data(now.year - 1)
    return {
        "from": earliest,
        "to": now.year if has_recent else None,
        "note": "estimated from coarse probes; gaps possible",
    }
