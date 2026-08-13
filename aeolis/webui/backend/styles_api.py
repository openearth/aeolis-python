"""User-defined colormap/style presets, shared across ALL projects.

A preset bundles everything about how a layer is drawn: colormap (a
trailing "!r" marks reversed), display mode (cells/dots), dot size,
opacity and optional fixed z-limits (null = auto from the data). The
store is a single JSON file in the per-user app dir; the frontend owns
the editing and sends the full list on every change.
"""

from aeolis.webui.backend import settings
from aeolis.webui.backend.httpd import route
from aeolis.webui.backend.util import load_json, save_json, send_error_json, send_json

# starter presets, written once when no styles file exists yet
DEFAULT_STYLES = [
    {"id": "elevation", "name": "Elevation", "cmap": "topo_dutch",
     "mode": "cells", "dotSize": 6, "opacity": 0.9, "min": -5, "max": 15},
    {"id": "bedchange", "name": "Bed level change", "cmap": "RdBu",
     "mode": "cells", "dotSize": 6, "opacity": 0.9, "min": -1, "max": 1},
    {"id": "vegetation", "name": "Vegetation", "cmap": "Greens",
     "mode": "cells", "dotSize": 6, "opacity": 0.9, "min": None, "max": None},
    {"id": "grayscale", "name": "Grayscale", "cmap": "gray",
     "mode": "cells", "dotSize": 6, "opacity": 0.9, "min": None, "max": None},
]


def _load_styles():
    data = load_json(settings.STYLES_FILE, default=None)
    if not isinstance(data, dict) or not isinstance(data.get("styles"), list):
        return [dict(s) for s in DEFAULT_STYLES]
    return data["styles"]


@route("GET", "/api/styles")
def _get_styles(handler, query, tail):
    send_json(handler, {"styles": _load_styles()})


@route("POST", "/api/styles")
def _set_styles(handler, body, tail):
    styles = (body or {}).get("styles")
    if not isinstance(styles, list):
        send_error_json(handler, "missing 'styles' list")
        return
    cleaned = []
    for s in styles:
        if not isinstance(s, dict) or not s.get("id") or not s.get("name"):
            send_error_json(handler, "every style needs an 'id' and a 'name'")
            return
        cleaned.append({
            "id": str(s["id"]), "name": str(s["name"]),
            "cmap": str(s.get("cmap") or "viridis"),
            "mode": "dots" if s.get("mode") == "dots" else "cells",
            "dotSize": float(s.get("dotSize") or 6),
            "opacity": float(s["opacity"]) if s.get("opacity") is not None else 0.9,
            "min": None if s.get("min") is None else float(s["min"]),
            "max": None if s.get("max") is None else float(s["max"]),
        })
    settings.APP_DIR.mkdir(parents=True, exist_ok=True)
    save_json(settings.STYLES_FILE, {"styles": cleaned})
    send_json(handler, {"ok": True, "count": len(cleaned)})
