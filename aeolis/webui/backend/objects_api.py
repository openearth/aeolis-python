"""Objects store API: persist user-drawn shapes (polygons, boxes,
transects) to ``gui/polygons.json`` in the project folder. The frontend
Objects module is the live authority; this simply loads/saves it."""

from aeolis.webui.backend import project
from aeolis.webui.backend.httpd import route
from aeolis.webui.backend.util import send_json


@route("GET", "/api/objects")
def _load(handler, query, tail):
    data = project.require().load_polygons()
    # stored format: {"objects": [...]} (legacy "features" tolerated)
    objects = data.get("objects", data.get("features", []))
    send_json(handler, {"objects": objects})


@route("POST", "/api/objects/save")
def _save(handler, body, tail):
    objects = body.get("objects", [])
    project.require().save_polygons({"objects": objects})
    send_json(handler, {"ok": True, "count": len(objects)})
