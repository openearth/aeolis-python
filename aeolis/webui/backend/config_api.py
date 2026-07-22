"""aeolis.txt load/save API for the Settings tab.

Reading and writing go through ``aeolis.inout`` so the GUI behaves
exactly like the model: values are merged over ``DEFAULT_CONFIG`` on
read, and only non-default values are written back.
"""

import os

import numpy as np

import aeolis.inout
from aeolis.constants import DEFAULT_CONFIG
from aeolis.webui.backend import project
from aeolis.webui.backend.httpd import route
from aeolis.webui.backend.util import send_error_json, send_json


def _jsonable(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    return value


def load_config(configfile):
    """Merged config values (defaults + file), JSON-safe, filenames as
    strings (files are not parsed into arrays)."""
    values = aeolis.inout.read_configfile(
        str(configfile), parse_files=False, load_defaults=True
    )
    return {k: _jsonable(v) for k, v in values.items()}


def explicit_keys(configfile):
    """Keys explicitly present in the config file text."""
    keys = []
    with open(configfile, "r") as fp:
        for line in fp:
            if "=" in line and not line.strip().startswith("%"):
                keys.append(line.split("=")[0].strip())
    return keys


@route("GET", "/api/config")
def _get(handler, query, tail):
    current = project.require()
    if not current.configfile.is_file():
        send_error_json(handler, f"config not found: {current.configfile}", 404)
        return
    send_json(handler, {
        "path": str(current.configfile),
        "values": load_config(current.configfile),
        "explicit": explicit_keys(current.configfile),
    })


@route("GET", "/api/config/raw")
def _raw(handler, query, tail):
    """The aeolis.txt file text as-is (for the config viewer popup)."""
    current = project.require()
    if not current.configfile.is_file():
        send_error_json(handler, f"config not found: {current.configfile}", 404)
        return
    text = current.configfile.read_text(encoding="utf-8", errors="replace")
    send_json(handler, {"path": str(current.configfile), "text": text})


def _prepare_values(values):
    """JSON payload -> config dict suitable for write_configfile."""
    prepared = {}
    for key, value in values.items():
        default = DEFAULT_CONFIG.get(key)
        if value is None:
            prepared[key] = None
        elif isinstance(default, bool):
            prepared[key] = bool(value)
        elif isinstance(value, list):
            prepared[key] = np.asarray(value) if value and isinstance(value[0], (int, float)) else value
        else:
            prepared[key] = value
    return prepared


@route("POST", "/api/config/save")
def _save(handler, body, tail):
    current = project.require()
    values = body.get("values")
    if not isinstance(values, dict):
        send_error_json(handler, "missing 'values'")
        return
    target = body.get("path") or str(current.configfile)
    target = os.path.abspath(target)
    aeolis.inout.write_configfile(target, _prepare_values(values))
    if body.get("path"):
        # "save as" re-anchors the project on the new file
        project.open_project(target)
    send_json(handler, {"ok": True, "path": target})
