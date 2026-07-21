"""Shared helpers for the AeoLiS web GUI backend."""

import json
import threading
from pathlib import Path

import numpy as np

# netCDF4/HDF5 is not thread-safe: every open/read/close of any netCDF
# file anywhere in the (threaded) backend must hold this lock.
NC_LOCK = threading.Lock()

MIME = {
    ".html": "text/html; charset=utf-8",
    ".js": "application/javascript; charset=utf-8",
    ".mjs": "application/javascript; charset=utf-8",
    ".css": "text/css; charset=utf-8",
    ".json": "application/json; charset=utf-8",
    ".png": "image/png",
    ".jpg": "image/jpeg",
    ".svg": "image/svg+xml",
    ".ico": "image/x-icon",
    ".map": "application/json",
    ".woff2": "font/woff2",
}


def content_type(path):
    return MIME.get(Path(path).suffix.lower(), "application/octet-stream")


def is_within(path, root):
    """True if *path* resolves to a location inside *root*."""
    try:
        Path(path).resolve().relative_to(Path(root).resolve())
        return True
    except ValueError:
        return False


def nanfilled(var):
    """Read a netCDF variable slice as float32 with masked/fill -> NaN."""
    arr = var[:]
    if np.ma.isMaskedArray(arr):
        arr = arr.filled(np.nan)
    arr = np.asarray(arr, dtype=np.float32)
    return arr


def send_json(handler, obj, status=200):
    body = json.dumps(obj).encode("utf-8")
    handler.send_response(status)
    handler.send_header("Content-Type", "application/json; charset=utf-8")
    handler.send_header("Content-Length", str(len(body)))
    handler.send_header("Cache-Control", "no-store")
    handler.end_headers()
    handler.wfile.write(body)


def send_error_json(handler, message, status=400):
    send_json(handler, {"error": str(message)}, status=status)


def send_bytes(handler, data, ctype="application/octet-stream", extra_headers=None):
    handler.send_response(200)
    handler.send_header("Content-Type", ctype)
    handler.send_header("Content-Length", str(len(data)))
    handler.send_header("Cache-Control", "no-store")
    for key, value in (extra_headers or {}).items():
        handler.send_header(key, value)
    handler.end_headers()
    handler.wfile.write(data)


def read_body_json(handler):
    length = int(handler.headers.get("Content-Length") or 0)
    if length == 0:
        return {}
    raw = handler.rfile.read(length)
    return json.loads(raw.decode("utf-8"))


def load_json(path, default=None):
    path = Path(path)
    if not path.is_file():
        return default
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return default


def save_json(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2), encoding="utf-8")
    tmp.replace(path)
