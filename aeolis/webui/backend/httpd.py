"""HTTP server for the AeoLiS web GUI.

Stdlib-only threaded HTTP server. API modules register their routes with
the :func:`route` decorator; static frontend files are served from
``aeolis/webui/web``. Binary payloads (grid meshes, netCDF field slabs)
are served as ``application/octet-stream``.
"""

from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import parse_qs, unquote, urlparse

from aeolis.webui.backend import settings
from aeolis.webui.backend.util import (
    content_type,
    is_within,
    read_body_json,
    send_error_json,
    send_json,
)

# --- route registry ----------------------------------------------------
# exact:  {("GET", "/api/ping"): fn}
# prefix: [("GET", "/api/job/", fn)]  (longest prefix wins)
_EXACT = {}
_PREFIX = []


def route(method, path, prefix=False):
    """Register ``fn(handler, query_or_body, tail)`` for an API path."""

    def _register(fn):
        if prefix:
            _PREFIX.append((method, path, fn))
            _PREFIX.sort(key=lambda item: len(item[1]), reverse=True)
        else:
            _EXACT[(method, path)] = fn
        return fn

    return _register


def _load_api_modules():
    """Import all API modules so their routes register."""
    from aeolis.webui.backend import project_api  # noqa: F401
    from aeolis.webui.backend import schema_api  # noqa: F401
    from aeolis.webui.backend import config_api  # noqa: F401
    from aeolis.webui.backend import grid_api  # noqa: F401
    from aeolis.webui.backend import objects_api  # noqa: F401
    from aeolis.webui.backend import domain_api  # noqa: F401
    from aeolis.webui.backend import conditions_api  # noqa: F401
    from aeolis.webui.backend import run_api  # noqa: F401
    from aeolis.webui.backend import output_api  # noqa: F401
    from aeolis.webui.backend import jobs_api  # noqa: F401
    from aeolis.webui.backend import dialogs  # noqa: F401


class GuiHandler(BaseHTTPRequestHandler):
    protocol_version = "HTTP/1.1"

    # --- logging: keep the console quiet -----------------------------
    def log_message(self, fmt, *args):  # noqa: A003
        pass

    # --- request handling --------------------------------------------
    def do_GET(self):  # noqa: N802
        parsed = urlparse(self.path)
        path = unquote(parsed.path)
        query = {k: v[-1] for k, v in parse_qs(parsed.query).items()}
        try:
            if path.startswith("/api/"):
                self._dispatch("GET", path, query)
            else:
                self._send_static(path)
        except (ConnectionAbortedError, ConnectionResetError, BrokenPipeError):
            pass
        except Exception as exc:  # noqa: BLE001 - reported to the client
            self._safe_error(exc)

    def do_POST(self):  # noqa: N802
        parsed = urlparse(self.path)
        path = unquote(parsed.path)
        try:
            body = read_body_json(self)
            self._dispatch("POST", path, body)
        except (ConnectionAbortedError, ConnectionResetError, BrokenPipeError):
            pass
        except Exception as exc:  # noqa: BLE001 - reported to the client
            self._safe_error(exc)

    def _dispatch(self, method, path, payload):
        fn = _EXACT.get((method, path))
        if fn is not None:
            fn(self, payload, "")
            return
        for pmethod, prefix, pfn in _PREFIX:
            if pmethod == method and path.startswith(prefix):
                pfn(self, payload, path[len(prefix):])
                return
        send_error_json(self, f"unknown endpoint {method} {path}", status=404)

    def _safe_error(self, exc):
        try:
            send_error_json(self, exc, status=500)
        except (ConnectionAbortedError, ConnectionResetError, BrokenPipeError, OSError):
            pass

    # --- static files -------------------------------------------------
    def _send_static(self, path):
        if path in ("", "/"):
            path = "/index.html"
        target = (settings.WEB_DIR / path.lstrip("/")).resolve()
        if not is_within(target, settings.WEB_DIR) or not target.is_file():
            send_error_json(self, f"not found: {path}", status=404)
            return
        data = target.read_bytes()
        self.send_response(200)
        self.send_header("Content-Type", content_type(target))
        self.send_header("Content-Length", str(len(data)))
        # no-store for app code so the edit->F5 dev loop always works;
        # vendored libs are fingerprinted by version and may cache.
        if "/vendor/" in path:
            self.send_header("Cache-Control", "max-age=86400")
        else:
            self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(data)


class QuietServer(ThreadingHTTPServer):
    daemon_threads = True
    # fail loudly instead of silently sharing a stale port (Windows)
    allow_reuse_address = False

    def handle_error(self, request, client_address):
        import sys
        exc = sys.exc_info()[1]
        if isinstance(exc, (ConnectionAbortedError, ConnectionResetError, BrokenPipeError)):
            return
        super().handle_error(request, client_address)


@route("GET", "/api/ping")
def _ping(handler, query, tail):
    send_json(handler, {"ok": True, "app": "aeolis-webui"})


def make_server(host, base_port=None, tries=1):
    """Bind the server on the first free port starting at *base_port*."""
    _load_api_modules()
    base_port = settings.BASE_PORT if base_port is None else base_port
    last_error = None
    for offset in range(max(1, tries)):
        try:
            return QuietServer((host, base_port + offset), GuiHandler)
        except OSError as exc:
            last_error = exc
    raise OSError(
        f"no free port in range {base_port}-{base_port + tries - 1}: {last_error}"
    )
