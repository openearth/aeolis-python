"""Launcher for the AeoLiS web GUI.

Starts the backend HTTP server on a background thread and presents the
frontend in a native desktop window (pywebview, titled "AeoLiS", no URL
bar). Falls back to the default browser when pywebview is not installed
or when ``--browser`` is requested.
"""

import faulthandler
import logging
import os
import threading
import webbrowser

from aeolis.webui.backend import settings

logger = logging.getLogger(__name__)


def _prepare_environment():
    # Allow the GUI process to read netCDF files that are still being
    # written by a running model (Windows/HDF5 file locking). The run
    # subprocess gets this variable stripped again in run_manager.
    os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
    faulthandler.enable()


def _start_server(port=None):
    from aeolis.webui.backend.httpd import make_server

    server = make_server(
        settings.HOST,
        base_port=port or settings.BASE_PORT,
        tries=1 if port else settings.PORT_TRIES,
    )
    thread = threading.Thread(
        target=server.serve_forever, daemon=True, name="aeolis-webui-httpd"
    )
    thread.start()
    return server


def _open_native_window(url):
    """Open the GUI in a native window via pywebview. Returns False if
    pywebview is unavailable (caller falls back to the browser)."""
    try:
        import webview
    except ImportError:
        return False

    width, height = settings.WINDOW_SIZE
    webview.create_window(
        settings.WINDOW_TITLE,
        url,
        width=width,
        height=height,
        min_size=(1024, 640),
        confirm_close=False,
        text_select=True,
    )
    webview.start()  # blocks until the window is closed
    return True


def launch(configfile=None, port=None, browser=False):
    """Start the AeoLiS web GUI.

    Parameters
    ----------
    configfile : str, optional
        aeolis.txt to open as project on startup.
    port : int, optional
        Fixed port; default scans from ``settings.BASE_PORT``.
    browser : bool
        Force a browser tab instead of the native window.
    """
    _prepare_environment()

    if configfile:
        from aeolis.webui.backend import project

        opened = project.open_project(configfile)
        if not opened.configfile.is_file():
            logger.warning("configuration file not found: %s", configfile)

    server = _start_server(port=port)
    host, bound_port = server.server_address[:2]
    url = f"http://{host}:{bound_port}/"
    print(f"AeoLiS GUI serving at {url}")

    try:
        if browser:
            webbrowser.open(url)
            _wait_forever()
        elif not _open_native_window(url):
            print(
                "pywebview not installed (pip install aeolis[webui]) - "
                "opening in the default browser instead."
            )
            webbrowser.open(url)
            _wait_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.shutdown()


def _wait_forever():
    """Keep the process alive while the server runs (browser mode)."""
    event = threading.Event()
    try:
        event.wait()
    except KeyboardInterrupt:
        raise


if __name__ == "__main__":
    launch()
