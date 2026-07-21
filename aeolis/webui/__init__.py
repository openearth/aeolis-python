"""AeoLiS web-based graphical user interface.

This subpackage implements the modern AeoLiS GUI: a zero-build web
application (vanilla JS + MapLibre GL) served by a Python standard
library HTTP server and presented in a native desktop window via
pywebview (with a plain-browser fallback).

It is fully independent from the legacy Tkinter GUI in ``aeolis.gui``.

Launch from the command line::

    aeolis webui [configfile] [--port PORT] [--browser]

or programmatically::

    from aeolis.webui.launcher import launch
    launch()

This file is part of AeoLiS, distributed under the GNU General Public
License v3. See LICENSE.txt for details.
"""

__all__ = ["launch"]

from aeolis.webui.launcher import launch
