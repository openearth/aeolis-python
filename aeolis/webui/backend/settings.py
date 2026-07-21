"""Constants and paths for the AeoLiS web GUI backend."""

from pathlib import Path

# Package directories
WEBUI_DIR = Path(__file__).resolve().parent.parent
WEB_DIR = WEBUI_DIR / "web"

# Per-user application directory (recent projects, global preferences)
APP_DIR = Path.home() / ".aeolis_webui"
RECENT_FILE = APP_DIR / "recent.json"

# Server
HOST = "127.0.0.1"
BASE_PORT = 8790
PORT_TRIES = 21

# Project folder convention: everything the GUI creates lives in a
# single "gui" subfolder next to aeolis.txt so a project stays portable.
GUI_DIRNAME = "gui"
RAWDATA_DIRNAME = "rawdata"
CACHE_DIRNAME = "cache"
PROJECT_FILE = "project.json"
POLYGONS_FILE = "polygons.json"
MANIFEST_FILE = "manifest.json"

# Window
WINDOW_TITLE = "AeoLiS"
WINDOW_SIZE = (1520, 920)
