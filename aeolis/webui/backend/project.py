"""Project (model folder) management for the AeoLiS web GUI.

A "project" is simply the folder that contains an ``aeolis.txt``
configuration file. All GUI-specific artifacts live in a single ``gui``
subfolder so the whole project remains one portable directory:

    <project>/
      aeolis.txt            model configuration
      *.grd, wind.txt, ...  model input files
      aeolis.nc, aeolis.log run output
      gui/
        project.json        GUI state (crs, layers, grid draft, ...)
        polygons.json       objects store (polygons, transects, boxes)
        rawdata/            downloaded/imported raw data + manifest.json
        cache/              derived binary caches (safe to delete)
"""

import threading
import time
from pathlib import Path

from aeolis.webui.backend import settings
from aeolis.webui.backend.util import load_json, save_json

_LOCK = threading.Lock()
_current = None


class Project:

    def __init__(self, configfile):
        self.configfile = Path(configfile).resolve()
        self.root = self.configfile.parent
        self.gui_dir = self.root / settings.GUI_DIRNAME
        self.rawdata_dir = self.gui_dir / settings.RAWDATA_DIRNAME
        self.cache_dir = self.gui_dir / settings.CACHE_DIRNAME
        self.ensure_structure()

    def ensure_structure(self):
        self.rawdata_dir.mkdir(parents=True, exist_ok=True)
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    # --- persisted GUI state ------------------------------------------

    @property
    def state_file(self):
        return self.gui_dir / settings.PROJECT_FILE

    @property
    def polygons_file(self):
        return self.gui_dir / settings.POLYGONS_FILE

    def load_state(self):
        return load_json(self.state_file, default={})

    def save_state(self, state):
        save_json(self.state_file, state)

    def load_polygons(self):
        return load_json(self.polygons_file, default={"features": []})

    def save_polygons(self, obj):
        save_json(self.polygons_file, obj)

    def info(self):
        return {
            "configfile": str(self.configfile),
            "root": str(self.root),
            "name": self.root.name,
            "exists": self.configfile.is_file(),
        }


def current():
    with _LOCK:
        return _current


def require():
    project = current()
    if project is None:
        raise RuntimeError("no project is open")
    return project


def open_project(path):
    """Open a project from an aeolis.txt path or a project folder."""
    global _current
    path = Path(path).expanduser().resolve()
    if path.is_dir():
        configfile = path / "aeolis.txt"
    else:
        configfile = path
    with _LOCK:
        _current = Project(configfile)
        project = _current
    _remember_recent(project)
    return project


def new_project(folder):
    """Create a fresh project (empty aeolis.txt) in *folder*."""
    folder = Path(folder).expanduser().resolve()
    folder.mkdir(parents=True, exist_ok=True)
    configfile = folder / "aeolis.txt"
    if not configfile.exists():
        from aeolis.inout import write_configfile
        write_configfile(str(configfile), None)
    return open_project(configfile)


def _remember_recent(project):
    recent = load_json(settings.RECENT_FILE, default=[])
    entry = {"configfile": str(project.configfile), "name": project.root.name}
    recent = [r for r in recent if r.get("configfile") != entry["configfile"]]
    entry["opened"] = time.strftime("%Y-%m-%d %H:%M:%S")
    recent.insert(0, entry)
    save_json(settings.RECENT_FILE, recent[:10])


def recent_projects():
    """Recent projects whose config file still exists (pruned)."""
    recent = load_json(settings.RECENT_FILE, default=[])
    kept = []
    for entry in recent:
        if Path(entry.get("configfile", "")).is_file():
            entry["exists"] = True
            kept.append(entry)
    if len(kept) != len(recent):
        save_json(settings.RECENT_FILE, kept)
    return kept[:10]
