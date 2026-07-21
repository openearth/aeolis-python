"""Model execution for the Run tab.

A small runner-backend abstraction so the GUI can later submit to the
Deltares HPC cluster (SSH + sbatch) through the same interface; today
there is one implementation: ``LocalRunner`` executes the model as a
subprocess of the GUI server.

The subprocess runs ``AeoLiSRunner`` directly (not the ``aeolis``
console script) so it works in any environment where the aeolis package
is importable. Progress is parsed from the model's own logging
(``print_progress`` lines); the log is kept in a ring buffer served
incrementally by offset.
"""

import os
import re
import subprocess
import sys
import threading
import time

# "010.5%   0:01:23 / 0:10:00 / 0:08:37 / 1.0"
PROGRESS_RE = re.compile(
    r"(\d{1,3}\.\d)%\s+(\S+)\s*/\s*(\S+)\s*/\s*(\S+)\s*/\s*([\d.]+)"
)
MAX_LOG_LINES = 8000


class RunnerBackend:
    """Interface for run backends (local subprocess, HPC later)."""

    id = "abstract"
    title = "Abstract"
    available = True
    note = ""

    def start(self, project):
        raise NotImplementedError

    def stop(self):
        raise NotImplementedError

    def status(self):
        raise NotImplementedError

    def log_tail(self, offset):
        raise NotImplementedError


class LocalRunner(RunnerBackend):

    id = "local"
    title = "This computer"

    def __init__(self):
        self._proc = None
        self._lock = threading.Lock()
        self._lines = []
        self._state = "idle"      # idle | running | finished | error | stopped
        self._progress = {}
        self._started = None
        self._exit_code = None
        self._configfile = None

    # --- lifecycle ----------------------------------------------------

    def start(self, project):
        with self._lock:
            if self._proc is not None and self._proc.poll() is None:
                raise RuntimeError("a simulation is already running")
            self._lines = []
            self._progress = {}
            self._state = "running"
            self._started = time.time()
            self._exit_code = None
            self._configfile = str(project.configfile)

            env = dict(os.environ)
            # the GUI sets this to read files being written; the model
            # itself must use normal locking
            env.pop("HDF5_USE_FILE_LOCKING", None)
            env["PYTHONUNBUFFERED"] = "1"

            code = (
                "from aeolis.model import AeoLiSRunner; "
                f"AeoLiSRunner(configfile=r'{project.configfile}').run()"
            )
            self._proc = subprocess.Popen(
                [sys.executable, "-u", "-c", code],
                cwd=str(project.root),
                env=env,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
                creationflags=getattr(subprocess, "CREATE_NEW_PROCESS_GROUP", 0),
            )
        threading.Thread(target=self._reader, daemon=True, name="aeolis-run-reader").start()

    def _reader(self):
        proc = self._proc
        buffer = ""
        while True:
            chunk = proc.stdout.read(1)
            if chunk == "" and proc.poll() is not None:
                break
            if chunk == "":
                continue
            if chunk in ("\n", "\r"):
                if buffer.strip():
                    self._append(buffer.rstrip())
                buffer = ""
            else:
                buffer += chunk
        if buffer.strip():
            self._append(buffer.rstrip())

        exit_code = proc.wait()
        with self._lock:
            self._exit_code = exit_code
            if self._state == "running":
                self._state = "finished" if exit_code == 0 else "error"
            if self._state == "finished":
                self._progress["percent"] = 100.0

    def _append(self, line):
        with self._lock:
            self._lines.append(line)
            if len(self._lines) > MAX_LOG_LINES:
                del self._lines[: len(self._lines) - MAX_LOG_LINES]
        match = PROGRESS_RE.search(line)
        if match:
            with self._lock:
                self._progress = {
                    "percent": float(match.group(1)),
                    "elapsed": match.group(2),
                    "total": match.group(3),
                    "remaining": match.group(4),
                    "avg_dt": float(match.group(5)),
                }

    def stop(self):
        with self._lock:
            proc = self._proc
            if proc is None or proc.poll() is not None:
                return False
            self._state = "stopped"
        if os.name == "nt":
            subprocess.run(
                ["taskkill", "/PID", str(proc.pid), "/T", "/F"],
                capture_output=True, check=False,
            )
        else:
            proc.terminate()
            try:
                proc.wait(timeout=5)
            except subprocess.TimeoutExpired:
                proc.kill()
        return True

    # --- reporting ----------------------------------------------------

    def status(self):
        with self._lock:
            running = self._proc is not None and self._proc.poll() is None
            return {
                "backend": self.id,
                "state": self._state,
                "running": running,
                "progress": dict(self._progress),
                "started": self._started,
                "wall_elapsed": time.time() - self._started if self._started else None,
                "exit_code": self._exit_code,
                "configfile": self._configfile,
                "log_length": len(self._lines),
            }

    def log_tail(self, offset):
        with self._lock:
            offset = max(0, min(int(offset), len(self._lines)))
            return {"offset": len(self._lines), "lines": self._lines[offset:]}


class HpcRunner(RunnerBackend):
    """Placeholder for the Deltares HPC (SLURM) backend.

    Planned flow: paramiko SSH to the login node, rsync/sftp the project
    folder, generate an sbatch script, submit with ``sbatch``, poll
    ``squeue``/``sacct`` by job id, and sync results back. The Run tab
    already exposes the backend selector so this can slot in without UI
    changes.
    """

    id = "hpc"
    title = "Deltares HPC (Deltareken)"
    available = False
    note = "Planned - runs will be submitted over SSH via sbatch."

    def start(self, project):
        raise RuntimeError("HPC backend not implemented yet")

    def stop(self):
        return False

    def status(self):
        return {"backend": self.id, "state": "unavailable"}

    def log_tail(self, offset):
        return {"offset": 0, "lines": []}


BACKENDS = {b.id: b for b in [LocalRunner(), HpcRunner()]}
active = BACKENDS["local"]


def select(backend_id):
    global active
    backend = BACKENDS.get(backend_id)
    if backend is None:
        raise ValueError(f"unknown backend '{backend_id}'")
    if not backend.available:
        raise RuntimeError(f"backend '{backend_id}' is not available yet")
    active = backend
    return backend
