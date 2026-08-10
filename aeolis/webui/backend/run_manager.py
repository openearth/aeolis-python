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
import posixpath
import re
import shlex
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


# ---------------------------------------------------------------------
# Deltares HYDRAX (HAL8, ex-h7) SLURM backend
#
# Deltareken uses SLURM; login node hal8.directory.intra with AD
# username+password (no MFA). The /p project filesystem is mounted both
# on the workstation and the cluster, so nothing is uploaded except the
# tiny job script: we SSH in, write <run_dir>/<job>.sh, sbatch it, then
# poll squeue/sacct and tail the job's .o<id> log (whose AeoLiS progress
# lines feed the same Run-tab progress bar as a local run).
# ---------------------------------------------------------------------

HPC_PARTITIONS = ["1vcpu", "4vcpu", "16vcpu", "24vcpu", "44vcpu", "60vcpu", "gpu", "test"]

DEFAULT_HPC_PROFILE = {
    "host": "hal8.directory.intra",
    "user": "",
    "job_name": "aeolis",
    "partition": "1vcpu",
    "ntasks": 1,
    "cpus_per_task": 1,
    "walltime": "5-00:00:00",         # days-hours:min:sec (mandatory on HYDRAX)
    "modules": ["miniforge/latest"],
    "conda_setup": "/opt/miniforge3/etc/profile.d/conda.sh",
    "env_path": "",                    # conda activate target (user-selected)
    "run_dir": "",                     # /p/... folder holding the config (user-selected)
    "config": "aeolis.txt",
    "mail_user": "",
    "extra_sbatch": [],
    "run_mode": "inplace",             # "inplace" (project already on /p) | "copy"
}


def local_to_linux(path):
    """Windows mount path -> cluster path, e.g. P:\\proj\\run -> /p/proj/run.
    The drive letter maps to a top-level /<letter> (the /p project drive)."""
    m = re.match(r"^([A-Za-z]):[\\/](.*)$", str(path or ""))
    if m:
        return "/" + m.group(1).lower() + "/" + m.group(2).replace("\\", "/")
    return str(path or "").replace("\\", "/")


def linux_to_local(path):
    """Cluster path -> Windows mount path, e.g. /p/proj/run -> P:\\proj\\run."""
    parts = str(path or "").strip("/").split("/")
    if parts and len(parts[0]) == 1 and parts[0].isalpha():
        return parts[0].upper() + ":\\" + "\\".join(parts[1:])
    return str(path or "")


def build_job_script(profile):
    """Render a SLURM batch script from a profile. Pure + deterministic
    so it can be previewed in the UI and unit-tested."""
    p = {**DEFAULT_HPC_PROFILE, **(profile or {})}
    name = p.get("job_name") or "aeolis"
    lines = [
        "#!/bin/bash",
        f"#SBATCH --job-name={name}",
        f"#SBATCH --output={name}.o%j",
        f"#SBATCH --partition={p['partition']}",
        f"#SBATCH --ntasks={int(p['ntasks'])}",
        f"#SBATCH --cpus-per-task={int(p['cpus_per_task'])}",
        f"#SBATCH --time={p['walltime']}",
    ]
    if p.get("mail_user"):
        lines.append(f"#SBATCH --mail-user={p['mail_user']}")
        lines.append("#SBATCH --mail-type=BEGIN,END,FAIL")
    for extra in p.get("extra_sbatch") or []:
        extra = str(extra).strip()
        if extra:
            lines.append(extra if extra.startswith("#SBATCH") else f"#SBATCH {extra}")
    lines.append("")
    for mod in p.get("modules") or []:
        if str(mod).strip():
            lines.append(f"module load {str(mod).strip()}")
    if p.get("conda_setup"):
        lines.append(f"source {p['conda_setup']}")
    if p.get("env_path"):
        lines.append(f"conda activate {p['env_path']}")
    if p.get("run_dir"):
        lines.append(f"cd {p['run_dir']}")
    lines.append(f"aeolis run ./{p.get('config') or 'aeolis.txt'}")
    return "\n".join(lines) + "\n"


def parse_job_id(text):
    """Extract the job id from sbatch's 'Submitted batch job 12345'."""
    match = re.search(r"Submitted batch job (\d+)", text or "")
    return match.group(1) if match else None


def parse_squeue(text):
    """Parse one row of `squeue -h -o '%T|%r|%M|%D|%P'` (state, reason,
    time, nodes, partition). Returns None if the job is no longer queued."""
    for line in (text or "").splitlines():
        parts = line.strip().split("|")
        if len(parts) >= 5 and parts[0]:
            return {"state": parts[0], "reason": parts[1], "time": parts[2],
                    "nodes": parts[3], "partition": parts[4]}
    return None


def parse_sacct(text, job_id):
    """Parse the main accounting row of `sacct -n -P -o
    JobID,State,ExitCode,Elapsed` for a finished job."""
    for line in (text or "").splitlines():
        parts = line.strip().split("|")
        if len(parts) >= 4 and parts[0] == str(job_id):
            return {"state": parts[1].split()[0], "exit_code": parts[2], "elapsed": parts[3]}
    return None


class HpcRunner(RunnerBackend):
    """Deltares HYDRAX (SLURM) backend over SSH (paramiko)."""

    id = "hpc"
    title = "Deltares HPC (HYDRAX)"

    POLL_SECONDS = 12

    def __init__(self):
        self._lock = threading.Lock()
        self.profile = dict(DEFAULT_HPC_PROFILE)
        self._password = None
        self._script = None          # optional raw-script override
        self._lines = []
        self._state = "idle"         # idle|submitting|queued|running|finished|error|stopped
        self._progress = {}
        self._started = None
        self._job_id = None
        self._remote = {}            # last squeue/sacct fields
        self._log_bytes = 0
        self._stop_flag = False
        self._poller = None

    # --- capability (paramiko may be absent) --------------------------

    @property
    def available(self):
        try:
            import paramiko  # noqa: F401
            return True
        except Exception:  # noqa: BLE001
            return False

    @property
    def note(self):
        return "" if self.available else "run: pip install paramiko"

    def configure(self, profile=None, password=None, script=None):
        if profile is not None:
            self.profile = {**DEFAULT_HPC_PROFILE, **profile}
        if password is not None:
            self._password = password
        self._script = script or None

    # --- SSH helpers --------------------------------------------------

    def _connect(self):
        import paramiko
        client = paramiko.SSHClient()
        client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        client.connect(
            self.profile["host"], username=self.profile["user"],
            password=self._password, look_for_keys=False, allow_agent=False,
            timeout=20,
        )
        return client

    @staticmethod
    def _run(client, cmd):
        stdin, stdout, stderr = client.exec_command(cmd, timeout=45)
        out = stdout.read().decode("utf-8", "replace")
        err = stderr.read().decode("utf-8", "replace")
        code = stdout.channel.recv_exit_status()
        return out, err, code

    @staticmethod
    def _copy_project(src, dest):
        """Copy the project inputs to ``dest`` (a mounted /p folder),
        skipping the GUI cache and previous run outputs."""
        import shutil
        src = Path(src)
        dest = Path(dest)
        dest.mkdir(parents=True, exist_ok=True)

        def _ignore(_dir, names):
            return [n for n in names
                    if n in ("gui", ".git") or n.endswith((".nc", ".log"))]

        shutil.copytree(src, dest, ignore=_ignore, dirs_exist_ok=True)

    @staticmethod
    def _write_remote(client, path, content):
        sftp = client.open_sftp()
        try:
            with sftp.open(path, "w") as fh:
                fh.write(content)
        finally:
            sftp.close()

    def _tail_remote(self, client, path):
        try:
            sftp = client.open_sftp()
            try:
                with sftp.open(path, "r") as fh:
                    fh.seek(self._log_bytes)
                    data = fh.read()
            finally:
                sftp.close()
        except IOError:
            return   # output file not created yet
        if not data:
            return
        self._log_bytes += len(data)
        for line in data.decode("utf-8", "replace").splitlines():
            if line.strip():
                self._ingest(line.rstrip())

    # --- lifecycle ----------------------------------------------------

    def start(self, project):
        if not self.available:
            raise RuntimeError("paramiko is not installed (pip install paramiko)")
        p = self.profile
        missing = [k for k in ("host", "user", "run_dir", "env_path") if not p.get(k)]
        if missing:
            raise RuntimeError("HPC settings incomplete: " + ", ".join(missing))
        if not self._password:
            raise RuntimeError("no password provided for the HPC connection")
        if self._poller and self._poller.is_alive():
            raise RuntimeError("a simulation is already submitted")

        script = self._script or build_job_script(p)
        with self._lock:
            self._lines = []
            self._progress = {}
            self._remote = {}
            self._job_id = None
            self._log_bytes = 0
            self._stop_flag = False
            self._state = "submitting"
            self._started = time.time()

        # optionally copy the project onto the (mounted) /p run dir first
        if p.get("run_mode") == "copy":
            dest = Path(linux_to_local(p["run_dir"]))
            self._ingest(f"Copying project to {dest} …")
            try:
                self._copy_project(project.root, dest)
            except Exception as exc:  # noqa: BLE001
                with self._lock:
                    self._state = "error"
                self._ingest(f"copy failed: {exc}")
                raise RuntimeError(f"copying the project failed: {exc}")

        self._ingest(f"Connecting to {p['user']}@{p['host']} …")
        try:
            client = self._connect()
        except Exception as exc:  # noqa: BLE001
            with self._lock:
                self._state = "error"
            self._ingest(f"connection failed: {exc}")
            raise RuntimeError(f"SSH connection failed: {exc}")

        try:
            name = p.get("job_name") or "aeolis"
            remote_sh = posixpath.join(p["run_dir"], f"{name}.sh")
            self._write_remote(client, remote_sh, script)
            self._ingest(f"Wrote {remote_sh}")
            out, err, code = self._run(
                client, f"cd {shlex.quote(p['run_dir'])} && sbatch {shlex.quote(name + '.sh')}")
            if code != 0:
                raise RuntimeError((err or out).strip() or "sbatch failed")
            job_id = parse_job_id(out)
            if not job_id:
                raise RuntimeError(f"could not read job id from: {out.strip()}")
            with self._lock:
                self._job_id = job_id
                self._state = "queued"
            self._ingest(f"Submitted batch job {job_id}")
        except Exception as exc:  # noqa: BLE001
            with self._lock:
                self._state = "error"
            self._ingest(str(exc))
            raise RuntimeError(str(exc))
        finally:
            client.close()

        self._poller = threading.Thread(target=self._poll_loop, daemon=True,
                                        name="aeolis-hpc-poll")
        self._poller.start()

    def _poll_loop(self):
        p = self.profile
        name = p.get("job_name") or "aeolis"
        out_file = posixpath.join(p["run_dir"], f"{name}.o{self._job_id}")
        while not self._stop_flag:
            client = None
            try:
                client = self._connect()
                self._tail_remote(client, out_file)
                q_out, _, _ = self._run(
                    client, f"squeue -j {self._job_id} -h -o '%T|%r|%M|%D|%P'")
                info = parse_squeue(q_out)
                if info:
                    with self._lock:
                        self._remote = info
                        self._state = ("running"
                                       if info["state"] in ("RUNNING", "COMPLETING")
                                       else "queued")
                else:
                    a_out, _, _ = self._run(
                        client, f"sacct -j {self._job_id} -n -P -o JobID,State,ExitCode,Elapsed")
                    fin = parse_sacct(a_out, self._job_id)
                    with self._lock:
                        if fin:
                            self._remote = fin
                            st = fin["state"]
                            self._state = ("finished" if st.startswith("COMPLETED")
                                           else "stopped" if st.startswith("CANCELLED")
                                           else "error")
                        else:
                            self._state = "finished"
                        if self._state == "finished":
                            self._progress["percent"] = 100.0
                    self._ingest(f"Job {self._job_id} {self._state}")
                    break
            except Exception as exc:  # noqa: BLE001
                self._ingest(f"poll error: {exc}")
            finally:
                if client is not None:
                    try:
                        client.close()
                    except Exception:  # noqa: BLE001
                        pass
            for _ in range(self.POLL_SECONDS):
                if self._stop_flag:
                    break
                time.sleep(1)

    def stop(self):
        self._stop_flag = True
        if not self._job_id:
            return False
        try:
            client = self._connect()
            try:
                self._run(client, f"scancel {self._job_id}")
            finally:
                client.close()
            with self._lock:
                self._state = "stopped"
            self._ingest(f"Cancelled job {self._job_id}")
            return True
        except Exception as exc:  # noqa: BLE001
            self._ingest(f"scancel failed: {exc}")
            return False

    # --- reporting ----------------------------------------------------

    def _ingest(self, line):
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

    def status(self):
        with self._lock:
            running = self._state in ("submitting", "queued", "running")
            return {
                "backend": self.id,
                "state": self._state,
                "running": running,
                "progress": dict(self._progress),
                "started": self._started,
                "wall_elapsed": time.time() - self._started if self._started else None,
                "exit_code": None,
                "configfile": self.profile.get("config"),
                "log_length": len(self._lines),
                "hpc": {"job_id": self._job_id, **self._remote},
            }

    def log_tail(self, offset):
        with self._lock:
            offset = max(0, min(int(offset), len(self._lines)))
            return {"offset": len(self._lines), "lines": self._lines[offset:]}


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
