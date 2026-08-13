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
from pathlib import Path

import aeolis.inout

# "010.5%   0:01:23 / 0:10:00 / 0:08:37 / 1.0"
# The time fields are datetime.timedelta reprs, which read "1 day,
# 0:08:37" beyond 24 h - early in a long run the remaining-time estimate
# almost always does, so the fields must be allowed to contain spaces.
PROGRESS_RE = re.compile(
    r"(\d{1,3}\.\d)%\s+(.+?)\s*/\s*(.+?)\s*/\s*(.+?)\s*/\s*([\d.]+)\s*$"
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


def linuxify_config(cfg_path):
    """Rewrite a COPIED aeolis config for the cluster: absolute Windows
    paths become their /<drive> mount form, relative backslash paths get
    forward slashes, and the text is forced to ASCII - the cluster opens
    the config with the Linux default (UTF-8), so a stray Windows-encoded
    character in a comment would abort the run with a UnicodeDecodeError.
    Comments (%) are preserved. The original project's config is never
    touched - only the /p copy."""
    cfg_path = Path(cfg_path)
    if not cfg_path.is_file():
        return 0
    changed = 0
    out_lines = []
    for line in aeolis.inout.read_config_lines(str(cfg_path)):
        line = line.rstrip("\n").rstrip("\r")
        clean = aeolis.inout.to_ascii(line)
        if clean != line:
            changed += 1
            line = clean
        if "=" in line and not line.lstrip().startswith("%"):
            head, rest = line.split("=", 1)
            val, comment = (rest.split("%", 1) + [None])[:2]
            v = val.strip()
            v2 = local_to_linux(v) if re.match(r"^[A-Za-z]:[\\/]", v) else v.replace("\\", "/")
            if v2 != v:
                changed += 1
                line = f"{head}= {v2}" + (f"   %{comment}" if comment is not None else "")
        out_lines.append(line)
    if changed:
        cfg_path.write_text("\n".join(out_lines) + "\n", encoding="utf-8")
    return changed


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


def parse_scontrol_stdout(text):
    """Extract the job's StdOut path from `scontrol show job` output."""
    match = re.search(r"StdOut=(\S+)", text or "")
    return match.group(1) if match else None


def record_hpc_job(project, entry):
    """Remember a submitted SLURM job in the project's GUI state so the
    run can be re-attached after the GUI was closed and reopened. Holds
    no secrets (the password is never persisted)."""
    state = project.load_state()
    jobs = [j for j in state.get("hpc_jobs") or []
            if str(j.get("job_id")) != str(entry.get("job_id"))]
    jobs.insert(0, entry)
    state["hpc_jobs"] = jobs[:10]
    project.save_state(state)


def recorded_hpc_jobs(project):
    """Jobs previously submitted from this project (newest first)."""
    return list(project.load_state().get("hpc_jobs") or [])


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
        self._out_file = None        # remote stdout path being tailed
        self._log_bytes = 0
        self._stop_flag = False
        self._poller = None
        self._gen = 0                # bumped when the monitor switches jobs;
                                     # a poll loop from an older generation exits

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

    def _copy_project(self, src, dest):
        """Copy the project inputs to ``dest`` (a mounted /p folder),
        skipping the GUI cache and previous run outputs. Logs progress
        so long copies to the P-drive don't look stalled."""
        import shutil
        src = Path(src)
        dest = Path(dest)
        dest.mkdir(parents=True, exist_ok=True)
        skip = ("gui", ".git")
        n_files = 0
        n_bytes = 0
        for root, dirs, files in os.walk(src):
            dirs[:] = [d for d in dirs if d not in skip]
            rel = Path(root).relative_to(src)
            (dest / rel).mkdir(parents=True, exist_ok=True)
            for fname in files:
                if self._stop_flag:
                    raise RuntimeError("cancelled")
                if fname.endswith((".nc", ".log")):
                    continue
                fsrc = Path(root) / fname
                shutil.copy2(fsrc, dest / rel / fname)
                n_files += 1
                n_bytes += fsrc.stat().st_size
                if n_files % 100 == 0:
                    self._ingest(f"  … {n_files} files copied ({n_bytes / 1e6:.1f} MB)")
        self._ingest(f"Copy done: {n_files} files, {n_bytes / 1e6:.1f} MB")

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

        with self._lock:
            self._lines = []
            self._progress = {}
            self._remote = {}
            self._job_id = None
            self._log_bytes = 0
            self._stop_flag = False
            self._state = "submitting"
            self._started = time.time()

        # the copy to /p and the SSH handshake can take minutes: run the
        # whole submission in the background so the UI gets its response
        # immediately and can follow every stage in the live log
        self._ingest("Submitting to HYDRAX — follow the stages below.")
        self._poller = threading.Thread(
            target=self._submit_and_poll, args=(project,), daemon=True,
            name="aeolis-hpc-submit")
        self._poller.start()

    def _submit_and_poll(self, project):
        p = self.profile
        gen = self._gen
        script = self._script or build_job_script(p)

        # optionally copy the project onto the (mounted) /p run dir first
        if p.get("run_mode") == "copy":
            dest = Path(linux_to_local(p["run_dir"]))
            self._ingest(f"[1/4] Copying project inputs to {dest} …")
            try:
                self._copy_project(project.root, dest)
                # the cluster cannot read C:\... references - convert the
                # COPY's config to /p-style paths (the original stays as-is)
                n = linuxify_config(dest / (p.get("config") or "aeolis.txt"))
                if n:
                    self._ingest(f"Converted {n} file reference(s) to /p form in the copied config")
            except Exception as exc:  # noqa: BLE001
                with self._lock:
                    self._state = "error"
                self._ingest(f"✖ copying the project failed: {exc}")
                return
        else:
            self._ingest(f"[1/4] Running in place — no copy needed ({p['run_dir']})")

        self._ingest(f"[2/4] Connecting to {p['user']}@{p['host']} (SSH) …")
        try:
            client = self._connect()
        except Exception as exc:  # noqa: BLE001
            with self._lock:
                self._state = "error"
            self._ingest(f"✖ SSH connection failed: {exc}")
            return
        self._ingest("Connected.")

        try:
            name = p.get("job_name") or "aeolis"
            remote_sh = posixpath.join(p["run_dir"], f"{name}.sh")
            self._ingest(f"[3/4] Writing job script {remote_sh} …")
            self._write_remote(client, remote_sh, script)
            self._ingest(f"[4/4] Submitting (sbatch {name}.sh) …")
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
            self._out_file = posixpath.join(p["run_dir"], f"{name}.o{job_id}")
            self._ingest(f"✔ Submitted batch job {job_id} — waiting for SLURM to schedule it.")
            try:
                # remembered per project so the job can be re-attached
                # after the GUI is closed and reopened
                record_hpc_job(project, {
                    "job_id": job_id, "job_name": name,
                    "run_dir": p["run_dir"], "partition": p.get("partition"),
                    "config": p.get("config"), "host": p.get("host"),
                    "user": p.get("user"), "submitted": time.time(),
                })
            except Exception:  # noqa: BLE001 - remembering must not kill the submit
                pass
            if p.get("run_mode") == "copy":
                # the run now writes its output in the /p copy - point the
                # Viewer there so it shows THIS run, not the stale local
                # file (the user can switch back in the Viewer anytime)
                try:
                    from aeolis.webui.backend import output_api
                    output_api.set_source(project, p["run_dir"], p.get("config"))
                    self._ingest("Viewer output now follows the run folder on /p "
                                 "(switch back anytime in the Viewer).")
                except Exception as exc:  # noqa: BLE001 - viewing is optional
                    self._ingest(f"could not point the Viewer at the run folder: {exc}")
        except Exception as exc:  # noqa: BLE001
            with self._lock:
                self._state = "error"
            self._ingest(f"✖ {exc}")
            return
        finally:
            client.close()

        self._poll_loop(gen)

    def _queue_position(self, client, partition):
        """1-based place of our job among the partition's pending jobs
        (None when it cannot be determined)."""
        try:
            out, _, code = self._run(
                client,
                f"squeue -p {shlex.quote(partition)} -t PD -h -o %i --sort=-p,i")
            if code != 0:
                return None, None
            ids = [ln.strip() for ln in out.splitlines() if ln.strip()]
            return ids.index(str(self._job_id)) + 1, len(ids)
        except Exception:  # noqa: BLE001 - incl. ValueError when not listed
            return None, None

    def attach(self, job_id, out_file=None):
        """(Re-)attach the monitor to a submitted SLURM job — after the
        GUI was reopened while the job kept running, or to switch the
        monitor to another job. Attaching to the job that is already
        being monitored is a no-op; attaching to a different one
        abandons the current poll loop (the job itself is untouched —
        unlike stop(), this never scancels) and follows the new job.
        The password must have been configured (it is never persisted)."""
        if not self.available:
            raise RuntimeError("paramiko is not installed (pip install paramiko)")
        if not self._password:
            raise RuntimeError("no password provided for the HPC connection")
        if self._poller and self._poller.is_alive():
            if str(job_id) == str(self._job_id):
                self._ingest(f"Already monitoring job {job_id}.")
                return
            self._gen += 1     # the old poll loop sees this and exits
        with self._lock:
            self._lines = []
            self._progress = {}
            self._remote = {}
            self._job_id = str(job_id)
            self._log_bytes = 0
            self._stop_flag = False
            self._state = "queued"
            self._started = time.time()
        self._out_file = out_file or None
        self._ingest(f"Reattaching to job {job_id} — fetching its log and status …")
        self._poller = threading.Thread(
            target=self._attach_and_poll, args=(self._gen,), daemon=True,
            name="aeolis-hpc-attach")
        self._poller.start()

    def _attach_and_poll(self, gen):
        if not self._out_file:
            # scontrol knows the exact stdout path while the job lives
            try:
                client = self._connect()
                try:
                    out, _, code = self._run(
                        client, f"scontrol show job {self._job_id}")
                    self._out_file = parse_scontrol_stdout(out) if code == 0 else None
                finally:
                    client.close()
            except Exception as exc:  # noqa: BLE001
                with self._lock:
                    self._state = "error"
                self._ingest(f"✖ SSH connection failed: {exc}")
                return
        if not self._out_file:
            p = self.profile
            name = p.get("job_name") or "aeolis"
            self._out_file = posixpath.join(p["run_dir"], f"{name}.o{self._job_id}")
        self._ingest(f"Following {self._out_file}")
        self._poll_loop(gen)

    def list_jobs(self, recorded=None):
        """The user's live SLURM jobs merged with this project's
        recorded submissions; recorded jobs that are no longer queued
        get their final state from sacct."""
        client = self._connect()
        try:
            out, err, code = self._run(
                client,
                f"squeue -u {shlex.quote(self.profile['user'])} -h -o '%i|%j|%T|%M|%P|%Z'")
            if code != 0:
                raise RuntimeError((err or out).strip() or "squeue failed")
            jobs, live = [], {}
            for line in out.splitlines():
                parts = line.strip().split("|")
                if len(parts) >= 6 and parts[0]:
                    job = {"job_id": parts[0], "name": parts[1], "state": parts[2],
                           "time": parts[3], "partition": parts[4],
                           "workdir": parts[5], "recorded": False}
                    live[parts[0]] = job
                    jobs.append(job)
            for entry in recorded or []:
                jid = str(entry.get("job_id"))
                if jid in live:
                    live[jid]["recorded"] = True
                    live[jid]["submitted"] = entry.get("submitted")
                    live[jid]["config"] = entry.get("config")
                else:
                    jobs.append({"job_id": jid, "name": entry.get("job_name"),
                                 "state": None, "time": None,
                                 "partition": entry.get("partition"),
                                 "workdir": entry.get("run_dir"), "recorded": True,
                                 "submitted": entry.get("submitted"),
                                 "config": entry.get("config")})
            finished = [j["job_id"] for j in jobs if j["state"] is None]
            if finished:
                a_out, _, a_code = self._run(
                    client,
                    f"sacct -j {','.join(finished)} -n -P -o JobID,State,ExitCode,Elapsed")
                if a_code == 0:
                    for job in jobs:
                        if job["state"] is None:
                            fin = parse_sacct(a_out, job["job_id"])
                            if fin:
                                job["state"] = fin["state"]
                                job["time"] = fin["elapsed"]
            return jobs
        finally:
            client.close()

    def _poll_loop(self, gen=None):
        if gen is None:
            gen = self._gen
        p = self.profile
        out_file = self._out_file
        job_id = self._job_id
        last_qmsg = None       # last queue-status log line (dedup)
        while not self._stop_flag and gen == self._gen:
            client = None
            try:
                client = self._connect()
                self._tail_remote(client, out_file)
                q_out, _, _ = self._run(
                    client, f"squeue -j {job_id} -h -o '%T|%r|%M|%D|%P'")
                info = parse_squeue(q_out)
                if gen != self._gen:
                    return     # the monitor switched jobs while we polled
                if info:
                    with self._lock:
                        prev = self._remote.get("state")
                        self._remote = info
                        self._state = ("running"
                                       if info["state"] in ("RUNNING", "COMPLETING")
                                       else "queued")
                    if info["state"] != prev:
                        if info["state"] == "RUNNING":
                            self._ingest(
                                f"▶ Job {self._job_id} is running "
                                f"({info['nodes']} node(s), partition {info['partition']}).")
                        else:
                            reason = (f" ({info['reason']})"
                                      if info.get("reason") not in (None, "", "None") else "")
                            self._ingest(f"Job {self._job_id}: {info['state']}{reason}")
                    if info["state"] == "PENDING":
                        pos, total = self._queue_position(
                            client, info.get("partition") or p["partition"])
                        if pos:
                            qmsg = (f"Waiting in queue: position {pos} of {total} "
                                    f"in partition {info.get('partition') or p['partition']}")
                            if qmsg != last_qmsg:
                                self._ingest(qmsg)
                                last_qmsg = qmsg
                else:
                    a_out, _, _ = self._run(
                        client, f"sacct -j {job_id} -n -P -o JobID,State,ExitCode,Elapsed")
                    fin = parse_sacct(a_out, job_id)
                    if gen != self._gen:
                        return
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
                if self._stop_flag or gen != self._gen:
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
