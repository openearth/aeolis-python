"""Background job registry for the AeoLiS web GUI.

Long-running work (downloads, interpolation, cache builds) runs on
daemon threads. The frontend polls ``/api/job/<id>`` until the job
reports ``done`` or ``error``.
"""

import threading
import traceback
import uuid

_LOCK = threading.Lock()
_jobs = {}


class Job:

    def __init__(self, label):
        self.id = uuid.uuid4().hex[:12]
        self.label = label
        self.status = "running"      # running | done | error
        self.progress = 0.0          # 0..1, -1 = indeterminate
        self.message = ""
        self.result = None
        self.error = None
        self.cancel_requested = False

    def update(self, progress=None, message=None):
        if progress is not None:
            self.progress = float(progress)
        if message is not None:
            self.message = str(message)

    def to_dict(self):
        return {
            "id": self.id,
            "label": self.label,
            "status": self.status,
            "progress": self.progress,
            "message": self.message,
            "result": self.result if self.status == "done" else None,
            "error": self.error,
        }


def start(label, target):
    """Run ``target(job)`` on a daemon thread; return the job id."""
    job = Job(label)
    with _LOCK:
        _jobs[job.id] = job

    def _runner():
        try:
            job.result = target(job)
            job.status = "done"
            job.progress = 1.0
        except Exception as exc:  # noqa: BLE001 - reported to the client
            job.status = "error"
            job.error = f"{exc}"
            job.message = traceback.format_exc(limit=8)

    thread = threading.Thread(target=_runner, daemon=True, name=f"job-{job.id}")
    thread.start()
    return job.id


def get(job_id):
    with _LOCK:
        return _jobs.get(job_id)


def cancel(job_id):
    with _LOCK:
        job = _jobs.get(job_id)
    if job is not None:
        job.cancel_requested = True
    return job
