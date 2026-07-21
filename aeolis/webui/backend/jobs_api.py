"""Background job polling routes."""

from aeolis.webui.backend import jobs
from aeolis.webui.backend.httpd import route
from aeolis.webui.backend.util import send_error_json, send_json


@route("GET", "/api/job/", prefix=True)
def _job_status(handler, query, job_id):
    job = jobs.get(job_id)
    if job is None:
        send_error_json(handler, f"unknown job {job_id}", status=404)
        return
    send_json(handler, job.to_dict())


@route("POST", "/api/job/cancel/", prefix=True)
def _job_cancel(handler, body, job_id):
    job = jobs.cancel(job_id)
    if job is None:
        send_error_json(handler, f"unknown job {job_id}", status=404)
        return
    send_json(handler, job.to_dict())
