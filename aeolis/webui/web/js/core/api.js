/* Backend REST helpers. JSON for control, ArrayBuffer for bulk data. */
"use strict";

const Api = (() => {

  async function get(path) {
    const res = await fetch(path);
    return _json(res);
  }

  async function post(path, body = {}) {
    const res = await fetch(path, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(body),
    });
    return _json(res);
  }

  async function _json(res) {
    let data = null;
    try { data = await res.json(); } catch { /* non-json error body */ }
    if (!res.ok) {
      const msg = (data && data.error) ? data.error : `${res.status} ${res.statusText}`;
      throw new Error(msg);
    }
    return data;
  }

  /* Binary GET; returns {buffer, headers}. Used for field slabs/meshes. */
  async function binary(path) {
    const res = await fetch(path);
    if (!res.ok) {
      let msg = `${res.status} ${res.statusText}`;
      try { msg = (await res.json()).error || msg; } catch { /* ignore */ }
      throw new Error(msg);
    }
    const buffer = await res.arrayBuffer();
    return { buffer, headers: res.headers };
  }

  /* Poll a background job until done/error. onProgress(job) optional. */
  async function waitJob(jobId, onProgress = null, intervalMs = 400) {
    for (;;) {
      const job = await get(`/api/job/${jobId}`);
      if (onProgress) onProgress(job);
      if (job.status === "done") return job.result;
      if (job.status === "error") throw new Error(job.error || "job failed");
      await new Promise((resolve) => setTimeout(resolve, intervalMs));
    }
  }

  /* Native OS dialogs (served by the backend). Returns path or null. */
  async function pickFile(options = {}) {
    const res = await post("/api/pickfile", options);
    return res.path;
  }
  async function pickFolder(options = {}) {
    const res = await post("/api/pickfolder", options);
    return res.path;
  }

  return { get, post, binary, waitJob, pickFile, pickFolder };
})();
