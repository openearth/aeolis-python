"""Tests for the AeoLiS web GUI backend (aeolis.webui)."""

import json
import threading
import urllib.error
import urllib.request
from pathlib import Path

import numpy as np
import pytest

from aeolis.constants import DEFAULT_CONFIG
from aeolis.webui.backend import grd_io
from aeolis.webui.backend.datasources import synthetic
from aeolis.webui.backend.run_manager import PROGRESS_RE
from aeolis.webui.backend.schema_api import build_schema


# ---------------------------------------------------------------------
# schema (dynamic settings form)
# ---------------------------------------------------------------------

class TestSchema:

    def test_all_default_config_keys_covered(self):
        from aeolis.webui.backend.schema_api import HIDDEN
        schema = build_schema()
        keys = {p["key"] for s in schema["sections"] for p in s["params"]}
        assert keys == set(DEFAULT_CONFIG) - HIDDEN
        assert "alfa" in HIDDEN

    def test_conditional_visibility_rules(self):
        schema = build_schema()
        by_key = {p["key"]: p for s in schema["sections"] for p in s["params"]}
        assert by_key["Ck"]["visible_if"] == {"key": "method_transport", "in": ["kawamura"]}
        assert by_key["veg_file"]["visible_if"]["in"] == ["duran"]
        assert by_key["hveg_file"]["visible_if"]["in"] == ["grass"]
        # nx/ny are managed by the Grid tab and hidden from the form
        assert "nx" not in by_key and "ny" not in by_key
        # flux factors only shown for 'flux' boundaries
        assert by_key["offshore_flux"]["visible_if"] == {"key": "boundary_offshore", "in": ["flux"]}
        assert by_key["tstop"]["time_tool"]
        # options verified against model source
        assert "sauermann" in by_key["method_transport"]["options"]
        assert by_key["method_moist_process"]["options"] == ["infiltration", "surf_moisture"]
        # masks link to the Domain tab, 4D bedcomp does not
        assert by_key["wave_mask"]["link"] == "domain"
        assert by_key["bedcomp_file"]["link"] is None

    def test_sections_and_metadata(self):
        schema = build_schema()
        by_key = {p["key"]: p for s in schema["sections"] for p in s["params"]}
        tstop = by_key["tstop"]
        assert tstop["type"] == "float"
        assert tstop["unit"] == "s"
        assert "End time" in tstop["desc"]
        # enum options extracted from '(circular, flux or constant)'
        assert by_key["boundary_offshore"]["options"] == ["circular", "flux", "constant"]
        # cross-tab links
        assert by_key["xgrid_file"]["link"] == "grid"
        assert by_key["bed_file"]["link"] == "domain"
        assert by_key["wind_file"]["link"] == "conditions"

    def test_json_serializable(self):
        json.dumps(build_schema())

    def test_section_order_and_disabled_rules(self):
        schema = build_schema()
        names = [s["name"] for s in schema["sections"]]
        assert names[:5] == ["Time settings", "Grid files (*.grd)",
                             "Domain files (*.grd)", "Timeseries",
                             "Output settings"]
        assert names[-1] == "Other"
        by_name = {s["name"]: s for s in schema["sections"]}
        assert by_name["Avalanching"]["enabled_if"] == {"key": "process_avalanche", "in": [True]}
        assert "any" in by_name["Moisture and groundwater"]["enabled_if"]

    def test_output_vars_catalog(self):
        from aeolis.webui.backend.schema_api import list_output_vars
        catalog = {v["name"]: v for v in list_output_vars()}
        assert "zb" in catalog and "Ct" in catalog
        assert catalog["zb"]["dims"][:2] == ["ny", "nx"]
        assert "Bed level" in catalog["zb"]["desc"]

    def test_output_vars_picker_replaces_output_types(self):
        schema = build_schema()
        by_key = {p["key"]: p for s in schema["sections"] for p in s["params"]}
        assert "output_types" not in by_key
        assert by_key["output_vars"]["picker"] == "output_vars"


class TestHpcRunner:

    def test_build_job_script_matches_shape(self):
        from aeolis.webui.backend import run_manager as rm
        s = rm.build_job_script({
            "user": "weste_bt", "partition": "4vcpu", "walltime": "20-00:00:00",
            "env_path": "/p/x/00_environments/env", "run_dir": "/p/x/03_simulations/run",
            "job_name": "NPZK77c", "config": "blowout.txt", "ntasks": 2,
        })
        assert s.startswith("#!/bin/bash\n")
        assert "#SBATCH --job-name=NPZK77c" in s
        assert "#SBATCH --output=NPZK77c.o%j" in s
        assert "#SBATCH --partition=4vcpu" in s
        assert "#SBATCH --ntasks=2" in s
        assert "#SBATCH --time=20-00:00:00" in s
        assert "module load miniforge/latest" in s
        assert "conda activate /p/x/00_environments/env" in s
        assert "cd /p/x/03_simulations/run" in s
        assert s.strip().endswith("aeolis run ./blowout.txt")

    def test_job_script_optional_mail(self):
        from aeolis.webui.backend import run_manager as rm
        s = rm.build_job_script({"mail_user": "a@deltares.nl"})
        assert "#SBATCH --mail-user=a@deltares.nl" in s
        assert "#SBATCH --mail-type=BEGIN,END,FAIL" in s
        assert "--mail-user" not in rm.build_job_script({})

    def test_build_job_script_venv(self):
        # env_kind "venv" activates by sourcing bin/activate and emits no
        # conda/module bootstrap; pre_lines land between cd and aeolis run
        from aeolis.webui.backend import run_manager as rm
        s = rm.build_job_script({
            "env_kind": "venv", "modules": [],
            "env_path": "/p/x/00_environments/aeolis_linux/",
            "run_dir": "/p/x/03_simulations/run",
            "pre_lines": ["mkdir -p output"],
        })
        assert "source /p/x/00_environments/aeolis_linux/bin/activate" in s
        assert "conda" not in s
        assert "module load" not in s
        cd_i = s.index("cd /p/x/03_simulations/run")
        mk_i = s.index("mkdir -p output")
        run_i = s.index("aeolis run ./aeolis.txt")
        assert cd_i < mk_i < run_i
        # default profile stays conda and is unchanged by the new keys
        assert "conda activate" not in rm.build_job_script({})  # no env_path set
        assert "module load miniforge/latest" in rm.build_job_script({})

    def test_parse_job_id(self):
        from aeolis.webui.backend import run_manager as rm
        assert rm.parse_job_id("Submitted batch job 215578") == "215578"
        assert rm.parse_job_id("nope") is None

    def test_path_conversion_roundtrip(self):
        from aeolis.webui.backend import run_manager as rm
        assert rm.local_to_linux(r"P:\proj\run") == "/p/proj/run"
        assert rm.local_to_linux("P:/proj/run") == "/p/proj/run"
        assert rm.linux_to_local("/p/proj/run") == r"P:\proj\run"
        # already-linux path is left alone
        assert rm.local_to_linux("/p/proj/run") == "/p/proj/run"

    def test_copy_project_skips_outputs_and_logs(self, tmp_path):
        # regression: _copy_project used Path without a module-level
        # import -> NameError "Path is not defined" on copy-mode submit
        from aeolis.webui.backend import run_manager as rm
        src = tmp_path / "src"
        (src / "gui").mkdir(parents=True)
        (src / "inputs").mkdir()
        (src / "aeolis.txt").write_text("x")
        (src / "inputs" / "z.grd").write_text("0")
        (src / "old.nc").write_bytes(b"x")
        (src / "run.log").write_text("x")
        (src / "gui" / "state.json").write_text("{}")
        dest = tmp_path / "dest"
        runner = rm.HpcRunner()
        runner._copy_project(src, dest)
        assert (dest / "aeolis.txt").is_file()
        assert (dest / "inputs" / "z.grd").is_file()
        assert not (dest / "old.nc").exists()
        assert not (dest / "run.log").exists()
        assert not (dest / "gui").exists()
        assert any("Copy done" in line for line in runner._lines)

    def test_parse_squeue_and_sacct(self):
        from aeolis.webui.backend import run_manager as rm
        assert rm.parse_squeue("PENDING|Priority|0:00|1|1vcpu")["reason"] == "Priority"
        assert rm.parse_squeue("RUNNING|None|1:23|1|4vcpu")["state"] == "RUNNING"
        assert rm.parse_squeue("") is None
        fin = rm.parse_sacct("215578|COMPLETED|0:0|00:03:36\n215578.batch|COMPLETED|0:0|00:03:36", "215578")
        assert fin["state"] == "COMPLETED" and fin["elapsed"] == "00:03:36"
        assert rm.parse_sacct("", "1") is None

    def test_parse_scontrol_stdout(self):
        from aeolis.webui.backend import run_manager as rm
        text = ("JobId=244359 JobName=aeolis\n"
                "   StdOut=/p/proj/run/aeolis.o244359\n"
                "   StdErr=/p/proj/run/aeolis.o244359\n")
        assert rm.parse_scontrol_stdout(text) == "/p/proj/run/aeolis.o244359"
        assert rm.parse_scontrol_stdout("") is None

    def test_record_hpc_jobs(self, tmp_path):
        """Submissions are remembered per project (newest first, capped,
        deduplicated on job id) so a run survives a GUI restart."""
        from aeolis.webui.backend import project as prj
        from aeolis.webui.backend import run_manager as rm
        p = prj.Project(tmp_path / "aeolis.txt")
        for i in range(12):
            rm.record_hpc_job(p, {"job_id": str(100 + i), "run_dir": "/p/x",
                                  "job_name": "aeolis", "submitted": float(i)})
        jobs = rm.recorded_hpc_jobs(p)
        assert len(jobs) == 10
        assert jobs[0]["job_id"] == "111" and jobs[-1]["job_id"] == "102"
        rm.record_hpc_job(p, {"job_id": "105", "run_dir": "/p/x"})
        jobs = rm.recorded_hpc_jobs(p)
        assert jobs[0]["job_id"] == "105"
        assert len(jobs) == 10
        assert [j["job_id"] for j in jobs].count("105") == 1

    def test_hpc_attach_finished_job(self, monkeypatch):
        """Reattaching to a job that already ended replays its state via
        scontrol/sacct and finishes cleanly (fake SSH transport)."""
        import io as _io
        from aeolis.webui.backend import run_manager as rm

        responses = {
            "scontrol show job": "JobId=999 StdOut=/p/x/aeolis.o999\n",
            "squeue -j": "\n",
            "sacct -j": "999|COMPLETED|0:0|01:00:00\n",
        }

        class _Chan:
            def recv_exit_status(self):
                return 0

        class _Stream(_io.BytesIO):
            channel = _Chan()

        class _FakeSSH:
            def exec_command(self, cmd, timeout=None):
                data = ""
                for key, val in responses.items():
                    if cmd.startswith(key):
                        data = val
                        break
                return None, _Stream(data.encode()), _io.BytesIO(b"")

            def open_sftp(self):
                raise IOError("no sftp in this test")

            def close(self):
                pass

        backend = rm.HpcRunner()
        backend.configure(profile={"user": "u", "host": "h", "run_dir": "/p/x",
                                   "env_path": "/p/env"}, password="pw")
        monkeypatch.setattr(backend, "_connect", lambda: _FakeSSH())
        backend.attach("999")
        backend._poller.join(timeout=10)
        assert not backend._poller.is_alive()
        status = backend.status()
        assert status["state"] == "finished"
        assert status["progress"]["percent"] == 100.0
        assert status["hpc"]["job_id"] == "999"
        assert backend._out_file == "/p/x/aeolis.o999"   # resolved via scontrol
        log = "\n".join(backend._lines)
        assert "Reattaching to job 999" in log

    def test_hpc_attach_requires_password_and_job(self, server_project):
        tmp_path, get, post = server_project
        with pytest.raises(urllib.error.HTTPError) as err:
            post("/api/run/hpc/attach", {"job_id": "1"})
        assert err.value.code == 400
        with pytest.raises(urllib.error.HTTPError) as err:
            post("/api/run/hpc/attach", {"password": "x"})
        assert err.value.code == 400

    def test_hpc_get_lists_recorded_jobs(self, server_project):
        tmp_path, get, post = server_project
        from aeolis.webui.backend import project as prj
        from aeolis.webui.backend import run_manager as rm
        rm.record_hpc_job(prj.current(), {"job_id": "244359", "run_dir": "/p/x",
                                          "job_name": "aeolis", "submitted": 1.0})
        cfg = get("/api/run/hpc")
        assert cfg["jobs"] and cfg["jobs"][0]["job_id"] == "244359"


class TestOutputSource:
    """Viewer output-source override: follow the run folder of an HPC
    submission (mounted on P:) instead of the project's own output."""

    def test_override_roundtrip(self, server_project):
        tmp_path, get, post = server_project
        from aeolis.webui.backend.run_manager import local_to_linux

        run_dir = tmp_path / "p_run"
        run_dir.mkdir()
        (run_dir / "aeolis.txt").write_text((tmp_path / "aeolis.txt").read_text())

        # no override: the output resolves inside the project
        src = get("/api/output/source")
        assert src["override"] is False
        assert Path(src["path"]).parent == tmp_path

        # point at the run folder, given in the cluster's /p path form
        src = post("/api/output/source",
                   {"dir": local_to_linux(str(run_dir)), "config": "aeolis.txt"})
        assert src["override"] is True
        assert Path(src["path"]).parent == run_dir
        assert src["exists"] is False           # nothing written yet

        # meta reports the override even while no output file exists
        meta = get("/api/output/meta")
        assert meta["exists"] is False and meta["override"] is True
        assert Path(meta["path"]).parent == run_dir

        # and back to the project's own file
        src = post("/api/output/source", {"dir": None})
        assert src["override"] is False
        assert Path(src["path"]).parent == tmp_path

    def test_override_rejects_missing_folder(self, server_project):
        tmp_path, get, post = server_project
        with pytest.raises(urllib.error.HTTPError) as err:
            post("/api/output/source", {"dir": str(tmp_path / "nope")})
        assert err.value.code == 404

    def test_override_is_per_project(self, server_project):
        """An override set for one project must not leak into another."""
        import aeolis.inout
        tmp_path, get, post = server_project
        run_dir = tmp_path / "p_run2"
        run_dir.mkdir()
        post("/api/output/source", {"dir": str(run_dir)})
        assert get("/api/output/source")["override"] is True
        other = tmp_path / "other"
        other.mkdir()
        aeolis.inout.write_configfile(str(other / "aeolis.txt"), None)
        post("/api/project/open", {"path": str(other / "aeolis.txt")})
        src = get("/api/output/source")
        assert src["override"] is False
        assert Path(src["path"]).parent == other
        post("/api/output/source", {"dir": None})   # leave no module state behind


class TestRawClean:

    def test_flag_sentinels(self):
        from aeolis.webui.backend.conditions_api import _flag_flaws
        col = np.array([1.0, 999.0, 2.0, -999.0, 3.0])
        flags, counts = _flag_flaws(col, {"sentinels": [999, -999]})
        assert flags.tolist() == [False, True, False, True, False]
        assert counts["= 999"] == 1 and counts["= -999"] == 1

    def test_flag_constant_runs(self):
        from aeolis.webui.backend.conditions_api import _flag_flaws
        col = np.array([1.0, 5.0, 5.0, 5.0, 5.0, 2.0, 3.0, 3.0])
        flags, _ = _flag_flaws(col, {"run_min": 4})
        assert flags.tolist() == [False, True, True, True, True, False, False, False]
        # restricted to a specific value: the run of 5s only
        flags, _ = _flag_flaws(col, {"run_min": 2, "run_value": 3})
        assert flags.tolist() == [False, False, False, False, False, False, True, True]

    def test_runs_not_bridged_by_nan(self):
        from aeolis.webui.backend.conditions_api import _flag_flaws
        col = np.array([0.0, 0.0, np.nan, 0.0, 0.0])
        flags, _ = _flag_flaws(col, {"run_min": 3})
        assert not flags.any()   # two 2-runs separated by NaN, no 3-run


class TestDownloadCaching:

    def test_bounds_tag_distinguishes_areas(self):
        from aeolis.webui.backend.datasources.rws_lidar import bounds_tag
        a = bounds_tag((70000, 445000, 72000, 447000))
        b = bounds_tag((80000, 445000, 82000, 447000))
        assert a != b
        assert bounds_tag((70000, 445000, 72000, 447000)) == a
        assert len(a) == 6


# ---------------------------------------------------------------------
# grid generation
# ---------------------------------------------------------------------

class TestGrdIO:

    def test_generate_shape_and_spacing(self):
        X, Y = grd_io.generate(x0=1000.0, y0=2000.0, dx=5.0, nx=12, ny=8, rotation_deg=0.0)
        assert X.shape == (9, 13)
        assert np.allclose(np.diff(X, axis=1), 5.0)
        assert np.allclose(np.diff(Y, axis=0), 5.0)

    def test_rotation_roundtrip(self):
        X, Y = grd_io.generate(75000, 450000, 2.5, 40, 30, 33.0)
        params = grd_io.derive_params(X, Y)
        assert params["uniform"]
        assert params["nx"] == 40 and params["ny"] == 30
        assert params["dx"] == pytest.approx(2.5, rel=1e-9)
        assert params["rotation"] == pytest.approx(33.0, rel=1e-9)
        assert params["x0"] == pytest.approx(75000)
        assert params["y0"] == pytest.approx(450000)

    def test_write_read_roundtrip(self, tmp_path):
        X, Y = grd_io.generate(0, 0, 1.0, 5, 4, 10.0)
        grd_io.write_grd(tmp_path / "x.grd", X)
        X2 = grd_io.read_grd(tmp_path / "x.grd")
        assert np.allclose(X, X2)

    def test_boundary_orientation(self):
        # offshore = column 0 (aeolis.inout.visualize_grid convention)
        X, Y = grd_io.generate(0, 0, 1.0, 10, 6, 0.0)
        geom = grd_io.geometry(X, Y)
        assert geom["boundaries"]["offshore"][0] == pytest.approx(0.0)
        assert geom["boundaries"]["onshore"][0] == pytest.approx(10.0)

    def test_shear_preview_extents(self):
        X, Y = grd_io.generate(0, 0, 1.0, 100, 50, 0.0)
        prev = grd_io.shear_preview(X, Y, 1.0, 1.0, 10.0, 270.0)
        # wind-aligned box must contain the grid + 2x buffer
        assert prev["length"] >= 100 + 20
        assert prev["width"] >= 50 + 20
        assert len(prev["ring"]) == 4

    def test_invalid_input(self):
        with pytest.raises(ValueError):
            grd_io.generate(0, 0, -1.0, 10, 10, 0.0)
        with pytest.raises(ValueError):
            grd_io.generate(0, 0, 1.0, 0, 10, 0.0)


# ---------------------------------------------------------------------
# synthetic boundary conditions
# ---------------------------------------------------------------------

class TestSynthetic:

    def test_time_axis(self):
        t = synthetic.time_axis(0, 86400, 3600)
        assert len(t) == 25 and t[0] == 0 and t[-1] == 86400

    def test_profiles(self):
        t = synthetic.time_axis(0, 4 * 3600, 3600)
        assert np.allclose(synthetic.profile(t, {"type": "constant", "value": 7}), 7)
        blocks = synthetic.profile(t, {"type": "blocks", "values": [1, 2], "block_duration": 7200})
        assert list(blocks) == [1, 1, 2, 2, 1]
        linear = synthetic.profile(t, {"type": "linear", "start": 0, "end": 8})
        assert linear[0] == 0 and linear[-1] == 8
        rot = synthetic.profile(t, {"type": "rotational", "start": 350, "rate": 10})
        assert rot[0] == pytest.approx(350) and rot[2] == pytest.approx(10)

    def test_harmonic_tide(self):
        data = synthetic.tide(0, 12.42 * 3600, 600,
                              {"type": "harmonic", "mean": 0.3, "amplitude": 1.0,
                               "period": 12.42 * 3600})
        assert data.shape[1] == 2
        assert data[:, 1].mean() == pytest.approx(0.3, abs=0.05)

    def test_wind_columns_and_clipping(self):
        data = synthetic.wind(0, 3600, 600,
                              {"type": "harmonic", "mean": 1, "amplitude": 5, "period": 1800},
                              {"type": "constant", "value": 270})
        assert data.shape[1] == 3
        assert data[:, 1].min() >= 0.0          # no negative wind speeds
        assert np.allclose(data[:, 2], 270.0)

    def test_file_format_matches_model(self, tmp_path):
        data = synthetic.wind(0, 3600, 600, {"type": "constant", "value": 10},
                              {"type": "constant", "value": 270})
        synthetic.write_series(tmp_path / "wind.txt", data)
        loaded = np.loadtxt(tmp_path / "wind.txt")
        assert loaded.shape == data.shape


class TestResolveTarget:
    """resolve_target must never leak a '..' segment into the config ref
    (the file browser can hand back .../input/../input/wind.txt)."""

    def _current(self, root):
        import types
        return types.SimpleNamespace(root=root)

    def test_collapses_dotdot_to_clean_relative(self, tmp_path):
        from aeolis.webui.backend.grid_api import resolve_target
        root = tmp_path / "input"
        root.mkdir()
        weird = str(root / ".." / "input" / "wind.txt")
        write_path, ref = resolve_target(self._current(root), weird, "wind.txt")
        assert ".." not in ref
        assert ref == "wind.txt"
        assert Path(write_path).name == "wind.txt"

    def test_outside_root_stays_absolute(self, tmp_path):
        from aeolis.webui.backend.grid_api import resolve_target
        root = tmp_path / "proj"
        root.mkdir()
        other = tmp_path / "elsewhere" / "wind.txt"
        write_path, ref = resolve_target(self._current(root), str(other), "wind.txt")
        assert Path(ref).is_absolute()

    def test_bare_name_joins_root(self, tmp_path):
        from aeolis.webui.backend.grid_api import resolve_target
        root = tmp_path / "proj"
        root.mkdir()
        write_path, ref = resolve_target(self._current(root), "wind.txt", "wind.txt")
        assert ref == "wind.txt"


class TestFillColumns:
    """The priority-layered column builder and its remaining-NaN policies."""

    def test_priority_and_fill_methods(self):
        from aeolis.webui.backend.conditions_api import _column_from_sources
        master = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
        # top source covers only 0..2 -> 3,4 stay NaN until the fill method
        top = [(np.array([0.0, 1.0, 2.0]), np.array([10.0, 11.0, 12.0]))]

        val = _column_from_sources(master, top, {"method": "value", "value": -1}, None)
        assert val[0] == 10 and val[2] == 12 and val[3] == -1 and val[4] == -1

        near = _column_from_sources(master, top, {"method": "nearest"}, None)
        assert near[3] == 12 and near[4] == 12

        lin = _column_from_sources(master, top, {"method": "linear"}, None)
        assert lin[3] == pytest.approx(12) and lin[4] == pytest.approx(12)

    def test_lower_priority_fills_gaps(self):
        from aeolis.webui.backend.conditions_api import _column_from_sources
        master = np.array([0.0, 1.0, 2.0, 3.0])
        top = (np.array([0.0, 1.0]), np.array([5.0, 6.0]))          # covers 0..1
        low = (np.array([2.0, 3.0]), np.array([70.0, 80.0]))        # covers 2..3
        out = _column_from_sources(master, [top, low], {"method": "value", "value": 0}, None)
        assert list(out) == [5.0, 6.0, 70.0, 80.0]

    def test_internal_gap_not_bridged_by_top_source(self):
        # a stretch of MISSING TIMESTAMPS in the top source must fall through
        # to the lower-priority source, not be linearly bridged by the top one
        from aeolis.webui.backend.conditions_api import _column_from_sources
        master = np.array([0.0, 1.0, 2.0, 5.0, 6.0, 10.0, 11.0, 12.0])
        top = (np.array([0.0, 1.0, 2.0, 10.0, 11.0, 12.0]),   # gap 2 -> 10
               np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0]))
        low = (np.array([4.0, 5.0, 6.0, 7.0]),
               np.array([9.0, 9.0, 9.0, 9.0]))
        out = _column_from_sources(master, [top, low], {"method": "value", "value": 0}, None)
        assert list(out) == [1.0, 1.0, 1.0, 9.0, 9.0, 1.0, 1.0, 1.0]

    def test_nan_gap_falls_through_to_lower_source(self):
        # ... and so must a stretch of NaN samples
        from aeolis.webui.backend.conditions_api import _column_from_sources
        master = np.arange(0.0, 7.0)
        top = (master.copy(),
               np.array([1.0, 1.0, np.nan, np.nan, np.nan, 1.0, 1.0]))
        low = (master.copy(), np.full(7, 9.0))
        out = _column_from_sources(master, [top, low], {"method": "value", "value": 0}, None)
        assert list(out) == [1.0, 1.0, 9.0, 9.0, 9.0, 1.0, 1.0]


# ---------------------------------------------------------------------
# run progress parsing
# ---------------------------------------------------------------------

class TestRunProgress:

    def test_progress_line(self):
        line = "010.5%       0:01:23 /    0:10:00 /       0:08:37 / 60.0"
        match = PROGRESS_RE.search(line)
        assert match
        assert float(match.group(1)) == 10.5
        assert match.group(4) == "0:08:37"
        assert float(match.group(5)) == 60.0

    def test_real_model_format(self):
        # exact format string used by AeoLiSRunner.print_progress
        line = "%05.1f%%  %12s / %10s / %14s / %0.1f" % (
            42.0, "0:05:00", "0:12:00", "0:07:00", 3456.0)
        match = PROGRESS_RE.search(line)
        assert match and float(match.group(1)) == 42.0

    def test_times_beyond_a_day(self):
        # timedelta reprs read "1 day, 0:54:00" past 24 h - early in a
        # long run the remaining estimate almost always does, which used
        # to keep the progress bar stuck at 0.0%
        from datetime import timedelta
        line = "%05.1f%%  %12s / %10s / %14s / %0.1f" % (
            0.4, timedelta(seconds=360), timedelta(seconds=90000),
            timedelta(seconds=89640), 36.0)
        assert "day" in line
        match = PROGRESS_RE.search(line)
        assert match
        assert float(match.group(1)) == 0.4
        assert match.group(2) == "0:06:00"
        assert match.group(3) == "1 day, 1:00:00"
        assert match.group(4) == "1 day, 0:54:00"
        assert float(match.group(5)) == 36.0


# ---------------------------------------------------------------------
# HTTP API integration (config + grid roundtrip on a temp project)
# ---------------------------------------------------------------------

@pytest.fixture()
def server_project(tmp_path, monkeypatch):
    import aeolis.inout
    from aeolis.webui.backend import settings
    from aeolis.webui.backend.httpd import make_server

    # keep test projects out of the user's real recent-projects list and
    # style edits out of the user's real colormap presets
    monkeypatch.setattr(settings, "RECENT_FILE", tmp_path / "recent.json")
    monkeypatch.setattr(settings, "STYLES_FILE", tmp_path / "styles.json")

    configfile = tmp_path / "aeolis.txt"
    aeolis.inout.write_configfile(str(configfile), None)

    server = make_server("127.0.0.1", base_port=9450, tries=30)
    host, port = server.server_address[:2]
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    base = f"http://{host}:{port}"

    def get(path):
        with urllib.request.urlopen(base + path) as response:
            return json.loads(response.read())

    def post(path, body):
        request = urllib.request.Request(
            base + path, data=json.dumps(body).encode(),
            headers={"Content-Type": "application/json"})
        with urllib.request.urlopen(request) as response:
            return json.loads(response.read())

    post("/api/project/open", {"path": str(configfile)})
    yield tmp_path, get, post
    server.shutdown()


def _wait_job(get, job_id, timeout=15.0):
    """Poll a background job until it finishes (or the timeout elapses)."""
    import time
    deadline = time.time() + timeout
    while time.time() < deadline:
        job = get(f"/api/job/{job_id}")
        if job["status"] in ("done", "error"):
            return job
        time.sleep(0.05)
    raise AssertionError(f"job {job_id} did not finish in {timeout}s")


class TestApi:

    def test_ping_and_project(self, server_project):
        tmp_path, get, post = server_project
        assert get("/api/ping")["ok"]
        assert (tmp_path / "gui" / "rawdata").is_dir()

    def test_config_roundtrip(self, server_project):
        tmp_path, get, post = server_project
        cfg = get("/api/config")
        values = cfg["values"]
        values["tstop"] = 999.0
        values["process_avalanche"] = True
        post("/api/config/save", {"values": values})

        import aeolis.inout
        parsed = aeolis.inout.read_configfile(str(tmp_path / "aeolis.txt"), parse_files=False)
        assert parsed["tstop"] == 999.0
        assert parsed["process_avalanche"] is True or parsed["process_avalanche"] == True  # noqa: E712

    def test_grid_save_and_model_parse(self, server_project):
        tmp_path, get, post = server_project
        result = post("/api/grid/save", {
            "x0": 100.0, "y0": 200.0, "dx": 2.0, "nx": 10, "ny": 6, "rotation": 20.0,
        })
        assert result["ok"]
        X = grd_io.read_grd(tmp_path / "x.grd")
        assert X.shape == (7, 11)

        grid = get("/api/grid")
        assert grid["exists"]
        assert grid["params"]["nx"] == 10
        assert grid["params"]["rotation"] == pytest.approx(20.0)

        import aeolis.inout
        parsed = aeolis.inout.read_configfile(str(tmp_path / "aeolis.txt"), parse_files=False)
        assert parsed["nx"] == 10 and parsed["ny"] == 6 and parsed["alfa"] == 0

    def test_synthetic_conditions(self, server_project):
        tmp_path, get, post = server_project
        result = post("/api/conditions/synthetic", {
            "kind": "wind", "dt": 3600, "tstart": 0, "tstop": 86400,
            "speed": {"type": "constant", "value": 12},
            "direction": {"type": "constant", "value": 300},
        })
        assert result["ok"] and result["rows"] == 25
        data = np.loadtxt(tmp_path / "wind.txt")
        assert data.shape == (25, 3)
        assert np.allclose(data[:, 1], 12.0)

    def test_domain_modify_with_init(self, server_project):
        tmp_path, get, post = server_project
        post("/api/grid/save", {
            "x0": 0.0, "y0": 0.0, "dx": 1.0, "nx": 8, "ny": 8, "rotation": 0.0,
        })
        result = post("/api/domain/modify", {
            "target": "bed", "op": "set", "value": 3.0,
            "indices": [0, 3, 0, 3], "init": 0.0,
        })
        assert result["cells"] == 16
        Z = np.loadtxt(tmp_path / "zb.grd")
        assert Z.shape == (9, 9)
        assert np.allclose(Z[0:4, 0:4], 3.0)
        assert np.allclose(Z[5:, 5:], 0.0)

    def test_run_checklist_catches_shape_mismatch(self, server_project):
        tmp_path, get, post = server_project
        post("/api/grid/save", {
            "x0": 0.0, "y0": 0.0, "dx": 1.0, "nx": 10, "ny": 10, "rotation": 0.0,
        })
        # wrong-shaped bed file
        np.savetxt(tmp_path / "zb_bad.grd", np.zeros((3, 4)))
        cfg = get("/api/config")
        values = cfg["values"]
        values["bed_file"] = "zb_bad.grd"
        post("/api/config/save", {"values": values})

        res = get("/api/run/checklist")
        assert not res["ready"]
        errors = [c["text"] for c in res["checks"] if c["level"] == "error"]
        assert any("does not match" in t for t in errors)

    def test_synthetic_raw_single_variable(self, server_project):
        tmp_path, get, post = server_project
        res = post("/api/conditions/synthetic_raw", {
            "variable": "wind_dir", "dt": 3600, "tstart": 0, "tstop": 86400,
            "segments": {"type": "segments",
                         "segments": [{"type": "constant", "value": 123}]},
        })
        assert res["ok"]
        entry = res["entry"]
        # a single-variable series carries its own column label, not the kind default
        assert entry["labels"] == ["direction [deg]"]
        raw = get(f"/api/conditions/raw_series?id={entry['id']}")
        cols = raw["series"]["columns"]
        assert len(cols) == 1
        assert all(abs(v - 123) < 1e-6 for v in cols[0])

    def test_fill_columns_body(self, server_project):
        tmp_path, get, post = server_project
        r = post("/api/conditions/synthetic_raw", {
            "kind": "wind", "dt": 3600, "tstart": 0, "tstop": 86400,
            "speed": {"type": "constant", "value": 8},
            "direction": {"type": "constant", "value": 210},
        })
        rid = r["entry"]["id"]
        res = post("/api/conditions/fill", {
            "kind": "wind",
            "columns": [
                {"sources": [{"id": rid, "column": 0}], "fill": {"method": "linear"}},
                {"sources": [{"id": rid, "column": 1}], "fill": {"method": "nearest"}},
            ],
        })
        assert res["ok"]
        data = np.loadtxt(tmp_path / "wind.txt")
        assert data.shape[1] == 3
        assert np.allclose(data[:, 1], 8) and np.allclose(data[:, 2], 210)

    def test_fill_secondary_station_fills_primary_gap(self, server_project):
        # primary station with a 6h hole (no samples at all) + a secondary
        # station: the hole must appear in the output, valued from the
        # secondary — the user's Stortemelk-west + Amelander-zeegat case
        tmp_path, get, post = server_project
        mk = lambda v: post("/api/conditions/synthetic_raw", {
            "kind": "tide", "dt": 3600, "tstart": 0, "tstop": 86400,
            "level": {"type": "constant", "value": v},
        })["entry"]
        prim, sec = mk(1.0), mk(9.0)
        # punch the hole into the primary's stored npz (samples 6h..12h gone)
        path = tmp_path / prim["path"]
        data = np.load(path)
        t, cols = data["t"], np.atleast_2d(data["cols"])
        if cols.shape[0] != t.size:
            cols = cols.T
        rel = t - t[0]                       # npz stores epoch seconds
        keep = (rel < 6 * 3600) | (rel > 12 * 3600)
        np.savez_compressed(path, t=t[keep], cols=cols[keep])

        res = post("/api/conditions/fill", {
            "kind": "tide",
            "columns": [{"sources": [{"id": prim["id"], "column": 0},
                                     {"id": sec["id"], "column": 0}],
                         "fill": {"method": "linear"}}],
        })
        assert res["ok"]
        out = np.atleast_2d(np.loadtxt(tmp_path / "tide.txt"))
        t_out, v_out = out[:, 0], out[:, 1]
        hole = (t_out >= 7 * 3600) & (t_out <= 11 * 3600)
        assert hole.any(), "gap timestamps missing from the output"
        assert np.allclose(v_out[hole], 9.0), "gap not filled by the secondary station"
        assert np.allclose(v_out[t_out < 6 * 3600], 1.0)

    def test_duplicate_include_outputs(self, server_project):
        tmp_path, get, post = server_project
        (tmp_path / "aeolis.nc").write_text("nc")
        (tmp_path / "wind.txt").write_text("0 1 2\n")
        # the copy destination must live OUTSIDE the project root
        parent = tmp_path.parent / f"dup_{tmp_path.name}"
        parent.mkdir()

        post("/api/project/duplicate",
             {"parent": str(parent), "name": "noout", "include_outputs": False})
        assert (parent / "noout" / "wind.txt").is_file()
        assert not (parent / "noout" / "aeolis.nc").is_file()

        post("/api/project/open", {"path": str(tmp_path / "aeolis.txt")})
        post("/api/project/duplicate",
             {"parent": str(parent), "name": "withnc", "include_outputs": True})
        assert (parent / "withnc" / "aeolis.nc").is_file()

    def test_duplicate_gather_external_input(self, server_project):
        tmp_path, get, post = server_project
        # an input file referenced from OUTSIDE the project root
        ext = tmp_path.parent / f"ext_{tmp_path.name}"
        ext.mkdir()
        grd_io.write_grd(ext / "zext.grd", np.ones((4, 4)))
        values = get("/api/config")["values"]
        values["bed_file"] = str(ext / "zext.grd")
        post("/api/config/save", {"values": values})

        parent = tmp_path.parent / f"dupg_{tmp_path.name}"
        parent.mkdir()
        post("/api/project/duplicate",
             {"parent": str(parent), "name": "g", "input_mode": "gather"})
        # gather copies the external file in and repoints the config locally
        assert (parent / "g" / "zext.grd").is_file()
        import aeolis.inout
        parsed = aeolis.inout.read_configfile(
            str(parent / "g" / "aeolis.txt"), parse_files=False)
        assert parsed["bed_file"] == "zext.grd"

    def test_duplicate_keep_external_input(self, server_project):
        tmp_path, get, post = server_project
        ext = tmp_path.parent / f"ext2_{tmp_path.name}"
        ext.mkdir()
        grd_io.write_grd(ext / "zext.grd", np.ones((4, 4)))
        extfile = str(ext / "zext.grd")
        values = get("/api/config")["values"]
        values["bed_file"] = extfile
        post("/api/config/save", {"values": values})

        parent = tmp_path.parent / f"dupk_{tmp_path.name}"
        parent.mkdir()
        post("/api/project/duplicate",
             {"parent": str(parent), "name": "k", "input_mode": "keep"})
        # keep leaves the original in place and stores an absolute link
        assert not (parent / "k" / "zext.grd").is_file()
        import aeolis.inout
        parsed = aeolis.inout.read_configfile(
            str(parent / "k" / "aeolis.txt"), parse_files=False)
        assert Path(parsed["bed_file"]) == Path(extfile)

    def test_domain_sample_rename_file(self, server_project):
        tmp_path, get, post = server_project
        xyz = tmp_path / "pts.xyz"
        xyz.write_text("0 0 1\n1 0 2\n0 1 3\n1 1 4\n")
        eid = post("/api/domain/import_xyz", {"path": str(xyz)})["entry"]["id"]
        old_path = get("/api/domain")["entries"][0]["path"]

        res = post("/api/domain/sample_rename",
                   {"id": eid, "name": "My Bed Points", "rename_file": True})
        entry = res["entry"]
        assert entry["label"] == "My Bed Points"
        # the file on disk is renamed to a slug and the id follows the new path
        assert entry["path"].endswith("My_Bed_Points.npz")
        assert (tmp_path / entry["path"]).is_file()
        assert not (tmp_path / old_path).is_file()
        assert entry["id"] != eid

    def test_link_external_references_in_place(self, server_project):
        import json
        tmp_path, get, post = server_project
        # build a second project on disk with one downloaded dataset
        ext = tmp_path.parent / "otherproj"
        ext_raw = ext / "gui" / "rawdata"
        ext_raw.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(ext_raw / "shared.npz",
                            x=np.array([0.0, 1.0]), y=np.array([0.0, 1.0]),
                            z=np.array([2.0, 3.0], dtype="float32"))
        (ext_raw / "manifest.json").write_text(json.dumps({"entries": [{
            "id": "abc123", "source": "xyz", "kind": "points",
            "path": "gui/rawdata/shared.npz", "label": "Shared points",
        }]}))

        scan = post("/api/domain/scan_external", {"path": str(ext)})
        assert len(scan["entries"]) == 1
        assert scan["entries"][0]["exists"] and not scan["entries"][0]["is_local"]

        res = post("/api/domain/link_external", {"path": str(ext), "ids": ["abc123"]})
        assert res["linked"] == 1
        entry = get("/api/domain")["entries"][0]
        assert entry["linked"] is True
        # the manifest points at the external file, not a local copy
        assert entry["path"] == str((ext_raw / "shared.npz").resolve())
        assert not (tmp_path / "gui" / "rawdata" / "shared.npz").exists()

        # deleting a linked entry must NOT remove the external source file
        post("/api/domain/forget", {"id": entry["id"], "delete_file": True})
        assert (ext_raw / "shared.npz").is_file()
        assert not get("/api/domain")["entries"]

    def test_link_external_bare_data_folder(self, server_project):
        # a shared data folder (e.g. 01_data/lidar on a network drive)
        # carries its manifest NEXT TO the files with plain filenames as
        # paths - no gui/rawdata nesting
        import json
        tmp_path, get, post = server_project
        ext = tmp_path.parent / "shared_data" / "lidar"
        ext.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(ext / "flat.npz",
                            x=np.array([0.0, 1.0]), y=np.array([0.0, 1.0]),
                            z=np.array([2.0, 3.0], dtype="float32"))
        (ext / "manifest.json").write_text(json.dumps({"entries": [{
            "id": "bare01", "source": "rws_lidar", "kind": "points",
            "path": "flat.npz", "label": "Bare-folder points",
        }]}))

        scan = post("/api/domain/scan_external", {"path": str(ext)})
        assert len(scan["entries"]) == 1
        e = scan["entries"][0]
        assert e["exists"] and not e["is_local"]
        assert e["abspath"] == str((ext / "flat.npz").resolve())

        res = post("/api/domain/link_external", {"path": str(ext), "ids": ["bare01"]})
        assert res["linked"] == 1
        entry = get("/api/domain")["entries"][0]
        assert entry["linked"] is True
        assert entry["path"] == str((ext / "flat.npz").resolve())

    def test_target_constant_and_modify(self, server_project):
        tmp_path, get, post = server_project
        post("/api/grid/save",
             {"x0": 0.0, "y0": 0.0, "dx": 1.0, "nx": 4, "ny": 4, "rotation": 0.0})
        r = post("/api/domain/target_constant", {"target": "bed", "value": 1})
        assert r["min"] == 1 and r["max"] == 1
        assert get("/api/domain")["targets"]["bed"]["has_draft"] is True
        # set a 2x2 index box to 0 (wave-mask style edit)
        r = post("/api/domain/target_modify", {"target": "bed", "op": "set", "value": 0,
                 "scope": {"type": "indices", "indices": [0, 1, 0, 1]}})
        assert r["cells"] == 4 and r["min"] == 0 and r["max"] == 1
        post("/api/domain/target_save", {"target": "bed"})
        assert (tmp_path / "zb.grd").is_file()

    def test_raster_derive_ndvi_threshold(self, server_project):
        import rasterio
        from rasterio.transform import from_origin
        tmp_path, get, post = server_project
        tif = tmp_path / "cir.tif"
        nir = np.array([[0.8, 0.1], [0.9, 0.2]], dtype="float32")
        red = np.array([[0.2, 0.3], [0.1, 0.4]], dtype="float32")
        with rasterio.open(tif, "w", driver="GTiff", height=2, width=2, count=2,
                           dtype="float32", crs="EPSG:28992",
                           transform=from_origin(0, 2, 1, 1)) as ds:
            ds.write(nir, 1)
            ds.write(red, 2)
        eid = post("/api/domain/import_tiff", {"path": str(tif)})["entry"]["id"]
        ent = next(e for e in get("/api/domain")["entries"] if e["id"] == eid)
        assert ent["bands"] == 2
        # NDVI = (b1-b2)/(b1+b2), classify > 0.4 -> 1 else 0
        r = post("/api/domain/raster_derive", {"id": eid, "expr": "(b1-b2)/(b1+b2)",
                 "threshold": {"op": ">", "x": 0.4, "then": 1, "else": 0}, "save_as": "veg"})
        assert r["min"] == 0 and r["max"] == 1
        assert "veg" in [e["label"] for e in get("/api/domain")["entries"]]

    def test_domain_target_draft_load_and_save(self, server_project):
        tmp_path, get, post = server_project
        post("/api/grid/save",
             {"x0": 0.0, "y0": 0.0, "dx": 1.0, "nx": 4, "ny": 4, "rotation": 0.0})
        grd_io.write_grd(tmp_path / "src.grd", np.full((5, 5), 7.0))

        r = post("/api/domain/target_load",
                 {"target": "bed", "path": str(tmp_path / "src.grd")})
        assert r["ok"] and r["draft"] and r["shape_ok"]
        # loaded as a draft: nothing written to disk yet
        assert get("/api/domain")["targets"]["bed"]["has_draft"] is True
        assert not (tmp_path / "zb.grd").is_file()

        s = post("/api/domain/target_save", {"target": "bed"})
        Z = grd_io.read_grd(tmp_path / s["file"])
        assert np.allclose(Z, 7.0)
        # after save the draft is cleared and the file exists
        ov = get("/api/domain")["targets"]["bed"]
        assert ov["has_draft"] is False and ov["exists"] is True

    def test_domain_interpolate_produces_draft(self, server_project):
        tmp_path, get, post = server_project
        post("/api/grid/save",
             {"x0": 0.0, "y0": 0.0, "dx": 1.0, "nx": 4, "ny": 4, "rotation": 0.0})
        xyz = tmp_path / "pts.xyz"
        lines = [f"{i} {j} {i + j}" for i in range(5) for j in range(5)]
        xyz.write_text("\n".join(lines) + "\n")
        eid = post("/api/domain/import_xyz", {"path": str(xyz)})["entry"]["id"]

        res = post("/api/domain/interpolate",
                   {"target": "bed", "layers": [eid], "extrapolate": True})
        job = _wait_job(get, res["job"])
        assert job["status"] == "done", job.get("error")
        assert job["result"]["draft"] is True
        # interpolation makes a draft, it does NOT write the .grd
        assert get("/api/domain")["targets"]["bed"]["has_draft"] is True
        assert not (tmp_path / "zb.grd").is_file()
        # committing the draft writes the file
        post("/api/domain/target_save", {"target": "bed"})
        assert (tmp_path / "zb.grd").is_file()

    def test_interpolate_copies_from_another_target(self, server_project):
        tmp_path, get, post = server_project
        post("/api/grid/save",
             {"x0": 0.0, "y0": 0.0, "dx": 1.0, "nx": 4, "ny": 4, "rotation": 0.0})
        # build a wave mask draft: 1 everywhere, a 2x2 box set to 0
        post("/api/domain/target_constant", {"target": "wave_mask", "value": 1})
        post("/api/domain/target_modify", {"target": "wave_mask", "op": "set", "value": 0,
             "scope": {"type": "indices", "indices": [0, 1, 0, 1]}})
        # tide_mask = direct copy of the wave mask via a "tgt:" source
        res = post("/api/domain/interpolate",
                   {"target": "tide_mask", "layers": ["tgt:wave_mask"]})
        job = _wait_job(get, res["job"])
        assert job["status"] == "done", job.get("error")
        assert job["result"]["min"] == 0 and job["result"]["max"] == 1
        post("/api/domain/target_save", {"target": "tide_mask"})
        Z = grd_io.read_grd(tmp_path / "tide_mask.grd")
        assert Z.shape == (5, 5)                      # nx cells -> nx+1 nodes
        assert Z[0, 0] == 0 and Z[4, 4] == 1 and Z.sum() == 21

    def test_save_as_rejects_glued_absolute_path(self, server_project):
        # a picker bug once produced "<root>\C:\other\file.grd"; the backend
        # must refuse such a path with a clear error instead of WinError 123
        tmp_path, get, post = server_project
        post("/api/grid/save",
             {"x0": 0.0, "y0": 0.0, "dx": 1.0, "nx": 4, "ny": 4, "rotation": 0.0})
        post("/api/domain/target_constant", {"target": "bed", "value": 1})
        glued = str(tmp_path) + "\\C:\\somewhere\\else\\zb.grd"
        with pytest.raises(urllib.error.HTTPError) as err:
            post("/api/domain/target_save_as", {"target": "bed", "path": glued})
        assert err.value.code in (400, 500)

    def test_era5_cached_summary_reuses_covering_parts(self, tmp_path):
        from datetime import datetime, timezone
        from aeolis.webui.backend.datasources import era5
        # AL00_a-style cache: 2011 part starts mid-year, 2012 is a full year
        (tmp_path / "era5_wind_5.50_53.50_20110601_20111231.nc").touch()
        (tmp_path / "era5_wind_5.50_53.50_20120101_20121231.nc").touch()
        d = lambda y, m, day: datetime(y, m, day, tzinfo=timezone.utc)  # noqa: E731
        # request starting inside the partial 2011 part: covered
        s = era5.cached_summary(5.5, 53.5, d(2011, 7, 1), d(2013, 12, 31), tmp_path)
        assert s == {"cached": [2011, 2012], "missing": [2013], "total": 3}
        # request starting BEFORE the partial part: 2011 must re-download
        s = era5.cached_summary(5.5, 53.5, d(2011, 1, 1), d(2012, 12, 31), tmp_path)
        assert s["cached"] == [2012] and s["missing"] == [2011]
        # a different cell shares nothing
        s = era5.cached_summary(5.75, 53.5, d(2012, 1, 1), d(2012, 12, 31), tmp_path)
        assert s["cached"] == [] and s["missing"] == [2012]

    def test_era5_cached_endpoint(self, server_project):
        tmp_path, get, post = server_project
        raw = tmp_path / "gui" / "rawdata"
        raw.mkdir(parents=True, exist_ok=True)
        (raw / "era5_wind_5.50_53.50_20200101_20201231.nc").touch()
        res = post("/api/conditions/era5_cached",
                   {"lon": 5.5, "lat": 53.5, "date0": "2020-01-01", "date1": "2021-12-31"})
        assert res == {"cached": [2020], "missing": [2021], "total": 2}

    def test_save_as_writes_portable_relative_refs(self, server_project):
        # default path format: relative with FORWARD slashes (runs on /p too)
        tmp_path, get, post = server_project
        post("/api/grid/save",
             {"x0": 0.0, "y0": 0.0, "dx": 1.0, "nx": 4, "ny": 4, "rotation": 0.0})
        post("/api/domain/target_constant", {"target": "bed", "value": 1})
        sub = tmp_path / "bathy"
        sub.mkdir()
        res = post("/api/domain/target_save_as",
                   {"target": "bed", "path": str(sub / "zb.grd")})
        assert res["file"] == "bathy/zb.grd"
        assert (sub / "zb.grd").is_file()

    def test_save_as_honours_absolute_path_format(self, server_project):
        tmp_path, get, post = server_project
        post("/api/grid/save",
             {"x0": 0.0, "y0": 0.0, "dx": 1.0, "nx": 4, "ny": 4, "rotation": 0.0})
        post("/api/domain/target_constant", {"target": "bed", "value": 1})
        state = get("/api/project/state") or {}
        state.setdefault("ui", {})["pathFormat"] = "absolute"
        post("/api/project/state", state)
        res = post("/api/domain/target_save_as",
                   {"target": "bed", "path": str(tmp_path / "zb2.grd")})
        assert res["file"] == str(tmp_path / "zb2.grd")

    def test_linuxify_config(self, tmp_path):
        from aeolis.webui.backend.run_manager import linuxify_config
        cfg = tmp_path / "aeolis.txt"
        cfg.write_text(
            "bed_file    = C:\\proj\\input\\z.grd % bed level\n"
            "wind_file   = sub\\wind.txt\n"
            "tide_file   = tide.txt   % local, untouched\n"
            "nx          = 100\n")
        n = linuxify_config(cfg)
        text = cfg.read_text()
        assert n == 2
        assert "/c/proj/input/z.grd" in text
        assert "sub/wind.txt" in text
        assert "% bed level" in text
        assert "tide_file   = tide.txt   % local, untouched" in text
        assert "nx          = 100" in text

    def test_linuxify_config_forces_ascii(self, tmp_path):
        """Regression: a cp1252 dash in a comment (from the "Lotka-Volterra"
        description in constants.py) made the cluster abort with
        UnicodeDecodeError, because Linux opens the config as UTF-8."""
        from aeolis.webui.backend.run_manager import linuxify_config
        import aeolis.inout
        cfg = tmp_path / "aeolis.txt"
        cfg.write_bytes(
            "alpha_comp  = 0    % [-] Lotka–Volterra competition\n"
            "kappa       = 0.41 % [-] Von Kármán constant\n"
            "nx          = 100\n".encode("cp1252"))
        n = linuxify_config(cfg)
        raw = cfg.read_bytes()
        assert n == 2                                   # two comments folded
        assert all(b < 128 for b in raw), raw
        raw.decode("utf-8")                             # what the cluster does
        text = cfg.read_text(encoding="utf-8")
        assert "Lotka-Volterra" in text and "Von Karman" in text
        assert aeolis.inout.read_configfile(
            str(cfg), parse_files=False)["nx"] == 100

    def test_write_configfile_is_ascii_utf8(self, tmp_path):
        """The config writer must never emit platform-encoded bytes: the
        file is written on Windows and read back on the Linux cluster."""
        import aeolis.inout
        from aeolis.constants import DEFAULT_CONFIG
        cfg = tmp_path / "aeolis.txt"
        values = DEFAULT_CONFIG.copy()
        values["alpha_comp"] = [0.5]        # non-default -> comment is written
        values["kappa"] = 0.42
        aeolis.inout.write_configfile(str(cfg), values)
        raw = cfg.read_bytes()
        assert all(b < 128 for b in raw), [b for b in raw if b > 127]
        text = raw.decode("utf-8")
        assert "Lotka-Volterra" in text and "Von Karman" in text

    def test_write_configfile_keeps_float_parameters_float(self, tmp_path):
        """Regression: whole numbers typed in the GUI were written as
        "G_h = 4" and read back as int, so grass.initialize's in-place
        `p['G_h'] /= 365.25*24*3600` raised UFuncTypeError on the cluster."""
        import aeolis.inout
        from aeolis.constants import DEFAULT_CONFIG
        cfg = tmp_path / "aeolis.txt"
        values = DEFAULT_CONFIG.copy()
        values.update({"G_h": 4, "G_c": 2, "dzb_opt_c": [1, 2],
                       "nx": 100, "nlayers": 40})
        aeolis.inout.write_configfile(str(cfg), values)
        p = aeolis.inout.read_configfile(str(cfg), parse_files=False)
        for key in ("G_h", "G_c"):
            assert isinstance(p[key], float), (key, p[key], type(p[key]))
        assert np.asarray(p["dzb_opt_c"]).dtype.kind == "f", p["dzb_opt_c"]
        # integer parameters stay integers
        assert isinstance(p["nx"], int) and isinstance(p["nlayers"], int)
        # and the values themselves are unchanged
        assert p["G_h"] == 4.0 and p["G_c"] == 2.0

    def test_read_configfile_falls_back_to_latin1(self, tmp_path):
        """Configs written by older versions on Windows must still load."""
        import aeolis.inout
        cfg = tmp_path / "aeolis.txt"
        cfg.write_bytes(
            "nx = 100 % [-] Lotka–Volterra\n".encode("cp1252"))
        assert aeolis.inout.read_configfile(
            str(cfg), parse_files=False)["nx"] == 100

    def test_unknown_route_404(self, server_project):
        tmp_path, get, post = server_project
        with pytest.raises(urllib.error.HTTPError) as err:
            get("/api/nonexistent")
        assert err.value.code == 404

    def test_resample_starts_at_input_start(self):
        """Resampled series must keep the exact start time of the input;
        stamping bins at their centres shifted a projected wind file by
        half an interval (first sample at t = 1800 instead of 0)."""
        from aeolis.webui.backend.conditions_api import _resample
        t = 1.3e9 + np.arange(0.0, 7200.0 + 1.0, 600.0)
        v = np.arange(t.size, dtype="float64")
        t_out, cols = _resample(t, [v], "tide", 3600.0)
        assert t_out[0] == t[0]
        assert np.allclose(np.diff(t_out), 3600.0)
        # bin means are unchanged by the stamping convention
        assert np.isclose(cols[0][0], v[:6].mean())

    def test_backup_creates_zip(self, server_project):
        import zipfile
        tmp_path, get, post = server_project
        (tmp_path / "z.grd").write_text("0 0\n1 1\n")
        (tmp_path / "aeolis.nc").write_bytes(b"fake output")
        (tmp_path / "gui" / "cache" / "scratch.bin").write_bytes(b"x")
        (tmp_path / "gui" / "rawdata" / "series.npz").write_bytes(b"raw")

        res = post("/api/project/backup",
                   {"include_rawdata": True, "include_outputs": False})
        job = _wait_job(get, res["job"])
        assert job["status"] == "done", job.get("error")
        zip_path = Path(job["result"]["file"])
        assert zip_path.is_file()
        assert zip_path.parent == tmp_path / "backups"
        names = set(zipfile.ZipFile(zip_path).namelist())
        assert "aeolis.txt" in names
        assert "z.grd" in names
        assert "gui/rawdata/series.npz" in names
        assert "aeolis.nc" not in names                      # outputs excluded
        assert not any(n.startswith("gui/cache/") for n in names)

        res = post("/api/project/backup",
                   {"include_rawdata": False, "include_outputs": True})
        job = _wait_job(get, res["job"])
        assert job["status"] == "done", job.get("error")
        names = set(zipfile.ZipFile(Path(job["result"]["file"])).namelist())
        assert "aeolis.nc" in names
        assert "gui/rawdata/series.npz" not in names
        # a second backup never archives an earlier one
        assert not any(n.startswith("backups") for n in names)

    def test_styles_seeded_and_roundtrip(self, server_project):
        tmp_path, get, post = server_project
        # no styles file yet -> seeded defaults
        names = {s["name"] for s in get("/api/styles")["styles"]}
        assert {"Elevation", "Bed level change", "Vegetation", "Grayscale"} <= names
        # full-replace roundtrip with normalization
        post("/api/styles", {"styles": [
            {"id": "x1", "name": "My style", "cmap": "viridis!r", "mode": "dots",
             "dotSize": 8, "opacity": 0.7, "min": None, "max": 12},
        ]})
        res = get("/api/styles")["styles"]
        assert res == [{"id": "x1", "name": "My style", "cmap": "viridis!r",
                        "mode": "dots", "dotSize": 8.0, "opacity": 0.7,
                        "min": None, "max": 12.0}]
        # invalid: every style needs id + name
        with pytest.raises(urllib.error.HTTPError):
            post("/api/styles", {"styles": [{"name": "no id"}]})
        # the valid store was not clobbered by the rejected write
        assert get("/api/styles")["styles"][0]["id"] == "x1"
