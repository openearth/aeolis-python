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


# ---------------------------------------------------------------------
# HTTP API integration (config + grid roundtrip on a temp project)
# ---------------------------------------------------------------------

@pytest.fixture()
def server_project(tmp_path, monkeypatch):
    import aeolis.inout
    from aeolis.webui.backend import settings
    from aeolis.webui.backend.httpd import make_server

    # keep test projects out of the user's real recent-projects list
    monkeypatch.setattr(settings, "RECENT_FILE", tmp_path / "recent.json")

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

    def test_unknown_route_404(self, server_project):
        tmp_path, get, post = server_project
        with pytest.raises(urllib.error.HTTPError) as err:
            get("/api/nonexistent")
        assert err.value.code == 404
