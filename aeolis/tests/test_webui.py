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
        schema = build_schema()
        keys = {p["key"] for s in schema["sections"] for p in s["params"]}
        assert keys == set(DEFAULT_CONFIG)

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
def server_project(tmp_path):
    import aeolis.inout
    from aeolis.webui.backend.httpd import make_server

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

    def test_unknown_route_404(self, server_project):
        tmp_path, get, post = server_project
        with pytest.raises(urllib.error.HTTPError) as err:
            get("/api/nonexistent")
        assert err.value.code == 404
