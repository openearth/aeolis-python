"""Pre-run validation for the Run tab.

Goes beyond file existence: shape consistency of every configured 2D
spatial file against the grid, timeseries column counts (mirroring
aeolis.inout.check_configuration), time coverage of the simulation
window, and vegetation-framework file requirements.

Each check: {"key", "level": "ok"|"warn"|"error", "text"}.
"""

import numpy as np

from aeolis.webui.backend import grd_io
from aeolis.webui.backend.config_api import load_config

# config key -> (required columns, label)
TIMESERIES = {
    "wind_file": (3, "time, speed, direction"),
    "tide_file": (2, "time, water level"),
    "wave_file": (3, "time, Hs, Tp"),
    "meteo_file": (6, "time + 5 meteo columns"),
}

# 2D spatial files checked against the grid shape
SPATIAL_2D = ["bed_file", "ne_file", "veg_file", "threshold_file", "fence_file",
              "supply_file", "wave_mask", "tide_mask", "runup_mask",
              "threshold_mask", "gw_mask", "vver_mask"]

# per-species vegetation files: (ny+1, (nx+1) * nspecies)
SPATIAL_SPECIES = ["hveg_file", "Nt_file"]


def run_checks(project):
    checks = []

    def add(key, level, text):
        checks.append({"key": key, "level": level, "ok": level != "error", "text": text})

    try:
        values = load_config(project.configfile)
    except Exception as exc:  # noqa: BLE001
        add("config", "error", f"configuration unreadable: {exc}")
        return checks

    # --- time settings ---
    tstart = float(values.get("tstart") or 0)
    tstop = float(values.get("tstop") or 0)
    duration = tstop - tstart
    add("time", "ok" if duration > 0 else "error",
        f"simulation duration {duration:g} s" + ("" if duration > 0 else " (tstop <= tstart)"))
    dt = float(values.get("dt") or 0)
    if dt <= 0:
        add("dt", "error", f"dt = {dt:g} s must be positive")
    out_t = values.get("output_times")
    if out_t and dt > 0 and float(out_t) < dt:
        add("output_times", "warn", f"output_times ({out_t:g} s) is smaller than dt ({dt:g} s)")

    # --- grid ---
    grid_shape = None
    xf, yf = values.get("xgrid_file"), values.get("ygrid_file")
    if not xf or not yf:
        add("grid", "error", "xgrid_file / ygrid_file not set (Grid tab)")
    else:
        try:
            X = grd_io.read_grd(project.root / str(xf))
            Y = grd_io.read_grd(project.root / str(yf))
            if X.shape != Y.shape:
                add("grid", "error", f"x/y grid shapes differ: {X.shape} vs {Y.shape}")
            else:
                grid_shape = X.shape
                add("grid", "ok",
                    f"grid {X.shape[1] - 1} x {X.shape[0] - 1} cells ({xf}, {yf})")
                params = grd_io.derive_params(X, Y)
                if params and not params["uniform"]:
                    add("grid_uniform", "warn",
                        "grid cells are not equidistant/square - AeoLiS assumes dx = dy")
        except FileNotFoundError:
            add("grid", "error", f"grid file missing: {xf} / {yf}")
        except (ValueError, OSError) as exc:
            add("grid", "error", f"grid files unreadable: {exc}")

    # --- 2D spatial files vs grid shape ---
    method_veg = values.get("method_vegetation")
    for key in SPATIAL_2D:
        filename = values.get(key)
        if not filename:
            required = key == "bed_file"
            if required:
                add(key, "error", f"{key} not set (required)")
            continue
        if key == "veg_file" and method_veg == "grass":
            add(key, "warn", f"{key} is set but method_vegetation='grass' uses hveg/Nt files")
            continue
        path = project.root / str(filename)
        if not path.is_file():
            add(key, "error", f"{key} = {filename} (file missing)")
            continue
        try:
            Z = grd_io.read_grd(path)
        except (ValueError, OSError) as exc:
            add(key, "error", f"{key} = {filename} unreadable: {exc}")
            continue
        if grid_shape and Z.shape != grid_shape:
            add(key, "error",
                f"{key} = {filename} shape {Z.shape[0]}x{Z.shape[1]} does not match "
                f"the grid {grid_shape[0]}x{grid_shape[1]} - re-interpolate (Domain tab)")
        else:
            add(key, "ok", f"{key} = {filename}")

    # --- per-species vegetation files ---
    nspecies = len(values.get("species_names") or []) or 1
    for key in SPATIAL_SPECIES:
        filename = values.get(key)
        if not filename:
            if method_veg == "grass" and values.get("process_vegetation"):
                add(key, "error", f"{key} required for method_vegetation='grass'")
            continue
        path = project.root / str(filename)
        if not path.is_file():
            add(key, "error", f"{key} = {filename} (file missing)")
            continue
        try:
            Z = grd_io.read_grd(path)
        except (ValueError, OSError) as exc:
            add(key, "error", f"{key} = {filename} unreadable: {exc}")
            continue
        if grid_shape:
            ny1, nx1 = grid_shape
            if Z.shape[0] != ny1 or Z.shape[1] % nx1 != 0:
                add(key, "error",
                    f"{key} = {filename} shape {Z.shape[0]}x{Z.shape[1]} does not fit "
                    f"the grid {ny1}x{nx1} (x nspecies)")
            elif Z.shape[1] // nx1 != nspecies:
                add(key, "warn",
                    f"{key} holds {Z.shape[1] // nx1} species but species_names lists {nspecies}")
            else:
                add(key, "ok", f"{key} = {filename}")

    # --- timeseries ---
    for key, (ncols, label) in TIMESERIES.items():
        filename = values.get(key)
        required = key == "wind_file"
        if not filename:
            if required:
                add(key, "error", f"{key} not set (required)")
            continue
        path = project.root / str(filename)
        if not path.is_file():
            add(key, "error", f"{key} = {filename} (file missing)")
            continue
        try:
            data = np.atleast_2d(np.loadtxt(path))
        except (ValueError, OSError) as exc:
            add(key, "error", f"{key} = {filename} unreadable: {exc}")
            continue
        if data.shape[1] < ncols:
            add(key, "error",
                f"{key} = {filename} has {data.shape[1]} columns, expected >= {ncols} ({label})")
            continue
        t0, t1 = float(data[0, 0]), float(data[-1, 0])
        if t0 > tstart or t1 < tstop:
            add(key, "warn",
                f"{key} covers t = {t0:g} ... {t1:g} s but the simulation runs "
                f"{tstart:g} ... {tstop:g} s (series will be repeated/extrapolated)")
        else:
            add(key, "ok", f"{key} = {filename} ({data.shape[0]} rows)")

    # --- fractions ---
    grain_size = values.get("grain_size") or []
    grain_dist = values.get("grain_dist") or []
    try:
        if len(grain_size) != len(np.atleast_1d(np.asarray(grain_dist, dtype=float))):
            add("fractions", "warn",
                f"grain_size has {len(grain_size)} fraction(s) but grain_dist has "
                f"{len(np.atleast_1d(np.asarray(grain_dist, dtype=float)))} entries")
    except (TypeError, ValueError):
        pass

    return checks
