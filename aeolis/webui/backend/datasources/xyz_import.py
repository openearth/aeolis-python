"""Custom sample file import (*.xyz and friends).

Accepts whitespace- or comma-separated text with at least three columns
(x, y, z). Stores a fast-loading .npz copy in the project rawdata folder
next to the manifest entry; the original file is left untouched.
"""

from datetime import datetime, timezone
from pathlib import Path

import numpy as np


def import_file(path, dest_dir, crs=None):
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"file not found: {path}")

    try:
        data = np.loadtxt(path)
    except ValueError:
        data = np.loadtxt(path, delimiter=",")
    data = np.atleast_2d(data)
    if data.shape[1] < 3:
        raise ValueError("expected at least 3 columns (x, y, z)")

    x, y, z = data[:, 0], data[:, 1], data[:, 2]
    valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    x, y, z = x[valid], y[valid], z[valid]
    if x.size == 0:
        raise ValueError("no valid samples in file")

    out_name = f"xyz_{path.stem}.npz"
    np.savez_compressed(dest_dir / out_name, x=x, y=y, z=z.astype("float32"))

    return {
        "source": "xyz",
        "kind": "points",
        "path": f"gui/rawdata/{out_name}",
        "origin": str(path),
        "crs": crs,
        "res": None,
        "bounds": [float(x.min()), float(y.min()), float(x.max()), float(y.max())],
        "label": f"XYZ {path.name} ({x.size:,} pts)",
        "npoints": int(x.size),
        "downloaded": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M"),
    }
