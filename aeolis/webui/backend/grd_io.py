"""Reading, writing and generating AeoLiS rectangular grids (.grd).

Grid convention (verified against ``aeolis.inout.visualize_grid`` and
``aeolis.model.AeoLiS.initialize``):

- ``x.grd``/``y.grd`` are ASCII matrices of shape (ny+1, nx+1); element
  [j, i] is the same physical point in every grid file.
- Columns (i) run cross-shore: column 0 is the OFFSHORE boundary, the
  last column is the ONSHORE boundary. Rows (j) run alongshore: the
  first and last rows are the LATERAL boundaries.
- Cells are square and equidistant (dx = dy); orientation is embedded
  in the coordinates themselves (config ``alfa`` stays 0).

Generation follows the documented protocol::

    X, Y = np.meshgrid(arange(nx+1)*dx, arange(ny+1)*dx)
    Xr = X*cos(t) - Y*sin(t) + x0
    Yr = X*sin(t) + Y*cos(t) + y0

with ``t`` the counter-clockwise rotation of the local (cross-shore)
x-axis w.r.t. the model-CRS east axis, and (x0, y0) the (0,0) corner,
i.e. the offshore/lateral-A corner.
"""

import math

import numpy as np


def generate(x0, y0, dx, nx, ny, rotation_deg):
    """Return (X, Y) coordinate matrices of shape (ny+1, nx+1)."""
    nx = int(nx)
    ny = int(ny)
    if nx < 1 or ny < 1:
        raise ValueError("nx and ny must be >= 1")
    if dx <= 0:
        raise ValueError("dx must be > 0")
    x = np.arange(nx + 1, dtype=float) * dx
    y = np.arange(ny + 1, dtype=float) * dx
    X, Y = np.meshgrid(x, y)
    t = math.radians(rotation_deg)
    Xr = X * math.cos(t) - Y * math.sin(t) + x0
    Yr = X * math.sin(t) + Y * math.cos(t) + y0
    return Xr, Yr


def write_grd(path, matrix):
    # 11 significant digits: RD-scale coordinates (~1e5-1e6 m) keep
    # sub-millimetre precision through the ASCII roundtrip
    np.savetxt(path, np.atleast_2d(matrix), fmt="%0.10e")


def read_grd(path):
    try:
        return np.atleast_2d(np.loadtxt(path))
    except ValueError:
        # mask files may hold complex values (real=multiplier, imag=offset)
        return np.atleast_2d(np.loadtxt(path, dtype=complex))


def derive_params(X, Y):
    """Recover {x0, y0, dx, nx, ny, rotation} from coordinate matrices.
    Returns None for non-equidistant/curvilinear grids."""
    ny1, nx1 = X.shape
    if nx1 < 2 or ny1 < 2:
        return None
    dxi = float(np.hypot(X[0, 1] - X[0, 0], Y[0, 1] - Y[0, 0]))
    rotation = math.degrees(math.atan2(Y[0, 1] - Y[0, 0], X[0, 1] - X[0, 0]))
    # equidistance check (tolerant)
    steps_x = np.hypot(np.diff(X, axis=1), np.diff(Y, axis=1))
    steps_y = np.hypot(np.diff(X, axis=0), np.diff(Y, axis=0))
    uniform = (
        np.allclose(steps_x, dxi, rtol=1e-3, atol=1e-6)
        and np.allclose(steps_y, dxi, rtol=1e-3, atol=1e-6)
    )
    return {
        "x0": float(X[0, 0]),
        "y0": float(Y[0, 0]),
        "dx": dxi,
        "nx": nx1 - 1,
        "ny": ny1 - 1,
        "rotation": rotation,
        "uniform": bool(uniform),
    }


def _mid(a, b):
    return [(a[0] + b[0]) / 2, (a[1] + b[1]) / 2]


def geometry(X, Y, max_lines=60):
    """Displayable geometry: outline ring, decimated grid lines and
    boundary midpoints/corners, all in model coordinates."""
    ny1, nx1 = X.shape

    def _pt(j, i):
        return [float(X[j, i]), float(Y[j, i])]

    outline = (
        [_pt(0, i) for i in range(nx1)]
        + [_pt(j, nx1 - 1) for j in range(1, ny1)]
        + [_pt(ny1 - 1, i) for i in range(nx1 - 2, -1, -1)]
        + [_pt(j, 0) for j in range(ny1 - 2, 0, -1)]
    )

    stride_i = max(1, (nx1 - 1) // max_lines)
    stride_j = max(1, (ny1 - 1) // max_lines)
    lines = []
    for i in range(0, nx1, stride_i):
        lines.append([_pt(j, i) for j in range(ny1)])
    for j in range(0, ny1, stride_j):
        lines.append([_pt(j, i) for i in range(nx1)])

    boundaries = {
        # matches inout.visualize_grid: column 0 = offshore
        "offshore": _mid(_pt(0, 0), _pt(ny1 - 1, 0)),
        "onshore": _mid(_pt(0, nx1 - 1), _pt(ny1 - 1, nx1 - 1)),
        "lateral_a": _mid(_pt(0, 0), _pt(0, nx1 - 1)),
        "lateral_b": _mid(_pt(ny1 - 1, 0), _pt(ny1 - 1, nx1 - 1)),
    }
    corners = {
        "c00": _pt(0, 0),
        "c0n": _pt(0, nx1 - 1),
        "cm0": _pt(ny1 - 1, 0),
        "cmn": _pt(ny1 - 1, nx1 - 1),
    }
    return {
        "outline": outline,
        "lines": lines,
        "boundaries": boundaries,
        "corners": corners,
        "shape": [ny1, nx1],
    }


def shear_preview(X, Y, dx_c, dy_c, buffer_width, udir):
    """Geometry of the rotating computational (shear) grid for a given
    wind direction, mirroring shear.WindShear.set_computational_grid:
    a rectangle centred on the grid centroid whose extents are the
    corner projections onto the wind-parallel/perpendicular axes plus
    2 x buffer_width."""
    x0, y0 = float(np.mean(X)), float(np.mean(Y))
    corners = np.array([
        [X[0, 0], Y[0, 0]], [X[-1, 0], Y[-1, 0]],
        [X[0, -1], Y[0, -1]], [X[-1, -1], Y[-1, -1]],
    ])
    # wind vector (nautical: direction wind comes from)
    t = math.radians(udir)
    par = np.array([-math.sin(t), -math.cos(t)])     # wind-parallel
    perp = np.array([-math.sin(t - math.pi / 2), -math.cos(t - math.pi / 2)])

    rel = corners - [x0, y0]
    proj_par = rel @ par
    proj_perp = rel @ perp
    length = float(proj_par.max() - proj_par.min()) + 2 * buffer_width
    width = float(proj_perp.max() - proj_perp.min()) + 2 * buffer_width

    half_l, half_w = length / 2, width / 2
    ring = []
    for sl, sw in [(-1, -1), (1, -1), (1, 1), (-1, 1)]:
        pt = np.array([x0, y0]) + par * (sl * half_l) + perp * (sw * half_w)
        ring.append([float(pt[0]), float(pt[1])])

    inner = []
    if buffer_width > 0:
        hl, hw = half_l - buffer_width, half_w - buffer_width
        for sl, sw in [(-1, -1), (1, -1), (1, 1), (-1, 1)]:
            pt = np.array([x0, y0]) + par * (sl * hl) + perp * (sw * hw)
            inner.append([float(pt[0]), float(pt[1])])

    n_cells = [max(1, int(round(length / dx_c))), max(1, int(round(width / dy_c)))]
    return {
        "ring": ring,
        "inner": inner,
        "center": [x0, y0],
        "length": length,
        "width": width,
        "n_cells": n_cells,
        "udir": udir,
    }
