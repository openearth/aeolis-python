"""Flow-separation envelope computed directly on the model grid.

The separating streamline model is the one of :mod:`aeolis.separation`
(streamline curvature limits, look-ahead for upward curvature, detachment
at steep lee slopes, re-attachment after about 6 times the obstacle
height). There it is applied to each row of a computational grid that is
rotated into the wind. Here it is applied on the model grid itself, as a
semi-Lagrangian march along the wind:

* The grid is flipped and/or transposed so that the wind blows towards +i
  (columns) with a cross component of at most 45 degrees. Moving one column
  downwind then moves ``t = tan(angle) <= 1`` rows, and the step length
  along the streamline is ``h = ds / cos(angle)``.
* Every point on a streamline therefore lies exactly on a grid column and
  only needs linear interpolation between two rows. The upstream state of
  a cell (streamline height and slope) is interpolated from the two cells
  in the previous column.
* For a wind along a grid axis (``t = 0``) no interpolation happens at all
  and the result is identical to running the 1D core on every row.

Unlike the original core, slopes, curvatures and curvature increments are
in physical units (per metre), so the bubble does not change with the grid
spacing or with the step length of oblique streamlines. At a step length of
1 m the arithmetic is the same as in the original core.

Plain numpy in, numpy out; usable outside AeoLiS:

>>> zsep = separation_envelope(zb, theta=np.deg2rad(20.), ds=1., dn=1.)
>>> hsep = zsep - zb
"""

import numpy as np
from numba import njit


def separation_envelope(zb, theta, ds, dn, look_dist=50., k_press_up=0.05, k_crit_down=0.18,
                        s_crit=0.18, s_leeside=0.25, periodic=(False, False), skip_flat=True):
    '''Height of the separating streamline (>= zb) for a wind towards ``theta``.

    Parameters
    ----------
    zb : ndarray (ny, nx)
        Bed level [m].
    theta : float
        Direction the wind blows towards, counter-clockwise from the s axis
        (columns) of the grid [rad].
    ds, dn : float
        Grid spacing along the columns (s) and rows (n) [m].
    look_dist, k_press_up, k_crit_down, s_crit, s_leeside : float
        Parameters of the streamline model (``sep_*`` in the configuration).
    periodic : (bool, bool)
        Periodic boundaries in n and s. Streamlines that leave the grid
        through a periodic boundary re-enter on the other side. The march
        itself always starts attached at the upwind edge.
    skip_flat : bool
        Only treat streamlines that cross a steep lee slope or a sharp bend,
        like the row test of the original 2D implementation. Elsewhere the
        envelope equals the bed.

    Returns
    -------
    zsep : ndarray (ny, nx)
    '''
    zb = np.asarray(zb, dtype=float)
    ex, ey = np.cos(theta), np.sin(theta)
    per_n, per_s = bool(periodic[0]), bool(periodic[1])

    # canonical frame: march along +i (axis 1) with 0 <= t <= 1
    transpose = abs(ey) * ds > abs(ex) * dn  # compare index-space components
    Z = zb.T if transpose else zb
    a, b = (ey, ex) if transpose else (ex, ey)          # along-march, cross-march components
    d_along, d_cross = (dn, ds) if transpose else (ds, dn)
    per_cross = per_s if transpose else per_n
    flip_i, flip_j = a < 0., b < 0.
    if flip_i:
        Z = Z[:, ::-1]
    if flip_j:
        Z = Z[::-1, :]
    a, b = abs(a), abs(b)
    t = (b * d_along) / (a * d_cross) if a > 0. else 0.     # rows per column
    t = min(t, 1.)
    if t < 1e-12:            # wind along a grid axis (cos(pi/2) is 6e-17, not 0)
        t = 0.
    elif t > 1. - 1e-12:
        t = 1.
    h = d_along / a if a > 0. else d_along                  # step length along the streamline

    zs = _march(np.ascontiguousarray(Z), t, h, look_dist, k_press_up, k_crit_down,
                s_crit, s_leeside, per_cross, skip_flat)

    if flip_j:
        zs = zs[::-1, :]
    if flip_i:
        zs = zs[:, ::-1]
    if transpose:
        zs = zs.T
    return np.ascontiguousarray(zs)


def shear_reduction(hsep, tau_sep=0.5, slope=0.2):
    '''Factor (0..1) that reduces the bed shear stress below the envelope.

    Same form as :meth:`aeolis.shear.WindShear.separation_shear` (Duran 2010).
    '''
    delta = 1. / (slope * tau_sep)
    return np.minimum(np.maximum(1. - delta * hsep, 0.), 1.)


@njit(cache=True)
def _sample(A, jf, i, per):
    '''Linear interpolation of A[:, i] at fractional row jf; NaN outside.'''
    ny = A.shape[0]
    if per:
        jf = jf % ny
    elif jf < 0. or jf > ny - 1:
        return np.nan
    j0 = int(np.floor(jf))
    w = jf - j0
    if w == 0.:
        return A[j0, i]
    j1 = j0 + 1
    if j1 == ny:
        j1 = 0  # only reachable when periodic
    return A[j0, i] * (1. - w) + A[j1, i] * w


@njit(cache=True)
def _steep_flags(Z, t, h, k_crit_down, s_crit, per):
    '''Cells on a streamline that crosses a steep lee slope or a sharp bend.'''
    ny, nx = Z.shape
    steep = np.zeros((ny, nx), dtype=np.bool_)
    for i in range(2, nx - 1):
        for j in range(ny):
            zp = _sample(Z, j - t, i - 1, per)
            zpp = _sample(Z, j - 2. * t, i - 2, per)
            if np.isnan(zp) or np.isnan(zpp):
                continue
            s = (Z[j, i] - zp) / h
            s_prev = (zp - zpp) / h
            if s < -s_crit or (s_prev - s) / h < -k_crit_down:
                steep[j, i] = True

    # a streamline is flagged when any of its cells is: propagate down- and upwind
    fwd = steep.copy()
    for i in range(1, nx):
        for j in range(ny):
            jf = j - t
            if per:
                jf = jf % ny
            elif jf < 0. or jf > ny - 1:
                continue
            j0 = int(np.floor(jf))
            w = jf - j0
            j1 = j0 + 1 if j0 + 1 < ny else 0
            if fwd[j0, i - 1] or (w > 0. and fwd[j1, i - 1]):
                fwd[j, i] = True
    bwd = steep.copy()
    for i in range(nx - 2, -1, -1):
        for j in range(ny):
            jf = j + t
            if per:
                jf = jf % ny
            elif jf < 0. or jf > ny - 1:
                continue
            j0 = int(np.floor(jf))
            w = jf - j0
            j1 = j0 + 1 if j0 + 1 < ny else 0
            if bwd[j0, i + 1] or (w > 0. and bwd[j1, i + 1]):
                bwd[j, i] = True
    return fwd | bwd


@njit(cache=True)
def _march(Z, t, h, look_dist, k_press_up, k_crit_down, s_crit, s_leeside, per, skip_flat):
    ny, nx = Z.shape
    zs = Z.copy()                 # streamline height
    ss = np.zeros((ny, nx))       # streamline slope arriving at the cell
    if skip_flat:
        flag = _steep_flags(Z, t, h, k_crit_down, s_crit, per)
    else:
        flag = np.ones((ny, nx), dtype=np.bool_)

    look_n = int(look_dist / h)
    k_press_down_base = 0.05
    r_LH = 6.
    s_L = 1. / r_LH

    for i in range(1, nx):
        for j in range(ny):
            zc = Z[j, i]
            jp = j - t                                  # upstream point P in column i-1
            zbP = _sample(Z, jp, i - 1, per)
            if np.isnan(zbP):
                zs[j, i] = zc
                continue
            zbPP = _sample(Z, jp - t, i - 2, per) if i >= 2 else np.nan
            if np.isnan(zbPP) or not flag[j, i]:
                # attached: inflow boundary or a streamline without separation
                zs[j, i] = zc
                ss[j, i] = (zc - zbP) / h
                continue

            zP = _sample(zs, jp, i - 1, per)
            s_str = _sample(ss, jp, i - 1, per)
            gap = zP - zbP
            s_bed = (zbP - zbPP) / h
            s_bed_next = (zc - zbP) / h
            ds_bed = (s_bed_next - s_bed) / h

            # ----- 1. upward curvature anticipation (bed ahead of P) -----
            k_req_max = 0.
            zprev = zP
            for m in range(1, look_n + 1):
                if i - 1 + m > nx - 1:
                    break
                zm = zc if m == 1 else _sample(Z, jp + m * t, i - 1 + m, per)
                if np.isnan(zm):
                    break
                dxj = m * h
                dzj = zm - zP
                sbj = (zm - zprev) / h
                zprev = zm
                s_req_height = dzj / dxj
                k_req_height = (s_req_height - s_str) / dxj
                k_req_slope = (sbj - s_str) / dxj
                z_est = zP + s_str * dxj + 0.5 * k_req_slope * dxj * dxj
                if z_est > zm:
                    k_req_slope = 0.
                d_req_height = np.sqrt(2 * max(0., dzj - s_str * dxj) / k_press_up) if k_req_height > 0 else 0.
                d_req_slope = (sbj - s_str) / k_press_up if k_req_slope > 0 else 0.
                d_req = max(d_req_height, d_req_slope)
                if d_req > dxj:
                    k_req_max = max(k_req_max, k_req_slope, k_req_height)

            if k_req_max >= k_press_up:
                v_z = s_str + k_req_max * h
            else:
                # ----- 2. downward curvature by detachment and re-attachment -----
                if gap < 1e-6 and (ds_bed < -k_crit_down or s_str < -s_crit):
                    L_curv = max((s_str + s_leeside) / k_press_down_base, h)
                    # same rounding as the original int(i + L_curv / dx), with P in column i-1
                    m_curv = int((i - 1) + L_curv / h) - (i - 1)
                    # last downward crossing of a 1:6 line with the bed; without a
                    # crossing, up to where the streamline leaves the grid
                    n_int = -1
                    zmax = zbP
                    zmin = zbP
                    m = 1
                    while i - 1 + m <= nx - 1:
                        zm = _sample(Z, jp + m * t, i - 1 + m, per)
                        if np.isnan(zm):
                            break
                        if m >= m_curv:
                            zj = zP - s_L * (m - m_curv) * h
                            zj_prev = zP - s_L * (m - m_curv - 1) * h
                            if zj < zm and zj_prev >= zm:
                                n_int = m
                        m += 1
                    m_last = m - 1
                    if n_int < 0:
                        n_int = m_last
                    for m in range(1, n_int + 1):
                        zm = _sample(Z, jp + m * t, i - 1 + m, per)
                        zmax = max(zmax, zm)
                        zmin = min(zmin, zm)
                    L_attach = (zmax - zmin) * r_LH + L_curv
                    k_press_down = 1.5 / L_attach
                else:
                    k_press_down = k_press_down_base

                if ds_bed < -k_crit_down or s_bed < -s_crit or gap > 0.:
                    f_converge = np.exp(1.4 * (max(s_str + s_leeside, 0.) ** 0.85)) - 1.
                    v_z = s_str - k_press_down * f_converge * h
                else:
                    v_z = s_bed_next

            # ----- 3. advance the streamline -----
            z_new = max(zP + v_z * h, zc)
            zs[j, i] = z_new
            ss[j, i] = (z_new - zP) / h
    return zs
