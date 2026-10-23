"""
Utility functions for vegetation_grass.
Contains long, technical, or reusable logic only.
"""

import numpy as np
import logging
from numba import njit, prange

from aeolis.utils import *

logger = logging.getLogger(__name__)

def ensure_grass_parameters(p):
    """
    Check the length of grass parameters and make iterable over species.
    """ 

    ns = p['nspecies'] # based on length of 'species_names'

    param_lists = [
        'd_tiller', 'r_stem', 'alpha_uw', 'alpha_Nt', 'alpha_0',
        'G_h', 'G_c', 'G_s', 'Hveg', 'phi_h',
        'Nt_max', 'Nt_cap', 'R_cov',
        'lmax_c', 'mu_c', 'alpha_s', 'nu_s',
        'gamma_h', 'dzb_tol_c', 'dzb_tol_s',
        'dzb_opt_h', 'dzb_opt_c', 'dzb_opt_s',
        'beta_veg', 'm_veg', 'c1_okin', 'bounce', 'alpha_lift', 'gamma_Nt_decay', 'pNt_zeta'
        ]

    for param in param_lists:
        if param in p:

            # Check if the parameter is converted to a string
            if isinstance(p[param], str):
                p[param] = [float(t) for t in p[param].replace(',', ' ').split()]

            p[param] = makeiterable(p[param]) # Coverts to array

            # If the parameter length does not match nspecies, extend it
            if len(p[param]) != ns:
                if len(p[param]) == 1:
                    p[param] = np.full(ns, p[param][0]) 
                    logger.info(f"Parameter '{param}' extended to match nspecies={ns}.")
                else:
                    logger.error(f"Parameter '{param}' length does not match nspecies={ns}.")
                    raise ValueError(f"Parameter '{param}' length does not match nspecies={ns}.")
                
        else:
            logger.error(f"Parameter '{param}' is missing for grass model.")
            raise KeyError(f"Parameter '{param}' is missing for grass model.")
        
    # Competition density defaults to the maximum tiller density, which also scales the shelter
    p['Nt_cap'] = np.where(p['Nt_cap'] > 0., p['Nt_cap'], p['Nt_max'])

    # Handle alpha_comp separately as it is a matrix (nspecies x nspecies)
    if 'alpha_comp' in p:

        # Check if (elements of) alpha_comp is converted to a string
        if isinstance(p['alpha_comp'], str):
            p['alpha_comp'] = [float(t) for t in p['alpha_comp'].replace(',', ' ').split()]
        if isinstance(p['alpha_comp'], ndarray) and p['alpha_comp'].dtype.type is np.str_:
            p['alpha_comp'] = [float(item) for sublist in p['alpha_comp'] for item in sublist.replace(',', ' ').split()]

        p['alpha_comp'] = makeiterable(p['alpha_comp'])

        if p['alpha_comp'].size == 1:
            p['alpha_comp'] = np.zeros((ns, ns))
            np.fill_diagonal(p['alpha_comp'], 1.0)

        elif p['alpha_comp'].size == ns * ns:
            p['alpha_comp'] = p['alpha_comp'].reshape((ns, ns))

        else:
            logger.error(f"Parameter 'alpha_comp' length must be {ns*ns} ({ns} x {ns}).")
            raise ValueError(f"Parameter 'alpha_comp' length must be {ns*ns} ({ns} x {ns}).")
    else:
        p['alpha_comp'] = np.zeros((ns, ns))

    # --- Check intraspecific competition ---------------------------------------
    for k in range(ns):
        if p['alpha_comp'][k, k] != 1.0:
            logger.error(f"alpha_comp[{k},{k}] = {p['alpha_comp'][k,k]:.3f}. Expected value is 1.0 (equal to intraspecific competition).")
            raise ValueError(f"alpha_comp[{k},{k}] = {p['alpha_comp'][k,k]:.3f}. Expected value is 1.0 (equal to intraspecific competition).")

            
    return p
            

def generate_grass_subgrid(x, y, veg_res_factor):
    """
    Generate a refined vegetation subgrid inside each main grid cell.
    """

    f = veg_res_factor

    # Local grid vectors (rotation-safe)
    ex_x, ex_y = np.gradient(x, axis=1), np.gradient(y, axis=1)
    ey_x, ey_y = np.gradient(x, axis=0), np.gradient(y, axis=0)

    # Subcell offsets inside parent cell
    offs = (np.arange(f) + 0.5) / f - 0.5

    # Expand base grid
    x_veg = np.repeat(np.repeat(x, f, axis=1), f, axis=0)
    y_veg = np.repeat(np.repeat(y, f, axis=1), f, axis=0)

    for j in range(f):
        for i in range(f):
            x_veg[j::f, i::f] += offs[i] * ex_x + offs[j] * ey_x
            y_veg[j::f, i::f] += offs[i] * ex_y + offs[j] * ey_y

    return x_veg, y_veg


def expand_to_subgrid(A, f):
    """Expand (ny, nx, nspecies) → (ny*f, nx*f, nspecies) by replication."""
    return np.repeat(np.repeat(A, f, axis=0), f, axis=1)


@njit(cache=True)
def expand_to_subgrid_linear(A, f):
    """
    Expand (ny, nx, nspecies) → (ny*f, nx*f, nspecies)
    using bilinear interpolation on a uniform grid.
    """

    ny, nx, ns = A.shape
    nyf = ny * f
    nxf = nx * f

    Aout = np.empty((nyf, nxf, ns))

    # Precompute subcell weights
    w = np.empty(f)
    for i in range(f):
        w[i] = (i + 0.5) / f

    # Interior cells (bilinear)
    for iy in range(ny - 1):
        for ix in range(nx - 1):

            v00 = A[iy,     ix    ]
            v10 = A[iy,     ix + 1]
            v01 = A[iy + 1, ix    ]
            v11 = A[iy + 1, ix + 1]

            for jy in range(f):
                wy = w[jy]
                iyf = iy * f + jy

                for jx in range(f):
                    wx = w[jx]
                    ixf = ix * f + jx

                    for k in range(ns):
                        Aout[iyf, ixf, k] = (
                            (1.0 - wy) * ((1.0 - wx) * v00[k] + wx * v10[k]) +
                            wy         * ((1.0 - wx) * v01[k] + wx * v11[k])
                        )

    # Last column: linear in y
    for iy in range(ny - 1):
        for jy in range(f):
            iyf = iy * f + jy
            wy = w[jy]
            for k in range(ns):
                Aout[iyf, (nx - 1) * f:, k] = (
                    (1.0 - wy) * A[iy,     nx - 1, k] +
                    wy         * A[iy + 1, nx - 1, k]
                )

    # Last row: linear in x
    for ix in range(nx - 1):
        for jx in range(f):
            ixf = ix * f + jx
            wx = w[jx]
            for k in range(ns):
                Aout[(ny - 1) * f:, ixf, k] = (
                    (1.0 - wx) * A[ny - 1, ix,     k] +
                    wx         * A[ny - 1, ix + 1, k]
                )

    # Bottom-right corner
    for k in range(ns):
        Aout[(ny - 1) * f:, (nx - 1) * f:, k] = A[ny - 1, nx - 1, k]

    return Aout



def aggregate_from_subgrid(A, f):
    """Aggregate (ny*f, nx*f, nspecies) → (ny, nx, nspecies) by averaging."""
    nyf, nxf, ns = A.shape
    ny, nx = nyf // f, nxf // f
    return A.reshape(ny, f, nx, f, ns).mean(axis=(1,3))


def smooth_burial(s, p):
    """
    Compute smoothed burial rate over a trailing window T_burial.
    Stores bed levels over time and computes:
        dzb_veg = (zb(t) - zb(t - T_burial)) / T_burial   [m/s]
    """

    t = p['_time']
    T = p['T_burial']

    # Initialize history on first call
    if 'zb_hist' not in s:
        s['zb_hist'] = [(t, s['zb'].copy())]
        return np.zeros_like(s['zb'])

    # Append current state
    s['zb_hist'].append((t, s['zb'].copy()))

    # Drop states older than window
    t_min = t - T
    while s['zb_hist'][0][0] < t_min:
        s['zb_hist'].pop(0)
        
    t0, zb0 = s['zb_hist'][0] # Oldest retained bed level

    # Burial rate [m/s]
    dzb_smooth = (s['zb'] - zb0) / max(t - t0, 1e-12)

    # Ramp up smoothing over first T_burial
    T_rat = (t - t0) / T
    if T_rat < 1.0:
        dzb_smooth *= T_rat
    
    return dzb_smooth


@njit(cache=True)
def neighbourhood_average(Nt, R_cov, dx, wrap=False):
    """
    Compute neighbourhood average of Nt within radius R_cov.

    With ``wrap`` the rows are periodic (circular lateral boundary): the
    neighbourhood continues across the first and last row.
    """

    # Determine radius in grid cells
    r = int(np.ceil(R_cov / dx))
    ny = Nt.shape[0]

    # Loop over grid and compute average
    Nt_avg = np.zeros_like(Nt)
    for i in range(Nt.shape[0]):
        for j in range(Nt.shape[1]):
            ssum = 0.0
            cnt = 0 

            # Define neighbourhood bounds
            if wrap:
                i0 = i - r
                i1 = i + r + 1
            else:
                i0 = max(i - r, 0)
                i1 = min(i + r + 1, ny)
            j0 = max(j - r, 0)
            j1 = min(j + r + 1, Nt.shape[1])

            # Loop over neighbourhood
            for ii in range(i0, i1):
                for jj in range(j0, j1):
                    if ((ii - i)**2 + (jj - j)**2) * dx**2 <= R_cov**2:
                        ssum += Nt[ii % ny, jj]
                        cnt += 1

            # Compute average
            if cnt > 0:
                Nt_avg[i, j] = ssum / cnt

    return Nt_avg


def build_clonal_kernel(s_min, s_max, mu, dx):
    """ 
    Construct clonal dispersal kernel using a truncated Pareto distribution.
    """

    # Determine kernel radius in grid cells
    r = int(s_max / dx)

    # Create distance grid in physical space
    offsets = np.arange(-r, r + 1) * dx
    KX, KY = np.meshgrid(offsets, offsets)
    dist = np.sqrt(KX**2 + KY**2)

    # Set centre cell distance for internal recruitment
    dist[r, r] = dx / np.sqrt(6)

    # Initialise kernel
    kernel = np.zeros_like(dist)

    # Apply truncated Pareto distribution
    mask = (dist >= s_min) & (dist <= s_max)
    kernel[mask] = (
        mu
        / (s_min**(1 - mu) - s_max**(1 - mu))
        * dist[mask] ** (-mu)
    )

    # Normalize kernel
    kernel /= kernel.sum()

    return kernel, r


@njit(cache=True)
def apply_clonal_kernel(S_c, kernel, wrap=False):
    """
    Apply clonal dispersal kernel to clonal production rate S_c.

    With ``wrap`` the rows are periodic (circular lateral boundary), so no
    production is lost across the first and last row.
    """

    # Kernel radius
    r = kernel.shape[0] // 2
    ny = S_c.shape[0]

    # Initialize output
    dNt_clonal = np.zeros_like(S_c)

    # Loop over grid cells
    for i in range(S_c.shape[0]):
        for j in range(S_c.shape[1]):
            ssum = 0.0

            # Loop over kernel
            for ki in range(kernel.shape[0]):
                for kj in range(kernel.shape[1]):
                    ii = i + ki - r
                    jj = j + kj - r

                    if wrap:
                        if 0 <= jj < S_c.shape[1]:
                            ssum += kernel[ki, kj] * S_c[ii % ny, jj]
                    elif 0 <= ii < ny and 0 <= jj < S_c.shape[1]:
                        ssum += kernel[ki, kj] * S_c[ii, jj]

            # Assign spreaded value
            dNt_clonal[i, j] = ssum

    return dNt_clonal

@njit(cache=True)
def sample_seed_germination(S_s, B_s, a_s, nu_s, dx, wrap=False):
    """
    Sample stochastic seed germination events from seed production rate S_s.

    With ``wrap`` the rows are periodic (circular lateral boundary): a seed
    that lands beyond the first or last row enters on the other side.
    """

    ny, nx = S_s.shape
    dNt_seed = np.zeros_like(S_s)

    lambda_total = S_s.sum()
    if lambda_total <= 0.0:
        return dNt_seed

    n_seeds = np.random.poisson(lambda_total)
    if n_seeds == 0:
        return dNt_seed

    # Flatten source strengths
    flat = S_s.ravel()

    # Build cumulative distribution
    cdf = np.empty(flat.size)
    csum = 0.0
    for i in range(flat.size):
        csum += flat[i]
        cdf[i] = csum

    for _ in range(n_seeds):

        # --- Weighted source selection (CDF sampling) ---
        u = np.random.rand() * cdf[-1]
        j = 0
        while cdf[j] < u:
            j += 1

        iy = j // nx
        ix = j - iy * nx

        # --- Sample jump distance (2Dt) ---
        u = np.random.rand()
        r = np.sqrt(a_s * (u ** (-1.0 / (nu_s - 1.0)) - 1.0))

        # --- Sample direction ---
        theta = 2.0 * np.pi * np.random.rand()

        dy = int(round((np.sin(theta) * r) / dx))
        dx_ = int(round((np.cos(theta) * r) / dx))

        yy = iy + dy
        xx = ix + dx_
        if wrap:
            yy = yy % ny

        if 0 <= yy < ny and 0 <= xx < nx:
            # Check burial probability at the landing site (yy, xx)
            if np.random.rand() <= B_s[yy, xx]: 
                dNt_seed[yy, xx] += 1

    return dNt_seed


def weight_parameters(p):
    """Per-species parameters of :func:`shear_reduction_weights`, as float arrays."""
    f = lambda k: np.asarray(p[k], dtype=np.float64)
    d = f('d_tiller')
    return (f('r_stem'), f('alpha_uw'), f('alpha_Nt'), f('alpha_0'), d, (d / 2.0)**2,
            f('Hveg'), f('Nt_max'), f('m_veg') * f('beta_veg'), f('c1_okin'))


@njit(cache=True)
def shear_reduction_weights(uw, Nt, hveg, r_stem, alpha_uw, alpha_Nt, alpha_0, d_tiller, d_sq,
                            Hveg, Nt_max, mbeta, c1):
    """Bending, main-grid vegetation metrics and the weighted means of the shear reduction.

    Per species k: fbend = clip(r + (1 - r)(a_uw uw + a_Nt Nt + a_0), 0, 1),
    hvegeff = hveg fbend, lamveg = Nt hvegeff d, rhoveg = Nt pi (d/2)^2,
    R0 = 1 / sqrt(1 + m beta lamveg) and the weight w = (hveg / Hveg)(Nt / Nt_max).
    Returns fbend, hvegeff, lamveg, rhoveg (per species) and the w-weighted
    means R0veg (1 where w sums to 0), L_decay = hvegeff / c1, hvegeff and Nt / Nt_max.
    Each expression is evaluated in the order of the numpy version
    (:func:`shear_reduction_weights_reference`), so the results are identical.
    """
    ny, nx, ns = Nt.shape
    fbend = np.empty((ny, nx, ns))
    hvegeff = np.empty((ny, nx, ns))
    lamveg = np.empty((ny, nx, ns))
    rhoveg = np.empty((ny, nx, ns))
    R0veg = np.ones((ny, nx))
    L_decay = np.zeros((ny, nx))
    hv = np.zeros((ny, nx))
    rNt = np.zeros((ny, nx))
    for i in range(ny):
        for j in range(nx):
            w_sum = 0.0
            n_R0 = 0.0
            n_L = 0.0
            n_h = 0.0
            n_rNt = 0.0
            for k in range(ns):
                nt = Nt[i, j, k]
                fb = r_stem[k] + (1.0 - r_stem[k]) * (alpha_uw[k] * uw[i, j] + alpha_Nt[k] * nt + alpha_0[k])
                fb = min(max(fb, 0.0), 1.0)
                he = hveg[i, j, k] * fb
                lam = nt * he * d_tiller[k]
                fbend[i, j, k] = fb
                hvegeff[i, j, k] = he
                lamveg[i, j, k] = lam
                rhoveg[i, j, k] = (nt * np.pi) * d_sq[k]
                w = (hveg[i, j, k] / Hveg[k]) * (nt / Nt_max[k])
                R0 = 1.0 / np.sqrt(1.0 + mbeta[k] * lam)
                w_sum += w
                n_R0 += w * R0
                n_L += w * (he / c1[k])
                n_h += w * he
                n_rNt += w * (nt / Nt_max[k])
            if w_sum != 0.0:
                R0veg[i, j] = n_R0 / w_sum
                L_decay[i, j] = n_L / w_sum
                hv[i, j] = n_h / w_sum
                rNt[i, j] = n_rNt / w_sum
    return fbend, hvegeff, lamveg, rhoveg, R0veg, L_decay, hv, rNt


def shear_reduction_weights_reference(uw, Nt, hveg, r_stem, alpha_uw, alpha_Nt, alpha_0, d_tiller, d_sq,
                                      Hveg, Nt_max, mbeta, c1):
    """The numpy form of :func:`shear_reduction_weights` (the code it replaced), for the tests."""
    ny, nx, ns = Nt.shape
    fbend = np.ones_like(hveg)
    for k in range(ns):
        fbend[:, :, k] = (r_stem[k] + (1.0 - r_stem[k]) * (alpha_uw[k] * uw + alpha_Nt[k] * Nt[:, :, k] + alpha_0[k]))
        fbend[:, :, k] = np.clip(fbend[:, :, k], 0.0, 1.0)
    hvegeff = hveg * fbend
    lamveg = Nt * hvegeff * d_tiller
    rhoveg = Nt * np.pi * d_sq
    w_sum, n_R0, n_L, n_h, n_rNt = (np.zeros((ny, nx)) for _ in range(5))
    R0veg, L_decay, hv, rNt = np.ones((ny, nx)), np.zeros((ny, nx)), np.zeros((ny, nx)), np.zeros((ny, nx))
    for k in range(ns):
        w = (hveg[:, :, k] / Hveg[k]) * (Nt[:, :, k] / Nt_max[k])
        R0_k = 1.0 / np.sqrt(1.0 + mbeta[k] * lamveg[:, :, k])
        w_sum += w
        n_R0 += w * R0_k
        n_L += w * (hvegeff[:, :, k] / c1[k])
        n_h += w * hvegeff[:, :, k]
        n_rNt += w * (Nt[:, :, k] / Nt_max[k])
    ix = w_sum != 0.0
    R0veg[ix] = n_R0[ix] / w_sum[ix]
    L_decay[ix] = n_L[ix] / w_sum[ix]
    hv[ix] = n_h[ix] / w_sum[ix]
    rNt[ix] = n_rNt[ix] / w_sum[ix]
    return fbend, hvegeff, lamveg, rhoveg, R0veg, L_decay, hv, rNt


@njit(cache=True)
def _okin_sources(R0, udir, L_decay, dx, wind_convention_id):
    """Wake direction, wake length and search radius of every vegetated source cell."""
    ny, nx = R0.shape
    deg2rad = np.pi / 180.0
    R_end = 0.999

    # source-cell properties (identical expressions to compute_okin_reduction)
    active = np.zeros((ny, nx), dtype=np.bool_)
    ux = np.zeros((ny, nx))
    uy = np.zeros((ny, nx))
    L_end = np.zeros((ny, nx))
    rad = np.zeros((ny, nx), dtype=np.int64)
    rmax = -1
    for iy in range(ny):
        for ix in range(nx):
            R0_ij = R0[iy, ix]
            if R0_ij >= 1.0:
                continue
            L_ij = L_decay[iy, ix]
            if L_ij <= 0.0:
                continue
            if wind_convention_id == 0:
                th = (270.0 - udir[iy, ix]) * deg2rad
            else:
                th = udir[iy, ix] * deg2rad
            ux[iy, ix] = np.cos(th)
            uy[iy, ix] = np.sin(th)
            L_end[iy, ix] = -L_ij * np.log((1.0 - R_end) / (1.0 - R0_ij))
            rad[iy, ix] = int(L_end[iy, ix] / dx) + 2
            active[iy, ix] = True
            if rad[iy, ix] > rmax:
                rmax = rad[iy, ix]
    return active, ux, uy, L_end, rad, rmax


@njit(parallel=True, cache=True)
def compute_okin_reduction_periodic(x, y, R0, udir, L_decay, hvegeff, rNt, wind_convention_id,
                                    row_dx, row_dy):
    """
    Okin leeside shear reduction on a grid whose rows repeat (circular lateral boundary).

    As :func:`compute_okin_reduction`, but a wake also reaches targets beyond
    the first and last row: the rows are continued periodically, row ``iv``
    being row ``iv % ny`` moved by ``iv // ny`` periods. The position of a
    source relative to a target is the column offset within the target's row
    plus a whole number of row steps (``row_dx``, ``row_dy``, the vector from
    one row to the next), so a single row with identical neighbours gives the
    same wake as any number of identical rows.
    """

    ny, nx = x.shape
    R = np.ones((ny, nx))
    hvegeff_zeta = np.zeros((ny, nx))
    rNt_zeta = np.zeros((ny, nx))
    Rveg_zeta = np.ones((ny, nx))

    dx = np.sqrt((x[0, 1] - x[0, 0])**2 + (y[0, 1] - y[0, 0])**2)
    active, ux, uy, L_end, rad, rmax = _okin_sources(R0, udir, L_decay, dx, wind_convention_id)

    if rmax < 0:
        return R, hvegeff_zeta, rNt_zeta, Rveg_zeta

    for jy in prange(ny):
        for jx in range(nx):
            Rt = 1.0
            ht = 0.0
            nt = 0.0
            zt = 1.0
            for iv in range(jy - rmax, jy + rmax + 1):
                iy = iv % ny
                ox = (jy - iv) * row_dx
                oy = (jy - iv) * row_dy
                for ix in range(max(0, jx - rmax), min(nx, jx + rmax + 1)):
                    if not active[iy, ix]:
                        continue
                    r = rad[iy, ix]
                    if jy < iv - r or jy > iv + r or jx < ix - r or jx > ix + r:
                        continue

                    rx = (x[jy, jx] - x[jy, ix]) + ox
                    ry = (y[jy, jx] - y[jy, ix]) + oy
                    s = rx * ux[iy, ix] + ry * uy[iy, ix]
                    if s <= 0.0 or s > L_end[iy, ix]:
                        continue
                    d_perp = np.sqrt((rx - s * ux[iy, ix])**2 + (ry - s * uy[iy, ix])**2)
                    if d_perp >= dx:
                        continue

                    R0_ij = R0[iy, ix]
                    R_s = 1.0 - (1.0 - R0_ij) * np.exp(-s / L_decay[iy, ix])
                    w = 1.0 - d_perp / dx
                    R_loc = 1.0 - w * (1.0 - R_s)

                    if R_loc < Rt:
                        Rt = R_loc
                        ht = hvegeff[iy, ix]
                        nt = rNt[iy, ix]
                        zt = (1 - R_loc) / (1 - R0_ij)
            R[jy, jx] = Rt
            hvegeff_zeta[jy, jx] = ht
            rNt_zeta[jy, jx] = nt
            Rveg_zeta[jy, jx] = zt

    return R, hvegeff_zeta, rNt_zeta, Rveg_zeta


@njit(cache=True)
def compute_okin_reduction(x, y, R0, udir, L_decay, hvegeff, rNt, wind_convention_id):
    """
    Compute Okin leeside shear reduction using an effective decay length.

    Every vegetated source cell reduces the shear in a wake downwind of it;
    the strongest reduction wins. The loop runs over the sources, each over
    its own wake. A parallel form that runs over the targets (e5df619) had to
    search the largest wake radius of the grid around every target: on the
    blowout it visited 2.8 times as many cells and took twice as long on one
    core, and was no faster on four.

    Parameters
    ----------
    x, y               : 2D arrays of cell-center coordinates [m]
    R0                 : 2D array of local Raupach reduction [-]
    udir               : 2D array of wind direction [deg]
    L_decay            : 2D array of Okin decay length h/c1 [m]
    wind_convention_id : int
                          0 = nautical (from North, clockwise)
                          1 = cartesian (0° = +x, CCW)

    Returns
    -------
    R : 2D array of shear reduction including leeside effects [-]
    hvegeff_zeta, rNt_zeta, Rveg_zeta : 2D arrays of the strongest wake per target
    """

    ny, nx = x.shape
    R = np.ones((ny, nx))

    hvegeff_zeta = np.zeros((ny, nx)) # Representative vegetation height for zeta adjustment
    rNt_zeta = np.zeros((ny, nx)) # Relative Nt for zeta adjustment in vegetation wake
    Rveg_zeta = np.ones((ny, nx)) # Relative shear reduction for zeta adjustment in vegetation wake

    # grid spacing (assumed uniform)
    dx = np.sqrt((x[0, 1] - x[0, 0])**2 + (y[0, 1] - y[0, 0])**2)

    deg2rad = np.pi / 180.0
    R_end = 0.999

    # Loop over source cells
    for iy in range(ny):
        for ix in range(nx):

            R0_ij = R0[iy, ix]
            hveg0 = hvegeff[iy, ix]  # Effective vegetation height at source cell
            rNt0 = rNt[iy, ix]       # rNt (=Nt/Nt_max) at source cell

            if R0_ij >= 1.0:
                continue

            L_ij = L_decay[iy, ix]
            if L_ij <= 0.0:
                continue

            # Local wind direction at source cell
            if wind_convention_id == 0:      # nautical
                th = (270.0 - udir[iy, ix]) * deg2rad
            else:                            # cartesian
                th = udir[iy, ix] * deg2rad

            ux = np.cos(th)
            uy = np.sin(th)

            # Max downwind distance (R = R_end)
            L_end = -L_ij * np.log((1.0 - R_end) / (1.0 - R0_ij))
            r = int(L_end / dx) + 2  # margin

            # Compute window bounds
            j0 = max(0, iy - r)
            j1 = min(ny, iy + r + 1)
            i0 = max(0, ix - r)
            i1 = min(nx, ix + r + 1)

            # Loop over target cells in window
            x0 = x[iy, ix]
            y0 = y[iy, ix]
            for jy in range(j0, j1):
                for jx in range(i0, i1):

                    rx = x[jy, jx] - x0
                    ry = y[jy, jx] - y0

                    # Along-wind distance
                    s = rx * ux + ry * uy
                    if s <= 0.0 or s > L_end:
                        continue

                    # Perpendicular distance
                    d_perp = np.sqrt((rx - s * ux)**2 + (ry - s * uy)**2)
                    if d_perp >= dx:
                        continue

                    # Possible debug visualization
                    # debug_okin_geometry(x, y, iy, ix, jy, jx, ux, uy, L_end, W)

                    # Okin along-wind reduction (using decay length)
                    R_s = 1.0 - (1.0 - R0_ij) * np.exp(-s / L_ij)

                    # Cross-wind triangular weighting (width = 2*dx)
                    w = 1.0 - d_perp / dx
                    R_loc = 1.0 - w * (1.0 - R_s)

                    # Strongest reduction wins
                    if R_loc < R[jy, jx]:
                        R[jy, jx] = R_loc
                        hvegeff_zeta[jy, jx] = hveg0
                        rNt_zeta[jy, jx] = rNt0 
                        Rveg_zeta[jy, jx] = (1 - R_loc) / (1 - R0_ij)

    return R, hvegeff_zeta, rNt_zeta, Rveg_zeta


def debug_okin_geometry(x, y, iy, ix, jy, jx, ux, uy, L_end, W):
    """
    Visualize Okin geometry for one source cell (iy,ix)
    and one target cell (jy,jx).
    """
    import matplotlib.pyplot as plt

    # Extract coordinates
    x0 = x[iy, ix]
    y0 = y[iy, ix]
    xt = x[jy, jx]
    yt = y[jy, jx]
    rx = xt - x0
    ry = yt - y0
    s = rx * ux + ry * uy

    # Projection point on ray
    xp = x0 + s * ux
    yp = y0 + s * uy
    d_perp = np.sqrt((rx - s * ux)**2 + (ry - s * uy)**2)

    # Plotting
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(x, y, s=10, c='lightgray', label='grid') # grid points
    ax.scatter(x0, y0, c='red', s=80, label='source') # source cell
    ax.scatter(xt, yt, c='blue', s=80, label='target') # target cell
    ax.plot([x0, x0 + L_end * ux],[y0, y0 + L_end * uy],'r--', lw=2, label='wind ray')
    ax.plot([xt, xp], [yt, yp], 'k:', lw=2, label='d_perp')
    ax.set_aspect('equal')
    ax.set_title(
        f"s = {s:.2f}, d_perp = {d_perp:.2f}\n"
        f"rx = {rx:.2f}, ry = {ry:.2f}"
    )
    ax.legend()
    plt.show()


