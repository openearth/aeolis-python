"""
Streamline-based separation bubble model for AeoLiS

This module provides a unified 1D/2D separation bubble computation based on
a streamline curvature model. The physics is contained in a single numba-
compiled 1D core. For 2D grids (FFT shear mode), the 1D streamline model
is applied to each wind-aligned row.

Public Entry Point:
    compute_separation(z_bed, dx, wind_sign, p)

Internal:
    _compute_separation_1d(...)
    _compute_separation_2d(...)
    _streamline_core_1d(...)   # numba auto-compiled

"""

import numpy as np
from numba import njit

# Wrapper for 1D vs 2D
def compute_separation(p, z_bed, dx, udir=0):

    # Get all parameters
    look_dist = p['sep_look_dist']
    k_press_up = p['sep_k_press_up']
    k_crit_down = p['sep_k_crit_down']
    s_crit = p['sep_s_crit']
    s_leeside = p['sep_s_leeside']

    # Determine size
    ny_, nx = z_bed.shape

    # Flip the bed for negative wind directions (only in 1D)
    if udir < 0 and ny_ == 0:
            z = z[::-1]

    # Run streamline computation once in 1D
    if ny_ == 0:
        zsep = _streamline_core_1d(
            z_bed, dx, look_dist,       # Base input
            k_press_up, k_crit_down,   # Curviture 
            s_crit, s_leeside)          # Leeside slopes

    
    else: # Run streamline computation over every transect in 2D (Numba-compiled over 1D)
        zsep = _compute_separation_2d(
            z_bed, dx, look_dist,
            k_press_up, k_crit_down,
            s_crit, s_leeside)

    # Flip back for the 1D case
    if udir < 0 and ny_ == 0:
            zsep = zsep[::-1]

    return zsep


@njit
def _compute_separation_2d(
    z_bed, dx, look_dist,
    k_press_up, k_crit_down,
    s_crit, s_leeside):

    ny_, nx = z_bed.shape
    zsep = np.zeros_like(z_bed)

    for j in range(ny_):

        row = z_bed[j, :]

        # Quick skip: slope + curvature checks
        skip = True
        for i in range(1, nx-1):
            s  = (row[i]   - row[i-1]) / dx
            s_prev = (row[i-1] - row[i-2]) / dx
            ds = s_prev - s  # curvature approx

            if (s < -s_crit) or (ds < -k_crit_down):
                skip = False
                break

        if skip:
            zsep[j, :] = row
        else:
            zsep[j, :] = _streamline_core_1d(
                row, dx, look_dist,
                k_press_up, k_crit_down,
                s_crit, s_leeside)

    return zsep


# Streamline core
@njit
def _streamline_core_1d(
    z_bed, dx, look_dist,
    k_press_up, k_crit_down,
    s_crit, s_leeside):

    # Initialize size and z (= streamline)
    N = z_bed.size
    z = z_bed.copy()

    # Loop over all points in windward direction
    for i in range(1, N - 1):

        # Compute slope of streamline
        s_str = (z[i] - z[i-1]) / dx
        gap = z[i] - z_bed[i]

        # Compute slope of bed
        s_bed = (z_bed[i] - z_bed[i-1]) / dx
        s_bed_next = (z_bed[i+1] - z_bed[i]) / dx
        ds_bed = s_bed_next - s_bed

        # Determine how far to look ahead for upward curviture
        look_n = int(look_dist / dx)
        i_end  = min(i + look_n, N - 1)

        # Initialize maximum required curviture (k_req_max) and the resulting slope (v_z)
        k_req_max = 0.0
        v_z = None
        k_press_down_base = 0.05

        # ----- 1. UPWARD CURVATURE ANTICIPATION -----

        # Start looking forward for upward curvature anticipation
        if i_end > i:
            for j in range(i + 1, i_end + 1):

                # Compute difference in distance, height and slope downwind (forward = j)
                dxj = (j - i) * dx          # Total distance between current (i) and forward (j)
                dzj = z[j] - z[i]           # Total height difference between current (i) and forward (j)
                sbj = (z[j] - z[j-1]) / dx  # Slope at forward (j)

                # Required slope (s) to close height difference from (i) to (j)
                s_req_height = dzj / dxj

                # Required curviture (k) to close slope difference for height (k_req_height) and slope (k_req_slope)
                k_req_height = (s_req_height - s_str) / dxj
                k_req_slope  = (sbj - s_str) / dxj

                # Prevent that the streamline "overshoots" the bed level due a too steep curviture
                z_est = z[i] + s_str * dxj + 0.5 * k_req_slope * dxj * dxj
                if z_est > z[j]:
                    k_req_slope = 0.0

                # Required distance to reach either height or slope
                d_req_height = np.sqrt(2 * max(0.0, dzj - s_str * dxj) / k_press_up) if k_req_height > 0 else 0.
                d_req_slope  = (sbj - s_str) / k_press_up if k_req_slope > 0 else 0.
                d_req = max(d_req_height, d_req_slope)

                # Check whether d_req is within reach (pass if dxj > d_req)
                # i.e., if d_req < dxj; it is not necessary to bend upward yet
                if d_req > dxj:
                    k_req_max = max(k_req_max, k_req_slope, k_req_height)

        # Apply curvature anticipation
        if k_req_max >= k_press_up: 
            v_z = s_str + k_req_max 

        # ----- 2. DOWNWARD CURVATURE BY DE- and RE-ATTACHMENT -----

        # Don't apply downward curvature if we are doing upward anticipation
        if v_z is None:

            # Check if we are at the first point of detachment
            if gap < 1e-6 and (ds_bed < -k_crit_down or s_str < -s_crit): 

                # The downward curvature (k_press_down) is based on an estimated attachment length
                # This attachment length depends on geometry (L_attach is roughly equal to 6 * H_hill)
                r_LH = 6.
                s_L = 1 / r_LH

                # An addition to the dettachment length is added based on the distance its takes to bend to leeslope
                L_curv = max((s_str + s_leeside) / k_press_down_base, dx)
                i_curv = int(i + L_curv / dx)

                # Loop forward going downward with s_L per m and track where it intersects the bed
                bed_intersect = []
                for j in range(i_curv, N):
                    zj = z[i] - s_L * (j - i_curv) * dx
                    zj_prev = z[i] - s_L * (j - i_curv - 1) * dx
                    if zj < z_bed[j] and zj_prev >= z_bed[j]:
                        bed_intersect.append(j)
                    
                # Double crossing could be caused due to a small hill on the leeslide
                # If the bed is crossed multiple times, the last one is used
                if bed_intersect:
                    n_intersect = bed_intersect[-1]  
                else:
                    n_intersect = N-1

                # Estimate the height of the feature by taking the max and min bed levels
                h_gap = np.max(z_bed[i:n_intersect+1])- np.min(z_bed[i:n_intersect+1])
                
                # Estimate attachment length and subsequent downward curvature
                L_attach = h_gap * r_LH + L_curv
                k_press_down = 1.5 / L_attach

            else:
                k_press_down = k_press_down_base
            
            # Apply downward curvature
            if ds_bed < -k_crit_down or s_bed < -s_crit or gap > 0.0:

                # An exponential term drives s → -s_leeside and bend steep upward slopes faster
                f_converge = np.exp(1.4*(max(s_str + s_leeside, 0)**0.85)) - 1
                v_z = s_str - k_press_down * f_converge

            # Remains when no anticipation or downward curvature is applied
            else:
                v_z = s_bed_next

        # ----- 3. ADVANCE STREAMLINE -----
        z[i+1] = max(z[i] + v_z * dx, z_bed[i+1])

    return z


def set_sep_params(p, s):
    """
    Set separation–bubble parameters based on perturbation–theory
    scaling (Kroy-type). Values are stored inside the main parameter
    dictionary 'p'. Returns updated p.

    Inputs
    ------
    p : dict
        Global model parameter dictionary.
    s : dict
        Shear fields (contains L, l, z0 via WindShear object).

    Notes
    ------
    Upward curvature (k_press_up) follows Kroy-type perturbation theory.
    Downward curvature is taken as a fraction of the upward limit.
    Other parameters remain simple fixed constants.
    """

    # Characteristic length scales
    f_min = 0.3         # [-] Minimum ratio of tau to tau_0
    # gamma_down = 0.4    # [-] Ratio of k_down to k_up

    # --- Upward curvature limit (perturbation theory) ---
    p['sep_k_press_up'] = (1 - f_min) * (2*np.pi / (np.sqrt(p['L']) * 40)) 

    # # --- Downward curvature limit (softer than upward) ---
    # p['sep_k_press_down'] = gamma_down * p['sep_k_press_up']

    # print(p['sep_k_press_down'])

    # # --- Curvature threshold for detachment trigger ---
    # p['sep_curv_limit_down'] = np.maximum(p['sep_s_crit'] / 0.5, 1.5 * p['sep_k_press_down'])

    return p
