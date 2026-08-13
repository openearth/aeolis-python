"""
Vegetation module for dune grasses (new framework).

This module is an alternative to aeolis.vegetation without overwriting it.
It describes the local growth and spreading of dune grasses,
including their effects on shear stress and sediment transport.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
from aeolis import grass_utils as gutils

def initialize(s, p):
    """
    Initialize vegetation state variables.
    Vegetation subgrid is prognostic; main grid is diagnostic.
    """

    f = p['veg_res_factor']

    # --- Make parameters iterable over species and check length -------------
    p = gutils.ensure_grass_parameters(p)

    # --- Convert yearly to secondly rates ------------------------------------
    for param in ['G_h', 'G_c', 'G_s',
                  'dzb_tol_c', 'dzb_tol_s',
                  'dzb_opt_h', 'dzb_opt_c', 'dzb_opt_s']:
        # integer-valued config entries arrive as int arrays; in-place
        # division would fail, so always convert to float first
        p[param] = np.asarray(p[param], dtype=np.float64) / (365.25 * 24.0 * 3600.0)

    # --- Read main-grid vegetation state variables --------------------------
    try:
        s['Nt'][:] = p['Nt_file'].reshape((p['ny']+1, p['nx']+1, p['nspecies']))
        s['hveg'][:] = p['hveg_file'].reshape((p['ny']+1, p['nx']+1, p['nspecies']))
    except Exception as e:
        raise RuntimeError("Shape mismatch for vegetation files.") from e
    
    # --- Initialize vegetation subgrid geometry -----------------------------
    s['x_vsub'], s['y_vsub'] = gutils.generate_grass_subgrid(s['x'], s['y'], f)

    # --- Compute resolution of vegetation subgrid (assumed uniform) ---------
    p['dx_veg'] = np.sqrt((s['x_vsub'][0, 1] - s['x_vsub'][0, 0])**2 
                          + (s['y_vsub'][0, 1] - s['y_vsub'][0, 0])**2)
    dx_main = np.sqrt((s['x'][0, 1] - s['x'][0, 0])**2 
                      + (s['y'][0, 1] - s['y'][0, 0])**2)
    print(f"Grass.py: Vegetation subgrid created with a resolution of {p['dx_veg']:.3f} m (main grid: {dx_main:.3f} m, factor: {f})")

    # --- One-time lift: main grid → vegetation subgrid ----------------------
    s['Nt_vsub']   = gutils.expand_to_subgrid(s['Nt'], f)
    s['hveg_vsub'] = gutils.expand_to_subgrid(s['hveg'], f)

    # --- Build kernel functions for spreading --------------------------------
    s['kernel_c'] = [None] * p['nspecies']
    s['radius_c'] = np.zeros(p['nspecies'], dtype=int)
    for k in range(p['nspecies']):
        s['kernel_c'][k], s['radius_c'][k] = gutils.build_clonal_kernel(
            0.001, p['lmax_c'][k], p['mu_c'][k], p['dx_veg'])
        
    # --- Initial main-grid vegetation metrics -------------------------------
    s['fbend'] = np.ones_like(s['hveg'])
    s['hvegeff'] = s['hveg'].copy()  # initial effective height
    s['lamveg'] = s['Nt'] * s['hvegeff'] * p['d_tiller']
    s['rhoveg'] = s['Nt'] * np.pi * (p['d_tiller'] / 2.0)**2

    return s, p


def update(s, p):

    """
    Main vegetation update.
    All dynamics occur on the vegetation subgrid.
    """

    # --- Time step and resolution factor ------------------------------------
    dt = p['dt_veg']
    f = p['veg_res_factor']

    # --- Burial smoothing (main grid → subgrid, diagnostic → prognostic) ----
    dzb_main = gutils.smooth_burial(s, p)
    s['dzbavg'][:] = dzb_main
    dzb_vsub = gutils.expand_to_subgrid_linear(dzb_main[:,:,None], p['veg_res_factor'])[:,:,0]

    # # --- Expand main-grid state variables to subgrid (for flooding) ---------
    zb_vsub = gutils.expand_to_subgrid_linear(s['zb'][:,:,None], p['veg_res_factor'])[:,:,0]
    TWL_vsub = gutils.expand_to_subgrid_linear(s['TWL'][:,:,None], p['veg_res_factor'])[:,:,0]

    # --- Neighbourhood-averaged densities (all species) ------------------------
    Nt_avg = np.zeros_like(s['Nt_vsub'])
    for l in range(p['nspecies']):
        Nt_avg[:, :, l] = gutils.neighbourhood_average(
            s['Nt_vsub'][:, :, l], p['R_cov'][l], p['dx_veg'])

    # --- Loop over species (subgrid physics) --------------------------------
    for k in range(p['nspecies']):

        Nt   = s['Nt_vsub'][:,:,k]
        hveg = s['hveg_vsub'][:,:,k]

        # --- Burial responses ----------------------------------------------
        B_h = - p['gamma_h'][k] * np.abs(dzb_vsub - p['dzb_opt_h'][k])                          # additive factor
        B_c = np.maximum(1.0 - np.abs(dzb_vsub - p['dzb_opt_c'][k]) / p['dzb_tol_c'][k], 0.0)   # multiplicative factor
        B_s = np.maximum(1.0 - np.abs(dzb_vsub - p['dzb_opt_s'][k]) / p['dzb_tol_s'][k], 0.0)   # multiplicative factor

        # --- Spreading ------------------------------------------------------
        # NOTE: dNt can only be negative when over-saturated (due to competition)
        dNt = spreading(k, Nt, hveg, Nt_avg, B_c, B_s, p, s)                        # [tillers/dt]
        s['Nt_vsub'][:,:,k] = np.maximum(Nt + dNt, 0.0)                             # [tillers]
               
        # --- Local growth ---------------------------------------------------
        dhveg_local = ( p['G_h'][k] * 
                       (1.0 - hveg / p['Hveg'][k])**p['phi_h'][k] + B_h ) * dt      # [m/dt] Local growth rate
        dhveg = dhveg_local - hveg / np.maximum(Nt, 1e-6) * dNt                     # [m/dt] Compensate for density changes
        s['hveg_vsub'][:,:,k] = np.clip(hveg + dhveg, 0.0, p['Hveg'][k])            # [m] Updated height
        s['hveg_vsub'][:,:,k][~(Nt > 0.0)] = 0.0                                    # [m] Zero height where no tillers

        # --- Mortality due to flooding --------------------------------------
        if p['process_tide']:
            hw = np.maximum(TWL_vsub - zb_vsub, 0.0)                                # [m] Local water depth
            hveg = np.maximum(s['hveg_vsub'][:, :, k], 1e-6)                        # [m] Avoid division by zero
            rel_hw = hw / hveg                                                      # [-] Relative water depth
            s['Nt_vsub'][:, :, k] *= np.maximum(0.0, 
                                                1.0 - (rel_hw * dt / p['T_flood'])) # [-] Decay due to flooding

        # --- Mortality due to height reduction ------------------------------
        # Mean tiller height is used; density loss represents implicit 
        # mortality of smaller tillers during burial or erosion.
        hveg = np.maximum(s['hveg_vsub'][:, :, k], 1e-6)
        rel_dhveg = np.clip(dhveg_local / hveg, -1.0, 0.0)                          # [-] Relative height decrease
        mortality_factor = np.clip(1.0 + p['gamma_Nt_decay'][k] * rel_dhveg, 0.0, 1.0) # [-] Mortality factor
        s['Nt_vsub'][:, :, k] *= mortality_factor
        s['Nt_vsub'][:, :, k] = np.maximum(s['Nt_vsub'][:, :, k], 0.0)

        # --- Complete die-off where vegetation height is zero ---------------
        ix_decayed = (s['hveg_vsub'][:, :, k] <= 0.0)
        s['Nt_vsub'][:, :, k][ix_decayed] = 0.0

    # --- Aggregate back to main grid (diagnostic only) ----------------------
    s['Nt']       = gutils.aggregate_from_subgrid(s['Nt_vsub'], f)
    s['hveg']     = gutils.aggregate_from_subgrid(s['hveg_vsub'], f)

    return s


def spreading(k, Nt, hveg, Nt_avg, B_c, B_s, p, s):
    """
    Spatial redistribution of vegetation:
    clonal expansion and seed dispersal.
    """

    # --- Competition-weighted relative densities (Option B) --------------------
    comp = np.zeros_like(Nt)
    for l in range(p['nspecies']):
        comp += p['alpha_comp'][k, l] * (Nt_avg[:, :, l] / p['Nt_max'][l])

    saturation = 1.0 - comp #  np.maximum(1.0 - comp, 0.0)

    maturity = np.clip(hveg / p['Hveg'][k], 0.0, 1.0)

    # --- Tiller production rates --------------------------------------------
    S_c = p['G_c'][k] * Nt * B_c * maturity #* saturation    # [tillers/s] clonal rate 
    S_s = p['G_s'][k] * Nt * maturity # B_s * saturation    # [tillers/s] seed rate CHECK THIS SATURATION (OTHERWISE NO SEED PRODUCTION FROM FULLY COVERED AREAS)
    
    S_c *= p['dt_veg']                                      # [tillers/dt]
    S_s *= p['dt_veg']                                      # [tillers/dt]

    # --- Clonal expansion ---------------------------------------------------
    Nt_clonal_new = np.zeros_like(Nt)
    dNt_clonal = gutils.apply_clonal_kernel(S_c, s['kernel_c'][k])

    dNt_clonal *= saturation  # Limit expansion in saturated areas (after kernel, so on target instead of source)
    ix_pos = dNt_clonal >= 0.0
    Nt_clonal_new[ix_pos] = np.random.poisson(dNt_clonal[ix_pos])
    ix_neg = dNt_clonal < 0.0
    Nt_clonal_new[ix_neg] = -np.random.poisson(-dNt_clonal[ix_neg])

    # --- Seed dispersal -----------------------------------------------------
    Nt_seed_new = gutils.sample_seed_germination(S_s, B_s, p['alpha_s'][k], 
                                                 p['nu_s'][k], p['dx_veg'])

    # --- Sum contributions --------------------------------------------------
    dNt = Nt_clonal_new + Nt_seed_new

    return dNt  


def compute_shear_reduction(s, p):
    """
    Compute vegetation-induced shear reduction.
    """

    # --- Vegetation bending -------------------------------------------------
    for k in range(p['nspecies']):
        s['fbend'][:,:,k] = (p['r_stem'][k] + (1.0 - p['r_stem'][k])
                       * (p['alpha_uw'][k] * s['uw']
                         + p['alpha_Nt'][k] * s['Nt'][:,:,k] 
                         + p['alpha_0'][k]))
        s['fbend'][:,:,k] = np.clip(s['fbend'][:,:,k], 0.0, 1.0)

    # --- Main-grid vegetation metrics ---------------------------------------
    s['hvegeff'] = s['hveg'] * s['fbend']
    s['lamveg'] = s['Nt'] * s['hvegeff'] * p['d_tiller']
    s['rhoveg'] = s['Nt'] * np.pi * (p['d_tiller'] / 2.0)**2

    # --- Weights for normalization ------------------------------------------
    w_sum = np.zeros_like(s['x'])
    w_num_R0 = np.zeros_like(s['x'])
    w_num_L = np.zeros_like(s['x'])
    w_num_h = np.zeros_like(s['x'])
    w_sum_rNt = np.zeros_like(s['x'])

    s['R0veg'] = np.ones_like(s['x'])
    L_decay = np.zeros_like(s['x'])
    hvegeff = np.zeros_like(s['x'])
    rNt = np.zeros_like(s['x'])

    # --- Wind convention ID (Numba) -----------------------------------------
    if p['wind_convention'] == 'nautical':
        udir_id = 0
    elif p['wind_convention'] == 'cartesian':
        udir_id = 1
    else:
        raise ValueError(f"Unknown wind_convention: {p['wind_convention']}")

    # --- Compute per-species local and leeside reduction --------------------
    for k in range(p['nspecies']):

        # Weighting function based on maturity and density
        maturity = s['hveg'][:, :, k] / p['Hveg'][k]
        density  = s['Nt'][:, :, k] / p['Nt_max'][k]
        w = maturity * density

        # Local Raupach reduction per species 
        R0_k = 1.0 / np.sqrt(1.0 + p['m_veg'][k] * p['beta_veg'][k] * s['lamveg'][:, :, k])

        # Accumulate weighted numerators for later normalization
        w_sum   += w
        w_num_R0 += w * R0_k
        w_num_L += w * (s['hvegeff'][:, :, k] / p['c1_okin'][k])
        w_num_h += w * s['hvegeff'][:, :, k]
        w_sum_rNt += w * (s['Nt'][:, :, k] / p['Nt_max'][k])

    # --- Final normalization (separate loop / block) ------------------------
    ix = w_sum != 0.0
    s['R0veg'][ix] = w_num_R0[ix] / w_sum[ix]
    L_decay[ix] = w_num_L[ix] / w_sum[ix]
    hvegeff[ix] = w_num_h[ix] / w_sum[ix]
    rNt[ix] = w_sum_rNt[ix] / w_sum[ix]

    # --- Compute leeside Okin reduction -------------------------------------
    if p['process_vegetation_leeside']:
        R_okin, hvegeff_zeta, rNt_zeta, Rveg_zeta = gutils.compute_okin_reduction(
            s['x'], s['y'], s['R0veg'], s['udir'], L_decay, hvegeff, rNt, udir_id)
        s['Rveg'] = np.minimum(s['R0veg'], R_okin)

        # Store representative vegetation metrics for zeta adjustment in vegetation wake
        s['hvegeff_zeta'] = np.maximum(hvegeff_zeta, hvegeff)
        s['rNt_zeta'] = np.maximum(rNt_zeta, rNt)

        s['Rveg_zeta'] = Rveg_zeta
        s['Rveg_zeta'][(s['Rveg'] == s['R0veg'])] = 1. 

    else:
        s['Rveg'] = s['R0veg'].copy()

    # --- Apply Gaussian filter to Rveg --------------------------------------
    if p['veg_sigma'] > 0.0:
        s['Rveg'] = ndimage.gaussian_filter(s['Rveg'], sigma=p['veg_sigma'])

    return s


def apply_shear_reduction(s, p):
    """
    Apply vegetation-induced shear reduction to wind shear.
    """

    ets = np.zeros(s['zb'].shape)
    etn = np.zeros(s['zb'].shape)

    ix = s['ustar'] != 0

    ets[ix] = s['ustars'][ix] / s['ustar'][ix]
    etn[ix] = s['ustarn'][ix] / s['ustar'][ix]

    s['ustar'] *= s['Rveg']
    s['ustars'] = s['ustar'] * ets
    s['ustarn'] = s['ustar'] * etn

    return s
