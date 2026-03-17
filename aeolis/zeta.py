"""
Module for computing bed–interaction factor zeta.
"""

import numpy as np
import logging
import matplotlib.pyplot as plt
from scipy.special import gamma
from scipy import ndimage

# initialize logger
import logging
logger = logging.getLogger(__name__)

def initialize(s, p):
    """
    Initialize zeta-related parameters.
    """

    s['zeta'][:] = p['zeta_base']   # Bed–interaction factor
    s['kzeta'][:] = 1.0    # Shape parameter for Weibull function (1 = exponential)
    s['Lzeta'][:] = 1.0  # Scale parameter for Weibull function

    return s


def compute_zeta(s, p):
    """
    Compute bed–interaction factor zeta. Order matters.
    """

    s['zeta'][:] = p['zeta_base']

    if p['th_nelayer']: 
        s = zeta_from_non_erodible(s, p)

    if p['process_moist']:
        s = zeta_from_moisture(s, p)

    if p['process_separation']:
        s = zeta_from_separation(s, p)

    if p['zeta_sheltering']:
        if p['th_sheltering']:
            s = zeta_from_sheltering(s, p)
        else:
            logger.warning("Zeta sheltering requested but sheltering threshold not set.")

    if p['process_vegetation']: # and p['method_vegetation'] == 'grass':
        s = zeta_from_vegetation(s, p)

    # --- Apply Gaussian filter to zeta --------------------------------------
    s['zeta'] = np.clip(s['zeta'], 0.004, 1.0)

    return s


def zeta_from_non_erodible(s, p):
    """Apply non-erodible layer effect on zeta.
    When the non-erodible layer is exposed (or close to the surface),
    bed interaction is assumed to be air-dominated and zeta is set to
    zero locally.
    """

    # Influence of non-erodible layer on bed interaction parameter zeta
    thuthlyr = 0.05
    ix = (s['zb'] <= s['zne'] + thuthlyr)

    # Air-dominated interaction when non-erodible layer is exposed
    s['zeta'][ix] = 0.  

    return s


def zeta_from_moisture(s, p):
    """Apply moisture effect on zeta.
    Moisture increases the fraction of transport that interacts with the bed.
    This increases zeta towards 1 for wet beds, while keeping the current
    baseline zeta for dry conditions.
    """

    # Parameters
    w = s['moist']
    zeta = s['zeta']
    p_zeta_moist = p['p_zeta_moist']

    # Compute zeta modification due to moisture
    s['zeta'] = zeta + (1.0 - zeta) * (w / 0.4)**p_zeta_moist

    return s

def zeta_from_separation(s, p):
    """Apply separation bubble effect on zeta.
    The separation bubble elevates the transport layer, reducing the
    fraction of sediment transport that interacts with the bed.
    This reduces zeta towards 0 for large separation bubbles, while
    keeping the current baseline zeta for no separation.
    """

    tau_sep, slope = 0.5, 0.2
    delta = 1. / (slope * tau_sep)

    zsepdelta = np.maximum(np.minimum(delta * s['hsep'], 1.), 0.)
    s['zeta'] = zsepdelta + (1. - zsepdelta) * s['zeta']

    return s


def zeta_from_sheltering(s, p):
   
    Rti = s['Rti'].reshape(s['zeta'].shape)
    s['zeta'] *= (1.0 / Rti)

    return s


def zeta_from_vegetation(s, p):
    """
    Compute vegetation-induced bed–interaction factor zeta.

    Vegetation modifies the fraction of sediment transport that interacts
    with the bed by elevating the transport layer. Only the part of the
    vertical transport distribution below the effective vegetation height
    contributes to bed-dominated transport; sediment above behaves as
    airborne transport.

    The vertical transport distribution is represented by a Weibull-type
    formulation with characteristic decay length L_d and shape parameter k.
    The computation is performed per vegetation species and subsequently
    normalised using species weights based on vegetation maturity and density.
    """

    # --- Reference decay length for flat-bed (L_d) --------------------------
    lambda_salt = 1.
    LD_FAC = 2. # FIX ????????
    Ld = LD_FAC * s['ustarAir']**2 / (p['g'] * lambda_salt)  # [m] 
    Ld = np.maximum(Ld, 1e-12)
    s['Lzeta'] = Ld # Store Lzeta for output

    # --- Empirical Weibull parameters ---------------------------------------
    a = p['a_weibull']
    b = p['b_weibull']

    # --- Defaults for cells without vegetation ------------------------------
    zeta_out = s['zeta'].copy()
    k_out = s['kzeta'].copy()
    L_out = s['Lzeta'].copy()

    # --- Per-species computation --------------------------------------------
    for ksp in range(p['nspecies']):

        # --- Species weighting based on maturity and density ----------------
        if p['method_vegetation'] == 'grass':
            hveg_eff = s['hvegeff'][:, :, ksp]  # Effective vegetation height for zeta adjustment
            hveg_eff_zeta = s['hvegeff_zeta'] # Effective vegetation height for zeta adjustment
        elif p['method_vegetation'] == 'duran':
            hveg_eff = s['hveg']
            hveg_eff_zeta = s['hveg']
        else:
            ValueError(f"Unknown vegetation method: {p['method_vegetation']}")

        # --- Upward lift of transport layer ---------------------------------
        h_sep = np.maximum(s['zsep'] - s['zb'], 0.0)            # [m]
        h_veg_lift = hveg_eff * p['alpha_lift'][ksp]            # [m] OR ACTUAL HVEG?
        h_lift = np.maximum(h_sep, h_veg_lift) + Ld             # [m]

        # --- Weibull shape parameter k (1 → 3) ------------------------------
        h_ratio = np.maximum(h_lift / Ld - 1.0, 0.0)
        k_zeta = 1.0 + 2.0 * (1.0 - np.exp(-a * h_ratio**b))    # [-]
        k_zeta = np.maximum(k_zeta, 1e-12)

        # --- Bed–interaction factor from Weibull CDF ------------------------
        Ld_eff = h_lift / gamma(1.0 + 1.0 / k_zeta)
        zeta_k = 1.0 - np.exp(-(hveg_eff_zeta / Ld_eff)**k_zeta)
        
        # --- Reduction due to bouncing (skimming) ---------------------------
        zeta_k *= (1.0 - p['bounce'][ksp])

        # --- Compensate for vegetation density ------------------------------
        # Prevents excessive (and spikey) reduction of zeta at low vegetation densities
        pNt_zeta = p['pNt_zeta'][ksp] # Exponent for density effect on zeta (lower values = less effect)
        if p['method_vegetation'] == 'grass':
            density_factor = s['rNt_zeta'] ** pNt_zeta
        elif p['method_vegetation'] == 'duran':
            density_factor = s['rhoveg'] ** pNt_zeta
        zeta_k = 1.0 - density_factor * (1.0 - zeta_k)

        # Adjust in vegetation wake (sheltering effect)
        zeta_k = 1.0 - s['Rveg_zeta'] * (1.0 - zeta_k)

        # --- Take the minimum zeta across species ---------------------------
        ix_zeta = (zeta_k <= zeta_out) 
        zeta_out[ix_zeta] = zeta_k[ix_zeta]
        k_out[ix_zeta] = k_zeta[ix_zeta]
        L_out[ix_zeta] = Ld_eff[ix_zeta]

    # --- Store results ------------------------------------------------------
    s['zeta'] = np.clip(zeta_out, 0.0, 1.0)
    s['kzeta'] = k_out
    s['Lzeta'] = L_out

    return s