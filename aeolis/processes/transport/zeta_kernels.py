'''Bed-interaction factor from vegetation as one fused kernel.

The numpy version (:func:`aeolis.zeta.zeta_from_vegetation`) takes about 25
full-grid transcendental passes per species. Here each cell is computed in
one visit, rows in parallel, with the same arithmetic in the same order, so
the results agree to round-off.
'''

import math
import numpy as np
from numba import njit, prange


@njit(cache=True, parallel=True, error_model='numpy')
def zeta_vegetation(zeta, kzeta, ustar_air, zsep, zb, hveg_eff, hveg_eff_zeta, density,
                    Rveg_zeta, alpha_lift, bounce, pNt_zeta, a, b, g):
    '''Minimum over species of the Weibull bed-interaction factor.

    Parameters
    ----------
    zeta, kzeta : (ny, nx)
        Current factor and shape parameter (the defaults for cells without vegetation).
    ustar_air : (ny, nx)
    zsep, zb : (ny, nx)
        Separation streamline and bed level.
    hveg_eff : (ny, nx, nspecies)
        Effective vegetation height per species.
    hveg_eff_zeta, density, Rveg_zeta : (ny, nx)
        Effective height for the factor, density measure (``rNt_zeta`` or
        ``rhoveg``) and wake reduction.
    alpha_lift, bounce, pNt_zeta : (nspecies,)
    a, b, g : float
        Weibull parameters and gravity.

    Returns
    -------
    zeta_out : (ny, nx)
        Clipped to [0, 1].
    k_out, L_out : (ny, nx)
    '''
    ny, nx, nsp = hveg_eff.shape
    zeta_out = np.empty((ny, nx))
    k_out = np.empty((ny, nx))
    L_out = np.empty((ny, nx))
    lambda_salt = 1.
    LD_FAC = 2.
    for j in prange(ny):
        for i in range(nx):
            Ld = LD_FAC * ustar_air[j, i] ** 2 / (g * lambda_salt)
            Ld = max(Ld, 1e-12)
            z_o = zeta[j, i]
            k_o = kzeta[j, i]
            L_o = Ld
            h_sep = max(zsep[j, i] - zb[j, i], 0.0)
            for k in range(nsp):
                h_veg_lift = hveg_eff[j, i, k] * alpha_lift[k]
                h_lift = max(h_sep, h_veg_lift) + Ld
                h_ratio = max(h_lift / Ld - 1.0, 0.0)
                k_zeta = 1.0 + 2.0 * (1.0 - np.exp(-a * h_ratio ** b))
                k_zeta = max(k_zeta, 1e-12)
                Ld_eff = h_lift / math.gamma(1.0 + 1.0 / k_zeta)
                zeta_k = 1.0 - np.exp(-(hveg_eff_zeta[j, i] / Ld_eff) ** k_zeta)
                zeta_k *= (1.0 - bounce[k])
                density_factor = density[j, i] ** pNt_zeta[k]
                zeta_k = 1.0 - density_factor * (1.0 - zeta_k)
                zeta_k = 1.0 - Rveg_zeta[j, i] * (1.0 - zeta_k)
                if zeta_k <= z_o:
                    z_o = zeta_k
                    k_o = k_zeta
                    L_o = Ld_eff
            zeta_out[j, i] = min(max(z_o, 0.0), 1.0)
            k_out[j, i] = k_o
            L_out[j, i] = L_o
    return zeta_out, k_out, L_out


def zeta_from_vegetation(s, p):
    '''State adapter with the contract of :func:`aeolis.zeta.zeta_from_vegetation`.'''
    nsp = p['nspecies']
    if p['method_vegetation'] == 'grass':
        hveg_eff = np.ascontiguousarray(s['hvegeff'][:, :, :nsp], dtype=float)
        hveg_eff_zeta = np.ascontiguousarray(s['hvegeff_zeta'], dtype=float)
        density = np.ascontiguousarray(s['rNt_zeta'], dtype=float)
    elif p['method_vegetation'] == 'duran':
        hveg_eff = np.ascontiguousarray(np.repeat(s['hveg'][:, :, np.newaxis], nsp, axis=2), dtype=float)
        hveg_eff_zeta = np.ascontiguousarray(s['hveg'], dtype=float)
        density = np.ascontiguousarray(s['rhoveg'], dtype=float)
    else:
        raise ValueError("Unknown vegetation method: %s" % p['method_vegetation'])
    par = [np.ascontiguousarray(np.asarray(p[k], dtype=float).ravel()[:nsp]) for k in ('alpha_lift', 'bounce', 'pNt_zeta')]
    zeta, kzeta, Lzeta = zeta_vegetation(
        np.ascontiguousarray(s['zeta'], dtype=float), np.ascontiguousarray(s['kzeta'], dtype=float),
        np.ascontiguousarray(s['ustarAir'], dtype=float), np.ascontiguousarray(s['zsep'], dtype=float),
        np.ascontiguousarray(s['zb'], dtype=float), hveg_eff, hveg_eff_zeta, density,
        np.ascontiguousarray(s['Rveg_zeta'], dtype=float), *par,
        float(p['a_weibull']), float(p['b_weibull']), float(np.asarray(p['g']).ravel()[0]))
    s['zeta'] = zeta
    s['kzeta'] = kzeta
    s['Lzeta'] = Lzeta
    return s
