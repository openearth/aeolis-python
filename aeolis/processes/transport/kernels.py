'''Grain speed (Duran 2007) and equilibrium concentration as fused kernels.

The numpy versions in :mod:`aeolis.transport` spread the same pointwise
arithmetic over about a hundred full-grid passes and forty temporaries per
step. Here every cell and fraction is computed in one visit, rows in
parallel. The arithmetic and its order are those of the numpy versions, so the
results agree to round-off (the transcendental functions of numba and numpy
differ in the last bit); the numpy versions stay as the reference in the
tests.

Method codes
------------
grain speed (``GRAINSPEED``): 0 duran_uniform, 1 duran_flat, 2 duran, 3 duran_full
transport (``TRANSPORT``): 0 bagnold, 1 bagnold_gs, 2 kawamura, 3 lettau, 4 dk,
5 sauermann, 6 vanrijn_strypsteen
'''

import numpy as np
from numba import njit, prange

GRAINSPEED = {'duran_uniform': 0, 'duran_flat': 1, 'duran': 2, 'duran_full': 3}
TRANSPORT = {'bagnold': 0, 'bagnold_gs': 1, 'kawamura': 2, 'lettau': 3, 'dk': 4,
             'sauermann': 5, 'vanrijn_strypsteen': 6}


@njit(cache=True, parallel=True, error_model='numpy')
def grainspeed_duran(ustar, ustars, ustarn, ustar0, uth0, uthST, uw, uws, uwn, dzs, dzn,
                     z0, z1, zm_fac, uf, alpha, kappa, method, uw00, uws00, uwn00,
                     ucap, maxiter, tol):
    '''Grain speed per cell and fraction (Duran 2007, eqs 1.60-1.62).

    Parameters
    ----------
    ustar, ustars, ustarn : (ny, nx)
        Shear velocity and components used for the grain speed direction
        (``ustarAir`` for mode 'normal', ``ustar`` for mode 'bed').
    ustar0 : (ny, nx)
        Undisturbed shear velocity.
    uth0 : (ny, nx, nf)
        Grain-size threshold (the ``uth`` of the numpy version).
    uthST : (ny, nx, nf)
        Full threshold, only for the SedTRAILS mask.
    uw, uws, uwn : (ny, nx)
        Wind speed and components (direction where ustar is zero).
    dzs, dzn : (ny, nx)
        Bed slopes along s and n.
    z0, zm_fac, uf, alpha : (nf,)
        Grain roughness d/20, ``c * tv`` (so zm = zm_fac * uth0), settling
        velocity and restitution coefficient per fraction.
    z1 : float
        Reference height.
    method : int
        See ``GRAINSPEED``.
    uw00, uws00, uwn00 : float
        Wind at cell (0, 0), used by ``duran_uniform``.
    ucap : float
        Hard cap on the grain speed magnitude.
    maxiter, tol : int, float
        Newton iteration of ``duran_full``.

    Returns
    -------
    u0, us, un, u, usST, unST : (ny, nx, nf)
    '''
    ny, nx = ustar.shape
    nf = uth0.shape[2]
    u0 = np.zeros((ny, nx, nf))
    us = np.zeros((ny, nx, nf))
    un = np.zeros((ny, nx, nf))
    u = np.zeros((ny, nx, nf))
    usST = np.zeros((ny, nx, nf))
    unST = np.zeros((ny, nx, nf))
    logz = np.empty(nf)
    sq2a = np.empty(nf)
    for f in range(nf):
        logz[f] = np.log(z1 / z0[f])
        sq2a[f] = np.sqrt(2. * alpha[f])

    for j in prange(ny):
        for i in range(nx):
            # direction: from the wind, overruled by the shear where it is non-zero
            ets = 1.0
            etn = 0.0
            if uw[j, i] > 0.:
                ets = uws[j, i] / uw[j, i]
                etn = uwn[j, i] / uw[j, i]
            if ustar[j, i] > 0.:
                ets = ustars[j, i] / ustar[j, i]
                etn = ustarn[j, i] / ustar[j, i]
            for f in range(nf):
                zm = zm_fac[f] * uth0[j, i, f]
                ueff0 = (uth0[j, i, f] / kappa) * (logz[f] + (z1 / zm) * (ustar0[j, i] / uth0[j, i, f] - 1))
                u0f = ueff0 - uf[f] / sq2a[f]
                u0[j, i, f] = u0f
                if method == 0:
                    if uw00 > 0.:
                        usf = u0f * uws00 / uw00
                        unf = u0f * uwn00 / uw00
                        uf_ = u0f
                    else:
                        usf = 0.0
                        unf = 0.0
                        uf_ = 0.0
                elif method == 1:
                    usf = u0f * ets
                    unf = u0f * etn
                    uf_ = u0f
                else:
                    Axs = ets + 2 * alpha[f] * dzs[j, i]
                    Axn = etn + 2 * alpha[f] * dzn[j, i]
                    Ax = np.hypot(Axs, Axn)
                    ueff = (uth0[j, i, f] / kappa) * (logz[f] + (z1 / zm) * (ustar[j, i] / uth0[j, i, f] - 1))
                    usf = (ueff - uf[f] / (sq2a[f] * Ax)) * ets - (sq2a[f] * uf[f] / Ax) * dzs[j, i]
                    unf = (ueff - uf[f] / (sq2a[f] * Ax)) * etn - (sq2a[f] * uf[f] / Ax) * dzn[j, i]
                    uf_ = np.hypot(usf, unf)
                    if method == 3:
                        # Newton-Raphson on eq 1.61 (as solve_duran)
                        u_c = usf + 1j * unf
                        ueff_c = ueff * (ets + 1j * etn)
                        dh_c = dzs[j, i] + 1j * dzn[j, i]
                        for k in range(maxiter):
                            u_mag = max(abs(u_c), 1e-6)
                            rel_u = ueff_c - u_c
                            drag = (rel_u * abs(rel_u)) / uf[f] ** 2
                            fric = u_c / (2 * alpha[f] * u_mag)
                            res = drag - fric - dh_c
                            if abs(res) < tol:
                                break
                            d_drag = -2 * abs(rel_u) / uf[f] ** 2
                            d_fric = -1 / (2 * alpha[f] * u_mag)
                            u_c -= res / (d_drag + d_fric)
                        usf = u_c.real
                        unf = u_c.imag
                        uf_ = abs(u_c)
                # hard cap on the magnitude
                if uf_ > ucap:
                    usf = usf * ucap / uf_
                    unf = unf * ucap / uf_
                    uf_ = ucap
                us[j, i, f] = usf
                un[j, i, f] = unf
                u[j, i, f] = uf_
                # SedTRAILS: no grain speed where the full threshold exceeds the shear velocity
                if uthST[j, i, f] > ustar[j, i]:
                    usST[j, i, f] = usf * 0.
                    unST[j, i, f] = unf * 0.
                else:
                    usST[j, i, f] = usf
                    unST[j, i, f] = unf
    return u0, us, un, u, usST, unST


@njit(cache=True, parallel=True, error_model='numpy')
def equilibrium_concentration(ustar, ustar0, ustar_air, u, uth, uthf, uth0, grain_size,
                              method, coef, rhoa, g, accfac, Cu, Cuf, Cu0, CuAir, CuBed):
    '''Equilibrium concentrations per cell and fraction.

    Writes ``Cu``, ``Cuf``, ``CuAir``, ``CuBed`` (zero where ``ustar`` or ``u``
    is zero) and ``Cu0`` (only where they are not; elsewhere it keeps its
    value), all scaled by ``accfac``, as the numpy version does.

    Parameters
    ----------
    ustar, ustar0, ustar_air : (ny, nx)
        Bed, undisturbed and air shear velocity.
    u : (ny, nx, nf)
        Grain speed magnitude.
    uth, uthf, uth0 : (ny, nx, nf)
        Thresholds.
    grain_size : (nf,)
    method : int
        See ``TRANSPORT``.
    coef : float
        The formula's coefficient (``Cb``, ``Ck``, ``Cl`` or ``Cdk``).
    '''
    ny, nx, nf = u.shape
    Dref = 0.000250
    alpha_sauermann = 0.35
    for j in prange(ny):
        for i in range(nx):
            us_ = ustar[j, i]
            u0_ = ustar0[j, i]
            ua_ = ustar_air[j, i]
            for f in range(nf):
                uf_ = u[j, i, f]
                if us_ != 0. and uf_ != 0.:
                    t = uth[j, i, f]
                    tf = uthf[j, i, f]
                    t0 = uth0[j, i, f]
                    if method == 0:
                        c = coef * rhoa / g
                        v = c * (us_ - t) ** 3 / uf_
                        vf = c * (us_ - tf) ** 3 / uf_
                        v0 = c * (u0_ - t0) ** 3 / uf_
                        va = c * (ua_ - t0) ** 3 / uf_
                    elif method == 1:
                        c = coef * np.sqrt(grain_size[f] / Dref) * rhoa / g
                        v = c * (us_ - t) ** 3 / uf_
                        vf = v
                        v0 = v
                        va = c * (ua_ - t0) ** 3 / uf_
                    elif method == 2:
                        c = coef * rhoa / g
                        v = c * (us_ + t) ** 2 * (us_ - t) / uf_
                        vf = c * (us_ + tf) ** 2 * (us_ - tf) / uf_
                        v0 = c * (u0_ + t0) ** 2 * (u0_ - t0) / uf_
                        va = c * (ua_ + t0) ** 2 * (ua_ - t0) / uf_
                    elif method == 3:
                        c = coef * rhoa / g
                        v = c * us_ ** 2 * (us_ - t) / uf_
                        vf = c * us_ ** 2 * (us_ - tf) / uf_
                        v0 = c * u0_ ** 2 * (u0_ - t0) / uf_
                        va = c * ua_ ** 2 * (ua_ - t0) / uf_
                    elif method == 4:
                        c = coef * rhoa / g
                        v = c * 0.8 * t * (us_ ** 2 - (0.8 * t) ** 2) / uf_
                        vf = c * 0.8 * tf * (us_ ** 2 - (0.8 * tf) ** 2) / uf_
                        v0 = c * 0.8 * t0 * (u0_ ** 2 - (0.8 * t0) ** 2) / uf_
                        va = c * 0.8 * t0 * (ua_ ** 2 - (0.8 * t0) ** 2) / uf_
                    elif method == 5:
                        c = 2. * alpha_sauermann * rhoa / g
                        v = c * (us_ ** 2 - t ** 2)
                        vf = c * (us_ ** 2 - tf ** 2)
                        v0 = c * (u0_ ** 2 - t0 ** 2)
                        va = c * (ua_ ** 2 - t0 ** 2)
                    else:
                        c = coef * rhoa / g
                        v = c * (us_ ** 3 - t ** 3) / uf_
                        vf = v
                        v0 = c * (u0_ ** 3 - t0 ** 3) / uf_
                        va = c * (ua_ ** 3 - t0 ** 3) / uf_
                    v = v if v > 0. else 0.
                    vf = vf if vf > 0. else 0.
                    v0 = v0 if v0 > 0. else 0.
                    va = va if va > 0. else 0.
                    Cu[j, i, f] = v * accfac
                    Cuf[j, i, f] = vf * accfac
                    Cu0[j, i, f] = v0 * accfac
                    CuAir[j, i, f] = va * accfac
                    CuBed[j, i, f] = v * accfac
                else:
                    Cu[j, i, f] = 0.
                    Cuf[j, i, f] = 0.
                    CuAir[j, i, f] = 0.
                    CuBed[j, i, f] = 0.
                    Cu0[j, i, f] = Cu0[j, i, f] * accfac


# ----------------------------------------------------------------------------
# adapters: model state in, model state out (same contracts as aeolis.transport)
# ----------------------------------------------------------------------------

def duran_parameters(p):
    '''Per-fraction constants of the Duran (2007) grain speed (as in :func:`aeolis.transport.duran_grainspeed`).'''
    nf = p['nfractions']
    d = np.asarray(p['grain_size'], dtype=float)
    g = np.repeat(p['g'], nf, axis=0)
    v = np.repeat(p['v'], nf, axis=0)
    srho = p['rhog'] / p['rhoa']
    r = 1.
    c = 14. / (1. + 1.4 * r)
    tv = (v / g ** 2) ** (1 / 3)
    lv = (v ** 2 / (p['Aa'] ** 2 * g * (srho - 1))) ** (1 / 3)
    z0 = d / 20.
    z1 = 35. * lv
    A, B = 0.95, 5.12
    alpha = 0.17 * d / lv
    Sstar = d / (4 * v) * np.sqrt(g * d * (srho - 1.))
    Cd = (4 / 3) * (A + np.sqrt(2 * alpha) * B / Sstar) ** 2
    uf = np.sqrt(4 / (3 * Cd) * (srho - 1.) * g * d)
    return dict(z0=np.ascontiguousarray(z0), z1=float(z1[0]), zm_fac=np.ascontiguousarray(c * tv),
                uf=np.ascontiguousarray(uf), alpha=np.ascontiguousarray(alpha))


def bed_slopes(s, p):
    '''Central-difference bed slopes along s and n (as in the numpy version).'''
    z, x, y = s['zb'], s['x'], s['y']
    dzs = np.zeros(z.shape)
    dzn = np.zeros(z.shape)
    dzs[:, 1:-1] = (z[:, 2:] - z[:, :-2]) / (x[:, 2:] - x[:, :-2])
    dzn[1:-1, :] = (z[:-2, :] - z[2:, :]) / (y[:-2, :] - y[2:, :])
    return dzs, dzn


def duran_grainspeed(s, p, mode='normal'):
    '''Grain speed from the model state; same returns as :func:`aeolis.transport.duran_grainspeed`.'''
    if mode == 'normal':
        ustar, ustars, ustarn = s['ustarAir'], s['ustarsAir'], s['ustarnAir']
    else:
        ustar, ustars, ustarn = s['ustar'], s['ustars'], s['ustarn']
    par = duran_parameters(p)
    dzs, dzn = bed_slopes(s, p)
    return grainspeed_duran(ustar, ustars, ustarn, s['ustar0'], s['uth0'], s['uth'],
                            s['uw'], s['uws'], s['uwn'], dzs, dzn,
                            par['z0'], par['z1'], par['zm_fac'], par['uf'], par['alpha'],
                            float(p['kappa']), GRAINSPEED[p['method_grainspeed']],
                            float(s['uw'][0, 0]), float(s['uws'][0, 0]), float(s['uwn'][0, 0]),
                            5., 100, 1e-2)


def equilibrium(s, p, u):
    '''Fill Cu, Cuf, Cu0, CuAir and CuBed of the state from the grain speed magnitude ``u``.'''
    method = p['method_transport'].lower()
    coef = {'bagnold': 'Cb', 'bagnold_gs': 'Cb', 'kawamura': 'Ck', 'lettau': 'Cl', 'dk': 'Cdk',
            'sauermann': None, 'vanrijn_strypsteen': 'Cb'}[method]
    coef = float(p[coef]) if coef else 0.
    shape = s['uth'].shape
    for k in ('Cu', 'Cuf', 'CuAir', 'CuBed'):
        s[k] = np.zeros(shape)
    s['Cu0'] = np.ascontiguousarray(s['Cu0'], dtype=float)
    equilibrium_concentration(s['ustar'], s['ustar0'], s['ustarAir'], u, s['uth'], s['uthf'], s['uth0'],
                              np.asarray(p['grain_size'], dtype=float), TRANSPORT[method], coef,
                              float(p['rhoa']), float(np.asarray(p['g']).ravel()[0]), float(p['accfac']),
                              s['Cu'], s['Cuf'], s['Cu0'], s['CuAir'], s['CuBed'])
    return s
