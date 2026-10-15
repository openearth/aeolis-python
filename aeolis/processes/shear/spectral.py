"""Linear shear-stress perturbation computed directly on the model grid.

The surface shear stress over low topography is

    tau = tau0 * (1 + dtau),   dtau_hat(k) = T(k) * h_hat(k)

with T the Weng et al. (1991) transfer function (as used by
:class:`aeolis.shear.WindShear`). The operator commutes with rotation, so
instead of building a computational grid aligned with the wind, the
wavenumbers of the *model grid* are rotated:

    k_par  =  kx cos(theta) + ky sin(theta)     (along the wind)
    k_perp = -kx sin(theta) + ky cos(theta)     (across the wind)

No secondary grid and no interpolation are needed. The convolution is made
linear (non-periodic) by padding the bed to about twice its size. The pad
continues the edge values of the bed (relative to its boundary mean) and is
blended smoothly to that mean in both directions, i.e. the topography is
continued flat outside the domain. Compared with the same operator on a
domain extended 3-6x beyond every open edge, this reproduces the stress on
the embryo-dune, foredune, blowout and Boschplaat beds to within 0.4-2 % of
the perturbation rms (see docs/overhaul/PROGRESS.md, "Why the shear fields
differ"). The pad uses the edge cells only, so a feature that does not reach
the edges cannot change how the edges are treated. A first version removed a
best-fit plane of the whole bed before padding; that deleted the real
speed-up over a finite ramp such as a beach rising to a dune and was removed
(2026-10-15). See docs/overhaul/SHEAR_ON_MAIN_GRID.md.

Three transfer functions are available (``transfer``): ``weng``, the
original form with one length scale L for all modes; ``weng_mode``, the same
Weng et al. (1991) coefficients with L/4 = lambda/4 for every mode; and
``fca2010``, the fit of Fourriere, Claudin & Andreotti (2010) for a rough bed.
The last two depend on k z0 only, so the configured L has no effect.

Everything here works on plain numpy arrays and can be used outside AeoLiS:

>>> sh = SpectralShear(ny=200, nx=300, ds=2., dn=2., L=100., l=1., z0=0.001)
>>> dtau_par, dtau_perp = sh.perturbation(zb, theta=np.deg2rad(20.))
"""

import numpy as np
import scipy.fft
from scipy.special import kve, lambertw

KAPPA = 0.41

#: transfer functions SpectralShear accepts
TRANSFERS = ('weng', 'weng_mode', 'fca2010')
SHORT_WAVES = ('none', 'cell_avg', 'floor')

#: Fourriere, Claudin & Andreotti (2010): coefficients of the empirical fit of
#: A(k z0) and B(k z0) for the rough mixing-length closure (arXiv:0708.2823v3,
#: section "Geometrically induced roughness"; not in the printed JFM paper)
FCA_A = (1.0702, 0.093069, 0.10838, 0.024835, 0.041603, 0.0010625)
FCA_B = (0.036989, 0.15765, 0.11518, 0.0020249, 0.0028725, 0.00053483)


def boundary_layer_scales(L, l, z0):
    '''Inner-layer height, middle-layer height and U(l)/U(hm) as in WindShear.

    ``L`` is the configured hill length scale; like the original
    implementation, L/4 is used (see the open question in the design note).
    '''
    L4 = L / 4.
    for _ in range(5):
        l = 2 * KAPPA ** 2 * L4 / np.log(l / z0)
    hm = 1.0
    for _ in range(5):
        hm = L4 / np.sqrt(np.log(hm / z0))
    ul = np.log(l / z0) / np.log(hm / z0)
    return L4, l, hm, ul


def weng_transfer(kpar, kperp, L, l, z0, kfloor=0.):
    '''Weng et al. (1991) transfer functions for along- and cross-wind stress.

    Parameters
    ----------
    kpar, kperp : ndarray
        Along-wind and cross-wind wavenumbers [rad/m] (any shape).
    L, l, z0 : float
        Hill length scale, inner-layer height (initial guess) and roughness [m].

    Returns
    -------
    T_par, T_perp : complex ndarray
        dtau_par_hat / h_hat and dtau_perp_hat / h_hat.

    Notes
    -----
    The along-wind function is regrouped so that it has no 1/k_par factor:

        T_par = 2/ul^2 * ( -k_par^2/k + (2 ln(l/z0) k_par^2/k + k) * R(sigma) )

    with R = sigma K1(2 sigma) / K0(2 sigma), sigma = sqrt(i L k_par z0 / l).
    R -> 0 for k_par -> 0, but only logarithmically: modes almost parallel to
    the wind keep a finite response down to extremely small k_par. Below
    ``kfloor`` (e.g. the fundamental wavenumber of the domain, i.e. modes
    longer than the domain along the wind) R is therefore evaluated at
    ``kfloor`` and ramped linearly to 0, which keeps T continuous in the
    wind direction. The ratio uses exponentially scaled Bessel functions.
    '''
    L4, l, hm, ul = boundary_layer_scales(L, l, z0)
    kpar = np.asarray(kpar, dtype=float)
    kperp = np.asarray(kperp, dtype=float)
    k = np.hypot(kpar, kperp)

    T_par = np.zeros(np.broadcast(kpar, kperp).shape, dtype=complex)
    T_perp = np.zeros_like(T_par)
    ok = (kpar != 0.) & (k != 0.)
    if not np.any(ok):
        return T_par, T_perp

    kp, kk, kq = np.broadcast_to(kpar, ok.shape)[ok], k[ok], np.broadcast_to(kperp, ok.shape)[ok]
    kr = np.where(np.abs(kp) < kfloor, np.copysign(kfloor, kp), kp)
    ramp = np.minimum(np.abs(kp) / kfloor, 1.) if kfloor > 0 else 1.
    sigma = np.sqrt(1j * L4 * kr * z0 / l)
    R = sigma * kve(1, 2. * sigma) / kve(0, 2. * sigma) * ramp
    c = 2. / ul ** 2
    T_par[ok] = c * (-kp ** 2 / kk + (2. * np.log(l / z0) * kp ** 2 / kk + kk) * R)
    s2 = 2. * np.sqrt(2.) * sigma
    # sigma K1(2 sqrt2 sigma) = sigma * kve(1, z) * exp(-z)
    T_perp[ok] = c * (kp * kq / kk) * s2 * kve(1, s2) * np.exp(-s2)
    return T_par, T_perp


def mode_coefficients(kz0, transfer):
    '''A + iB of a bed wave along the wind as a function of k z0.

    For a single mode along the wind (k > 0) the stress perturbation is
    dtau_hat = (A + iB) k h_hat: A is the response in phase with the bed (the
    stress is highest on the crest), B the upwind shift of the maximum.

    ``weng_mode``: the Weng et al. (1991) function of :func:`weng_transfer`
    with L/4 = lambda/4 = pi/(2k) for every mode (Charru et al. 2013, eq. 9),
    with the inner- and middle-layer heights solved exactly.
    ``fca2010``: the fit of Fourriere, Claudin & Andreotti (2010),
    A = 2 + P_a(R)/Q_a(R), B = P_b(R)/Q_b(R), R = ln(2 pi / k z0); A -> 2 and
    B -> 0 for very long waves (potential flow), but only as 1/R.
    '''
    kz0 = np.asarray(kz0, dtype=float)
    if transfer == 'fca2010':
        R = np.log(2. * np.pi / kz0)
        a, b = FCA_A, FCA_B
        A = 2. + (a[0] + a[1] * R + a[2] * R ** 2 + a[3] * R ** 3) / (1. + a[4] * R ** 2 + a[5] * R ** 4)
        B = (b[0] + b[1] * R + b[2] * R ** 2 + b[3] * R ** 3) / (1. + b[4] * R ** 2 + b[5] * R ** 4)
        return A + 1j * B
    if transfer == 'weng_mode':
        L4 = np.pi / (2. * kz0)                                 # lambda/4 in units of z0
        lnl = lambertw(2. * KAPPA ** 2 * L4).real               # l ln(l/z0) = 2 kappa^2 lambda/4
        lnhm = 0.5 * lambertw(2. * L4 ** 2).real                # hm ln^1/2(hm/z0) = lambda/4
        sigma = np.sqrt(1j * L4 * kz0 * np.exp(-lnl))           # sqrt(i (lambda/4) k z0 / l)
        R = sigma * kve(1, 2. * sigma) / kve(0, 2. * sigma)
        return 2. * (lnhm / lnl) ** 2 * (-1. + (1. + 2. * lnl) * R)
    raise ValueError('mode_coefficients: unknown transfer %r' % (transfer,))


class SpectralShear:
    '''Shear perturbation on a uniform rectilinear grid (no secondary grid).

    Parameters
    ----------
    ny, nx : int
        Number of grid points in n (rows) and s (columns).
    ds, dn : float
        Grid spacing in s and n [m].
    L, l, z0 : float
        Hill length scale, inner-layer height and roughness length [m].
    pad_factor : float
        Padded FFT size relative to the grid (>= 2 for a non-circular result;
        3 halves the remaining boundary error at 2.25x the FFT cost).
    blend_length : float, optional
        Distance over which the padded bed relaxes to the reference level
        [m]; default ``max(L/2, 10 grid cells)``, limited by the pad.
    filter_frac : float
        Isotropic low-pass cut-off as a fraction of the Nyquist wavenumber of
        the coarser grid direction; 0 or ``None`` (the default) disables it.
        The filter removes the whole stress response of the shortest
        resolved waves, including the in-phase part (stress highest on the
        crest) that erodes the crest and flattens a ripple of two to three
        cells. With 0.6 such ripples neither grow nor decay and grid-scale
        bed noise accumulates into stripes; without it the model flattens
        them as on a fine grid (docs/overhaul/SHEAR_ON_MAIN_GRID.md).
    short_waves : str
        Treatment of the shortest resolved waves, on top of ``filter_frac``.
        ``none`` (default): full response. ``cell_avg``: the stress averaged
        over a grid cell, F = sinc(kx ds/2) sinc(ky dn/2), which is how the
        grid represents the bed itself (no parameter; keeps 64 % of the
        response at two cells along an axis, 41 % on the diagonal).
        ``floor``: F = f + (1 - f) exp(-(k/kc)^8) with kc at four cells of the
        coarser direction and f = ``short_floor``. Both lower the cell-scale u*
        noise from a rough bed but keep part of the in-phase response that
        flattens two-cell ripples, unlike ``filter_frac``, which removes it.
    short_floor : float
        Fraction f of the response kept for waves shorter than four cells
        with ``short_waves = floor``.
    transfer : str
        ``weng`` (default): Weng et al. (1991) with the configured L for all
        modes. ``weng_mode`` and ``fca2010``: A(k z0) + iB(k z0) per mode
        (see :func:`mode_coefficients`), applied in the usual
        two-dimensional form T_par = (A + iB) k_par^2/k and
        T_perp = (A + iB) k_par k_perp/k, with B changing sign with k_par.
        Here L and l have no effect.
    periodic : (bool, bool)
        Treat the n (rows) and/or s (columns) direction as periodic, e.g. for
        circular boundary conditions. A periodic direction is neither padded
        nor tapered, so the FFT gives the exact solution for a
        topography that repeats itself. For a semi-1D transect with circular
        lateral boundaries this is the solution for an infinitely long ridge.
    '''

    def __init__(self, ny, nx, ds, dn, L, l, z0, pad_factor=2., blend_length=None,
                 filter_frac=0., ntable=4096, periodic=(False, False), transfer='weng',
                 short_waves='none', short_floor=0.5):
        if short_waves not in SHORT_WAVES:
            raise ValueError('SpectralShear: unknown short_waves %r, use one of %s'
                             % (short_waves, ', '.join(SHORT_WAVES)))
        if transfer not in TRANSFERS:
            raise ValueError('SpectralShear: unknown transfer %r, use one of %s' % (transfer, ', '.join(TRANSFERS)))
        self.kind = transfer
        self.ny, self.nx, self.ds, self.dn = int(ny), int(nx), float(ds), float(dn)
        self.L, self.l, self.z0 = float(L), float(l), float(z0)
        self.periodic = (bool(periodic[0]), bool(periodic[1]))

        def padded(n, per):
            if per:
                return n
            return scipy.fft.next_fast_len(max(int(np.ceil(pad_factor * n)), n + 2), real=True)
        self.NY = padded(ny, self.periodic[0])
        self.NX = padded(nx, self.periodic[1])
        self.oy = (self.NY - ny) // 2          # offset of the model grid inside the pad
        self.ox = (self.NX - nx) // 2

        # separable taper: 1 on the model grid, extrusion zone, cosine blend to 0
        if blend_length is None:
            blend_length = max(L / 2., 10. * max(ds, dn))
        self.wy = np.ones(ny) if self.periodic[0] else self._taper(self.NY, self.oy, ny, dn, blend_length)
        self.wx = np.ones(nx) if self.periodic[1] else self._taper(self.NX, self.ox, nx, ds, blend_length)

        # wavenumbers of the padded grid (rfft along s)
        self.kx = 2. * np.pi * scipy.fft.rfftfreq(self.NX, ds)[np.newaxis, :]
        self.ky = 2. * np.pi * scipy.fft.fftfreq(self.NY, dn)[:, np.newaxis]
        k = np.hypot(self.kx, self.ky)
        if filter_frac:
            kc = filter_frac * np.pi / max(ds, dn)
            self.lowpass = np.exp(-(k / kc) ** 8)
        else:
            self.lowpass = np.ones_like(k)
        if short_waves == 'cell_avg':
            self.lowpass = self.lowpass * (np.sinc(self.kx * ds / (2. * np.pi)) * np.sinc(self.ky * dn / (2. * np.pi)))
        elif short_waves == 'floor':
            kc = 0.5 * np.pi / max(ds, dn)
            self.lowpass = self.lowpass * (short_floor + (1. - short_floor) * np.exp(-(k / kc) ** 8))

        # 1D tables of the k_par-dependent Bessel factors (exact symmetry for k_par < 0);
        # below the fundamental wavenumber of the padded domain R is ramped to 0
        self.kfloor = 2. * np.pi / max(self.NX * ds, self.NY * dn)
        kmax = np.sqrt(2.) * np.pi / min(ds, dn) * 1.01
        self._tk = np.geomspace(self.kfloor, kmax, ntable)
        L4, l_, hm, ul = boundary_layer_scales(L, l, z0)
        sigma = np.sqrt(1j * L4 * self._tk * z0 / l_)
        self._R = sigma * kve(1, 2. * sigma) / kve(0, 2. * sigma)
        s2 = 2. * np.sqrt(2.) * sigma
        self._S = s2 * kve(1, s2) * np.exp(-s2)
        self._c = 2. / ul ** 2
        self._loglz0 = 2. * np.log(l_ / z0)
        self._ltk = np.log(self._tk)
        if transfer != 'weng':
            self._AB = mode_coefficients(self._tk * z0, transfer)
        self._cache = {}

    @staticmethod
    def _taper(N, offset, n, d, blend_length):
        w = np.zeros(N)
        pad_lo, pad_hi = offset, N - offset - n
        nb = int(min(round(blend_length / d), min(pad_lo, pad_hi)))
        ne = max(min(pad_lo, pad_hi) - nb, 0) // 2        # short extrusion before the blend
        w[offset - ne:offset + n + ne] = 1.
        ramp = 0.5 * (1. + np.cos(np.linspace(0., np.pi, nb + 2)[1:-1]))
        if nb > 0:
            w[offset - ne - nb:offset - ne] = ramp[::-1]
            w[offset + n + ne:offset + n + ne + nb] = ramp
        return w

    def _table(self, kpar, table, ramp=False):
        a = np.abs(kpar)
        la = np.log(np.maximum(a, self._tk[0]))
        re = np.interp(la, self._ltk, table.real)
        im = np.interp(la, self._ltk, table.imag)
        v = re + 1j * np.where(kpar < 0., -im, im)      # T(-k) = conj T(k)
        if ramp:
            v = v * np.minimum(a / self.kfloor, 1.)
        return v

    def transfer(self, theta):
        '''Transfer functions (T_par, T_perp) on the padded rfft grid for wind angle theta.'''
        key = round(float(theta), 10)
        if key in self._cache:
            return self._cache[key]
        c, s = np.cos(theta), np.sin(theta)
        kpar = self.kx * c + self.ky * s
        kperp = -self.kx * s + self.ky * c
        k = np.hypot(self.kx, self.ky)
        ksafe = np.where(k == 0., 1., k)
        if self.kind == 'weng':
            R = self._table(kpar, self._R, ramp=True)
            S = self._table(kpar, self._S)
            T_par = self._c * (-kpar ** 2 / ksafe + (self._loglz0 * kpar ** 2 / ksafe + k) * R)
            T_perp = self._c * (kpar * kperp / ksafe) * S
        else:
            lk = np.log(np.maximum(k, self._tk[0]))
            AB = (np.interp(lk, self._ltk, self._AB.real)
                  + 1j * np.sign(kpar) * np.interp(lk, self._ltk, self._AB.imag))   # T(-k) = conj T(k)
            T_par = AB * kpar ** 2 / ksafe
            T_perp = AB * kpar * kperp / ksafe
        T_par[k == 0.] = 0.
        T_perp[k == 0.] = 0.
        out = (T_par * self.lowpass, T_perp * self.lowpass)
        if len(self._cache) > 64:
            self._cache.pop(next(iter(self._cache)))
        self._cache[key] = out
        return out

    def pad(self, zb):
        '''Bed on the padded grid: edge extrusion blended to the boundary mean.'''
        zb = np.asarray(zb, dtype=float)
        edges = []
        if not self.periodic[0]:
            edges += [zb[0, :], zb[-1, :]]
        if not self.periodic[1]:
            edges += [zb[:, 0], zb[:, -1]]
        if not edges:
            return zb
        zref = np.concatenate(edges).mean()
        h = np.pad(zb - zref, ((self.oy, self.NY - self.ny - self.oy),
                               (self.ox, self.NX - self.nx - self.ox)), mode='edge')
        return h * self.wy[:, np.newaxis] * self.wx[np.newaxis, :]

    def perturbation(self, zb, theta, workers=-1):
        '''Relative stress perturbation along and across the wind on the model grid.

        Parameters
        ----------
        zb : ndarray (ny, nx)
            Bed level (including any separation envelope) [m].
        theta : float
            Direction the wind blows towards, counter-clockwise from the s
            axis of the grid [rad].

        Returns
        -------
        dtau_par, dtau_perp : ndarray (ny, nx)
        '''
        H = scipy.fft.rfft2(self.pad(zb), workers=workers)
        T_par, T_perp = self.transfer(theta)
        sl = (slice(self.oy, self.oy + self.ny), slice(self.ox, self.ox + self.nx))
        dpar = scipy.fft.irfft2(T_par * H, s=(self.NY, self.NX), workers=workers)[sl]
        dperp = scipy.fft.irfft2(T_perp * H, s=(self.NY, self.NX), workers=workers)[sl]
        return dpar, dperp

    def stress(self, zb, taus0, taun0, workers=-1):
        '''Perturbed stress components on the grid axes from the flat-bed stress.

        ``taus0``, ``taun0`` are the (uniform) flat-bed stress components in
        the grid frame. The along-wind stress is limited to >= 0 (no reversed
        flow in linear theory).
        '''
        tau0 = np.hypot(taus0, taun0)
        if tau0 == 0.:
            return np.zeros((self.ny, self.nx)), np.zeros((self.ny, self.nx))
        theta = np.arctan2(taun0, taus0)
        dpar, dperp = self.perturbation(zb, theta, workers)
        tpar = np.maximum(tau0 * (1. + dpar), 0.)
        tperp = tau0 * dperp
        c, s = np.cos(theta), np.sin(theta)
        return c * tpar - s * tperp, s * tpar + c * tperp
