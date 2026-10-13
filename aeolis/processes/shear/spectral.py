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
linear (non-periodic) by padding the bed to about twice its size. A
best-fit plane is removed first: a uniform slope produces no perturbation
in linear theory, and removing it makes the edges of the residual much
smaller. The pad continues the edge values of the residual and is blended
smoothly to its boundary mean in both directions. See
docs/overhaul/SHEAR_ON_MAIN_GRID.md.

Everything here works on plain numpy arrays and can be used outside AeoLiS:

>>> sh = SpectralShear(ny=200, nx=300, ds=2., dn=2., L=100., l=1., z0=0.001)
>>> dtau_par, dtau_perp = sh.perturbation(zb, theta=np.deg2rad(20.))
"""

import numpy as np
import scipy.fft
from scipy.special import kve

KAPPA = 0.41


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
        Padded FFT size relative to the grid (>= 2 for a non-circular result).
    blend_length : float, optional
        Distance over which the padded bed relaxes to the reference level
        [m]; default ``max(L/2, 10 grid cells)``, limited by the pad.
    filter_frac : float
        Isotropic low-pass cut-off as a fraction of the Nyquist wavenumber of
        the coarser grid direction (``None`` to disable).
    detrend : bool
        Remove the best-fit plane of the bed before padding (default). On the
        blowout case this reduces the boundary error against a much larger
        domain 6-9x compared to padding the bed itself.
    '''

    def __init__(self, ny, nx, ds, dn, L, l, z0, pad_factor=2., blend_length=None,
                 filter_frac=0.6, ntable=4096, detrend=True):
        self.ny, self.nx, self.ds, self.dn = int(ny), int(nx), float(ds), float(dn)
        self.L, self.l, self.z0 = float(L), float(l), float(z0)

        self.NY = scipy.fft.next_fast_len(max(int(np.ceil(pad_factor * ny)), ny + 2), real=True)
        self.NX = scipy.fft.next_fast_len(max(int(np.ceil(pad_factor * nx)), nx + 2), real=True)
        self.oy = (self.NY - ny) // 2          # offset of the model grid inside the pad
        self.ox = (self.NX - nx) // 2

        # separable taper: 1 on the model grid, extrusion zone, cosine blend to 0
        if blend_length is None:
            blend_length = max(L / 2., 10. * max(ds, dn))
        self.wy = self._taper(self.NY, self.oy, ny, dn, blend_length)
        self.wx = self._taper(self.NX, self.ox, nx, ds, blend_length)

        # wavenumbers of the padded grid (rfft along s)
        self.kx = 2. * np.pi * scipy.fft.rfftfreq(self.NX, ds)[np.newaxis, :]
        self.ky = 2. * np.pi * scipy.fft.fftfreq(self.NY, dn)[:, np.newaxis]
        k = np.hypot(self.kx, self.ky)
        if filter_frac:
            kc = filter_frac * np.pi / max(ds, dn)
            self.lowpass = np.exp(-(k / kc) ** 8)
        else:
            self.lowpass = np.ones_like(k)

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
        self._cache = {}

        # least-squares plane fit z ~ c0 + c1 * i + c2 * j as a precomputed projection
        self.detrend = detrend
        j, i = np.mgrid[0:ny, 0:nx]
        A = np.c_[np.ones(ny * nx), i.ravel(), j.ravel()].astype(float)
        self._plane_A = A
        self._plane_pinv = np.linalg.pinv(A)

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
        R = self._table(kpar, self._R, ramp=True)
        S = self._table(kpar, self._S)
        T_par = self._c * (-kpar ** 2 / ksafe + (self._loglz0 * kpar ** 2 / ksafe + k) * R)
        T_perp = self._c * (kpar * kperp / ksafe) * S
        T_par[k == 0.] = 0.
        T_perp[k == 0.] = 0.
        out = (T_par * self.lowpass, T_perp * self.lowpass)
        if len(self._cache) > 64:
            self._cache.pop(next(iter(self._cache)))
        self._cache[key] = out
        return out

    def pad(self, zb):
        '''Bed on the padded grid: plane removed, edge extrusion blended to the boundary mean.'''
        zb = np.asarray(zb, dtype=float)
        if self.detrend:
            zb = zb - (self._plane_A @ (self._plane_pinv @ zb.ravel())).reshape(zb.shape)
        ring = np.concatenate([zb[0, :], zb[-1, :], zb[:, 0], zb[:, -1]])
        zref = ring.mean()
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
