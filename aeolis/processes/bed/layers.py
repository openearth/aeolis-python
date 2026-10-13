"""Per-column kernels for the layered sediment bed.

The bed in every grid cell is a column of ``nl`` layers with constant total
mass; each layer holds a mass per sediment fraction (``nf``). These kernels
update one column at a time, so they only touch the memory of that column
and can run in parallel over cells.

They reproduce the original vectorised implementation in
:func:`aeolis.bed.update` and :func:`aeolis.bed.mixtoplayer` operation by
operation (same arithmetic in the same order), so results are bitwise
identical. Keep it that way when editing: every change of the arithmetic
here is a change of model results.

Array shapes (flattened grid, ``nc = (ny+1)*(nx+1)`` cells):

* ``mass``   (nc, nl, nf) [kg/m2]  - updated in place
* ``pickup`` (nc, nf)     [kg/m2]  - net pickup this step (>0 erosion); updated in place
* ``thlyr``  (nc, nl)     [m]      - layer thickness
"""

import math

import numpy as np
from numba import njit, prange


@njit(cache=True, error_model='numpy')
def _floor_divide(a, b):
    '''numpy.floor_divide for two floats (numpy's npy_floor_divide).'''
    if b == 0.:
        return a / b if a != 0. else np.nan
    mod = np.fmod(a, b)
    div = (a - mod) / b
    if mod != 0.:
        if (b < 0.) != (mod < 0.):
            div -= 1.0
    if div != 0.:
        floordiv = math.floor(div)
        if div - floordiv > 0.5:
            floordiv += 1.0
    else:
        floordiv = math.copysign(0., a / b)
    return floordiv


@njit(cache=True)
def _pairwise_block(a, start, n):
    '''numpy's pairwise summation of ``a[start:start+n]`` for ``n <= 128``.'''
    if n < 8:
        res = 0.
        for i in range(n):
            res += a[start + i]
        return res
    r0 = a[start]
    r1 = a[start + 1]
    r2 = a[start + 2]
    r3 = a[start + 3]
    r4 = a[start + 4]
    r5 = a[start + 5]
    r6 = a[start + 6]
    r7 = a[start + 7]
    i = 8
    while i < n - (n % 8):
        r0 += a[start + i]
        r1 += a[start + i + 1]
        r2 += a[start + i + 2]
        r3 += a[start + i + 3]
        r4 += a[start + i + 4]
        r5 += a[start + i + 5]
        r6 += a[start + i + 6]
        r7 += a[start + i + 7]
        i += 8
    res = ((r0 + r1) + (r2 + r3)) + ((r4 + r5) + (r6 + r7))
    while i < n:
        res += a[start + i]
        i += 1
    return res


@njit(cache=True)
def _split(n):
    n2 = n // 2
    return n2 - n2 % 8


@njit(cache=True)
def _pairwise(a, start, n):
    '''numpy's pairwise summation of ``a[start:start+n]``.

    numpy halves blocks larger than 128 recursively; numba cannot cache
    functions that call recursive functions, so the first three levels of
    that recursion are written out (exact for n <= 1024).
    '''
    if n <= 128:
        return _pairwise_block(a, start, n)
    if n > 1024:
        raise ValueError('pairwise summation emulated up to 1024 values')
    h = _split(n)
    return _pairwise_l2(a, start, h) + _pairwise_l2(a, start + h, n - h)


@njit(cache=True)
def _pairwise_l2(a, start, n):
    if n <= 128:
        return _pairwise_block(a, start, n)
    h = _split(n)
    return _pairwise_l3(a, start, h) + _pairwise_l3(a, start + h, n - h)


@njit(cache=True)
def _pairwise_l3(a, start, n):
    if n <= 128:
        return _pairwise_block(a, start, n)
    h = _split(n)
    return _pairwise_block(a, start, h) + _pairwise_block(a, start + h, n - h)


@njit(cache=True)
def _sum(x):
    '''Sum of a 1D array, bitwise identical to ``numpy.sum`` along a contiguous axis.'''
    return _pairwise(x, 0, x.shape[0])


@njit(cache=True, parallel=True, error_model='numpy')
def update_layers(mass, pickup, gd_bottom, max_error):
    '''Move sediment between bed layers after erosion and deposition.

    Equivalent to the layer bookkeeping in :func:`aeolis.bed.update` for a
    uniform bottom-layer grain size distribution. Sediment removed from (or
    added to) the top layer is replenished from (or pushed to) the layers
    below, keeping the mass per layer constant. The bottom layer exchanges
    sediment with an infinite source of composition ``gd_bottom``.

    Parameters
    ----------
    mass : ndarray (nc, nl, nf)
        Bed mass per layer and fraction; updated in place.
    pickup : ndarray (nc, nf)
        Net sediment pickup per fraction (positive is erosion); updated in
        place exactly as the original routine does.
    gd_bottom : ndarray (nf,)
        Normalised grain size distribution of the sediment source below the
        bottom layer.
    max_error : float
        Negative masses smaller than this in magnitude are set to zero.

    Returns
    -------
    dm : ndarray (nc,)
        Total mass exchanged between the layers (positive is deposition).
    burial : ndarray (nc,)
        Total mass pushed out of the bottom layer by large deposition.
    '''

    nc, nl, nf = mass.shape
    dm_out = np.zeros(nc)
    burial_out = np.zeros(nc)

    for c in prange(nc):
        m = mass[c]
        pk = pickup[c]

        # no pickup in any fraction: the bookkeeping below is the identity
        # apart from the tiny-negative clean-up, so only that is done (most
        # cells of a large grid have no transport in a given step)
        active = False
        for f in range(nf):
            if pk[f] != 0.:
                active = True
        if not active:
            for l in range(nl):
                for f in range(nf):
                    if m[l, f] < 0. and m[l, f] > -max_error:
                        m[l, f] = 0.
            continue

        # net exchange; erosion/deposition classification uses the value
        # *before* the negative-mass corrections, as in the original
        dm = -_sum(pk)
        is_ero = dm < 0.
        is_dep = dm > 0.

        burial = np.zeros(nf)

        # --- case 1: net deposition cell with some erosional fractions
        if is_dep:
            erog = 0.
            for f in range(nf):
                erog += max(0., pk[f])
            if erog > 0.:
                dep = np.empty(nf)
                for f in range(nf):
                    dep[f] = -min(0., pk[f])
                depsum = _sum(dep)
                for f in range(nf):
                    d = dep[f] / depsum if depsum != 0. else 0.
                    ddep = erog * d
                    ero = max(0., pk[f])
                    pk[f] = -dep[f] + ddep
                    m[0, f] -= ero - ddep
                dm = -_sum(pk)

        # --- case 2: deposition larger than the mass of the top layer
        mx = _sum(m[0])
        n = _floor_divide(dm, mx)
        if n > 0.:
            pksum = _sum(pk)
            d2 = np.empty(nf)
            for f in range(nf):
                d2[f] = pk[f] / pksum if pksum != 0. else 0.
            for i in range(nl):
                if not n > i:
                    break
                for f in range(nf):
                    burial[f] += m[nl - 1, f]
                for k in range(nl - 1, i, -1):      # shift layers below i down by one
                    for f in range(nf):
                        m[k, f] = m[k - 1, f]
                for f in range(nf):
                    m[i, f] = mx * d2[f]
                    pk[f] += m[i, f]
            dm = -_sum(pk)

        # --- distribution of each layer over the fractions
        d = np.empty((nl, nf))
        for l in range(nl):
            ref = _sum(m[l])
            for f in range(nf):
                d[l, f] = m[l, f] / ref if ref != 0. else 0.

        # --- remove pickup from the top layer and rearrange layers
        for f in range(nf):
            m[0, f] -= pk[f]
        if is_ero:
            for i in range(1, nl):
                for f in range(nf):
                    m[i - 1, f] -= dm * d[i, f]
                    m[i, f] += dm * d[i, f]
            for f in range(nf):
                m[nl - 1, f] -= dm * gd_bottom[f]
        elif is_dep:
            for i in range(1, nl):
                for f in range(nf):
                    m[i - 1, f] -= dm * d[i - 1, f]
                    m[i, f] += dm * d[i - 1, f]
            for f in range(nf):
                m[nl - 1, f] -= dm * d[nl - 1, f]

        # --- remove tiny negatives
        for l in range(nl):
            for f in range(nf):
                if m[l, f] < 0. and m[l, f] > -max_error:
                    m[l, f] = 0.

        dm_out[c] = dm
        burial_out[c] = _sum(burial)

    return dm_out, burial_out


@njit(cache=True, parallel=True, error_model='numpy')
def mix_toplayers(mass, thlyr, dod):
    '''Mix the grain size distribution over the layers above the depth of disturbance.

    Equivalent to :func:`aeolis.bed.mixtoplayer`. Only cells with a
    positive depth of disturbance are changed (for all other cells the
    original routine reproduces the top layer exactly).

    Parameters
    ----------
    mass : ndarray (nc, nl, nf)
        Bed mass per layer and fraction; updated in place.
    thlyr : ndarray (nc, nl)
        Layer thickness [m].
    dod : ndarray (nc,)
        Depth of disturbance [m] (``facDOD * Hsmix``).
    '''

    nc, nl, nf = mass.shape
    for c in prange(nc):
        if not dod[c] > 0.:
            continue
        m = mass[c]

        f_mix = min(1., _sum(thlyr[c]) / dod[c])

        # layers above the depth of disturbance (the top layer always mixes)
        nmix = 0
        cum = 0.
        mixed = np.zeros(nl, dtype=np.bool_)
        for l in range(nl):
            cum += thlyr[c, l]
            if cum <= dod[c] or l == 0:
                mixed[l] = True
                nmix += 1

        for f in range(nf):
            tot = 0.
            for l in range(nl):
                if mixed[l]:
                    tot += m[l, f]
            mean = tot / nmix
            for l in range(nl):
                if mixed[l]:
                    m[l, f] = mean * f_mix + m[l, f] * (1. - f_mix)
