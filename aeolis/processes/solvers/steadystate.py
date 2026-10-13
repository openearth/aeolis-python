'''Steady-state advection-pickup solve in upwind order.

The steady-state solver of AeoLiS solves, per sediment fraction, the upwind
finite-volume balance of airborne sediment

    sum_out |u_f| L_f Ct  -  sum_in |u_f| L_f Ct_upwind  =  (w Cu(Ct) - Ct) A / T

with the pickup limited by the mass in the top bed layer. Every cell depends
only on its upwind neighbours, so the cells form a directed graph with an edge
from a cell to each neighbour it feeds. The legacy solver visits the grid in a
fixed order per wind quadrant and repeats the whole sweep until nothing changes
more than ``max_error``. Wherever a streamline switches quadrant, the downwind
cell reads a stale value and the information lags one sweep, so a wind that
meanders around a grid axis needs tens of sweeps.

Here the cells are visited in topological order of the dependency graph
instead. Strongly connected components (recirculation cells, a few per grid)
are iterated locally until converged; everything else is solved exactly in one
visit. Periodic boundaries are cut: the wrap neighbours are ghost values that
are refreshed per pass, and passes repeat until the wrap has converged.

The per-cell arithmetic is that of the quadrant kernels in
:mod:`aeolis.advection`, generalised to any combination of face signs per
fraction: inflow and outflow are decided per face and per fraction, and n-faces
carry the face length ``ds`` and s-faces ``dn`` for every quadrant.

All arrays are on the ghosted grid (one ghost cell on each side).
'''

import numpy as np
from numba import njit, prange


# ----------------------------------------------------------------------------
# boundary pass and fraction weights (shared with the legacy quadrant sweeps)
# ----------------------------------------------------------------------------

def _bc_code(name):
    '''Boundary condition name -> integer code for compiled code.

    As in :func:`aeolis.advection.apply_boundary`: 'circular' and 'flux' are
    special, every other name (e.g. 'constant', 'gradient') copies the
    neighbouring value.
    '''
    return {'circular': 0, 'flux': 1}.get(name, 2)


@njit(cache=True)
def _apply_boundary_nb(Ct, Cu, ufs, ufn, off_bc, on_bc, lat_bc,
                       offshore_flux, onshore_flux, lateral_flux):
    '''Compiled equivalent of :func:`aeolis.advection.apply_boundary` (same order of updates).'''
    ny, nx, nf = Ct.shape

    if off_bc == 0:   # circular in s (both offshore and onshore)
        for j in range(ny):
            for f in range(nf):
                Ct[j, 0, f] = Ct[j, nx - 2, f]
                Ct[j, nx - 1, f] = Ct[j, 1, f]
                Cu[j, 0, f] = Cu[j, nx - 2, f]
                Cu[j, nx - 1, f] = Cu[j, 1, f]
    else:
        for j in range(ny):
            for f in range(nf):
                v_in = offshore_flux * Cu[j, 1, f] if off_bc == 1 else Ct[j, 1, f]
                Ct[j, 0, f] = v_in if ufs[j, 1, f] >= 0 else Ct[j, 1, f]
        for j in range(ny):
            for f in range(nf):
                v_in = onshore_flux * Cu[j, nx - 2, f] if on_bc == 1 else Ct[j, nx - 2, f]
                Ct[j, nx - 1, f] = v_in if ufs[j, nx - 1, f] < 0 else Ct[j, nx - 2, f]
        for j in range(ny):
            for f in range(nf):
                Cu[j, 0, f] = Cu[j, 1, f]
                Cu[j, nx - 1, f] = Cu[j, nx - 2, f]

    if lat_bc == 0:
        for i in range(nx):
            for f in range(nf):
                Ct[0, i, f] = Ct[ny - 2, i, f]
                Ct[ny - 1, i, f] = Ct[1, i, f]
                Cu[0, i, f] = Cu[ny - 2, i, f]
                Cu[ny - 1, i, f] = Cu[1, i, f]
    else:
        for i in range(nx):
            for f in range(nf):
                v_in = lateral_flux * Cu[1, i, f] if lat_bc == 1 else Ct[1, i, f]
                Ct[0, i, f] = v_in if ufn[1, i, f] >= 0 else Ct[1, i, f]
        for i in range(nx):
            for f in range(nf):
                v_in = lateral_flux * Cu[ny - 2, i, f] if lat_bc == 1 else Ct[ny - 2, i, f]
                Ct[ny - 1, i, f] = v_in if ufn[ny - 1, i, f] < 0 else Ct[ny - 2, i, f]
        for i in range(nx):
            for f in range(nf):
                Cu[0, i, f] = Cu[1, i, f]
                Cu[ny - 1, i, f] = Cu[ny - 2, i, f]


@njit(cache=True)
def update_weights(Ct, Cu, mass, bi):
    ''' Update weights based on current Ct, Cu and bed mass.
    '''
    ny, nx, nf = Ct.shape
    w = np.zeros((ny, nx, nf))

    # loop over all cells
    for iy in range(ny):
        for ix in range(nx):

            # air contribution
            sum_air = 0.0
            for f in range(nf):
                ca = Cu[iy, ix, f]
                if ca != 0.0:
                    w[iy, ix, f] = Ct[iy, ix, f] / ca
                    sum_air += w[iy, ix, f]
                else:
                    w[iy, ix, f] = 0.0

            # bed contribution
            sum_bed = 0.0
            for f in range(nf):
                sum_bed += mass[iy, ix, 0, f]
            scale_bed = 1.0 - min(1.0, (1.0 - bi) * sum_air)

            for f in range(nf):
                w[iy, ix, f] = (
                    (1.0 - bi) * w[iy, ix, f] +
                    (mass[iy, ix, 0, f] / sum_bed if sum_bed > 0.0 else 0.0)
                    * scale_bed
                )

            # normalize weights
            s = 0.0
            for f in range(nf):
                s += w[iy, ix, f]

            if s > 0.0:
                for f in range(nf):
                    w[iy, ix, f] /= s

    return w


# ----------------------------------------------------------------------------
# ordering
# ----------------------------------------------------------------------------

@njit(cache=True)
def _tarjan(ufs, ufn, ny, nx):
    '''Strongly connected components of the physical cells (Tarjan, iterative).

    Edge c -> d when the flow crosses the face shared by c and d from c to d
    (periodic wraps are not edges). Returns the component id of every physical
    cell in row-major numbering and the component sizes. Tarjan emits the
    components in reverse topological order: a component comes out after every
    component downstream of it.
    '''
    W = nx - 2
    N = (ny - 2) * W
    index = np.full(N, -1, np.int32)
    low = np.zeros(N, np.int32)
    onstack = np.zeros(N, np.bool_)
    comp = np.full(N, -1, np.int32)
    stack = np.empty(N, np.int32)
    cs_node = np.empty(N, np.int32)     # explicit call stack: node ...
    cs_edge = np.empty(N, np.int32)     # ... and next edge to try
    sp = 0
    idx = 0
    ncomp = 0
    for root in range(N):
        if index[root] >= 0:
            continue
        depth = 0
        cs_node[0] = root
        cs_edge[0] = 0
        index[root] = idx
        low[root] = idx
        idx += 1
        stack[sp] = root
        sp += 1
        onstack[root] = True
        while depth >= 0:
            v = cs_node[depth]
            e = cs_edge[depth]
            n = v // W + 1
            s = v % W + 1
            pushed = False
            while e < 4:
                m = -1
                t = -1
                if e == 0:
                    if ufs[n, s + 1] > 0 and s + 1 <= nx - 2:
                        m = n
                        t = s + 1
                elif e == 1:
                    if ufs[n, s] < 0 and s - 1 >= 1:
                        m = n
                        t = s - 1
                elif e == 2:
                    if ufn[n + 1, s] > 0 and n + 1 <= ny - 2:
                        m = n + 1
                        t = s
                else:
                    if ufn[n, s] < 0 and n - 1 >= 1:
                        m = n - 1
                        t = s
                e += 1
                if m < 0:
                    continue
                u = (m - 1) * W + (t - 1)
                if index[u] < 0:
                    cs_edge[depth] = e
                    index[u] = idx
                    low[u] = idx
                    idx += 1
                    stack[sp] = u
                    sp += 1
                    onstack[u] = True
                    depth += 1
                    cs_node[depth] = u
                    cs_edge[depth] = 0
                    pushed = True
                    break
                elif onstack[u]:
                    if index[u] < low[v]:
                        low[v] = index[u]
            if pushed:
                continue
            if low[v] == index[v]:
                while True:
                    sp -= 1
                    u = stack[sp]
                    onstack[u] = False
                    comp[u] = ncomp
                    if u == v:
                        break
                ncomp += 1
            depth -= 1
            if depth >= 0:
                p = cs_node[depth]
                if low[v] < low[p]:
                    low[p] = low[v]
    sizes = np.zeros(ncomp, np.int64)
    for i in range(N):
        sizes[comp[i]] += 1
    return comp, sizes


@njit(cache=True)
def _heap_push(hk, hv, size, key, val):
    i = size
    hk[i] = key
    hv[i] = val
    while i > 0:
        p = (i - 1) // 2
        if hk[p] <= hk[i]:
            break
        hk[p], hk[i] = hk[i], hk[p]
        hv[p], hv[i] = hv[i], hv[p]
        i = p
    return size + 1


@njit(cache=True)
def _heap_pop(hk, hv, size):
    val = hv[0]
    size -= 1
    if size > 0:
        hk[0] = hk[size]
        hv[0] = hv[size]
        i = 0
        while True:
            l = 2 * i + 1
            r = l + 1
            m = i
            if l < size and hk[l] < hk[m]:
                m = l
            if r < size and hk[r] < hk[m]:
                m = r
            if m == i:
                break
            hk[m], hk[i] = hk[i], hk[m]
            hv[m], hv[i] = hv[i], hv[m]
            i = m
    return val, size


@njit(cache=True)
def ordering_plan(ufs, ufn):
    '''Visiting plan for one fraction.

    The strongly connected components of the upwind dependency graph are the
    nodes of a directed acyclic graph. They are released in topological order
    (Kahn's algorithm) with a priority queue keyed on the position of the
    component in a row-major scan aligned with the dominant flow direction, so
    the visiting order follows the memory layout wherever the dependencies
    allow it. Recirculation cells (components of more than one cell) are
    visited as one block, in scan order.

    Parameters
    ----------
    ufs : ndarray (ny, nx+1)
        Face velocities on the s-faces of the ghosted grid.
    ufn : ndarray (ny+1, nx)
        Face velocities on the n-faces of the ghosted grid.

    Returns
    -------
    order : ndarray (N,) int32
        Physical cells (row-major index over the (ny-2, nx-2) interior),
        grouped by strongly connected component, components in upwind
        topological order.
    cstart : ndarray (ncomp+1,) int64
        Component c occupies ``order[cstart[c]:cstart[c+1]]``.
    '''
    ny = ufs.shape[0]
    nx = ufn.shape[1]
    W = nx - 2
    H = ny - 2
    N = H * W
    comp, sizes = _tarjan(ufs, ufn, ny, nx)
    ncomp = sizes.size

    # scan frame: rows in the direction of the mean n-velocity, cells within a
    # row in the direction of the mean s-velocity
    sum_s = 0.0
    sum_n = 0.0
    for n in range(1, ny - 1):
        for s in range(1, nx - 1):
            sum_s += ufs[n, s]
            sum_n += ufn[n, s]
    flip_s = sum_s < 0.0
    flip_n = sum_n < 0.0

    # scan key per cell, component key = smallest member key, members grouped
    key = np.empty(N, np.int64)
    ckey = np.full(ncomp, N, np.int64)
    for v in range(N):
        n = v // W
        s = v % W
        kn = (H - 1 - n) if flip_n else n
        ks = (W - 1 - s) if flip_s else s
        key[v] = kn * W + ks
        if key[v] < ckey[comp[v]]:
            ckey[comp[v]] = key[v]
    moff = np.zeros(ncomp + 1, np.int64)
    for c in range(ncomp):
        moff[c + 1] = moff[c] + sizes[c]
    fill = moff[:-1].copy()
    members = np.empty(N, np.int32)
    for v in range(N):
        members[fill[comp[v]]] = v
        fill[comp[v]] += 1
    # members of a recirculation block in scan order (insertion sort, blocks are small)
    for c in range(ncomp):
        a = moff[c]
        b = moff[c + 1]
        if b - a > 1:
            for i in range(a + 1, b):
                v = members[i]
                j = i - 1
                while j >= a and key[members[j]] > key[v]:
                    members[j + 1] = members[j]
                    j -= 1
                members[j + 1] = v

    # in-degree of every component: edges from other components
    indeg = np.zeros(ncomp, np.int64)
    for v in range(N):
        n = v // W + 1
        s = v % W + 1
        c = comp[v]
        if ufs[n, s + 1] > 0 and s + 1 <= nx - 2:
            u = comp[v + 1]
            if u != c:
                indeg[u] += 1
        if ufs[n, s] < 0 and s - 1 >= 1:
            u = comp[v - 1]
            if u != c:
                indeg[u] += 1
        if ufn[n + 1, s] > 0 and n + 1 <= ny - 2:
            u = comp[v + W]
            if u != c:
                indeg[u] += 1
        if ufn[n, s] < 0 and n - 1 >= 1:
            u = comp[v - W]
            if u != c:
                indeg[u] += 1

    # Kahn with a priority queue on the scan key
    hk = np.empty(ncomp, np.int64)
    hv = np.empty(ncomp, np.int64)
    hsize = 0
    for c in range(ncomp):
        if indeg[c] == 0:
            hsize = _heap_push(hk, hv, hsize, ckey[c], c)
    order = np.empty(N, np.int32)
    cstart = np.empty(ncomp + 1, np.int64)
    pos = 0
    cc = 0
    while hsize > 0:
        c, hsize = _heap_pop(hk, hv, hsize)
        cstart[cc] = pos
        cc += 1
        for i in range(moff[c], moff[c + 1]):
            v = members[i]
            order[pos] = v
            pos += 1
            n = v // W + 1
            s = v % W + 1
            if ufs[n, s + 1] > 0 and s + 1 <= nx - 2:
                u = comp[v + 1]
                if u != c:
                    indeg[u] -= 1
                    if indeg[u] == 0:
                        hsize = _heap_push(hk, hv, hsize, ckey[u], u)
            if ufs[n, s] < 0 and s - 1 >= 1:
                u = comp[v - 1]
                if u != c:
                    indeg[u] -= 1
                    if indeg[u] == 0:
                        hsize = _heap_push(hk, hv, hsize, ckey[u], u)
            if ufn[n + 1, s] > 0 and n + 1 <= ny - 2:
                u = comp[v + W]
                if u != c:
                    indeg[u] -= 1
                    if indeg[u] == 0:
                        hsize = _heap_push(hk, hv, hsize, ckey[u], u)
            if ufn[n, s] < 0 and n - 1 >= 1:
                u = comp[v - W]
                if u != c:
                    indeg[u] -= 1
                    if indeg[u] == 0:
                        hsize = _heap_push(hk, hv, hsize, ckey[u], u)
    cstart[cc] = pos
    return order, cstart[:cc + 1]


@njit(cache=True, parallel=True)
def ordering_plans(ufs, ufn):
    '''Visiting plans of all fractions, built in parallel.

    Parameters
    ----------
    ufs : ndarray (ny, nx+1, nf), ufn : ndarray (ny+1, nx, nf)
        Face velocities on the ghosted grid.

    Returns
    -------
    orders : ndarray (nf, N) int32
    cstarts : ndarray (nf, N+1) int64
        Row f holds ``ncomps[f] + 1`` valid entries.
    ncomps : ndarray (nf,) int64
    '''
    ny, nx, nf = ufs.shape[0], ufn.shape[1], ufs.shape[2]
    N = (ny - 2) * (nx - 2)
    orders = np.empty((nf, N), np.int32)
    cstarts = np.zeros((nf, N + 1), np.int64)
    ncomps = np.zeros(nf, np.int64)
    for f in prange(nf):
        o, c = ordering_plan(np.ascontiguousarray(ufs[:, :, f]), np.ascontiguousarray(ufn[:, :, f]))
        orders[f, :] = o
        cstarts[f, :c.size] = c
        ncomps[f] = c.size - 1
    return orders, cstarts, ncomps


# ----------------------------------------------------------------------------
# per-cell update
# ----------------------------------------------------------------------------

@njit(cache=True, inline='always')
def _solve_cell(n, s, f, Ct, Cu_air, Cu_bed, zeta, mass, pickup, pickup0, w,
                dt, Ts, ds, dn, ufs, ufn):
    '''Quadrant-kernel arithmetic for one cell and fraction, any face signs.

    Returns the absolute change of Ct in the cell. Inlined into the pass: a
    call with this many array arguments costs more than the arithmetic.
    '''
    A = ds[n, s] * dn[n, s]
    wf = w[n, s, f]

    if Cu_air[n, s, f] > 0 and Ct[n, s, f] > 0:
        a = (1 - zeta[n, s]) * (Cu_air[n, s, f] - Cu_bed[n, s, f]) / Cu_air[n, s, f]
        b = Cu_bed[n, s, f]
    else:
        a = 0.0
        b = Cu_bed[n, s, f]

    # inflow faces feed N, outflow faces add to D; n-faces have length ds,
    # s-faces have length dn
    N = 0.0
    D = A / Ts
    if ufn[n, s, f] > 0:
        N += Ct[n - 1, s, f] * ufn[n, s, f] * ds[n, s]
    else:
        D += -ufn[n, s, f] * ds[n, s]
    if ufs[n, s, f] > 0:
        N += Ct[n, s - 1, f] * ufs[n, s, f] * dn[n, s]
    else:
        D += -ufs[n, s, f] * dn[n, s]
    if ufn[n + 1, s, f] < 0:
        N += -Ct[n + 1, s, f] * ufn[n + 1, s, f] * ds[n, s]
    else:
        D += ufn[n + 1, s, f] * ds[n, s]
    if ufs[n, s + 1, f] < 0:
        N += -Ct[n, s + 1, f] * ufs[n, s + 1, f] * dn[n, s]
    else:
        D += ufs[n, s + 1, f] * dn[n, s]

    wa = wf * a
    if wa > 0.999:
        wa = 0.999
    D_eff = D - wa * A / Ts
    if D_eff < 1e-12:
        D_eff = 1e-12

    Ct_new = (N + wf * b * A / Ts) / D_eff

    max_C = mass[n, s, 0, f] / A + 1.0
    if Ct_new > max_C:
        Ct_new = max_C

    p = (wf * (b + a * Ct_new) - Ct_new) * dt / Ts
    pickup0[n, s, f] = abs(p)

    if p > mass[n, s, 0, f]:
        p_lim = mass[n, s, 0, f]
        Ct_new = (N + p_lim * A / dt) / D_eff
        p = p_lim - (1.0 - wa) * Ct_new * dt / Ts

    d = abs(Ct_new - Ct[n, s, f])
    Ct[n, s, f] = Ct_new
    pickup[n, s, f] = p
    return d


@njit(cache=True)
def _ordered_pass(order, cstart, f, tol, max_inner, Ct, Cu_air, Cu_bed, zeta, mass,
                  pickup, pickup0, w, dt, Ts, ds, dn, ufs, ufn, inner):
    '''One pass over all physical cells of fraction f in plan order.

    Components of one cell are solved once; larger components are Gauss-Seidel
    iterated until their change is <= tol (or max_inner iterations). ``inner``
    receives the iteration count per cell. Returns the maximum change in the
    cells' first visit of this pass.

    The cell arithmetic is written out here (the same as :func:`_solve_cell`):
    a call per cell with this many array arguments costs several times the
    arithmetic itself, even when numba inlines it.
    '''
    W = ufn.shape[1] - 2
    dmax = 0.0
    for c in range(cstart.size - 1):
        a = cstart[c]
        b = cstart[c + 1]
        single = b - a == 1
        k = 0
        while True:
            dm = 0.0
            for i in range(a, b):
                v = order[i]
                n = v // W + 1
                s = v % W + 1

                A = ds[n, s] * dn[n, s]
                wf = w[n, s, f]
                if Cu_air[n, s, f] > 0 and Ct[n, s, f] > 0:
                    aa = (1 - zeta[n, s]) * (Cu_air[n, s, f] - Cu_bed[n, s, f]) / Cu_air[n, s, f]
                    bb = Cu_bed[n, s, f]
                else:
                    aa = 0.0
                    bb = Cu_bed[n, s, f]
                N = 0.0
                D = A / Ts
                if ufn[n, s, f] > 0:
                    N += Ct[n - 1, s, f] * ufn[n, s, f] * ds[n, s]
                else:
                    D += -ufn[n, s, f] * ds[n, s]
                if ufs[n, s, f] > 0:
                    N += Ct[n, s - 1, f] * ufs[n, s, f] * dn[n, s]
                else:
                    D += -ufs[n, s, f] * dn[n, s]
                if ufn[n + 1, s, f] < 0:
                    N += -Ct[n + 1, s, f] * ufn[n + 1, s, f] * ds[n, s]
                else:
                    D += ufn[n + 1, s, f] * ds[n, s]
                if ufs[n, s + 1, f] < 0:
                    N += -Ct[n, s + 1, f] * ufs[n, s + 1, f] * dn[n, s]
                else:
                    D += ufs[n, s + 1, f] * dn[n, s]
                wa = wf * aa
                if wa > 0.999:
                    wa = 0.999
                D_eff = D - wa * A / Ts
                if D_eff < 1e-12:
                    D_eff = 1e-12
                Ct_new = (N + wf * bb * A / Ts) / D_eff
                max_C = mass[n, s, 0, f] / A + 1.0
                if Ct_new > max_C:
                    Ct_new = max_C
                p = (wf * (bb + aa * Ct_new) - Ct_new) * dt / Ts
                pickup0[n, s, f] = abs(p)
                if p > mass[n, s, 0, f]:
                    p_lim = mass[n, s, 0, f]
                    Ct_new = (N + p_lim * A / dt) / D_eff
                    p = p_lim - (1.0 - wa) * Ct_new * dt / Ts
                d = abs(Ct_new - Ct[n, s, f])
                Ct[n, s, f] = Ct_new
                pickup[n, s, f] = p
                if d > dm:
                    dm = d
            k += 1
            if k == 1 and dm > dmax:
                dmax = dm
            if single or dm <= tol or k >= max_inner:
                break
        for i in range(a, b):
            v = order[i]
            inner[v // W + 1, v % W + 1, f] = k
    return dmax


# ----------------------------------------------------------------------------
# driver
# ----------------------------------------------------------------------------

@njit(cache=True)
def solve_ordered(Ct, Ct_last, Cu_air, Cu_bed, zeta, mass, pickup, pickup0, w,
                  dt, Ts, ds, dn, ufs, ufn, orders, cstarts, ncomps,
                  off_bc, on_bc, lat_bc,
                  offshore_flux, onshore_flux, lateral_flux,
                  nf_weights, bi,
                  max_iter, max_error, iters, inner, history):
    '''Steady-state solve in plan order; the loop structure of the legacy driver.

    ``orders`` (nf, N) and ``cstarts`` (nf, maxcomp+1) hold the plan of every
    fraction (``ncomps[f]`` components). With ``nf_weights`` the fraction
    weights are refreshed from Ct before every pass (bed interaction ``bi``).
    Passes repeat until no cell changes by more than ``max_error`` (one pass
    plus a confirming pass for a grid without periodic boundaries). Returns the
    number of passes and whether ``max_iter`` was hit.
    '''
    ny, nx, nf = Ct.shape
    Cu = np.empty_like(Cu_air)
    k = 0
    stopped = False
    inner_tol = max_error * 1e-3    # components converge well below the outer tolerance
    # a second pass can only change something when a value of the first pass
    # is fed back into it: a ghost that is not a constant (circular wraps and
    # 'constant'/'gradient' copies of the boundary cell itself; only 'flux'
    # ghosts are constants), Ct-dependent weights, or the clip
    lagged = (off_bc != 1 or on_bc != 1 or lat_bc != 1 or (nf_weights and bi < 1.0))
    cmax = np.empty(nf)

    while True:
        Ct_last[:] = Ct

        _apply_boundary_nb(Ct, Cu_air, ufs, ufn, off_bc, on_bc, lat_bc,
                           offshore_flux, onshore_flux, lateral_flux)

        # physical limits: 0 <= Ct <= max(Cu_air) per fraction
        for f in range(nf):
            if k == 0:
                cm = -np.inf
                has_nan = False
                for j in range(ny):
                    for i in range(nx):
                        v = Cu_air[j, i, f]
                        if np.isnan(v):
                            has_nan = True
                        elif v > cm:
                            cm = v
                cmax[f] = np.nan if has_nan else cm
            cm = cmax[f]
            for j in range(ny):
                for i in range(nx):
                    x = Ct[j, i, f]
                    if not np.isnan(x):
                        x = x if x > 0.0 else 0.0
                        x = x if x < cm else cm
                    Ct[j, i, f] = x

        if nf_weights:
            for j in range(ny):
                for i in range(nx):
                    for f in range(nf):
                        v = Cu_air[j, i, f]
                        Cu[j, i, f] = v if (v >= 0.0 or np.isnan(v)) else 0.0
            w[:] = update_weights(Ct, Cu, mass, bi)

        for f in range(nf):
            _ordered_pass(orders[f], cstarts[f, :ncomps[f] + 1], f, inner_tol, 1000,
                          Ct, Cu_air, Cu_bed, zeta, mass, pickup, pickup0, w,
                          dt, Ts, ds, dn, ufs, ufn, inner)
        k += 1

        notconverged = False
        in_range = True
        dmax = -1.0
        dval = 0.0
        for j in range(ny):
            for i in range(nx):
                for f in range(nf):
                    x = Ct[j, i, f]
                    d = x - Ct_last[j, i, f]
                    ad = abs(d)
                    if ad > max_error:
                        iters[j, i, f] = k
                        notconverged = True
                    if ad > dmax:
                        dmax = ad
                        dval = d
                    if not (x >= 0.0 and x <= cmax[f]):
                        in_range = False
        history[min(k - 1, history.shape[0] - 1)] = dval

        if k > max_iter:
            stopped = True
            break
        if not notconverged:
            break
        if not lagged and in_range:
            # nothing of this pass feeds back into the next one: the next pass
            # would reproduce these values exactly, so it is not run and the
            # reported last change is zero, as it would have been
            Ct_last[:] = Ct
            break

    return k, stopped
