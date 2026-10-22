'''One row for a transect of identical rows (configuration key ``transect_1d``).

A cross-shore transect runs as a grid of identical rows with a circular
lateral boundary: AeoLiS turns a 1D input into three rows, and setup tools
often write three rows themselves. When every row holds the same state, a
row whose neighbours are copies of itself is the same problem as the whole
grid, so all rows but one are repeated work. With ``transect_1d = T`` the
model checks at its first time step whether the grid is such a transect and,
if so, computes one row. To everything outside the time step (output files,
the callback, ``get_var``/``set_var``, restart files) the model keeps all
its rows.

What the one row solves is the transect that is uniform alongshore. The
grid of identical rows gives the same, with three exceptions that are
properties of the many-row grid, not of the transect:

- The steady-state solver stops at ``max_error``. With an alongshore wind
  component the rows are coupled, and the rows of a many-row grid then
  differ within that tolerance (Gauss-Seidel order, lateral ghosts one
  iteration behind). The one row has no lateral coupling to iterate.
- The grass model draws its random numbers per cell, so the rows of a
  many-row grid are different realisations. The one row is one realisation.
- Round-off: the FFT of identical rows is exact only to round-off.

Kept on one row, unchanged from the many-row grid: circular lateral
boundaries in the solver (the ghost rows are the row itself), avalanching
(no slope across identical rows), the bed-slope terms, the shear on the model
grid (periodic along the rows), and the grass kernels and leeside wake, which
wrap across the lateral boundary (``grass_wrap_lateral``).

If the callback or ``set_var`` makes the rows differ, the model goes back to
all rows for the rest of the run and says so in the log.
'''

import logging
import numpy as np

logger = logging.getLogger('aeolis.model')      # the logger the runner writes to aeolis.log

#: coordinates: these differ between rows by a constant row step
GEOMETRY = ('x', 'y', 'x_vsub', 'y_vsub')


def _equal_rows(v, n):
    '''True if the n blocks of rows along axis 0 of v are identical (NaN equals NaN).'''
    b = v.reshape((n, v.shape[0] // n) + v.shape[1:])
    first = b[:1]
    if v.dtype.kind == 'f':
        return bool(np.all((b == first) | (np.isnan(b) & np.isnan(first))))
    return bool(np.all(b == first))


def _row_shifted(v, n):
    '''True if row r of v is row 0 shifted by r times one row step.'''
    if v.shape[0] != n or n < 2:
        return True
    step = v[1] - v[0]
    r = np.arange(n).reshape((n,) + (1,) * (v.ndim - 1))
    return bool(np.allclose(v, v[:1] + r * step[np.newaxis], rtol=0., atol=1e-9 * max(1., np.abs(v).max())))


class OneRow:
    '''The model state of a transect on one row, and the full grid it stands for.'''

    def __init__(self, model):
        s, p = model.s, model.p
        self.n = p['ny'] + 1                                    # rows of the full grid
        self.nx = p['nx'] + 1
        self.f = int(p['veg_res_factor']) if p['method_vegetation'] == 'grass' else 1
        self.flat = set()                                       # 1D arrays over all cells (row-major)
        self.full = {}                                          # coordinates on the full grid
        for k, v in s.items():
            if self.layout(k, v, self.n) == 'flat':
                self.flat.add(k)

    # --- shapes ---------------------------------------------------------------

    def layout(self, k, v, rows=1):
        '''How an array spans ``rows`` rows: 'row', 'sub' (vegetation subgrid), 'flat' or None.'''
        if not isinstance(v, np.ndarray):
            return None
        if v.ndim >= 2 and v.shape[0] == rows and v.shape[1] == self.nx:
            return 'row'
        if v.ndim >= 2 and self.f > 1 and v.shape[0] == rows * self.f and v.shape[1] == self.nx * self.f:
            return 'sub'
        if v.ndim == 1 and v.size == rows * self.nx and (rows > 1 or k in self.flat):
            return 'flat'
        return None

    def one(self, v, kind):
        '''The first row of a full-grid array.'''
        return v[:{'row': 1, 'sub': self.f, 'flat': self.nx}[kind]].copy()

    def expand(self, k, v):
        '''A one-row value on the full grid (a new array); anything else unchanged.'''
        if k in self.full:
            return self.full[k].copy()
        kind = self.layout(k, v)
        if kind is None:
            return v
        if kind == 'flat':
            return np.tile(v, self.n)
        return np.tile(v, (self.n,) + (1,) * (v.ndim - 1))

    def identical(self, k, v):
        '''True if a full-grid array has identical rows (coordinates: rows a row step apart).'''
        if k in GEOMETRY:
            return _row_shifted(v, v.shape[0])
        return _equal_rows(v, self.n)

    def expand_state(self, s):
        for k, v in list(s.items()):
            e = self.expand(k, v)
            if e is not v:
                dict.__setitem__(s, k, e)
        if isinstance(s.get('zb_hist'), list):
            dict.__setitem__(s, 'zb_hist', [(t, np.tile(z, (self.n, 1))) for t, z in s['zb_hist']])
        if hasattr(s.get('shear'), 'with_rows'):
            dict.__setitem__(s, 'shear', s['shear'].with_rows(self.n))

    # --- collapse and expand the model ------------------------------------

    def collapse(self, model):
        s, p = model.s, model.p
        for k, v in list(s.items()):
            kind = self.layout(k, v, self.n)
            if kind is None:
                continue
            if k in GEOMETRY:
                self.full[k] = v.copy()
            dict.__setitem__(s, k, self.one(v, kind))
        if isinstance(s.get('zb_hist'), list):
            dict.__setitem__(s, 'zb_hist', [(t, z[:1].copy()) for t, z in s['zb_hist']])
        if hasattr(s.get('shear'), 'with_rows'):
            dict.__setitem__(s, 'shear', s['shear'].with_rows(1))
        p['ny'] = 0
        p['ny_vsub'] = 0
        p['_one_row'] = True
        model.l = model.s.copy()

    def restore(self, model, given=None):
        '''Back to the full grid; ``given`` holds full-grid values that replace the expansion.'''
        s, p = model.s, model.p
        self.expand_state(s)
        for k, v in (given or {}).items():
            dict.__setitem__(s, k, v)
        p['ny'] = self.n - 1
        p['ny_vsub'] = (self.n - 1) * self.f
        p['_one_row'] = False
        model.l = model.s.copy()

    def full_state(self, model):
        '''Copies of (s, l, p) on the full grid, e.g. for a restart file.'''
        s = model.s.copy()
        self.expand_state(s)
        l = {k: self.expand(k, v) for k, v in model.l.items()}
        p = dict(model.p, ny=self.n - 1, ny_vsub=(self.n - 1) * self.f, _one_row=False)
        return s, l, p

    # --- the full grid as seen from outside -------------------------------

    def get(self, model, k):
        return self.expand(k, model.s[k])

    def shape(self, model, k):
        v = model.s[k]
        if self.layout(k, v) is None:
            return np.shape(v)
        return (v.shape[0] * self.n,) + v.shape[1:]

    def put(self, model, values, assigned=()):
        '''Take full-grid values back. Returns False (and restores the full grid) if any has rows that differ.'''
        differ, ones = [], {}
        for k, v in values.items():
            kind = self.layout(k, v, self.n)
            if kind is None:
                if k in model.s and self.layout(k, model.s[k]) is not None:
                    differ.append(k)            # a one-row variable replaced by something else
                elif k in assigned:
                    model.s[k] = v
                continue
            if self.identical(k, v):
                ones[k] = (kind, v)
            else:
                differ.append(k)
        if differ:
            logger.warning('transect_1d: %s no longer has identical rows; the model computes all %d rows '
                           'from here on' % (', '.join(sorted(differ)), self.n))
            self.restore(model, values)
            model._one_row = False
            return False
        for k, (kind, v) in ones.items():
            if kind == 'flat':
                self.flat.add(k)
            if k in GEOMETRY:
                self.full[k] = v.copy()
            one = self.one(v, kind)
            if k in assigned:
                model.s[k] = one                 # respects immutable state variables
            else:
                dict.__setitem__(model.s, k, one)
        return True

    def call(self, model, callback):
        '''Call ``callback(model)`` with the model on its full grid.'''
        view = FullRows(model.s, self)
        model.s, model._one_row = view, False        # get_var/set_var inside the callback use the view
        try:
            callback(model)
        finally:
            model.s, model._one_row = view.base, self
        self.put(model, view.touched, view.assigned)


class FullRows(dict):
    '''The model state as the full grid, expanded on first access (for the callback).'''

    def __init__(self, base, onerow):
        super().__init__(base)
        self.base, self.r = base, onerow
        self.touched, self.assigned = {}, set()

    def __getitem__(self, k):
        if k in self.touched:
            return self.touched[k]
        v = dict.__getitem__(self, k)
        if self.r.layout(k, v) is not None:
            v = self.r.expand(k, v)
            self.touched[k] = v
        return v

    def __setitem__(self, k, v):
        dict.__setitem__(self, k, v)
        self.touched[k] = v
        self.assigned.add(k)

    def get(self, k, default=None):
        return self[k] if k in self else default

    def items(self):
        return [(k, self[k]) for k in self.keys()]

    def values(self):
        return [self[k] for k in self.keys()]

    def copy(self):
        return {k: self[k] for k in self.keys()}


def reason(model):
    '''Why the model cannot run as one row, or None if it can.'''
    s, p = model.s, model.p
    n = p['ny'] + 1
    if n < 2:
        return 'the grid has one row already'
    if p['boundary_lateral'] != 'circular':
        return 'the lateral boundary is %s, not circular' % p['boundary_lateral']
    if model.solver != 'steadystate':
        return 'the solver is %s; one row is implemented for steadystate' % model.solver
    if p['process_shear'] and p['method_shear'] not in ('fft_grid', 'spectral', '1Dstacks'):
        return 'method_shear = %s interpolates from a rotated grid, so its rows are not identical' % p['method_shear']
    if p['process_vegetation'] and p['method_vegetation'] != 'grass':
        return 'method_vegetation = %s spreads vegetation without a periodic lateral boundary' % p['method_vegetation']
    if p['process_vegetation'] and not p['grass_wrap_lateral']:
        return 'grass_wrap_lateral = F: the grass kernels stop at the lateral boundary'
    for key in ('process_fences', 'process_dune_erosion'):
        if p.get(key):
            return '%s = T is not checked for one row' % key
    for key in ('bedcomp_file', 'supply_file', 'threshold_file'):
        if isinstance(p.get(key), np.ndarray):
            return '%s is read during the run on the full grid' % key
    probe = OneRow(model)
    for k, v in s.items():
        if probe.layout(k, v, n) is not None and not probe.identical(k, v):
            return '%s differs between rows' % k
    hist = s.get('zb_hist')
    if isinstance(hist, list) and not all(_equal_rows(z, n) for _, z in hist):
        return 'the burial history differs between rows'
    return None


def start(model):
    '''Collapse the model to one row if ``transect_1d`` is set and the grid allows it.'''
    if not model.p.get('transect_1d'):
        return False
    why = reason(model)
    if why is not None:
        logger.info('transect_1d: computing all rows, because %s' % why)
        return False
    r = OneRow(model)
    r.collapse(model)
    logger.info('transect_1d: %d identical rows, computing one' % r.n)
    return r
