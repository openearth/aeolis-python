import numpy as np

# package modules
from utils import *


def initialize(s, p):
    '''Initialize bathymetry and bed composition

    Initialized bathymetry, computes cell sizes and orientation, bed
    layer thickness and bed composition.

    Parameters
    ----------
    s : dict
        Spatial grids
    p : dict
        Model configuration parameters

    Returns
    -------
    dict
        Spatial grids

    '''
    
    # get model dimensions
    ny = p['ny']
    nl = p['nlayers']
    nf = p['nfractions']

    # initialize x-dimension
    s['x'][:,:] = p['xgrid_file']
    s['ds'][:,1:] = np.diff(s['x'], axis=1)
    s['ds'][:,0] = s['ds'][:,1]

    # initialize y-dimension
    if ny == 0:
        s['y'][:,:] = 0.
        s['dn'][:,:] = 1.
        s['alfa'][:,:] = 0.
    else:
        s['y'][:,:] = p['ygrid_file']
        s['dn'][1:,:] = np.diff(s['y'], axis=0)
        s['dn'][0,:] = s['dn'][1,:]

        s['alfa'][1:-1,:] = np.arctan2(s['x'][2:,:] - s['x'][:-2,:],
                                       s['y'][2:,:] - s['y'][:-2,:])
        s['alfa'][0,:] = s['alfa'][1,:]
        s['alfa'][-1,:] = s['alfa'][-2,:]

    # compute cell areas
    s['dsdn'][:,:] = s['ds'] * s['dn']
    s['dsdni'][:,:] = 1. / s['dsdn']

    # initialize bathymetry
    s['zb'][:,:] = p['bed_file']

    # initialize bed layers
    s['thlyr'][:,:,:] = p['layer_thickness']

    # initialize bed composition
    if p['bedcomp_file'] is None:
        gs = makeiterable(p['grain_dist'])
        gs = gs / np.sum(gs)
        for i in range(nl):
            for j in range(nf):
                s['mass'][:,:,i,j] = p['rhop'] * p['porosity'] \
                                     * s['thlyr'][:,:,i] * gs[j]
    else:
        s['mass'][:,:,:,:] = p['bedcomp_file'].reshape(s['mass'].shape)                

    return s


def update(s, p):
    '''Update bathymetry and bed composition

    Update bed composition by moving sediment fractions between bed
    layers. The total mass in a single bed layer does not change as
    sediment removed from a layer is repleted with sediment from
    underlying layers. Similarly, excess sediment added in a layer is
    moved to underlying layers in order to keep the layer mass
    constant. The lowest bed layer exchanges sediment with an infinite
    sediment source that follows the original grain size distribution
    as defined in the model configuration file by ``grain_size`` and
    ``grain_dist``. The bathymetry is updated following the
    cummulative erosion/deposition over the fractions if ``bedupdate``
    is ``True``.

    Parameters
    ----------
    s : dict
        Spatial grids
    p : dict
        Model configuration parameters

    Returns
    -------
    dict
        Spatial grids

    '''

    nx = p['nx']
    ny = p['ny']
    nl = p['nlayers']
    nf = p['nfractions']

    # get total erosion
    pickup = s['pickup'].reshape((-1,nf))
    ero = np.sum(pickup, axis=1, keepdims=True).repeat(nf, axis=1)

    # get erosion and deposition cells
    ix_ero = ero[:,0] > 0.
    ix_dep = ero[:,0] < 0.

    # get sediment distribution in each cell
    m = s['mass'].reshape((-1,nl,nf))
    d = normalize(m, axis=2)

    # move mass among layers
    m[ix_ero,0,:] -= pickup[ix_ero,:]
    for i in range(1,nl):
        m[ix_ero,i-1,:] += ero[ix_ero,:] * d[ix_ero,i,:]
        m[ix_ero,i,  :] -= ero[ix_ero,:] * d[ix_ero,i,:]
        m[ix_dep,i-1,:] += ero[ix_dep,:] * d[ix_dep,i-1,:]
        m[ix_dep,i,  :] -= ero[ix_dep,:] * d[ix_dep,i-1,:]
    m[ix_dep,0,:] -= pickup[ix_dep,:]
    m[ix_dep,-1,:] += ero[ix_dep,:] * d[ix_dep,-1,:]
    m[ix_ero,-1,:] += ero[ix_ero,:] * makeiterable(p['grain_dist'])[np.newaxis,:].repeat(np.sum(ix_ero), axis=0)

    s['mass'] = m.reshape((ny+1,nx+1,nl,nf))

    # update bathy
    if p['bedupdate']:
        s['zb'] -= ero[:,0].reshape((ny+1,nx+1)) / (p['rhop'] * p['porosity'])

    return s


def mixtoplayer(s, p):
    '''Mix grain size distribution in top layer of the bed

    Simulates mixing of the top layers of the bed by wave action. The
    wave action is represented by a local wave height maximized by a
    maximum wave hieght over depth ratio ``gamma``. The mixing depth
    is a fraction of the local wave height indicated by
    ``facDOD``. The mixing depth is used to compute the number of bed
    layers that should be included in the mixing. The grain size
    distribution in these layers is then replaced by the average grain
    size distribution over these layers.

    Parameters
    ----------
    s : dict
        Spatial grids
    p : dict
        Model configuration parameters

    Returns
    -------
    dict
        Spatial grids

    '''

    if p['mixtoplayer']:
        
        # get model dimensions
        nl = p['nlayers']
        nf = p['nfractions']

        # compute depth of disturbence for each cell and repeat for each layer
        DOD = np.minimum(s['Hs'], (s['zs'] - s['zb']) * p['gamma']) * p['facDOD']
        DOD = DOD[:,:,np.newaxis].repeat(nl, axis=2)
        
        # determine what layers are above the depth of disturbance
        ix = (s['thlyr'].cumsum(axis=2) <= DOD) & (DOD > 0.)
        ix = ix[:,:,:,np.newaxis].repeat(nf, axis=3)
        
        # average mass over layers
        if np.any(ix):
            ix[:,:,0,:] = True # at least mix the top layer
            mass = s['mass'].copy()
            mass[~ix] = np.nan
            
            s['mass'][ix] = np.nanmean(mass, axis=2)[:,:,np.newaxis,:].repeat(nl, axis=2)[ix]
        
    return s
