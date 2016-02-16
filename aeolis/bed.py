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

    # determine net erosion
    pickup = s['pickup'].reshape((-1,nf))
    dm = np.sum(pickup, axis=1, keepdims=True).repeat(nf, axis=1)
    
    # get erosion and deposition cells
    ix_ero = dm[:,0] > 0.
    ix_dep = dm[:,0] < 0.
    
    # get sediment distribution in each cell
    m = s['mass'].reshape((-1,nl,nf))
    
    ## START: HANDLING OF HIGH EXPLOSIVE CELLS
    
    if True:
    
            # determine erosion and deposition fractions per cell
            ero =  np.maximum(0., pickup)
            dep = -np.minimum(0., pickup)

            # determine gross erosion
            erog = np.sum(ero, axis=1, keepdims=True).repeat(nf, axis=1)

            # determine net deposition cells with some erosional fractions
            ix = ix_dep & (erog[:,0] > 0)

            # determine distribution of deposition fractions
            if np.any(ix):
                d = normalize(dep, axis=1)
                ddep = erog[ix,:] * d[ix,:]
                pickup[ix,:] = -dep[ix,:] + ddep
                dm[ix,:] = np.sum(pickup[ix,:], axis=1, keepdims=True).repeat(nf, axis=1)
                m[ix,0,:] += -ero[ix,:] + ddep

    ## END
    
    
    # erode and deposit in echt cell (sierd has moved this step above the normalization.)
 #   m[:,0,:] -= pickup[:,:]
    # it can happen that pickup exceeds the available mass. Then set the mass to 0.    
    #m[m<0]=0
    #m[ix_dep,0,:] -= pickup[ix_dep,:] # here does the mass get negative
    
    # now detmine normalized mass (weighing factors)
    d = normalize(m, axis=2)
    m[:,0,:] -= pickup

    # move mass among layers
    
    if m.min()<0:
        print 'before pick up'
        print m.min()
   
    
    if m.min()<0:
        print 'after pick up'
        print m.min()
    
    for i in range(1,nl):
        m[ix_ero,i-1,:] += dm[ix_ero,:] * d[ix_ero,i,:]
        if (m.min()<0): #&((np.abs(m.min()))<1e-6) :
            m[m<0]=0
            print 'in loop sand1'
            print m.min()
            print i

        m[ix_ero,i,  :] -= dm[ix_ero,:] * d[ix_ero,i,:]
        if (m.min()<0): #&((np.abs(m.min()))<1e-6) :
            print 'in loop sand2'
            m[m<0]=0
            print m.min()
            print i
            
            
            
        m[ix_dep,i-1,:] += dm[ix_dep,:] * d[ix_dep,i-1,:]
        if m.min()<0:
            print 'in loop sand3'
            print m.min()
            print i

        m[ix_dep,i,  :] -= dm[ix_dep,:] * d[ix_dep,i-1,:]
      

    
    if m.min()<0:
        print 'after loop sand'
        print m.min()
    
    
    if m.min()<0:
        print 'after loop sand2'
        print m.min()
        
    m[ix_dep,-1,:] += dm[ix_dep,:] * d[ix_dep,-1,:]
    
    if m.min()<0:
        print 'after loop sand3'
        print m.min()
        
    m[ix_ero,-1,:] += dm[ix_ero,:] * makeiterable(p['grain_dist'])[np.newaxis,:].repeat(np.sum(ix_ero), axis=0)
    
    
        
    if m.min()<0:
        print 'after moving sand'
        print m.min()

    s['mass'] = m.reshape((ny+1,nx+1,nl,nf))

    # update bathy
    if p['bedupdate']:
        s['zb'] -= dm[:,0].reshape((ny+1,nx+1)) / (p['rhop'] * p['porosity'])

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
