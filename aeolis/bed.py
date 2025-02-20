'''This file is part of AeoLiS.
   
AeoLiS is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
   
AeoLiS is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
   
You should have received a copy of the GNU General Public License
along with AeoLiS.  If not, see <http://www.gnu.org/licenses/>.
   
AeoLiS  Copyright (C) 2015 Bas Hoonhout

bas.hoonhout@deltares.nl         b.m.hoonhout@tudelft.nl
Deltares                         Delft University of Technology
Unit of Hydraulic Engineering    Faculty of Civil Engineering and Geosciences
Boussinesqweg 1                  Stevinweg 1
2629 HVDelft                     2628CN Delft
The Netherlands                  The Netherlands

'''


from __future__ import absolute_import, division

import logging
import numpy as np
import aeolis.gridparams
from matplotlib import pyplot as plt
from numba import njit

# package modules
from aeolis.utils import *
#import matplotlib.pyplot as plt


# initialize logger
logger = logging.getLogger(__name__)


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
    nx = p['nx']
    nl = p['nlayers']
    nf = p['nfractions']

    # initialize bathymetry
    s['zb'][:,:] = p['bed_file']
    s['zb0'][:,:] = p['bed_file']
    s['zne'][:,:] = p['ne_file']

    #initialize thickness of erodable or dry top layer
    s['zdry'][:,:] = 0.05
    
    # initialize bed layers
    s['thlyr'][:,:,:] = p['layer_thickness']

    # initialize bed composition
    if isinstance(p['grain_dist'], str):
            logger.log_and_raise('Grain size file not recognized as array, check file path and whether all values have been filled in.', exc=ValueError) 

    if p['bedcomp_file'] is None and p['grain_dist'].ndim == 1 and p['grain_dist'].dtype == 'float64' or p['grain_dist'].dtype == 'int': 
        # Both float and int are included as options for the grain dist to make sure there is no error when grain_dist is filled in as 1 instead of 1.0. 
        for i in range(nl):
            gs = makeiterable(p['grain_dist'])
            gs = gs / np.sum(gs)
            for j in range(nf):
                s['mass'][:,:,i,j] = p['rhog'] * (1. - p['porosity']) \
                                     * s['thlyr'][:,:,i] * gs[j]
    elif p['bedcomp_file'] is None and p['grain_dist'].ndim > 1: #allows simple cases with layering, txt file containing distribution per fraction per column and layers in the rows.
        if nl != p['grain_dist'].shape[0]:
            logger.log_and_raise('Grain size distribution not assigned for each layer, not enough rows for the number of layers', exc=ValueError)
        for i in range(nl):
            gs = makeiterable(p['grain_dist'][i,:])
            gs = gs / np.sum(gs)
            for j in range(nf):
                s['mass'][:,:,i,j] = p['rhog'] * (1. - p['porosity']) \
                                     * s['thlyr'][:,:,i] * gs[j]
    else:
        s['mass'][:,:,:,:] = p['bedcomp_file'].reshape(s['mass'].shape)                

    # initialize masks
    for k, v in p.items():
        if k.endswith('_mask'):
            if v is None:
                s[k] = 1.
            else:
                s[k] = v.reshape(s['zb'].shape)

    # initialize threshold
    if p['threshold_file'] is not None:
        s['uth'] = p['threshold_file'][:,:,np.newaxis].repeat(nf, axis=-1)

    # initialize sand and cobble layers for composite beaches
    s['zsand'][:,:] = p['zsand_file']
    s['dcob'][:,:] = p['dcob_file']
    s['dsandcover'][:,:] = p['dsandcover_file']
    s['doverlap'][:,:] = p['doverlap_file']
    if p['process_bedupdate_comp'] is True and (p['zsand_file'] is None or p['dcob_file'] is None or p['dsandcover_file'] is None or p['doverlap_file'] is None):
        logger.log_and_raise('Process bedupdate for composite beaches is turned on but no initial sand and cobble locations are provided. Please'
                             'provide a zsand_file, dcob_file, dsandcover_file and doverlap_file', exc=ValueError) 

    return s


def mixtoplayer(s, p):
    '''Mix grain size distribution in top layer of the bed.

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

    if p['process_mixtoplayer']:
        
        # get model dimensions
        nx = p['nx']+1
        ny = p['ny']+1
        nl = p['nlayers']
        nf = p['nfractions']

        # compute depth of disturbance for each cell and repeat for each layer
        DOD = p['facDOD'] * s['Hsmix']

        # compute ratio total layer thickness and depth of disturbance 
        ix = DOD > 0.
        f = np.ones(DOD.shape)
        f[ix] = np.minimum(1., s['thlyr'].sum(axis=2)[ix] / DOD[ix])

        # correct shapes
        DOD = DOD[:,:,np.newaxis].repeat(nl, axis=2)
        f = f[:,:,np.newaxis].repeat(nl, axis=2)

        # determine what layers are above the depth of disturbance
        ix = (s['thlyr'].cumsum(axis=2) <= DOD) & (DOD > 0.)
        ix = ix[:,:,:,np.newaxis].repeat(nf, axis=3)
        f = f[:,:,:,np.newaxis].repeat(nf, axis=3)
        
        # average mass over layers
        if np.any(ix):
            ix[:,:,0,:] = True # at least mix the top layer
            mass = s['mass'].copy()
            mass[~ix] = np.nan
            
            # gd = normalize(p['grain_dist']) * p['rhog'] * (1. - p['porosity'])
            # gd = gd.reshape((1,1,1,-1)).repeat(ny, axis=0) \
                                       # .repeat(nx, axis=1) \
                                       # .repeat(nl, axis=2)

            mass1 = np.nanmean(mass, axis=2, keepdims=True).repeat(nl, axis=2)
            # mass2 = gd * s['thlyr'][:,:,:,np.newaxis].repeat(nf, axis=-1)
            mass = mass1 * f + mass * (1. - f)

        
            s['mass'][ix] = mass[ix]
            
    return s


def wet_bed_reset(s, p):
    ''' Text



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

    if p['process_wet_bed_reset']:
        
        Tbedreset = p['dt_opt'] / p['Tbedreset']
        
        ix = s['zs'] > (s['zb'] + 0.01)
        s['zb'][ix] += (s['zb0'][ix] - s['zb'][ix]) * Tbedreset
            
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

    # determine total mass that should be exchanged between layers
    dm = -np.sum(pickup, axis=-1, keepdims=True).repeat(nf, axis=-1)
    
    # get erosion and deposition cells
    ix_ero = dm[:,0] < 0.
    ix_dep = dm[:,0] > 0.
    
    # reshape mass matrix
    m = s['mass'].reshape((-1,nl,nf))

    # negative mass may occur in case of deposition due to numerics,
    # which should be prevented
    m, dm, pickup = prevent_negative_mass(m, dm, pickup)
    
    # determine weighing factors
    d = normalize(m, axis=2)
    
    # move mass among layers
    m[:,0,:] -= pickup
    m = arrange_layers(m,dm,d,nl,ix_ero,ix_dep)
    
    # this is replaced by arrange_layers and speed up using numba
    # for i in range(1,nl):
    #     m[ix_ero,i-1,:] -= dm[ix_ero,:] * d[ix_ero,i,:]
    #     m[ix_ero,i,  :] += dm[ix_ero,:] * d[ix_ero,i,:]
    #     m[ix_dep,i-1,:] -= dm[ix_dep,:] * d[ix_dep,i-1,:]
    #     m[ix_dep,i,  :] += dm[ix_dep,:] * d[ix_dep,i-1,:]
    #m[ix_dep,-1,:] -= dm[ix_dep,:] * d[ix_dep,-1,:]

    if p['grain_dist'].ndim == 2: 
        m[ix_ero,-1,:] -= dm[ix_ero,:] * normalize(p['grain_dist'][-1,:])[np.newaxis,:].repeat(np.sum(ix_ero), axis=0)
    elif type(p['bedcomp_file']) == np.ndarray:
        gs = p['bedcomp_file'].reshape((-1,nl,nf))
        m[ix_ero,-1,:] -= dm[ix_ero,:] * normalize(gs[ix_ero,-1, :], axis=1)
    else:
        m[ix_ero,-1,:] -= dm[ix_ero,:] * normalize(p['grain_dist'])[np.newaxis,:].repeat(np.sum(ix_ero), axis=0)
    # remove tiny negatives
    m = prevent_tiny_negatives(m, p['max_error'])

    # warn if not all negatives are gone
    if m.min() < 0:
        logger.warning(format_log('Negative mass',
                                  nrcells=np.sum(np.any(m<0., axis=-1)),
                                  minvalue=m.min(),
                                  minwind=s['uw'].min(),
                                  time=p['_time']))
    
    # reshape mass matrix
    s['mass'] = m.reshape((ny+1,nx+1,nl,nf))

    # update bathy
    if p['process_bedupdate']:

        dz = dm[:, 0].reshape((ny + 1, nx + 1)) / (p['rhog'] * (1. - p['porosity']))
        # s['dzb'] = dm[:, 0].reshape((ny + 1, nx + 1))
        s['dzb'] = dz.copy()

        if p['process_bedupdate_comp']: 
            s = update_composite(s, p)
        else:
            s['zb'] += dz

        # redistribute sediment from inactive zone to marine interaction zone
        if p['process_tide']:
            s['zs'] += dz #???
    
    return s


def prevent_negative_mass(m, dm, pickup):
    '''Handle situations in which negative mass may occur due to numerics

    Negative mass may occur by moving sediment to lower layers down to
    accomodate deposition of sediments. In particular two cases are
    important:

    #. A net deposition cell has some erosional fractions.

       In this case the top layer mass is reduced according to the
       existing sediment distribution in the layer to accomodate
       deposition of fresh sediment. If the erosional fraction is
       subtracted afterwards, negative values may occur. Therefore the
       erosional fractions are subtracted from the top layer
       beforehand in this function. An equal mass of deposition
       fractions is added to the top layer in order to keep the total
       layer mass constant. Subsequently, the distribution of the
       sediment to be moved to lower layers is determined and the
       remaining deposits are accomodated.

    #. Deposition is larger than the total mass in a layer.

       In this case a non-uniform distribution in the bed may also
       lead to negative values as the abundant fractions are reduced
       disproportionally as sediment is moved to lower layers to
       accomodate the deposits. This function fills the top layers
       entirely with fresh deposits and moves the existing sediment
       down such that the remaining deposits have a total mass less
       than the total bed layer mass. Only the remaining deposits are
       fed to the routine that moves sediment through the layers.

    Parameters
    ----------
    m : np.ndarray
        Sediment mass in bed (nx*ny, nl, nf)
    dm : np.ndarray
        Total sediment mass exchanged between layers (nx*ny, nf)
    pickup : np.ndarray
        Sediment pickup (nx*ny, nf)

    Returns
    -------
    np.ndarray
        Sediment mass in bed (nx*ny, nl, nf)
    np.ndarray
        Total sediment mass exchanged between layers (nx*ny, nf)
    np.ndarray
        Sediment pickup (nx*ny, nf)

    Note
    ----
    The situations handled in this function can also be prevented by
    reducing the time step, increasing the layer mass or increasing
    the adaptation time scale.

    '''

    nl = m.shape[1]
    nf = m.shape[2]

    ###
    ### case #1: deposition cells with some erosional fractions
    ###
    
    ix_dep = dm[:,0] > 0.
    
    # determine erosion and deposition fractions per cell
    ero =  np.maximum(0., pickup)
    dep = -np.minimum(0., pickup)

    # determine gross erosion
    erog = np.sum(ero, axis=1, keepdims=True).repeat(nf, axis=1)

    # determine net deposition cells with some erosional fractions
    ix = ix_dep & (erog[:,0] > 0)

    # remove erosional fractions from pickup and remove an equal mass
    # of accretive fractions from the pickup, adapt sediment exchange
    # mass and bed composition accordingly
    if np.any(ix):
        d = normalize(dep, axis=1)
        ddep = erog[ix,:] * d[ix,:]
        pickup[ix,:] = -dep[ix,:] + ddep
        dm[ix,:] = -np.sum(pickup[ix,:], axis=-1, keepdims=True).repeat(nf, axis=-1)
        m[ix,0,:] -= ero[ix,:] - ddep # FIXME: do not use deposition in normalization

    ###
    ### case #2: deposition cells with deposition larger than the mass present in the top layer
    ###

    mx = m[:,0,:].sum(axis=-1, keepdims=True)

    # determine deposition in terms of layer mass (round down)
    n = dm[:,:1] // mx

    # determine if deposition is larger than a sinle layer mass
    if np.any(n > 0):

        # determine distribution of deposition
        d = normalize(pickup, axis=1)

        # walk through layers from top to bottom
        for i in range(nl):

            ix = (n > i).flatten()
            if not np.any(ix):
                break

            # move all sediment below current layer down one layer
            m[ix,(i+1):,:] = m[ix,i:-1,:]

            # fill current layer with deposited sediment
            m[ix,i,:] = mx[ix,:].repeat(nf, axis=1) * d[ix,:]

            # remove deposited sediment from pickup
            pickup[ix,:] -= m[ix,i,:]

        # discard any remaining deposits at locations where all layers
        # are filled with fresh deposits
        ix = (dm[:,:1] > mx).flatten()
        if np.any(ix):
            pickup[ix,:] = 0.

        # recompute sediment exchange mass
        dm[ix,:] = -np.sum(pickup[ix,:], axis=-1, keepdims=True).repeat(nf, axis=-1)

    return m, dm, pickup


def average_change(l, s, p):

    #Compute bed level change with previous time step [m/timestep]
    s['dzb'] = s['zb'] - l['zb']

    # Collect time steps
    s['dzbyear'] = s['dzb'] * (3600. * 24. * 365.25) / (p['dt_opt'] * p['accfac'])
    n = (p['dt_opt'] * p['accfac']) / p['avg_time']
    s['dzbavg'] = n*s['dzbyear']+(1-n)*l['dzbavg']
    
    # Calculate average bed level change as input for vegetation growth [m/year]
    # s['dzbveg'] = s['dzbavg'].copy()
    s['dzbveg'] = s['dzbyear'].copy()
    
    
    return s

@njit
def arrange_layers(m,dm,d,nl,ix_ero,ix_dep):
    '''Arranges mass redistrubution between layers. 
    This function is called in the bed.update fucntion to speed up code using numba
    
    

    Parameters
    ----------
    m       :   array
                mass in layers
    dm      :   array
                total mass exchanged between layers derrived from pickup
    d       :   array
                normalized mass in layers
    nl      :   int
                number of layers
    ix_dep  :   array
                cells for deposition
    ix_ero  :   array
                cells for erosion

   
    Returns
    -------
    m
     
    '''
    for i in range(1,nl):
        m[ix_ero,i-1,:] -= dm[ix_ero,:] * d[ix_ero,i,:]
        m[ix_ero,i,  :] += dm[ix_ero,:] * d[ix_ero,i,:]
        m[ix_dep,i-1,:] -= dm[ix_dep,:] * d[ix_dep,i-1,:]
        m[ix_dep,i,  :] += dm[ix_dep,:] * d[ix_dep,i-1,:]
    m[ix_dep,-1,:] -= dm[ix_dep,:] * d[ix_dep,-1,:]

    return m

def update_composite(s, p):
    '''Update bed and sand level and cobble infilling for composite beaches

    Update XX by YYY
    layers. 
    Initial sand level, cobble thickness, overlap and sandcover are defined 
    in the model configuration file by ``zsand``, ``dcob``, ``doverlap`` and
    ``dsandcover``. The bathymetry is updated if the sand cover increases.

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

    dz = s['dzb']
    por = p['porosity'] # Assumes that cobble porosity is the same as sand porosity

    # FIND LOCATIONS
    # with deposition
    ix_depo = dz >0
    # with erosion
    ix_ero = dz < 0
    # with 0 bed level change
    ix_none = dz == 0.

    # where no cobble is present, i.e. pure sand
    ix_nocob = s['dcob'] == 0.
    # where open cobbles are present
    ix_cob = s['dcob'] > 0.
    # where sand cover is present
    ix_sc = s['dsandcover'] > 0.
    # where no sand cover is present
    ix_nosc = s['dsandcover'] == 0.
    # where elevation change is bigger than sand cover, these are locations that can accommodate erosion
    ix_sc_accom = dz + s['dsandcover'] >= 0.
    # where elevation change is smaller than sand cover, these are locations that cannot accommodate erosion
    ix_sc_noaccom = dz + s['dsandcover'] < 0.
    # where open cobbles are present at the top
    ix_opencob = s['dcob']-s['doverlap'] > 0.
    # where elevation change is smaller than pore space in open cobble, these are locations that can accommodate deposition
    ix_oc_accom = (s['dcob']-s['doverlap']) - dz/(1-por) >= 0.
    # where elevation change is bigger than pore space in open cobble, these are locations that can't accommodate deposition
    ix_oc_noaccom = (s['dcob']-s['doverlap']) - dz/(1-por) < 0.
    # where cobbles are fully filled with sand, so doverlap == dcob
    ix_fullcob = s['dcob']==s['doverlap']
    # where filled cobbles are present at the top, dcob > 0 & doverlap == dcob & dsandcover == 0 
    ix_fillcob = (ix_cob & ix_fullcob & ix_nosc)

    # where no bed level change occurs and cobble is present
    ix_step0 = ix_none & ix_cob
    # where deposition takes place and sand cover is present
    ix_step1 = ix_depo & ix_sc
    # where deposition takes place and filled cobble is present  
    ix_step2 = ix_depo & ix_cob & ix_fillcob
    # where deposition takes place and open cobble is present with enough accommodation space
    ix_step3 = ix_depo & ix_cob & ix_opencob & ix_oc_accom
    # where deposition takes place and open cobble is present without enough accommodation space
    ix_step4 = ix_depo & ix_cob & ix_opencob & ix_oc_noaccom
    # where erosion takes place and open cobble is present
    ix_step5 = ix_ero & ix_cob & ix_opencob
    # where erosion takes place and filled cobble is present
    ix_step6 = ix_ero & ix_cob & ix_fillcob
    # where erosion takes place, sand cover is present and there is enough accommodation space
    ix_step7 = ix_ero & ix_sc & ix_sc_accom
    # where erosion takes place, sand cover is present and there is not enough accommodation space
    ix_step8 = ix_ero & ix_sc & ix_sc_noaccom

    # Check if all locations are assigned to a single step
    ix_steps = np.vstack((ix_nocob[1,:], ix_step0[1,:], ix_step1[1,:], ix_step2[1,:], ix_step3[1,:], ix_step4[1,:], ix_step5[1,:], ix_step6[1,:], ix_step7[1,:], ix_step8[1,:]))
    check = np.sum(ix_steps, axis=0)
    if np.sum(s['dcob'] - s['doverlap'] < 0) > 0:
        print('Cobble thickness is smaller than overlap thickness')
    if np.sum(check!=1) > 0:
        args = np.argwhere(check!=1)
        nocob_check = ix_nocob[1,args]
        step0_check = ix_step0[1,args]
        step1_check = ix_step1[1,args]
        step2_check = ix_step2[1,args]
        step3_check = ix_step3[1,args]
        step4_check = ix_step4[1,args]
        step5_check = ix_step5[1,args]
        step6_check = ix_step6[1,args]
        step7_check = ix_step7[1,args]
        step8_check = ix_step8[1,args]
        print(check)

    # if no cobble is present, 
    # change sand level. This is independent of erosion or deposition
    s['zsand'][ix_nocob] += dz[ix_nocob]

    # if no bed level change occurs and cobble is present,
    # keep all variables the same.
    s['zsand'][ix_step0] = s['zsand'][ix_step0]

    # if deposition takes place and sand cover is present 
    # change sand cover and sand level. 
    s['dsandcover'][ix_step1] += dz[ix_step1]
    s['zsand'][ix_step1] += dz[ix_step1]

    # if deposition takes place and filled cobble is present 
    # change sand cover and sand level. 
    s['dsandcover'][ix_step2] += dz[ix_step2]
    s['zsand'][ix_step2] += dz[ix_step2]

    # if deposition takes place and open cobble is present with enough accommodation space, 
    # change sand level and overlap. 
    s['zsand'][ix_step3] += dz[ix_step3]/(1-por)
    s['doverlap'][ix_step3] += dz[ix_step3]/(1-por)

    if np.sum(s['dcob'] - s['doverlap'] < 0) > 0:
        print('Cobble thickness is smaller than overlap thickness')

    # if deposition takes place and open cobble is present without enough accommodation space, 
    # change sand cover, sand level and overlap. 
    s['dsandcover'][ix_step4] += dz[ix_step4] - (s['dcob'][ix_step4] - s['doverlap'][ix_step4])*(1-por)
    s['zsand'][ix_step4] += (s['dcob'][ix_step4] - s['doverlap'][ix_step4]) + s['dsandcover'][ix_step4]
    s['doverlap'][ix_step4] = s['dcob'][ix_step4]

    if np.sum(s['dcob'] - s['doverlap'] < 0) > 0:
        np.argwhere(s['dcob'] - s['doverlap'] < 0)
        print('Cobble thickness is smaller than overlap thickness')

    # if erosion takes place and open cobble is present, 
    # change sand level and overlap. 
    s['zsand'][ix_step5] += dz[ix_step5]/(1-por)
    s['doverlap'][ix_step5] += dz[ix_step5]/(1-por)

    if np.sum(s['dcob'] - s['doverlap'] < 0) > 0:
        print('Cobble thickness is smaller than overlap thickness')

    # if erosion takes place and filled cobble is present, 
    # change sand level and overlap. 
    s['zsand'][ix_step6] += dz[ix_step6]/(1-por)
    s['doverlap'][ix_step6] += dz[ix_step6]/(1-por)

    # if erosion takes place, sand cover is present and there is enough accommodation space, 
    # change sand cover and sand level. 
    s['dsandcover'][ix_step7] += dz[ix_step7]
    s['zsand'][ix_step7] += dz[ix_step7]

    # print(s['dcob'][0,380], s['doverlap'][0,380], s['dsandcover'][0,380], s['zsand'][0,380], dz[0,380])

    # if erosion takes place, sand cover is present and there is not enough accommodation space, 
    # change overlap, sand level and sand cover. 
    s['doverlap'][ix_step8] += (dz[ix_step8] + s['dsandcover'][ix_step8])/(1-por)
    s['zsand'][ix_step8] += (dz[ix_step8] + s['dsandcover'][ix_step8])/(1-por) - s['dsandcover'][ix_step8]
    s['dsandcover'][ix_step8] = 0

    # print(s['dcob'][0,380], s['doverlap'][0,380], s['dsandcover'][0,380], s['zsand'][0,380], dz[0,380])

    if np.sum(s['dcob'] - s['doverlap'] < 0) > 0:
        arg = np.argwhere(s['dcob'] - s['doverlap'] < 0)
        print('Cobble thickness is smaller than overlap thickness')
    # for now, assume not enough erosion takes place to halt erosion through changing the roughness. Requires small enough time step.

    # For all locations calculate the new bed level
    s['zb'] = s['zsand'] + (s['dcob']-s['doverlap'])

    return s

