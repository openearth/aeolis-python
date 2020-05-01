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
import matplotlib.pyplot as plt

# package modules
from aeolis.utils import *


# initialize logger
logger = logging.getLogger(__name__)

def saltationvelocity(s, p):
    '''Define saltation velocity u [m/s]

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
    
    # Initialize arrays
    
    nf = p['nfractions']
    
    uw  = s['uw']
    uws = s['uws']
    uwn = s['uwn']
        
    ustar  = s['ustar']
    ustars = s['ustars']
    ustarn = s['ustarn']
    
    us = np.zeros(ustar.shape)
    un = np.zeros(ustar.shape)
    u  = np.zeros(ustar.shape)

    # u with direction of perturbation theory
    
    ix = ustar != 0
    
    us[ix] += 1 * ustars[ix] / ustar[ix]
    un[ix] += 1 * ustarn[ix] / ustar[ix]
            
    u[ix] = np.hypot(us[ix], un[ix])
    
    # u under the sep bubble
    
    ix = ustar == 0
    
    us[ix] += 0.2 * uws[ix] / uw[ix]
    un[ix] += 0.2 * uwn[ix] / uw[ix]
            
    u[ix] = np.hypot(us[ix], un[ix])
                       
    s['us'] = us[:,:,np.newaxis].repeat(nf, axis=2)
    s['un'] = un[:,:,np.newaxis].repeat(nf, axis=2)
    s['u']  = u[:,:,np.newaxis].repeat(nf, axis=2)
        
    return s


def equilibrium(s, p):
    '''Compute equilibrium sediment concentration following Bagnold (1937)

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

    if p['process_transport']:
                
        nf = p['nfractions']
        
        us = s['us']
        un = s['un']   
        u  = s['u']
        
        ustar  = s['ustar'][:,:,np.newaxis].repeat(nf, axis=2)
        uth    = s['uth']
        uthf   = s['uthf']

        rhoa   = p['rhoa'] 
        g      = p['g']
        
        s['Cu']  = np.zeros(uth.shape)
        s['Cuf'] = np.zeros(uth.shape)
    
                
        ix = (ustar != 0.)*(u != 0.)
        
        if p['method_transport'].lower() == 'bagnold':
            s['Cu'][ix]  = np.maximum(0., p['Cb'] * rhoa / g * (ustar[ix] - uth[ix])**3 / u[ix])
            s['Cuf'][ix] = np.maximum(0., p['Cb'] * rhoa / g * (ustar[ix] - uthf[ix])**3 / u[ix])
        
        elif p['method_transport'].lower() == 'kawamura':
            s['Cu'][ix]  = np.maximum(0., p['Ck'] * rhoa / g * (ustar[ix] - uth[ix])**2 * (ustar[ix] + uth[ix]) / u[ix])
            s['Cuf'][ix] = np.maximum(0, p['Ck'] * rhoa / g * (ustar[ix] - uthf[ix])**2 * (ustar[ix] + uthf[ix]) / u[ix])
        
        elif p['method_transport'].lower() == 'lettau':
            s['Cu'][ix]  = np.maximum(0., p['Cl'] * rhoa / g * ustar[ix]**2 * (ustar[ix] - uth[ix]) / u[ix])
            s['Cuf'][ix] = np.maximum(0., p['Cl'] * rhoa / g * ustar[ix]**2 * (ustar[ix] - uthf[ix]) / u[ix])

        elif p['method_transport'].lower() == 'dk':
            s['Cu'][ix]  = np.maximum(0., p['Cdk'] * rhoa / g * uth[ix] * (ustar[ix]**2 - uth[ix]**2) / u[ix])
            s['Cuf'][ix] = np.maximum(0., p['Cdk'] * rhoa / g * uthf[ix] * (ustar[ix]**2 - uthf[ix]**2) / u[ix])
        
        else:
            logger.log_and_raise('Unknown transport formulation [%s]' % method, exc=ValueError)   
                                       
    s['Cu']  *= p['accfac']
    s['Cuf'] *= p['accfac']
    
    return s


def compute_weights(s, p):
    '''Compute weights for sediment fractions

    Multi-fraction sediment transport needs to weigh the transport of
    each sediment fraction to prevent the sediment transport to
    increase with an increasing number of sediment fractions. The
    weighing is not uniform over all sediment fractions, but depends
    on the sediment availibility in the air and the bed and the bed
    interaction parameter ``bi``.

    Parameters
    ----------
    s : dict
        Spatial grids
    p : dict
        Model configuration parameters

    Returns
    -------
    numpy.ndarray
        Array with weights for each sediment fraction

    '''

    w_air = normalize(s['Ct'], s['Cu'])
    w_bed = normalize(s['mass'][:,:,0,:], axis=2)

    w = (1. - p['bi']) * w_air \
        + (1. - np.minimum(1., (1. - p['bi']) * np.sum(w_air, axis=2, keepdims=True))) * w_bed
    w = normalize(w, axis=2)
    
    return w, w_air, w_bed


def renormalize_weights(w, ix):
    '''Renormalizes weights for sediment fractions

    Renormalizes weights for sediment fractions such that the sum of
    all weights is unity. To ensure that the erosion of specific
    fractions does not exceed the sediment availibility in the bed,
    the normalization only modifies the weights with index equal or
    larger than ``ix``.

    Parameters
    ----------
    w : numpy.ndarray
        Array with weights for each sediment fraction
    ix : int
        Minimum index to be modified

    Returns
    -------
    numpy.ndarray
        Array with weights for each sediment fraction

    '''
    
    f = np.sum(w[:,:,:ix], axis=2, keepdims=True)
    w[:,:,ix:] = normalize(w[:,:,ix:], axis=2) * (1. - f)

    # normalize in case of supply-limitation
    # use uniform distribution in case of no supply
    w = normalize(w, axis=2, fill=1./w.shape[2])

    return w
