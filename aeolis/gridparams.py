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

# package modules
from aeolis.utils import *

# initialize logger
logger = logging.getLogger(__name__)

def initialize(s, p):
    '''EXPLAIN WHAT HAPPENS IN THIS FUNCTION?
    
    
    
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
    
    
    # initialize x-dimensions
    s['x'][:,:] = p['xgrid_file']
    
    # # Initializing all other arrays
    # s['xz'] = np.zeros(np.shape(s['x']))
    # s['xu'] = np.zeros(np.shape(s['x']))
    # s['xv'] = np.zeros(np.shape(s['x']))
    # s['xc'] = np.zeros(np.shape(s['x']))
    
    # s['yz'] = np.zeros(np.shape(s['x']))
    # s['yu'] = np.zeros(np.shape(s['x']))
    # s['yv'] = np.zeros(np.shape(s['x']))
    # s['yc'] = np.zeros(np.shape(s['x']))
    
    # s['dsz'] = np.zeros(np.shape(s['x']))
    # s['dsu'] = np.zeros(np.shape(s['x']))
    # s['dsv'] = np.zeros(np.shape(s['x']))
    # s['dsc'] = np.zeros(np.shape(s['x']))
    
    # s['dnz'] = np.zeros(np.shape(s['x']))
    # s['dnu'] = np.zeros(np.shape(s['x']))
    # s['dnv'] = np.zeros(np.shape(s['x']))
    # s['dnc'] = np.zeros(np.shape(s['x']))

    # s['dsdnz'] = np.zeros(np.shape(s['x']))
    # s['dsdnzi'] = np.zeros(np.shape(s['x']))
    
    # s['alfaz'] = np.zeros(np.shape(s['x']))
    # s['alfau'] = np.zeros(np.shape(s['x']))
    # s['alfav'] = np.zeros(np.shape(s['x']))
    
    # # World coordinates of z-points
    # s['xz'][:,:] = s['x'][:,:]
    
    # # World coordinates of u-points
    # s['xu'][:,1:] = 0.5 * (s['xz'][:,:-1] + s['xz'][:,1:])
    # s['xu'][:,0]  = 1.5 *  s['xz'][:,0] - 0.5 * s['xz'][:,1]
    
    # # World coordinates of v-points
    # s['xv'][1:,:] = 0.5 * (s['xz'][:-1,:] + s['xz'][1:,:])
    # s['xv'][0,:]  = 1.5 *  s['xz'][0,:]   - 0.5 * s['xz'][1,:]
    
    # # World coordinates of c-points
    # s['xc'][1:,1:] = 0.25 *(s['xz'][:-1,:-1] + s['xz'][:-1,1:] + s['xz'][1:,:-1] + s['xz'][1:,1:])
    # s['xc'][1:,0]  = 0.5 * (s['xu'][:-1,0]   + s['xu'][1:,0])
    # s['xc'][0,1:]  = 0.5 * (s['xv'][0,:-1]   + s['xv'][0,1:])
    # s['xc'][0,0]   = s['xu'][0,0]
    
    # # initialize y-dimension
    ny = p['ny']
    
    if ny == 0:
        s['y'][:,:] = 0.
    #     s['yz'][:,:] = 0.
    #     s['yu'][:,:] = 0.
    #     s['yv'][:,:] = 0.
    #     s['dnz'][:,:] = 1.
    #     s['dnu'][:,:] = 1.
    #     s['dnv'][:,:] = 1.
    #     s['dnc'][:,:] = 1.
    #     s['alfaz'][:,:] = 0.
    else:
        # initialize y-dimensions
        s['y'][:,:] = p['ygrid_file']
        
    #     # World coordinates of z-points
    #     s['yz'][:,:] = s['y'][:,:] # Different from XBeach
        
    #     # World coordinates of u-points
    #     s['yu'][:,1:] = 0.5 * (s['yz'][:,:-1] + s['yz'][:,1:])
    #     s['yu'][:,0]  = 1.5 *  s['yz'][:,0]   - 0.5 * s['yz'][:,1]
        
    #     # World coordinates of v-points
    #     s['yv'][1:,:] = 0.5 * (s['yz'][:-1,:] + s['yz'][1:,:])
    #     s['yv'][0,:]  = 1.5 *  s['yz'][0,:]   - 0.5 * s['yz'][1,:]
        
    #     # World coordinates of c-points
    #     s['yc'][1:,1:] = 0.25 *(s['yz'][:-1,:-1] + s['yz'][:-1,1:] + s['yz'][1:,:-1] + s['yz'][1:,1:])
    #     s['yc'][0,1:]  = 0.5 * (s['yv'][0,:-1]  + s['yv'][0,1:])
    #     s['yc'][1:,0]  = 0.5 * (s['yu'][:-1,0]  + s['yu'][1:,0])
    #     s['yc'][0,0]   = s['yv'][0,0]
        
    #     # Distances in n-direction
    #     s['dnz'][:-1,:] = ((s['yv'][:-1,:]-s['yv'][1:,:])**2.+(s['xv'][:-1,:]-s['xv'][1:,:])**2.)**0.5
    #     s['dnu'][1:,:] = ((s['xc'][:-1,:]-s['xc'][1:,:])**2.+(s['yc'][:-1,:]-s['yc'][1:,:])**2.)**0.5
    #     s['dnv'][1:,:] = ((s['xz'][:-1,:]-s['xz'][1:,:])**2.+(s['yz'][:-1,:]-s['yz'][1:,:])**2.)**0.5
    #     s['dnc'][1:,:] = ((s['xu'][:-1,:]-s['xu'][1:,:])**2.+(s['yu'][:-1,:]-s['yu'][1:,:])**2.)**0.5
        
    #     s['dnz'][-1,:] = s['dnz'][-2,:] 
    #     s['dnu'][0,:] = s['dnu'][1,:]
    #     s['dnv'][0,:] = s['dnv'][1,:]
    #     s['dnc'][0,:] = s['dnc'][1,:]
    
    # # Distances in s-direction
    # s['dsz'][:,:-1] = ((s['xu'][:,:-1]-s['xu'][:,1:])**2.+(s['yu'][:,:-1]-s['yu'][:,1:])**2.)**0.5
    # s['dsu'][:,1:] = ((s['xz'][:,:-1]-s['xz'][:,1:])**2.+(s['yz'][:,:-1]-s['yz'][:,1:])**2.)**0.5
    # s['dsv'][:,1:] = ((s['xc'][:,:-1]-s['xc'][:,1:])**2.+(s['yc'][:,:-1]-s['yc'][:,1:])**2.)**0.5
    # s['dsc'][:,1:] = ((s['xv'][:,:-1]-s['xv'][:,1:])**2.+(s['yv'][:,:-1]-s['yv'][:,1:])**2.)**0.5
    
    # s['dsz'][:,-1] = s['dsz'][:,-2] 
    # s['dsu'][:,0] = s['dsu'][:,1]
    # s['dsv'][:,0] = s['dsv'][:,1]
    # s['dsc'][:,0] = s['dsc'][:,1]

    # # Cell areas
    # s['dsdnz'][:-1,:-1] = (0.5*(s['dsv'][:-1,:-1]+s['dsv'][1:,:-1])) * (0.5*(s['dnu'][:-1,:-1]+s['dnu'][:-1,1:]))
    # s['dsdnz'][:-1,-1] = s['dsdnz'][:-1,-2]
    # s['dsdnz'][-1,:] = s['dsdnz'][-2,:]
    # s['dsdnzi'][:,:] = 1. / s['dsdnz']
    
    # # Alfaz, grid orientation in z-points
    # s['alfaz'][:-1,:] = np.arctan2(s['yu'][1:,:] - s['yu'][:-1,:], s['xu'][1:,:] - s['xu'][:-1,:])
    # s['alfaz'][-1,:] = s['alfaz'][-2,:]
    
    # # Alfau, grid orientation in u-points
    # s['alfau'][1:,:] = np.arctan2(s['yz'][1:,:] - s['yz'][:-1,:], s['xz'][1:,:] - s['xz'][:-1,:])
    # s['alfau'][0,:] = s['alfau'][1,:]
    
    # # Alfav, grid orientation in v-points
    # s['alfav'][:-1,:] = np.arctan2(s['yc'][1:,:] - s['yc'][:-1,:], s['xc'][1:,:] - s['xc'][:-1,:])
    # s['alfav'][-1,:] = s['alfav'][-2,:]
    
    # First compute angle with horizontal
    dx = s['x'][0,1] - s['x'][0,0]
    dy = s['y'][0,1] - s['y'][0,0]

    if dx == 0.:
        p['alpha'] = 90.
    else:
        p['alpha'] = np.rad2deg(np.arctan(dy/dx))

    if dx <= 0 and dy <= 0:
        p['alpha'] += 180.
        
    # Rotate grids to allign with horizontal
    xr, yr = rotate(s['x'], s['y'], p['alpha'], origin=(np.mean(s['x']), np.mean(s['y'])))

    # initialize y-dimension
    if ny == 0:
        s['dn'][:,:] = 1.
        s['ds'][:, 1:] = np.diff(s['x'], axis=1)
        s['ds'][:, 0] = s['ds'][:, 1]
    else:
        s['dn'][:,:] = ((yr[0,1]-yr[0,0])**2.+(xr[0,1]-xr[0,0])**2.)**0.5
        s['ds'][:,:] = ((xr[1,0]-xr[0,0])**2.+(yr[1,0]-yr[0,0])**2.)**0.5
         
    # compute cell areas
    s['dsdn'][:,:] = s['ds'] * s['dn']
    s['dsdni'][:,:] = 1. / s['dsdn']

    if ny > 0:
        dx_test = s['x'][0,1] - s['x'][0,0]
        dy_test = s['y'][1,0] - s['y'][0,0]

        if (dx_test <= 0.) or (dy_test <= 0.):
            logger.warn(format_log('WARNING: After rotation to the horizontal orientation, both x and y should be ascending. Otherwise he solver might produce false results. It is recommended to use the following function: create_grd (see https://github.com/openearth/aeolis-python/blob/AEOLIS_V2/tools/setup/setup_tools.py)'))

    return s, p