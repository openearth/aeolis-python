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

from numba import njit

# initialize logger
logger = logging.getLogger(__name__)

def angele_of_repose(s,p):
    '''Determine the dynamic and static angle of repose.
    
    Both the critical dynamic and static angle of repose are spatial varying
    and depend on surface moisture content and roots of present vegetation
    and .... 
        
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
        
    # comment Lisa: dependence on moisture content is not yet implemented 
    # Can we do something with theta dependent on vegetation cover (larger rhoveg = larger theta?)    
        
    # theta_stat = p['theta_stat']
    theta_dyn  = p['theta_dyn']
    
    # s['theta_stat'] = theta_stat
    s['theta_dyn'] = theta_dyn
        
    return s


def avalanche(s, p):
    '''Avalanching occurs if bed slopes exceed critical slopes.
    
    Simulates the process of avalanching that is triggered by the exceedence
    of a critical static slope ``theta_stat`` by the bed slope. The iteration
    stops if the bed slope does not exceed the dynamic critical slope
    ``theta_dyn``.
    
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

    if p['process_avalanche']:
        nx = p['nx'] + 1
        ny = p['ny'] + 1

        # parameters - only dynamic angle used in loop for now. 
        # Static angle can be used for more complex criterions in later
        # tan_stat = np.tan(np.deg2rad(s['theta_stat']))
        tan_dyn = np.tan(np.deg2rad(s['theta_dyn']))

        E = 0.1

        max_iter_ava = p['max_iter_ava']

        # call the njit-compiled loop for performance
        zb, grad_h = avalanche_loop(
            s['zb'].copy(), s['zne'], s['ds'], s['dn'], nx, ny, E, max_iter_ava, tan_dyn
            )
   
        # Ensure water level is up-to-date with bed level
        s['zb'] = zb
        s['gradh'] = grad_h
        s['zs'] = s['SWL']
        ix = (s['zb'] > s['zs'])
        s['zs'][ix] = s['zb'][ix]

    return s

@njit(cache=True)
def avalanche_loop(zb, zne, ds, dn, nx, ny, E, max_iter_ava, tan_dyn):
    # Rewritten to use explicit loops and avoid numpy boolean indexing
    # Allocate temporaries once
    grad_h_down = np.zeros((ny, nx, 2))
    flux_down = np.zeros((ny, nx, 2))
    slope_diff = np.zeros((ny, nx))
    grad_h = np.zeros((ny, nx))
    for it in range(max_iter_ava):
        # Reset temporaries to zero
        grad_h_down.fill(0)
        flux_down.fill(0)
        slope_diff.fill(0)
        grad_h.fill(0)

        # first calculate the downslope gradients to see if there is avalanching
        # Compute downslope gradients grad_h_down (ny,nx,2), grad_h (ny,nx), and max_grad_h
        # Initialize
        max_grad_h = 0.0

        # Directions: 0 => +X, 1 => +Y
        for i in range(ny):
            for j in range(nx):
                # disable avalanching where zne >= zb
                if zne[i, j] >= zb[i, j]:
                    continue
                else:
                    center = zb[i, j]
                    # +X direction
                    g0 = 0.0
                    # Handle boundaries: set gradient to zero at edges
                    if j == 0 or j == nx - 1:
                        grad_h_down[i, j, 0] = 0.0
                    else:
                        right = zb[i, j + 1]
                        left = zb[i, j - 1]
                        if not ((right > center) and (left > center)):
                            if right > left:
                                g0 = left - center
                            else:
                                g0 = center - right
                        grad_h_down[i, j, 0] = g0 / ds[i, j]

                    # +Y direction
                    g1 = 0.0
                    if i == 0 or i == ny - 1:
                        grad_h_down[i, j, 1] = 0.0
                    else:
                        down = zb[i + 1, j]
                        up = zb[i - 1, j]
                        if not ((down > center) and (up > center)):
                            if down > up:
                                g1 = up - center
                            else:
                                g1 = center - down
                        grad_h_down[i, j, 1] = g1 / dn[i, j]

                # gradient magnitude and maximum
                gh2 = grad_h_down[i, j, 0] * grad_h_down[i, j, 0] + grad_h_down[i, j, 1] * grad_h_down[i, j, 1]
                gh = np.sqrt(gh2)
                grad_h[i, j] = gh
                # derive maximum slope
                if gh > max_grad_h:
                    max_grad_h = gh

        # ok now the gradients are calculated
        # these are max_grad_h, grad_h, grad_h_down
        # check for stopping criterion
        if max_grad_h < tan_dyn:
            break       
        
        # we continue to compute fluxes and update zb
        
        # compute grad_h_nonerod and slope_diff per cell using explicit loops
        for i in range(ny):
            for j in range(nx):
                if grad_h[i, j] > tan_dyn: 
                    slope_diff[i, j] = np.tanh(grad_h[i, j]) - np.tanh(0.9 * tan_dyn)
                    flux_down[i, j, 0] = slope_diff[i, j] * grad_h_down[i, j, 0]# / grad_h[i, j]
                    flux_down[i, j, 1] = slope_diff[i, j] * grad_h_down[i, j, 1]# / grad_h[i, j]
 
        # Build q_in and q_out from 2-component flux representation
        f_x = flux_down[:, :, 0]
        f_y = flux_down[:, :, 1]

        # preserve sign: compute positive outgoing components per direction
        out_east = np.maximum(f_x, 0.0)
        out_south = np.maximum(f_y, 0.0)

        # average with neighbor contributions at faces (keeps sign info)
        out_north = np.maximum(-f_y, 0.0)
        out_west  = np.maximum(-f_x, 0.0)

        q_out = out_east + out_west + out_south + out_north

        inc_west = np.zeros_like(f_x)
        # from west neighbor (positive f_x of west cell) with periodic wrap
        inc_west[:, 1:] = np.maximum(f_x[:, :-1], 0.0)
        inc_west[:, 0] = np.maximum(f_x[:, -1], 0.0)

        inc_east = np.zeros_like(f_x)
        # from east neighbor (negative f_x of east cell) with periodic wrap
        inc_east[:, :-1] = np.maximum(-f_x[:, 1:], 0.0)
        inc_east[:, -1] = np.maximum(-f_x[:, 0], 0.0)

        inc_north = np.zeros_like(f_y)
        # from north neighbor (positive f_y of north cell) with periodic wrap
        inc_north[1:, :] = np.maximum(f_y[:-1, :], 0.0)
        inc_north[0, :] = np.maximum(f_y[-1, :], 0.0)

        inc_south = np.zeros_like(f_y)
        # from south neighbor (negative f_y of south cell) with periodic wrap
        inc_south[:-1, :] = np.maximum(-f_y[1:, :], 0.0)
        inc_south[-1, :] = np.maximum(-f_y[0, :], 0.0)

        q_in = (inc_west + inc_east + inc_north + inc_south)

        # update bed level without non-erodible layer       
        zb += E * (q_in - q_out)

    return zb, grad_h