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
        Spatial grids'''
        
    # comment Lisa: dependence on moisture content is not yet implemented 
    # Can we do something with theta dependent on vegetation cover (larger rhoveg = larger theta?)    
        
    theta_stat = p['theta_stat']
    theta_dyn  = p['theta_dyn']
    
    s['theta_stat'] = theta_stat
    s['theta_dyn'] = theta_dyn

    
    #s['theta_stat'] += theta_stat
    #s['theta_dyn'] += theta_dyn
        
        
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

        nx = p['nx']+1
        ny = p['ny']+1

        #parameters

        tan_stat = np.tan(np.deg2rad(s['theta_stat']))
        tan_dyn = np.tan(np.deg2rad(s['theta_dyn']))
        


        E = 0.2

        grad_h_down = np.zeros((ny,nx,4))
        flux_down = np.zeros((ny,nx,4))
        slope_diff = np.zeros((ny,nx))
        grad_h = np.zeros((ny,nx))

        max_iter_ava = p['max_iter_ava']
        
        max_grad_h, grad_h, grad_h_down = calc_gradients(s['zb'], nx, ny, s['ds'], s['dn'], s[ 'zne'])
        

        
        s['gradh'] = grad_h.copy()

        initiate_avalanche = (max_grad_h > tan_stat) 

        if initiate_avalanche:

            for i in range(0,max_iter_ava):

                grad_h_down *= 0.
                flux_down *= 0.
                slope_diff *= 0.
                grad_h *= 0.

                max_grad_h, grad_h, grad_h_down = calc_gradients(s['zb'], nx, ny, s['ds'], s['dn'], s[ 'zne'])

                if max_grad_h < tan_dyn:
                    break

                # Calculation of flux

                grad_h_nonerod = (s['zb'] - s['zne']) / s['ds'] # HAS TO BE ADJUSTED!    
				
                ix = np.logical_and(grad_h > tan_dyn, grad_h_nonerod > 0)
                slope_diff[ix] = np.tanh(grad_h[ix]) - np.tanh(0.9*tan_dyn)    
                
                ix = grad_h_nonerod < grad_h - tan_dyn 
                slope_diff[ix] = np.tanh(grad_h_nonerod[ix])				                    

                ix = grad_h != 0
                
                if ny == 1:
                    #1D interpretation
                    flux_down[:,:,0][ix] = slope_diff[ix] * grad_h_down[:,:,0][ix] / grad_h[ix]
                    flux_down[:,:,2][ix] = slope_diff[ix] * grad_h_down[:,:,2][ix] / grad_h[ix]
                    
                    # Calculation of change in bed level
                    
                    q_in = np.zeros((ny,nx))
                    
                    q_out = 0.5*np.abs(flux_down[:,:,0]) + 0.5*np.abs(flux_down[:,:,2])
                    
                    q_in[0,1:-1] =   0.5*(np.maximum(flux_down[0,:-2,0],0.) \
                                        - np.minimum(flux_down[0,2:,0],0.) \
                                        + np.maximum(flux_down[0,2:,2],0.) \
                                        - np.minimum(flux_down[0,:-2,2],0.))
                else:
                    # 2D interpretation
                    flux_down[:,:,0][ix] = slope_diff[ix] * grad_h_down[:,:,0][ix] / grad_h[ix]
                    flux_down[:,:,1][ix] = slope_diff[ix] * grad_h_down[:,:,1][ix] / grad_h[ix]
                    flux_down[:,:,2][ix] = slope_diff[ix] * grad_h_down[:,:,2][ix] / grad_h[ix]
                    flux_down[:,:,3][ix] = slope_diff[ix] * grad_h_down[:,:,3][ix] / grad_h[ix]

                    # Calculation of change in bed level

                    q_in = np.zeros((ny,nx))

                    q_out = 0.5*np.abs(flux_down[:,:,0]) + 0.5* np.abs(flux_down[:,:,1]) + 0.5*np.abs(flux_down[:,:,2]) + 0.5* np.abs(flux_down[:,:,3])

                    q_in[1:-1,1:-1] =   0.5*(np.maximum(flux_down[1:-1,:-2,0],0.) \
                                        - np.minimum(flux_down[1:-1,2:,0],0.) \
                                        + np.maximum(flux_down[:-2,1:-1,1],0.) \
                                        - np.minimum(flux_down[2:,1:-1,1],0.) \

                                        + np.maximum(flux_down[1:-1,2:,2],0.) \
                                        - np.minimum(flux_down[1:-1,:-2,2],0.) \
                                        + np.maximum(flux_down[2:,1:-1,3],0.) \
                                        - np.minimum(flux_down[:-2,1:-1,3],0.))

                s['zb'] += E * (q_in - q_out)

    return s	    


def calc_gradients(zb, nx, ny, ds, dn, zne):
    '''Calculates the downslope gradients in the bed that are needed for
    avalanching module
 
    Parameters
    ----------
        
        
    Returns
    -------
    np.ndarray
        Downslope gradients in 4 different directions (nx*ny, 4)
    '''
    
    grad_h_down = np.zeros((ny,nx,4))

    # Calculation of slope (positive x-direction)
    grad_h_down[:,1:-1,0] = zb[:,1:-1] - zb[:,2:] 
    ix = zb[:,2:] > zb[:,:-2]
    grad_h_down[:,1:-1,0][ix] = - (zb[:,1:-1][ix] - zb[:,:-2][ix])    
    ix = np.logical_and(zb[:,2:]>zb[:,1:-1], zb[:,:-2]>zb[:,1:-1])
    grad_h_down[:,1:-1,0][ix] = 0.

    # Calculation of slope (positive y-direction)
    grad_h_down[1:-1,:,1] = zb[1:-1,:] - zb[2:,:]    
    ix = zb[2:,:] > zb[:-2,:]
    grad_h_down[1:-1,:,1][ix] = - (zb[1:-1,:][ix] - zb[:-2,:][ix])    
    ix = np.logical_and(zb[2:,:]>zb[1:-1,:], zb[:-2,:]>zb[1:-1,:])
    grad_h_down[1:-1,:,1][ix] = 0.

    # Calculation of slope (negative x-direction)
    grad_h_down[:,1:-1,2] = zb[:,1:-1] - zb[:,:-2]    
    ix = zb[:,:-2] > zb[:,2:]
    grad_h_down[:,1:-1,2][ix] = - (zb[:,1:-1][ix] - zb[:,2:][ix])    
    ix = np.logical_and(zb[:,:-2]>zb[:,1:-1], zb[:,2:]>zb[:,1:-1])
    grad_h_down[:,1:-1,2][ix] = 0.

    # Calculation of slope (negative y-direction)
    grad_h_down[1:-1,:,3] = zb[1:-1,:] - zb[:-2,:]    
    ix = zb[:-2,:] > zb[2:,:]
    grad_h_down[1:-1,:,3][ix] = - (zb[1:-1,:][ix] - zb[2:,:][ix])    
    ix = np.logical_and(zb[:-2,:]>zb[1:-1,:], zb[2:,:]>zb[1:-1,:])
    grad_h_down[1:-1,:,3][ix] = 0.
  
    if ny == 1:
        #1D interpretation
        grad_h_down[:,0,:]  = 0
        grad_h_down[:,-1,:] = 0

    else:
        # 2D interpretation
        grad_h_down[:,0,:]  = 0
        grad_h_down[:,-1,:] = 0
        grad_h_down[0,:,:]  = 0
        grad_h_down[-1,:,:] = 0
        
    grad_h_down[:,:,0] /= ds
    grad_h_down[:,:,1] /= dn
    grad_h_down[:,:,2] /= ds
    grad_h_down[:,:,3] /= dn

    grad_h2 = 0.5*grad_h_down[:,:,0]**2 + 0.5*grad_h_down[:,:,1]**2 + 0.5*grad_h_down[:,:,2]**2 + 0.5*grad_h_down[:,:,3]**2

    if 0: #Sierd_com; to be changed in future release
        ix = zb < zne + 0.005
        grad_h2[ix] = 0.

    grad_h = np.sqrt(grad_h2)

    max_grad_h = np.max(grad_h)

    return max_grad_h, grad_h, grad_h_down
    

