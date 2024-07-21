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
from scipy import optimize


# package modules
from aeolis.utils import *


# initialize logger
logger = logging.getLogger(__name__)

def duran_grainspeed(s, p):
    '''Compute grain speed according to Duran 2007 (p. 42)

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
    
    # Nr of fractions and grainsizes
    nf = p['nfractions']
    d = p['grain_size']
    
    # Copy basic variables
    z = s['zb']
    x = s['x']
    y = s['y']
    ny = p['ny']
    
    # Shear velocity and threshold
    ustar = s['ustar']
    ustars = s['ustars']
    ustarn = s['ustarn']
    ustar0 = s['ustar0']
    uth = s['uth0'] # uth0 or uth???
    uth0 = s['uth0'] 

    # Wind input for filling up ets/ets (udir), where ustar == 0
    uw = s['uw']
    uws = s['uws'] 
    uwn = s['uwn'] 
    
    # Other settings
    rhog = p['rhog']
    rhoa = p['rhoa']
    srho = rhog/rhoa
    kappa = p['kappa']

    # Add dimension for grain fractions
    g = np.repeat(p['g'], nf, axis = 0)
    v = np.repeat(p['v'], nf, axis = 0)

    # Initiate arrays
    ets = np.zeros(uth.shape)
    etn = np.zeros(uth.shape)
    ueff = np.zeros(uth.shape)
    ueff0 = np.zeros(uth.shape)
    u0 = np.zeros(uth.shape)
    us = np.zeros(uth.shape)
    un = np.zeros(uth.shape)
    u_approx = np.zeros(uth.shape)
    us_approx = np.zeros(uth.shape)
    un_approx = np.zeros(uth.shape)
    u  = np.zeros(uth.shape)

    dzs = np.zeros(z.shape)
    dzn = np.zeros(z.shape)

    # Extend 2D to 3D arrays for grain fractions
    ustar0 = np.repeat(ustar0[:,:,np.newaxis], nf, axis=2)
    ustar = np.repeat(ustar[:,:,np.newaxis], nf, axis=2)
    ustars = np.repeat(ustars[:,:,np.newaxis], nf, axis=2)
    ustarn = np.repeat(ustarn[:,:,np.newaxis], nf, axis=2)
    uw = np.repeat(uw[:,:,np.newaxis], nf, axis=2)
    uws = np.repeat(uws[:,:,np.newaxis], nf, axis=2)
    uwn = np.repeat(uwn[:,:,np.newaxis], nf, axis=2)
        
    # The formulations below are obtained from the PhD Thesis by Duran 
    # https://elib.uni-stuttgart.de/handle/11682/4810
    r       = 1.                                            # Ratio z/zm, assumption given in text p. 33
    c       = 14./(1.+1.4*r)                                # Given below eq. 1.27, p.32
    tv      = (v/g**2)**(1/3)                               # Timescale, eq 1.28, p.32, approx. 5.38 ms        
    lv      = (v**2/(p['Aa']**2*g*(srho-1)))**(1/3)         # Lengthscale, eq 1.45, p.36

    # Heights
    zm      = c * uth * tv                                  # Characteristic height of the saltation layer, eq 1.27, p.32, approx. 20 mm
    z0      = d/20.                                         # Grain based roughness layer, eq 1.29, p.32, approx 10 mu m
    z1      = 35. * lv                                      # reference height +- 35 mm, eq 1.46
  
    # Drag coefficient
    A       = 0.95                                          # Below eq 1.44, p.35
    B       = 5.12                                          # Below eq 1.44, p.35
    alpha   = 0.17 * d / lv                                 # Effective restitution coefficient, eq 1.47, p.36, approx. 0.42
    Sstar   = d/(4*v)*np.sqrt(g*d*(srho-1.))                # Fluid-sediment paramter, below eq 1.44, p.35
    Cd      = (4/3)*(A+np.sqrt(2*alpha)*B/Sstar)**2         # Drag coefficient, eq 1.44, p.35
    uf      = np.sqrt(4/(3*Cd)*(srho-1.)*g*d)               # Grain settling velocity - Jimnez and Madsen, 2003
    
    # Efficient wind velocity
    ueff    = (uth0 / kappa) * (np.log(z1 / z0))            # Effective wind velocity over a flat bed, 
    ueff0   = (uth0 / kappa) * (np.log(z1 / z0))

    # Compute Surface gradient
    dzs[:,1:-1] = (z[:,2:]-z[:,:-2])/(x[:,2:]-x[:,:-2])
    dzn[1:-1,:] = (z[:-2,:]-z[2:,:])/(y[:-2,:]-y[2:,:])
    if ny > 0:
        dzs[:,[0, -1]] = dzs[:,[0, -1]]
        dzn[[0, -1],:] = dzn[[0, -1],:]
    
    # Extend dimension of dzs/dzn (becomes dhs)
    dhs = np.repeat(dzs[:,:,np.newaxis], nf, axis = 2)
    dhn = np.repeat(dzn[:,:,np.newaxis], nf, axis = 2)
    
    # Compute ets/etn (wind direction)
    ix = (ustar > 0.)
    ets = uws / uw                                          # s-component of wind (where ustar == 0)
    etn = uwn / uw                                          # n-component of wind (where ustar == 0)
    ets[ix] = ustars[ix] / ustar[ix]                        # s-component of ustar
    etn[ix] = ustarn[ix] / ustar[ix]                        # n-component of ustar

    # Compute A parameter (Below Figure 1.13, p. 43)
    Axs = ets + 2*alpha*dhs
    Axn = etn + 2*alpha*dhn
    Ax = np.hypot(Axs, Axn)

    # Start looping over fractions
    for i in range(nf):  

        # Compute effective wind velocity, eq 1.60 p.42
        ueff[:,:,i] = (uth[:,:,i] / kappa) * (np.log(z1[i] / z0[i]) + (z1[i]/zm[:,:,i]) * (ustar[:,:,i]/uth[:,:,i]-1)) 
        ueff0[:,:,i] = (uth0[:,:,i] / kappa) * (np.log(z1[i] / z0[i]) + (z1[i]/zm[:,:,i]) * (ustar0[:,:,i]/uth0[:,:,i]-1))

        # Compute grainspeed over a flat bed (if dhs and dhn = 0)
        u0[:,:,i] = (ueff0[:,:,i] - uf[i] / (np.sqrt(2 * alpha[i])))
        
        # Compute grain speed: First approximation (eq 1.62) in case of gentle slopes
        us_approx[:,:,i] = (ueff[:,:,i] - uf[i] / (np.sqrt(2. * alpha[i]) * Ax[:,:,i])) * ets[:,:,i] \
                        - (np.sqrt(2*alpha[i]) * uf[i] / Ax[:,:,i]) * dhs[:,:,i]  
        
        un_approx[:,:,i] = (ueff[:,:,i] - uf[i] / (np.sqrt(2. * alpha[i]) * Ax[:,:,i])) * etn[:,:,i] \
                        - (np.sqrt(2*alpha[i]) * uf[i] / Ax[:,:,i]) * dhn[:,:,i] 
        
        u_approx[:,:,i] = np.hypot(us_approx[:,:,i], un_approx[:,:,i])

        # If 'regular' duran method is chosen, u_approx is the final solution
        if p['method_grainspeed'] == 'duran':
            us[:,:,i] = us_approx[:,:,i]
            un[:,:,i] = un_approx[:,:,i]
            u[:,:,i] = u_approx[:,:,i]

        # When duran_full is chosen the full formulation (eq 1.61) will be solved
        elif p['method_grainspeed'] == 'duran_full':

            # Transform into complex numbers
            u_approx_i = us_approx[:,:,i] + un_approx[:,:,i] * 1j
            veff_i = ueff[:,:,i] * ets[:,:,i] + ueff[:,:,i] * etn[:,:,i] * 1j
            dh_i = dhs[:,:,i] + dhn[:,:,i] * 1j
            uf_i = uf[i]
            alpha_i = alpha[i]

            # Solver van eq 1.61
            def solve_u(u_i: complex, veff_i: complex, uf_i: float, alpha_i: float, dh_i: complex) -> complex:
                return (veff_i - u_i) * np.abs(veff_i - u_i) / (uf_i ** 2) - u_i / (2 * alpha_i * np.abs(u_i)) - dh_i
            u_i = optimize.newton(solve_u, u_approx_i, maxiter=20, tol=0.05, args=(veff_i, uf_i, alpha_i, dh_i)) 

            # Transform back into components
            us[:,:,i] = np.real(u_i)
            un[:,:,i] = np.imag(u_i)
            u[:,:,i]= np.abs(u_i)

        else:
            logger.error('Grainspeed method not found!')
        
    return u0, us, un, u



def constant_grainspeed(s, p):
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
    
    uw  = s['uw'].copy() 
    uws = s['uws'].copy()
    uwn = s['uwn'].copy()
        
    ustar  = s['ustar'].copy()
    ustars = s['ustars'].copy()
    ustarn = s['ustarn'].copy()
    
    us = np.zeros(ustar.shape)
    un = np.zeros(ustar.shape)
    u  = np.zeros(ustar.shape)

    # u with direction of perturbation theory
    uspeed = 1. # m/s
    ix = ustar != 0
    us[ix] = uspeed * ustars[ix] / ustar[ix]
    un[ix] = uspeed * ustarn[ix] / ustar[ix]

    # u under the sep bubble
    sepspeed = 1.0 #m/s
    ix = (ustar == 0.)*(uw != 0.)
    us[ix] = sepspeed * uws[ix] / uw[ix]
    un[ix] = sepspeed * uwn[ix] / uw[ix]

    u = np.hypot(us, un)
                       
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
        nx = p['nx']      
        
        # u via grainvelocity:
        
        if p['method_grainspeed']=='duran' or p['method_grainspeed']=='duran_full':
            #the syntax inside grainspeed needs to be cleaned up
            u0, us, un, u = duran_grainspeed(s,p)
            s['u0'] = u0
            s['us'] = us
            s['un'] = un
            s['u']  = u
            
        elif p['method_grainspeed']=='windspeed':
            s['u0'] = s['uw'][:,:,np.newaxis].repeat(nf, axis=2)
            s['us'] = s['uws'][:,:,np.newaxis].repeat(nf, axis=2)
            s['un'] = s['uwn'][:,:,np.newaxis].repeat(nf, axis=2)
            s['u']  = s['uw'][:,:,np.newaxis].repeat(nf, axis=2) 
            u = s['u']
            
        elif p['method_grainspeed']=='constant':
            s = constant_grainspeed(s,p)
            u0 = s['u']
            us = s['us']
            un = s['un']
            u = s['u']
     
        
        ustar  = s['ustar'][:,:,np.newaxis].repeat(nf, axis=2)
        ustar0 = s['ustar0'][:,:,np.newaxis].repeat(nf, axis=2)
        
        uth    = s['uth']
        uthf   = s['uthf']
        uth0 = s['uth0']

        rhoa   = p['rhoa'] 
        g      = p['g']
        
        s['Cu']  = np.zeros(uth.shape)
        s['Cuf'] = np.zeros(uth.shape)
    
                
        ix = (ustar != 0.)*(u != 0.)
        
        if p['method_transport'].lower() == 'bagnold':
            s['Cu'][ix]  = np.maximum(0., p['Cb'] * rhoa / g * (ustar[ix] - uth[ix])**3 / u[ix])
            s['Cuf'][ix] = np.maximum(0., p['Cb'] * rhoa / g * (ustar[ix] - uthf[ix])**3 / u[ix])
            
            s['Cu0'][ix] = np.maximum(0., p['Cb'] * rhoa / g * (ustar0[ix] - uth0[ix])**3 / u[ix])
            
        elif p['method_transport'].lower() == 'bagnold_gs':
            Dref = 0.000250
            d = p['grain_size'][np.newaxis,np.newaxis,:].repeat(nx+1, axis=1)
            s['Cu'][ix]  = np.maximum(0., p['Cb'] * np.sqrt(d[ix]/Dref) * rhoa / g * (ustar[ix] - uth[ix])**3 / u[ix])
            s['Cuf'][ix] = np.maximum(0., p['Cb'] * np.sqrt(d[ix]/Dref) * rhoa / g * (ustar[ix] - uth[ix])**3 / u[ix])
            
            s['Cu0'][ix] = np.maximum(0., p['Cb'] * np.sqrt(d[ix]/Dref) * rhoa / g * (ustar[ix] - uth[ix])**3 / u[ix])
        
        elif p['method_transport'].lower() == 'kawamura':
            s['Cu'][ix]  = np.maximum(0., p['Ck'] * rhoa / g * (ustar[ix] + uth[ix])**2 * (ustar[ix] - uth[ix]) / u[ix])
            s['Cuf'][ix] = np.maximum(0, p['Ck'] * rhoa / g * (ustar[ix] + uthf[ix])**2 * (ustar[ix] - uthf[ix]) / u[ix])
        
        elif p['method_transport'].lower() == 'lettau':
            s['Cu'][ix]  = np.maximum(0., p['Cl'] * rhoa / g * ustar[ix]**2 * (ustar[ix] - uth[ix]) / u[ix])
            s['Cuf'][ix] = np.maximum(0., p['Cl'] * rhoa / g * ustar[ix]**2 * (ustar[ix] - uthf[ix]) / u[ix])

        elif p['method_transport'].lower() == 'dk':
            s['Cu'][ix]  = np.maximum(0., p['Cdk'] * rhoa / g * 0.8*uth[ix] * (ustar[ix]**2 - (0.8*uth[ix])**2) / u[ix])
            s['Cuf'][ix] = np.maximum(0., p['Cdk'] * rhoa / g * 0.8*uthf[ix] * (ustar[ix]**2 - (0.8*uthf[ix])**2) / u[ix])
            
            s['Cu0'][ix]  = np.maximum(0., p['Cdk'] * rhoa / g * 0.8*uth0[ix] * (ustar0[ix]**2 - (0.8*uth0[ix])**2) / u[ix])
         
        elif p['method_transport'].lower() == 'sauermann':
            alpha_sauermann = 0.35
            s['Cu'][ix]  = np.maximum(0., 2.* alpha_sauermann * rhoa / g * (ustar[ix]**2 - uth[ix]**2))
            s['Cuf'][ix] = np.maximum(0., 2.* alpha_sauermann * rhoa / g * (ustar[ix]**2 - uthf[ix]**2))
             
            s['Cu0'][ix]  = np.maximum(0., 2.* alpha_sauermann * rhoa / g * (ustar0[ix]**2 - uth0[ix]**2))
         
        elif p['method_transport'].lower() == 'vanrijn_strypsteen':
            s['Cu'][ix]  = np.maximum(0., p['Cb'] * rhoa / g * ((ustar[ix])**3 - (uth[ix])**3) / u[ix])
            s['Cuf'][ix] = np.maximum(0., p['Cb'] * rhoa / g * ((ustar[ix])**3 - (uth[ix])**3) / u[ix])
            
            s['Cu0'][ix] = np.maximum(0., p['Cb'] * rhoa / g * ((ustar0[ix])**3 - (uth0[ix])**3) / u[ix])
        
        else:
            logger.log_and_raise('Unknown transport formulation [%s]' % p['method_transport'], exc=ValueError)   

                     
    s['Cu']  *= p['accfac']
    s['Cuf'] *= p['accfac']
    s['Cu0'] *= p['accfac']
    
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

