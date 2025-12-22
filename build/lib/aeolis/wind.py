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

import numpy as np
import logging
import operator
#import scipy.special
#import scipy.interpolate
from scipy import ndimage, misc
import matplotlib.pyplot as plt

# package modules
import aeolis.shear
from aeolis.utils import *


# initialize logger
logger = logging.getLogger(__name__)


def initialize(s, p):
    '''Initialize wind model

    '''

    # apply wind direction convention
    if isarray(p['wind_file']):
        if p['wind_convention'] == 'nautical':

            #fix issue associated with longshore winds/divide by zero
            ifix = p['wind_file'][:, 2] == 0.
            p['wind_file'][ifix, 2] = 0.01

        elif p['wind_convention'] == 'cartesian':
            #fix issue associated with longshore winds/divide by zero
            ifix = p['wind_file'][:, 2] == 270.
            p['wind_file'][ifix, 2] = 270.01

            p['wind_file'][:,2] = 270.0 - p['wind_file'][:,2]

        else:
            logger.log_and_raise('Unknown convention: %s' 
                                 % p['wind_convention'], exc=ValueError)

    # initialize wind shear model (z0 according to Duran much smaller)
    # Otherwise no Barchan
    z0    = calculate_z0(p, s)
    
    if p['process_shear']:
        if p['ny'] > 0:
            s['shear'] = aeolis.shear.WindShear(s['x'], s['y'], s['zb'],
                                                dx=p['dx'], dy=p['dy'],
                                                L=p['L'], l=p['l'], z0=z0,
                                                buffer_width=p['buffer_width'])
        else:
            s['shear'] = np.zeros(s['x'].shape)

    return s
    
def interpolate(s, p, t):
    '''Interpolate wind velocity and direction to current time step

    Interpolates the wind time series for velocity and direction to
    the current time step. The cosine and sine of the direction angle
    are interpolated separately to prevent zero-crossing errors. The
    wind velocity is decomposed in two grid components based on the
    orientation of each individual grid cell. In case of a
    one-dimensional model only a single positive component is used.

    Parameters
    ----------
    s : dict
        Spatial grids
    p : dict
        Model configuration parameters
    t : float
        Current time

    Returns
    -------
    dict
        Spatial grids

    '''
        
    if p['process_wind'] and p['wind_file'] is not None:

        uw_t = p['wind_file'][:,0]
        uw_s = p['wind_file'][:,1]
        uw_d = p['wind_file'][:,2] / 180. * np.pi

        s['uw'][:,:] = interp_circular(t, uw_t, uw_s)
        s['udir'][:,:] = np.arctan2(interp_circular(t, uw_t, np.sin(uw_d)),
                                    interp_circular(t, uw_t, np.cos(uw_d))) * 180. / np.pi

    s['uws'] = - s['uw'] * np.sin((-p['alfa'] + s['udir']) / 180. * np.pi)        # alfa [deg] is real world grid cell orientation (clockwise)
    s['uwn'] = - s['uw'] * np.cos((-p['alfa'] + s['udir']) / 180. * np.pi)

          
    s['uw'] = np.abs(s['uw'])
    
    # Compute wind shear velocity
    kappa = p['kappa']
    z     = p['z']
    z0    = calculate_z0(p, s)                                                                                                           
    
    s['ustars'] = s['uws'] * kappa / np.log(z/z0)
    s['ustarn'] = s['uwn'] * kappa / np.log(z/z0) 
    s['ustar']  = np.hypot(s['ustars'], s['ustarn'])
    
    s = velocity_stress(s,p)
        
    s['ustar0'] = s['ustar'].copy()
    s['ustars0'] = s['ustar'].copy()
    s['ustarn0'] = s['ustar'].copy()
        
    s['tau0'] = s['tau'].copy()
    s['taus0'] = s['taus'].copy()
    s['taun0'] = s['taun'].copy()
    
    return s
    
def calculate_z0(p, s):
    '''Calculate z0 according to chosen roughness method

    The z0 is required for the calculation of the shear velocity. Here, z0
    is calculated based on a user-defined method. The constant method defines 
    the value of z0 as equal to k (z0 = ks). This was implemented to ensure 
    backward compatibility and does not follow the definition of Nikuradse 
    (z0 = k / 30). For following the definition of Nikuradse use the method 
    constant_nikuradse. The mean_grainsize_initial method uses the intial
    mean grain size ascribed to the bed (grain_dist and grain_size in the 
    input file) to calculate the z0. The median_grainsize_adaptive bases the 
    z0 on the median grain size (D50) in the surface layer in every time step. 
    The resulting z0 is variable accross the domain (x,y). The 
    strypsteen_vanrijn method is based on the roughness calculation in their 
    paper. 

    Parameters
    ----------
    s : dict
        Spatial grids
    p : dict
        Model configuration parameters

    Returns
    -------
    array
        z0

    '''
    if p['method_roughness'] == 'constant':
        z0    = p['k']  # Here, the ks (roughness length) is equal to the z0, this method is implemented to assure backward compatibility. Note, this does not follow the definition of z0 = ks /30 by Nikuradse    
    if p['method_roughness'] == 'constant_nikuradse':
        z0    = p['k'] / 30   # This equaion follows the definition of the bed roughness as introduced by Nikuradse
    if p['method_roughness'] == 'mean_grainsize_initial': #(based on Nikuradse and Bagnold, 1941), can only be applied in case with uniform grain size and is most applicable to a flat bed
        z0    = np.sum(p['grain_size']*p['grain_dist']) / 30.
    if p['method_roughness'] == 'mean_grainsize_adaptive': # makes Nikuradse roughness method variable through time and space depending on grain size variations
        z0    = calc_mean_grain_size(p, s) / 30.
    if p['method_roughness'] == 'median_grainsize_adaptive': # based on Sherman and Greenwood, 1982 - only appropriate for naturally occurring grain size distribution
        d50 = calc_grain_size(p, s, 50)
        z0 = 2*d50 / 30.
    if p['method_roughness'] == 'vanrijn_strypsteen': # based on van Rijn and Strypsteen, 2019; Strypsteen et al., 2021
        if len(p['grain_dist']) == 1:  # if one grainsize is used the d90 is calculated with the d50 
            d50 = p['grain_size']
            d90 = 2*d50
        else:
            d50 = calc_grain_size(p, s, 50) #calculate d50 and d90 per cell.
            d90 = calc_grain_size(p, s, 90)
        
        ustar_grain_stat = p['kappa'] * (s['uw'] / np.log(30*p['z']/d90))
        
        ustar_th_B = 0.1 * np.sqrt((p['rhog'] - p['rhoa']) / p['rhoa'] * p['g'] * d50) # Note that Aa could be filled in in the spot of 0.1
        
        T = (np.square(ustar_grain_stat) - np.square(ustar_th_B))/np.square(ustar_th_B) # T represents different phases of the transport related to the saltation layer and ripple formation
        #T[T < 0] = 0
        
        alpha1 = 15
        alpha2 = 1
        gamma_r = 1 + 1/T
        z0    = (d90 + alpha1 * gamma_r * d50 * np.power(T, alpha2)) / 30
        
        # print("z0 in wind.py = ", np.mean(z0))
    return z0


def shear(s,p):
    
    # Compute shear velocity field (including separation)

    if 'shear' in s.keys() and p['process_shear'] and p['ny'] > 0:
        
        s['shear'](x=s['x'], y=s['y'], z=s['zb'],
                   taux=s['taus'], tauy=s['taun'],
                   u0=s['uw'][0,0], udir=s['udir'][0,0],
                   process_separation = p['process_separation'],
                   c = p['c_b'],
                   mu_b = p['mu_b'],
                   taus0 = s['taus0'][0,0], taun0 = s['taun0'][0,0],
                   sep_filter_iterations=p['sep_filter_iterations'],
                   zsep_y_filter=p['zsep_y_filter'])

        s['taus'], s['taun'] = s['shear'].get_shear()
        s['tau'] = np.hypot(s['taus'], s['taun'])
        
        s = stress_velocity(s,p)
                               
        # Returns separation surface     
        if p['process_separation']:
            s['hsep'] = s['shear'].get_separation()
            s['zsep'] = s['hsep'] + s['zb']
        

    elif p['process_shear'] and p['ny'] == 0: #NTC - Added in 1D only capabilities
        s = compute_shear1d(s, p)
        s = stress_velocity(s, p)

        if p['process_separation']:
            zsep = separation1d(s, p)
            s['zsep'] = zsep
            s['hsep'] = s['zsep'] - s['zb']
            tau_sep = 0.5
            slope = 0.2  # according to Durán 2010 (Sauermann 2001: c = 0.25 for 14 degrees)
            delta = 1. / (slope * tau_sep)
            zsepdelta = np.minimum(np.maximum(1. - delta * s['hsep'], 0.), 1.)
            s['taus'] *= zsepdelta
            s['taun'] *= zsepdelta
            s = stress_velocity(s, p)

    # if p['process_nelayer']:
    # if p['th_nelayer']:

    #     ustar = s['ustar'].copy()
    #     ustars = s['ustars'].copy()
    #     ustarn = s['ustarn'].copy()
            
    #     s['zne'][:,:] = p['ne_file']
            
    #     ix = s['zb'] <= s['zne']
    #     s['ustar'][ix] = np.maximum(0., s['ustar'][ix] - (s['zne'][ix]-s['zb'][ix])* (1/p['layer_thickness']) * s['ustar'][ix])
        
    #     ix = ustar != 0.
    #     s['ustars'][ix] = s['ustar'][ix] * (ustars[ix] / ustar[ix])
    #     s['ustarn'][ix] = s['ustar'][ix] * (ustarn[ix] / ustar[ix])


    return s

def velocity_stress(s, p):

    s['tau'] = p['rhoa'] * s['ustar'] ** 2

    ix = s['ustar'] > 0.
    s['taus'][ix] = s['tau'][ix]*s['ustars'][ix]/s['ustar'][ix]
    s['taun'][ix] = s['tau'][ix]*s['ustarn'][ix]/s['ustar'][ix]
    s['tau'] = np.hypot(s['taus'], s['taun'])

    ix = s['ustar'] == 0.
    s['taus'][ix] = 0.
    s['taun'][ix] = 0.
    s['tau'][ix] = 0.

    return s

def stress_velocity(s, p):

    s['ustar'] = np.sqrt(s['tau'] / p['rhoa'])

    ix = s['tau'] > 0.
    s['ustars'][ix] = s['ustar'][ix] * s['taus'][ix] / s['tau'][ix]
    s['ustarn'][ix] = s['ustar'][ix] * s['taun'][ix] / s['tau'][ix]

    ix = s['tau'] == 0.
    s['ustar'][ix] = 0.
    s['ustars'][ix] = 0.
    s['ustarn'][ix] = 0.

    return s



def compute_shear1d(s, p):
    '''Compute wind shear perturbation for given free-flow wind
    speed on computational grid. based on same implementation in Duna'''

    tau = s['tau'].copy()
    taus = s['taus'].copy()
    taun = s['taun'].copy()
    ets = np.zeros(s['tau'].shape)
    etn = np.zeros(s['tau'].shape)
    ix = tau != 0
    ets[ix] = taus[ix] / tau[ix]
    etn[ix] = taun[ix] / tau[ix]

    x = s['x'][0,:]
    zb = s['zb'][0,:]

    #Bart: check for negative wind direction
    if np.sum(taus) < 0:
        x = np.flip(x)
        zb = np.flip(zb)


    dzbdx = np.zeros(x.shape)
    tau_over_tau0 = np.zeros(x.shape)
    dx = x[1] - x[0]
    dx = np.abs(dx)
    dzbdx[1:-1] = (zb[2:] - zb[0:-2]) / 2 / dx
    nx = x.size - 1
    alfa = 3
    beta = 1
    for i in range(nx + 1):
        integ = 0
        startval = i - nx
        endval = i - 1
        for j in np.arange(startval, endval + 1):
            if j != 0:
                integ = integ + dzbdx[i - j] / (j * np.pi)
        tau_over_tau0[i] = alfa * (integ + beta * dzbdx[i]) + 1
        tau_over_tau0[i] = np.maximum(tau_over_tau0[i], 0.1)

    #should double check this - but i think this is right. duna is in u10, so slightly different

    #Bart: check for negative wind direction
    if np.sum(taus) < 0:
        tau_over_tau0 = np.flip(tau_over_tau0)

    s['tau'] = tau * tau_over_tau0
    s['taus'] = s['tau'] * ets
    s['taun'] = s['tau'] * etn

    return s


def separation1d(s, p):
    # Initialize grid and bed dimensions

    #load relevant input
    x = s['x'][0,:]
    #x = s['x']
    z = s['zb'][0,:]
    dx = p['dx']
    dy = dx
    c = p['c_b']
    mu_b = p['mu_b']
    nx = np.size(z)
    udir = s['udir'][0][0]

    #make the grids 2d to utilize same code as in the shear module
    ny = 3
    #z = np.matlib.repmat(z, ny, 1)
    z = np.tile(z, [ny, 1])

    if udir < 360:
        udir = udir + 360

    if udir > 360:
        udir = udir - 360


    if udir > 180 and udir < 360:
        udir = np.abs(udir-270)
        dx = dx / np.cos(udir * np.pi / 180)
        dy = dx
        direction = 1
    elif udir == 180:
        dx = 0.0001
        direction = 1
    elif udir == 360:
        dx = 0.0001
        direction = 1
    else:
        udir = np.abs(udir-90)
        dx = dx / np.cos(udir * np.pi / 180)
        dy = dx
        direction = 2

    x = np.tile(x, [ny, 1])

    if direction == 2:
        z = np.flip(z, 1)

    #y = np.matrix.transpose(np.tile(y, [ny, 1]))

    # Initialize arrays
    dzx = np.zeros(z.shape)
    dzdx0 = np.zeros(z.shape)
    dzdx1 = np.zeros(z.shape)

    stall = np.zeros(z.shape)
    bubble = np.zeros(z.shape)

    k = np.array(range(0, nx))

    zsep = z.copy()  # total separation bubble

    zsep0 = np.zeros(z.shape)  # zero-order separation bubble surface
    zsep1 = np.zeros(z.shape)  # first-oder separation bubble surface

    zfft = np.zeros((ny, nx), dtype=complex)

    # Compute bed slope angle in x-dir
    dzx[:, :-1] = np.rad2deg(np.arctan((z[:, 1:] - z[:, :-1]) / dx))
    dzx[:, 0] = dzx[:, 1]
    dzx[:, -1] = dzx[:, -2]

    # Determine location of separation bubbles
    '''Separation bubble exist if bed slope angle (lee side) 
    is larger than max angle that wind stream lines can 
    follow behind an obstacle (mu_b = ..)'''

    stall += np.logical_and(abs(dzx) > mu_b, dzx < 0.)

    stall[:, 1:-1] += np.logical_and(stall[:, 1:-1] == 0, stall[:, :-2] > 0., stall[:, 2:] > 0.)

    # Define separation bubble
    bubble[:, :-1] = np.logical_and(stall[:, :-1] == 0., stall[:, 1:] > 0.)

    # Shift bubble back to x0: start of separation bubble
    p = 2
    bubble[:, :-p] = bubble[:, p:]
    bubble[:, :p] = 0

    bubble = bubble.astype(int)

    # Count separation bubbles
    n = np.sum(bubble)
    bubble_n = np.asarray(np.where(bubble == True)).T

    # Walk through all separation bubbles and determine polynoms
    for k in range(0, n):

        i = bubble_n[k, 1]
        j = bubble_n[k, 0]

        ix_neg = (dzx[j, i + 5:] >= 0)  # i + 5??

        if np.sum(ix_neg) == 0:
            zbrink = z[j, i]  # z level of brink at z(x0)
        else:
            zbrink = z[j, i] - z[j, i + 5 + np.where(ix_neg)[0][0]]

        # Zero order polynom
        dzdx0 = (z[j, i - 1] - z[j, i - 2]) / dx

        # if dzdx0 > 0.1:
        #    dzdx0 = 0.1

        a = dzdx0 / c

        ls = np.minimum(np.maximum((3. * zbrink / (2. * c) * (1. + a / 4. + a ** 2 / 8.)), 0.1), 200.)
        a2 = -3 * zbrink / ls ** 2 - 2 * dzdx0 / ls
        a3 = 2 * zbrink / ls ** 3 + dzdx0 / ls ** 2
        i_max = min(i + int(ls / dx), int(nx - 1))
        xs = x[j, i:i_max] - x[j, i]
        zsep0[j, i:i_max] = (a3 * xs ** 3 + a2 * xs ** 2 + dzdx0 * xs + z[j, i])

        # Zero order filter
        Cut = 1.5
        dk = 2.0 * np.pi / (np.max(x))
        zfft[j, :] = np.fft.fft(zsep0[j, :])
        zfft[j, :] *= np.exp(-(dk * k * dx) ** 2 / (2. * Cut ** 2))
        zsep0[j, :] = np.real(np.fft.ifft(zfft[j, :]))

        # First order polynom
        dzdx1 = (zsep0[j, i - 1] - zsep0[j, i - 2]) / dx
        a = dzdx1 / c
        ls = np.minimum(np.maximum((3. * z[j, i] / (2. * c) * (1. + a / 4. + a ** 2 / 8.)), 0.1), 200.)
        a2 = -3 * z[j, i] / ls ** 2 - 2 * dzdx1 / ls
        a3 = 2 * z[j, i] / ls ** 3 + dzdx1 / ls ** 2
        i_max1 = min(i + int(ls / dx), int(nx - 1))
        xs1 = x[j, i:i_max1] - x[j, i]

        # Combine Seperation Bubble
        zsep1[j, i:i_max1] = (a3 * xs1 ** 3 + a2 * xs1 ** 2 + dzdx1 * xs1 + z[j, i])
        zsep[j, i:i_max] = np.maximum(zsep1[j, i:i_max], z[j, i:i_max])

    # Smooth surface of separation bubbles over y direction
    zsep = ndimage.gaussian_filter1d(zsep, sigma=0.2, axis=0)
    ilow = zsep < z
    zsep[ilow] = z[ilow]

    #remove the 2d aspect of results
    zsepout = zsep[1,:]

    if direction == 2:
        zsepout = np.flip(zsepout)

    return zsepout