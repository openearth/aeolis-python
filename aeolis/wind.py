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

from scipy import ndimage
from scipy.optimize import root_scalar


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
            if p['method_shear'] == 'fft':
                s['shear'] = aeolis.shear.WindShear(s['x'], s['y'], s['zb'],
                                                    dx=p['dx'], dy=p['dy'],
                                                    L=p['L'], l=p['l'], z0=z0,
                                                    buffer_width=p['buffer_width'])
            else:
                s['shear'] = aeolis.rotation.rotationClass(s['x'], s['y'], s['zb'],
                                                          dx=p['dx'], dy=p['dy'],
                                                          buffer_width=100)
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
        # defining the wind inputs the same as the timestep speeds up the simulation significantly 
        if (np.any(p['wind_file'][:,0]==t)):
            s['uw'][:,:] = p['wind_file'][p['wind_file'][:,0]==t,1][0]      # this extra bracket is needed to accound for messy input files
            s['udir'][:,:] = p['wind_file'][p['wind_file'][:,0]==t,2][0] 
        
        # alternatively, wind inputs are interpolated based on a circular interpolation.
        # this is more time expensive
        else:
            uw_t = p['wind_file'][:,0]
            uw_s = p['wind_file'][:,1]
            uw_d = p['wind_file'][:,2] / 180. * np.pi

            s['uw'][:,:] = interp_circular_nearest(t, uw_t, uw_s)
            
            s['udir'][:,:] = np.arctan2(interp_circular_nearest(t, uw_t, np.sin(uw_d)),
                                        interp_circular_nearest(t, uw_t, np.cos(uw_d))) * 180. / np.pi


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
    s['ustars0'] = s['ustars'].copy()
    s['ustarn0'] = s['ustarn'].copy()
        
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
    return z0


def shear(s,p):
    
    # Compute shear velocity field (including separation)

    if 'shear' in s.keys() and p['process_shear'] and p['ny'] > 0 and p['method_shear'] == 'duna2d':

        shear_params = {'alfa': 3, 'beta': 1, 'c': p['c_b'], 'mu_b': p['mu_b'], 'sep_filter_iterations': p['sep_filter_iterations'], 'zsep_y_filter': p['zsep_y_filter'], 'process_separation': p['process_separation'], 'tau_sep': 0.5, 'slope': 0.2, 'rhoa': p['rhoa'], 'shear_type': p['method_shear']}

        s['shear'](x=s['x'], y=s['y'], z=s['zb'], udir=s['udir'][0, 0], ustar=s['ustar'], tau=s['tau'],  hveg=s['hveg'], type_flag = 3, params=shear_params)

        s['tau'] = s['shear'].get_tau()
        s['taus'] = - s['tau'] * np.sin((-p['alfa'] + s['udir']) / 180. * np.pi)
        s['taun'] = - s['tau'] * np.cos((-p['alfa'] + s['udir']) / 180. * np.pi)

        s = stress_velocity(s, p)

    elif 'shear' in s.keys() and p['process_shear'] and p['ny'] > 0 and p['method_shear'] == 'quasi2d':

        z0 = calculate_z0(p, s)

        shear_params = {'L': p['L'], 'l': p['l'], 'kappa': p['kappa'], 'viscNu': 1.5e-5, 'z0': z0, 'rho': p['rhoa'], 'sep_angle': 14, 'shear_type': p['method_shear'], 'process_separation': p['process_separation'], 'c': p['c_b'], 'mu_b': p['mu_b'], 'sep_filter_iterations': p['sep_filter_iterations'], 'zsep_y_filter': p['zsep_y_filter'], 'tau_sep': 0.5, 'slope': 0.2}

        s['shear'](x=s['x'], y=s['y'], z=s['zb'], udir=s['udir'][0, 0], ustar=s['ustar'], tau=s['tau'],  hveg=s['hveg'], type_flag = 3, params=shear_params)
        s['tau'] = s['shear'].get_tau()

        s['taus'] = - s['tau'] * np.sin((-p['alfa'] + s['udir']) / 180. * np.pi)
        s['taun'] = - s['tau'] * np.cos((-p['alfa'] + s['udir']) / 180. * np.pi)
        s = stress_velocity(s, p)

    elif 'shear' in s.keys() and p['process_shear'] and p['ny'] > 0:
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
    dx = s['ds'][0,0] # p['dx']
    dy = dx
    c = p['c_b']
    mu_b = p['mu_b']
    nx = np.size(z)
    udir = s['udir'][0][0]

    #make the grids 2d to utilize same code as in the shear module
    ny = 3
    #z = np.matlib.repmat(z, ny, 1)
    z = np.tile(z, [ny, 1])

    if udir < 0:
        udir = udir + 360

    if udir > 360:
        udir = udir - 360


    if udir > 180 and udir < 360:
        udir = np.abs(udir-270)
        dx = dx / np.cos(udir * np.pi / 180)
        dy = dx
        direction = 1
        idir = 1
    elif udir == 180:
        dx = 0.0001
        direction = 1
        idir = 1
    elif udir == 360:
        dx = 0.0001
        direction = 1
        idir = 1
    else:
        udir = np.abs(udir-90)
        dx = dx / np.cos(udir * np.pi / 180)
        dy = dx
        direction = 2
        idir = -1

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
    p = 1
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

        ix_neg = (dzx[j, i+idir*5:] >= 0) 

        if np.sum(ix_neg) == 0:
            zbrink = z[j, i]  # z level of brink at z(x0)
        else:
            zbrink = z[j, i] - z[j,i+idir*5+idir*np.where(ix_neg)[0][0]]
            
        # Zero order polynom
        dzdx0 = (z[j,i] - z[j,i-3]) / (3.*dx)

        a = dzdx0 / c
        ls = np.minimum(np.maximum((3.*zbrink/(2.*c) * (1. + a/4. + a**2/8.)), 0.1), 200.)
        
        a2 = -3 * zbrink/ls**2 - 2 * dzdx0 / ls
        a3 =  2 * zbrink/ls**3 +     dzdx0 / ls**2
        
        i_max = min(i+int(ls/dx)+1,int(nx-1))

        if idir == 1:
            xs = x[j,i:i_max] - x[j,i]
        else:
            xs = -(x[j,i:i_max] - x[j,i])

        zsep0[j, i:i_max] = (a3 * xs ** 3 + a2 * xs ** 2 + dzdx0 * xs + z[j, i])

        # # Zero order filter
        # Cut = 1.5
        # dk = 2.0 * np.pi / (np.max(x))
        # zfft[j, :] = np.fft.fft(zsep0[j, :])
        # zfft[j, :] *= np.exp(-(dk * k * dx) ** 2 / (2. * Cut ** 2))
        # zsep0[j, :] = np.real(np.fft.ifft(zfft[j, :]))

        # # First order polynom
        # dzdx1 = (zsep0[j, i - 1] - zsep0[j, i - 2]) / dx
        # a = dzdx1 / c
        # ls = np.minimum(np.maximum((3. * z[j, i] / (2. * c) * (1. + a / 4. + a ** 2 / 8.)), 0.1), 200.)
        # a2 = -3 * z[j, i] / ls ** 2 - 2 * dzdx1 / ls
        # a3 = 2 * z[j, i] / ls ** 3 + dzdx1 / ls ** 2
        # i_max1 = min(i + int(ls / dx), int(nx - 1))
        # xs1 = x[j, i:i_max1] - x[j, i]

        # # Combine Seperation Bubble
        # zsep1[j, i:i_max1] = (a3 * xs1 ** 3 + a2 * xs1 ** 2 + dzdx1 * xs1 + z[j, i])
        # zsep[j, i:i_max] = np.maximum(zsep1[j, i:i_max], z[j, i:i_max])
        zsep[j, i:i_max] = np.maximum(zsep0[j, i:i_max], z[j, i:i_max])

    # Smooth surface of separation bubbles over y direction
    # zsep = ndimage.gaussian_filter1d(zsep, sigma=0.2, axis=0)
    ilow = zsep < z
    zsep[ilow] = z[ilow]

    #remove the 2d aspect of results
    zsepout = zsep[1,:]

    if direction == 2:
        zsepout = np.flip(zsepout)

    return zsepout



class Kroy2002():

    def __init__(self,xgrid, z, z0, dune_length, rho, mu, tau, sep=False, sepInd=None, sep_ang=14, kappa=0.4):
        self.xgrid = xgrid
        self.z0 = z0
        self.z = z
        self.L = 0.5 * dune_length
        self.Phi = self.L / z0
        self.rho = rho
        self.mu = mu
        self._kappa = kappa
        self.tau_inf = tau
        self.tau_p = None
        self.tau_tot = None
        self.sep_slope = np.tan(np.radians(sep_ang))

        self.dz = np.zeros(self.xgrid.shape)
        ## check for equal spacing
        dx = self.xgrid[1:] - self.xgrid[0:-1]
        if np.all(np.abs(dx[0] - dx) <= 1e-6):
            dx = dx[0]
        else:
            raise Exception('Unequally spaced dune profile data is not supported!')

        ## use second order central differences except on ends
        self.dz[0] = (-3.0 * self.z[0] + 4.0 * self.z[1] - self.z[2]) / (2.0 * dx)
        self.dz[-1] = (3.0 * self.z[-1] - 4.0 * self.z[-2] + self.z[-3]) / (2.0 * dx)
        self.dz[1:-1] = (self.z[2:] - self.z[0:-2]) / (2.0 * dx)

        if sep and sepInd is None:
            raise Exception('Must provide and index for the location of separation if "sep=True"')
        self._sep = sep
        self._sepInd = sepInd

    def calc_tau_p(self):
        # self._Phi = L/z0

        # check for separation (max slope of -14 deg.)
        #if self.dz.min() < -self.sep_slope or self._sep:
        #    print('Separation possible based on dune profile!')
        #    self._build_envelope()
        #else:
        #    self.env = self.z
        self.env = self.z

        # determine alpha and beta using eq. 13 from Kroy 03/2002
        def phi_func(x):
            # implicit equation for determining phi
            val = x - (2.0 * self._kappa ** 2 * self.Phi / np.log(x))
            deriv1 = 1 + (2.0 * self._kappa ** 2 * self.Phi / (x * np.log(x) ** 2))
            # deriv2 = (2.0*self._kappa**2*Phi) * (np.log(x) + 2.0)/(x**2 * np.log(x)**3)
            return val, deriv1  # , deriv2

        phi_sol = root_scalar(phi_func,
                              method='newton',
                              x0=2.0 * self._kappa ** 2 * self.Phi / np.log(2.0 * self._kappa ** 2 * self.Phi),
                              fprime=True)  # , fprime2=True)

        if (phi_sol.converged):
            phi = phi_sol.root
        else:
            exit('Implicit equation for phi did not converge! Exiting!')

        fact = 1 + np.log(phi) + 2.0 * np.log(0.5 * np.pi) + 4.0 * np.euler_gamma
        self._alpha = ((np.log(self.Phi ** 2 / np.log(self.Phi)) ** 2) / (2.0 * np.log(phi) ** 3)) * fact
        self._beta = np.pi / fact

        # calculate shear stress perturbation (2D) using eq. 16a from Kroy 03/2002
        x = self.xgrid
        L = self.L
        fft_out = np.fft.rfft(self.env)
        freq = 2.0 * np.pi * np.fft.rfftfreq(x.shape[0], x[1] - x[0])
        fft_tau = self._alpha * (freq + 1j * self._beta * np.abs(freq)) * fft_out
        self.tau_p = np.fft.irfft(fft_tau)
        #self.tau_p_ana = 2.0 * self._alpha * (
        #            np.pi ** (-0.5) - (x / L) * (self._beta + erfi(x / L)) * np.exp(-(x / L) ** 2))


        #dx = x[1]-x[0]
        #x2 = x[1:-1]-dx/2

        #print(np.size(x))
        #print(np.size(x2))
        #print(np.size(self.tau_p))

       # tau_p_interp = np.interp(x, x2, self.tau_p)

        #func = np.interp1d(x2, self.tau_p,  axis=0,  # interpolate along columns
        #                bounds_error=False,
        #                kind='linear')

        #tau_p_interp = func(x)

        #self.tau_tot = self.tau_inf * (1 + tau_p_interp)

    def calc_init_y(self, yPlus):
        # calculate smallest needed y for desired y+
        max_loc_u_star = np.sqrt(self.tau_tot.max() / self.rho)
        first_y = (self.mu / self.rho) * yPlus / max_loc_u_star
        return 2.0 * first_y  # y+ taken as centroid of cell, want to return total cell thickness

    def S(self, x):
        s = np.ones(x.shape)
        s[x >= 0.5 * np.pi] = 0.0
        s[x <= 0.5 * np.pi] = 0.0
        return s

    def _build_envelope(self):
        def interp_end_vals(xr):
            """Interpolates the elevation and slope at the reattachment point given by
            xr"""

            # get index of value just less than xr
            #print(self.xgrid.max(), xr)
            xrInd = np.argwhere(self.xgrid > xr)[0][0] - 1

            # interp the values
            fact = ((xr - self.xgrid[xrInd]) / (self.xgrid[xrInd + 1] - self.xgrid[xrInd]))
            zr = self.z[xrInd] + fact * (self.z[xrInd + 1] - self.z[xrInd])
            dzr = self.dz[xrInd] + fact * (self.dz[xrInd + 1] - self.dz[xrInd])

            return zr, dzr

        def solve_coef(xd, xr, zd, dzd, zr, dzr):
            """Solves for the coefficients of a cubic polynomial given the value (zd, zr) and
            derivative (dzd, dzr) at either end (xd, xr) """
            # setup the system
            A = np.array(
                [[xd ** 3, xd ** 2, xd, 1.0],
                 [3.0 * xd ** 2, 2.0 * xd, 1.0, 0.0],
                 [xr ** 3, xr ** 2, xr, 1.0],
                 [3.0 * xr ** 2, 2.0 * xr, 1.0, 0.0]]
            )
            b = [zd, dzd, zr, dzr]

            # solve the system
            return np.linalg.solve(A, b)

        def sepLengthRoot(xr):
            """Function for determining the separation bubble length which results in a
            cubic polynomial with max slope equal to the limit slope"""
            # interpolate endpoint values
            zr, dzr = interp_end_vals(xr)

            # solve for the polynomial
            poly = solve_coef(xSepVal, xr, zd, dzd, zr, dzr)

            # calc the max derivative
            deriv_max = -(poly[1] ** 2 / (3.0 * poly[0])) + poly[2]

            # return the difference between the deriv and the limit
            return deriv_max - (-self.sep_slope)

        # location of potential separation
        if self._sepInd is None:
            xSepInd = np.argwhere(self.dz < -self.sep_slope)[0][0] - 1
        else:
            xSepInd = self._sepInd

        # useful values at separation point
        xSepVal = self.xgrid[xSepInd]
        zd = self.z[xSepInd]
        # back-difference for left-side derivative at separation point
        dzd = (self.z[xSepInd] - self.z[xSepInd - 1]) / (self.xgrid[xSepInd] - self.xgrid[xSepInd - 1])

        # bubble length initial guess (Kroy2002 Eq. 29)
        nu = dzd / self.sep_slope
        Lb = 1.5 * (zd / self.sep_slope) * (
                1.0 + 0.25 * nu + 0.125 * nu ** 2)

        print('Predicted separation point: {:.4f}'.format(xSepVal))
        print('Initial reattachment point: {:.4f}'.format(xSepVal + Lb))

        # solve for actual reattachment length
        # root_sol = root_scalar(sepLengthRoot,
        #                        method='secant',
        #                        x0=(xSepVal + Lb),
        #                        x1=(xSepVal + Lb + 1e-6))
        root_sol = root_scalar(sepLengthRoot,
                               method='toms748',
                               bracket=[xSepVal + 0.1, self.xgrid[-5]])
        xr_fin = root_sol.root

        print('\nSeparation Bubble Solve:')
        print('\tInitial Separation Length: {:}'.format(Lb))
        print('\tConverged: {:}'.format(root_sol.converged))
        print('\tReason: {:}'.format(root_sol.flag))
        print('\tSeparation Length: {:}'.format(xr_fin - xSepVal))
        print('\tReattachment point: {:}'.format(xr_fin))

        # fake it
        # xr_fin = xSepVal + Lb
        # solve for the final polynomial coefficients
        zr, dzr = interp_end_vals(xr_fin)
        poly = solve_coef(xSepVal, xr_fin, zd, dzd, zr, dzr)

        # set the envelope
        xrInd = np.argwhere(self.xgrid > xr_fin)[0][0] - 1
        bub = lambda x: poly[0] * x ** 3 + poly[1] * x ** 2 + poly[2] * x + poly[3]

        self.env = self.z.copy()
        self.env[xSepInd:xrInd + 1] = bub(self.xgrid[xSepInd:xrInd + 1])

def compute_shear_perturbation_1D(x, z, tau, params):
    '''Compute wind shear perturbation for given free-flow wind
    speed on computational grid. based on same implementation in Duna'''

    ny, nx = x.shape
    x = x[0,:]
    dx = x[1] - x[0]
    dx = np.abs(dx)
    alfa = params['alfa']
    beta = params['beta']
    nx = nx-1
    tau_update = tau.copy()


    for igridy in range(ny):
        zb = z[igridy, :]
        dzbdx = np.zeros(x.size)

        dzbdx[1:-1] = (zb[2:] - zb[0:-2]) / 2 / dx
        tau_over_tau0 = np.zeros(x.size)

        for i in range(nx + 1):
            integ = 0
            startval = i - nx
            endval = i - 1
            for j in np.arange(startval, endval + 1):
                if j != 0:
                    integ = integ + dzbdx[i - j] / (j * np.pi)
            tau_over_tau0[i] = alfa * (integ + beta * dzbdx[i]) + 1
            tau_over_tau0[i] = np.maximum(tau_over_tau0[i], 0.1)

        tau_update[igridy,:] = tau[igridy,:] * tau_over_tau0

    return tau_update

def separation_quasi2d(x,y,z,tau,dx,dy, params):

    if params['process_separation']:

        c = params['c']
        mu_b = params['mu_b']
        sep_filter_iterations = params['sep_filter_iterations']
        zsep_y_filter = params['zsep_y_filter']

        # Initialize grid and bed dimensions
        nx = len(z[1])
        ny = len(z[0])
        dx = dx
        dy = dy
        dzx = np.zeros(z.shape)
        dzdx0 = np.zeros(z.shape)
        dzdx1 = np.zeros(z.shape)
        stall = np.zeros(z.shape)
        bubble = np.zeros(z.shape)
        k = np.array(range(0, nx))

        zsep = np.zeros(z.shape)  # total separation bubble
        zsep_new = np.zeros(z.shape)  # first-oder separation bubble surface

        zfft = np.zeros((ny, nx), dtype=complex)

        # Compute bed slope angle in x-dir
        dzx[:, :-2] = np.rad2deg(np.arctan((z[:, 2:] - z[:, :-2]) / (2. * dx)))
        dzx[:, -2] = dzx[:, -3]
        dzx[:, -1] = dzx[:, -2]

        # Determine location of separation bubbles
        '''Separation bubble exist if bed slope angle (lee side) 
        is larger than max angle that wind stream lines can 
        follow behind an obstacle (mu_b = ..)'''

        stall += np.logical_and(abs(dzx) > mu_b, dzx < 0.)

        stall[:, 1:-1] += np.logical_and(stall[:, 1:-1] == 0, stall[:, :-2] > 0., stall[:, 2:] > 0.)

        # Define separation bubble
        bubble[:, :-1] = (stall[:, :-1] == 0.) * (stall[:, 1:] > 0.)

        # Better solution for cleaner separation bubble, but no working Barchan dune (yet)
        p = 1
        bubble[:, p:] = bubble[:, :-p]
        bubble[:, -p:] = 0

        bubble = bubble.astype(int)

        # Count separation bubbles
        n = np.sum(bubble)
        bubble_n = np.asarray(np.where(bubble == True)).T

        # Walk through all separation bubbles and determine polynoms
        j = 9999
        for k in range(0, n):

            i = bubble_n[k, 1]
            j = bubble_n[k, 0]

            # Bart: check for negative wind direction
            #if np.sum(gc['taux']) >= 0:
            idir = 1
            #else:
            #    idir = -1

            ix_neg = (dzx[j, i + idir * 5:] >= 0)  # i + 5??

            if np.sum(ix_neg) == 0:
                zbrink = z[j, i]  # z level of brink at z(x0)
            else:
                zbrink = z[j, i] - z[j, i + idir * 5 + idir * np.where(ix_neg)[0][0]]

            # Better solution and cleaner separation bubble, but no working Barchan dune (yet)
            dzdx0 = (z[j, i] - z[j, i - 3]) / (3. * dx)

            a = dzdx0 / c
            ls = np.minimum(np.maximum((3. * zbrink / (2. * c) * (1. + a / 4. + a ** 2 / 8.)), 0.1), 200.)

            a2 = -3 * zbrink / ls ** 2 - 2 * dzdx0 / ls
            a3 = 2 * zbrink / ls ** 3 + dzdx0 / ls ** 2

            i_max = min(i + int(ls / dx) + 1, int(nx - 1))

            if idir == 1:
                xs = x[j, i:i_max] - x[j, i]
            else:
                xs = -(x[j, i:i_max] - x[j, i])

            zsep_new[j, i:i_max] = (a3 * xs ** 3 + a2 * xs ** 2 + dzdx0 * xs + z[j, i])

            # Choose maximum of bedlevel, previous zseps and new zseps
            zsep[j, :] = np.maximum.reduce([z[j, :], zsep[j, :], zsep_new[j, :]])

            for filter_iter in range(sep_filter_iterations):

                zsep_new = np.zeros(zsep.shape)

                Cut = 1.5
                dk = 2.0 * np.pi / (np.max(x))
                zfft[j, :] = np.fft.fft(zsep[j, :])
                zfft[j, :] *= np.exp(-(dk * k * dx) ** 2 / (2. * Cut ** 2))
                zsep_fft = np.real(np.fft.ifft(zfft[j, :]))

                if np.sum(ix_neg) == 0:
                    zbrink = zsep_fft[i]
                else:
                    zbrink = zsep_fft[i] - zsep_fft[i + idir * 5 + idir * np.where(ix_neg)[0][0]]

                # First order polynom
                dzdx1 = (zsep_fft[i] - zsep_fft[i - 3]) / (3. * dx)

                a = dzdx1 / c

                ls = np.minimum(np.maximum((3. * zbrink / (2. * c) * (1. + a / 4. + a ** 2 / 8.)), 0.1), 200.)

                a2 = -3 * zbrink / ls ** 2 - 2 * dzdx1 / ls
                a3 = 2 * zbrink / ls ** 3 + dzdx1 / ls ** 2

                i_max1 = min(i + idir * int(ls / dx), int(nx - 1))

                if idir == 1:
                    xs1 = x[j, i:i_max1] - x[j, i]
                else:
                    xs1 = -(x[j, i:i_max1] - x[j, i])

                zsep_new[j, i:i_max1] = (a3 * xs1 ** 3 + a2 * xs1 ** 2 + dzdx1 * xs1 + zbrink)

                # Pick the maximum seperation bubble hieght at all locations
                zsep[j, :] = np.maximum.reduce([z[j, :], zsep[j, :], zsep_new[j, :]])

        # Smooth surface of separation bubbles over y direction
        if zsep_y_filter:
            zsep = ndimage.gaussian_filter1d(zsep, sigma=0.2, axis=0)

        # Correct for any seperation bubbles that are below the bed surface following smoothing
        ilow = zsep < z
        zsep[ilow] = z[ilow]

        hsep = zsep - z
        tau_sep = params['tau_sep']
        slope = params['slope']
        delta = 1. / (slope * tau_sep)

        zsepdelta = np.minimum(np.maximum(1. - delta * hsep, 0.), 1.)

        tau = tau * zsepdelta

    else:

        zsep = np.zeros(z.shape)
        hsep = np.zeros(z.shape)
        tau = tau

    return zsep, hsep, tau

def compute_shear_perturbation_kroy1D(x, y, z, tau, params):
    '''Compute wind shear perturbation for given free-flow wind
    speed on computational grid. based on same implementation in Duna'''

    ny, nx = x.shape
    x = x[0,:]
    dx = x[1] - x[0]
    dx = np.abs(dx)

    x2 = x[1:]-dx/2

    L = params['L']
    l = params['l']
    kappa = params['kappa']
    mu = params['rho']*params['viscNu']
    rho = params['rho']
    sep_angle = params['sep_angle'] #14
    z0 = params['z0']

    tau_update = tau.copy()
    ix = np.isnan(tau_update)
    tau_update[ix] = 0

    for igridy in range(ny):
        zb = z[igridy, :]
        tau_in = tau[igridy, :]
        kroy = Kroy2002(x, zb, z0, L, rho, mu, tau_in, sep=False, sepInd=None, sep_ang=sep_angle, kappa=kappa)
        kroy.calc_tau_p()
        #tau_kroy = kroy.tau_tot
        tau_p_kroy = kroy.tau_p
        try:
            tau_p_interp = np.interp(x, x2, tau_p_kroy)
        except:
            tau_p_interp = tau_p_kroy

        ix = np.isnan(tau_p_interp)
        tau_p_interp[ix] = 0
        tau_kroy = tau_in * (1 + tau_p_interp)

        ix = np.isinf(tau_kroy)
        tau_kroy[ix] = tau_in[ix]

        tau_update[igridy, :] = tau_kroy

    #fix bad values
    ix = np.isinf(tau_update)
    tau_update[ix] = tau[ix]
    ix = np.isnan(tau_update)
    tau_update[ix] = tau[ix]
    ix = tau_update < 0
    tau_update[ix] = tau[ix]

    xdata = x
    ydata = y
    tau_data = tau_update

    return tau_update