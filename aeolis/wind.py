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
from numba import njit

# package modules
import aeolis.shear
from aeolis.processes.shear.spectral import SpectralShear
from aeolis.processes.shear.separation import separation_envelope, shear_reduction
from aeolis.utils import *

# FFT on the model grid (no rotated grid, no interpolation); 'spectral' is the
# name it had until 2026-10-15 and stays accepted
FFT_GRID = ('fft_grid', 'spectral')
from aeolis.separation import compute_separation, set_sep_params

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

        # # Get separation parameters based on Perturbation theory settings
        # if p['process_separation']:# and p['sep_auto_tune']:
        #     p = set_sep_params(p, s)

        if p['ny'] > 0:
            if p['method_shear'] == 'fft':
                s['shear'] = aeolis.shear.WindShear(s['x'], s['y'], s['zb'],
                                                    dx=p['dx'], dy=p['dy'],
                                                    L=p['L'], l=p['l'], z0=z0,
                                                    buffer_width=p['buffer_width'])
            elif p['method_shear'] in FFT_GRID:
                if p['method_shear'] == 'spectral':
                    logger.info('method_shear = spectral is now called fft_grid; both names work')
                # grid spacing along the columns (s) and rows (n), from the coordinates
                ds = np.hypot(s['x'][0, 1] - s['x'][0, 0], s['y'][0, 1] - s['y'][0, 0])
                dn = np.hypot(s['x'][1, 0] - s['x'][0, 0], s['y'][1, 0] - s['y'][0, 0])
                z0_eff = float(np.exp(np.mean(np.log(z0)))) if np.size(z0) > 1 else float(z0)
                if p.get('shear_detrend'):
                    logger.warning('shear_detrend was removed on 2026-10-15 and is ignored: removing the '
                                   'best-fit plane deleted the real speed-up over a beach ramp')
                # circular boundaries: the topography repeats itself, so no padding
                periodic = (p['boundary_lateral'] == 'circular',
                            p['boundary_offshore'] == 'circular' and p['boundary_onshore'] == 'circular')
                s['shear'] = SpectralShear(s['zb'].shape[0], s['zb'].shape[1], ds, dn,
                                           L=p['L'], l=p['l'], z0=z0_eff, periodic=periodic,
                                           pad_factor=p['shear_pad_factor'],
                                           filter_frac=p['shear_filter_frac'],
                                           transfer=p['shear_transfer'],
                                           short_waves=p['shear_short_waves'],
                                           short_floor=p['shear_short_floor'])
            elif p['method_shear'] != '1Dstacks':
                logger.log_and_raise('Unknown shear method [%s], use fft, fft_grid or 1Dstacks'
                                     % p['method_shear'], exc=ValueError)
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


def shear(s, p):
    """Compute shear stress and shear velocity fields."""

    # --- Shear disabled: copy initial values for Air ------------------------
    if not p['process_shear']:
        s['tauAir']    = s['tau'].copy()
        s['tausAir']   = s['taus'].copy()
        s['taunAir']   = s['taun'].copy()
        s['ustarAir']  = s['ustar'].copy()
        s['ustarsAir'] = s['ustars'].copy()
        s['ustarnAir'] = s['ustarn'].copy()
        return s

    # --- 1D and stacked-1D (independent rows) --------------------------------
    if p['ny'] == 0 or p['method_shear'] == '1Dstacks':

        # Compute separation bubble
        if p['process_separation']:
            dx = s['ds'][0, 0]
            s['zsep'] = compute_separation(p, s['zb'], dx, udir=s['udir'][0, 0])
            s['hsep'] = s['zsep'] - s['zb']
        else:
            s['zsep'] = s['zb'].copy()
            s['hsep'] = np.zeros_like(s['zb'])

        # Compute shear perturbation
        s = compute_shear1d(s, p)
        s['tauAir'] = s['tau'].copy()
        s['tausAir'] = s['taus'].copy()
        s['taunAir'] = s['taun'].copy()

        # Reduce stress within separation bubble
        if p['process_separation']:
            tau_sep, slope = 0.5, 0.2
            delta = 1. / (slope * tau_sep)
            zsepdelta = np.minimum(np.maximum(1. - delta * s['hsep'], 0.), 1.)

            s['taus'] *= zsepdelta
            s['taun'] *= zsepdelta
            s['tau']  = np.hypot(s['taus'], s['taun'])

        return stress_velocity(s, p)

    # --- Fully 2D shear (FFT-based) -----------------------------------------
    elif p['method_shear'] == 'fft':
        
        s['shear'](
            p=p, x=s['x'], y=s['y'], z=s['zb'],
            u0=s['uw'][0, 0], udir=s['udir'][0, 0],
            taux=s['taus'], tauy=s['taun'],
            taus0=s['taus0'][0, 0], taun0=s['taun0'][0, 0]
        )

        # Get shear stress and velocity
        s['taus'], s['taun'], s['tausAir'], s['taunAir'] = s['shear'].get_shear()
        s['tau'] = np.hypot(s['taus'], s['taun'])
        s['tauAir'] = np.hypot(s['tausAir'], s['taunAir'])

        # Get separation bubble
        if p['process_separation']:
            s['hsep'] = s['shear'].get_separation()
            s['zsep'] = s['hsep'] + s['zb']
        else:
            s['zsep'] = s['zb'].copy()
            s['hsep'] = np.zeros_like(s['zb'])

        s = stress_velocity(s, p)

    # --- Fully 2D shear on the model grid (fft_grid, no secondary grid) -----
    elif p['method_shear'] in FFT_GRID:

        sh = s['shear']
        taus0, taun0 = s['taus0'][0, 0], s['taun0'][0, 0]

        # separation envelope along the wind, on the model grid
        if p['process_separation'] and (taus0 != 0. or taun0 != 0.):
            s['zsep'] = separation_envelope(
                s['zb'], np.arctan2(taun0, taus0), sh.ds, sh.dn,
                look_dist=p['sep_look_dist'], k_press_up=p['sep_k_press_up'],
                k_crit_down=p['sep_k_crit_down'], s_crit=p['sep_s_crit'],
                s_leeside=p['sep_s_leeside'], periodic=sh.periodic)
        else:
            s['zsep'] = s['zb'].copy()
        s['hsep'] = s['zsep'] - s['zb']

        # perturbation of the envelope; the air stress is not reduced in the bubble
        s['taus'], s['taun'] = sh.stress(s['zsep'], taus0, taun0)
        s['tausAir'] = s['taus'].copy()
        s['taunAir'] = s['taun'].copy()
        s['tauAir'] = np.hypot(s['tausAir'], s['taunAir'])
        if p['process_separation']:
            fac = shear_reduction(s['hsep'])
            s['taus'] *= fac
            s['taun'] *= fac
        s['tau'] = np.hypot(s['taus'], s['taun'])

        s = stress_velocity(s, p)

    else:
        logger.log_and_raise('Unknown shear method: %s' 
                             % p['method_shear'], exc=ValueError)

    return s


@njit(cache=True)
def _tau_factor_stacked(zb, dx, alfa, beta):
    """
    Kroy/Duna perturbation factor tau/tau0 along x, applied per row (transect).
    """
    ny, nx = zb.shape

    # Compute bed slope
    dzbdx = np.zeros((ny, nx))
    dzbdx[:, 1:-1] = (zb[:, 2:] - zb[:, :-2]) / (2.0 * dx)

    # Compute perturbation factor
    fac = np.zeros((ny, nx))
    nxm1 = nx - 1

    for i in range(nxm1 + 1):
        # integral term per transect
        integ = np.zeros(ny)
        for j in range(i - nxm1, i):
            if j != 0:
                integ += dzbdx[:, i - j] / (j * np.pi)

        fac[:, i] = alfa * (integ + beta * dzbdx[:, i]) + 1.0

        # lower bound
        for k in range(ny):
            if fac[k, i] < 0.1:
                fac[k, i] = 0.1

    return fac


def compute_shear1d(s, p):
    """
    1D shear perturbation (and stacked-1D over rows if ny>0).

    Uses zsep as an optional separation envelope: zb_eff = max(zb, zsep).
    Does not apply the separation *stress reduction*; that stays in shear().
    """

    tau0  = s['tau'].copy()
    taus0 = s['taus'].copy()
    taun0 = s['taun'].copy()

    # Unit direction of stress to reconstruct vector after updating magnitude
    ets = np.zeros_like(tau0)
    etn = np.zeros_like(tau0)
    ix = tau0 != 0.0
    ets[ix] = taus0[ix] / tau0[ix]
    etn[ix] = taun0[ix] / tau0[ix]

    # Effective bed used for perturbation
    zb = s['zb']
    if ('zsep' in s) and (s['zsep'] is not None):
        zb = np.maximum(zb, s['zsep'])

    # Wind sign handling (kept consistent with existing convention)
    flip = np.sum(taus0) < 0.0
    zb_in = np.flip(zb, axis=1) if flip else zb

    # Compute perturbation factor (Kroy et al. 2002 coefficients; fixed for now)
    shear_alfa, shear_beta = 3., 1.
    fac = _tau_factor_stacked(zb_in, s['ds'][0, 0], shear_alfa, shear_beta)
    if flip:
        fac = np.flip(fac, axis=1)

    # Update shear stress and velocity
    s['tau']  = tau0 * fac
    s['taus'] = s['tau'] * ets
    s['taun'] = s['tau'] * etn

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

    # --- Shear for airborne sediment ----------------------------------------
    s['ustarAir'] = np.sqrt(s['tauAir'] / p['rhoa'])

    ix = s['tauAir'] > 0.
    s['ustarsAir'][ix] = s['ustarAir'][ix] * s['tausAir'][ix] / s['tauAir'][ix]
    s['ustarnAir'][ix] = s['ustarAir'][ix] * s['taunAir'][ix] / s['tauAir'][ix]

    ix = s['tauAir'] == 0.
    s['ustarAir'][ix] = 0.
    s['ustarsAir'][ix] = 0.
    s['ustarnAir'][ix] = 0.

    return s