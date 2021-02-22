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


def interpolate(s, p, t):
    '''Interpolate hydrodynamic and meteorological conditions to current time step

    Interpolates the hydrodynamic and meteorological time series to
    the current time step, if available. Meteorological parameters are
    stored as dictionary rather than a single value.

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

    if p['process_tide']:
        if p['tide_file'] is not None:
            s['zs'][:,:] = interp_circular(t, 
                                           p['tide_file'][:,0], 
                                           p['tide_file'][:,1])
        else:
            s['zs'][:,:] = 0.

        # apply complex mask
        s['zs'] = apply_mask(s['zs'], s['tide_mask'])

    if p['process_wave'] and p['wave_file'] is not None:

        # determine water depth
        h = np.maximum(0., s['zs'] - s['zb'])
    
        s['Hs'][:,:] = interp_circular(t,
                                       p['wave_file'][:,0],
                                       p['wave_file'][:,1])

        # apply complex mask
        s['Hs'] = apply_mask(s['Hs'], s['wave_mask'])

        # add wave runup
        if p['process_runup']:
            ix = s['Hs'] > 0.
            R = p['xi'] * s['Hs']
            s['zs'][ix] += R[ix] * (1. - np.minimum(1., h[ix] * p['gamma'] / s['Hs'][ix]))
        
        # maximize wave height by depth ratio ``gamma``
        s['Hs'] = np.minimum(h * p['gamma'], s['Hs'])
        
    if p['process_meteo'] and  p['meteo_file'] is not None:

        m = interp_array(t,
                         p['meteo_file'][:,0],
                         p['meteo_file'][:,1:], circular=True)

        # Symbols according to KNMI files: T = temperature [oC], RH =
        # precipitation [mm/hr], U = relative humidity [%], Q = global
        # radiation [J/m2], P = air pressure [kPa]
        s['meteo'] = dict(zip(('T','RH','U','Q','P') , m))

    # ensure compatibility with XBeach: zs >= zb
    s['zs'] = np.maximum(s['zs'], s['zb'])

    return s


def update(s, p, dt):
    '''Update soil moisture content

    Updates soil moisture content in all cells. The soil moisure
    content in flooded cells is set to the porosity. The soil moisture
    content in non-flooded cells is decreased by simulating
    infiltration using an exponential decay function with a rate ``F``
    following:

    .. math::

       M = M \\cdot e^{-F \\cdot \\Delta t}

    Cells are considered flooded if the water depth is larger than
    ``eps``. Aditionally, is meteorological conditions are provided
    the soil moisture content is decreased by simulating evaporation
    following the Penman equation:

    .. math::

        M = M - \\frac{m \\cdot R + \\gamma \cdot 6.43 \cdot (1 + 0.536 \cdot u]) \cdot \\delta}{\\lambda \cdot (m + \\gamma)}

    Parameters
    ----------
    s : dict
        Spatial grids
    p : dict
        Model configuration parameters
    dt : float
        Current time step

    Returns
    -------
    dict
        Spatial grids

    '''

    # infiltration using Darcy
    if p['process_moist']:
        F1 = -np.log(.5) / p['Tdry']
        ix = s['zs'] - s['zb'] > p['eps']
        s['moist'][ ix,0] = p['porosity']
        s['moist'][~ix,0] *= np.exp(-F1 * dt)

    # salinitation
    if p['process_salt']:
        F2 = -np.log(.5) / p['Tsalt']
        s['salt'][ ix,0] = 1.
        s['salt'][~ix,0] *= np.exp(-F2 * dt)

    # evaporation using Penman
    if 'meteo' in s:

        met = s['meteo']
        
        # convert radiation from J/m2 to MJ/m2/day
        rad = met['Q'] / 1e6 / dt * 3600. * 24.

        l = 2.26
        m = vaporation_pressure_slope(met['T']) # [kPa/K]
        delta = saturation_pressure(met['T']) * (1. - met['RH']) # [kPa]
        gamma = (p['cpair'] * met['P']) / (.622 * l) # [kPa/K]
        u2 = .174 / np.log10(p['z'] / 2.) * s['uw'] # [m/s]
   
        evo = np.maximum(0., (m * rad \
                              + gamma * 6.43 * (1. + 0.536 * u2) * delta) \
                         / (l * (m + gamma)))

        # convert evaporation from mm/day to m/s
        evo = evo / 24. / 3600. / 1000.

        # convert precipitation from mm/hr to m/s
        pcp = met['RH'] / 3600. / 1000.

        # update moisture and salt content
        s['moist'][:,:,0] = np.maximum(0.,
                                       np.minimum(p['porosity'],
                                                  s['moist'][:,:,0] \
                                                  + (pcp - evo) * dt \
                                                  / p['layer_thickness'] / p['nlayers']))
        s['salt'][:,:,0] = np.minimum(1., s['salt'][:,:,0] + pcp * dt / p['layer_thickness'])

    return s


def vaporation_pressure_slope(T):
    '''Compute vaporation pressure slope based on air temperature

    Parameters
    ----------
    T : float
        Air temperature in degrees Celcius

    Returns
    -------
    float
        Vaporation pressure slope

    '''
    
    # Tetens, 1930; Murray, 1967
    s = 4098. * 0.6108 * np.exp((17.27 * T) / (T - 237.3)) / (T + 237.3)**2 # [kPa/K]

    return s


def saturation_pressure(T):
    '''Compute saturation pressure based on air temperature

    Parameters
    ----------
    T : float
        Air temperature in degrees Celcius

    Returns
    -------
    float
        Saturation pressure

    '''

    TK = T + 273.15
    A  = -1.88e4
    B  = -13.1
    C  = -1.5e-2
    D  =  8e-7
    E  = -1.69e-11
    F  =  6.456
    vp = np.exp(A/TK + B + C*TK + D*TK**2 + E*TK**3 + F*np.log(TK)) # [kPa]

    return vp
