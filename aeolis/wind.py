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


import numpy as np
import operator

# package modules
from utils import *


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
        
    uw_t = p['wind_file'][:,0]
    uw_s = p['wind_file'][:,1]
    uw_d = p['wind_file'][:,2] / 180. * np.pi

    s['uw'][:,:] = interp_circular(t, uw_t, uw_s)
    s['udir'][:,:] = np.arctan2(np.interp(t, uw_t, np.sin(uw_d)),
                                np.interp(t, uw_t, np.cos(uw_d))) * 180. / np.pi

    s['uws'] = s['uw'] * np.cos(s['alfa'] + s['udir'] / 180. * np.pi)
    s['uwn'] = s['uw'] * np.sin(s['alfa'] + s['udir'] / 180. * np.pi)

    if p['ny'] == 0:
        s['uw'] = s['uws']
        s['uwn'][:,:] = 0.

    s['uw'] = np.abs(s['uw'])

    s = compute_shear(s, p)

    return s


def compute_shear(s, p):
    '''Compute shear velocity from wind velocity following Law of the Wall

    Parameters
    ----------
    uw : numpy.ndarray
        Spatial wind field
    p : dict
        Model configuration parameters

    Returns
    -------
    dict
        Spatial grids

    '''

    alpha = .174 / np.log10(p['z'] / p['k'])

    s['tau'] = alpha * s['uw']
    s['taun'] = alpha * s['uwn']
    s['taus'] = alpha * s['uws']

    return s
