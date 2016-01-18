import numpy as np
import operator
import matplotlib.pyplot as plt

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

    s['uw'][:,:] = np.interp(t, uw_t, uw_s)
    s['udir'][:,:] = np.arctan2(np.interp(t, uw_t, np.sin(uw_d)),
                                np.interp(t, uw_t, np.cos(uw_d))) * 180. / np.pi

    s['uws'] = s['uw'] * np.cos(s['alfa'] + s['udir'] / 180. * np.pi)
    s['uwn'] = s['uw'] * np.sin(s['alfa'] + s['udir'] / 180. * np.pi)

    if p['ny'] == 0:
       s['uw'] = s['uws']
       s['uwn'][:,:] = 0.

    s['uw'] = np.abs(s['uw'])

    return s
