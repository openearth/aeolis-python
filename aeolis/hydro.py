import numpy as np

# package modules
from utils import *


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

    if p['tide_file'] is not None:
    
        s['zs'][:,:] = interp_circular(t,
                                       p['tide_file'][:,0],
                                       p['tide_file'][:,1])

        # ensure compatibility with XBeach: zs >= zb
        s['zs'] = np.maximum(s['zs'], s['zb'])

    if p['wave_file'] is not None:
    
        s['Hs'][:,:] = interp_circular(t,
                                       p['wave_file'][:,0],
                                       p['wave_file'][:,1])

        # maximize wave height by depth ratio ``gamma``
        s['Hs'] = np.minimum((s['zs'] - s['zb']) * p['gamma'], s['Hs'])
        
    if p['meteo_file'] is not None:

        m = interp_array(t,
                         p['meteo_file'][:,0],
                         p['meteo_file'][:,1:], circular=True)

        s['meteo'] = dict(zip(('S','T','R','s','p','l') , m))

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

    F = -np.log(.05) / p['Tdry']
    
    # infiltration using Darcy
    ix = s['zs'] - s['zb'] > p['eps']
    s['moist'][ ix,:] = p['porosity']
    s['moist'][~ix,:] *= np.exp(-F * dt)

    # evaporation using Penman
    if s.has_key('meteo'):

        m = s['meteo']

        # conversion from J/m2 to MJ/m2/day
        rad = m['S'] / 1e6 / dt * 3600. * 24.
        
        m = vaporation_pressure_slope(m['T']) # [kPa/K]
        delta = saturation_pressure(m['T']) * (1. - m['R']) # [kPa]
        gamma = (m['s'] * m['p']) / (.622 * m['l']) # [kPa/K]
   
        evo = np.maximum(0., (m * rad \
                              + gamma * 6.43 * (1. + 0.536 * s['uw']) * delta) \
                         / (m['l'] * (m + gamma)))

        # conversion from mm/day to m/s
        evo = evo / 24. / 3600. / 1000.

        s['moist'] = np.maximum(0., s['moist'] - evo * dt / p['layer_thickness'])

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
