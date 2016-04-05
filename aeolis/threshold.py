import logging
import numpy as np


def compute(s, p):
    '''Compute wind velocity threshold based on bed surface properties

    Computes wind velocity threshold based on grain size fractions,
    bed slope, soil moisture content, air humidity and the presence of
    roughness elements. All bed surface properties increase the
    current wind velocity threshold, except for the grain size
    fractions. Therefore, the computation is initialized by the grain
    size fractions and subsequently altered by the other bed surface
    properties.

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

    See Also
    --------
    :func:`~aeolis.threshold.compute_grainsize`
    :func:`~aeolis.threshold.compute_bedslope`
    :func:`~aeolis.threshold.compute_moisture`
    :func:`~aeolis.threshold.compute_humidity`
    :func:`~aeolis.threshold.compute_roughness`

    '''

    if p['th_grainsize']:
        s = compute_grainsize(s, p)
    if p['th_bedslope']:
        s = compute_bedslope(s, p)
    if p['th_moisture']:
        s = compute_moisture(s, p)
    if p['th_humidity']:
        s = compute_humidity(s, p)
    if p['th_roughness']:
        s = compute_roughness(s, p)

    return s


def compute_grainsize(s, p):
    '''Compute wind velocity threshold based on grain size fractions following Bagnold (1937)

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

    s['uth'][:,:,:] = 1.
    s['uth'][:,:,:] *= p['A'] * np.sqrt(((p['rhop'] - p['rhoa']) * p['g'] * p['grain_size']) / p['rhop'])
    return s


def compute_bedslope(s, p):
    '''Modify wind velocity threshold based on bed slopes following Dyer (1986)

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

    return s


def compute_moisture(s, p):
    '''Modify wind velocity threshold based on soil moisture content following Belly (1964) or Hotta (1984)

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

    nf = p['nfractions']
    
    # convert from volumetric content (percentage of volume) to
    # geotechnical mass content (percentage of dry mass)
    mg = (s['moist'][:,:,:1] * p['rhow'] / (p['rhop'] * (1. - p['porosity']))).repeat(nf, axis=2)
    ix = mg > 0.005
    
    if p['method_moist'] == 'belly_johnson':
        s['uth'][ix] *= np.maximum(1., 1.8+0.6*np.log10(mg[ix]))
    elif p['method_moist'] == 'hotta':
        s['uth'][ix] += 7.5 * mg[ix]
    else:
        raise ValuerError('Unknown moisture formulation [%s]' % p['method_moist'])

    # should be .04 according to Pye and Tsoar
    # should be .64 according to Delgado-Fernandez (10% vol.)
    ix = mg > 0.064
    s['uth'][ix] = np.inf
    
    return s


def compute_humidity(s, p):
    '''Modify wind velocity threshold based on air humidity following Arens (1996)

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

    return s


def compute_roughness(s, p):
    '''Modify wind velocity threshold based on the presence of roughness elements following Raupach (1993)

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

    # TODO: these are shape dependent constants and should be
    # configurable
    m = 1.
    sigma = 3.
    beta = 45.

    # TODO: now only the largest fraction is taken into account and
    # assumed to be non-eridible, this should be configurable
    lmb = s['mass'][:,:,0,-1:] / s['mass'][:,:,0,:].sum(axis=-1, keepdims=True)
    lmb = lmb.repeat(p['nfractions'], axis=2)
    
    s['uth'] *= np.sqrt((1. - m * sigma * lmb) * (1 + m * beta * lmb))
    
    return s


