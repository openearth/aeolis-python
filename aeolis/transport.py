import numpy as np

# package modules
from utils import *


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
    
    nf = p['nfractions']
    uw = s['uw'][:,:,np.newaxis].repeat(nf, axis=2)
    ix = uw != 0.
    
    alpha = (0.174 / np.log10(p['z0'] / p['k']))**3

    s['Cu'] = np.zeros(uw.shape)
    s['Cu'][ix] = np.maximum(0., alpha * p['Cb'] * p['rhoa'] / p['g'] \
                             * (uw[ix] - s['uth'][ix])**3 / uw[ix])

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

    ix = s['Cu'] > 0.
    w_air = np.zeros(s['Cu'].shape)
    w_air[ix] = s['Ct'][ix] / s['Cu'][ix]

    w_bed = s['mass'][:,:,0,:] / np.sum(s['mass'][:,:,0,:], axis=2, keepdims=True)

    w = (1. - p['bi']) * w_air + w_bed * \
        (1. - np.minimum(1., (1. - p['bi']) * np.sum(w_air, axis=2, keepdims=True)))
    w = normalize(w, axis=2)
    
    return w


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
