import numpy as np


def isiterable(x):
    '''Check if variable is iterable'''
    
    if type(x) is str:
        return False
    try:
        _ = [i for i in x]
    except:
        return False
    return True

                
def makeiterable(x):
    '''Ensure that variable is iterable'''
    
    if not isiterable(x):
        x = np.asarray([x])
    return x


def isarray(x):
    '''Check if variable is an array'''
    
    if type(x) is str:
        return False
    if hasattr(x, '__getitem__'):
        return True
    else:
        return False


def interp_array(x, xp, fp, **kwargs):
    '''Interpolate multiple time series at once

    Parameters
    ----------
    x : array_like
        The x-coordinates of the interpolated values.
    xp : 1-D sequence of floats
        The x-coordinates of the data points, must be increasing.
    fp : 2-D sequence of floats
        The y-coordinates of the data points, same length as ``xp``.
    kwargs : dict
        Keyword options to the numpy.interp function

    Returns
    -------
    ndarray
        The interpolated values, same length as second dimension of ``fp``.

    '''
    
    f = np.zeros((1,fp.shape[1]))
    for i in range(fp.shape[1]):
        f[i] = np.interp(x, xp, fp[:,i], **kwargs)
    return f


def normalize(x, axis=0, fill=0.):
    '''Normalize array

    Normalizes an array to make it sum to unity over a specific
    axis. The procedure is safe for dimensions that sum to zero. These
    dimensions return the ``fill`` value instead.

    Parameters
    ----------
    x : array_like
        The array to be normalized
    axis : int, optional
        The normalization axis (default: 0)
    fill : float, optional
        The return value for all-zero dimensions (default: 0.)

    '''

    s = np.sum(x, axis=axis, keepdims=True).repeat(x.shape[axis], axis=axis)
    ix = s != 0.
    y = np.zeros(x.shape) + fill
    y[ix] = x[ix] / s[ix]
    return y
