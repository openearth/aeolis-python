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


import re
import numpy as np


def isiterable(x):
    '''Check if variable is iterable'''

    if isinstance(x, str):
        return False
    try:
        _ = [i for i in x]
    except:
        return False
    return True

                
def makeiterable(x):
    '''Ensure that variable is iterable'''
    
    if not isiterable(x):
        if x is None:
            x = np.asarray([])
        else:
            x = np.asarray([x])
    return np.asarray(x)


def isarray(x):
    '''Check if variable is an array'''
    
    if isinstance(x, str):
        return False
    if hasattr(x, '__getitem__'):
        return True
    else:
        return False


def interp_array(x, xp, fp, circular=False, **kwargs):
    '''Interpolate multiple time series at once

    Parameters
    ----------
    x : array_like
        The x-coordinates of the interpolated values.
    xp : 1-D sequence of floats
        The x-coordinates of the data points, must be increasing.
    fp : 2-D sequence of floats
        The y-coordinates of the data points, same length as ``xp``.
    circular : bool
        Use the :func:`interp_circular` function rather than the 
        :func:`numpy.interp` function.
    kwargs : dict
        Keyword options to the :func:`numpy.interp` function

    Returns
    -------
    ndarray
        The interpolated values, same length as second dimension of ``fp``.

    '''

    f = np.zeros((fp.shape[1],))
    for i in range(fp.shape[1]):
        if circular:
            f[i] = interp_circular(x, xp, fp[:,i], **kwargs)
        else:
            f[i] = np.interp(x, xp, fp[:,i], **kwargs)
    return f


def interp_circular(x, xp, fp, **kwargs):
    '''One-dimensional linear interpolation.

    Returns the one-dimensional piecewise linear interpolant to a
    function with given values at discrete data-points. Values beyond
    the limits of ``x`` are interpolated in circular manner. For
    example, a value of ``x > x.max()`` evaluates as ``f(x-x.max())``
    assuming that ``x.max() - x < x.max()``.

    Parameters
    ----------
    x : array_like
        The x-coordinates of the interpolated values.
    xp : 1-D sequence of floats
        The x-coordinates of the data points, must be increasing.
    fp : 1-D sequence of floats
        The y-coordinates of the data points, same length as ``xp``.
    kwargs : dict
        Keyword options to the :func:`numpy.interp` function

    Returns
    -------
    y : {float, ndarray}
        The interpolated values, same shape as ``x``.

    Raises
    ------
    ValueError
        If ``xp`` and ``fp`` have different length

    '''
    
    xmin = xp.min()
    xmax = xp.max()
    xrng = xmax - xmin
    
    x = xmin + np.mod(x - xmax - 1., xrng + 1.)
    return np.interp(x, xp, fp, **kwargs)


def normalize(x, ref=None, axis=0, fill=0.):
    '''Normalize array

    Normalizes an array to make it sum to unity over a specific
    axis. The procedure is safe for dimensions that sum to zero. These
    dimensions return the ``fill`` value instead.

    Parameters
    ----------
    x : array_like
        The array to be normalized
    ref : array_like, optional
        Alternative normalization reference, if not specified, the sum of x is used
    axis : int, optional
        The normalization axis (default: 0)
    fill : float, optional
        The return value for all-zero dimensions (default: 0.)

    '''

    x = makeiterable(x)
    if ref is None:
        ref = np.sum(x, axis=axis, keepdims=True).repeat(x.shape[axis], axis=axis)
    ix = ref != 0.
    y = np.zeros(x.shape) + fill
    y[ix] = x[ix] / ref[ix]
    return y


def prevent_tiny_negatives(x, max_error=1e-10, replacement=0.):
    '''Replace tiny negative values in array
    
    Parameters
    ----------
    x : np.ndarray
        Array with potential tiny negative values
    max_error : float
        Maximum absolute value to be replaced
    replacement : float
        Replacement value
        
    Returns
    -------
    np.ndarray
        Array with tiny negative values removed
        
    '''
    
    ix = (x < 0.) & (x > -max_error)
    x[ix] = replacement
    
    return x

                           
def print_value(val, fill='<novalue>'):
    '''Construct a string representation from an arbitrary value

    Parameters
    ----------
    val : misc
        Value to be represented as string
    fill : str, optional
        String representation used in case no value is given

    Returns
    -------
    str
        String representation of value

    '''

    if isiterable(val):
        return ' '.join([print_value(x) for x in val])
    elif val is None:
        return fill
    elif isinstance(val, bool):
        return 'T' if val else 'F'
    elif isinstance(val, int):
        return '%d' % val
    elif isinstance(val, float):
        if val < 1.:
            return '%0.6f' % val
        else:
            return '%0.2f' % val
    else:
        return str(val)


def format_log(msg, ncolumns=2, **props):
    '''Format log message into columns
    
    Prints log message and additional data into a column format
    that fits into a 70 character terminal.
    
    Parameters
    ----------
    msg : str
        Main log message
    ncolumns : int
        Number of columns
    props : key/value pairs
        Properties to print in column format
        
    Returns
    -------
    str
        Formatted log message
        
    Note
    ----
    Properties names starting with ``min``, ``max`` or ``nr`` are
    respectively replaced by ``min.``, ``max.`` or ``#``.

    '''
            
    fmt = []
    fmt.append(msg)

    i = 0
    fmt.append('')
    for k, v in sorted(props.items()):
        
        if i == ncolumns:
            fmt.append('')
            i = 0
            
        k = re.sub('^min', 'min. ', k)
        k = re.sub('^max', 'max. ', k)
        k = re.sub('^nr', '# ', k)
    
        fmt[-1] += '%-15s: %-10s ' % (k.ljust(15, '.'),
                                      print_value(v))
        i += 1
            
    return fmt
    

def apply_mask(arr, mask):
    '''Apply complex mask

    The real part of the complex mask is multiplied with the input
    array. Subsequently the imaginary part is added and the result
    returned.

    The shape of the mask is assumed to match the first few dimensions
    of the input array. If the input array is larger than the mask,
    the mask is repeated for any additional dimensions.

    Parameters
    ----------
    arr : numpy.ndarray
        Array or matrix to which the mask needs to be applied
    mask : numpy.ndarray
        Array or matrix with complex mask values

    Returns
    -------
    arr : numpy.ndarray
        Array or matrix to which the mask is applied

    '''

    # repeat mask to match shape of input array
    mask = np.asarray(mask)
    shp = arr.shape[mask.ndim:]
    for d in shp:
        mask = mask.reshape(mask.shape + tuple([1])).repeat(d, axis=-1)

    # apply mask
    arr *= np.real(mask)
    arr += np.imag(mask)

    return arr
