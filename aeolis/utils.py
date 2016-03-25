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
    fmt.append('%s\n         ' % msg)

    i = 0
    for k, v in sorted(props.iteritems()):
        k = re.sub('^min', 'min. ', k)
        k = re.sub('^max', 'max. ', k)
        k = re.sub('^nr', '# ', k)
    
        fmt.append('%-15s: %-10s ' % (k.ljust(15, '.'),
                                      print_value(v)))
        i += 1

        if i == ncolumns:
            fmt.append('\n         ')
            i = 0
            
    fmt.append('\n')
            
    return ''.join(fmt)
    
