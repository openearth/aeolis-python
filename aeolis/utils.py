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
#import numba
#import numba_scipy
import scipy
from numpy import ndarray
from numba import njit


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


def interp_array(x: ndarray, xp: ndarray, 
                 fp: ndarray, circular: bool=False, **kwargs: dict) -> ndarray:
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


def interp_circular(x: ndarray, xp: ndarray, fp: ndarray, **kwargs) -> ndarray:
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


def interp_circular_nearest(x: ndarray, xp: ndarray, fp: ndarray) -> ndarray:
    '''One-dimensional linear interpolation to nearest neighbor.

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

    x = xmin + np.mod(x - xmin, xrng)
  
    return fp[xp<=x][-1]

def normalize(x: ndarray, ref:ndarray = None, axis: int =0, fill: float =0.):
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
    axis : optional
        The normalization axis (default: 0)
    fill : optional
        The return value for all-zero dimensions (default: 0.)

    '''

    x = makeiterable(x)
    if ref is None:
        ref = np.sum(x, axis=axis, keepdims=True).repeat(x.shape[axis], axis=axis)
    ix = ref != 0.
    y = np.zeros(x.shape) + fill
    y[ix] = x[ix] / ref[ix]
    return y


def prevent_tiny_negatives(x: ndarray, max_error: float =1e-10, replacement: float =0.) -> ndarray:
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


def rotate(x, y, alpha, origin=(0,0)):
    '''Rotate a matrix over given angle around given origin'''
    
    xr = x - origin[0]
    yr = y - origin[1]
    
    a = alpha / 180. * np.pi
    
    R = np.asmatrix([[np.cos(a), -np.sin(a)],
                     [np.sin(a),  np.cos(a)]])
    
    xy = np.concatenate((xr.reshape((-1,1)), 
                         yr.reshape((-1,1))), axis=1) * R
                     
    return (np.asarray(xy[:,0].reshape(x.shape) + origin[0]),
            np.asarray(xy[:,1].reshape(y.shape) + origin[1]))

# @numba.njit 
def sc_kv(v, z):
    return scipy.special.kv(v, z)

def calc_grain_size(p, s, percent):
    '''Calculate grain size characteristics based on mass in each fraction

    Calculate grain size distribution for each cell based on weight 
    distribution over the fractions. Interpolates to the requested percentage 
    in the grain size distribution. For example, percent=50 will result 
    in calculation of the D50. Calculation is only executed for the top layer


    Parameters
    ----------
    s : dict
        Spatial grids
    p : dict
        Model configuration parameters
    percent : float
        Requested percentage in grain size dsitribution

    Returns
    -------
    array
        grain size per grid cell

    '''
    from scipy.interpolate import interp1d
    mass = s['mass'][:,:,0,:] # only select upper, surface layer because that is the relevant layer for transport
    D = np.zeros((mass.shape[0], mass.shape[1]))
    for yi in range(mass.shape[0]):
        for xi in range(mass.shape[1]):
                diameters = np.insert(p['grain_size'], 0, 0)
                cummulative_weights = np.cumsum(np.insert(mass[yi, xi,:], 0, 0))
                percentages = 100*cummulative_weights/np.max(cummulative_weights)
                f_linear = interp1d(list(percentages), diameters, fill_value='extrapolate') # get interpolation function
                
                # Retrieve grain size characteristics based on interpolation
                D[yi, xi] = f_linear(percent)
    return D

def calc_mean_grain_size(p, s):
    '''Calculate mean grain size based on mass in each fraction

    Calculate mean grain size for each cell based on weight distribution 
    over the fractions. Calculation is only executed for the top layer.


    Parameters
    ----------
    s : dict
        Spatial grids
    p : dict
        Model configuration parameters
    percent : float
        Requested percentage in grain size dsitribution

    Returns
    -------
    array
        mean grain size per grid cell

    '''
    mass = s['mass'][:,:,0,:] # only select upper, surface layer because that is the relevant layer for transport
    D_mean = np.zeros((mass.shape[0], mass.shape[1]))
    for yi in range(mass.shape[0]):
        for xi in range(mass.shape[1]):
                diameters = p['grain_size']
                weights = mass[yi, xi,:]/ np.sum(mass[yi, xi,:])
                
                # Retrieve mean grain size based on weight of mass
                D_mean[yi, xi] = np.sum(diameters*weights)
    return D_mean


@njit
def sweep(Cu, mass, dt, Ts, ds, dn, us, un):

    Ct = np.zeros(Cu.shape)
    pickup = np.zeros(Cu.shape)
    i=0
    k=0
    q=0

   
    # determine quadrants this code is currently only compatible with spacially non-varying winds
    # in future code the 4 loops could be made paralell and wind domains spatially varying.    


    if  np.all(un[:,:,0]>=0) & np.all(us[:,:,0]>=0): 
        q=1   
    elif np.all(un[:,:,0]>=0) & np.all(us[:,:,0]<0): 
        q=2
    elif np.all(un[:,:,0]<0) & np.all(us[:,:,0]<0): 
        q=3
    elif np.all(un[:,:,0]<0) & np.all(us[:,:,0]>=0): 
        q=4

    if q==0:
        raise NotImplementedError('Spatially varying input detected please choose other solver')
        # print('Model will crash now')
        # input("Press Enter to continue...")
        # crash
        

    

    # The while loop accounts for the circular boundary.
    # The zero difference on the circular boundary will likely only work in theoretical cases.
    # It works for now, it could be replaced by a max_error value instead.  
    # This first loop is valid for the first quadrant
    if q==1:
        while k==0 or np.any(np.abs(Ct[0,:]-Ct[-1,:])!=0):
            #circular boundary this loop stops if lateral boundaries are equal.
            Ct[0,:]=Ct[-1,:]  
            for n in range(1,Ct.shape[0]):
                for s in range(1,Ct.shape[1]):
                    # sweep from [0,0] corner through domain in positive direction. 
                    Ct[n,s,i] = (+ Ct[n-1,s,i] * un[n,s,0] * ds[n,s] \
                                    + Ct[n,s-1,i] * us[n,s,0] * dn[n,s] \
                                    + Cu[n,s,i] * ds[n,s] * dn [n,s] / Ts  ) \
                                    / ((us[n,s,0] * dn[n,s]) + (un[n,s,0] * ds[n,s]) + (ds[n,s] * dn [n,s] / Ts) )
                    #calculate pickup                
                    pickup[n,s,i] = (Cu[n,s,i]-Ct[n,s,i]) * dt/Ts
                    #check for supply limitations and re-iterate concentration to account for supply limitations
                    if pickup[n,s,i]>mass[n,s,0,i]:
                        pickup[n,s,i] = mass[n,s,0,i]              
                        Ct[n,s,i] = (+ Ct[n-1,s,i] * un[n,s,0] * ds[n,s] \
                                        + Ct[n,s-1,i] * us[n,s,0] * dn[n,s] \
                                        + pickup[n,s,i] * ds[n,s] * dn [n,s] / dt ) \
                                        / ((us[n,s,0] * dn[n,s]) + (un[n,s,0] * ds[n,s]))
            #count the amount of iterations
            k+=1
            #print(k)        
        
    if q==2:
        while k==0 or np.any(np.abs(Ct[0,:]-Ct[-1,:])!=0):
            #circular boundary this loop stops if lateral boundaries are equal.
            Ct[0,:]=Ct[-1,:] 
            for n in range(1,Ct.shape[0]):
                #print(n)
                for s in range(Ct.shape[1]-2,-1,-1):
                    # sweep from [-1,0] corner through domain in positive direction. 
                    Ct[n,s,i] = (+ Ct[n-1,s,i] * un[n,s,0] * ds[n,s] \
                                    - Ct[n,s+1,i] * us[n,s,0] * dn[n,s] \
                                    + Cu[n,s,i] * ds[n,s] * dn [n,s] / Ts  ) \
                                    / (( - us[n,s,0] * dn[n,s]) + (un[n,s,0] * ds[n,s]) + (ds[n,s] * dn [n,s] / Ts) )
                    #calculate pickup                
                    pickup[n,s,i] = (Cu[n,s,i]-Ct[n,s,i]) * dt/Ts
                    #check for supply limitations and re-iterate concentration to account for supply limitations
                    if pickup[n,s,i]>mass[n,s,0,i]:
                        pickup[n,s,i] = mass[n,s,0,i]              
                        Ct[n,s,i] = (+ Ct[n-1,s,i] * un[n,s,0] * ds[n,s] \
                                        - Ct[n,s+1,i] * us[n,s,0] * dn[n,s] \
                                        + pickup[n,s,i] * ds[n,s] * dn [n,s] / dt ) \
                                        / (( - us[n,s,0] * dn[n,s]) + (un[n,s,0] * ds[n,s]))
            #count the amount of iterations
            k+=1
            #print(k)

        
    if q==3:
        while k==0 or np.any(np.abs(Ct[0,:]-Ct[-1,:])!=0):
            #circular boundary this loop stops if lateral boundaries are equal.
            Ct[-1,:]=Ct[0,:] 
            for n in range(Ct.shape[0]-2,-1,-1):
                #print(n)
                for s in range(Ct.shape[1]-2,-1,-1):
                    # sweep from [-1,0] corner through domain in positive direction. 
                    Ct[n,s,i] = (- Ct[n+1,s,i] * un[n,s,0] * ds[n,s] \
                                    - Ct[n,s+1,i] * us[n,s,0] * dn[n,s] \
                                    + Cu[n,s,i] * ds[n,s] * dn [n,s] / Ts  ) \
                                    / (( - us[n,s,0] * dn[n,s]) + ( - un[n,s,0] * ds[n,s]) + (ds[n,s] * dn [n,s] / Ts) )
                    #calculate pickup                
                    pickup[n,s,i] = (Cu[n,s,i]-Ct[n,s,i]) * dt/Ts
                    #check for supply limitations and re-iterate concentration to account for supply limitations
                    if pickup[n,s,i]>mass[n,s,0,i]:
                        pickup[n,s,i] = mass[n,s,0,i]              
                        Ct[n,s,i] = (- Ct[n+1,s,i] * un[n,s,0] * ds[n,s] \
                                        - Ct[n,s+1,i] * us[n,s,0] * dn[n,s] \
                                        + pickup[n,s,i] * ds[n,s] * dn [n,s] / dt ) \
                                        / (( - us[n,s,0] * dn[n,s]) + ( - un[n,s,0] * ds[n,s]))
            #count the amount of iterations
            k+=1
            #print(k)

        
    if q==4:
        while k==0 or np.any(np.abs(Ct[0,:]-Ct[-1,:])!=0):
            #circular boundary this loop stops if lateral boundaries are equal.
            Ct[-1,:]=Ct[0,:] 
            for n in range(Ct.shape[0]-2,-1,-1):
                #print(n)
                for s in range(1,Ct.shape[1]):
                    # sweep from [-1,0] corner through domain in positive direction. 
                    Ct[n,s,i] = (- Ct[n+1,s,i] * un[n,s,0] * ds[n,s] \
                                    + Ct[n,s-1,i] * us[n,s,0] * dn[n,s] \
                                    + Cu[n,s,i] * ds[n,s] * dn [n,s] / Ts  ) \
                                    / ((us[n,s,0] * dn[n,s]) + ( - un[n,s,0] * ds[n,s]) + (ds[n,s] * dn [n,s] / Ts) )
                    #calculate pickup                
                    pickup[n,s,i] = (Cu[n,s,i]-Ct[n,s,i]) * dt/Ts
                    #check for supply limitations and re-iterate concentration to account for supply limitations
                    if pickup[n,s,i]>mass[n,s,0,i]:
                        pickup[n,s,i] = mass[n,s,0,i]              
                        Ct[n,s,i] = (- Ct[n+1,s,i] * un[n,s,0] * ds[n,s] \
                                        + Ct[n,s-1,i] * us[n,s,0] * dn[n,s] \
                                        + pickup[n,s,i] * ds[n,s] * dn [n,s] / dt ) \
                                        / ((us[n,s,0] * dn[n,s]) + (- un[n,s,0] * ds[n,s]))
            #count the amount of iterations
            k+=1
            #print(k)


    return Ct, pickup

@njit
def sweep2(Ct, Cu, mass, dt, Ts, ds, dn, us, un):

    pickup = np.zeros(Cu.shape)
    i=0
    k=0

    # Are the lateral boundary conditions circular?
    circ_lateral = False
    if Ct[0,1,0]==-1:
        circ_lateral = True
        Ct[0,:,0] = 0                
        Ct[-1,:,0] = 0

    circ_offshore = False
    if Ct[1,0,0]==-1:
        circ_offshore = True
        Ct[:,0,0] = 0                
        Ct[:,-1,0] = 0

    recirc_offshore = False
    if Ct[1,0,0]==-2:
        recirc_offshore = True
        Ct[:,0,0] = 0                
        Ct[:,-1,0] = 0
    
    
    ufs = np.zeros((np.shape(us)[0], np.shape(us)[1]+1, np.shape(us)[2]))
    ufn = np.zeros((np.shape(un)[0]+1, np.shape(un)[1], np.shape(un)[2]))
    
    # define fluxes
    ufs[:,1:-1, :] = 0.5*us[:,:-1, :] + 0.5*us[:,1:, :]
    ufn[1:-1,:, :] = 0.5*un[:-1,:, :] + 0.5*un[1:,:, :]
  
    #boundary values
    ufs[:,0, :]  = us[:,0, :]
    ufs[:,-1, :] = us[:,-1, :]
   
    ufn[0,:, :]  = un[0,:, :]
    ufn[-1,:, :] = un[-1,:, :]

    Ct_last = Ct.copy()
    while k==0 or np.any(np.abs(Ct[:,:,i]-Ct_last[:,:,i])>1e-10):
        Ct_last = Ct.copy()

        # lateral boundaries circular
        if circ_lateral:
            Ct[0,:,0],Ct[-1,:,0] = Ct[-1,:,0],Ct[0,:,0]

        if circ_offshore:
            Ct[:,0,0],Ct[:,-1,0] = Ct[:,-1,0],Ct[:,0,0]

        if recirc_offshore:
            # print(Ct[:,1,0])
            # print(Ct[:,-2,0]) 
            Ct[:,0,0],Ct[:,-1,0] = np.average(Ct[:,-2,0]),np.average(Ct[:,1,0])
            # print(Ct[:,0,0])
            # print(Ct[:,-1,0]) 

        # make an array with a bolean operator. This keeps track of considerd cells. We start with all False (not considered)
        q = np.zeros(Cu.shape[:2]) 

    ########################################################################################
        # in this sweeping algorithm we sweep over the 4 quadrants
        # assuming that most cells have no converging/divering charactersitics.
        # In the last quadrant we take converging and diverging cells into account. 

        # The First quadrant  
        for n in range(1,Ct.shape[0]-1):
            for s in range(1,Ct.shape[1]-1):
                if (not q[n,s]) and (ufn[n,s,0]>=0) and (ufs[n,s,0]>=0) and (ufn[n+1,s,0]>=0) and (ufs[n,s+1,0]>=0):
                    Ct[n,s,i] = (+ (Ct[n-1,s,i] * ufn[n,s,0] * ds[n,s]) \
                                    + (Ct[n,s-1,i] * ufs[n,s,0] * dn[n,s]) \
                                    + Cu[n,s,i] * ds[n,s] * dn [n,s] / Ts  ) \
                                    / ( + (ufn[n+1,s,0] * ds[n,s]) \
                                    + (ufs[n,s+1,0] * dn[n,s]) \
                                    + (ds[n,s] * dn [n,s] / Ts) )
                    #calculate pickup                
                    pickup[n,s,i] = (Cu[n,s,i]-Ct[n,s,i]) * dt/Ts
                    #check for supply limitations and re-iterate concentration to account for supply limitations
                    if pickup[n,s,i]>mass[n,s,0,i]:
                        pickup[n,s,i] = mass[n,s,0,i]              
                        Ct[n,s,i] = (+ (Ct[n-1,s,i] * ufn[n,s,0] * ds[n,s]) \
                                    + (Ct[n,s-1,i] * ufs[n,s,0] * dn[n,s]) \
                                        + pickup[n,s,i] * ds[n,s] * dn [n,s] / dt ) \
                                        / ( + (ufn[n+1,s,0] * ds[n,s]) \
                                        + (ufs[n,s+1,0] * dn[n,s]) \
                                        + (ds[n,s] * dn [n,s] / Ts) )
                    q[n,s]=1
        # The second quadrant
        for n in range(1,Ct.shape[0]):
            for s in range(Ct.shape[1]-2,-1,-1):  
                if (not q[n,s]) and (ufn[n,s,0]>=0) and (ufs[n,s,0]<=0) and (ufn[n+1,s,0]>=0) and (ufs[n,s+1,0]<=0):
                    Ct[n,s,i] = (+ (Ct[n-1,s,i] * ufn[n,s,0] * ds[n,s]) \
                                    + ( -Ct[n,s+1,i] * ufs[n,s+1,0] * dn[n,s]) \
                                    + Cu[n,s,i] * ds[n,s] * dn [n,s] / Ts  ) \
                                    / ( + (ufn[n+1,s,0] * ds[n,s]) \
                                    + (-ufs[n,s,0] * dn[n,s]) \
                                    + (ds[n,s] * dn [n,s] / Ts) )
                    #calculate pickup                
                    pickup[n,s,i] = (Cu[n,s,i]-Ct[n,s,i]) * dt/Ts
                    #check for supply limitations and re-iterate concentration to account for supply limitations
                    if pickup[n,s,i]>mass[n,s,0,i]:
                        pickup[n,s,i] = mass[n,s,0,i]              
                        Ct[n,s,i] = (+ (Ct[n-1,s,i] * ufn[n,s,0] * ds[n,s]) \
                                        + ( -Ct[n,s+1,i] * ufs[n,s+1,0] * dn[n,s]) \
                                        + pickup[n,s,i] * ds[n,s] * dn [n,s] / dt ) \
                                        / ( + (ufn[n+1,s,0] * ds[n,s]) \
                                        + (-ufs[n,s,0] * dn[n,s]) \
                                        + (ds[n,s] * dn [n,s] / Ts) )  
                    q[n,s]=2
        # The third quadrant
        for n in range(Ct.shape[0]-2,-1,-1):
            for s in range(Ct.shape[1]-2,-1,-1):
                if (not q[n,s]) and (ufn[n,s,0]<=0) and (ufs[n,s,0]<=0) and (ufn[n+1,s,0]<=0) and (ufs[n,s+1,0]<=0):
                    Ct[n,s,i] = (+ ( -Ct[n+1,s,i] * ufn[n+1,s,0] * dn[n,s]) \
                                    + ( -Ct[n,s+1,i] * ufs[n,s+1,0] * dn[n,s]) \
                                    + Cu[n,s,i] * ds[n,s] * dn [n,s] / Ts  ) \
                                    / ( + (-ufn[n,s,0] * dn[n,s]) \
                                    + (-ufs[n,s,0] * dn[n,s]) \
                                    + (ds[n,s] * dn [n,s] / Ts) )
                    #calculate pickup                
                    pickup[n,s,i] = (Cu[n,s,i]-Ct[n,s,i]) * dt/Ts
                    #check for supply limitations and re-iterate concentration to account for supply limitations
                    if pickup[n,s,i]>mass[n,s,0,i]:
                        pickup[n,s,i] = mass[n,s,0,i]              
                        Ct[n,s,i] = (+ ( -Ct[n+1,s,i] * ufn[n+1,s,0] * dn[n,s]) \
                                        + ( -Ct[n,s+1,i] * ufs[n,s+1,0] * dn[n,s]) \
                                        + pickup[n,s,i] * ds[n,s] * dn [n,s] / dt ) \
                                        / ( + (-ufn[n,s,0] * dn[n,s]) \
                                        + (-ufs[n,s,0] * dn[n,s]) \
                                        + (ds[n,s] * dn [n,s] / Ts) )   
                    q[n,s]=3  
        # The fourth guadrant including all remainnig unadressed cells  
        for n in range(Ct.shape[0]-2,-1,-1):
            for s in range(1,Ct.shape[1]-1): 
                if (not q[n,s]):
                    if (ufn[n,s,0]<=0) and (ufs[n,s,0]>=0) and (ufn[n+1,s,0]<=0) and (ufs[n,s+1,0]>=0): 
                        # this is the fourth quadrant
                        Ct[n,s,i] = (+ (Ct[n,s-1,i] * ufs[n,s,0] * dn[n,s]) \
                                        + ( -Ct[n+1,s,i] * ufn[n+1,s,0] * dn[n,s]) \
                                        + Cu[n,s,i] * ds[n,s] * dn [n,s] / Ts  ) \
                                        / ( + (ufs[n,s+1,0] * dn[n,s]) \
                                        + (-ufn[n,s,0] * dn[n,s]) \
                                        + (ds[n,s] * dn [n,s] / Ts) )
                        #calculate pickup                
                        pickup[n,s,i] = (Cu[n,s,i]-Ct[n,s,i]) * dt/Ts
                        #check for supply limitations and re-iterate concentration to account for supply limitations
                        if pickup[n,s,i]>mass[n,s,0,i]:
                            pickup[n,s,i] = mass[n,s,0,i]              
                            Ct[n,s,i] = (+ (Ct[n,s-1,i] * ufs[n,s,0] * dn[n,s]) \
                                            + ( -Ct[n+1,s,i] * ufn[n+1,s,0] * dn[n,s]) \
                                            + pickup[n,s,i] * ds[n,s] * dn [n,s] / dt ) \
                                            / ( + (ufs[n,s+1,0] * dn[n,s]) \
                                            + (-ufn[n,s,0] * dn[n,s]) \
                                            + (ds[n,s] * dn [n,s] / Ts) )        
                        q[n,s]=4
                    else:
                        # This is where we apply a generic stencil where all posible directions on the grid boundaries are solved for.
                        # all remaining cells will be calculated for.
                        Ct[n,s,i] = (+ (ufn[n,s,0]>0) * (Ct[n-1,s,i] * ufn[n,s,0] * ds[n,s]) \
                                        + (ufs[n,s,0]>0) * (Ct[n,s-1,i] * ufs[n,s,0] * dn[n,s]) \
                                        + (ufn[n+1,s,0]<0) * ( -Ct[n+1,s,i] * ufn[n+1,s,0] * dn[n,s]) \
                                        + (ufs[n,s+1,0]<0) * ( -Ct[n,s+1,i] * ufs[n,s+1,0] * dn[n,s]) \
                                        + Cu[n,s,i] * ds[n,s] * dn [n,s] / Ts  ) \
                                        / ( + (ufn[n+1,s,0]>0) * (ufn[n+1,s,0] * ds[n,s]) \
                                        + (ufs[n,s+1,0]>0) * (ufs[n,s+1,0] * dn[n,s]) \
                                        + (ufn[n,s,0]<0) * (-ufn[n,s,0] * dn[n,s]) \
                                        + (ufs[n,s,0]<0) * (-ufs[n,s,0] * dn[n,s]) \
                                        + (ds[n,s] * dn [n,s] / Ts) )
                        #calculate pickup                
                        pickup[n,s,i] = (Cu[n,s,i]-Ct[n,s,i]) * dt/Ts
                        #check for supply limitations and re-iterate concentration to account for supply limitations
                        if pickup[n,s,i]>mass[n,s,0,i]:
                            pickup[n,s,i] = mass[n,s,0,i]              
                            Ct[n,s,i] = (+ (ufn[n,s,0]>0) * (Ct[n-1,s,i] * ufn[n,s,0] * ds[n,s]) \
                                            + (ufs[n,s,0]>0) * (Ct[n,s-1,i] * ufs[n,s,0] * dn[n,s]) \
                                            + (ufn[n+1,s,0]<0) * ( -Ct[n+1,s,i] * ufn[n+1,s,0] * dn[n,s]) \
                                            + (ufs[n,s+1,0]<0) * ( -Ct[n,s+1,i] * ufs[n,s+1,0] * dn[n,s]) \
                                            + pickup[n,s,i] * ds[n,s] * dn [n,s] / dt ) \
                                            / ( + (ufn[n+1,s,0]>0) * (ufn[n+1,s,0] * ds[n,s]) \
                                            + (ufs[n,s+1,0]>0) * (ufs[n,s+1,0] * dn[n,s]) \
                                            + (ufn[n,s,0]<0) * (-ufn[n,s,0] * dn[n,s]) \
                                            + (ufs[n,s,0]<0) * (-ufs[n,s,0] * dn[n,s]) \
                                            + (ds[n,s] * dn [n,s] / Ts) )   
                        q[n,s]=5


        k+=1
    
    # print("q1 = " + str(np.sum(q==1)) + "     q2 = " + str(np.sum(q==2)) \
    #       + "     q3 = " + str(np.sum(q==3)) + "     q4 = " + str(np.sum(q==4)) \
    #         + "     q5 = " + str(np.sum(q==5)))
    return Ct, pickup

@njit
def sweep3(Ct, Cu, mass, dt, Ts, ds, dn, us, un):

    pickup = np.zeros(Cu.shape)
    i=0
    k=0

    # Are the lateral boundary conditions circular?
    circ_lateral = False
    if Ct[0,1,0]==-1:
        circ_lateral = True
        Ct[0,:,0] = 0                
        Ct[-1,:,0] = 0

    circ_offshore = False
    if Ct[1,0,0]==-1:
        circ_offshore = True
        Ct[:,0,0] = 0                
        Ct[:,-1,0] = 0

    recirc_offshore = False
    if Ct[1,0,0]==-2:
        recirc_offshore = True
        Ct[:,0,0] = 0                
        Ct[:,-1,0] = 0
    
    
    ufs = np.zeros((np.shape(us)[0], np.shape(us)[1]+1, np.shape(us)[2]))
    ufn = np.zeros((np.shape(un)[0]+1, np.shape(un)[1], np.shape(un)[2]))
    
    # define fluxes
    ufs[:,1:-1, :] = 0.5*us[:,:-1, :] + 0.5*us[:,1:, :]
    ufn[1:-1,:, :] = 0.5*un[:-1,:, :] + 0.5*un[1:,:, :]

    # print(ufs[5,:,0])
    
    #boundary values
    ufs[:,0, :]  = us[:,0, :]
    ufs[:,-1, :] = us[:,-1, :]
   
    ufn[0,:, :]  = un[0,:, :]
    ufn[-1,:, :] = un[-1,:, :]

    Ct_last = Ct.copy()
    while k==0 or np.any(np.abs(Ct[:,:,i]-Ct_last[:,:,i])>1e-10):
    # while k==0 or np.any(np.abs(Ct[:,:,i]-Ct_last[:,:,i])!=0):
        Ct_last = Ct.copy()

        # lateral boundaries circular
        if circ_lateral:
            Ct[0,:,0],Ct[-1,:,0] = Ct[-1,:,0].copy(),Ct[0,:,0].copy()
        if circ_offshore:
            Ct[:,0,0],Ct[:,-1,0] = Ct[:,-1,0].copy(),Ct[:,0,0].copy()

        if recirc_offshore:
            # print(Ct[:,1,0])
            # print(Ct[:,-2,0]) 
            Ct[:,0,0],Ct[:,-1,0] = np.average(Ct[:,-2,0]),np.average(Ct[:,1,0])
            # print(Ct[:,0,0])
            # print(Ct[:,-1,0]) 

        # make an array with a bolean operator. This keeps track of considerd cells. We start with all False (not considered)
        q = np.zeros(Cu.shape[:2]) 

    ########################################################################################
        # in this sweeping algorithm we sweep over the 4 quadrants
        # assuming that most cells have no converging/divering charactersitics.
        # In the last quadrant we take converging and diverging cells into account. 

        # The First quadrant  
        for n in range(1,Ct.shape[0]):
            for s in range(1,Ct.shape[1]):
                if (not q[n,s]) and (ufn[n,s,0]>=0) and (ufs[n,s,0]>=0) and (ufn[n+1,s,0]>=0) and (ufs[n,s+1,0]>=0):
                    Ct[n,s,i] = (+ (Ct[n-1,s,i] * ufn[n,s,0] * ds[n,s]) \
                                    + (Ct[n,s-1,i] * ufs[n,s,0] * dn[n,s]) \
                                    + Cu[n,s,i] * ds[n,s] * dn [n,s] / Ts  ) \
                                    / ( + (ufn[n+1,s,0] * ds[n,s]) \
                                    + (ufs[n,s+1,0] * dn[n,s]) \
                                    + (ds[n,s] * dn [n,s] / Ts) )
                    #calculate pickup                
                    pickup[n,s,i] = (Cu[n,s,i]-Ct[n,s,i]) * dt/Ts
                    #check for supply limitations and re-iterate concentration to account for supply limitations
                    if pickup[n,s,i]>mass[n,s,0,i]:
                        pickup[n,s,i] = mass[n,s,0,i]              
                        Ct[n,s,i] = (+ (Ct[n-1,s,i] * ufn[n,s,0] * ds[n,s]) \
                                    + (Ct[n,s-1,i] * ufs[n,s,0] * dn[n,s]) \
                                        + pickup[n,s,i] * ds[n,s] * dn [n,s] / dt ) \
                                        / (+(ufn[n+1,s,0] * ds[n,s]) \
                                        + (ufs[n,s+1,0] * dn[n,s]))
                    q[n,s]=1
        # The second quadrant
        for n in range(1,Ct.shape[0]):
            for s in range(Ct.shape[1]-2,-1,-1):  
                if (not q[n,s]) and (ufn[n,s,0]>=0) and (ufs[n,s,0]<=0) and (ufn[n+1,s,0]>=0) and (ufs[n,s+1,0]<=0):
                    Ct[n,s,i] = (+ (Ct[n-1,s,i] * ufn[n,s,0] * ds[n,s]) \
                                    + ( -Ct[n,s+1,i] * ufs[n,s+1,0] * dn[n,s]) \
                                    + Cu[n,s,i] * ds[n,s] * dn [n,s] / Ts  ) \
                                    / ( + (ufn[n+1,s,0] * ds[n,s]) \
                                    + (-ufs[n,s,0] * dn[n,s]) \
                                    + (ds[n,s] * dn [n,s] / Ts) )
                    #calculate pickup                
                    pickup[n,s,i] = (Cu[n,s,i]-Ct[n,s,i]) * dt/Ts
                    #check for supply limitations and re-iterate concentration to account for supply limitations
                    if pickup[n,s,i]>mass[n,s,0,i]:
                        pickup[n,s,i] = mass[n,s,0,i]              
                        Ct[n,s,i] = (+ (Ct[n-1,s,i] * ufn[n,s,0] * ds[n,s]) \
                                        + ( -Ct[n,s+1,i] * ufs[n,s+1,0] * dn[n,s]) \
                                        + pickup[n,s,i] * ds[n,s] * dn [n,s] / dt ) \
                                        / ( + (ufn[n+1,s,0] * ds[n,s]) \
                                        + (-ufs[n,s,0] * dn[n,s]))  
                    q[n,s]=2
        # The third quadrant
        for n in range(Ct.shape[0]-2,-1,-1):
            for s in range(Ct.shape[1]-2,-1,-1):
                if (not q[n,s]) and (ufn[n,s,0]<=0) and (ufs[n,s,0]<=0) and (ufn[n+1,s,0]<=0) and (ufs[n,s+1,0]<=0):
                    Ct[n,s,i] = (+ ( -Ct[n+1,s,i] * ufn[n+1,s,0] * dn[n,s]) \
                                    + ( -Ct[n,s+1,i] * ufs[n,s+1,0] * dn[n,s]) \
                                    + Cu[n,s,i] * ds[n,s] * dn [n,s] / Ts  ) \
                                    / ( + (-ufn[n,s,0] * dn[n,s]) \
                                    + (-ufs[n,s,0] * dn[n,s]) \
                                    + (ds[n,s] * dn [n,s] / Ts) )
                    #calculate pickup                
                    pickup[n,s,i] = (Cu[n,s,i]-Ct[n,s,i]) * dt/Ts
                    #check for supply limitations and re-iterate concentration to account for supply limitations
                    if pickup[n,s,i]>mass[n,s,0,i]:
                        pickup[n,s,i] = mass[n,s,0,i]              
                        Ct[n,s,i] = (+ ( -Ct[n+1,s,i] * ufn[n+1,s,0] * dn[n,s]) \
                                        + ( -Ct[n,s+1,i] * ufs[n,s+1,0] * dn[n,s]) \
                                        + pickup[n,s,i] * ds[n,s] * dn [n,s] / dt ) \
                                        / ( + (-ufn[n,s,0] * dn[n,s]) \
                                        + (-ufs[n,s,0] * dn[n,s]))   
                    q[n,s]=3  
        # The fourth guadrant including all remainnig unadressed cells  
        for n in range(Ct.shape[0]-2,-1,-1):
            for s in range(1,Ct.shape[1]): 
                if (not q[n,s]):
                    if (ufn[n,s,0]<=0) and (ufs[n,s,0]>=0) and (ufn[n+1,s,0]<=0) and (ufs[n,s+1,0]>=0): 
                        # this is the fourth quadrant
                        Ct[n,s,i] = (+ (Ct[n,s-1,i] * ufs[n,s,0] * dn[n,s]) \
                                        + ( -Ct[n+1,s,i] * ufn[n+1,s,0] * dn[n,s]) \
                                        + Cu[n,s,i] * ds[n,s] * dn [n,s] / Ts  ) \
                                        / ( + (ufs[n,s+1,0] * dn[n,s]) \
                                        + (-ufn[n,s,0] * dn[n,s]) \
                                        + (ds[n,s] * dn [n,s] / Ts) )
                        #calculate pickup                
                        pickup[n,s,i] = (Cu[n,s,i]-Ct[n,s,i]) * dt/Ts
                        #check for supply limitations and re-iterate concentration to account for supply limitations
                        if pickup[n,s,i]>mass[n,s,0,i]:
                            pickup[n,s,i] = mass[n,s,0,i]              
                            Ct[n,s,i] = (+ (Ct[n,s-1,i] * ufs[n,s,0] * dn[n,s]) \
                                            + ( -Ct[n+1,s,i] * ufn[n+1,s,0] * dn[n,s]) \
                                            + pickup[n,s,i] * ds[n,s] * dn [n,s] / dt ) \
                                            / ( + (ufs[n,s+1,0] * dn[n,s]) \
                                            + (-ufn[n,s,0] * dn[n,s]))        
                        q[n,s]=4
                    else:
                        if (not n==0) and (not s==Ct.shape[1]-1):
                            # This is where we apply a generic stencil where all posible directions on the cell boundaries are solved for.
                            # all remaining cells will be calculated for and q=5 is assigned. 
                            # this stencil is nested in the q4 loop which is the final quadrant.
                            # grid boundaries are filtered in both if statements.
                            Ct[n,s,i] = (+ (ufn[n,s,0]>0) * (Ct[n-1,s,i] * ufn[n,s,0] * ds[n,s]) \
                                            + (ufs[n,s,0]>0) * (Ct[n,s-1,i] * ufs[n,s,0] * dn[n,s]) \
                                            + (ufn[n+1,s,0]<0) * ( -Ct[n+1,s,i] * ufn[n+1,s,0] * dn[n,s]) \
                                            + (ufs[n,s+1,0]<0) * ( -Ct[n,s+1,i] * ufs[n,s+1,0] * dn[n,s]) \
                                            + Cu[n,s,i] * ds[n,s] * dn [n,s] / Ts  ) \
                                            / ( + (ufn[n+1,s,0]>0) * (ufn[n+1,s,0] * ds[n,s]) \
                                            + (ufs[n,s+1,0]>0) * (ufs[n,s+1,0] * dn[n,s]) \
                                            + (ufn[n,s,0]<0) * (-ufn[n,s,0] * dn[n,s]) \
                                            + (ufs[n,s,0]<0) * (-ufs[n,s,0] * dn[n,s]) \
                                            + (ds[n,s] * dn [n,s] / Ts) )
                            #calculate pickup                
                            pickup[n,s,i] = (Cu[n,s,i]-Ct[n,s,i]) * dt/Ts
                            #check for supply limitations and re-iterate concentration to account for supply limitations
                            if pickup[n,s,i]>mass[n,s,0,i]:
                                pickup[n,s,i] = mass[n,s,0,i]              
                                Ct[n,s,i] = (+ (ufn[n,s,0]>0) * (Ct[n-1,s,i] * ufn[n,s,0] * ds[n,s]) \
                                                + (ufs[n,s,0]>0) * (Ct[n,s-1,i] * ufs[n,s,0] * dn[n,s]) \
                                                + (ufn[n+1,s,0]<0) * ( -Ct[n+1,s,i] * ufn[n+1,s,0] * dn[n,s]) \
                                                + (ufs[n,s+1,0]<0) * ( -Ct[n,s+1,i] * ufs[n,s+1,0] * dn[n,s]) \
                                                + pickup[n,s,i] * ds[n,s] * dn [n,s] / dt ) \
                                                / ( + (ufn[n+1,s,0]>0) * (ufn[n+1,s,0] * ds[n,s]) \
                                                + (ufs[n,s+1,0]>0) * (ufs[n,s+1,0] * dn[n,s]) \
                                                + (ufn[n,s,0]<0) * (-ufn[n,s,0] * dn[n,s]) \
                                                + (ufs[n,s,0]<0) * (-ufs[n,s,0] * dn[n,s]))   
                            q[n,s]=5


        k+=1
    
    # print("q1 = " + str(np.sum(q==1)) + "     q2 = " + str(np.sum(q==2)) \
    #       + "     q3 = " + str(np.sum(q==3)) + "     q4 = " + str(np.sum(q==4)) \
    #         + "     q5 = " + str(np.sum(q==5)))
    return Ct, pickup
