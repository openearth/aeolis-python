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

# Note: @njit(cache=True) is intentionally not used here.
# This function acts as an orchestrator, delegating work to Numba-compiled helper functions.
# Decorating the orchestrator itself with njit provides no performance benefit,
# since most of the computation is already handled by optimized Numba functions.
def sweep(Ct, Cu, mass, dt, Ts, ds, dn, us, un, w):


    pickup = np.zeros(Cu.shape)
    i=0
    k=0

    nf = np.shape(Ct)[2]

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
    
    # define velocity at cell faces
    ufs[:,1:-1, :] = 0.5*us[:,:-1, :] + 0.5*us[:,1:, :]
    ufn[1:-1,:, :] = 0.5*un[:-1,:, :] + 0.5*un[1:,:, :]

    # print(ufs[5,:,0])

    # set empty boundary values, extending the velocities at the boundaries
    ufs[:,0, :]  = ufs[:,1, :]
    ufs[:,-1, :] = ufs[:,-2, :]
   
    ufn[0,:, :]  = ufn[1,:, :]
    ufn[-1,:, :] = ufn[-2,:, :]
    
    # Lets take the average of the top and bottom and left/right boundary cells
    # apply the average to the boundary cells
    # this ensures that the inflow at one side is equal to the outflow at the other side

    ufs[:,0,:]  = (ufs[:,0,:]+ufs[:,-1,:])/2
    ufs[:,-1,:] = ufs[:,0,:]     
    ufs[0,:,:]  = (ufs[0,:,:]+ufs[-1,:,:])/2
    ufs[-1,:,:] = ufs[0,:,:]     
    
    ufn[:,0,:]  = (ufn[:,0,:]+ufn[:,-1,:])/2
    ufn[:,-1,:] = ufn[:,0,:]     
    ufn[0,:,:]  = (ufn[0,:,:]+ufn[-1,:,:])/2
    ufn[-1,:,:] = ufn[0,:,:] 

    # now make sure that there is no gradients at the boundaries
    # ufs[:,1,:]  = ufs[:,0,:]
    # ufs[:,-2,:] = ufs[:,-1,:]
    # ufs[1,:,:]  = ufs[0,:,:]
    # ufs[-2,:,:] = ufs[-1,:,:]

    # ufn[:,1,:]  = ufn[:,0,:]
    # ufn[:,-2,:] = ufn[:,-1,:]
    # ufn[1,:,:]  = ufn[0,:,:]
    # ufn[-2,:,:] = ufn[-1,:,:]

    # ufn[:,:,:] = ufn[-2,:,:]

    # also correct for the potential gradients at the boundary cells in the equilibrium concentrations
    Cu[:,0,:]  = Cu[:,1,:]
    Cu[:,-1,:] = Cu[:,-2,:]
    Cu[0,:,:]  = Cu[1,:,:]
    Cu[-1,:,:] = Cu[-2,:,:]
    
    # #boundary values
    # ufs[:,0, :]  = us[:,0, :]
    # ufs[:,-1, :] = us[:,-1, :]
   
    # ufn[0,:, :]  = un[0,:, :]
    # ufn[-1,:, :] = un[-1,:, :]

    Ct_last = Ct.copy()
    while k==0 or np.any(np.abs(Ct[:,:,i]-Ct_last[:,:,i])>1e-10):
    # while k==0 or np.any(np.abs(Ct[:,:,i]-Ct_last[:,:,i])!=0):
        Ct_last = Ct.copy()

        # lateral boundaries circular
        if circ_lateral:
            Ct[0,:,0],Ct[-1,:,0] = Ct[-1,:,0].copy(),Ct[0,:,0].copy()
            # pickup[0,:,0],pickup[-1,:,0] = pickup[-1,:,0].copy(),pickup[0,:,0].copy()
        if circ_offshore:
            Ct[:,0,0],Ct[:,-1,0] = Ct[:,-1,0].copy(),Ct[:,0,0].copy()
            # pickup[:,0,0],pickup[:,-1,0] = pickup[:,-1,0].copy(),pickup[:,0,0].copy()

        if recirc_offshore:
            Ct[:,0,0],Ct[:,-1,0] = np.mean(Ct[:,-2,0]), np.mean(Ct[:,1,0])

        # Track visited  cells and quadrant classification
        visited = np.zeros(Cu.shape[:2], dtype=bool)
        quad = np.zeros(Cu.shape[:2], dtype=np.uint8)

    ########################################################################################
        # in this sweeping algorithm we sweep over the 4 quadrants
        # assuming that most cells have no converging/divering charactersitics.
        # In the last quadrant we take converging and diverging cells into account. 

        # The First quadrant (Numba-optimized)
        _solve_quadrant1(Ct, Cu, mass, pickup, dt, Ts, ds, dn, ufs, ufn, w, visited, quad, nf)
        
        # The second quadrant (Numba-optimized)
        _solve_quadrant2(Ct, Cu, mass, pickup, dt, Ts, ds, dn, ufs, ufn, w, visited, quad, nf)
        
        # The third quadrant (Numba-optimized)
        _solve_quadrant3(Ct, Cu, mass, pickup, dt, Ts, ds, dn, ufs, ufn, w, visited, quad, nf)
        
        # The fourth quadrant (Numba-optimized)
        _solve_quadrant4(Ct, Cu, mass, pickup, dt, Ts, ds, dn, ufs, ufn, w, visited, quad, nf)
        
        # Generic stencil for remaining cells including boundaries (Numba-optimized)
        _solve_generic_stencil(Ct, Cu, mass, pickup, dt, Ts, ds, dn, ufs, ufn, w, visited, quad, nf)

        # check the boundaries of the pickup matrix for unvisited cells
        # print(np.shape(visited[0,:]==False))
        pickup[0,:,0] = pickup[1,:,0].copy() 
        pickup[-1,:,0] = pickup[-2,:,0].copy() 
                
        k+=1
    
    # # plot Ct
    # import matplotlib.pyplot as plt
    # plt.imshow(quad[:10,:10], origin='lower')
    # # plt.colorbar()
    # plt.title('Concentration after %d sweeps' % k)
    # plt.show()
    # plt.imshow(Ct[:50,:50], origin='lower')
    # # plt.colorbar()
    # plt.title('Concentration after %d sweeps' % k)
    # plt.show()
    # plt.plot(pickup[0,:,0])
    # plt.plot(pickup[-1,:,0])
    # plt.show()

        # print(k)


    # print("q1 = " + str(np.sum(q==1)) + "     q2 = " + str(np.sum(q==2)) \
    #       + "     q3 = " + str(np.sum(q==3)) + "     q4 = " + str(np.sum(q==4)) \
    #         + "     q5 = " + str(np.sum(q==5)))
    # print("pickup deviation percentage = " + str(pickup.sum()/pickup[pickup>0].sum()*100) + " %")
    # print("pickup deviation percentage = " + str(pickup[1,:,0].sum()/pickup[1,pickup[1,:,0]>0,0].sum()*100) + " %")
    # print("pickup maximum = " + str(pickup.max()) + " mass max = " + str(mass.max()))
    # print("pickup minimum = " + str(pickup.min()))
    # print("pickup average = " + str(pickup.mean()))
    # print("number of cells for pickup maximum = " + str((pickup == mass.max()).sum()))
                                                #  pickup[1,:,0].sum()/pickup[1,pickup[1,:,0]<0,0].sum()

    return Ct, pickup


@njit(cache=True)
def _solve_quadrant1(Ct, Cu, mass, pickup, dt, Ts, ds, dn, ufs, ufn, w, visited, quad, nf):
    """Solve first quadrant (positive flow in both directions) with Numba optimization."""
    for n in range(1, Ct.shape[0]):
        for s in range(1, Ct.shape[1]):
            if (
                (not visited[n, s])
                and (ufn[n, s, 0] >= 0)
                and (ufs[n, s, 0] >= 0)
                and (ufn[n + 1, s, 0] >= 0)
                and (ufs[n, s + 1, 0] >= 0)
            ):
                
                # Compute concentration for all fractions
                for f in range(nf):
                    num = (Ct[n - 1, s, f] * ufn[n, s, f] * ds[n, s] + 
                           Ct[n, s - 1, f] * ufs[n, s, f] * dn[n, s] + 
                           w[n, s, f] * Cu[n, s, f] * ds[n, s] * dn[n, s] / Ts)
                    
                    den = (ufn[n + 1, s, f] * ds[n, s] + 
                           ufs[n, s + 1, f] * dn[n, s] + 
                           ds[n, s] * dn[n, s] / Ts)
                    
                    Ct[n, s, f] = num / den
                    
                    # Calculate pickup
                    pickup[n, s, f] = (w[n, s, f] * Cu[n, s, f] - Ct[n, s, f]) * dt / Ts
                    
                    # Check for supply limitations and re-iterate
                    if pickup[n, s, f] > mass[n, s, 0, f]:
                        pickup[n, s, f] = mass[n, s, 0, f]
                        
                        num_limited = (Ct[n - 1, s, f] * ufn[n, s, f] * ds[n, s] + 
                                      Ct[n, s - 1, f] * ufs[n, s, f] * dn[n, s] + 
                                      pickup[n, s, f] * ds[n, s] * dn[n, s] / dt)
                        
                        den_limited = (ufn[n + 1, s, f] * ds[n, s] + 
                                      ufs[n, s + 1, f] * dn[n, s])
                        
                        Ct[n, s, f] = num_limited / den_limited
                
                visited[n, s] = True
                quad[n, s] = 1


@njit(cache=True)
def _solve_quadrant2(Ct, Cu, mass, pickup, dt, Ts, ds, dn, ufs, ufn, w, visited, quad, nf):
    """Solve second quadrant (positive n-flow, negative s-flow) with Numba optimization."""
    for n in range(1, Ct.shape[0]):
        for s in range(Ct.shape[1] - 2, -1, -1):
            if (
                (not visited[n, s])
                and (ufn[n, s, 0] >= 0)
                and (ufs[n, s, 0] <= 0)
                and (ufn[n + 1, s, 0] >= 0)
                and (ufs[n, s + 1, 0] <= 0)
            ):
                
                # Compute concentration for all fractions
                for f in range(nf):
                    num = (Ct[n - 1, s, f] * ufn[n, s, f] * ds[n, s] + 
                           -Ct[n, s + 1, f] * ufs[n, s + 1, f] * dn[n, s] + 
                           w[n, s, f] * Cu[n, s, f] * ds[n, s] * dn[n, s] / Ts)
                    
                    den = (ufn[n + 1, s, f] * ds[n, s] + 
                           -ufs[n, s, f] * dn[n, s] + 
                           ds[n, s] * dn[n, s] / Ts)
                    
                    Ct[n, s, f] = num / den
                    
                    # Calculate pickup
                    pickup[n, s, f] = (w[n, s, f] * Cu[n, s, f] - Ct[n, s, f]) * dt / Ts
                    
                    # Check for supply limitations and re-iterate
                    if pickup[n, s, f] > mass[n, s, 0, f]:
                        pickup[n, s, f] = mass[n, s, 0, f]
                        
                        num_limited = (Ct[n - 1, s, f] * ufn[n, s, f] * ds[n, s] + 
                                      -Ct[n, s + 1, f] * ufs[n, s + 1, f] * dn[n, s] + 
                                      pickup[n, s, f] * ds[n, s] * dn[n, s] / dt)
                        
                        den_limited = (ufn[n + 1, s, f] * ds[n, s] + 
                                      -ufs[n, s, f] * dn[n, s])
                        
                        Ct[n, s, f] = num_limited / den_limited
                
                visited[n, s] = True
                quad[n, s] = 2


@njit(cache=True)
def _solve_quadrant3(Ct, Cu, mass, pickup, dt, Ts, ds, dn, ufs, ufn, w, visited, quad, nf):
    """Solve third quadrant (negative flow in both directions) with Numba optimization."""
    for n in range(Ct.shape[0] - 2, -1, -1):
        for s in range(Ct.shape[1] - 2, -1, -1):
            if (
                (not visited[n, s])
                and (ufn[n, s, 0] <= 0)
                and (ufs[n, s, 0] <= 0)
                and (ufn[n + 1, s, 0] <= 0)
                and (ufs[n, s + 1, 0] <= 0)
            ):
                
                # Compute concentration for all fractions
                for f in range(nf):
                    num = (-Ct[n + 1, s, f] * ufn[n + 1, s, f] * dn[n, s] + 
                           -Ct[n, s + 1, f] * ufs[n, s + 1, f] * dn[n, s] + 
                           w[n, s, f] * Cu[n, s, f] * ds[n, s] * dn[n, s] / Ts)
                    
                    den = (-ufn[n, s, f] * dn[n, s] + 
                           -ufs[n, s, f] * dn[n, s] + 
                           ds[n, s] * dn[n, s] / Ts)
                    
                    Ct[n, s, f] = num / den
                    
                    # Calculate pickup
                    pickup[n, s, f] = (w[n, s, f] * Cu[n, s, f] - Ct[n, s, f]) * dt / Ts
                    
                    # Check for supply limitations and re-iterate
                    if pickup[n, s, f] > mass[n, s, 0, f]:
                        pickup[n, s, f] = mass[n, s, 0, f]
                        
                        num_limited = (-Ct[n + 1, s, f] * ufn[n + 1, s, f] * dn[n, s] + 
                                      -Ct[n, s + 1, f] * ufs[n, s + 1, f] * dn[n, s] + 
                                      pickup[n, s, f] * ds[n, s] * dn[n, s] / dt)
                        
                        den_limited = (-ufn[n, s, f] * dn[n, s] + 
                                      -ufs[n, s, f] * dn[n, s])
                        
                        Ct[n, s, f] = num_limited / den_limited
                
                visited[n, s] = True
                quad[n, s] = 3


@njit(cache=True)
def _solve_quadrant4(Ct, Cu, mass, pickup, dt, Ts, ds, dn, ufs, ufn, w, visited, quad, nf):
    """Solve fourth quadrant (negative n-flow, positive s-flow) with Numba optimization."""
    for n in range(Ct.shape[0] - 2, -1, -1):
        for s in range(1, Ct.shape[1]):
            if (
                (not visited[n, s])
                and (ufn[n, s, 0] <= 0)
                and (ufs[n, s, 0] >= 0)
                and (ufn[n + 1, s, 0] <= 0)
                and (ufs[n, s + 1, 0] >= 0)
            ):
                
                # Compute concentration for all fractions
                for f in range(nf):
                    num = (Ct[n, s - 1, f] * ufs[n, s, f] * dn[n, s] + 
                           -Ct[n + 1, s, f] * ufn[n + 1, s, f] * dn[n, s] + 
                           w[n, s, f] * Cu[n, s, f] * ds[n, s] * dn[n, s] / Ts)
                    
                    den = (ufs[n, s + 1, f] * dn[n, s] + 
                           -ufn[n, s, f] * dn[n, s] + 
                           ds[n, s] * dn[n, s] / Ts)
                    
                    Ct[n, s, f] = num / den
                    
                    # Calculate pickup
                    pickup[n, s, f] = (w[n, s, f] * Cu[n, s, f] - Ct[n, s, f]) * dt / Ts
                    
                    # Check for supply limitations and re-iterate
                    if pickup[n, s, f] > mass[n, s, 0, f]:
                        pickup[n, s, f] = mass[n, s, 0, f]
                        
                        num_limited = (Ct[n, s - 1, f] * ufs[n, s, f] * dn[n, s] + 
                                      -Ct[n + 1, s, f] * ufn[n + 1, s, f] * dn[n, s] + 
                                      pickup[n, s, f] * ds[n, s] * dn[n, s] / dt)
                        
                        den_limited = (ufs[n, s + 1, f] * dn[n, s] + 
                                      -ufn[n, s, f] * dn[n, s])
                        
                        Ct[n, s, f] = num_limited / den_limited
                
                visited[n, s] = True
                quad[n, s] = 4


@njit(cache=True)
def _solve_generic_stencil(Ct, Cu, mass, pickup, dt, Ts, ds, dn, ufs, ufn, w, visited, quad, nf):
    """Solve remaining cells with generic stencil using conditionals (Numba-optimized)."""
    for n in range(Ct.shape[0] - 2, -1, -1):
        for s in range(1, Ct.shape[1]):
            if (not visited[n, s]) and (n != 0) and (s != Ct.shape[1] - 1):
                # Apply generic stencil with conditionals instead of boolean multiplication
                for f in range(nf):
                    # Initialize with source term
                    num = w[n, s, f] * Cu[n, s, f] * ds[n, s] * dn[n, s] / Ts
                    den = ds[n, s] * dn[n, s] / Ts
                    
                    # Add flux contributions conditionally
                    if ufn[n, s, 0] > 0:
                        num += Ct[n - 1, s, f] * ufn[n, s, f] * ds[n, s]
                    
                    if ufs[n, s, 0] > 0:
                        num += Ct[n, s - 1, f] * ufs[n, s, f] * dn[n, s]
                    
                    if ufn[n + 1, s, 0] < 0:
                        num += -Ct[n + 1, s, f] * ufn[n + 1, s, f] * dn[n, s]
                    elif ufn[n + 1, s, 0] > 0:
                        den += ufn[n + 1, s, f] * ds[n, s]
                    
                    if ufs[n, s + 1, 0] < 0:
                        num += -Ct[n, s + 1, f] * ufs[n, s + 1, f] * dn[n, s]
                    elif ufs[n, s + 1, 0] > 0:
                        den += ufs[n, s + 1, f] * dn[n, s]
                    
                    if ufn[n, s, 0] < 0:
                        den += -ufn[n, s, f] * dn[n, s]
                    
                    if ufs[n, s, 0] < 0:
                        den += -ufs[n, s, f] * dn[n, s]
                    
                    Ct[n, s, f] = num / den
                    
                    # Calculate pickup
                    pickup[n, s, f] = (w[n, s, f] * Cu[n, s, f] - Ct[n, s, f]) * dt / Ts
                    
                    # Check for supply limitations and re-iterate
                    if pickup[n, s, f] > mass[n, s, 0, f]:
                        pickup[n, s, f] = mass[n, s, 0, f]
                        
                        # Recompute with limited pickup
                        num_lim = pickup[n, s, f] * ds[n, s] * dn[n, s] / dt
                        den_lim = 0.0
                        
                        if ufn[n, s, 0] > 0:
                            num_lim += Ct[n - 1, s, f] * ufn[n, s, f] * ds[n, s]
                        
                        if ufs[n, s, 0] > 0:
                            num_lim += Ct[n, s - 1, f] * ufs[n, s, f] * dn[n, s]
                        
                        if ufn[n + 1, s, 0] < 0:
                            num_lim += -Ct[n + 1, s, f] * ufn[n + 1, s, f] * dn[n, s]
                        elif ufn[n + 1, s, 0] > 0:
                            den_lim += ufn[n + 1, s, f] * ds[n, s]
                        
                        if ufs[n, s + 1, 0] < 0:
                            num_lim += -Ct[n, s + 1, f] * ufs[n, s + 1, f] * dn[n, s]
                        elif ufs[n, s + 1, 0] > 0:
                            den_lim += ufs[n, s + 1, f] * dn[n, s]
                        
                        if ufn[n, s, 0] < 0:
                            den_lim += -ufn[n, s, f] * dn[n, s]
                        
                        if ufs[n, s, 0] < 0:
                            den_lim += -ufs[n, s, f] * dn[n, s]
                        
                        Ct[n, s, f] = num_lim / den_lim
                
                visited[n, s] = True
                quad[n, s] = 5

