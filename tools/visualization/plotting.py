import os
import glob
import netCDF4
import numpy as np
from matplotlib import cm
from scipy import special
import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)

from scipy.ndimage import gaussian_filter, uniform_filter



def plot2d(ncfile, itimes, param='zb', cmap=plt.cm.jet, clim='auto'):
    
    '''Generic function for plotting 2D parameters from AeoLiS netcdf results

    Parameters
    ----------
    ncfile : str
        filename of the netcdf-file containing the results
    itimes : integer or array of integers
        integer or list of integers describing the timestep-indices in the netcdf-file for plotting
    param:
        string  describing the name of the parameter for plotting
    cmap: colormap
        matplotlib colormap used for the plot
    clim: list or string
        if 'auto': the min and max values from the chosen parameter are used.
        else the min and max value can be given manually for the colorlimits of the colormap
    quiver

    Returns
    -------
    figs/ axs?
        ...

    Examples
    --------
    >>> ploted(parse_value('T'))
        bool

    '''
    
    if len(itimes) == 1:
        itimes = [itimes]
        
    for it in itimes:
    
        with netCDF4.Dataset(ncfile, 'r') as ds:
        
            # get spatial dimensions and bed levels
            x = ds.variables['x'][:,:]
            y = ds.variables['y'][:,:]
            val = ds.variables[param][it, :, :]
            v = -ds.variables[parameter_quivy][time_index, ::d, ::d]
    
    return


    