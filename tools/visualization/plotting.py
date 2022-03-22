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



def plot2d(ncfile, itimes=0, ifrac=0, param='zb', cmap=plt.cm.jet, clim='auto', delta=False, itimeDelta=0):
    
    '''Generic function for plotting 2D parameters from AeoLiS netcdf results

    Parameters
    ----------
    ncfile : str
        filename of the netcdf-file containing the results
    itimes : integer or array of integers
        index or list of indices describing the timestep-indices in the netcdf-file for plotting
    ifrac : integer
        index of fraction, only used when the results has multiple dimensions ( >=4 )
    param :
        string  describing the name of the parameter for plotting
    cmap : colormap
        matplotlib colormap used for the plot
    clim : list or string
        if 'auto': the min and max values from the chosen parameter are used.
        else the min and max value can be given manually for the colorlimits of the colormap
    delta : Bool (True or False)
        Boolean for indicating plotting difference between timesteps or the actual value
    itimeDelta : integer
        index for the starting timestep in the netcdf-file when merging, only used when delta == Bool
    quiver?

    Returns
    -------
    figs/ axs?
        ...

    Examples
    --------
    >>> ploted(parse_value('T'))
        bool

    '''
    
    if type(itimes) == int:
        itimes = [itimes]
        
    for it in itimes:
    
        with netCDF4.Dataset(ncfile, 'r') as ds:
        
            # get spatial dimensions and bed levels
            x = ds.variables['x'][:,:]
            y = ds.variables['y'][:,:]
            
            # Taking delta or acutal value 
            if delta == True:
                val = ds.variables[param][itimeDelta, ::] - ds.variables[param][it, ::]
            else:
                val = ds.variables[param][it, ::]
                
            # Selecting fractions if included
            if len(val.shape) == 3:
                val = val[:, :, ifrac]
            
            # Determining clim
            if clim == 'auto':
                vmin = np.min(val)
                vmax = np.max(val)
            else:
                vmin = clim[0]
                vmax = clim[1]
                
            # Plotting
            fig, ax = plt.subplots()
            pc = ax.pcolormesh(x, y, val, cmap=cmap, vmin=vmin, vmax=vmax)
            plt.colorbar(pc)
    
    return fig, ax


    