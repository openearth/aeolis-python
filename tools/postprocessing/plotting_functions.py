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



def plot1d(ncfile, itransects=0, itimes=0, ifrac=0, params='zb'):
    
    '''Generic function for plotting 1D parameters from AeoLiS netcdf results

    Parameters
    ----------
    ncfile : str
        filename of the netcdf-file containing the results
    itransects : integer or array of integers
        index or list of indices describing the transect-indices (y) in the netcdf-file for plotting
    itimes : integer or array of integers
        index or list of indices describing the timestep-indices in the netcdf-file for plotting
    ifrac : integer
        index of fraction, only used when the results has multiple dimensions ( >=4 )
    params :
        string or list of strings describing the name of the parameter for plotting'
        the parameters that can be plotted are defined in the AeoLiS input file.

    Returns
    -------
    figs/ ax
        ...

    Examples
    --------
    >>> plot1d(ncfile, itransects=0, itimes=0, ifrac=0, params='zb')            # plot the bed level for the first transect and timestep 0
    >>> plot1d(ncfile, itransects=0, itimes=-1, ifrac=0, params=['zb', 'Hs'])   # plot the bed level and wave height for the first transect and last timestep

    '''
    
    if type(itimes) == int:
        itimes = [itimes]
        
    if type(itransects) == int:
        itransects = [itransects]
        
    if type(params) == str:
        params = [params]
    
    # Plotting
    fig, ax = plt.subplots()
                    
    with netCDF4.Dataset(ncfile, 'r') as ds:
        
        for it in itimes:
            for iy in itransects:
                for param in params:
 
                    # get spatial dimensions and bed levels
                    x = ds.variables['x'][iy, :]
                    
                    # Taking the values 
                    val_full = ds.variables[param][it, ::]
                        
                    # Selecting fractions if included
                    if len(val_full.shape) == 3:
                        val = val_full[iy, :, ifrac]
                    else: 
                        val = val_full[iy, :]
 
                    ax.plot(x, val)
        
    return fig, ax



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

    Returns
    -------
    figs/ ax
        ...

    Examples
    --------
    >>> plot2d(ncfile, itimes=0, ifrac=0, param='zb')            # plot the bed level for True first timestep
    >>> plot2d(ncfile, itimes=-1, ifrac=0, param='zb', cmap=plt.cm.jet, clim='auto', delta=False, itimeDelta=0)   # plot the bed level difference for the first timestep

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
            ax.set_aspect('equal')
            plt.colorbar(pc)
    
    return fig, ax
