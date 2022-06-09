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
    
    '''Generic function for plotting 2D parameters from AeoLiS netcdf results

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
        string or list of strings describing the name of the parameter for plotting

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
            ax.set_aspect('equal')
            plt.colorbar(pc)
    
    return fig, ax


def plot3d(ncfile, itimes, scalez=1., vangle=90., hangle=0.):
    
    '''Generic function for plotting 2D parameters from AeoLiS netcdf results

    Parameters
    ----------
    ncfile : str
        filename of the netcdf-file containing the results
    itimes : integer or array of integers
        index or list of indices describing the timestep-indices in the netcdf-file for plotting

    Returns
    -------
    figs/ axs?
        ...

    Examples
    --------
    >>> ploted(parse_value('T'))
        bool

    '''
    
    rgb=255.
    sand = [230/rgb,230/rgb,210/rgb]
    blue = [0/rgb, 166/rgb, 214/rgb]
    sand = [230/rgb,230/rgb,210/rgb]
    vegetation = [120/rgb,180/rgb,50/rgb]
    white = [255/rgb,255/rgb,255/rgb]
    gray = [100/rgb,100/rgb,100/rgb]
    
        
    if type(itimes) == int:
        itimes = [itimes]
        
    for it in itimes:
    
        with netCDF4.Dataset(ncfile, 'r') as ds:
            
            # get spatial dimensions and bed levels
            x = ds.variables['x'][:,:]
            y = ds.variables['y'][:,:]
            z = ds.variables['zb'][it,:,:]
            
            fig = plt.figure()
            ax = plt.gca(projection='3d')
    
            x_scale = 1.
            y_scale = 1. # np.max(x)/np.max(y)
            z_scale = scalez*np.max(z)/np.max(x)
            
            scale=np.diag([x_scale, y_scale, z_scale, 1.0])
            scale=scale*(1.0/scale.max())
            scale[3,3]=1.0
            
            def short_proj():
                return np.dot(Axes3D.get_proj(ax), scale)
            
            ax.get_proj = short_proj   
            
            # lay-out of axes
            ax.set_ylabel('alongshore distance [m]')
            ax.set_xlabel('cross-shore distance [m]')
            ax.axes.get_xaxis().set_ticks([])
            ax.axes.get_yaxis().set_ticks([])
            ax.axes.get_xaxis().set_visible(False)
            ax.axes.get_yaxis().set_visible(False)
            
            # perspective
            ax.view_init(vangle, hangle)
            ax.set_proj_type(proj_type='persp')
            
            # make the panes transparent           
            for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
                axis.set_pane_color((0/rgb,166/rgb,214/rgb,0.0))
                axis._axinfo["grid"]['color'] =  (1,1,1,0) 
                axis.line.set_color((1.0, 1.0, 1.0, 0.0)) 
            
            ax.set_axis_off()
            pos_z = z.copy()
            neg_z = z.copy()
            # neg_z[(neg_z > 0.2)] = np.nan
            pos_z[(pos_z < 0.2)] = np.nan
            
            # ax.plot_surface(y, x, neg_z[:,:], color=white, linewidth=0, antialiased=False, shade=False, alpha=1.0, rstride=1, cstride=1)
            ax.plot_surface(y, x, pos_z[:,:], color=sand, linewidth=0, antialiased=False, shade=True, alpha=1.0, rstride=1, cstride=1)          
            
            # plt.savefig(r'c:\Users\weste_bt\AeoLiS\Examples\Parabolic figures\plottoptime' + str('{:03}'.format(time_index)) + '.png', dpi=200)
            # plt.show()
            # plt.close(fig)
            
    return fig, ax


    