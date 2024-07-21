import os
import glob
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib.animation import FuncAnimation


# import colormap
colors = np.loadtxt('cmap_sandmotor.txt')
cmap_sandmotor = matplotlib.colors.ListedColormap(colors)


def plot_bathymetry(ncfile, change=False, figsize=(10,5), time_index=-1, ax=None):
    '''Plot bathymetries or the erosion/deposition from an AeoLiS result file
    
    Parameters
    ----------
    ncfile : str
      Path to AeoLiS result file
    change : bool
      Plot bathymetric change rather than actual bathymetry
    figsize : 2-tuple, optional
      Dimensions of resulting figure
    time_index : int, optional
      Index of time dimension to plot [default: -1]
    ax : matplotlib.axes.SubplotAxis, optional
      Axis used for plotting, if not given a new figure is created
      
    Returns
    -------
    ax : matplotlib.axes.SubplotAxis
      plot axes objects
      
    '''
    
    with netCDF4.Dataset(ncfile, 'r') as ds:
        
        # get spatial dimensions and bed levels
        x = ds.variables['x'][:,:]
        y = ds.variables['y'][:,:]
        zb = ds.variables['zb'][...]
        pickup = ds.variables['pickup_sum'][...]
        
        # create figure
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        
        # plot bed levels and bed level change
        if not change:
            p = ax.pcolormesh(y, x, zb[time_index,:,:], cmap=cmap_sandmotor, vmin=-7.5, vmax=7.5)
            cb = fig.colorbar(p, shrink=.7)
            cb.set_label('bed level [m]')
        else:
            pickup = pickup.sum(axis=-1)             # sum over fractions
            pickup = np.cumsum(pickup, axis=0)       # cummulative sum in time
            dz = -pickup / (2650. * .6)               # convert from kg/m2 to m3/m2
            p = ax.pcolormesh(y, x, dz[time_index,:,:], cmap='bwr_r', vmin=-2, vmax=2)
            cb = fig.colorbar(p, shrink=.7)
            cb.set_label('bed level change [m]')

        ax.contour(y, x, zb[0,:,:], [0.], colors=['k'])
            
        ax.set_aspect('equal', adjustable='box')
        ax.invert_yaxis()
            
        ax.set_xlabel('alongshore distance [m]')
        ax.set_ylabel('cross-shore distance [m]')

    return ax


def plot_erosion(ncfile, figsize=(10,4), ax=None):
    '''Plot erosion in time from an AeoLiS result file
    
    Parameters
    ----------
    ncfile : str
      Path to AeoLiS result file
    figsize : 2-tuple, optional
      Dimensions of resulting figure
    ax : matplotlib.axes.SubplotAxis, optional
      Axis used for plotting, if not given a new figure is created
      
    Returns
    -------
    ax : matplotlib.axes.SubplotAxis
      plot axes object
      
    '''

    name = os.path.split(ncfile)[1]
    
    with netCDF4.Dataset(ncfile, 'r') as ds:

        # get spatial dimensions, time and pickup
        x = ds.variables['x'][:,:]
        y = ds.variables['y'][:,:]
        t = ds.variables['time'][:]
        pickup = ds.variables['pickup_sum'][...]
        
        dx = np.diff(x[0,:])[0] # assume equidistant grid
        dy = np.diff(y[:,0])[0] # assume equidistant grid
        
        # create figure
        if ax is None:
            fig, ax = plt.subplots(1, 1, figsize=figsize)
        
        # convert pickup to volumes
        pickup = pickup.sum(axis=-1)             # sum over fractions
        pickup = np.cumsum(pickup, axis=0)       # cummulative sum in time
        pickup = pickup * dx * dy / (2650. * .6) # convert from kg/m2 to m3
        
        # determine erosion and deposition and integrat in space
        E = np.maximum(0., pickup).sum(axis=-1).sum(axis=-1) * 1e-6
        D = np.minimum(0., pickup).sum(axis=-1).sum(axis=-1) * 1e-6
        
        # plot erosion in time
        dt = netCDF4.num2date(t, ds.variables['time'].units,only_use_cftime_datetimes=False)
        ax.plot(dt, E, label=name)
        
        # set titles and labels
        ax.set_xlabel('time')
        ax.set_ylabel('erosion volume [$\mathdefault{Mm^3}$]')
        ax.set_title(name)
        
        ax.grid(True)
            
    return ax


def plot_erosion_multi(ncfiles, ax=None, **kwargs):
    '''Plot erosion in time from multiple AeoLiS result files
    
    Parameters
    ----------
    ncfiles : str
      Path to AeoLiS result files, may contain wildcards (*)
    
    See ``plot_erosion`` for options and return values.
    
    See Also
    --------
    plot_erosion
    
    '''
    print(glob.glob(ncfiles))
    for ncfile in glob.glob(ncfiles):      
        ax = plot_erosion(ncfile, ax=ax, **kwargs)
    ax.legend(loc='upper left')
    ax.set_title('')
            
    return ax


def plot_coverage(ncfile, figsize=(10,4), ax=None):
    '''Plot coverage of non-erodible elements in top layer of bed from an AeoLiS result file
    
    Parameters
    ----------
    ncfile : str
      Path to AeoLiS result file
    figsize : 2-tuple, optional
      Dimensions of resulting figure
    ax : matplotlib.axes.SubplotAxis, optional
      Axis used for plotting, if not given a new figure is created
      
    Returns
    -------
    ax : matplotlib.axes.SubplotAxis
      plot axes object
      
    '''

    name = os.path.split(ncfile)[1]

    with netCDF4.Dataset(ncfile, 'r') as ds:
        
        # get spatial dimensions and bed mass
        x = ds.variables['x'][:,:]
        y = ds.variables['y'][:,:]
        zb = ds.variables['zb'][0,:,:]
        mass = ds.variables['mass'][-1,:,:,:,:]

        # compute coverage of largest fraction
        if mass.ndim < 4 or mass.shape[-1] < 2:
            coverage = np.ones(x.shape)
        else:
            coverage = mass[:,:,0,-1] / mass[:,:,0,:].sum(axis=-1)
        
        # create figure
        if ax is None:
            fig, ax = plt.subplots(1, 1, figsize=figsize)
        
        # plot bed levels and bed level change
        p = ax.pcolormesh(y, x, coverage, cmap='Reds', vmin=0, vmax=1)
        ax.contour(y, x, zb, [0.], colors=['k'])
        
        # create colorbars
        cb = fig.colorbar(p, shrink=.7)
        cb.set_label('coverage [%]')
    
        ax.set_aspect('equal', adjustable='box')
        ax.invert_yaxis()
        ax.set_xlabel('alongshore distance [m]')
        ax.set_ylabel('cross-shore distance [m]')
        ax.set_title(name)

    return ax


def create_animation(ncfile, figsize=(10,4), ext='mp4', nframes=None):
    '''Create animation from an AeoLiS result file

    Parameters
    ----------
    ncfile : str
      Path to AeoLiS result file
    figsize : 2-tuple, optional
      Dimensions of resulting video
    ext : str, optional
      Extension of resulting video file [default: mp4]
    nframes : int, optional
      Number of frames to include in video
      
    Returns
    -------
    videofile : str
      Path to video file
      
    '''


    fig, ax = plt.subplots(figsize=figsize)
    ax.invert_yaxis()
    
    with netCDF4.Dataset(ncfile, 'r') as ds:
        t = ds.variables['time'][:]
        x = ds.variables['x'][:,:]
        y = ds.variables['y'][:,:]
        zb = ds.variables['zb'][:,:,:]
        units = ds.variables['time'].units
        
    def update(i):
        ax.clear()
        p = ax.pcolormesh(y, x, zb[i,:,:], cmap=cmap_sandmotor, vmin=-7.5, vmax=7.5)
        ax.set_aspect('equal', adjustable='box')
        ax.invert_yaxis()
        ax.set_xlabel('alongshore [m]')
        ax.set_ylabel('cross-shore [m]')
        ax.set_title(netCDF4.num2date(t[i], units))
        return p

    
    if nframes is None:
        nframes = zb.shape[0]
    nframes = int(np.minimum(zb.shape[0], nframes))
    
    p = update(0)
    cb = fig.colorbar(p, shrink=.7)
    cb.set_label('bed level [m]')

    videofile = '%s.%s' % (os.path.splitext(ncfile)[0], ext)
    anim = FuncAnimation(fig, update, frames=nframes, interval=30)
    anim.save(videofile, dpi=150, writer='ffmpeg')

    return videofile
