import matplotlib
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
import os
import glob
import netCDF4
import numpy as np


cmap = plt.cm.jet
rgb=255.
#colormap = [[230/rgb, 230/rgb, 210/rgb],
#            [120/rgb,   180/rgb, 50/rgb]]

colormapveg = {'red': ((0.,  230/rgb, 230/rgb),
                   (1.,  120/rgb, 120/rgb)),

         'green': ((0.,  230/rgb, 230/rgb),
                   (1., 180/rgb, 180/rgb)),

         'blue':  ((0.,  210/rgb, 210/rgb),
                   (1.,  50/rgb, 50/rgb))}
         
colormap2 = {'red': ((0.,  112/rgb, 122/rgb),
                         (0.5,  230/rgb, 230/rgb),
                       (1.,  183/rgb, 183/rgb)),
    
             'green': ((0.,  146/rgb, 146/rgb),
                         (0.5,  230/rgb, 230/rgb),
                       (1.,  223/rgb, 223/rgb)),
    
             'blue':  ((0.,  190/rgb, 190/rgb),
                         (0.5,  210/rgb, 210/rgb),
                       (1.,  140/rgb, 140/rgb))}

plt.register_cmap(name='Vegetation', data = colormap2)

colormapwater = {'red': ((0.,  230/rgb, 230/rgb),
                   (1.,  0/rgb, 0/rgb)),

         'green': ((0.,  230/rgb, 230/rgb),
                   (1., 105/rgb, 105/rgb)),

         'blue':  ((0.,  210/rgb, 210/rgb),
                   (1.,  148/rgb, 148/rgb))}

plt.register_cmap(name='Water', data = colormapwater)


ncfile = 'C:/Users/Aeolis/dune development/parabolic/001/aeolis.nc'

nx = 300
ny = 200

blue = [0/rgb, 166/rgb, 214/rgb]
sand = [230/rgb,230/rgb,210/rgb]
vegetation = [120/rgb,180/rgb,50/rgb]
white = [255/rgb,255/rgb,255/rgb]
gray = [100/rgb,100/rgb,100/rgb]


def plot_bathymetry_3d_profile(ncfile, figsize=(10,5), time_index=-1, ax=None, skip=5, length=1., scalez=1., veg='yes', wind='no', sep='no', vertical_angle = 90, horizontal_angle = 0.):
    
    with netCDF4.Dataset(ncfile, 'r') as ds:
        
        # get spatial dimensions and bed levels
        x = ds.variables['x'][:,:]
        y = ds.variables['y'][:,:]
        z = ds.variables['zb'][time_index,:,:]
        zs = ds.variables['zs'][30,:,:]
        if veg == 'yes':
            rhoveg = ds.variables['rhoveg'][time_index,:,:]
        else:
            rhoveg = 0.
            
        
    fig = plt.figure(figsize=figsize)
    ax = fig.gca(projection='3d')

    x_scale = 1.
    y_scale = 1.
    z_scale = .2
    
    scale=np.diag([x_scale, y_scale, z_scale, 1.0])
    scale=scale*(1.0/scale.max())
    scale[3,3]=1.0

    def short_proj():
        return np.dot(Axes3D.get_proj(ax), scale)

    ax.get_proj=short_proj

    # lay-out of axes
    ax.set_ylabel('alongshore distance [m]')
    ax.set_xlabel('cross-shore distance [m]')
    ax.axes.get_xaxis().set_ticks([])
    ax.axes.get_yaxis().set_ticks([])
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
#        ax.zaxis.set_visible(False)
    
    # perspective
    ax.view_init(vertical_angle,horizontal_angle)
    
    # make the panes transparent
    ax.xaxis.set_pane_color((0/rgb,166/rgb,214/rgb,0.0))
    ax.yaxis.set_pane_color((0/rgb,166/rgb,214/rgb,0.0))
    ax.zaxis.set_pane_color((0/rgb,166/rgb,214/rgb,0.0))
    # make the grid lines transparent
    ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0) 
    # Get rid of the spines                         
    ax.xaxis.line.set_color((1.0, 1.0, 1.0, 0.0)) 
    ax.yaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.line.set_color((1.0, 1.0, 1.0, 0.0)) 

    ax.set_axis_off()
        
    colorfac = 0.5-np.minimum(np.maximum((zs-z)*4, 0.),0.5)+0.5*rhoveg
        
    ax.plot_surface(y, x, z, facecolors=cm.get_cmap('Vegetation')(colorfac), linewidth=0, antialiased=False, shade=True, alpha=1.0, rstride=1, cstride=1)
    
    #plt.show()
    plt.savefig( '/Users/Simulations/parabolic/001_' + str('{:03}'.format(time_index)) + '.png', dpi=200)
    plt.close(fig)
        
    return rhoveg


def plot_bathymetry_3d(ncfile, figsize=(10,5), time_index=-1, ax=None, skip=5, length=1., scalez=1., veg='yes', wind='no', sep='no', vertical_angle = 90, horizontal_angle = 0.):
    
    with netCDF4.Dataset(ncfile, 'r') as ds:
        
        # get spatial dimensions and bed levels
        x = ds.variables['x'][:,:]
        y = ds.variables['y'][:,:]
        z = ds.variables['zb'][time_index,:,:]
        zs = -1. #ds.variables['zs'][time_index,:,:]
        if veg == 'yes':
            rhoveg = ds.variables['rhoveg'][time_index,:,:]
        else:
            rhoveg = 0.
#        zshear = ds.variables['zshear'][...]
    
#        uth = ds.variables['uth'][1,:,:]
#        uth_min = np.min(uth)
        
#        colors = np.empty(x.shape)
#        colors = np.repeat(colors[:,:,np.newaxis], 3, axis = 2)
#        
#        for j in range(y.shape[0]):
#            for i in range(x.shape[1]):
#                if z[time_index,j,i]<0.1 and veg=='yes':
#                    colors[j, i, :] = vegetation
#                else:
#                    colors[j, i, :] = sand

        # create figure
#        plt.ioff()
        
        fig = plt.figure(figsize=figsize)
        ax = plt.gca(projection='3d')

        x_scale = 1.
        y_scale = 1. #np.max(x)/np.max(y)
        z_scale = scalez*np.max(z)/np.max(x)
        
        scale=np.diag([x_scale, y_scale, z_scale, 1.0])
        scale=scale*(1.0/scale.max())
        scale[3,3]=1.0
        
        def short_proj():
          return np.dot(Axes3D.get_proj(ax), scale)
        
        ax.get_proj=short_proj

        # plot bed levels and bed level change        
        
        if wind == 'yes':
#        
            x_wind = x[::skip,::skip]
            y_wind = y[::skip,::skip]
            z_wind = z[::skip,::skip] + 2.
            
            u = ds.variables['ustars'][::skip,::skip]
            v = ds.variables['ustarn'][::skip,::skip]
            w = v#*0
#
            ax.quiver(x_wind, y_wind, z_wind[:,:], u[:,:], v[:,:], w[:,:], color=white, length=length, normalize=True)
        
        # lay-out of axes
        ax.set_ylabel('alongshore distance [m]')
        ax.set_xlabel('cross-shore distance [m]')
        ax.axes.get_xaxis().set_ticks([])
        ax.axes.get_yaxis().set_ticks([])
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
#        ax.zaxis.set_visible(False)
        
        # perspective
#        ax.view_init(90,0)
        ax.view_init(vertical_angle,horizontal_angle)
        
        # make the panes transparent
        ax.xaxis.set_pane_color((0/rgb,166/rgb,214/rgb,0.0))
        ax.yaxis.set_pane_color((0/rgb,166/rgb,214/rgb,0.0))
        ax.zaxis.set_pane_color((0/rgb,166/rgb,214/rgb,0.0))
        # make the grid lines transparent
        ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
        ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
        ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0) 
        # Get rid of the spines                         
        ax.xaxis.line.set_color((1.0, 1.0, 1.0, 0.0)) 
        ax.yaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
        ax.zaxis.line.set_color((1.0, 1.0, 1.0, 0.0)) 
        
#        ax.vert_exag = 10.
        ax.set_axis_off()
        pos_z = z.copy()
#        pos_veg = zs.copy() + 0.1
        neg_z = z.copy()
        
#        pos_z[pos_z <= 0.2] = np.nan
#        pos_veg[rhoveg < 0.1] = np.nan
        neg_z[neg_z > 0.05] = np.nan
        
#        colorfac = 0.5-np.minimum(np.maximum((zs-z)*4, 0.),0.5)+0.5*rhoveg
        colorfac = 0.5
        
        ax.plot_surface(y, x, pos_z[:,:], color=sand, linewidth=0, antialiased=False, shade=True, alpha=1.0, rstride=1, cstride=1)

        
#        if veg == 'yes':
#            ax.plot_surface(y, x, pos_z[:,:], facecolors=cm.get_cmap('Vegetation')(1-rhoveg), linewidth=0, antialiased=False, shade=True, alpha=1.0, rstride=1, cstride=1)
#        else:   
#            ax.plot_surface(y, x, pos_z[:,:], color = sand, linewidth=0, antialiased=False, shade=True, rstride=1, cstride=1)
#            
        ax.plot_surface(y, x, neg_z[:,:], color=white, linewidth=0, antialiased=False, shade=False, alpha=1.0, rstride=1, cstride=1)
        
#        ax.plot(y[35,:],x[35,:],zshear[time_index,35,:]+2.,color = white)
#        ax.plot(y[50,:],x[50,:],zshear[time_index,50,:]+2.,color = white)
#        ax.plot(y[65,:],x[65,:],zshear[time_index,65,:]+2.,color = white)
        
#        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
#        plt.close(fig)
        plt.savefig('/Users/Simulations/parabolic/001_' + str(time_index) + '.png', dpi=200)
        plt.close(fig)
        
    return ax


def plot_bathymetry_3d_save(ncfile, figsize=(10,5), time_index=-1, ax=None, skip=5, length=1., scalez=1., veg='yes', wind='no', sep='no', vertical_angle = 90, horizontal_angle = 0.):
    
    with netCDF4.Dataset(ncfile, 'r') as ds:
        
        # get spatial dimensions and bed levels
        x = ds.variables['x'][:,:]
        y = ds.variables['y'][:,:]
        z = ds.variables['zb'][time_index,:,:]
        
        fig = plt.figure(figsize=figsize)
        ax = plt.gca(projection='3d')

        x_scale = 1.
        y_scale = 1. #np.max(x)/np.max(y)
        z_scale = scalez*np.max(z)/np.max(x)
        
        scale=np.diag([x_scale, y_scale, z_scale, 1.0])
        scale=scale*(1.0/scale.max())
        scale[3,3]=1.0
        
        def short_proj():
          return np.dot(Axes3D.get_proj(ax), scale)
        
        ax.get_proj=short_proj   
        
        # lay-out of axes
        ax.set_ylabel('alongshore distance [m]')
        ax.set_xlabel('cross-shore distance [m]')
        ax.axes.get_xaxis().set_ticks([])
        ax.axes.get_yaxis().set_ticks([])
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        
        # perspective
        ax.view_init(vertical_angle,horizontal_angle)
        
        # make the panes transparent
        ax.xaxis.set_pane_color((0/rgb,166/rgb,214/rgb,0.0))
        ax.yaxis.set_pane_color((0/rgb,166/rgb,214/rgb,0.0))
        ax.zaxis.set_pane_color((0/rgb,166/rgb,214/rgb,0.0))
        # make the grid lines transparent
        ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
        ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
        ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0) 
        # Get rid of the spines                         
        ax.xaxis.line.set_color((1.0, 1.0, 1.0, 0.0)) 
        ax.yaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
        ax.zaxis.line.set_color((1.0, 1.0, 1.0, 0.0)) 
        
        ax.set_axis_off()
        pos_z = z.copy()
        neg_z = z.copy()
        neg_z[neg_z > 0.2] = np.nan
        
        ax.plot_surface(y, x, pos_z[:,:], color=sand, linewidth=0, antialiased=False, shade=True, alpha=1.0, rstride=1, cstride=1)          
        ax.plot_surface(y, x, neg_z[:,:], color=white, linewidth=0, antialiased=False, shade=False, alpha=1.0, rstride=1, cstride=1)
        
        #plt.savefig('/Users/Simulations/parabolic/001_' + str('{:03}'.format(time_index)) + '.png', dpi=200)
        
        plt.close(fig)
        
    return ax


def plot_bathymetry_3d_animate(ncfile, tstart = 0, tend = 260, tskip = 10 ):

    for i in range(tstart,tend):
        if i % tskip == 0:
            plot_bathymetry_3d_profile(ncfile, figsize=(20, ny/nx*20), time_index=i, ax=None, skip=1,
                                       length=10., scalez=4., veg='yes', wind='yes', sep='yes', vertical_angle=60, horizontal_angle=-10)
            plt.close()

#plot_bathymetry_3d('params13.nc', figsize=(15,8), time_index=480, ax=None, skip=16, length=0.5, scalez=4., veg='yes', wind='no', sep='no')
#
# plot_bathymetry_3d(ncfile, figsize=(20, ny/nx*20), time_index=-1, ax=None, skip=10,
#                    length=10., scalez=4., veg='yes', wind='no', sep='no', vertical_angle=90, horizontal_angle=0)
# plt.show()
#plot_bathymetry_3d('params13.nc', figsize=(15,8), time_index=480, ax=None, skip=16, length=0.5, scalez=4., veg='yes', wind='no', sep='no')
# plot_bathymetry_3d_save(ncfile, figsize=(20, ny/nx*20), time_index=-1, ax=None, skip=10,
#                    length=10., scalez=4., veg='no', wind='no', sep='no', vertical_angle=90, horizontal_angle=0)
plot_bathymetry_3d_profile(ncfile, figsize=(20, ny/nx*20), time_index=45, ax=None, skip=10,
                   length=10., scalez=4., veg='yes', wind='yes', sep='yes', vertical_angle=60, horizontal_angle=-10)

#plot_bathymetry_3d_animate(ncfile, tstart=0, tend=50, tskip=10)
