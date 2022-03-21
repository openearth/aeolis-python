

import numpy as np
import logging
import netCDF4
from netCDF4 import Dataset
import matplotlib.pyplot as plt

# initialize logger
logger = logging.getLogger(__name__)

ncfile = 'C:/Users/Aeolis/dune development/parabolic/001/aeolis.nc'


def dz(outputfile1, ix):
    with netCDF4.Dataset(outputfile1, 'r') as ds:
        
        x = ds.variables['x'][:,50:]
        y = ds.variables['y'][:,50:]
        zb1 = ds.variables['zb'][ix,:,50:]
                
    dz = zb1 - z_2017
    
    a = np.maximum(-np.min(dz), np.max(dz))
    print (a)
    
    plt.figure(figsize=(12,10))
    mesh = plt.pcolormesh(x, y, dz, cmap='RdBu')
    mesh.set_clim(vmin=-2,vmax=2)
    bar = plt.colorbar(mesh)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('dzb')
    

def plot(outputfile, ix):
    
    with netCDF4.Dataset(outputfile, 'r') as ds:
        
        x = ds.variables['x'][:,:]
        y = ds.variables['y'][:,:]
        zb = ds.variables['zb'][ix,:,:]
        
        vmin = np.min(ds.variables['zb'][1,:,:])
        vmax = np.max(ds.variables['zb'][1,:,:])
    
        plt.figure(figsize=(10,8))
        mesh = plt.pcolormesh(x, y, zb, cmap='copper_r')
        mesh.set_clim(vmin,vmax)
        bar = plt.colorbar(mesh)
        bar.set_label('zb [m]')
        plt.xlabel('x [m]')
        plt.ylabel('y [m]')
        

def save_zb(outputfile,  ix):
    
    with netCDF4.Dataset(outputfile, 'r') as ds:
        

        zb = ds.variables['zb'][ix,:,:]
        
        
def Volume(outputfile, ix):
    
    with netCDF4.Dataset(outputfile, 'r') as ds:
        V0 =  np.sum(ds.variables['zb'][0,:,:])
        dV = np.sum(ds.variables['zb'][ix,:,:]) - np.sum(ds.variables['zb'][0,:,:])
        
        return V0, dV

def bedlevelchange(outputfile, ix):

    with netCDF4.Dataset(outputfile, 'r') as ds:

        x = ds.variables['x'][:,:] 
        y = ds.variables['y'][:,:]
        dzb = ds.variables['zb'][ix,:,:]-ds.variables['zb'][0,:,:]
        
        a = np.maximum(-np.min(dzb), np.max(dzb))
        
        fig, ax = plt.subplots(figsize=(10,8))
        mesh = ax.pcolormesh(x, y, dzb, cmap='RdBu');
        #mesh.set_clim(vmin=-a,vmax=a)
        mesh.set_clim(vmin=-4,vmax=4)
        bar = plt.colorbar(mesh)        
        bar.set_label('dzb [m]')
        #ax.set_xlim(0, 200)
        ax.set_xlabel('x [m]');
        ax.set_ylabel('y [m]');
        ax.set_title('Bed level change AeoLiS');
        
        return fig, ax

def plot_profiles(outputfile, y, ix):
    
    with netCDF4.Dataset(outputfile, 'r') as ds:
    
        nf = 0

        x = ds.variables['x'][y,:]
        
        zb0 = ds.variables['zb'][0,y,:]
        zb = ds.variables['zb'][ix,y,:]
        zsep = ds.variables['zsep'][ix,y,:]
        zs = ds.variables['zs'][ix,y,:]

        ustar0 = ds.variables['ustar0'][ix,y,:]
        ustar = ds.variables['ustar'][ix,y,:]
        uth = ds.variables['uth'][ix,y,:,nf]
    
        vegfac = ds.variables['vegfac'][ix,y,:]
        rhoveg = ds.variables['rhoveg'][ix,y,:]
        
        Cu = ds.variables['Cu'][ix,y,:,nf] * 1000
        Ct = ds.variables['Ct'][ix,y,:,nf] * 1000

    
        plt.figure(figsize=(12,10))
        plt.subplot(4,1,1)
        plt.plot(x, zs, 'cyan', label='zs')
        plt.plot(x, zb0, 'grey', label='zb0')
        plt.plot(x, zsep, 'lightblue', label='zsep')
        plt.plot(x, zb, 'orange', label='zb')
        plt.xlim(0, 225)
        #plt.ylim(-0.03, 0.1)
        plt.xlabel('x [m]');
        plt.ylabel('z [m]');
        plt.title('Bed level')
        plt.legend()
        
        plt.subplot(4,1,2)
        plt.plot(x, vegfac, 'palegreen', label='vegfac')
        plt.plot(x, rhoveg, 'seagreen', label='rhoveg')
        plt.xlim(0, 225)
        plt.ylim(-0.1, 1.1)
        plt.xlabel('x [m]');
        plt.ylabel('value');
        plt.title('Vegetation')
        plt.legend()
        
        plt.subplot(4,1,3)
        plt.plot(x, ustar0, 'pink', label='ustar0')
        plt.plot(x, ustar, 'coral', label='ustar')
        plt.plot(x, uth, 'lightgrey', label='uth')
        plt.xlim(0, 225)
        plt.xlabel('x [m]');
        plt.ylabel('ustar [m/s]');
        plt.title('Shear velocity profile')
        plt.legend()
        
        plt.subplot(4,1,4)
        plt.plot(x, Cu, 'grey', label='Cu');
        plt.plot(x, Ct, 'crimson', label='Ct')
        plt.xlim(0, 225)
        #plt.ylim(0,0.01)
        plt.xlabel('x [m]');
        plt.ylabel('C [g/m^2]');
        plt.title('Sediment concentration')
        plt.legend()
        
        plt.subplots_adjust(hspace=0.5)
        plt.show()
        

def flux(outputfile,ix):

    with netCDF4.Dataset(outputfile, 'r') as ds:
        
        x  = ds.variables['x'][:,:]
        y  = ds.variables['y'][:,:]
        zb = ds.variables['zb'][ix,:,:]
        Cu = ds.variables['Cu'][ix,:,:,0] * 1000
        Ct = ds.variables['Ct'][ix,:,:,0] * 1000
        
            
        n = 10
            
        vmin = 0
        vmax = 1. # np.max(Ct)
        
        fig, ax = plt.subplots(figsize=(10,8))
        mesh = ax.pcolormesh(x, y, Ct, cmap='pink_r');
        mesh.set_clim(vmin,vmax)
        bar = plt.colorbar(mesh)
        bar.set_label('Ct')
        ax.contour(x, y, zb, n, colors='black')
        ax.set_xlabel('x [m]');
        ax.set_ylabel('y [m]');
        ax.set_title('Sediment concentration');
            
        return fig, ax    

def vegetation(outputfile,ix):

    with netCDF4.Dataset(outputfile, 'r') as ds:
        
        x  = ds.variables['x'][:,:]
        y  = ds.variables['y'][:,:]
        zb = ds.variables['zb'][ix,:,:]
        rhoveg = ds.variables['rhoveg'][ix,:,:]
        vegfac = ds.variables['vegfac'][ix,:,:]
            
        n = 10
            
        vmin = 0
        vmax = 1.0
        
        fig, ax = plt.subplots(figsize=(10,8))
        mesh = ax.pcolormesh(x, y, rhoveg, cmap='YlGn');
        mesh.set_clim(vmin,vmax)
        bar = plt.colorbar(mesh)
        bar.set_label('rhoveg')
        ax.contour(x, y, zb, n, colors='black')
        ax.set_xlabel('x [m]');
        ax.set_ylabel('y [m]');
        #ax.set_title('Vegetation cover');
            
        return fig, ax
    
def dveg(outputfile,ix):
    with netCDF4.Dataset(outputfile, 'r') as ds:
        
        x  = ds.variables['x'][:,:]
        y  = ds.variables['y'][:,:]
        zb = ds.variables['zb'][ix,:,:]
        drhoveg = ds.variables['rhoveg'][ix,:,:] - ds.variables['rhoveg'][0,:,:]
            
        n = 10
            
        vmin = -1.0
        vmax = 1.0
        
        fig, ax = plt.subplots(figsize=(10,8))
        mesh = ax.pcolormesh(x, y, drhoveg, cmap='RdYlGn');
        mesh.set_clim(vmin,vmax)
        bar = plt.colorbar(mesh)
        bar.set_label('dveg')
        ax.contour(x, y, zb, n, colors='black')
        ax.set_xlabel('x [m]');
        ax.set_ylabel('y [m]');
        ax.set_title('Change in vegetation cover');
            
        return fig, ax    

def shear(outputfile,ix):

    with netCDF4.Dataset(outputfile, 'r') as ds:
        
        x  = ds.variables['x'][:,:]
        y  = ds.variables['y'][:,:]
        zb = ds.variables['zb'][ix,:,:]
        ustar0 = ds.variables['ustar'][ix,:,:]
        ustars = ds.variables['ustars'][ix,:,:]
        ustarn = ds.variables['ustarn'][ix,:,:]
        ustar = ds.variables['ustar'][ix,:,:]
            
        d = 10
        n = 10
            
        #vmin = 0
        #vmax = 1.0
        
        fig, ax = plt.subplots(figsize=(10,8))
        mesh = ax.pcolormesh(x, y, ustar, cmap='pink_r');
       #mesh.set_clim(vmin,vmax)
        bar = plt.colorbar(mesh)
        bar.set_label('zb')
        ax.contour(x, y, zb, n, colors='grey')
        plt.quiver(x[::d, ::d], y[::d, ::d], ustars[::d, ::d]/ustar0[::d, ::d], ustarn[::d, ::d]/ustar0[::d, ::d], color='black')
        ax.set_xlabel('x [m]');
        ax.set_ylabel('y [m]');
        ax.set_title('Shear velocities');
            
        return fig, ax   

def u(outputfile, ix):
    
    with netCDF4.Dataset(outputfile, 'r') as ds:
        
        x  = ds.variables['x'][:,:]
        y  = ds.variables['y'][:,:]
        zb = ds.variables['zb'][ix,:,:]

        #u0 = ds.variables['u0'][ix,:,:,0]
        us = ds.variables['us'][ix,:,:,0]
        un = ds.variables['un'][ix,:,:,0]
        u = ds.variables['u'][ix,:,:,0]
            
        d = 10
            
        #vmin = 0
        #vmax = 1.0
        
        fig, ax = plt.subplots(figsize=(12,7))
        mesh = ax.pcolormesh(x, y, zb, cmap='copper_r');
       #mesh.set_clim(vmin,vmax)
        bar = plt.colorbar(mesh)
        bar.set_label('zb')
        plt.quiver(x[::d, ::d], y[::d, ::d], us[::d, ::d], un[::d, ::d], color='white')
        ax.set_xlabel('x [m]');
        ax.set_ylabel('y [m]');
        ax.set_title('Grain velocities');
            
        return fig, ax
    
def wind(outputfile,ix):

    with netCDF4.Dataset(outputfile, 'r') as ds:
        
        x  = ds.variables['x'][:,:]
        y  = ds.variables['y'][:,:]
        zb = ds.variables['zb'][ix,:,:]
        uw = ds.variables['uw'][ix,:,:]
        uws = ds.variables['uws'][ix,:,:]
        uwn = ds.variables['uwn'][ix,:,:]

            
        d = 10
            
        #vmin = 0
        #vmax = 1.0
        
        fig, ax = plt.subplots(figsize=(10,8))
        mesh = ax.pcolormesh(x, y, zb, cmap='copper_r');
       #mesh.set_clim(vmin,vmax)
        bar = plt.colorbar(mesh)
        bar.set_label('zb')
        plt.quiver(x[::d, ::d], y[::d, ::d], uws[::d, ::d], uwn[::d, ::d], color='white')
        ax.set_xlabel('x [m]');
        ax.set_ylabel('y [m]');
        ax.set_title('Wind velocities');
            
        return fig, ax   