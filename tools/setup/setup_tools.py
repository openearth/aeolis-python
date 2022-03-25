import numpy as np
import os
import matplotlib.pyplot as plt


def create_grd(path, nx, ny, x0=0., y0=0., dx=1., dy=1., alfa=270., plotBool=False):
    '''
    Creates all required *.grd files for an AeoLiS simulation.
    
    z.grd, ne.grd, veg.grd are generated with zeros,
    so flat bed at z=0, flat non-erodible layer at z=0 and no vegetation 
    
    Parameters
    ----------
    path : str
        path of the folder where the *.grd files will be stored
    x0 : float
        x-coordinate of the base corner of the grid
    y0 : float
        y-coordinate of the base corner of the grid
    nx : int
        number of cells in x-direction
    ny : int
        number of cells in y-direction
    dx : float
        cellsize in x-direction
    dy : float
        cellsize in y-direction
    alfa : float
        rotation of the grid in nautical degrees
    '''
    
    
    # Convert from nautical to cartesian
    angle = np.deg2rad(270.0 - alfa)
    
    # Create arrays
    x_arr = np.arange(x0, x0 + dx*(nx+1), dx)
    y_arr = np.arange(y0, y0 + dy*(ny+1), dy)
    
    # Create meshgrid
    x, y = np.meshgrid(x_arr, y_arr)
    
    # Create z, ne and veg
    z = np.zeros(np.shape(x))
    ne = np.zeros(np.shape(x))
    veg = np.zeros(np.shape(x))
    
    # Rotate meshgrid
    xr = x0 + np.cos(angle) * (x - x0) - np.sin(angle) * (y - y0)
    yr = y0 + np.sin(angle) * (x - x0) + np.cos(angle) * (y - y0)
    
    # Saving files
    np.savetxt(os.path.join(path, 'x.grd'), xr)
    np.savetxt(os.path.join(path, 'y.grd'), yr)
    np.savetxt(os.path.join(path, 'z.grd'), z)
    np.savetxt(os.path.join(path, 'ne.grd'), ne)
    np.savetxt(os.path.join(path, 'veg.grd'), veg)
    
    # Plotting
    if plotBool:
        fig, ax = plt.subplots()
        ax.pcolormesh(xr, yr, z)
        ax.set_aspect('equal')
    
    
    return



def create_z_cone(path, nx, ny, height, diameter, midx, midy, dx=1., dy=1.):
    
    '''
    Creates a z.grd file with a cone shape on top of a flat bed 
    
    Parameters
    ----------
    path : str
        path of the folder where the *.grd files will be stored
    nx : integer
        description
    ny : integer
        description
    height : int
        description
    diameter : int
        description
    midx : float
        description
    midy : float
        description
    dx : float
        description
    dy : float
        description
    '''
    
    z=np.zeros((ny+1, nx+1))

    
    for j in range(0,nx+1):  
        for i in range(0,ny+1):
            if np.sqrt((midy-j*dy)**2+(midx-i*dx)**2)<=diameter/2:
                z[i,j]+=(1-(np.sqrt((midx-i*dx)**2+(midy-j*dy)**2)/(diameter/2)))*height
                
                
    np.savetxt(os.path.join(path, 'z.grd'), z)