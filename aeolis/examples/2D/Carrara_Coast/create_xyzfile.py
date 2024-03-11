# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 14:21:52 2019

@author: tspak
"""

import numpy as np
import matplotlib.pyplot as plt

"""
create spatial grid files (x.txt, y.txt, z.txt and ne_grid.txt) for AeoLiS

input parameters:
            dx (float): resolution in cross-shore direction (m)
            dy (float): resolution in alongshore direction (m)
            lx (float): length of domain in cross-shore direction (m) (must be a multiple of dx)
            ly (float): length of domain in alongshore direction (m) (must be a multiple of dy)
    
         z_top (float): height of the top of the nourishment (m + NAP)
         slope (float): slope of the nourishment (m/m)
    
    ly_seaside (float): length of the nourishment at the seaside boundary (m) (must be a multiple of dy)
   ly_landside (float): length of the nourishment at the landside boundary (m) (must be a multiple of dy)
     lx_street (float): width of the (non-erodible) street (m) (must be a multiple of dx)
lx_nourishment (float): width of the nourishment (m) (must be a multiple of dx)
"""
### input 
dx = 10 #m
dy = 10 #m
lx = 550 #m
ly = 1700 #m
z_top = 2 #m + NAP
slope = 0.02 #m/m
ly_seaside = 800 #m
ly_landside = 1500 #m
lx_street = 50 #m
lx_nourishment = 300 #m
### input

# define numerical parameters
nx = int(lx/dx) + 1                         #number of cross-shore grid points
ny = int(ly/dy) + 1                         #number of alongshore grid points
ny_seaside = int(ly_seaside/dy) + 1         #number of alongshore grid points of the nourishment at the seaside boundary
ny_landside = int(ly_landside/dy) + 1       #number of alongshore grid points of the nourishment at the landside boundary 
ny_notop_seaside = int((ny-ny_seaside)/2)   #number of alongshore grid points next to the nourishment at the seaside boundary
ny_notop_landside = int((ny-ny_landside)/2) #number of alongshore grid points next to the nourishment at the landside boundary 
nx_street = int(lx_street/dx)               #number of cross-shore grid points of the street
nx_nourishment = int(lx_nourishment/dx)     #number of cross-shore grid points of the nourishment

# create 2D grids for X, Y and Z
x = np.linspace(0, lx, nx)
y = np.linspace(0, ly, ny)
Y, X = np.meshgrid(y, x)
Z = np.zeros(X.shape)

# set values values to z_top at the middle of the nourishment
Z[:nx_street+nx_nourishment, ny_notop_seaside:ny_notop_seaside+ny_seaside] = z_top

# gradually decrease depth at the middle of the nourishment
values = -slope * (X - (nx_street+nx_nourishment)*dx) + z_top
ix = [nx_street+nx_nourishment, len(X)]
iy = [ny_notop_seaside,ny_notop_seaside+ny_seaside]
Z[ix[0]:ix[1], iy[0]:iy[1]] = values[ix[0]:ix[1], iy[0]:iy[1]]

# set values values to z_top at the street
Z[:nx_street, :] = z_top

# gradually decrease depth next to the nourishment
values = -slope * (X - nx_street*dx) + z_top
Z[nx_street:, :ny_notop_landside] = values[nx_street:, :ny_notop_landside]
Z[nx_street:, -ny_notop_landside:] = values[nx_street:, -ny_notop_landside:]

# gradually decrease depth at the sides of the nourishment
# create a temporal domain for the sides of the nourishment
grid = np.indices((nx-nx_street,ny_notop_seaside-ny_notop_landside))
# calculate slope of the side of the nourishment
dydx = nx_nourishment / (ny_notop_seaside-ny_notop_landside)
ix = grid[0] - grid[1]*dydx
values = -slope * (grid[0] - grid[1]*dydx) * dx + z_top
# round values to a multiple of (slope*dx)
roundup_value = 1/(slope*dx)
values = np.round(values * roundup_value) / roundup_value
values[values > z_top] = z_top

Z[nx_street:, ny_notop_landside:ny_notop_seaside] = values
Z[nx_street:, -ny_notop_seaside:-ny_notop_landside] = values[:,::-1]

# create non-erodible grid cells at the street 
ne_file = np.ones(Z.T.shape) * -999
ne_file[:,-nx_street:] = z_top

# plot the grid
plt.figure()
plt.pcolormesh(Z)
plt.colorbar()
plt.axis('equal')
plt.grid()

# save textfiles containing the grid and the ne_file
np.savetxt('x.txt', X.T)
np.savetxt('y.txt', Y.T)
np.savetxt('z.txt', Z[::-1,:].T)
np.savetxt('street.txt', ne_file)