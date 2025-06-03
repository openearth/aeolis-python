# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 17:58:28 2022

@author: cijzendoornvan
"""

# Load packages
import sys 
import netCDF4
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Introduce functions
def calculateD50(t, y, x, layers, fractions, mass):
    # Calculate D50 for each cell based on weight distribution over the fractions
    d50 = np.zeros((len(t), len(y[:,0]), len(x[0,:]), len(layers)))
       
    for idx, t_i in enumerate(t):    
        for yi, y_id in enumerate(y[:,0]):
            for xi, x_id in enumerate(x[0,:]):
                for ni, n in enumerate(layers):
                    d50[idx, yi, xi, ni] = np.asarray([np.average(fractions, weights=mass[idx, yi, xi, ni, :])])*1e6
    return d50

def animate(frame_num):
    frac1.set_data((x[0,:], pickup[frame_num,0,:,0]))
    frac2.set_data((x[0,:], pickup[frame_num,0,:,1]))
    
    bed.set_data(x[0,:], zb[frame_num,0,:])
    
    newdata = d50_all[frame_num,0,:,:].T
    bedplot.set_array((newdata.flatten()))
    
    time.set_text('time = ' + str(t[frame_num]) + ' s')
    
    wind.set_xdata((t[frame_num]))
    flux.set_xdata((t[frame_num]))
    
    return frac1, frac2, bedplot, time, wind, flux,

def prep_visualize():
    # Prep visualisation
    S = 14
    M = 18
    L = 20
    
    plt.rc('font', size=S)          # controls default text sizes
    plt.rc('axes', titlesize=S)     # fontsize of the axes title
    plt.rc('axes', labelsize=M)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=S)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=S)    # fontsize of the tick labels
    plt.rc('legend', fontsize=S)    # legend fontsize
    plt.rc('figure', titlesize=L)  # fontsize of the figure title
    
      
# Initialize directories
model_directory = "C:/Users/cijzendoornvan/OneDrive - Delft University of Technology/Documents/DuneForce/AEOLIS/aeolis-python/examples/grainsizevariations"
aeolis_directory = r"C:\Users\cijzendoornvan\OneDrive - Delft University of Technology\Documents\DuneForce\AEOLIS\aeolis-python\aeolis"
sys.path.append(aeolis_directory)
from console_debug import *

case1 = "aeolis_horizontalgradient1"
case2 = "aeolis_horizontalgradient2"

#%% Run models (note that the output file is saved, so unless input parameters are changed, this part of the code does not have to be run again)
# Run aeolis with case 1
case1_dir = model_directory + "/" + case1 + ".txt" 
aeolis_debug(case1_dir)

# Run aeolis with case 2
case2_dir = model_directory + "/" + case2 + ".txt"        
aeolis_debug(case2_dir)

#%% Visualization
prep_visualize()

#%% Visualisation of case 1 with coarse layer in top layer

# Load parameters resulting from model run
with netCDF4.Dataset(os.path.join(model_directory, case1+'.nc'), 'r') as ds:
    t = ds.variables['time'][:]# / 3600.
    uw = ds.variables['uw'][...]
    x = ds.variables['x'][:,:]
    y = ds.variables['y'][:,:]
    zb = ds.variables['zb'][:]
    pickup = ds.variables['pickup'][...]
    mass = ds.variables['mass'][:,:,:,:,:]
    qs = ds.variables['qs'][:]
    layers = ds.variables['layers'][:]
    fractions = ds.variables['fractions'][:]
          
# Calculate d50
d50_all = calculateD50(t, y, x, layers, fractions, mass)
    
# Set plotting parameters
vmin = 250#np.min(d50_all[:,0,:,:]) 
vmax = 500#np.max(d50_all[:,0,:,:])

frac1 = str(int(fractions[0]*1e6)) + ' $\mu$m'
frac2 = str(int(fractions[1]*1e6)) + ' $\mu$m'

# Calculate cumulative flux of both fractions
qs_sum = qs[:, 0, :, 0] + qs[:, 0, :, 1]
qs_cumsum = np.cumsum(qs_sum[:,-2])
                
# Create animation
fig, ax = plt.subplots(5, 1, figsize=(8,10))
fig.subplots_adjust(left = 0.16, right=0.75, hspace=1.1)

qs1, = ax[1].plot(t[:], np.cumsum(qs[:, 0, -2, 0]), 'goldenrod', label=frac1)
qs2, = ax[1].plot(t[:], np.cumsum(qs[:, 0, -2, 1]), 'maroon', label=frac2)
qs3, = ax[1].plot(t[:], qs_cumsum[:], 'dimgrey', label='summed')

ax[1].set_title('Sediment flux per fraction', fontsize=20)
ax[1].legend(bbox_to_anchor=(1.02, 1.03))
#ax[1].set_ylim(-np.max(qs_sum)*0.2, np.max(qs_sum)*1.1)
ax[1].set_xlabel('Time (s)')
ax[1].set_xlim(0, 600)
ax[1].set_ylabel('flux\n(kg/m/s)')
flux = ax[1].axvline(0, color='k', lw=2)

frac1, = ax[2].plot(x[0,:], pickup[0,0,:,0], 'goldenrod', label = frac1)
frac2, = ax[2].plot(x[0,:], pickup[0,0,:,1], 'maroon', label = frac2)

ax[2].set_title('Sediment pickup per fraction', fontsize=20)
#ax[1].legend()
ax[2].set_xlim(-1, 99)
ax[2].set_ylim(-np.max(pickup)*0.2, np.max(pickup)*1.1)
ax[2].set_ylabel('pickup\n(kg/m/s)')

bed, = ax[3].plot(x[0,:], zb[0,0,:], 'k')

ax[3].set_title('Bed level change', fontsize=20)
ax[3].set_xlim(-1, 99)
ax[3].set_ylim(np.min(zb)*1.6, np.min(zb)*-0.2)
ax[3].set_ylabel('zb (m)')

bedplot = ax[4].pcolor(d50_all[0,0,:,:].T, edgecolors='k', linewidths = 0.05, cmap='hot_r', vmin=vmin, vmax=vmax)
ax[4].invert_yaxis()
ax[4].set_xlim(-1, 99)
cax = fig.add_axes([0.78, 0.11, 0.03, 0.15])
fig.colorbar(bedplot, cax=cax, orientation='vertical')
cax.set_ylabel('gs ($\mu$m)', fontsize = 18)
ax[4].set_title('Average grain size in the bed', fontsize=20)
ax[4].set_ylabel('bed\nlayers')
ax[4].set_xlabel('Distance (m)')

axtext = fig.add_axes([0.75,0.83,0.1,0.15])
axtext.axis("off")
time = axtext.text(0.5,0.5, str(0), ha="left", va="top")

ax[0].plot(t[:], uw[:,0,0], 'royalblue')
ax[0].set_ylabel('wind speed\n(m/s)')
ax[0].set_xlabel('Time (s)')
ax[0].set_xlim(0, 600)
ax[0].set_ylim(0, 30)
wind = ax[0].axvline(0, color='k', lw=2)

anim = animation.FuncAnimation(fig, animate, frames=len(t)-1, interval=30)
writergif = animation.PillowWriter(fps=30) 
anim.save(case1+'_overview.gif', writer=writergif)
plt.show()

#%% Visualisation of case 2 with coarse layer in fourth layer

# Load parameters resulting from model run
with netCDF4.Dataset(os.path.join(model_directory, case2+'.nc'), 'r') as ds:
    t = ds.variables['time'][:]# / 3600.
    uw = ds.variables['uw'][...]
    x = ds.variables['x'][:,:]
    y = ds.variables['y'][:,:]
    zb = ds.variables['zb'][:]
    pickup = ds.variables['pickup'][...]
    mass = ds.variables['mass'][:,:,:,:,:]
    qs = ds.variables['qs'][:]
    layers = ds.variables['layers'][:]
    fractions = ds.variables['fractions'][:]
          
# Calculate d50
d50_all = calculateD50(t, y, x, layers, fractions, mass)
    
# Set plotting parameters
vmin = 250#np.min(d50_all[:,0,:,:]) 
vmax = 500#np.max(d50_all[:,0,:,:])

frac1 = str(int(fractions[0]*1e6)) + ' $\mu$m'
frac2 = str(int(fractions[1]*1e6)) + ' $\mu$m'

# Calculate cumulative flux of both fractions
qs_sum = qs[:, 0, :, 0] + qs[:, 0, :, 1]
qs_cumsum = np.cumsum(qs_sum[:,-2])
                
# Create animation
fig, ax = plt.subplots(5, 1, figsize=(8,10))
fig.subplots_adjust(left = 0.16, right=0.75, hspace=1.1)

qs1, = ax[1].plot(t[:], np.cumsum(qs[:, 0, -2, 0]), 'goldenrod', label=frac1)
qs2, = ax[1].plot(t[:], np.cumsum(qs[:, 0, -2, 1]), 'maroon', label=frac2)
qs3, = ax[1].plot(t[:], qs_cumsum[:], 'dimgrey', label='summed')

ax[1].set_title('Sediment flux per fraction', fontsize=20)
ax[1].legend(bbox_to_anchor=(1.02, 1.03))
#ax[1].set_ylim(-np.max(qs_sum)*0.2, np.max(qs_sum)*1.1)
ax[1].set_xlabel('Time (s)')
ax[1].set_xlim(0, 600)
ax[1].set_ylabel('flux\n(kg/m/s)')
flux = ax[1].axvline(0, color='k', lw=2)

frac1, = ax[2].plot(x[0,:], pickup[0,0,:,0], 'goldenrod', label = frac1)
frac2, = ax[2].plot(x[0,:], pickup[0,0,:,1], 'maroon', label = frac2)

ax[2].set_title('Sediment pickup per fraction', fontsize=20)
#ax[1].legend()
ax[2].set_xlim(-1, 99)
ax[2].set_ylim(-np.max(pickup)*0.2, np.max(pickup)*1.1)
ax[2].set_ylabel('pickup\n(kg/m/s)')

bed, = ax[3].plot(x[0,:], zb[0,0,:], 'k')

ax[3].set_title('Bed level change', fontsize=20)
ax[3].set_xlim(-1, 99)
ax[3].set_ylim(np.min(zb)*1.6, np.min(zb)*-0.2)
ax[3].set_ylabel('zb (m)')

bedplot = ax[4].pcolor(d50_all[0,0,:,:].T, edgecolors='k', linewidths = 0.05, cmap='hot_r', vmin=vmin, vmax=vmax)
ax[4].invert_yaxis()
ax[4].set_xlim(-1, 99)
cax = fig.add_axes([0.78, 0.11, 0.03, 0.15])
fig.colorbar(bedplot, cax=cax, orientation='vertical')
cax.set_ylabel('gs ($\mu$m)', fontsize = 18)
ax[4].set_title('Average grain size in the bed', fontsize=20)
ax[4].set_ylabel('bed\nlayers')
ax[4].set_xlabel('Distance (m)')

axtext = fig.add_axes([0.75,0.83,0.1,0.15])
axtext.axis("off")
time = axtext.text(0.5,0.5, str(0), ha="left", va="top")

ax[0].plot(t[:], uw[:,0,0], 'royalblue')
ax[0].set_ylabel('wind speed\n(m/s)')
ax[0].set_xlabel('Time (s)')
ax[0].set_xlim(0, 600)
ax[0].set_ylim(0, 30)
wind = ax[0].axvline(0, color='k', lw=2)

anim = animation.FuncAnimation(fig, animate, frames=len(t)-1, interval=30)
writergif = animation.PillowWriter(fps=30) 
anim.save(case2+'_overview.gif', writer=writergif)
plt.show()

