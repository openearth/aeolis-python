# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 11:24:34 2020

@author: tspak
"""

import os
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

"""
plot AeoLiS model results

input parameters:
 filedir (str): directory of the .nc file containing the AeoLiS model result
filename (str): fiilename of the .nc file containing the AeoLiS model result
            
"""

# Specifify path default is current working directory
filedir = os.getcwd()
# Specifify filename for output
filename = 'aeolis_v2noth_m.nc'

filepath = os.path.join(filedir, filename)
data = nc.Dataset(filepath)
t = data['time']
x = data['x']
y = data['y']
zb = data['zb']
dzb = zb[-1] - zb[0]
#uw = data['uw']
uws = data['uws']
#uth = data['uth']
qs = data['qs']
w = data['w']
w_bed = data['w_bed']
ny, nx = zb[0].shape
Qs = np.multiply(qs[:], w[:]).sum(axis=3) #summate sediment discharge of all sediment fractions

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
p1 = ax1.plot(t[:]/86400., Qs[:,int(ny/2),-2], 'r', label='Q [kg/ms]')[0]
p2 = ax2.plot(t[:]/86400., uws[:,int(ny/2),-1], 'b', label='uw [m/s]')[0]
ax1.set_ylim(-12/10000, 40/10000)
ax2.set_ylim(-12, 40)
ax1.set_xlabel('t [days]')
ax1.set_ylabel('Q [kg/ms]')
ax2.set_ylabel('uw [m/s]')
plt.legend([p1, p2], ['Q [kg/ms]', 'uw [m/s]'])
plt.grid(which='both', axis='both')
plt.title(filename)

plt.figure()
plt.plot(x[int(ny/2),:], dzb[int(ny/2),:])
plt.ylim(-0.0015, 0.0001)
plt.xlabel('x [m]')
plt.ylabel('$\Delta$ z [m]')
plt.grid()
plt.title(filename)

plt.figure()
plt.pcolormesh(x[:], y[:], zb[0])#, vmin=-0.3, vmax=0.3, cmap='RdBu_r')
plt.colorbar()
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.grid()
plt.title(filename)

#print(data['tau'])
#print(data['ustar'])


#print(data['tau'][:,85,10])
#print(data['ustar'][:,85,10])
#print(data['uw'][:,85,10])
#print(data['Cu'][50,85,:,7])
##print(data['Cu'][50,85,:,7])
#tau = (.41 / np.log(10/ 0.001)) * data['uw'][:,85,10]
#print(tau)
print(data['uth'][50,85,40,:])
#print(dzb[int(ny/2),:].shape)

data.close()