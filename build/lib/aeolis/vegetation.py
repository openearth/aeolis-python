'''This file is part of AeoLiS.
   
AeoLiS is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
   
AeoLiS is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
   
You should have received a copy of the GNU General Public License
along with AeoLiS.  If not, see <http://www.gnu.org/licenses/>.
   
AeoLiS  Copyright (C) 2015 Bas Hoonhout

bas.hoonhout@deltares.nl         b.m.hoonhout@tudelft.nl
Deltares                         Delft University of Technology
Unit of Hydraulic Engineering    Faculty of Civil Engineering and Geosciences
Boussinesqweg 1                  Stevinweg 1
2629 HVDelft                     2628CN Delft
The Netherlands                  The Netherlands

'''
from __future__ import absolute_import, division

import logging
from scipy import ndimage, misc
import numpy as np
import matplotlib.pyplot as plt

# package modules
import aeolis.wind
#from aeolis.utils import *

# initialize logger
logger = logging.getLogger(__name__)

def initialize (s,p):
    
    if p['veg_file'] is not None:
        s['rhoveg'][:, :] = p['veg_file']

        if np.isnan(s['rhoveg'][0, 0]):
            s['rhoveg'][:,:] = 0.

        ix = s['rhoveg'] < 0.03
        s['rhoveg'][ix] *= 0.
        s['hveg'][:,:] = p['hveg_max']*np.sqrt(s['rhoveg'])

        s['germinate'][:,:] = (s['rhoveg']>0)
        s['lateral'][:,:] = 0.

    return s


def vegshear(s, p):
    
    ustar  = s['ustar'].copy()
    ustars = s['ustars'].copy()
    ustarn = s['ustarn'].copy()
    
    ets = np.zeros(s['zb'].shape)
    etn = np.zeros(s['zb'].shape)
    
    ix = ustar != 0
    
    ets[ix] = ustars[ix] / ustar[ix]
    etn[ix] = ustarn[ix] / ustar[ix]
    

    # Raupach, 1993
    roughness = p['gamma_vegshear']
    
    vegfac = 1. / np.sqrt(1. + roughness * s['rhoveg'])
    
    # Smoothen the change in vegfac between surrounding cells following a gaussian distribution filter 

    s['vegfac'] = ndimage.gaussian_filter(vegfac, sigma=p['veg_sigma'])
    
    
    # Apply reduction factor of vegetation to the total shear stress 
    
    s['ustar']  *= s['vegfac']
    s['ustars'] = s['ustar'] * ets
    s['ustarn'] = s['ustar'] * etn
    
    return s


def germinate (s,p):
    ny = p['ny']
    s['germinate'][:, :] = (s['rhoveg'] > 0.)
    
    # time [year]
    n = (365.25*24.*3600. / (p['dt'] * p['accfac']))
    
    # Germination
    
    p_germinate_year = p['germinate']                                
    p_germinate_dt = 1-(1-p_germinate_year)**(1./n)
    germination = np.random.random((s['germinate'].shape))
    
    s['germinate'] += (s['dzbveg'] >= 0.) * (germination <= p_germinate_dt)
    s['germinate'] = np.minimum(s['germinate'], 1.)


    # Lateral expension
    if ny > 1:
        dx = s['ds'][2,2]
    else:
        dx = p['dx']

    p_lateral_year = p['lateral']  
    p_lateral_dt = 1-(1-p_lateral_year)**(1./n)
    p_lateral_cell = 1 - (1-p_lateral_dt)**(1./dx)
    
    drhoveg = np.zeros((p['ny']+1, p['nx']+1, 4))
    
    drhoveg[:,1:,0] = np.maximum((s['rhoveg'][:,:-1]-s['rhoveg'][:,1:]) / s['ds'][:,1:], 0.)    # positive x-direction
    drhoveg[:,:-1,1] = np.maximum((s['rhoveg'][:,1:]-s['rhoveg'][:,:-1]) / s['ds'][:,:-1], 0.)  # negative x-direction
    drhoveg[1:,:,2] = np.maximum((s['rhoveg'][:-1,:]-s['rhoveg'][1:,:]) / s['dn'][1:,:], 0.)    # positive y-direction
    drhoveg[:-1,:,3] = np.maximum((s['rhoveg'][1:,:]-s['rhoveg'][:-1,:]) / s['dn'][:-1,:], 0.)  # negative y-direction
    
    lat_veg = drhoveg > 0.
    
    s['drhoveg'] = np.sum(lat_veg[:,:,:], 2)
    
    p_lateral = p_lateral_cell * s['drhoveg']
    
    s['lateral'] += (germination <= p_lateral)
    s['lateral'] = np.minimum(s['lateral'], 1.)

    return s

def grow (s, p): #DURAN 2006
    
    ix = np.logical_or(s['germinate'] != 0., s['lateral'] != 0.) * ( p['V_ver'] > 0.)
                                                    

    # Reduction of vegetation growth due to sediment burial
    s['dhveg'][ix] = p['V_ver'] * (1 - s['hveg'][ix] / p['hveg_max']) - np.abs(s['dzbveg'][ix]-p['dzb_opt']) * p['veg_gamma']  # m/year


    # if p['_time'] > p['dzb_interval']:
        # plt.pcolormesh(s['x'],s['y'], s['dhveg'], vmin=-1, vmax=1, cmap='YlGn')
        # plt.gca().invert_yaxis()
        # bar = plt.colorbar()
        # bar.set_label('rhoveg')
        # plt.xlabel('x [m]')
        # plt.ylabel('y [m]')
        # plt.title('Vegetation cover')
        # plt.show()

    # Adding growth
    s['hveg'] += s['dhveg']*(p['dt'] * p['accfac']) / (365.25*24.*3600.)

    # Compute the density

    s['hveg'] = np.maximum(np.minimum(s['hveg'], p['hveg_max']), 0.)
    s['rhoveg'] = (s['hveg']/p['hveg_max'])**2
    
    
    # Plot has to vegetate again after dying

    s['germinate'] *= (s['rhoveg']!=0.)
    s['lateral'] *= (s['rhoveg']!=0.)

    # Dying of vegetation due to hydrodynamics (Dynamic Vegetation Limit)

    s['rhoveg']     *= (s['zb'] +0.01 >= s['zs'])
    s['hveg']       *= (s['zb'] +0.01 >= s['zs'])
    s['germinate']  *= (s['zb'] +0.01 >= s['zs'])
    s['lateral']    *= (s['zb'] +0.01 >= s['zs'])
    
    return s