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
import math
#import matplotlib.pyplot as plt
from aeolis.wind import *

# package modules
import aeolis.wind
#from aeolis.utils import *

# initialize logger
logger = logging.getLogger(__name__)

def initialize (s,p):
    '''Initialise vegetation based on vegetation file.
              
    Parameters
    ----------
    s : dict
        Spatial grids
    p : dict
        Model configuration parameters
        
    Returns
    -------
    dict
        Spatial grids
        
    '''
    
    if p['veg_file'] is not None:
        s['rhoveg'][:, :] = p['veg_file']

        if np.isnan(s['rhoveg'][0, 0]):
            s['rhoveg'][:,:] = 0.

    ix = s['rhoveg'] < 0
    s['rhoveg'][ix] *= 0.
    s['hveg'][:,:] = p['hveg_max']*np.sqrt(s['rhoveg'])

    # Change dtype of array to bool for germinate and lateral
    s['vegetated'] = s['vegetated'].astype('bool') 
    s['germinate'] = s['germinate'].astype('bool') 
    s['lateral'] = s['lateral'].astype('bool') 

    # Fill these
    s['vegetated'] = (s['rhoveg']>0)
    s['germinate'][:,:] = False
    s['lateral'][:,:] = False

    return s


def vegshear(s, p):
    if p['vegshear_type'] == 'okin' and p['ny'] == 0:
        s = vegshear_okin(s, p)
    else:
        s = vegshear_raupach(s, p)

    s = velocity_stress(s,p)

    return s


def germinate(s,p):
    ny = p['ny']
    
    # time [year]
    n = (365.25*24.*3600. / (p['dt_opt'] * p['accfac']))

    # Compute dhveg_max to determine a minimum threshold of rhoveg to be considered "vegetated"
    dhveg_max = p['V_ver'] * (p['dt_opt'] * p['accfac']) / (365.25*24.*3600.)
    drhoveg_max = (dhveg_max/p['hveg_max'])**2

    # Determine which cells are already germinated before
    # s['germinate'][:, :] = False
    s['germinate'][:, :] = (s['rhoveg'] > 0.5 * drhoveg_max)

    # Germination (convert from /year to /timestep)
    # p_germinate_year = p['germinate']                                
    p_germinate_dt = 1-(1-p['germinate']  )**(1./n)

    # random_prob = np.zeros((s['germinate'].shape))
    random_prob = np.random.random((s['germinate'].shape))
    
    # Germinate new cells (only if no erosion, dzb>=0) and add to vegetated
    s['germinate'] = (s['dzbveg'] >= 0.) * (random_prob <= p_germinate_dt)
    s['vegetated'] = np.logical_or(s['germinate'], s['vegetated'])
    # s['germinate'] = np.minimum(s['germinate'], 1.)

    # Lateral expansion
    if ny > 1:
        dx = s['ds'][2,2]
    else:
        dx = p['dx']

    # p_lateral_year = p['lateral']  
    # Lateral propagation (convert from /year to /timestep)
    p_lateral_dt = 1-(1-p['lateral'] )**(1./n)
    p_lateral_cell = 1 - (1-p_lateral_dt)**(1./dx)
    
    drhoveg = np.zeros((p['ny']+1, p['nx']+1, 4))
    
    drhoveg[:,1:,0] = np.maximum((s['rhoveg'][:,:-1]-s['rhoveg'][:,1:]) / s['ds'][:,1:], 0.)    # positive x-direction
    drhoveg[:,:-1,1] = np.maximum((s['rhoveg'][:,1:]-s['rhoveg'][:,:-1]) / s['ds'][:,:-1], 0.)  # negative x-direction
    drhoveg[1:,:,2] = np.maximum((s['rhoveg'][:-1,:]-s['rhoveg'][1:,:]) / s['dn'][1:,:], 0.)    # positive y-direction
    drhoveg[:-1,:,3] = np.maximum((s['rhoveg'][1:,:]-s['rhoveg'][:-1,:]) / s['dn'][:-1,:], 0.)  # negative y-direction
    
    lat_veg = drhoveg > 0.
    
    s['drhoveg'] = np.sum(lat_veg[:,:,:], 2)
    
    p_lateral = p_lateral_cell * s['drhoveg']
    
    s['vegetated'] = np.logical_or((random_prob <= p_lateral), s['vegetated'])
    # s['lateral'] = np.minimum(s['lateral'], 1.)

    return s


def grow (s, p): #DURAN 2006
    
    # ix = np.logical_or(s['germinate'] != 0., s['lateral'] != 0.) * ( p['V_ver'] > 0.)
    # ix = np.logical_or(s['germinate'], s['lateral']) * ( p['V_ver'] > 0.)
    ix = s['vegetated'] * ( p['V_ver'] > 0.)

    # Vegetation growth
    vertical_growth = np.zeros(np.shape(s['hveg']))
    vertical_growth[ix] = p['V_ver'] * (1 - s['hveg'][ix] / p['hveg_max'])
    vertical_growth = apply_mask(vertical_growth, s['vver_mask'])

    # Reduction of vegetation growth due to sediment burial
    dhveg_burial = np.abs(s['dzbveg']-p['dzb_opt']) * p['veg_gamma']

    # Compute change in vegetation height
    s['dhveg'][ix] = vertical_growth[ix] - dhveg_burial[ix]  # m/year

    # Adding height and convert to vegetation density (rhoveg)
    if p['veggrowth_type'] == 'orig': #based primarily on vegetation height
        s['hveg'] += s['dhveg']*(p['dt_opt'] * p['accfac']) / (365.25*24.*3600.)
        s['hveg'] = np.maximum(np.minimum(s['hveg'], p['hveg_max']), 0.)
        s['rhoveg'] = (s['hveg']/p['hveg_max'])**2

    else:
        t_veg = p['t_veg']/365
        v_gam = p['v_gam']
        rhoveg_max = p['rhoveg_max']
        ix2 = s['rhoveg'] > rhoveg_max
        s['rhoveg'][ix2] = rhoveg_max
        ixzero = s['rhoveg'] <= 0
        if p['V_ver'] > 0:
            s['drhoveg'][ix] = (rhoveg_max - s['rhoveg'][ix])/t_veg - (v_gam/p['hveg_max'])*np.abs(s['dzbveg'][ix] - p['dzb_opt'])*p['veg_gamma']
        else:
            s['drhoveg'][ix] = 0
        s['rhoveg'] += s['drhoveg']*(p['dt']*p['accfac'])/(365.25 * 24 *3600)
        irem = s['rhoveg'] < 0
        s['rhoveg'][irem] = 0
        s['rhoveg'][ixzero] = 0 #here only grow vegetation that already existed
        #now convert back to height for Okin or wherever else needed
        s['hveg'][:,:] = p['hveg_max']*np.sqrt(s['rhoveg'])

    # Plot has to vegetate again after dying
    dhveg_max = p['V_ver'] * (p['dt_opt'] * p['accfac']) / (365.25*24.*3600.)
    drhoveg_max = (dhveg_max/p['hveg_max'])**2
    s['vegetated'] *= (s['rhoveg'] >= 0.5 * drhoveg_max)
    # s['germinate'] *= (s['rhoveg']!=0.)
    # s['lateral'] *= (s['rhoveg']!=0.)

    # Dying of vegetation due to hydrodynamics (Dynamic Vegetation Limit)
    # if p['process_tide']:
    #     s['rhoveg']     *= (s['zb'] +0.01 >= s['zs'])
    #     s['hveg']       *= (s['zb'] +0.01 >= s['zs'])
    #     s['germinate']  *= (s['zb'] +0.01 >= s['zs'])
    #     s['lateral']    *= (s['zb'] +0.01 >= s['zs'])

    if p['process_tide']:
        ix_flooded = (s['zb'] < s['TWL'])
        s['rhoveg'][ix_flooded]     = 0. 
        s['hveg'][ix_flooded]       = 0.
        s['vegetated'][ix_flooded]  = False
        # s['lateral'][ix_flooded]    = False

    ix = s['zb'] < p['veg_min_elevation']
    s['rhoveg'][ix] = 0
    s['hveg'][ix] = 0
    s['vegetated'][ix] = False
    # s['lateral'][ix] = False

    return s


def vegshear_okin(s, p):
    #Approach to calculate shear reduction in the lee of plants using the general approach of:
    #Okin (2008), JGR, A new model of wind erosion in the presence of vegetation
    #Note that implementation only works in 1D currently

    #Initialize shear variables and other grid parameters
    ustar = s['ustar'].copy()
    ustars = s['ustars'].copy()
    ustarn = s['ustarn'].copy()
    ets = np.zeros(s['zb'].shape)
    etn = np.zeros(s['zb'].shape)
    ix = ustar != 0
    ets[ix] = ustars[ix] / ustar[ix]
    etn[ix] = ustarn[ix] / ustar[ix]
    udir = s['udir'][0,0] + 180

    x = s['x'][0,:]
    zp = s['hveg'][0,:]
    red = np.zeros(x.shape)
    red_all = np.zeros(x.shape)
    nx = x.size
    c1 = p['okin_c1_veg']
    intercept = p['okin_initialred_veg']

    if udir < 0:
        udir = udir + 360

    if udir > 360:
        udir = udir - 360

    #Calculate shear reduction by looking through all cells that have plants present and looking downwind of those features
    for igrid in range(nx):

        if zp[igrid] > 0:         # only look at cells with a roughness element
            mult = np.ones(x.shape)
            h = zp[igrid] #vegetation height at the appropriate cell

            if udir >= 180 and udir <= 360:
                xrel = -(x - x[igrid])
            else:
                xrel = x - x[igrid]

            for igrid2 in range(nx):

                if xrel[igrid2] >= 0 and xrel[igrid2]/h < 20:

                    # apply okin model
                    mult[igrid2] = intercept + (1 - intercept) * (1 - np.exp(-xrel[igrid2] * c1 / h))

            red = 1 - mult

            # fix potential issues for summation
            ix = red < 0.00001
            red[ix] = 0
            ix = red > 1
            red[ix] = 1
            ix = xrel < 0
            red[ix] = 0

            # combine all reductions between plants
            red_all = red_all + red

    # cant have more than 100% reduction
    ix = red_all > 1
    red_all[ix] = 1

    #update shear velocity according to Okin (note does not operate on shear stress)
    mult_all = 1 - red_all
    ustarveg = s['ustar'][0,:] * mult_all
    ix = ustarveg < 0.01
    ustarveg[ix] = 0.01 #some small number so transport code doesnt crash

    s['ustar'][0,:] = ustarveg
    s['ustars'][0,:] = s['ustar'][0,:] * ets[0,:]
    s['ustarn'][0,:] = s['ustar'][0,:] * etn[0,:]

    return s


def vegshear_raupach(s, p):
    ustar = s['ustar'].copy()
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

    s['ustar'] *= s['vegfac']
    s['ustars'] = s['ustar'] * ets
    s['ustarn'] = s['ustar'] * etn

    return s
