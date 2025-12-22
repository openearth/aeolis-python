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
        if p['process_vegflex']:
            print('VefFlex = T -> Veg is not initialized in def initialize')
        else:
            print('VefFlex = F -> Veg is initialized in def initialize')
            s['rhoveg'][:, :] = p['veg_file']
    
            if np.isnan(s['rhoveg'][0, 0]):
                s['rhoveg'][:,:] = 0.
    
            ix = s['rhoveg'] < 0
            s['rhoveg'][ix] *= 0.
            s['hveg'][:,:] = p['hveg_max']*np.sqrt(s['rhoveg'])
    
            s['germinate'][:,:] = (s['rhoveg']>0)
            s['lateral'][:,:] = 0.

    return s

def VegPositionWrite(s,p):
    # Check if 1D or 2D
    if p['ny'] == 0:
        print('VegPositionWrite1D')
        # Computation of the DownScaling factor
        # -------------------------------------
        ResolutionVeg = 0.001
        Resolution = np.around(p['xgrid_file'][1] - p['xgrid_file'][0], decimals=2)
        DownScale = np.round(Resolution / ResolutionVeg)
        nbCell = DownScale * DownScale
        # Computation of the surface according to the method
        # ---------------------------------------------------
        # No Wind Speed: Ammophila arenaria = 11411.19 mm² | Ammophila breviliglata = 22225.37 mm² | Leymus mollis = 35952.50 mm²
        VegSurface = ((p['stem_diam'] * p['stem_height']) + (
                    (p['leaf_width'] * p['leaf_height_adjusted']) * p['leafToStem'])) * 1000000  # mm²
        NbPlante = np.ceil((p['veg_file']*nbCell)/VegSurface)
        
        # Find index (x) of cell with vegetation
        #-----------------------------------------
        indices = np.where(NbPlante > 0)
        posX = indices[0]  # Index in x
    
        # Repeat the positions as many times as the number of plants in each cell
        #------------------------------------------------------------------------
        PlanteRepeat = NbPlante[indices].astype(np.int64)  # The cell values represent the number of plants
        VegPosXUpS = np.repeat(posX, PlanteRepeat)
        
        # Conversion of VegPos depending on the resolution to the 1mm resolution
        #-----------------------------------------------------------------------
        VegPosX = (p['xgrid_file'][VegPosXUpS] - p['xgrid_file'][0]) / ResolutionVeg
        
        # Saving the results
        #-------------------
        p['veg_posx'] = VegPosX 
        
    else:
        print('VegPositionWrite2D')
        # Computation of the DownScaling factor
        # -------------------------------------
        ResolutionVeg = 0.001
        Resolution = np.around(p['xgrid_file'][0, 1] - p['xgrid_file'][0, 0], decimals=1)
        DownScale = np.round(Resolution / ResolutionVeg)
        nbCell = DownScale * DownScale
        # Computation of the surface according to the method
        # ---------------------------------------------------
        # No Wind Speed: Ammophila arenaria = 11411.19 mm² | Ammophila breviliglata = 22225.37 mm² | Leymus mollis = 35952.50 mm²
        VegSurface = ((p['stem_diam'] * p['stem_height']) + (
                    (p['leaf_width'] * p['leaf_height_adjusted']) * p['leafToStem'])) * 1000000  # mm²
        NbPlante = np.ceil((p['veg_file']*nbCell)/VegSurface)
        
        # Find index (x,y) of cell with vegetation
        #-----------------------------------------
        indices = np.where(NbPlante > 0)
        posX = indices[1]  # Index in x
        posY = indices[0]  # Index in y
    
        # Repeat the positions as many times as the number of plants in each cell
        #------------------------------------------------------------------------
        PlanteRepeat = NbPlante[indices].astype(np.int64)  # The cell values represent the number of plants
        VegPosXUpS = np.repeat(posX, PlanteRepeat)
        VegPosYUpS = np.repeat(posY, PlanteRepeat)
        
        # Conversion of VegPos depending on the resolution to the 1mm resolution
        #-----------------------------------------------------------------------
        VegPosX = (p['xgrid_file'][0,VegPosXUpS] - p['xgrid_file'][0,0]) / ResolutionVeg
        VegPosY = (p['ygrid_file'][VegPosYUpS,0] - p['ygrid_file'][0,0]) / ResolutionVeg
        
        # Saving the results
        #-------------------
        p['veg_posx'] = VegPosX
        p['veg_posy'] = VegPosY
    
    return s


    

def VegFlexure(s, p):
    if p['ny'] == 0:
        print('VegUpdate1D')    
        if p['vegnumb_file'].any():
            print('vegnum_file is here')
            print("wind = ",s['uw'][0,0])
        
            
    #     # Check if vegetation positions are in aeolis.txt or not
    #     if np.isscalar(p['veg_posx']) and p['veg_posx'] == 0:
    #         s = VegPositionWrite(s,p)  
        # print("wind = ",s['uw'][0,0])
        # Adjust leaf height dependig on wind speed
        # ------------------------------------------
            if s['uw'][0,0] ==0:
                p['leaf_height_adjusted'] = p['leaf_height']
            else:
                Flexure = p['flexure_a'] * s['uw'][0,0]+ p['flexure_b'] * p['veg_density'] + p['flexure_c']
                p['leaf_height_adjusted'] = p['leaf_height'] * (Flexure/100)
        
            # Check if the ajusted height is higher to initial height or inf to 0
            # --------------------------------------------------------------------
            if p['leaf_height_adjusted'] > p['leaf_height']:
                p['leaf_height_adjusted'] = p['leaf_height']
            elif p['leaf_height_adjusted'] < 0:
                p['leaf_height_adjusted'] = 0;
        
            # Computation of the surface according to the method
            # ---------------------------------------------------
            VegSurface = ((p['stem_diam'] * p['stem_height']) + (
                        (p['leaf_width'] * p['leaf_height_adjusted']) * p['leafToStem'])) * 1000000  # mm²
            # VegSurfaceRayon = np.sqrt(VegSurface / np.pi) # mm
        
            # Computation of the DownScaling factor
            # -------------------------------------
            ResolutionVeg = 0.001
            Resolution = np.around(p['xgrid_file'][1] - p['xgrid_file'][0], decimals=2)
            DownScale = np.round(Resolution / ResolutionVeg)
            nbCell = DownScale * DownScale
           
            # Initialisation of the VegGridCount
            #-----------------------------------
            veg1D = np.zeros(p['xgrid_file'].shape)
            veg1D = ((p['vegnumb_file']*VegSurface)/nbCell)/100
            veg1D[veg1D > 1] = 1
            
            s['rhovegadjusted'] = veg1D
            
        
            if np.isnan(s['rhoveg'][0, 0]):
                s['rhoveg'][:,:] = 0.
        
            #ix = s['rhoveg'] < 0.03
            #s['rhoveg'][ix] *= 0.
            
            s['hvegadjusted'][:,:] = p['hveg_max']*np.sqrt(s['rhovegadjusted'])
        
            s['germinate'][:,:] = (s['rhovegadjusted']>0)
            s['lateral'][:,:] = 0.
            
            return s
     
        else:
            print('vegnum_file 1D not here')
    
    # # Conversion of VegPos depending on the resolution
    #     #------------------------------------------------- 
    #     VegPosXUpS = np.floor((p['xgrid_file'][0] + (p['veg_posx'] * ResolutionVeg)) / Resolution) * Resolution
    #     #VegPosYUpS = np.floor((p['ygrid_file'][0,0] + (p['veg_posy'] * ResolutionVeg)) / Resolution) * Resolution
        
    #     # Initialisation of the VegGridCount
    #     #-----------------------------------
    #     veg1D = np.zeros(p['xgrid_file'].shape)
        
    #     # Iterate over plant positions and update the grid
    #     #-------------------------------------------------
    #     for i in range(len(VegPosXUpS)):
    #         # Calculate the cell index corresponding to the plant position
    #         indX = int((VegPosXUpS[i] - p['xgrid_file'][0]) / Resolution)
    #         #indY = int((VegPosYUpS[i] - p['ygrid_file'][0,0]) / Resolution)
    #         # Update the corresponding cell value
    #         veg1D[indX] += 1
        
    #     veg1D = (veg1D*VegSurface)/nbCell
    #     veg1D[veg1D > 1] = 1
      
    #     s['rhovegadjusted'][:, :] = veg1D
        
    
    #     if np.isnan(s['rhoveg'][0, 0]):
    #         s['rhoveg'][:,:] = 0.
    
    #     #ix = s['rhoveg'] < 0.03
    #     #s['rhoveg'][ix] *= 0.
        
    #     s['hvegadjusted'][:,:] = p['hveg_max']*np.sqrt(s['rhovegadjusted'])
    
    #     s['germinate'][:,:] = (s['rhovegadjusted']>0)
    #     s['lateral'][:,:] = 0.
        
    #     return s
    
    elif p['ny'] > 0:
        print('VegUpdate2D')
        # Check if vegetation numbers are in aeolis.txt or not
        if p['vegnumb_file'].any():
            print('vegnum_file is here')
            print("wind = ",s['uw'][0,0])
            # Adjust leaf height dependig on wind speed
            # ------------------------------------------
            if s['uw'][0,0] == 0:
                p['leaf_height_adjusted'] = p['leaf_height']
            else:
                Flexure = p['flexure_a'] * s['uw'][0,0]+ p['flexure_b'] * p['veg_density'] + p['flexure_c']
                p['leaf_height_adjusted'] = p['leaf_height'] * (Flexure/100)
        
            # Check if the ajusted height is higher to initial height or inf to 0
            # --------------------------------------------------------------------
            if p['leaf_height_adjusted'] > p['leaf_height']:
                p['leaf_height_adjusted'] = p['leaf_height']
            elif p['leaf_height_adjusted'] < 0:
                p['leaf_height_adjusted'] = 0;
        
            # Computation of the surface according to the method
            # ---------------------------------------------------
            VegSurface = ((p['stem_diam'] * p['stem_height']) + (
                        (p['leaf_width'] * p['leaf_height_adjusted']) * p['leafToStem'])) * 1000000  # mm²
            
            # Computation of the DownScaling factor
            # -------------------------------------
            ResolutionVeg = 0.001
            Resolution = np.around(p['xgrid_file'][0, 1] - p['xgrid_file'][0, 0], decimals=1)
            DownScale = np.round(Resolution / ResolutionVeg)
            nbCell = DownScale * DownScale
            
            # Initialisation of the VegGridCount
            #-----------------------------------
            veg2D = np.zeros(p['xgrid_file'].shape)
            veg2D = ((p['vegnumb_file']*VegSurface)/nbCell)/100
            veg2D[veg2D > 1] = 1
            
            s['rhovegadjusted'][:, :] = veg2D
            
        
            if np.isnan(s['rhovegadjusted'][0, 0]):
                s['rhovegadjusted'][:,:] = 0.
        
            #ix = s['rhoveg'] < 0.03
            #s['rhoveg'][ix] *= 0.
            
            s['hvegadjusted'][:,:] = p['hveg_max']*np.sqrt(s['rhovegadjusted'])
        
            s['germinate'][:,:] = (s['rhovegadjusted']>0)
            s['lateral'][:,:] = 0.
            
            return s
                            
        else:
            print('vegnum_file 2D not here')
                       
            
            # # Check if vegetation positions are in aeolis.txt or not
            # if (p['veg_posx'] == 0).any() and (p['veg_posy'] == 0).any():
            #     s = VegPositionWrite(s,p)
            # print("wind = ",s['uw'][0,0])
            # # Adjust leaf height dependig on wind speed
            # # ------------------------------------------
            # if s['uw'][0,0] == 0:
            #     p['leaf_height_adjusted'] = p['leaf_height']
            # else:
            #     Flexure = p['flexure_a'] * s['uw'][0,0]+ p['flexure_b'] * p['veg_density'] + p['flexure_c']
            #     p['leaf_height_adjusted'] = p['leaf_height'] * (Flexure/100)
        
            # # Check if the ajusted height is higher to initial height or inf to 0
            # # --------------------------------------------------------------------
            # if p['leaf_height_adjusted'] > p['leaf_height']:
            #     p['leaf_height_adjusted'] = p['leaf_height']
            # elif p['leaf_height_adjusted'] < 0:
            #     p['leaf_height_adjusted'] = 0;
        
            # # Computation of the surface according to the method
            # # ---------------------------------------------------
            # VegSurface = ((p['stem_diam'] * p['stem_height']) + (
            #             (p['leaf_width'] * p['leaf_height_adjusted']) * p['leafToStem'])) * 1000000  # mm²
            # # VegSurfaceRayon = np.sqrt(VegSurface / np.pi) # mm
        
            # # Computation of the DownScaling factor
            # # -------------------------------------
            # ResolutionVeg = 0.001
            # Resolution = np.around(p['xgrid_file'][0, 1] - p['xgrid_file'][0, 0], decimals=1)
            # DownScale = np.round(Resolution / ResolutionVeg)
            # nbCell = DownScale * DownScale
        
            # # Conversion of VegPos depending on the resolution
            # #------------------------------------------------- 
            # VegPosXUpS = np.floor((p['xgrid_file'][0,0] + (p['veg_posx'] * ResolutionVeg)) / Resolution) * Resolution
            # VegPosYUpS = np.floor((p['ygrid_file'][0,0] + (p['veg_posy'] * ResolutionVeg)) / Resolution) * Resolution
            
            # # Initialisation of the VegGridCount
            # #-----------------------------------
            # veg2D = np.zeros(p['xgrid_file'].shape)
            
            # # Iterate over plant positions and update the grid
            # #-------------------------------------------------
            # for i in range(len(VegPosXUpS)):
            #     # Calculate the cell index corresponding to the plant position
            #     indX = int((VegPosXUpS[i] - p['xgrid_file'][0,0]) / Resolution)
            #     indY = int((VegPosYUpS[i] - p['ygrid_file'][0,0]) / Resolution)
            #     # Update the corresponding cell value
            #     veg2D[indY, indX] += 1
            
            # veg2D = (veg2D*VegSurface)/nbCell
            # veg2D[veg2D > 1] = 1
            
            # s['rhovegadjusted'][:, :] = veg2D
            
        
            # if np.isnan(s['rhovegadjusted'][0, 0]):
            #     s['rhovegadjusted'][:,:] = 0.
        
            # #ix = s['rhoveg'] < 0.03
            # #s['rhoveg'][ix] *= 0.
            
            # s['hvegadjusted'][:,:] = p['hveg_max']*np.sqrt(s['rhovegadjusted'])
        
            # s['germinate'][:,:] = (s['rhovegadjusted']>0)
            # s['lateral'][:,:] = 0.
            
            return s



def vegshear(s, p):
    if p['vegshear_type'] == 'okin':
        s = vegshear_okin(s, p)
    else:
        s = vegshear_raupach(s, p)

    s = velocity_stress(s,p)

    return s

def adjustment_length(s, p):
    #print('VegLag')
    if p['process_veglag']:
        #print('VegLag in if')
        RhoVeg = p['veg_density'] 
        if p['process_vegflex']: 
            s = VegFlexure(s, p)
            CanopyHeight = p['leaf_height_adjusted'] + p['stem_height']
            zp = s['hvegadjusted']
        else:
            CanopyHeight = p['h_canopy']
            zp = s['hveg']
        StretchHeight = p['h_stretch']
        StemDiam = p['stem_diam']
        LeafWidth = p['leaf_width']
        LeafThickness = p['leaf_thick']
        alpha = p['leafToStem']
        u = s['uw']  # chagne to 'uw' np.asarray([14.])
        theta = s['udir'][0, 0]
        theta = np.radians(theta)

        # Intialize other grid parameters
        x = s['x'][0, :]
        dx = s['x'][0, 1] - s['x'][0, 0]

        # Wind parameters
        u_10 = u[0, 0]
        zref = p['z']
        kappa = p['kappa']
        d = p['grain_size']
        nu = p['v']
        z0 = p['k']  # 2 * d[0]/30.
        z0 = float(z0)
        # find windspeed at canopy
        ustar = u_10 * kappa / np.log(zref / z0)

        # Calcualte canopy length
        posVeg = np.where(zp != 0)
        if len(posVeg[0]) > 0:
            posVegStart = np.min(posVeg[1])
            posVegEnd = np.max(posVeg[1])
            PatchLength = (posVegEnd - posVegStart) * dx
            PatchLength = np.abs(PatchLength)
        else:
            PatchLength = 0

        # Approximate solid volume fraciton and lateral cover
        SVF = alpha * RhoVeg * LeafThickness * LeafWidth * (StretchHeight / CanopyHeight)

        # Lateral cover
        LC = alpha * RhoVeg * LeafWidth * StretchHeight  # lateral cover approximation

        # Calculate element drag coefficent
        U_canopy = (ustar / kappa) * (np.log(CanopyHeight / z0))
        Reynolds = U_canopy * StemDiam / nu
        Cd = 1 + 10 * Reynolds ** (-2 / 3)
        CanopyDragLength = (2 * (1 - SVF) * CanopyHeight) / (Cd * LC)

        # Use empirical relation to determine deposition lag
        lag = .95 * CanopyDragLength - 0.09  # with HespData
        # lag = 4.48 * CanopyDragLength  - 1.29 # Hesp + John data

        # Correct lag is out of range
        # if lag > PatchLength:
        #     lag = PatchLength
        # elif lag < 0:
        #     lag = 0

    else:
        lag = 0
    print("lag = ",lag)
    return lag

def germinate(s,p):
    ny = p['ny']
    
    # time [year]
    n = (365.25*24.*3600. / (p['dt_opt'] * p['accfac']))
    
    # Determine which cells are already germinated before
    s['germinate'][:, :] = (s['rhoveg'] > 0.)
    
    # Germination
    p_germinate_year = p['germinate']                                
    p_germinate_dt = 1-(1-p_germinate_year)**(1./n)
    germination = np.random.random((s['germinate'].shape))
    
    # Germinate new cells
    germinate_new = (s['dzbveg'] >= 0.) * (germination <= p_germinate_dt)
    s['germinate'] += germinate_new.astype(float)
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

    # Adding growth
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
    s['germinate'] *= (s['rhoveg']!=0.)
    s['lateral'] *= (s['rhoveg']!=0.)

    # Dying of vegetation due to hydrodynamics (Dynamic Vegetation Limit)
    if p['process_tide']:
        s['rhoveg']     *= (s['zb'] +0.01 >= s['zs'])
        s['hveg']       *= (s['zb'] +0.01 >= s['zs'])
        s['germinate']  *= (s['zb'] +0.01 >= s['zs'])
        s['lateral']    *= (s['zb'] +0.01 >= s['zs'])

    ix = s['zb'] < p['veg_min_elevation']
    s['rhoveg'][ix] = 0
    s['hveg'][ix] = 0
    s['germinate'][ix] = 0
    s['lateral'][ix] = 0

    return s

def vegshear_okin(s, p):
    # Okin 1D
    if p['ny'] == 0:
        print('Okin1D')
        # Approach to calculate shear reduction in the lee of plants using the general approach of:
        # Okin (2008), JGR, A new model of wind erosion in the presence of vegetation
        # Note that implementation only works in 1D currently
        # print('Okin1D')
        # Initialize shear variables and other grid parameters
        ustar = s['ustar'].copy()
        ustars = s['ustars'].copy()
        ustarn = s['ustarn'].copy()
        ets = np.zeros(s['zb'].shape)
        etn = np.zeros(s['zb'].shape)
        ix = ustar != 0
        ets[ix] = ustars[ix] / ustar[ix]
        etn[ix] = ustarn[ix] / ustar[ix]
        udir = s['udir'][0, 0] + 180
    
        # intialize other grid parameters
        x = s['x'][0, :]
        dx = s['x'][0, 1] - s['x'][0, 0]
        if p['process_vegflex']:
            s = VegFlexure(s, p)
            zp = s['hvegadjusted'][0,:]
        else:
            zp = s['hveg'][0, :]
    
        red = np.zeros(x.shape)
        red_all = np.zeros(x.shape)
        nx = x.size
    
        # okin model defaults - hardcoded for now
        c1 = p['okin_c1_veg']
        if p['process_okinR0var']:
            print('OkinR0var')
            intercept = p['okin_initialred_vegLow'] + (p['okin_initialred_vegHigh']-p['okin_initialred_vegLow']) * ((s['uw'][0,0] - 5.97) / (9.58 - 5.97))
        else:            
            intercept = p['okin_initialred_veg']
        
        print("R0 = ",intercept)
    
        # adjustment length modification variables
        lag = adjustment_length(s, p)
    
        # Find index of vegetation
        zplag = zp.copy()
        posVeg = np.where(zplag != 0)[0]
    
        # If vegetation is present in zplag -> apply the lag
        if len(posVeg) > 0:
            posVegStart = posVeg[0]
            posVegEnd = posVeg[-1]
    
            # Orientation of the lag length depending on wind direction
            if udir >= 180 and udir <= 360:
                zplag[int(posVegEnd) - int(lag / dx): int(posVegEnd)] = 0
            else:
                zplag[int(posVegStart): int(posVegStart) + int(lag / dx)] = 0
    
        for igrid in range(nx):
    
            # only look at cells with a roughness element
            if zplag[igrid] > 0:
                # local parameters
                mult = np.ones(x.shape)
                h = zplag[igrid]
    
                if udir >= 180 and udir <= 360:
                    xrel = -(x - x[igrid])
                    # print('loop1')
                    # for igrid2 in range(nx-1):
                    # if edge[igrid2] == -1 and edge[igrid2-1] == 0:
                    #  deadzone[igrid2 - int(lag/dx) + 1:igrid2 + 1:1] = 1
                else:
                    xrel = x - x[igrid]
                    # for igrid2 in range(nx - 1):
                    # if edge[igrid2] == 1 and edge[igrid2 + 1] == 0:
                    #  realedge = igrid2 + 1
                    #  deadzone[realedge:realedge + int(lag/dx):1] = 1
    
                for igrid2 in range(nx):
    
                    if xrel[igrid2] >= 0 and xrel[igrid2] / h < 20:  # and deadzone[igrid2] !=1:
    
                        # apply okin model
                        mult[igrid2] = intercept + (1 - intercept) * (1 - math.exp(-xrel[igrid2] * c1 / h))
                        # print('Okin runing')
                    # else:
                    # print('No okin')
    
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
    
        # convert to a multiple
        mult_all = 1 - red_all
        ix = mult_all < 0.001
        mult_all[ix] = 0.001
    
        s['ustar'][0, :] = s['ustar'][0, :] * mult_all
        s['ustars'][0, :] = s['ustar'][0, :] * ets[0, :]
        s['ustarn'][0, :] = s['ustar'][0, :] * etn[0, :]
        
        return s
    
    elif p['ny'] > 1:
        print('Okin2D')
        # Okin 2D
        # Approach to calculate shear reduction in the lee of plants using the general approach of:
        # Okin (2008), JGR, A new model of wind erosion in the presence of vegetation
        # print('Okin2D')
        # Initialize shear variables and other grid parameters
        ustar = s['ustar'].copy()
        ustars = s['ustars'].copy()
        ustarn = s['ustarn'].copy()
        ets = np.zeros(s['zb'].shape)
        etn = np.zeros(s['zb'].shape)
        ix = ustar != 0
        ets[ix] = ustars[ix] / ustar[ix]
        etn[ix] = ustarn[ix] / ustar[ix]
        udir = s['udir'][0, 0] + 180

        # Intialize other grid parameters
        x = s['x'][:, :]
        dx = s['x'][0, 1] - s['x'][0, 0]
        nx = x[1, :].size
        y = s['y'][:, :]
        dy = s['y'][1, 0] - s['y'][0, 0]
        ny = y[:, 1].size

        if p['process_vegflex']:
            s = VegFlexure(s, p)
            zp = s['hvegadjusted'][:,:]
        else:
            zp = s['hveg'][:, :]
        red = np.zeros(x.shape)
        red_all = np.zeros(x.shape)
        mult_all = np.zeros(x.shape)

        # Okin model defaults
        c1 = p['okin_c1_veg']
        if p['process_okinR0var']:
            print('OkinR0var')
            intercept = p['okin_initialred_vegLow'] + (p['okin_initialred_vegHigh']-p['okin_initialred_vegLow']) * ((s['uw'][0,0] - 5.97) / (9.58 - 5.97))
        else:            
            intercept = p['okin_initialred_veg']
        
        print("R0 = ",intercept)
        # Adjustment length modification variables
        lag = adjustment_length(s, p)

        if udir < 360:
            udir = udir + 360

        if udir > 360:
            udir = udir - 360
        # START OF THE OKIN MODEL LOOP
        # Calculate shear reduction by looking through all cells that have plants present
        # and looking downwind of those features
        # For each cross-shore transect
        for jgrid in range(ny):
            # Create zplag where the lag max will be applied
            zplag = zp.copy()[jgrid]
            # Find index of vegetation
            posVeg = np.where(zplag != 0)[0]

            # If vegetation is present in zplag -> apply the lag
            if len(posVeg) > 0:
                posVegStart = posVeg[0]
                posVegEnd = posVeg[-1]
                # Orientation of the lag length depending on wind direction
                if udir >= 180 and udir <= 360:
                    zplag[int(posVegEnd) - int(lag / dx): int(posVegEnd)] = 0
                else:
                    zplag[int(posVegStart): int(posVegStart) + int(lag / dx)] = 0

                # For each cell of one cross-hore transect
                for igrid in range(nx):

                    if zplag[igrid] > 0:
                        # only look at cells with a roughness element
                        mult = np.ones(nx)
                        h = zplag[igrid]  # vegetation height at the appropriate cell

                        if udir >= 180 and udir <= 360:
                            xrel = -(x[jgrid, :] - x[jgrid, igrid])
                        else:
                            xrel = x[jgrid, :] - x[jgrid, igrid]

                        for igrid2 in range(nx):

                            if xrel[igrid2] >= 0 and xrel[igrid2] / h < 20:
                                # apply okin model
                                mult[igrid2] = intercept + (1 - intercept) * (1 - math.exp(-xrel[igrid2] * c1 / h))

                        red[jgrid, :] = 1 - mult

                        # fix potential issues for summation
                        ix = red[jgrid, :] < 0.00001
                        red[jgrid, ix] = 0
                        ix = red[jgrid, :] > 1
                        red[jgrid, ix] = 1
                        ix = xrel < 0
                        red[jgrid, ix] = 0

                        # combine all reductions between plants
                        red_all[jgrid, :] = red_all[jgrid, :] + red[jgrid, :]

                # cant have more than 100% reduction
                ix = red_all[jgrid, :] > 1
                red_all[jgrid, ix] = 1

                # update shear velocity according to Okin (note does not operate on shear stress)
                mult_all[jgrid, :] = 1 - red_all[jgrid, :]
                ustarveg = s['ustar'][jgrid, :] * mult_all[jgrid, :]
                ix = ustarveg < 0.01
                ustarveg[ix] = 0.01  # some small number so transport code doesnt crash

                s['ustar'][jgrid, :] = ustarveg
                s['ustars'][jgrid, :] = s['ustar'][jgrid, :] * ets[jgrid, :]
                s['ustarn'][jgrid, :] = s['ustar'][jgrid, :] * etn[jgrid, :]

        return s

# def vegshear_okin2d(s, p):
#     # Approach to calculate shear reduction in the lee of plants using the general approach of:
#     # Okin (2008), JGR, A new model of wind erosion in the presence of vegetation
#     # print('Okin2D')
#     # Initialize shear variables and other grid parameters
#     ustar = s['ustar'].copy()
#     ustars = s['ustars'].copy()
#     ustarn = s['ustarn'].copy()
#     ets = np.zeros(s['zb'].shape)
#     etn = np.zeros(s['zb'].shape)
#     ix = ustar != 0
#     ets[ix] = ustars[ix] / ustar[ix]
#     etn[ix] = ustarn[ix] / ustar[ix]
#     udir = s['udir'][0, 0] + 180

#     # Intialize other grid parameters
#     x = s['x'][:, :]
#     dx = s['x'][0, 1] - s['x'][0, 0]
#     nx = x[1, :].size
#     y = s['y'][:, :]
#     dy = s['y'][1, 0] - s['y'][0, 0]
#     ny = y[:, 1].size

#     if p['process_vegflex']:
#         s,LeafHeightAdjusted = VegFlexure(s, p)
#         zp = s['hvegadjusted'][:,:]
#     else:
#         zp = s['hveg'][:, :]
#     red = np.zeros(x.shape)
#     red_all = np.zeros(x.shape)
#     mult_all = np.zeros(x.shape)

#     # Okin model defaults
#     c1 = p['okin_c1_veg']
#     intercept = p['okin_initialred_veg']

#     # Adjustment length modification variables
#     lag = adjustment_length(s, p, LeafHeightAdjusted)

#     if udir < 360:
#         udir = udir + 360

#     if udir > 360:
#         udir = udir - 360
#     # START OF THE OKIN MODEL LOOP
#     # Calculate shear reduction by looking through all cells that have plants present
#     # and looking downwind of those features
#     # For each cross-shore transect
#     for jgrid in range(ny):
#         # Create zplag where the lag max will be applied
#         zplag = zp.copy()[jgrid]
#         # Find index of vegetation
#         posVeg = np.where(zplag != 0)[0]

#         # If vegetation is present in zplag -> apply the lag
#         if len(posVeg) > 0:
#             posVegStart = posVeg[0]
#             posVegEnd = posVeg[-1]
#             # Orientation of the lag length depending on wind direction
#             if udir >= 180 and udir <= 360:
#                 zplag[int(posVegEnd) - int(lag / dx): int(posVegEnd)] = 0
#             else:
#                 zplag[int(posVegStart): int(posVegStart) + int(lag / dx)] = 0

#             # For each cell of one cross-hore transect
#             for igrid in range(nx):

#                 if zplag[igrid] > 0:
#                     # only look at cells with a roughness element
#                     mult = np.ones(nx)
#                     h = zplag[igrid]  # vegetation height at the appropriate cell

#                     if udir >= 180 and udir <= 360:
#                         xrel = -(x[jgrid, :] - x[jgrid, igrid])
#                     else:
#                         xrel = x[jgrid, :] - x[jgrid, igrid]

#                     for igrid2 in range(nx):

#                         if xrel[igrid2] >= 0 and xrel[igrid2] / h < 20:
#                             # apply okin model
#                             mult[igrid2] = intercept + (1 - intercept) * (1 - math.exp(-xrel[igrid2] * c1 / h))

#                     red[jgrid, :] = 1 - mult

#                     # fix potential issues for summation
#                     ix = red[jgrid, :] < 0.00001
#                     red[jgrid, ix] = 0
#                     ix = red[jgrid, :] > 1
#                     red[jgrid, ix] = 1
#                     ix = xrel < 0
#                     red[jgrid, ix] = 0

#                     # combine all reductions between plants
#                     red_all[jgrid, :] = red_all[jgrid, :] + red[jgrid, :]

#             # cant have more than 100% reduction
#             ix = red_all[jgrid, :] > 1
#             red_all[jgrid, ix] = 1

#             # update shear velocity according to Okin (note does not operate on shear stress)
#             mult_all[jgrid, :] = 1 - red_all[jgrid, :]
#             ustarveg = s['ustar'][jgrid, :] * mult_all[jgrid, :]
#             ix = ustarveg < 0.01
#             ustarveg[ix] = 0.01  # some small number so transport code doesnt crash

#             s['ustar'][jgrid, :] = ustarveg
#             s['ustars'][jgrid, :] = s['ustar'][jgrid, :] * ets[jgrid, :]
#             s['ustarn'][jgrid, :] = s['ustar'][jgrid, :] * etn[jgrid, :]

#     return s
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

