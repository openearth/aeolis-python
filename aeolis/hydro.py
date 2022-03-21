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
#from numba import jit
import logging
import numpy as np
import sympy as sym

# package modules
from aeolis.utils import *


# initialize logger
logger = logging.getLogger(__name__)


def interpolate(s, p, t):
    '''Interpolate hydrodynamic and meteorological conditions to current time step

    Interpolates the hydrodynamic and meteorological time series to
    the current time step, if available. Meteorological parameters are
    stored as dictionary rather than a single value.

    Parameters
    ----------
    s : dict
        Spatial grids
    p : dict
        Model configuration parameters
    t : float
        Current time

    Returns
    -------
    dict
        Spatial grids

    '''

    if p['process_tide']:
        if p['tide_file'] is not None:
            s['zs'][:,:] = interp_circular(t,
                                           p['tide_file'][:,0],
                                           p['tide_file'][:,1])
            s['swl'][:,:] = interp_circular(t,
                                           p['tide_file'][:,0],
                                           p['tide_file'][:,1])
        else:
            s['zs'][:,:] = 0.
            s['swl'][:,:] = 0.

        # apply complex mask
        s['zs'] = apply_mask(s['zs'], s['tide_mask'])
        
    if p['process_wave'] and p['wave_file'] is not None:

        # determine water depth
        h = np.maximum(0., s['zs'] - s['zb'])
    
        s['Hs'][:,:] = interp_circular(t,
                                       p['wave_file'][:,0],
                                       p['wave_file'][:,1])

        # apply complex mask
        s['Hs'] = apply_mask(s['Hs'], s['wave_mask'])

        # add wave runup
        if p['process_runup']:
            ix = s['Hs'] > 0.
            R = p['xi'] * s['Hs']
            s['zs'][ix] += R[ix] * (1. - np.minimum(1., h[ix] * p['gamma'] / s['Hs'][ix]))
        
    if p['process_moist'] and p['method_moist_process'].lower() == 'surf_moisture' and p['meteo_file'] is not None:  ## här lägg till process

        m = interp_array(t,
                         p['meteo_file'][:,0],
                         p['meteo_file'][:,1:], circular=True)

        #Meteorological parameters (Symbols according to KNMI, units according to the Penman's equation)
        # T: Temperature, Degrees Celsius
        # Q : Global radiation, MJ/m2/d
        # RH : Precipitation, mm/h
        # P : Atmospheric pressure, kPa
        # U: Relative humidity, %
        
        s['meteo'] = dict(zip(('T','Q','RH','P','U') , m))
        

    # ensure compatibility with XBeach: zs >= zb
#    s['zs'] = np.maximum(s['zs'], s['zb'])

    return s


def update(s, p, dt, t):
    '''Update soil moisture content

    Updates soil moisture content in all cells. The soil moisure
    content in flooded cells is set to the porosity. The soil moisture
    content in non-flooded cells is decreased by simulating
    infiltration using an exponential decay function with a rate ``F``
    following:

    .. math::

       M = M \\cdot e^{-F \\cdot \\Delta t}

    Cells are considered flooded if the water depth is larger than
    ``eps``. Aditionally, is meteorological conditions are provided
    the soil moisture content is decreased by simulating evaporation
    following the Penman equation:

    .. math::

        M = M - \\frac{m \\cdot R + \\gamma \cdot 6.43 \cdot (1 + 0.536 \cdot u]) \cdot \\delta}{\\lambda \cdot (m + \\gamma)}

    Parameters
    ----------
    s : dict
        Spatial grids
    p : dict
        Model configuration parameters
    dt : float
        Current time step

    Returns
    -------
    dict
        Spatial grids

    '''
    # Groundwater level Boussinesq (1D CS-transects)
    if p['process_groundwater']:
        

        #Specify wetting or drying conditions in previous timestep
        s['wetting'] = s['gw'] > s['gw_prev']
        #Save groundwater from previous timestep
        s['gw_prev'] = s['gw']        
        
        #Decrease timestep for GW computations
        dt_gw = int(dt / p['tfac_gw'])
        for i in range(int(dt / dt_gw)):
            t_gw = t + i * dt_gw
            interpolate(s,p,t_gw)
        
            #Compute setup and runup
            setup=np.max(s['swl'] + 0.35 * p['xi'] * s['Hs']) #Stockdon et al (2006)
            
            runup=np.max(s['zs'])

            #Initialize GW levels
            if t==0:
                s['gw'][:,:]=p['in_gw']
                s['gw_prev'] = s['gw']
                s['wetting'] = s['gw'] > s['gw_prev']

                
            #Define index of shoreline location
            shl_ix =np.argmax(s['zb'] > setup,axis=1) - 1
            
            #Define index of runup limit
            runup_ix =np.argmax(s['zb'] > runup,axis=1) - 1
        
            #Runge-Kutta timestepping
            f1 = Boussinesq(s['gw'],s,p,setup, shl_ix)
            f2 = Boussinesq(s['gw'] + dt_gw / 2 * f1,s,p,setup, shl_ix)
            f3 = Boussinesq(s['gw'] + dt_gw / 2 * f2,s,p,setup, shl_ix)
            f4 = Boussinesq(s['gw'] + dt_gw * f3,s,p,setup, shl_ix)
            
            #Update groundwater level
            s['gw'] = s['gw'] + dt_gw / 6 * (f1 + 2 * f2 + 2 * f3 + f4)
        
            #Add infiltration from wave runup according to Nielsen (1990)
            if p['process_wave']:              
                #Compute f(x) = distribution of infiltrated water
                fx=np.zeros(s['gw'].shape)
                fx_ix=np.zeros_like(shl_ix)
                for i in range(len(s['gw'][:,0])):
                    #Define index of peak f(x)
                    fx_ix[i] = (shl_ix[i]) + (2/3 * (runup_ix[i] - shl_ix[i]))
                    #Compute f(X)
                    fx[i,shl_ix[i]:fx_ix[i]] = (s['x'][i,shl_ix[i]:fx_ix[i]] - s['x'][i,shl_ix[i]]) / (2 / 3 * (s['x'][i,runup_ix[i]] - s['x'][i,shl_ix[i]]))
                    fx[i,fx_ix[i]+1:runup_ix[i]] = 3 - (s['x'][i,fx_ix[i]+1:runup_ix[i]]- s['x'][i,shl_ix[i]])  / (1 / 3 * (s['x'][i,runup_ix[i]] - s['x'][i,shl_ix[i]]))
                    fx[i,fx_ix[i]]=1
                    
                # Update groundwater level with overheight due to runup
                s['gw'] = s['gw'] + p['Cl_gw'] * fx * p['K_gw'] / p['ne_gw']
            


            # #Set GW level to the setup level in submerged cells
            # for i in range(len(s['gw'][:,0])):
            #     for j in range(shl_ix[i]):                        
            #         s['gw'][i,j] = setup

            # Apply GW complex mask
            s['gw'] = apply_mask(s['gw'], s['gw_mask'])
            
                
            # Do not allow GW levels above ground level
            s['gw']=np.minimum(s['gw'], s['zb'])

            # Define cells below setup level
            ixg=s['zb'] < setup
            
            # Set gw level to setup level in cells below setup level
            
            s['gw'][ixg]=setup            



    # Compute surface moisture with infiltration method using Darcy
    if p['process_moist']:
        if p['method_moist_process'].lower() == 'infiltration':
            F1 = -np.log(.5) / p['Tdry']
            ix = s['zs'] - s['zb'] > p['eps']
            s['moist'][ ix] = p['porosity']
            s['moist'][~ix] *= np.exp(-F1 * dt)
            s['moist'][:,:] = np.maximum(0.,np.minimum(p['porosity'],\
                                                  s['moist'][:,:]))
        
        # Compute surface moisture accounting for runup, capillary rise and precipitation/evaporation
        elif p['method_moist_process'].lower() == 'surf_moisture':
            if p['process_groundwater'] is None :
                logger.log_and_raise('process_groundwater is not activated, the groundwater level is not computed within the program but set constant at 0 m', exc=ValueError)
                
            #Infiltration
            F1 = -np.log(.5) / p['Tdry']
            s['moist'] += np.minimum(0,(s['moist']-p['fc'])*(np.exp(-F1*dt)-1))
            
            #If the cell is flooded (runup) in this timestep, assume satiation
            ix = s['zs'] > s['zb']
            s['moist'][ix] = p['satd_moist']

            
            #Update surface moisture with respect to evaporation, condensation, and precipitation
            met = s['meteo']
            evo = evaporation(s,p,met)
            evo = evo / 24. / 3600. / 1000. # convert evaporation from mm/day to m/s
            pcp = met['RH'] / 3600. / 1000. # convert precipitation from mm/hr to m/s
            s['moist'][~ix] = np.maximum(s['moist'][~ix] + (pcp - evo[~ix]) * dt / p['thick_moist'], p['resw_moist'])
            s['moist'][~ix] = np.minimum(s['moist'][~ix],p['satd_moist'])
                        
            #Compute surface moisture due to capillary processes (van Genuchten and Mualem II)
            
            #Compute distance from gw table to the soil surface
            h=np.maximum(0,(s['zb'] - s['gw']) * 100) #h in cm to match convention of alfa (cm-1)
            
            if p['process_scanning']:
            #Initialize value of surface moisture due to capillary rise
                if t == 0:
                    s['moist_swr'] = p['resw_moist'] + (p['satw_moist'] - p['resw_moist']) \
                          / (1 + abs(p['alfaw_moist'] * h) ** p['nw_moist']) ** p['mw_moist']
                    s['h_delta'][:,:]=np.maximum(0,(s['zb'] - s['gw'])*100)
                    s['scan_w'][:,:] == False
                    s['scan_d'][:,:] == False
                else:
                    #Compute h_delta
                    for i in range(len(s['wetting'][:,0])):
                        for j in range(len(s['wetting'][0,:])):
                            #Compute h delta on the main drying and wetting curve
                            if s['scan_w'][i,j] == False and s['wetting'][i,j] == True and s['gw'][i,j] < s['gw_prev'][i,j] or s['scan_d'][i,j] == False and s['wetting'][i,j] == False and s['gw'][i,j] > s['gw_prev'][i,j]:
                                s['h_delta'][i,j]=np.maximum(0,(s['zb'][i,j] - s['gw'][i,j])*100)
                            #Compute h_delta if there is a reversal on the wetting scanning curve
                            if s['scan_w'][i,j] == True and s['wetting'][i,j] == True and s['gw'][i,j] < s['gw_prev'][i,j]:
                                #Solve hdelta from drying scanning curve for which moist(h) on drying scanning curve equals moist(h) on wetting scanning curve 
                                #intermediate solution:
                                w_hdelta_int = np.minimum((s['scan_w_moist'][i,j] - s['w_h'][i,j]) * (p['satd_moist'] - s['w_h'][i,j]) / (s['d_h'][i,j] - s['w_h'][i,j]) + s['w_h'][i,j], p['satw_moist'])
                                #Solve hdelta from wetting curve
                                s['h_delta'][i,j] = np.maximum( 1 / p['alfaw_moist'] * (((p['satw_moist'] - p['resw_moist']) \
                                  / np.maximum((w_hdelta_int - p['resw_moist']),0.00001)) ** (1 / p['mw_moist']) - 1) ** (1 / p['nw_moist']),0)
                            #Compute h_delta if there is a reversal on the drying scanning curve
                            if s['scan_d'][i,j] == True and s['wetting'][i,j] == False and s['gw'][i,j] > s['gw_prev'][i,j]:
                                #Solve hdelta from wetting scanning curve for which moist(h) on wetting scanning curve equals moist(h) on drying scanning curve
                               
                                #Simple iteration method
                                hdelta=0 #initialize hdelta
                                it_hdelta(hdelta,s,p,i,j)
                                s['h_delta'][i,j] = hdelta
                    
                    #Compute moisture of h for the wetting curve
                    s['w_h'] = p['resw_moist'] + (p['satw_moist'] - p['resw_moist']) \
                                      / (1 + abs(p['alfaw_moist'] * h) ** p['nw_moist']) ** p['mw_moist']
                    #Compute moisture of h_delta for the wetting curve
                    s['w_hdelta'] = p['resw_moist'] + (p['satw_moist'] - p['resw_moist']) \
                                      / (1 + abs(p['alfaw_moist'] * s['h_delta']) ** p['nw_moist']) ** p['mw_moist']
                    #Compute moisture of h for the drying curve
                    s['d_h'] = p['resd_moist'] + (p['satd_moist'] - p['resd_moist']) \
                                      / (1 + abs(p['alfad_moist'] * h) ** p['nd_moist']) ** p['md_moist']
                    #Compute moisture of h_delta for the drying curve
                    s['d_hdelta'] = p['resd_moist'] + (p['satd_moist'] - p['resd_moist']) \
                                      / (1 + abs(p['alfad_moist'] * s['h_delta']) ** p['nd_moist']) ** p['md_moist']
                    #Compute moisture content with the wetting scanning curve
                    s['scan_w_moist'] = np.maximum(np.minimum(s['w_h'] + (p['satw_moist'] - s['w_h']) / np.maximum(p['satw_moist'] - s['w_hdelta'],0.0001) \
                                                     * (s['d_hdelta'] - s['w_hdelta']),s['d_h']),s['w_h'])                
                    #Compute moisture content with the drying scanning curve
                    s['scan_d_moist'] = np.maximum(np.minimum(s['w_h'] + (s['w_hdelta'] - s['w_h']) / np.maximum(p['satd_moist'] - s['w_h'],0.0001) \
                                                      * (s['d_h'] - s['w_h']), s['d_h']),s['w_h'])
                    
                    #Select SWR curve to compute moisture content due to capillary processes
                    for i in range(len(s['wetting'][:,0])):
                        for j in range(len(s['wetting'][0,:])):
                            #Wetting conditions main curve
                            if s['gw'][i,j] >= s['gw_prev'][i,j] and s['wetting'][i,j] == True and s['scan_w'][i,j] == False:
                                s['moist_swr'][i,j]=s['w_h'][i,j]
                                s['scan_w'][i,j] = False
                                s['scan_d'][i,j] = False
                            #wetting conditions, timestep of reversal - move onto wetting scanning curve
                            elif s['gw'][i,j] >= s['gw_prev'][i,j] and s['wetting'][i,j] == False:
                                s['moist_swr'][i,j] = s['scan_w_moist'][i,j]
                                s['scan_w'][i,j] = s['scan_w_moist'] [i,j] > s['w_h'][i,j]
                                s['scan_d'][i,j] = False
                            #wetting conditions - followed a wetting scanning curve in previous timestep - continue following scanning curve unless main curve is reached
                            elif s['gw'][i,j] >= s['gw_prev'][i,j] and s['wetting'][i,j] == True and s['scan_w'][i,j] == True:
                                s['moist_swr'][i,j] = s['scan_w_moist'][i,j]
                                s['scan_w'][i,j] = s['scan_w_moist'] [i,j] > s['w_h'][i,j]
                                s['scan_d'][i,j] = False
                            #Drying conditions main curve
                            elif s['gw'][i,j] < s['gw_prev'][i,j] and s['wetting'][i,j] == False and s['scan_d'][i,j] == False:
                                s['moist_swr'][i,j]=s['d_h'][i,j]
                                s['scan_d'][i,j] = False
                                s['scan_w'][i,j] = False
                            #Drying conditions, timestep of reversal - move onto a drying scanning curve
                            elif s['gw'][i,j] < s['gw_prev'][i,j] and s['wetting'][i,j] == True:
                                s['moist_swr'][i,j] = s['scan_d_moist'][i,j]
                                s['scan_d'][i,j] = s['scan_d_moist'] [i,j] < s['d_h'][i,j]
                                s['scan_w'][i,j] = False
                            #Drying conditions - followed a drying scanning curve in previous timestep - continue following scanning curve unless main curve is reached
                            elif s['gw'][i,j] < s['gw_prev'][i,j] and s['wetting'][i,j] == False and s['scan_d'][i,j] == True:
                                s['moist_swr'][i,j] = s['scan_d_moist'][i,j]
                                s['scan_d'][i,j] = s['scan_d_moist'] [i,j] < s['d_h'][i,j]
                                s['scan_w'][i,j] = False
                
            else:
                ixw = s['wetting'] == True
                s['moist_swr'][ixw] = p['resw_moist'] + (p['satw_moist'] - p['resw_moist']) \
                                      / (1 + abs(p['alfaw_moist'] * h[ixw]) ** p['nw_moist']) ** p['mw_moist']
                s['moist_swr'][~ixw] = p['resd_moist'] + (p['satd_moist'] - p['resd_moist']) \
                                      / (1 + abs(p['alfad_moist'] * h[~ixw]) ** p['nd_moist']) ** p['md_moist']
                
            #Update surface moisture with respect to capillary processes
            s['moist'] = np.minimum(np.maximum(s['moist'],s['moist_swr']),p['satd_moist'])
            
            
            # if t > 0:
            #     q= 111
            #     w=1
            #     print(s['moist'][w,q],s['moist_swr'][p,q],s['wetting'][p,q],s['scan_d'][p,q],s['scan_w'][p,q],s['scan_d_moist'][p,q],s['scan_w_moist'][p,q], w_h[p,q],d_h[p,q],w_hdelta[p,q],d_hdelta[p,q],h[p,q],s['h_delta'][p,q])

        

        
        else:
            logger.log_and_raise('Unknown moisture process formulation [%s]' % p['method_moist_process'], exc=ValueError)
        

    # salinitation
    if p['process_salt']:
        met = s['meteo']
        F2 = -np.log(.5) / p['Tsalt']
        s['salt'][ ix,0] = 1.
        s['salt'][~ix,0] *= np.exp(-F2 * dt)
        pcp = met['RH'] / 3600. / 1000. # convert precipitation from mm/hr to m/s
        s['salt'][:,:,0] = np.minimum(1., s['salt'][:,:,0] + pcp * dt / p['layer_thickness'])

    return s




def Boussinesq (GW,s,p,setup,shl_ix):
    '''
    Add description
    
    '''

    #Define seaward boundary gw=setup
    GW[:,shl_ix] = setup
    GW[:,shl_ix-1] = setup
    
    if p['boundary_gw'].lower() == 'no_flow':
        #Define landward boundary dgw/dx=0
        GW[:,-1] = GW[:,-4] 
        GW[:,-2] = GW[:,-4] 
        GW[:,-3] = GW[:,-4]
    
    elif p['boundary_gw'].lower() == 'static':
       #Define landward boundary 
        GW[:,-1] = p['GW_stat']
        GW[:,-2] = p['GW_stat']
        GW[:,-3] = p['GW_stat'] 
        
    else:
         logger.log_and_raise('Unknown landward groundwater boundary condition' % p['boundary_gw'], exc=ValueError)
    
    # #Set GW levels to ground level within seepage face
    # ixs = np.argmin(GW + 0.001 >= s['zb'],axis=1)
    # for i in range(len(ixs)):
    #     if shl_ix[i] < ixs[i] - 1:
    #         GW[i,shl_ix[i]:ixs[i]-1] = s['zb'][i,shl_ix[i]:ixs[i]-1]
    
    #Compute groundwater level change dGW/dt (Boussinesq equation)
    dGW = np.zeros(s['gw'].shape)
    a = np.zeros(s['gw'].shape)
    b = np.zeros(s['gw'].shape)
    c = np.zeros(s['gw'].shape)

    

    for i in range(len(a[:,0])):
        #if shl_ix[i] > -1 & shl_ix[i] < len(a[0,:]) - 3:
        if shl_ix[i] < len(a[0,:]) - 3:   
            a[i,shl_ix[i]:-2]=(GW[i,shl_ix[i]+1:-1] - 2 * GW[i,shl_ix[i]:-2] + GW[i,shl_ix[i]-1:-3]) / s['ds'][i,shl_ix[i]:-2] ** 2
            b[i,shl_ix[i]:-2]=(GW[i,shl_ix[i]:-2] * (GW[i,shl_ix[i]+1:-1] + GW[i,shl_ix[i]-1:-3])) / s['ds'][i,shl_ix[i]:-2]                  
            c[i,shl_ix[i]+1:-3]=(b[i,shl_ix[i]+2:-2]-b[i,shl_ix[i]:-4])/s['ds'][i,shl_ix[i]+1:-3]
            dGW[i,shl_ix[i]+1:-3]=p['K_gw'] / p['ne_gw'] * (p['D_gw'] * a[i,shl_ix[i]+1:-3] + c[i,shl_ix[i]+1:-3])
    return dGW

#@jit
def it_hdelta(hdelta,s,p,i,j):
    '''
    Add description
    
    '''
    F_hdelta =1
    while F_hdelta > 0.01:           
        hdelta = hdelta + 0.01
        w_hdelta = (p['resw_moist'] + (p['satw_moist'] - p['resw_moist']) / (1 + abs(p['alfaw_moist'] * hdelta) ** p['nw_moist']) ** p['mw_moist'])
        d_hdelta = (p['resd_moist'] + (p['satd_moist'] - p['resd_moist']) / (1 + abs(p['alfad_moist'] * hdelta) ** p['nd_moist']) ** p['md_moist'])
        F_hdelta = s['w_h'][i,j] + (p['satw_moist'] - s['w_h'][i,j]) / np.maximum(p['satw_moist'] - w_hdelta,0.0001) * (d_hdelta - w_hdelta) - s['scan_d_moist'][i,j]

    return hdelta



def evaporation(s,p,met):
    '''Compute evaporation according to the Penman equation (Shuttleworth, 1993)

    Parameters
    ----------
    s : dict
        Spatial grids
    p : dict
        Model configuration parameters
    met : dict
        meteorologial parameters
        T: Temperature, degrees Celsius
        Q : Global radiation, MJ/m2/d
        P : Atmospheric pressure, kPa
        U: Relative humidity, %

    Returns
    -------
    float
        Evaporation (mm/day)

    '''
    l = 2.26 #latent heat of vaporization of water (MJ/kg)
    m = vaporation_pressure_slope(met['T']) # [kPa/K]
    delta = saturation_pressure(met['T']) * (1. - met['U'] / 100) # vapor pressure deficit [kPa]
    gamma = (p['cpair'] * met['P']) / (.622 * l) # [kPa/K]
    u2 = .174 / np.log10(p['z'] / 2.) * s['uw'] # [m/s]
    evo =(m * met['Q'] + 6.43 * gamma * delta * (1. + 0.86 * u2)) \
          / (l * (m + gamma))
    return evo




def vaporation_pressure_slope(T):
    '''Compute vaporation pressure slope based on air temperature

    Parameters
    ----------
    T : float
        Air temperature in degrees Celcius

    Returns
    -------
    float
        Vaporation pressure slope

    '''
    
    # Tetens, 1930; Murray, 1967
    s = 4098. * saturation_pressure(T) / (T + 237.3)**2 # [kPa/K]

    return s


def saturation_pressure(T):
    '''Compute saturation pressure based on air temperature, Tetens equation

    Parameters
    ----------
    T : float
        Air temperature in degrees Celcius

    Returns
    -------
    float
        Saturation pressure

    '''


    vp = 0.6108 * np.exp(17.27 * T / (T + 237.3)) # [kPa]

    return vp
