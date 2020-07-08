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
import numpy as np

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
        else:
            s['zs'][:,:] = 0.

        # apply complex mask
        s['zs'] = apply_mask(s['zs'], s['tide_mask'])
        
        #define SWL
        s['swl'] = s['zs']

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
        
        # maximize wave height by depth ratio ``gamma``
        s['Hs'] = np.minimum(h * p['gamma'], s['Hs'])
        
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
    s['zs'] = np.maximum(s['zs'], s['zb'])

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
        
        #Decrease timestep for GW computations
        dt_gw = int(dt / p['tfac_gw'])
        for i in range(int(dt / dt_gw)):
            t_gw = t + i * dt_gw
            interpolate(s,p,t_gw)
        
            #Compute setup
            setup=np.average(s['swl']) + 0.35 * p['xi'] #Stockdon et al (2006) CH:Should we include a function to compute the Irribaren number?

            #Initialize GW levels and h_delta (GW depth when changing from wetting/drying)
            if t==0:
                s['gw'][:,:]=setup
                s['h_delta'][:,:]=np.maximum(0,(s['zb'] - s['gw'])*100)
                
            #Define index of shoreline location
            shl_ix = np.argmax(s['zb'] > setup,axis=1) - 1

        
            #Runge-Kutta timestepping
            f1 = Boussinesq(s['gw'],s,p,setup, shl_ix)
            f2 = Boussinesq(s['gw'] + dt_gw / 2 * f1,s,p,setup, shl_ix)
            f3 = Boussinesq(s['gw'] + dt_gw / 2 * f2,s,p,setup, shl_ix)
            f4 = Boussinesq(s['gw'] + dt_gw * f3,s,p,setup, shl_ix)
            
            #Update groundwater level
            s['gw'] = s['gw'] + dt_gw / 6 * (f1 + 2 * f2 + 2 * f3 + f4)
        
            #Add infiltration from wave runup according to Nielsen (1990)
            if p['process_wave']:
                #Define index of runup limit
                runup_ix = np.argmax(s['zb'] >= s['zs'],axis=1) - 1
                
                #Compute f(x) = distribution of infiltrated water
                fx=np.zeros(s['gw'].shape)
                fx_ix=np.zeros_like(runup_ix)
                for i in range(len(s['gw'][:,0])):
                    #Define index of peak f(x)
                    fx_ix[i] = (shl_ix[i]) + (2/3 * (runup_ix[i] - shl_ix[i]))
                    #Compute f(X)
                    fx[i,shl_ix[i]:fx_ix[i]] = (s['x'][i,shl_ix[i]:fx_ix[i]] - s['x'][i,shl_ix[i]]) / (2 / 3 * (s['x'][i,runup_ix[i]] - s['x'][i,shl_ix[i]]))
                    fx[i,fx_ix[i]+1:runup_ix[i]] = 3 - (s['x'][i,fx_ix[i]+1:runup_ix[i]]- s['x'][i,shl_ix[i]])  / (1 / 3 * (s['x'][i,runup_ix[i]] - s['x'][i,shl_ix[i]]))
                    

                # Update groundwater level with overheight due to runup
                s['gw'] = s['gw'] + p['Cl_gw'] * fx * p['K_gw'] / p['ne_gw']
            
            # Define cells below setup level
            ixg=s['zb'] < setup
            
            #Set GW level to the setup level in submerged cells
            s['gw'][ixg] = setup 
                
            # Do not allow GW levels above ground level in areas that are not submerged (=above setup level)
            s['gw'][~ixg]=np.minimum(s['gw'][~ixg], s['zb'][~ixg])
            
        #Update h_delta if there is a reversal between wetting/drying conditions
        for i in range(len(s['wetting'][:,0])):
            for j in range(len(s['wetting'][0,:])):
                if (s['wetting'][i,j] == True and s['gw'][i,j] < s['gw_prev'][i,j]) or (s['wetting'][i,j] == False and s['gw'][i,j] > s['gw_prev'][i,j]):
                    s['h_delta'][i,j]=np.maximum(0,(s['zb'][i,j] - s['gw'][i,j])*100)

        #Specify wetting or drying conditions
        s['wetting'] = s['gw'] > s['gw_prev']+0.00001
            
        #Save groundwater level for next timestep
        s['gw_prev'] = s['gw']
        


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
            if p['process_moist'] is None :
                logger.log_and_raise('process_groundwater is not activated, the groundwater level is not computed within the program but set constant at 0 m', exc=ValueError)
                
            #Cells that were flooded in previous time step (by rain or runup) drains to field capacity 
            s['moist'] = np.minimum(s['moist'],p['fc'])
            
            #If the cell is flooded (runup) in this timestep, assume satiation
            ix = s['zs'] > s['zb']
            s['moist'][ix] = p['sat_moist']

            
            #Update surface moisture with respect to evaporation, condensation, and precipitation
            met = s['meteo']
            evo = evaporation(s,p,met)
            evo = evo / 24. / 3600. / 1000. # convert evaporation from mm/day to m/s
            pcp = met['RH'] / 3600. / 1000. # convert precipitation from mm/hr to m/s
            s['moist'][~ix] = np.maximum(s['moist'][~ix] + (pcp - evo[~ix]) * dt / p['thick_moist'], p['res_moist'])
            s['moist'][~ix] = np.minimum(s['moist'][~ix],p['sat_moist'])
                        
            #Compute surface moisture due to capillary processes (van Genuchten and Mualem II)
            h=np.maximum(0,(s['zb'] - s['gw']) * 100) #h in cm to match convention of alfa (cm-1)
            
            #Initialize value of surface moisture due to capillary rise
            if t == 0:
                s['moist_swr'] = p['res_moist'] + (p['sat_moist'] - p['res_moist']) \
                      / (1 + abs(p['alfaw_moist'] * h) ** p['n_moist']) ** (1 - 1 / p['n_moist'])
            else:
                #Compute moisture of h for the wetting curve
                w_h = p['res_moist'] + (p['sat_moist'] - p['res_moist']) \
                                  / (1 + abs(p['alfaw_moist'] * h) ** p['n_moist']) ** (1 - 1 / p['n_moist'])
                #Compute moisture of h_delta for the wetting curve
                w_hdelta = p['res_moist'] + (p['sat_moist'] - p['res_moist']) \
                                  / (1 + abs(p['alfaw_moist'] * s['h_delta']) ** p['n_moist']) ** (1 - 1 / p['n_moist'])
                #Compute moisture of h for the drying curve
                d_h = p['res_moist'] + (p['sat_moist'] - p['res_moist']) \
                                  / (1 + abs(p['alfad_moist'] * h) ** p['n_moist']) ** (1 - 1 / p['n_moist'])
                #Compute moisture of h_delta for the drying curve
                d_hdelta = p['res_moist'] + (p['sat_moist'] - p['res_moist']) \
                                  / (1 + abs(p['alfad_moist'] * s['h_delta']) ** p['n_moist']) ** (1 - 1 / p['n_moist'])
                                  
                ixw = s['wetting'] == True    
                
                #Wetting (select largest of scanning curve and main wetting curve)
                s['moist_swr'][ixw] = np.maximum(w_h[ixw] + (p['sat_moist'] - w_h[ixw]) / np.maximum(p['sat_moist'] - w_hdelta[ixw],0.0001) \
                                                 * (d_hdelta[ixw] - w_hdelta[ixw]),w_h[ixw])
                #Drying (select smallest of scanning curve and main drying curve)
                s['moist_swr'][~ixw] = np.minimum(w_h[~ixw] + (w_hdelta[~ixw] - w_h[~ixw]) / np.maximum(p['sat_moist'] - w_h[~ixw],0.0001) \
                                                  * (d_h[~ixw] - w_h[~ixw]), d_h[~ixw])

            
            #Update surface moisture with respect to capillary processes
            s['moist'] = np.maximum(s['moist'],s['moist_swr'])
        

        
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
    
    #Define landward boundary dgw/dx=0
    GW[:,-1] = GW[:,-3]
    GW[:,-2] = GW[:,-3]
    
    # Set GW levels to ground level within seepage face
    ixs = np.argmin(GW + 0.001 >= s['zb'],axis=1)
    for i in range(len(ixs)):
        if shl_ix[i] < ixs[i] - 1:
            GW[i,shl_ix[i]:ixs[i]-1] = s['zb'][i,shl_ix[i]:ixs[i]-1]
    
    #Compute groundwater level change dGW/dt (Boussinesq equation)
    dGW = np.zeros(s['gw'].shape)
    a = np.zeros(s['gw'].shape)
    b = np.zeros(s['gw'].shape)
    c = np.zeros(s['gw'].shape)

    for i in range(len(a[:,0])):
        a[i,shl_ix[i]:-2]=(GW[i,shl_ix[i]+1:-1] - 2 * GW[i,shl_ix[i]:-2] + GW[i,shl_ix[i]-1:-3]) / s['ds'][i,shl_ix[i]:-2] ** 2
        b[i,shl_ix[i]:-2]=(GW[i,shl_ix[i]:-2] * (GW[i,shl_ix[i]+1:-1] + GW[i,shl_ix[i]-1:-3])) / s['ds'][i,shl_ix[i]:-2]                  
        c[i,shl_ix[i]+1:-3]=(b[i,shl_ix[i]+2:-2]-b[i,shl_ix[i]:-4])/s['ds'][i,shl_ix[i]+1:-3]
        dGW[i,shl_ix[i]+1:-3]=p['K_gw'] / p['ne_gw'] * (p['D_gw'] * a[i,shl_ix[i]+1:-3] + c[i,shl_ix[i]+1:-3])
    return dGW

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
