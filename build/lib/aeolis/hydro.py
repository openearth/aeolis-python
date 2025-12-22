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
from matplotlib import pyplot as plt
from numba import njit
from scipy.interpolate import NearestNDInterpolator

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
        # Check if SWL or zs are not provided by some external model
        # In that case, skip initialization
        if ('zs' not in p['external_vars']) :
            if p['tide_file'] is not None:
                s['SWL'][:,:] = interp_circular(t,
                                            p['tide_file'][:,0],
                                            p['tide_file'][:,1])
            else:
                s['SWL'][:,:] = 0.

            # apply complex mask
            s['SWL'] = apply_mask(s['SWL'], s['tide_mask'])

        # External model input:
        elif ('zs' in p['external_vars']):
            

            s['SWL'] = s['zs'][:]
            s['SWL'] = apply_mask(s['SWL'], s['tide_mask'])

            # The dry points have to be filtered out, to prevent issues with run-up calculation later
            iwet = s['zs'] - s['zb'] > 2. * p['eps']
            s['SWL'][~iwet] = np.NaN
            mask = np.where(~np.isnan(s['SWL']))
            interp = NearestNDInterpolator(np.transpose(mask), s['SWL'][mask])
            s['SWL'] = interp( * np.indices(s['SWL'].shape))

            # fig, ax = plt.subplots()
            # pc = plt.pcolormesh(s['x'], s['y'], s['SWL'])#, vmin=1, vmax=1.3)
            # ax.set_aspect('equal')
            # fig.colorbar(pc, ax=ax)
            # plt.show()

            print('!Be carefull, according to current implementation of importing waterlevel from Flexible Mesh, SWL is equal to DSWL = zs!')
            logger.warning('!Be carefull, according to current implementation of importing waterlevel from Flexible Mesh, SWL is equal to DSWL = zs!')


    else:
        s['SWL'] = s['zb'] * 0.

    # Check if Hs or Tp are not provided by some external model
    # In that case, skip initialization
    if ('Hs' not in p['external_vars']) and ('Tp' not in p['external_vars']):

        if p['process_wave'] and p['wave_file'] is not None:

            # First compute wave height, than run-up + set-up and finally wave height including set-up for mixing

            # determine water depth
            h = np.maximum(0., s['SWL'] - s['zb'])
        
            s['Hs'][:,:] = interp_circular(t,
                                        p['wave_file'][:,0],
                                        p['wave_file'][:,1])

            s['Tp'][:,:] = interp_circular(t,
                                        p['wave_file'][:,0],
                                        p['wave_file'][:,2])

            s['Hs'] = np.minimum(h * p['gamma'], s['Hs'])

            # apply complex mask
            s['Hs'] = apply_mask(s['Hs'], s['wave_mask'])
            s['Tp'] = apply_mask(s['Tp'], s['wave_mask'])

        else:
            s['Hs'] = s['zb'] * 0.
            s['Tp'] = s['zb'] * 0.

    # apply complex mask (also for external model input)
    else:
        s['Hs'] = apply_mask(s['Hs'], s['wave_mask'])
        s['Tp'] = apply_mask(s['Tp'], s['wave_mask'])


    if p['process_runup']:
        ny = p['ny']

        if ('Hs' not in p['external_vars']):

            for iy in range(ny + 1):  # do this computation seperately on every y for now so alongshore variable wave runup can be added in the future
            
                hs = s['Hs'][iy][0]
                tp = s['Tp'][iy][0]
                wl = s['SWL'][iy][0]

                eta, sigma_s, R = calc_runup_stockdon(hs, tp, p['beach_slope'])
                s['R'][iy][:] = R
                s['eta'][iy][:] = eta
                s['sigma_s'][iy][:] = sigma_s
                
                if hasattr(s['runup_mask'], "__len__"):
                    s['eta'][iy][:] = apply_mask(s['eta'][iy][:], s['runup_mask'][iy][:])
                    s['R'][iy][:] = apply_mask(s['R'][iy][:], s['runup_mask'][iy][:])

                s['TWL'][iy][:] = s['SWL'][iy][:]  + s['R'][iy][:]
                s['DSWL'][iy][:] = s['SWL'][iy][:] + s['eta'][iy][:]            # Was s['zs'] before

        if ('Hs' in p['external_vars']):

            eta, sigma_s, R = calc_runup_stockdon(s['Hs'], s['Tp'], p['beach_slope'])
            s['R'][:] = R

            if hasattr(s['runup_mask'], "__len__"):
                s['eta'] = apply_mask(s['eta'], s['runup_mask'])
                s['R'] = apply_mask(s['R'], s['runup_mask'])

            s['TWL'][:] = s['SWL'][:]  + s['R'][:]
            s['DSWL'][:] = s['SWL'][:] # + s['eta'][:]       # DSWL is actually provided by FM (?)


    if p['process_wave'] and p['wave_file'] is not None:

        h_mix = np.maximum(0., s['TWL'] - s['zb'])

        s['Hsmix'][:,:] = interp_circular(t,
                                       p['wave_file'][:,0],
                                       p['wave_file'][:,1])

        s['Hsmix'] = np.minimum(h_mix * p['gamma'], s['Hsmix'])

        # apply complex mask
        s['Hsmix'] = apply_mask(s['Hsmix'], s['wave_mask'])

        
    if p['process_moist'] and p['method_moist_process'].lower() == 'surf_moisture' and p['meteo_file'] is not None: 

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


    # Ensure compatibility with XBeach: zs >= zb
    s['zs'] = s['SWL'].copy()
    ix = (s['zb'] > s['zs'])
    s['zs'][ix] = s['zb'][ix]

    return s


def update(s, p, dt,t):
    '''Update soil moisture content

    Updates soil moisture content in all cells. The soil moisure
    content is computed either with the infiltration-method or 
    surface_moist method. The infiltration method accounts for surface moisture
    as a function of runup and the subsequent infiltration and evaporation.
    The surface_moist method takes into account the effect of wave runup, 
    precipitation, evaporation, infiltration, and capillary rise from the 
    groundwater table. 
    
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
        
        #Initialize GW levels
        if t == p['tstart']:
            s['gw'][:,:] = p['in_gw']
            s['gw_prev'] = s['gw']
            s['wetting'] = s['gw'] > s['gw_prev']       
            
        #Specify wetting or drying conditions in previous timestep
        s['wetting'] = s['gw'] > s['gw_prev']
        #Save groundwater from previous timestep
        s['gw_prev'] = s['gw']        
        
        #Decrease timestep for GW computations
        dt_gw = int(dt / p['tfac_gw'])
        for i in range(int(dt / dt_gw)):
            t_gw = t + i * dt_gw
            interpolate(s,p,t_gw)
                
            #Define index of shoreline location
            shl_ix =np.argmax(s['zb'] > s['DSWL'],axis=1) - 1
            

            #Define index of runup limit
            runup_ix =np.argmax(s['zb'] > s['TWL'],axis=1) - 1
            
            # Landward boundary condition
            if p['boundary_gw'].lower() == 'no_flow':
                #Define landward boundary dgw/dx=0
                bound = 0
            
            elif p['boundary_gw'].lower() == 'static':
               #Define landward boundary 
                bound = 1
                
            else:
                 logger.log_and_raise('Unknown landward groundwater boundary condition' % p['boundary_gw'], exc=ValueError)

        
            #Runge-Kutta timestepping
            f1 = Boussinesq(s['gw'],s['DSWL'], s['ds'], p['GW_stat'], p['K_gw'], p['ne_gw'], p['D_gw'],shl_ix, bound,s['zb'],p['process_seepage_face'])
            f2 = Boussinesq(s['gw'] + dt_gw / 2 * f1,s['DSWL'], s['ds'], p['GW_stat'], p['K_gw'], p['ne_gw'], p['D_gw'], shl_ix, bound,s['zb'],p['process_seepage_face'])
            f3 = Boussinesq(s['gw'] + dt_gw / 2 * f2,s['DSWL'], s['ds'], p['GW_stat'], p['K_gw'], p['ne_gw'], p['D_gw'], shl_ix, bound,s['zb'],p['process_seepage_face'])
            f4 = Boussinesq(s['gw'] + dt_gw * f3,s['DSWL'], s['ds'], p['GW_stat'], p['K_gw'], p['ne_gw'], p['D_gw'], shl_ix, bound,s['zb'],p['process_seepage_face'])
            
            #Update groundwater level
            s['gw'] = s['gw'] + dt_gw / 6 * (f1 + 2 * f2 + 2 * f3 + f4)
        
            #Add infiltration from wave runup according to Nielsen (1990)
            if p['process_wave']:              
                #Compute f(x) = distribution of infiltrated water
                fx=np.zeros(s['gw'].shape)
                fx_ix=np.zeros_like(shl_ix)
                runup_overheight_distr(fx, fx_ix, shl_ix, runup_ix, s['x'])
                    
                # Update groundwater level with overheight due to runup
                s['gw'] = s['gw'] + p['Cl_gw'] * fx

            # Apply GW complex mask
            s['gw'] = apply_mask(s['gw'], s['gw_mask'])
            
                
            # Do not allow GW levels above ground level
            s['gw']=np.minimum(s['gw'], s['zb'])

            # Define cells below setup level
            ixg=s['zb'] < s['DSWL']
            
            # Set gw level to setup level in cells below setup level
            
            s['gw'][ixg]=s['DSWL'][ixg]            



    # Compute surface moisture with infiltration method using Darcy
    if p['process_moist']:
        if p['method_moist_process'].lower() == 'infiltration':
            F1 = -np.log(.5) / p['Tdry']
            ix = s['TWL'] - s['zb'] > p['eps']
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
            ix = s['TWL'] > s['zb']
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
                    s['h_delta'] = hdelta(s['wetting'], s['scan_w'], s['gw'],s['gw_prev'],s['scan_d'],s['h_delta'],s['zb'],s['scan_w_moist'],s['w_h'],p['satd_moist'],s['d_h'],p['satw_moist'],p['alfaw_moist'],p['resw_moist'],p['mw_moist'],p['nw_moist'])
                    
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
                    s['moist_swr'], s['scan_d'], s['scan_w'] = SWR_curve(s['wetting'],s['gw'],s['gw_prev'],s['scan_w'],s['moist_swr'],s['w_h'],s['scan_d'],s['scan_w_moist'],s['d_h'],s['scan_d_moist'])
                
            else:
                ixw = s['wetting'] == True
                s['moist_swr'][ixw] = p['resw_moist'] + (p['satw_moist'] - p['resw_moist']) \
                                      / (1 + abs(p['alfaw_moist'] * h[ixw]) ** p['nw_moist']) ** p['mw_moist']
                s['moist_swr'][~ixw] = p['resd_moist'] + (p['satd_moist'] - p['resd_moist']) \
                                      / (1 + abs(p['alfad_moist'] * h[~ixw]) ** p['nd_moist']) ** p['md_moist']
                
            #Update surface moisture with respect to capillary processes
            s['moist'] = np.minimum(np.maximum(s['moist'],s['moist_swr']),p['satd_moist'])
            

        
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



@njit
def Boussinesq (GW, DSWL, ds, GW_stat, K_gw, ne_gw, D_gw,shl_ix, bound,zb,process_seepage_face):
    '''
    Add description
    
    '''

    #Define seaward boundary gw=setup
    for i in range (len(GW[:,0])):
        GW[i,shl_ix[i]] = DSWL[i,shl_ix[i]]
        GW[i,shl_ix[i-1]] = DSWL[i,shl_ix[i-1]]
    # GW[:,shl_ix] = s['DSWL'][:,shl_ix]
    # GW[:,shl_ix-1] = s['DSWL'][:,shl_ix-1]

    if bound == 0:
        #Define landward boundary dgw/dx=0
        GW[:,-1] = GW[:,-4] 
        GW[:,-2] = GW[:,-4] 
        GW[:,-3] = GW[:,-4]
    
    elif bound == 1:
       #Define landward boundary 
        GW[:,-1] = GW_stat
        GW[:,-2] = GW_stat
        GW[:,-3] = GW_stat 
        
    
    #Set GW levels to ground level within seepage face
    if process_seepage_face:
        ixs = np.argmin(GW + 0.05 >= zb,axis=1)
        for i in range(len(ixs)):
            if shl_ix[i] < ixs[i] - 1:
                GW[i,shl_ix[i]:ixs[i]-1] = zb[i,shl_ix[i]:ixs[i]-1]
    
    #Compute groundwater level change dGW/dt (Boussinesq equation)
    dGW = np.zeros(GW.shape)
    a = np.zeros(GW.shape)
    b = np.zeros(GW.shape)
    c = np.zeros(GW.shape)
    



    for i in range(len(a[:,0])):
        if shl_ix[i] < len(a[0,:]) - 3:   
            a[i,shl_ix[i]:-2]=(GW[i,shl_ix[i]+1:-1] - 2 * GW[i,shl_ix[i]:-2] + GW[i,shl_ix[i]-1:-3]) / ds[i,shl_ix[i]:-2] ** 2
            b[i,shl_ix[i]:-2]=(GW[i,shl_ix[i]:-2] * (GW[i,shl_ix[i]+1:-1] + GW[i,shl_ix[i]-1:-3])) / ds[i,shl_ix[i]:-2]                  
            c[i,shl_ix[i]+1:-3]=(b[i,shl_ix[i]+2:-2]-b[i,shl_ix[i]:-4])/ds[i,shl_ix[i]+1:-3]
            dGW[i,shl_ix[i]+1:-3]=K_gw / ne_gw * (D_gw * a[i,shl_ix[i]+1:-3] + c[i,shl_ix[i]+1:-3])

    return dGW

@njit
def runup_overheight_distr(fx, fx_ix,shl_ix,runup_ix, x):
    '''
    Add description
    
    '''
    for i in range(len(fx[:,0])):
        #Define index of peak f(x)
        fx_ix[i] = (shl_ix[i]) + (2/3 * (runup_ix[i] - shl_ix[i]))
        #Compute f(X)
        fx[i,shl_ix[i]:fx_ix[i]] = (x[i,shl_ix[i]:fx_ix[i]] - x[i,shl_ix[i]]) / (2 / 3 * (x[i,runup_ix[i]] - x[i,shl_ix[i]]))
        fx[i,fx_ix[i]+1:runup_ix[i]] = 3 - (x[i,fx_ix[i]+1:runup_ix[i]]- x[i,shl_ix[i]])  / (1 / 3 * (x[i,runup_ix[i]] - x[i,shl_ix[i]]))
        fx[i,fx_ix[i]]=1
    return fx

@njit
def hdelta(wetting, scan_w, gw,gw_prev,scan_d,h_delta,zb,scan_w_moist,w_h,satd_moist,d_h,satw_moist,alfaw_moist,resw_moist,mw_moist,nw_moist):
    for i in range(len(wetting[:,0])):
        for j in range(len(wetting[0,:])):
            #Compute h delta on the main drying and wetting curve
            if scan_w[i,j] == False and wetting[i,j] == True and gw[i,j] < gw_prev[i,j] or scan_d[i,j] == False and wetting[i,j] == False and gw[i,j] > gw_prev[i,j]:
                h_delta[i,j]=np.maximum(0,(zb[i,j] - gw[i,j])*100)
            #Compute h_delta if there is a reversal on the wetting scanning curve
            if scan_w[i,j] == True and wetting[i,j] == True and gw[i,j] < gw_prev[i,j]:
                #Solve hdelta from drying scanning curve for which moist(h) on drying scanning curve equals moist(h) on wetting scanning curve 
                #intermediate solution:
                w_hdelta_int = np.minimum((scan_w_moist[i,j] - w_h[i,j]) * (satd_moist - w_h[i,j]) / (d_h[i,j] - w_h[i,j]) + w_h[i,j], satw_moist)
                #Solve hdelta from wetting curve
                h_delta[i,j] = np.maximum( 1 / alfaw_moist * (((satw_moist - resw_moist) \
                    / np.maximum((w_hdelta_int - resw_moist),0.00001)) ** (1 / mw_moist) - 1) ** (1 / nw_moist),0)
            #Compute h_delta if there is a reversal on the drying scanning curve
            if scan_d[i,j] == True and wetting[i,j] == False and gw[i,j] > gw_prev[i,j]:
            #Solve hdelta from wetting scanning curve for which moist(h) on wetting scanning curve equals moist(h) on drying scanning curve
                #Simple iteration method
                hdelta_it=0 #initialize hdelta
                F_hdelta =1
                while F_hdelta > 0.01:           
                    hdelta_it = hdelta_it + 0.01
                    w_hdelta = (resw_moist + (satw_moist - resw_moist) / (1 + np.abs(alfaw_moist * hdelta_it) ** nw_moist) ** mw_moist)
                    d_hdelta = (resd_moist + (satd_moist - resd_moist) / (1 + np.abs(alfad_moist * hdelta_it) ** nd_moist) ** md_moist)
                    F_hdelta = w_h + (satw_moist - w_h) / np.maximum(satw_moist - w_hdelta,0.0001) * (d_hdelta - w_hdelta) - scan_d_moist
                
                hdelta[i,j] = hdelta_it
                
    return hdelta

@njit
def SWR_curve(wetting,gw,gw_prev,scan_w,moist_swr,w_h,scan_d,scan_w_moist,d_h,scan_d_moist):
    for i in range(len(wetting[:,0])):
        for j in range(len(wetting[0,:])):
            #Wetting conditions main curve
            if gw[i,j] >= gw_prev[i,j] and wetting[i,j] == True and scan_w[i,j] == False:
                moist_swr[i,j]=w_h[i,j]
                scan_w[i,j] = False
                scan_d[i,j] = False
            #wetting conditions, timestep of reversal - move onto wetting scanning curve
            elif gw[i,j] >= gw_prev[i,j] and wetting[i,j] == False:
                moist_swr[i,j] = scan_w_moist[i,j]
                scan_w[i,j] = scan_w_moist[i,j] > w_h[i,j]
                scan_d[i,j] = False
            #wetting conditions - followed a wetting scanning curve in previous timestep - continue following scanning curve unless main curve is reached
            elif gw[i,j] >= gw_prev[i,j] and wetting[i,j] == True and scan_w[i,j] == True:
                moist_swr[i,j] = scan_w_moist[i,j]
                scan_w[i,j] = scan_w_moist[i,j] > w_h[i,j]
                scan_d[i,j] = False
            #Drying conditions main curve
            elif gw[i,j] < gw_prev[i,j] and wetting[i,j] == False and scan_d[i,j] == False:
                moist_swr[i,j]=d_h[i,j]
                scan_d[i,j] = False
                scan_w[i,j] = False
            #Drying conditions, timestep of reversal - move onto a drying scanning curve
            elif gw[i,j] < gw_prev[i,j] and wetting[i,j] == True:
                moist_swr[i,j] = scan_d_moist[i,j]
                scan_d[i,j] = scan_d_moist[i,j] < d_h[i,j]
                scan_w[i,j] = False
            #Drying conditions - followed a drying scanning curve in previous timestep - continue following scanning curve unless main curve is reached
            elif gw[i,j] < gw_prev[i,j] and wetting[i,j] == False and scan_d[i,j] == True:
                moist_swr[i,j] = scan_d_moist[i,j]
                scan_d[i,j] = scan_d_moist[i,j] < d_h[i,j]
                scan_w[i,j] = False
                                
    return moist_swr, scan_d, scan_w

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

def calc_runup_stockdon(Ho, Tp, beta):
    """
    Calculate runup according to /Stockdon et al 2006.
    """
    
    if hasattr(Ho, "__len__"):

        R = np.zeros(np.shape(Ho))
        sigma_s = np.zeros(np.shape(Ho))
        eta = np.zeros(np.shape(Ho))

        Lo = 9.81 * Tp * Tp / (2 * np.pi) #wavelength
        iribarren = beta / (Ho / Lo) ** (0.5) #irribarren number

        i_iri = (Ho > 0) * (iribarren < 0.3)
        R[i_iri] = 0.043 * np.sqrt(Ho[i_iri] * Lo[i_iri]) #formula for dissipative conditions
        sigma_s[i_iri] = 0.046 * np.sqrt(Ho[i_iri] * Lo[i_iri]) /2
        eta[i_iri] = R[i_iri] - sigma_s[i_iri]

        i_iri = (Ho > 0) * (iribarren > 0.3)
        nsigma = 2  # nsigma=1 for R16% and nsigma=2 for R2%
        eta[i_iri] = 0.35 * beta * np.sqrt(Ho[i_iri] * Lo[i_iri])
        sigma_s[i_iri] = np.sqrt(Ho[i_iri] * Lo[i_iri] * (0.563 * (beta * beta) + 0.0004)) * nsigma / 2 / 2
        R[i_iri] = 1.1 * (eta[i_iri] + sigma_s[i_iri]) #result for non-dissipative conditions

    else:
        if Ho > 0 and Tp > 0 and beta > 0:
            Lo = 9.81 * Tp * Tp / (2 * np.pi) #wavelength
            iribarren = beta / (Ho / Lo) ** (0.5) #irribarren number

            if iribarren < 0.3:
                R = 0.043 * np.sqrt(Ho * Lo) #formula for dissipative conditions
                sigma_s = 0.046 * np.sqrt(Ho * Lo) /2
                eta = R - sigma_s
            else:
                nsigma = 2  # nsigma=1 for R16% and nsigma=2 for R2%
                Lo = 9.81 * Tp * Tp /(2 * np.pi)
                eta = 0.35 * beta * np.sqrt(Ho * Lo)
                sigma_s = np.sqrt(Ho * Lo * (0.563 * (beta * beta) + 0.0004)) * nsigma / 2 / 2
                R = 1.1 * (eta + sigma_s) #result for non-dissipative conditions
        else:
            R = 0
            sigma_s = 0
            eta = 0

    return eta, sigma_s, R











