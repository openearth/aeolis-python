import math
from xml.parsers.expat import model

import numpy as np 

def update(self):

    veg_mortality_fix(self.get_var('zb'),self.get_var('TWL'), self)

    scr = self.p['shoreline_change_rate']
    zb = wet_supply(self.get_var('x'), self.get_var('zb'), self.p['beach_slope'], scr, self.p['dune_toe_elevation'], self)
    self.s['zb'] = zb

def veg_mortality_fix(zb, TWL, self):
    #Fixes mortality bug within the model, only killing vegetation seaward of the TWL and profile intersection'
    veg_min_elevation = self.p['veg_min_elevation']
    process_tide = self.get_var('process_tide')

    rhoveg = self.get_var('rhoveg')
    hveg = self.get_var('hveg')
    vegetated = self.get_var('vegetated')

    if process_tide:

        elev_dry = zb>=TWL

        try:
            if elev_dry[0][0] == True and TWL[0][0] < veg_min_elevation:
                #exception if TWL is below bathymetry & no intersection 
                limit = 0 
            else:
                limit = int(np.where(np.diff(elev_dry))[1][0]) # finds most seaward intersection of TWL and zb 

            ix_flooded1 = (zb[:,:limit] < TWL[:,:limit]) # identifies flodded area before limit
            
            rest = np.full((3, (len(TWL[0])-limit)), False, dtype=bool) # creates a False array to fill in the rest of the profile 
            
            ix_flooded = np.concatenate((ix_flooded1, rest), axis=1) # adds the arrays together 
        
        except:
            ix_flooded = (zb < TWL)
            rhoveg[ix_flooded]     = 0. 
            hveg[ix_flooded]       = 0.
            vegetated[ix_flooded]  = False
            # s['lateral'][ix_flooded]    = False
    self.set_var('rhoveg', rhoveg)
    self.set_var('hveg', hveg)
    self.set_var('vegetated', vegetated) 
    # return rhoveg, hveg, vegetated


def wet_supply(x_0, zb, beach_slope, shoreline_change_rate, dune_toe_elevation, self):

    ''' Increase elevation of beach topography.

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
    process_wet_supply = self.p['process_wet_supply']
    method_wet_supply = self.p['method_wet_supply']
    zshoreline = self.p['zshoreline']
    # xshoreline = self.p['xshoreline']
    #IN SOURCE CODE!!!!!!!!!!!!
    # Original wet-bed-reset function; basic resetting of the bed when inundated
    # if p['process_wet_supply'] or p['process_wet_bed_reset']:

    #     if p['method_wet_supply'] == 'wet_bed_reset':            
    #         Tbedreset = p['dt_opt'] / p['Tbedreset'] # []s
            
    #         ix = s['TWL'] > (s['zb'])
    #         s['zb'][ix] += (s['zb0'][ix] - s['zb'][ix]) * Tbedreset
    
    if process_wet_supply:

        if method_wet_supply == 'vertical_beach_growth':
            beach_inc = shoreline_change_rate*math.cos((math.pi/2)-math.atan(beach_slope))
            vrate = (beach_inc*(1/365.25/24/3600))*self.p['dt'] #(m/timestep)
            ny, nx = zb.shape  

            for iy in range(ny):
                x_all = x_0[iy,:]
                zb_all = zb[iy,:]               
                xi =  (zb_all < dune_toe_elevation)
                beach_z = zb_all[xi]
                b = beach_z[0] + vrate
                x = x_all[xi] 
                slope = beach_slope
                new_beach = slope*x + b 
                zb[iy,xi] = new_beach    

        if method_wet_supply == 'constant_SCR_constant_tanB':

            beach_inc = shoreline_change_rate*math.cos((math.pi/2)-math.atan(beach_slope))
            vrate = (beach_inc/(365.25*24*3600))*self.p['dt'] #(m/timestep)
            ny, nx = zb.shape  

            for iy in range(ny):
                x_all = x_0[iy,:]
                zb_all = zb[iy,:]

                xi = zb_all < dune_toe_elevation
                beach_z = zb_all[xi]
                x = x_all[xi]

                xi3 = np.where(beach_z > zshoreline)
                xi3 = xi3[0][0]
                b = beach_z[xi3] + vrate

                new_temp_beach = beach_slope*(x-x[xi3]) + b

                xi2 = new_temp_beach <= np.min(zb_all)
                new_temp_beach[xi2] = np.min(zb_all)

                zb[iy,xi]= new_temp_beach

        if method_wet_supply == 'constant_SCR_variable_tanB':
            beach_inc = shoreline_change_rate*math.cos((math.pi/2)-math.atan(beach_slope))
            vrate = (beach_inc/(365.25*24*3600))*self.p['dt'] #(m/timestep)
            hrate = (shoreline_change_rate/(365.25*24*3600))*self.p['dt'] #(m/timestep)
            ny, nx = zb.shape  

            for iy in range(ny):
                # print('x=' + str(x))
                x_all = x_0[iy,:]
                zb_all = zb[iy,:]

                xi = zb_all <= dune_toe_elevation
                beach_z = zb_all[xi]

                x = x_all[xi]

                xi3 = np.where(beach_z > zshoreline)
                xi3 = xi3[0][0]
                beach_x = x-x[xi3]

                xy1 = ((np.min(x[xi3])),np.min(beach_z[xi3]))
                xy2 = (np.max(x), np.max(beach_z))

                new_slope = (xy2[1]-xy1[1])/(xy2[0]-(xy1[0]))
                b = np.min(beach_z[xi3]) + vrate
                new_temp_beach = new_slope*(beach_x) + b

                xi2 = new_temp_beach <= np.min(zb_all)
                new_temp_beach[xi2] = np.min(zb_all)

                xi4 = new_temp_beach > dune_toe_elevation
                new_temp_beach[xi4] = dune_toe_elevation
                zb[iy,xi]= new_temp_beach
    return zb
