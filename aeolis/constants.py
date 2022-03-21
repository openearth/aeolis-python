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


#: Aeolis model state variables
INITIAL_STATE = {
    ('ny', 'nx') : (
        'uw',                               # [m/s] Wind velocity
        'uws',                              # [m/s] Component of wind velocity in x-direction
        'uwn',                              # [m/s] Component of wind velocity in y-direction
        'tau',                              # [m/s] Wind shear velocity
        'taus',                             # [m/s] Component of wind shear velocity in x-direction
        'taun',                             # [m/s] Component of wind shear velocity in y-direction
        'dtaus',                            # [-] Component of the wind shear perturbation in x-direction
        'dtaun',                            # [-] Component of the wind shear perturbation in y-direction
        'udir',                             # [rad] Wind direction
        'zs',                               # [m] Total water level above reference
        'swl',                              # [m] Still water level above reference
        'Hs',                               # [m] Wave height
        'zne',                              # NEW! [m] Non-erodible layer
    ),
}

MODEL_STATE = {
    ('ny', 'nx') : (
        'x',                                # [m] Real-world x-coordinate of grid cell center      
        'y',                                # [m] Real-world y-coordinate of grid cell center
# Gridparams 
        'xz',                               # [m] Real-world x-coordinate of grid cell center
        'xu',                               # [m] Real-world x-coordinates of u-points
        'xv',                               # [m] Real-world x-coordinates of v-points
        'xc',                               # [m] Real-world x-oordinates of c-points
        'yz',                               # [m] Real-world y-coordinate of grid cell center
        'yu',                               # [m] Real-world y-coordinates of u-points
        'yv',                               # [m] Real-world y-coordinates of v-points
        'yc',                               # [m] Real-world y-coordinates of c-points
        'ds',                               # [m] Real-world grid cell size in x-direction
        'dn',                               # [m] Real-world grid cell size in y-direction
        'dnz',                              # [m] Distances in n-direction
        'dnu',                              # [m] Distances in n-direction
        'dnv',                              # [m] Distances in n-direction
        'dnc',                              # [m] Distances in n-direction
        'dsz',                              # [m] Distances in s-direction
        'dsu',                              # [m] Distances in s-direction
        'dsv',                              # [m] Distances in s-direction
        'dsc',                              # [m] Distances in s-direction
        'dsdnz',                            # [m^2] Real-world grid cell surface area
        'dsdnzi',                           # [m^-2] Inverse of real-world grid cell surface area   
        'dsdn',                             # [m^2] Real-world grid cell surface area
        'dsdni',                            # [m^-2] Inverse of real-world grid cell surface area   
        'alfaz',                            # [rad] Real-world grid cell orientation around z
        'alfau',                            # [rad] Real-world grid cell orientation around u
        'alfav',                            # [rad] Real-world grid cell orientation around v    
        'alfa',                             # [rad] Real-world grid cell orientation #Sierd_comm in later releases this needs a revision 
        'zb',                               # [m] Bed level above reference
        'S',                                # [-] Level of saturation
        'moist',                            # [-] Moisture content (volumetric)
        'moist_swr',                        # [-] Moisture content soil water retention relationship (volumetric)
        'h_delta',                          # [-] Suction at reversal between wetting/drying conditions
        'gw',                               # [m] Groundwater level above reference
        'gw_prev',                          # [m] Groundwater level above reference in previous timestep
        'wetting',                          # [bool] Flag indicating wetting or drying of soil profile
        'scan_w',                           # [bool] Flag indicating that the moisture is calculated on the wetting scanning curve
        'scan_d',                           # [bool] Flag indicating that the moisture is calculated on the drying scanning curve
        'scan_w_moist',                     # [-] Moisture content (volumetric) computed on the wetting scanning curve
        'scan_d_moist',                     # [-] Moisture content (volumetric) computed on the drying scanning curve
        'w_h',                              # [-] Moisture content (volumetric) computed on the main wetting curve
        'd_h',                              # [-] Moisture content (volumetric) computed on the main drying curve
        'w_hdelta',                         # [-] Moisture content (volumetric) computed on the main wetting curve for hdelta
        'd_hdelta'                          # [-] Moisture content (volumetric) computed on the main drying curve for hdelta
    ), 
    ('ny','nx','nfractions') : (
        'Cu',                               # [kg/m^2] Equilibrium sediment concentration integrated over saltation height
        'Cuf',                              # [kg/m^2] Equilibrium sediment concentration integrated over saltation height, assuming the fluid shear velocity threshold
        'Ct',                               # [kg/m^2] Instantaneous sediment concentration integrated over saltation height
        'q',                                # [kg/m/s] Instantaneous sediment flux
        'qs',                               # [kg/m/s] Instantaneous sediment flux in x-direction
        'qn',                               # [kg/m/s] Instantaneous sediment flux in y-direction
        'pickup',                           # [kg/m^2] Sediment entrainment
        'w',                                # [-] Weights of sediment fractions
        'w_init',                           # [-] Initial guess for ``w''
        'w_air',                            # [-] Weights of sediment fractions based on grain size distribution in the air
        'w_bed',                            # [-] Weights of sediment fractions based on grain size distribution in the bed
        'uth',                              # [m/s] Shear velocity threshold
        'uthf',                             # [m/s] Fluid shear velocity threshold
    ),
    ('ny','nx','nlayers') : (
        'thlyr',                            # [m] Bed composition layer thickness
        'salt',                             # [-] Salt content
    ),
    ('ny','nx','nlayers','nfractions') : (
        'mass',                             # [kg/m^2] Sediment mass in bed
    )
}


#: AeoLiS model default configuration
DEFAULT_CONFIG = {
    'process_wind'                  : True,               # Enable the process of wind
    'process_shear'                 : False,              # Enable the process of wind shear
    'process_tide'                  : True,               # Enable the process of tides
    'process_wave'                  : True,               # Enable the process of waves
    'process_runup'                 : True,               # Enable the process of wave runup
    'process_moist'                 : True,               # Enable the process of moist
    'process_mixtoplayer'           : True,               # Enable the process of mixing
    'process_threshold'             : True,               # Enable the process of threshold
    'process_transport'             : True,               # Enable the process of transport
    'process_bedupdate'             : True,               # Enable the process of bed updating
    'process_salt'                  : False,              # Enable the process of salt
    'process_avalanche'             : True,               # NEW! Enable the process of avalanching
    'process_groundwater'           : False,              # Enable the process of groundwater
    'process_scanning'              : True,               # Enable the process of scanning curves
    'th_grainsize'                  : True,               # Enable wind velocity threshold based on grainsize
    'th_bedslope'                   : False,              # Enable wind velocity threshold based on bedslope
    'th_moisture'                   : True,               # Enable wind velocity threshold based on moisture
    'th_humidity'                   : False,              # Enable wind velocity threshold based on humidity
    'th_salt'                       : False,              # Enable wind velocity threshold based on salt
    'th_roughness'                  : True,               # Enable wind velocity threshold based on roughness
    'xgrid_file'                    : None,               # Filename of ASCII file with x-coordinates of grid cells
    'ygrid_file'                    : None,               # Filename of ASCII file with y-coordinates of grid cells
    'bed_file'                      : None,               # Filename of ASCII file with bed level heights of grid cells
    'wind_file'                     : None,               # Filename of ASCII file with time series of wind velocity and direction
    'tide_file'                     : None,               # Filename of ASCII file with time series of water levels
    'wave_file'                     : None,               # Filename of ASCII file with time series of wave heights
    'meteo_file'                    : None,               # Filename of ASCII file with time series of meteorlogical conditions
    'bedcomp_file'                  : None,               # Filename of ASCII file with initial bed composition
    'threshold_file'                : None,               # Filename of ASCII file with shear velocity threshold
    'ne_file'                       : None,               # NEW! Filename of ASCII file with non-erodible layer
    'gw_file'                       : None,               # Filename of ASCII file with time series of groundwater level in grid cells
    'wave_mask'                     : None,               # Filename of ASCII file with mask for wave height
    'tide_mask'                     : None,               # Filename of ASCII file with mask for tidal elevation
    'threshold_mask'                : None,               # Filename of ASCII file with mask for the shear velocity threshold
    'gw_mask'                       : None,               # Filename of ASCII file with mask for the groundwater level
    'dt'                            : 60.,                # [s] Time step size
    'CFL'                           : 1.,                 # [-] CFL number to determine time step in explicit scheme
    'accfac'                        : 1.,                 # [-] Numerical acceleration factor
    'tstart'                        : 0.,                 # [s] Start time of simulation
    'tstop'                         : 3600.,              # [s] End time of simulation
    'restart'                       : None,               # [s] Interval for which to write restart files
    'output_times'                  : 60.,                # [s] Output interval in seconds of simulation time
    'output_file'                   : None,               # Filename of netCDF4 output file
    'output_vars'                   : ['zb', 'zs',
                                       'Ct', 'Cu',
                                       'uw', 'uth',
                                       'mass', 'pickup'], # Names of spatial grids to be included in output
    'output_types'                  : [],                 # Names of statistical parameters to be included in output (avg, sum, var, min or max)
    'grain_size'                    : [225e-6],           # [m] Average grain size of each sediment fraction
    'grain_dist'                    : [1.],               # [-] Initial distribution of sediment fractions
    'nlayers'                       : 3,                  # [-] Number of bed layers
    'layer_thickness'               : .01,                # [m] Thickness of bed layers
    'g'                             : 9.81,               # [m/s^2] Gravitational constant
    'rhoa'                          : 1.25,               # [kg/m^3] Air density
    'rhop'                          : 2650.,              # [kg/m^3] Grain density
    'rhow'                          : 1025.,              # [kg/m^3] Water density
    'porosity'                      : .4,                 # [-] Sediment porosity
    'A'                             : .085,               # [-] Constant in formulation for wind velocity threshold based on grain size
    'z'                             : 10.,                # [m] Measurement height of wind velocity
    'h'                             : None,               # [m] Representative height of saltation layer
    'k'                             : 0.01,               # [m] Bed roughness
    'Cb'                            : 1.5,                # [-] Constant in formulation for equilibrium sediment concentration
    'm'                             : .5,                 # [-] Factor to account for difference between average and maximum shear stress
    'sigma'                         : 4.2,                # [-] Ratio between basal area and frontal area of roughness elements
    'beta'                          : 130.,               # [-] Ratio between drag coefficient of roughness elements and bare surface
    'bi'                            : 1.,                 # [-] Bed interaction factor
    'T'                             : 1.,                 # [s] Adaptation time scale in advection equation
    'Tdry'                          : 3600.*1.5,          # [s] Adaptation time scale for soil drying
    'Tsalt'                         : 3600.*24.*30.,      # [s] Adaptation time scale for salinitation
    'eps'                           : 1e-3,               # [m] Minimum water depth to consider a cell "flooded"
    'gamma'                         : .5,                 # [-] Maximum wave height over depth ratio
    'xi'                            : .3,                 # [-] Surf similarity parameter
    'facDOD'                        : .1,                 # [-] Ratio between depth of disturbance and local wave height
    'csalt'                         : 35e-3,              # [-] Maximum salt concentration in bed surface layer
    'cpair'                         : 1.013e-3,           # [MJ/kg/oC] Specific heat capacity of moist air
    'Mcr_stat'                      : 34.,                # [-] Critical static slope for avalanching
    'Mcr_dyn'                       : 33.,                # [-] Critical dynamic slope for avalanching
    'fc'                            : 0.11,               # [-] Moisture content at field capacity (volumetric) (Schmutz, 2014)
    'resw_moist'                    : 0.01,              # [-] Residual soil moisture content (volumetric) (Schmutz, 2014)
    'satw_moist'                    : 0.35,               # [-] Satiated soil moisture content (volumetric) (Schmutz, 2014)
    'resd_moist'                    : 0.01,              # [-] Residual soil moisture content (volumetric) (Schmutz, 2014)
    'satd_moist'                    : 0.5,               # [-] Satiated soil moisture content (volumetric) (Schmutz, 2014)
    'nw_moist'                      : 2.3,                # [-] Pore-size distribution index in the soil water retention function (Schmutz, 2014)
    'nd_moist'                      : 4.5,                # [-] Pore-size distribution index in the soil water retention function (Schmutz, 2014)
    'mw_moist'                      : 0.57,                # [-] 1-1/n
    'md_moist'                      : 0.42,                # [-] 1-1/n
    'alfaw_moist'                   : -0.070,              # [cm^-1] Inverse of the air-entry value for a wetting branch of the soil water retention function (Schmutz, 2014)
    'alfad_moist'                   : -0.035,              # [cm^-1] Inverse of the air-entry value for a drying branch of the soil water retention function (Schmutz, 2014)
    'thick_moist'                   : 0.002,              # [m] Thickness of surface moisture soil layer
    'K_gw'                          : 0.00078,           # [m/s] Hydraulic conductivity (Schmutz, 2014)
    'ne_gw'                         : 0.3,                # [-] Effective porosity
    'D_gw'                          : 12,                  # [m] Aquifer depth
    'tfac_gw'                       : 10,                 # [-] Reduction factor for time step in ground water calculations
    'Cl_gw'                         : 270,              # [-] Runup infiltration coefficient
    'in_gw'                         : 1,                  # [m] Initial groundwater level
    'GW_stat'                       : 1,                  # [m] Landward static groundwater boundary (if static boundary is defined)                 
    'scheme'                        : 'euler_backward',   # Name of numerical scheme (euler_forward, euler_backward or crank_nicolson)
    'boundary_lateral'              : 'circular',         # Name of lateral boundary conditions (circular, noflux)
    'boundary_offshore'             : 'noflux',           # Name of offshore boundary conditions (gradient, noflux, constant, uniform)
    'boundary_offshore_flux'        : 0.,                 # Constant offshore boundary flux
    'boundary_onshore'              : 'gradient',         # Name of onshore boundary conditions (gradient, noflux, constant, uniform)
    'boundary_onshore_flux'         : 0.,                 # Constant onshore boundary flux
    'boundary_gw'                   : 'no_flow',          # Landward groundwater boundary, dGw/dx = 0 (or 'static')
    'method_moist_threshold'        : 'belly_johnson',    # Name of method to compute wind velocity threshold based on soil moisture content
    'method_moist_process'          : 'infiltration',     # Name of method to compute soil moisture content(infiltration or surface_moisture)
    'method_transport'              : 'bagnold',          # Name of method to compute equilibrium sediment transport rate
    'max_error'                     : 1e-6,               # [-] Maximum error at which to quit iterative solution in implicit numerical schemes
    'max_iter'                      : 1000,               # [-] Maximum number of iterations at which to quit iterative solution in implicit numerical schemes
    'refdate'                       : '1970-01-01 00:00', # [-] Reference datetime in netCDF output
    'callback'                      : None,               # Reference to callback function (e.g. example/callback.py':callback)
    'wind_convention'               : 'cartesian',        # Convention used for the wind direction in the input files
}



#: Merge initial and model state
MODEL_STATE.update({
    (k, MODEL_STATE[k] + INITIAL_STATE[k])
    for k in set(MODEL_STATE).intersection(INITIAL_STATE)
})
