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
        'tau',                              # [N/m^2] Wind shear stress
        'taus',                             # [N/m^2] Component of wind shear stress in x-direction
        'taun',                             # [N/m^2] Component of wind shear stress in y-direction
        'dtaus',                            # [-] Component of the wind shear perturbation in x-direction
        'dtaun',                            # [-] Component of the wind shear perturbation in y-direction
        'taus_u',                      #NEW # [N/m^2] Saved direction of wind shear stress in x-direction
        'taun_u',                      #NEW # [N/m^2] Saved direction of wind shear stress in y-direction
        'udir',                             # [rad] Wind direction
        'zs',                               # [m] Water level above reference
        'Hs',                               # [m] Wave height
        'zne',                         #NEW # [m] Non-erodible layer
    ),
}

MODEL_STATE = {
    ('ny', 'nx') : (
        'x',                                # [m] Real-world x-coordinate of grid cell center
        'y',                                # [m] Real-world y-coordinate of grid cell center
        'ds',                               # [m] Real-world grid cell size in x-direction
        'dn',                               # [m] Real-world grid cell size in y-direction
        'dsdn',                             # [m^2] Real-world grid cell surface area
        'dsdni',                            # [m^-2] Inverse of real-world grid cell surface area
#        'alfa',                             # [rad] Real-world grid cell orientation #Sierd_comm in later releases this needs a revision 
        'zb',                               # [m] Bed level above reference
        'zb0',                        # NEW # [m] Initial bed level above reference
        'zdry',                       # NEW # [m]
        'dzdry',                      # NEW # [m]
        'dzb',                        # NEW # [m/dt] Bed level change per time step (computed after avalanching!)
        'dzbyear',                    # NEW # [m/yer] Bed level change tranlated to m/y
        'dzbavg',                     # NEW # [m/year] Bed level change averaged over collected time steps
        'S',                                # [-] Level of saturation
        'ustar',                      # NEW # [m/s] Shear velocity by wind
        'ustars',                     # NEW # [m/s] Component of shear velocity in x-direction by wind
        'ustarn',                     # NEW # [m/s] Component of shear velocity in y-direction by wind
        'ustar0',                     # NEW # [m/s] Initial shear velocity (without perturbation)
        'zsep',                       # NEW # [m] Z level of polynomial that defines the separation bubble
        'hsep',                       # NEW # [m] Height of separation bubbel = difference between z-level of zsep and of the bed level zb
        'theta_stat',                 # NEW # [degrees] Updated, spatially varying static angle of repose
        'theta_dyn',                  # NEW # [degrees] Updated, spatially varying dynamic angle of repose
        'rhoveg',                     # NEW # [-] Vegetation cover
        'drhoveg',                    # NEW # Change in vegetation cover
        'hveg',                       # NEW # [m] height of vegetation
        'dhveg',                      # NEW # [m] Difference in vegetation height per time step
        'dzbveg',                     # NEW # [m] Bed level change used for calculation of vegetation growth
        'germinate',                  # NEW # 
        'lateral',                    # NEW #
        'dxrhoveg',                   # NEW #
        'vegfac',                     # NEW # [] Vegetation factor 
    ),
    ('ny','nx','nfractions') : (
        'Cu',                               # [kg/m^2] Equilibrium sediment concentration integrated over saltation height
        'Cuf',                              # [kg/m^2] Equilibrium sediment concentration integrated over saltation height, assuming the fluid shear velocity threshold
        'Cu0',                              # [kg/m^2] Flat bad equilibrium sediment concentration integrated over saltation height
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
        'uth0',                             # [m/s] Shear velocity threshold based on grainsize only (aerodynamic entrainment)
        'u',                                # [m/s] Mean horizontal saltation velocity in saturated state
        'us',                               # [m/s] Component of the saltation velocity in x-direction
        'un',                               # [m/s] Component of the saltation velocity in y-direction
        'u0',
    ),
    ('ny','nx','nlayers') : (
        'thlyr',                            # [m] Bed composition layer thickness
        'moist',                            # [-] Moisure content
        'salt',                             # [-] Salt content
    ),
    ('ny','nx','nlayers','nfractions') : (
        'mass',                             # [kg/m^2] Sediment mass in bed
    ),
}


#: AeoLiS model default configuration
DEFAULT_CONFIG = {
    'process_wind'                  : True,               # Enable the process of wind
    'process_shear'                 : True,               # Enable the process of wind shear
    'process_tide'                  : True,               # Enable the process of tides
    'process_wave'                  : True,               # Enable the process of waves
    'process_runup'                 : True,               # Enable the process of wave runup
    'process_moist'                 : True,               # Enable the process of moist
    'process_mixtoplayer'           : True,               # Enable the process of mixing
    'process_threshold'             : True,               # Enable the process of threshold
    'process_transport'             : True,               # Enable the process of transport
    'process_bedupdate'             : True,               # Enable the process of bed updating
    'process_meteo'                 : False,              # Enable the process of meteo
    'process_salt'                  : False,              # Enable the process of salt
    'process_humidity'              : False,              # Enable the process of humidity
    'process_avalanche'             : True,         # NEW # Enable the process of avalanching
    'process_inertia'               : False,        # NEW 
    'process_separation'            : True,         # NEW # Enable the including of separation bubble
    'process_nelayer'               : False,        # NEW # Enable a non-erodible layer
    'process_vegetation'            : False,        # NEW # Enable the process of vegetation 
    'th_grainsize'                  : True,               # Enable wind velocity threshold based on grainsize
    'th_bedslope'                   : False,              # Enable wind velocity threshold based on bedslope
    'th_moisture'                   : True,               # Enable wind velocity threshold based on moisture
    'th_drylayer'                   : False,        # NEW # Enable threshold based on drying of layer
    'th_humidity'                   : False,              # Enable wind velocity threshold based on humidity
    'th_salt'                       : False,              # Enable wind velocity threshold based on salt
    'th_roughness'                  : True,               # Enable wind velocity threshold based on roughness
    'th_nelayer'                    : False,        # NEW # Enable wind velocity threshold based on a non-erodible layer
    'xgrid_file'                    : None,               # Filename of ASCII file with x-coordinates of grid cells
    'ygrid_file'                    : None,               # Filename of ASCII file with y-coordinates of grid cells
    'bed_file'                      : None,               # Filename of ASCII file with bed level heights of grid cells
    'wind_file'                     : None,               # Filename of ASCII file with time series of wind velocity and direction
    'tide_file'                     : None,               # Filename of ASCII file with time series of water levels
    'wave_file'                     : None,               # Filename of ASCII file with time series of wave heights
    'meteo_file'                    : None,               # Filename of ASCII file with time series of meteorlogical conditions
    'bedcomp_file'                  : None,               # Filename of ASCII file with initial bed composition
    'threshold_file'                : None,               # Filename of ASCII file with shear velocity threshold
    'ne_file'                       : None,         # NEW # Filename of ASCII file with non-erodible layer
    'veg_file'                      : None,         # NEW # Filename of ASCII file with initial vegetation density
    'wave_mask'                     : None,               # Filename of ASCII file with mask for wave height
    'tide_mask'                     : None,               # Filename of ASCII file with mask for tidal elevation
    'threshold_mask'                : None,               # Filename of ASCII file with mask for the shear velocity threshold
    'nx'                            : 0,                  # [-] Number of grid cells in x-dimension
    'ny'                            : 0,                  # [-] Number of grid cells in y-dimension
    'dt'                            : 60.,                # [s] Time step size
    'dx'                            : 1.,
    'dy'                            : 1.,
    'CFL'                           : 1.,                 # [-] CFL number to determine time step in explicit scheme
    'accfac'                        : 1.,                 # [-] Numerical acceleration factor
    'tstart'                        : 0.,                 # [s] Start time of simulation
    'tstop'                         : 3600.,              # [s] End time of simulation
    'restart'                       : None,               # [s] Interval for which to write restart files
    'dzb_interval'                  : 86400,        # NEW # [s] Interval used for calcuation of vegetation growth
    'output_times'                  : 60.,                # [s] Output interval in seconds of simulation time
    'output_file'                   : None,               # Filename of netCDF4 output file
    'output_vars'                   : ['zb', 'zs',
                                       'Ct', 'Cu',
                                       'uw', 'udir', 
                                       'uth', 'mass'
                                       'pickup', 'w'],    # Names of spatial grids to be included in output
    'output_types'                  : [],                 # Names of statistical parameters to be included in output (avg, sum, var, min or max)
    'grain_size'                    : [225e-6],           # [m] Average grain size of each sediment fraction
    'grain_dist'                    : [1.],               # [-] Initial distribution of sediment fractions
    'nfractions'                    : 1,                  # [-] Number of sediment fractions
    'nlayers'                       : 3,                  # [-] Number of bed layers
    'layer_thickness'               : .01,                # [m] Thickness of bed layers
    'g'                             : 9.81,               # [m/s^2] Gravitational constant
    'v'                             : 0.000015,           # [m^2/s] Air viscosity  
    'rhoa'                          : 1.225,              # [kg/m^3] Air density
    'rhog'                          : 2650.,              # [kg/m^3] Grain density
    'rhow'                          : 1025.,              # [kg/m^3] Water density
    'porosity'                      : .4,                 # [-] Sediment porosity
    'Aa'                            : .085,               # [-] Constant in formulation for wind velocity threshold based on grain size
    'z'                             : 10.,                # [m] Measurement height of wind velocity
    'h'                             : None,               # [m] Representative height of saltation layer
    'k'                             : 0.001,              # [m] Bed roughness
    'L'                             : 100.,          #NEW # [m] Typical length scale of dune feature (perturbation)
    'l'                             : 10.,           #NEW # [m] Inner layer height (perturbation)
    'c_b'                           : 0.2,          # NEW # [-] Slope at the leeside of the separation bubble # c = 0.2 according to Durán 2010 (Sauermann 2001: c = 0.25 for 14 degrees)
    'mu_b'                          : 30,           # NEW # [deg] Minimum required slope for the start of flow separation
    'Cb'                            : 1.5,                # [-] Constant in bagnold formulation for equilibrium sediment concentration
    'Ck'                            : 2.78,               # [-] Constant in kawamura formulation for equilibrium sediment concentration
    'Cl'                            : 6.7,                # [-] Constant in lettau formulation for equilibrium sediment concentration
    'Cdk'                           : 5.,                 # [-] Constant in DK formulation for equilibrium sediment concentration
    'm'                             : 0.5,                # [-] Factor to account for difference between average and maximum shear stress
#    'alpha'                         : 0.4,                # [-] Relation of vertical component of ejection velocity and horizontal velocity difference between impact and ejection 
    'kappa'                         : 0.41,               # [-] Von Kármán constant
    'sigma'                         : 4.2,                # [-] Ratio between basal area and frontal area of roughness elements
    'beta'                          : 130.,               # [-] Ratio between drag coefficient of roughness elements and bare surface
    'bi'                            : 1.,                 # [-] Bed interaction factor
    'T'                             : 1.,                 # [s] Adaptation time scale in advection equation
    'Tdry'                          : 3600.*1.5,          # [s] Adaptation time scale for soil drying
    'Tsalt'                         : 3600.*24.*30.,      # [s] Adaptation time scale for salinitation
    'Tswash'                        : 30.,                # [s] 
    'eps'                           : 1e-3,               # [m] Minimum water depth to consider a cell "flooded"
    'gamma'                         : .5,                 # [-] Maximum wave height over depth ratio
    'xi'                            : .3,                 # [-] Surf similarity parameter
    'facDOD'                        : .1,                 # [-] Ratio between depth of disturbance and local wave height
    'csalt'                         : 35e-3,              # [-] Maximum salt concentration in bed surface layer
    'cpair'                         : 1.0035e-3,          # [MJ/kg/oC] Specific heat capacity air
    'theta_dyn'                     : 33.,          # NEW # [degrees] Initial Dynamic angle of repose, critical dynamic slope for avalanching 
    'theta_stat'                    : 34.,          # NEW # [degrees] Initial Static angle of repose, critical static slope for avalanching
    'avg_time'                      : 86400.,       # NEW # [s] Indication of the time period over which the bed level change is averaged for vegetation growth
    'gamma_vegshear'                : 16.,          # NEW # [-] Roughness factor for the shear stress reduction by vegetation
    'hveg_max'                      : 1.,           # NEW # [m] Max height of vegetation
    'dzb_opt'                       : 0.,           # NEW # [m/year] Sediment burial for optimal growth
    'V_ver'                         : 0.,           # NEW # [m/year] Vertical growth 
    'V_lat'                         : 0.,           # NEW # [m/year] Lateral growth
    'germinate'                     : 0.,           # NEW # [1/year] Possibility of germination per year
    'lateral'                       : 0.,           # NEW # [1/year] Posibility of lateral expension per year
    'veg_gamma'                     : 1.,           # NEW # [-] Constant on influence of sediment burial
    'veg_sigma'                     : 0.8,          # NEW # [-] Sigma in gaussian distrubtion of vegetation cover filter  
    'sedimentinput'                 : 0.,           # NEW # [-] Constant boundary sediment influx (only used in solve_pieter)
    'scheme'                        : 'euler_backward',   # Name of numerical scheme (euler_forward, euler_backward or crank_nicolson)
    'boundary_lateral'              : 'circular',         # Name of lateral boundary conditions (circular, constant ==noflux)
    'boundary_offshore'             : 'constant',         # Name of offshore boundary conditions (flux, constant, uniform, gradient)
    'boundary_onshore'              : 'gradient',         # Name of onshore boundary conditions (flux, constant, uniform, gradient)
    'offshore_flux'                 : 0.,           # NEW # [-] Factor to determine offshore boundary flux as a function of Q0 (= 1 for saturated flux , = 0 for noflux)
    'constant_offshore_flux'        : 0.,           # NEW # [kg/m/s] Constant input flux at offshore boundary
    'onshore_flux'                  : 0.,           # NEW # [-] Factor to determine onshore boundary flux as a function of Q0 (= 1 for saturated flux , = 0 for noflux)
    'constant_onshore_flux'         : 0.,           # NEW # [kg/m/s] Constant input flux at offshore boundary
    'lateral_flux'                  : 0.,           # NEW # [-] Factor to determine lateral boundary flux as a function of Q0 (= 1 for saturated flux , = 0 for noflux)
    'method_moist'                  : 'belly_johnson',    # Name of method to compute wind velocity threshold based on soil moisture content
    'method_transport'              : 'bagnold',          # Name of method to compute equilibrium sediment transport rate
    'max_error'                     : 1e-6,               # [-] Maximum error at which to quit iterative solution in implicit numerical schemes
    'max_iter'                      : 1000,               # [-] Maximum number of iterations at which to quit iterative solution in implicit numerical schemes
    'max_iter_ava'                  : 1000,               # [-] Maximum number of iterations at which to quit iterative solution in avalanching calculation
    'refdate'                       : '2020-01-01 00:00', # [-] Reference datetime in netCDF output
    'callback'                      : None,               # Reference to callback function (e.g. example/callback.py':callback)
    'wind_convention'               : 'nautical',         # Convention used for the wind direction in the input files
    'alfa'                         : 0,                   # [deg] Real-world grid cell orientation wrt the North (clockwise)
    'solver'                        : 'trunk',      # NEW # Choose the solver to be used (steadystate / trunk / pieter)
}

#: Required model configuration parameters
REQUIRED_CONFIG = ['nx', 'ny']

#: Merge initial and model state
MODEL_STATE.update({
    (k, MODEL_STATE[k] + INITIAL_STATE[k])
    for k in set(MODEL_STATE).intersection(INITIAL_STATE)
})
