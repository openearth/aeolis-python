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

# Question; what are the different purposes of INITIAL_STATE and MODEL_STATE?

#: Aeolis model state variables
INITIAL_STATE = {
    ('ny', 'nx') : (
        'uw',                               # [m/s] Wind velocity
        'uws',                              # [m/s] Component of wind velocity in x-direction
        'uwn',                              # [m/s] Component of wind velocity in y-direction
        
        'tau',                              # [N/m^2] Wind shear stress
        'taus',                             # [N/m^2] Component of wind shear stress in x-direction
        'taun',                             # [N/m^2] Component of wind shear stress in y-direction
        'tau0',                             # [N/m^2] Wind shear stress over a flat bed
        'taus0',                            # [N/m^2] Component of wind shear stress in x-direction over a flat bed
        'taun0',                            # [N/m^2] Component of wind shear stress in y-direction over a flat bed
        'taus_u',                           # [N/m^2] Saved direction of wind shear stress in x-direction
        'taun_u',                           # [N/m^2] Saved direction of wind shear stress in y-direction
        'dtaus',                            # [-] Component of the wind shear perturbation in x-direction
        'dtaun',                            # [-] Component of the wind shear perturbation in y-direction

        'ustar',                            # [m/s] Wind shear velocity
        'ustars',                           # [m/s] Component of wind shear velocity in x-direction
        'ustarn',                           # [m/s] Component of wind shear velocity in y-direction
        'ustar0',                           # [m/s] Wind shear velocity over a flat bed
        'ustars0',                          # [m/s] Component of wind shear velocity in x-direction over a flat bed
        'ustarn0',                          # [m/s] Component of wind shear velocity in y-direction over a flat bed

        'udir',                             # [rad] Wind direction
        'zs',                               # [m] Water level above reference (or equal to zb if zb > zs)
        'SWL',                              # [m] Still water level above reference
        'Hs',                               # [m] Wave height
        'Hsmix',                            # [m] Wave height for mixing (including setup, TWL)
        'Tp',                               # [s] Wave period for wave runup calculations
        'zne',                              # [m] Non-erodible layer
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
        'zs',                               # [m] Water level above reference
        'zne',                              # [m] Height above reference of the non-erodible layer
        'zb0',                              # [m] Initial bed level above reference
        'zdry',                             # [m]
        'dzdry',                            # [m]
        'dzb',                              # [m/dt] Bed level change per time step (computed after avalanching!)
        'dzbyear',                          # [m/yr] Bed level change translated to m/y
        'dzbavg',                           # [m/year] Bed level change averaged over collected time steps
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
        'd_hdelta',                         # [-] Moisture content (volumetric) computed on the main drying curve for hdelta
        'ustar',                            # [m/s] Shear velocity by wind
        'ustars',                           # [m/s] Component of shear velocity in x-direction by wind
        'ustarn',                           # [m/s] Component of shear velocity in y-direction by wind
        'ustar0',                           # [m/s] Initial shear velocity (without perturbation)
        'zsep',                             # [m] Z level of polynomial that defines the separation bubble
        'hsep',                             # [m] Height of separation bubble = difference between z-level of zsep and of the bed level zb
        'theta_dyn',                        # [degrees] spatially varying dynamic angle of repose for avalanching
        'rhoveg',                           # [-] Vegetation cover
        'drhoveg',                          # Change in vegetation cover
        'hveg',                             # [m] height of vegetation
        'dhveg',                            # [m] Difference in vegetation height per time step
        'dzbveg',                           # [m] Bed level change used for calculation of vegetation growth
        'germinate',                        # [bool] Newly vegetated due to germination (or establishment) 
        'lateral',                          # [bool] Newly vegetated due to lateral propagation 
        'vegetated',                        # [bool] Vegetated, determines if vegetation growth or burial is allowed
        'vegfac',                           # Vegetation factor to modify shear stress by according to Raupach 1993
        'fence_height',                     # Fence height
        'R',                                # [m] wave runup
        'eta',                              # [m] wave setup
        'sigma_s',                          # [m] swash
        'TWL',                              # [m] Total Water Level above reference (SWL + Run-up)
        'SWL',                              # [m] Still Water Level above reference
        'DSWL',                             # [m] Dynamic Still water level above reference (SWL + Set-up)
        'Rti',                              # [-] Factor taking into account sheltering by roughness elements
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
        'usST',                            # [NEW] [m/s] Component of the saltation velocity in x-direction for SedTRAILS
        'unST',                            # [NEW] [m/s] Component of the saltation velocity in y-direction for SedTRAILS
        'u0',
        'masstop',                          # [kg/m^2] Sediment mass in bed toplayer, only stored for output
    ),
    ('ny','nx','nlayers') : (
        'thlyr',                            # [m] Bed composition layer thickness
        'salt',                             # [-] Salt content
    ),
    ('ny','nx','nlayers','nfractions') : (
        'mass',                             # [kg/m^2] Sediment mass in bed
    ),
}


#: AeoLiS model default configuration
DEFAULT_CONFIG = {

    # --- Grid files (convention *.grd) --- #
    'xgrid_file'                    : None,               # Filename of ASCII file with x-coordinates of grid cells
    'ygrid_file'                    : None,               # Filename of ASCII file with y-coordinates of grid cells
    'bed_file'                      : None,               # Filename of ASCII file with bed level heights of grid cells
    'ne_file'                       : None,               # Filename of ASCII file with non-erodible layer
    'veg_file'                      : None,               # Filename of ASCII file with initial vegetation density

    # --- Model, grid and time settings --- #
    'tstart'                        : 0.,                 # [s] Start time of simulation
    'tstop'                         : 3600.,              # [s] End time of simulation
    'dt'                            : 60.,                # [s] Time step size
    'restart'                       : None,               # [s] Interval for which to write restart files
    'refdate'                       : '2020-01-01 00:00', # [-] Reference datetime in netCDF output
    'callback'                      : None,               # Reference to callback function (e.g. example/callback.py':callback)
    'wind_convention'               : 'nautical',         # Convention used for the wind direction in the input files (cartesian or nautical)
    'alfa'                          : 0,                  # [deg] Real-world grid cell orientation wrt the North (clockwise)
    'nx'                            : 0,                  # [-] Number of grid cells in x-dimension
    'ny'                            : 0,                  # [-] Number of grid cells in y-dimension

    # --- Input Timeseries --- #
    'wind_file'                     : None,               # Filename of ASCII file with time series of wind velocity and direction
    'tide_file'                     : None,               # Filename of ASCII file with time series of water levels
    'wave_file'                     : None,               # Filename of ASCII file with time series of wave heights
    'meteo_file'                    : None,               # Filename of ASCII file with time series of meteorlogical conditions

    # --- Boundary conditions --- #
    'boundary_lateral'              : 'flux',             # Name of lateral boundary conditions (circular, flux or constant)
    'boundary_offshore'             : 'flux',             # Name of offshore boundary conditions (circular, flux or constant)
    'boundary_onshore'              : 'flux',             # Name of onshore boundary conditions (circular, flux or constant)
    'offshore_flux'                 : 1.,                 # [-] Factor to determine offshore boundary flux as a function of Cu (= 1 for saturated, = 0 for noflux)
    'onshore_flux'                  : 1.,                 # [-] Factor to determine onshore boundary flux as a function of Cu (= 1 for saturated, = 0 for noflux)
    'lateral_flux'                  : 1.,                 # [-] Factor to determine lateral boundary flux as a function of Cu (= 1 for saturated, = 0 for noflux)

    # --- Sediment fractions and layers --- #
    'grain_size'                    : [225e-6],           # [m] Average grain size of each sediment fraction
    'grain_dist'                    : [1.],               # [-] Initial distribution of sediment fractions
    'nlayers'                       : 3,                  # [-] Number of bed layers
    'layer_thickness'               : .01,                # [m] Thickness of bed layers

    # --- Output (and coupling) settings --- #
    'visualization'                 : False,              # Boolean for visualization of model interpretation before and just after initialization
    'output_times'                  : 60.,                # [s] Output interval in seconds of simulation time
    'output_file'                   : None,               # Filename of netCDF4 output file
    'output_types'                  : None,               # Names of statistical parameters to be included in output (avg, sum, var, min or max)
    'output_vars'                   : ['zb', 'zs',
                                       'Ct', 'Cu',
                                       'uw', 'udir', 
                                       'uth', 'mass'
                                       'pickup', 'w'],    # Names of spatial grids to be included in output
    'external_vars'                 : None,               # Names of variables that are overwritten by an external (coupling) model, i.e. CoCoNuT
    'output_sedtrails'              : False,              # NEW! [T/F] Boolean to see whether additional output for SedTRAILS should be generated
    'nfraction_sedtrails'           : 0,                  # [-] Index of selected fraction for SedTRAILS (0 if only one fraction)
    
    # --- Process Booleans (True/False) --- #
    'process_wind'                  : True,               # Enable the process of wind
    'process_transport'             : True,               # Enable the process of transport
    'process_bedupdate'             : True,               # Enable the process of bed updating
    'process_threshold'             : True,               # Enable the process of threshold
    'process_avalanche'             : False,              # Enable the process of avalanching
    'process_shear'                 : False,              # Enable the process of wind shear
    'process_tide'                  : False,              # Enable the process of tides
    'process_wave'                  : False,              # Enable the process of waves
    'process_runup'                 : False,              # Enable the process of wave runup
    'process_moist'                 : False,              # Enable the process of moist
    'process_mixtoplayer'           : False,              # Enable the process of mixing 
    'process_wet_bed_reset'         : False,              # Enable the process of bed-reset in the intertidal zone
    'process_meteo'                 : False,              # Enable the process of meteo
    'process_salt'                  : False,              # Enable the process of salt
    'process_humidity'              : False,              # Enable the process of humidity
    'process_groundwater'           : False,              # Enable the process of groundwater
    'process_scanning'              : False,              # Enable the process of scanning curves
    'process_inertia'               : False,              # NEW
    'process_separation'            : False,              # Enable the including of separation bubble
    'process_vegetation'            : False,              # Enable the process of vegetation
    'process_vegetation_leeside'    : False,              # Enable the process of leeside vegetation effects on shear stress
    'process_fences'                : False,              # Enable the process of sand fencing
    'process_dune_erosion'          : False,              # Enable the process of wave-driven dune erosion
    'process_seepage_face'          : False,              # Enable the process of groundwater seepage (NB. only applicable to positive beach slopes)
    'process_bedinteraction'        : False,              # Enable the process of bed interaction in the advection equation
    
    # --- Threshold Booleans (True/False) --- #
    'th_grainsize'                  : True,               # Enable wind velocity threshold based on grainsize
    'th_bedslope'                   : False,              # Enable wind velocity threshold based on bedslope
    'th_moisture'                   : False,              # Enable wind velocity threshold based on moisture
    'th_drylayer'                   : False,              # Enable threshold based on drying of layer
    'th_humidity'                   : False,              # Enable wind velocity threshold based on humidity
    'th_salt'                       : False,              # Enable wind velocity threshold based on salt
    'th_sheltering'                 : False,              # Enable wind velocity threshold based on sheltering by roughness elements
    'th_nelayer'                    : False,              # Enable wind velocity threshold based on a non-erodible layer
    
    # --- Other spatial files / masks --- #
    'bedcomp_file'                  : None,               # Filename of ASCII file with initial bed composition
    'threshold_file'                : None,               # Filename of ASCII file with shear velocity threshold
    'fence_file'                    : None,               # Filename of ASCII file with sand fence location/height (above the bed)
    'supply_file'                   : None,               # Filename of ASCII file with a manual definition of sediment supply (mainly used in academic cases)
    'wave_mask'                     : None,               # Filename of ASCII file with mask for wave height
    'tide_mask'                     : None,               # Filename of ASCII file with mask for tidal elevation
    'runup_mask'                    : None,               # Filename of ASCII file with mask for run-up
    'threshold_mask'                : None,               # Filename of ASCII file with mask for the shear velocity threshold
    'gw_mask'                       : None,               # Filename of ASCII file with mask for the groundwater level
    'vver_mask'                     : None,               # Filename of ASCII file with mask for the vertical vegetation growth   

    # --- Nummerical Solver --- #
    'T'                             : 1.,                 # [s] Adaptation time scale in advection equation
    'CFL'                           : 1.,                 # [-] CFL number to determine time step in explicit scheme
    'accfac'                        : 1.,                 # [-] Numerical acceleration factor
    'max_bedlevel_change'           : 999.,               # [m] Maximum bedlevel change after one timestep. Next timestep dt will be modified (use 999. if not used)
    'max_error'                     : 1e-8,               # [-] Maximum error at which to quit iterative solution in implicit numerical schemes
    'max_iter'                      : 1000,               # [-] Maximum number of iterations at which to quit iterative solution in implicit numerical schemes
    'solver'                        : 'steadystate',      # Name of the solver (steadystate, euler_backward, euler_forward)

    # --- General physical constants and model parameters --- #
    'method_roughness'              : 'constant',         # Name of method to compute the roughness height z0, note that here the z0 = k
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
    'kappa'                         : 0.41,               # [-] Von Kármán constant

    # --- Shear / Perturbation / Topographic steering --- #
    'method_shear'                  : 'fft',              # Name of method to compute topographic effects on wind shear stress (fft, quasi2d, duna2d (experimental))
    'dx'                            : 1.,
    'dy'                            : 1.,
    'L'                             : 100.,               # [m] Typical length scale of dune feature (perturbation)
    'l'                             : 10.,                # [m] Inner layer height (perturbation)

    # --- Flow separation bubble (OLD) --- #
    'buffer_width'                  : 10,                 # [m] Width of the bufferzone around the rotational grid for wind perturbation
    'sep_filter_iterations'         : 0,                  # [-] Number of filtering iterations on the sep-bubble (0 = no filtering)
    'zsep_y_filter'                 : False,              # [-] Boolean for turning on/off the filtering of the separation bubble in y-direction

    # --- Sediment transport formulations --- #
    'method_transport'              : 'bagnold',          # Name of method to compute equilibrium sediment transport rate
    'method_grainspeed'             : 'windspeed',        # Name of method to assume/compute grainspeed (windspeed, duran, constant)
    'Cb'                            : 1.5,                # [-] Constant in bagnold formulation for equilibrium sediment concentration
    'Ck'                            : 2.78,               # [-] Constant in kawamura formulation for equilibrium sediment concentration
    'Cl'                            : 6.7,                # [-] Constant in lettau formulation for equilibrium sediment concentration
    'Cdk'                           : 5.,                 # [-] Constant in DK formulation for equilibrium sediment concentration
    'sigma'                         : 4.2,                # [-] Ratio between basal area and frontal area of roughness elements
    'beta'                          : 130.,               # [-] Ratio between drag coefficient of roughness elements and bare surface
    'bi'                            : 1.,                 # [-] Bed interaction factor for sediment fractions

    # --- Bed update parameters --- #
    'Tbedreset'                     : 86400.,             # [s] 
    
    # --- Moisture parameters --- #
    'method_moist_threshold'        : 'belly_johnson',    # Name of method to compute wind velocity threshold based on soil moisture content
    'method_moist_process'          : 'infiltration',     # Name of method to compute soil moisture content(infiltration or surface_moisture)
    'Tdry'                          : 3600.*1.5,          # [s] Adaptation time scale for soil drying

    # --- Moisture / Groundwater (Hallin) --- #
    'boundary_gw'                   : 'no_flow',          # Landward groundwater boundary, dGw/dx = 0 (or 'static')
    'fc'                            : 0.11,               # [-] Moisture content at field capacity (volumetric)
    'w1_5'                          : 0.02,               # [-] Moisture content at wilting point (gravimetric)
    'resw_moist'                    : 0.01,               # [-] Residual soil moisture content (volumetric) 
    'satw_moist'                    : 0.35,               # [-] Satiated soil moisture content (volumetric)
    'resd_moist'                    : 0.01,               # [-] Residual soil moisture content (volumetric) 
    'satd_moist'                    : 0.5,                # [-] Satiated soil moisture content (volumetric) 
    'nw_moist'                      : 2.3,                # [-] Pore-size distribution index in the soil water retention function
    'nd_moist'                      : 4.5,                # [-] Pore-size distribution index in the soil water retention function 
    'mw_moist'                      : 0.57,               # [-] m, van Genucthen param (can be approximated as 1-1/n)
    'md_moist'                      : 0.42,               # [-] m, van Genucthen param (can be approximated as 1-1/n)
    'alfaw_moist'                   : -0.070,             # [cm^-1] Inverse of the air-entry value for a wetting branch of the soil water retention function (Schmutz, 2014)
    'alfad_moist'                   : -0.035,             # [cm^-1] Inverse of the air-entry value for a drying branch of the soil water retention function (Schmutz, 2014)
    'thick_moist'                   : 0.002,              # [m] Thickness of surface moisture soil layer
    'K_gw'                          : 0.00078,            # [m/s] Hydraulic conductivity (Schmutz, 2014)
    'ne_gw'                         : 0.3,                # [-] Effective porosity
    'D_gw'                          : 12,                 # [m] Aquifer depth
    'tfac_gw'                       : 10,                 # [-] Reduction factor for time step in ground water calculations
    'Cl_gw'                         : 0.7,                # [m] Groundwater overheight due to runup
    'in_gw'                         : 0,                  # [m] Initial groundwater level
    'GW_stat'                       : 1,                  # [m] Landward static groundwater boundary (if static boundary is defined)
    'max_moist'                     : 10.,           # NEWCH      # [%] Moisture content (volumetric in percent) above which the threshold shear velocity is set to infinity (no transport, default value Delgado-Fernandez, 2010)
    'max_moist'                     : 10.,                # [%] Moisture content (volumetric in percent) above which the threshold shear velocity is set to infinity (no transport, default value Delgado-Fernandez, 2010)
    
    # --- Avalanching parameters --- #
    'theta_dyn'                     : 33.,                # [degrees] Initial Dynamic angle of repose, critical dynamic slope for avalanching
    'theta_stat'                    : 34.,                # [degrees] Initial Static angle of repose, critical static slope for avalanching
    'max_iter_ava'                  : 1000,               # [-] Maximum number of iterations at which to quit iterative solution in avalanching calculation

    # --- Hydro and waves --- #
    'eps'                           : 1e-3,               # [m] Minimum water depth to consider a cell "flooded"
    'gamma'                         : .5,                 # [-] Maximum wave height over depth ratio
    'xi'                            : .3,                 # [-] Surf similarity parameter
    'facDOD'                        : .1,                 # [-] Ratio between depth of disturbance and local wave height

    # --- Vegetation (OLD) --- #
    'method_vegetation'             : 'duran',            # Name of method to compute vegetation: duran (original) or grass (new framework)
    'avg_time'                      : 86400.,             # [s] Indication of the time period over which the bed level change is averaged for vegetation growth
    'gamma_vegshear'                : 16.,                # [-] Roughness factor for the shear stress reduction by vegetation
    'hveg_max'                      : 1.,                 # [m] Max height of vegetation
    'dzb_opt'                       : 0.,                 # [m/year] Sediment burial for optimal growth
    'V_ver'                         : 0.,                 # [m/year] Vertical growth potential
    'germinate'                     : 0.,                 # [1/year] Possibility of germination per year
    'lateral'                       : 0.,                 # [1/year] Posibility of lateral expension per year
    'veg_gamma'                     : 1.,                 # [-] Constant on influence of sediment burial
    'veg_sigma'                     : 0.,                   # [-] Sigma in gaussian distrubtion of vegetation cover filter
    'vegshear_type'                 : 'raupach',          # Choose the Raupach grid based solver (1D or 2D) or the Okin approach (1D only)
    'okin_c1_veg'                   : 0.48,               #x/h spatial reduction factor in Okin model for use with vegetation
    'okin_c1_fence'                 : 0.48,               #x/h spatial reduction factor in Okin model for use with sand fence module
    'okin_initialred_veg'           : 0.32,               #initial shear reduction factor in Okin model for use with vegetation
    'okin_initialred_fence'         : 0.32,               #initial shear reduction factor in Okin model for use with sand fence module
    'veggrowth_type'                : 'orig',             #'orig', 'duranmoore14'
    'rhoveg_max'                    : 0.5,                #maximum vegetation density, only used in duran and moore 14 formulation
    't_veg'                         : 3,                  #time scale of vegetation growth (days), only used in duran and moore 14 formulation
    'v_gam'                         : 1,                  # only used in duran and moore 14 formulation

    # --- Dune erosion parameters --- #
    'dune_toe_elevation'            : 3,                  # Choose dune toe elevation, only used in the PH12 dune erosion solver
    'beach_slope'                   : 0.1,                # Define the beach slope, only used in the PH12 dune erosion solver
    'veg_min_elevation'             : -10.,               # Minimum elevation (m) where vegetation can grow; default -10 disables restriction.

    # --- Bed interaction in advection equation (new process) --- #
    'zeta_base'                     : 1.0,                # [-] Base value for bed interaction parameter in advection equation 
    'zeta_sheltering'               : False,              # [-] Include sheltering effect of roughness elements on bed interaction parameter
    'p_zeta_moist'                  : 0.8,                # [-] Exponent parameter for computing zeta from moisture
    'a_weibull'                     : 1.0,                # [-] Shape parameter k of Weibull function for bed interaction parameter zeta
    'b_weibull'                     : 0.5,                # [m] Scale parameter lambda of Weibull function for bed interaction parameter zeta
    'bounce'                        : [0.75],              # [-] Fraction of sediment skimming over vegetation canopy (species-specific)
    'alpha_lift'                    : 0.2,                # [-] Vegetation-induced upward lift (0-1) of transport-layer centroid

    # --- Grass vegetation model (new vegetation framework) --- #
    'method_vegetation'             : 'duran',       # ['duran' | 'grass'] Vegetation formulation
    'veg_res_factor'                : 5,             # [-] Vegetation subgrid refinement factor (dx_veg = dx / factor)
    'dt_veg'                        : 86400.,        # [s] Time step for vegetation growth calculations
    'species_names'                 : ['marram'],    # [-] Name(s) of vegetation species
    'hveg_file'                     : None,          # Filename of ASCII file with initial vegetation height (shape: ny * nx * nspecies)
    'Nt_file'                       : None,          # Filename of ASCII file with initial tiller density (shape: ny * nx * nspecies)

    'd_tiller'                      : [0.006],       # [m] Mean tiller diameter
    'r_stem'                        : [0.2],         # [-] Fraction of rigid (non-bending) stem height
    'alpha_uw'                      : [-0.0412],     # [s/m] Wind-speed sensitivity of vegetation bending
    'alpha_Nt'                      : [1.95e-4],     # [m^2] Tiller-density sensitivity of vegetation bending
    'alpha_0'                       : [0.9445],      # [-] Baseline bending factor (no wind, sparse vegetation)

    'G_h'                           : [1.0],         # [m/yr] Intrinsic vertical vegetation growth rate
    'G_c'                           : [2.5],         # [tillers/tiller/yr] Intrinsic clonal tiller production rate
    'G_s'                           : [0.01],        # [tillers/tiller/yr] Intrinsic seedling establishment rate
    'Hveg'                          : [0.8],         # [m] Maximum attainable vegetation height
    'phi_h'                         : [1.0],         # [-] Saturation exponent for height growth

    'Nt_max'                        : [900.0],       # [1/m^2] Maximum attainable tiller density
    'R_cov'                         : [1.2],         # [m] Radius for neighbourhood density averaging

    'lmax_c'                        : [0.9],         # [m] Maximum clonal dispersal distance
    'mu_c'                          : [2.5],         # [-] Shape parameter of clonal dispersal kernel
    'alpha_s'                       : [4.0],         # [m^2] Scale parameter of seed dispersal kernel
    'nu_s'                          : [2.5],         # [-] Tail-heaviness of seed dispersal kernel

    'T_burial'                      : 86400.*30.,    # [s] Time scale for sediment burial effect on vegetation growth (replaces avg_time)
    'gamma_h'                       : [1.0],         # [-] Sensitivity of vertical growth to burial (1 / dzb_tol_h)
    'dzb_tol_c'                     : [1.0],         # [m/yr] Tolerance burial range for clonal expansion
    'dzb_tol_s'                     : [0.1],         # [m/yr] Tolerance burial range for seed establishment
    'dzb_opt_h'                     : [0.5],         # [m/yr] Optimal burial rate for vertical growth
    'dzb_opt_c'                     : [0.5],         # [m/yr] Optimal burial rate for clonal expansion
    'dzb_opt_s'                     : [0.025],       # [m/yr] Optimal burial rate for seed establishment

    'beta_veg'                      : [120.0],       # [-] Vegetation momentum-extraction efficiency (Raupach)
    'm_veg'                         : [0.4],         # [-] Shear non-uniformity correction factor
    'c1_okin'                       : [0.48],        # [-] Downwind decay coefficient in Okin shear reduction

    'veg_sigma'                     : 0.,            # [-] Sigma in gaussian distrubtion of vegetation cover filter
    'zeta_sigma'                    : 0.,            # [-] Standard deviation for smoothing vegetation bed interaction parameter

    'alpha_comp'                    : [0.],          # [-] Lotka–Volterra competition coefficients
                                                     #      shape: nspecies * nspecies (flattened)
                                                     #      alpha_comp[k,l] = effect of species l on species k
                                    
    'T_flood'                       : 7200.,         # [s] Time scale for vegetation flood stress mortality (half-life under constant inundation)
    'gamma_Nt_decay'                : 0.,            # [-] Sensitivity of tiller density decay to relative reduction in hveg


    # --- Separation bubble parameters --- #
    'sep_look_dist'                 : 50.,           # [m] Flow separation: Look-ahead distance for upward curvature anticipation
    'sep_k_press_up'                : 0.05,          # [-] Flow separation: Press-up curvature 
    'sep_k_crit_down'               : 0.18,          # [1/m] Flow separation: Maximum downward curvature
    'sep_s_crit'                    : 0.18,          # [-] Flow separation: Critical bed slope below which reattachment is forced
    'sep_s_leeside'                 : 0.25,          # [-] Maximum downward leeside slope of the streamline

    # --- Other --- #
    'Tsalt'                         : 3600.*24.*30.,      # [s] Adaptation time scale for salinitation
    'csalt'                         : 35e-3,              # [-] Maximum salt concentration in bed surface layer
    'cpair'                         : 1.0035e-3,          # [MJ/kg/oC] Specific heat capacity air
}

REQUIRED_CONFIG = ['nx', 'ny']

#: Merge initial and model state
MODEL_STATE.update({
    (k, MODEL_STATE[k] + INITIAL_STATE[k])
    for k in set(MODEL_STATE).intersection(INITIAL_STATE)
})
