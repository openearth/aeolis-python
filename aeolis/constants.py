import numpy as np


#: AeoLiS model default configuration
DEFAULT_CONFIG = dict(
    
th_grainsize    = True,               # Enable wind velocity threshold based on grainsize
th_bedslope     = False,              # Enable wind velocity threshold baded on bed slopes
th_moisture     = True,               # Enable wind velocity threshold baded on soil moisture content
th_humidity     = False,              # Enable wind velocity threshold baded on air humidity
th_roughness    = False,              # Enable wind velocity threshold baded on bed roughness elements
mixtoplayer     = True,               # Enable mixing of top bed layer by waves
bedupdate       = True,               # Enable bed update
evaporation     = True,               # Enable simulation of evaporation
xgrid_file      = None,               # Filename of ASCII file with x-coordinates of grid cells
ygrid_file      = None,               # Filename of ASCII file with y-coordinates of grid cells
bed_file        = None,               # Filename of ASCII file with bed level heights of grid cells
wind_file       = None,               # Filename of ASCII file with time series of wind velocity and direction
tide_file       = None,               # Filename of ASCII file with time series of water levels
wave_file       = None,               # Filename of ASCII file with time series of wave heights
meteo_file      = None,               # Filename of ASCII file with time series of meteorlogical conditions
bedcomp_file    = None,               # Filename of ASCII file with initial bed composition
nx              = 0,                  # [-] Number of grid cells in x-dimension
ny              = 0,                  # [-] Number of grid cells in y-dimension
dt              = 60.,                # [s] Time step size
tstart          = 0.,                 # [s] Start time of simulation
tstop           = 3600.,              # [s] End time of simulation
output_times    = 60.,                # [s] Output interval in seconds of simulation time
output_file     = 'aeolis.nc',        # Filename of netCDF4 output file
output_vars     = ['zb', 'zs',
                   'Ct', 'Cu',
                   'uw', 'uth',
                   'mass', 'pickup'], # Names of spatial grids to be included in output
output_types    = [],                 # Names of statistical parameters to be included in output (avg, sum, var, min or max)
grain_size      = [225e-6],           # [m] Average grain size of each sediment fraction
grain_dist      = [1.],               # [-] Initial distribution of sediment fractions
nfractions      = 1,                  # [-] Number of sediment fractions
nlayers         = 3,                  # [-] Number of bed layers
layer_thickness = .01,                # [m] Thickness of bed layers
g               = 9.81,               # [m/s^2] Gravitational constant
rhoa            = 1.25,               # [kg/m^3] Air density
rhop            = 2650.,              # [kg/m^3] Grain density
rhow            = 1025.,              # [kg/m^3] Water density
porosity        = .4,                 # [-] Sediment porosity
A               = 100.,               # [-] Constant in formulation for wind velocity threshold based on grain size
z0              = 1.,                 # [m] Measurement height of wind velocity
k               = 0.01,               # [m] Bed roughness
Cb              = 1.5,                # [-] Constant in formulation for equilibrium sediment concentration
bi              = 1.,                 # [-] Bed interaction factor
T               = 1.,                 # [s] Adaptation time scale in advection equation
F               = 1e-4,               # [-] Soil drying rate
eps             = 1e-3,               # [m] Minimum water depth to consider a cell "flooded"
gamma           = .5,                 # [-] Maximum wave height over depth ratio
facDOD          = .1,                 # [-] Ratio between depth of disturbance and local wave height
scheme          = 'crank_nicolson',   # Name of numerical scheme (euler_forward, euler_backward or crank_nicolson)
method_moist    = 'belly_johnson',    # Name of method to compute wind velocity threshold based on soil moisture content
max_error       = 1e-6,               # [-] Maximum error at which to quit iterative solution in implicit numerical schemes
max_iter        = 1000,               # [-] Maximum number of iterations at which to quit iterative solution in implicit numerical schemes

)


#: Required model configuration parameters
REQUIRED_CONFIG = ['nx', 'ny']
