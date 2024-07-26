Model setup
=================
INTRODUCTORY TEXT

Model input
-----------

The computational grid and boundary conditions for AeoLiS are specified through external
input files called by the model parameter file aeolis.txt.  The computational grid is defined
with an x grid, y grid, and z grid.  Boundary conditions for wind, wave, and tides
are also specified with external text files.  A list of additional grid and boundary
files can be found in the table below.  Each file is further defined below.

.. list-table:: 
   :widths: 25 50
   :header-rows: 1

   * - Input File
     - File Description
   * - aeolis.txt
     - File containing parameter definitions
   * - x.grd
     - File containing cross-shore grid
   * - y.grd
     - File containing alongshore grid (can be all zeros for 1D cases)
   * - z.grd
     - File containing topography and bathymetry data
   * - veg.grd
     - File containing initial vegetation density   
   * - mass.txt
     - File containing sediment mass data when using space varying grain size distribution  
   * - wind.txt
     - File containing wind speed and direction data
   * - tide.txt
     - File containing water elevation data
   * - wave.txt
     - File containing wave height and period data
   * - meteo.txt
     - File containing meteorological time series data

aeolis.txt
^^^^^^^^^^

This is the parameter file for AeoLiS that defines the model processes and boundary conditions.
Parameters in the file are specified by various keywords; each keyword has a pre-defined
default value that will be used if it is not direclty specified in aeolis.txt (a list of default parameter
values can be found in the Default settings tab on the left).  Among the keywords
in aeolis.txt are the keywords to define the external computational grid files (xgrid_file, ygrid_file,
and bed_file) and external boundary condition files (tide_file, wave_file, wind_file, etc.).  
The different physical processes in AeoLiS can be turned on and off by changing the 
process keywords in aeolis.txt to T (True) and F (False).  Example aeolis.txt parameters files can be
found in the examples folder on the AeoLiS GitHub.

x.grd
^^^^^

The x.grd file defines the computational grid in the cross-shore direction defined in meters.  
In a 1-dimensional (1D) case, the file contains a single column of cross-shore locations 
starting at zero for a location of choice.  In a 2-dimesional (2D) case, the file
contains multiple columns (cross-shore positions) and rows (alongshore positions) 
where each value corresponds to a specific location in the 2D grid.  The file can be renamed 
and is referenced from the parameters file with the xgrid_file keyword.

y.grd
^^^^^

This file defines the computational grid in the alongshore direction.  In a 1D case,
y.grd will contain a single column of zeros.  In a 2D case, similar to the x.grd file, 
y.grd has multiple columns (cross-shore positions) and rows (alongshore positions)
where each row, column position corresponds to a specific location in the 2D gird.
x.grd and y.grd will always be the same size regardless of running a 1D or 2D simulation.
As with the x.grd file, this file can be renamed and is referenced from the parameters file with the
keyword: ygrid_file.    

z.grd
^^^^^

The z.grd file provides the model with the elevation information for the computational 
grid defined in x.grd and y.grd.  Similar to x.grd and y.grd, when running 
AeoLis in 1D the file contains a single column with the number of rows equal 
to the number of rows in x.grd and y.grd.  In 2D cases, z.grd has multiple columns and 
rows of equal size to x.grd and y.grd.  Elevation values in the file should be defined such that
positive is up and negative is down.  The file can be renamed and is referenced from the 
parameters file with the keyword: bed_file.

veg.grd
^^^^^^^

The veg.grd file is an optional grid providing initial vegetation coverage (density) at each position in the 
model domain defined in x.grd and y.grd.  Similar to the grid files, if simulations are in 
2D there will be multiple columns for each cross-shore position (x) and multiple rows for 
each alongshore position (y).  The format of a 1D vegetation grid file can be seen below 
where each red dots represent vegetation cover at each cross-shore position. 

.. _fig-veg-inputs:

.. figure:: /images/vegetation_text_file.jpeg
   :alt: vegetation input format
   :width: 200px
   :align: center
   
   File format for a 1D AeoLis vegetation grid.  Each red dot is the vegetation density at a specific location in the computational grid.

mass.txt
^^^^^^^^

The mass.txt file allows users to specify variations in grain size distribution in both 
horizontal and vertical directions.  If the grain size distribution is constant throughout
the model domain, multifraction sediment transport is possilbe without this file.  The file contains
the mass of each sediment fraction in each grid cell and bed layer. The file is formatted such that each
row corresponds to a specific location in the computational domain and the columns are grouped 
by bed layers and each individual column represents a single sediment fraction present in the model
domain.  An infinite number of sediment fractions can be defined in the model; however, it should be 
noted the more sediment fractions present the longer the simulation time and larger the output files.    

In a 1D case, the text file will have dimensions of number of cross-shore locations (x) by number 
of sediment fractions times the number of bed layers.  For example if you have 200 cross-shore positions
in your model domain and 4 different sediment fractions with 3 bed layers, your mass.txt file 
will contain a matrix of 200 rows by 12 columns.  An example of a 1D mass.txt file can be seen below 
where each red dot represents a sediment fraction mass at a specific location in the model domain.

.. _fig-mass-inputs-1D:

.. figure:: /images/mass_text_file_graphic.jpeg
   :alt: mass file format 1D
   :width: 550px
   :align: center
   
   File format for a 1D AeoLis mass for spatially variable grain size distributions.  Each red dot is the mass for each sediment fraction
   at each location in the computational grid (x, y, bed layer).
 
In a 2D case, the mass.txt file will have dimensions of number of cross-shore positions (x)
times the number of alongshore positions (x) by number of sediment fractions times the number of
bed layers.  The file will be formatted such that the columns are grouped by bed layer with all available
sediment fractions present in each bed layer and rows are grouped by alongshore position with all 
cross-shore prositions given for each alongshore position.  An visual example of a 2D mass.txt input
file for AeoLis can be seen below.

.. _fig-mass-inputs-2D:

.. figure:: /images/mass_text_file_2D.jpeg
   :alt: mass file format 2D
   :width: 550px
   :align: center
   
   File format for a 2D AeoLis mass file for spatially variable grain size distributions.  Each red dot is the mass for each sediment fraction
   at each location in the computational grid (x, y, bed layer).

wind.txt
^^^^^^^^

The wind.txt file provides the model with wind boundary conditions and is formatted similar to 
the tide.txt and wave.txt files.  The first column is time in seconds from 
start, the second column is wind speed, and the third column is wind direction.  The wind directions
can be specified in either nautical or cartesian convention (specified in aeolis.txt with keyword: wind_convention).  
The format of this file can be seen below were each of the red dots represents a data value of time, wind speed, 
or wind direction.  As AeoLiS is an aeolian sediment transport model, the wind boundary conditions are of particular
importance.      

.. _fig-wind-inputs:

.. figure:: /images/wind_text_file_graphic.jpeg
   :alt: wind input format
   :width: 300px
   :align: center
   
   File format for wind boundary conditions file for AeoLis input.

tide.txt
^^^^^^^^

The tide.txt file contains the water elevation data for the duration of the 
simulation.  It is formatted such that the first column is time in seconds and 
the second column is the water elevation data at each time step.  An example of 
the file format can be seen below where each red dot represents a data value for 
time or water elevation.

.. _fig-tide-inputs:

.. figure:: /images/tide_text_file.jpeg
   :alt: tide input format
   :width: 300px
   :align: center
   
   File format for the water elevation conditions file for AeoLis input.
   
wave.txt
^^^^^^^^

The wave.txt file provides the model with wave data used in AeoLiS for runup calculations.  
The file is formatted similar to tide.txt but has three columns instead of two.  
Here, the first column is time in seconds, the second column is wave height, 
and the third column is the wave period.  The format of this file can be seen 
below where each red dot represents 
a data value. 

.. _fig-wave-inputs:

.. figure:: /images/wave_text_file_graphic.jpeg
   :alt: wave input format
   :width: 300px
   :align: center
   
   File format for the wave conditions file for AeoLis input.

meteo.txt
^^^^^^^^^

The meteo.txt file contains meteorological data used to simulate surface moisture in the model domain (see Simulation of surface moisture 
in Model description on for surface moisture implementation in AeoLiS).  This file is formatted similar to the other environmental boundary
condition files (wind, wave, and tide) such that it contains a time series of environmental data read into AeoLiS through keyword specification. 
The keywords required to process surface moisture with evaporation and infiltration are process_moist = True, method_moist_process = surf_moisture, 
th_moisture = True, and meteo_file = meteo.txt (or name of file containing meteorological data).  An example of the meteo.txt file can be seen in the 
figure below where each red dot represents a time series data value.  The first column contains time (s), the second column is temperature (degrees C),
the thrid column is precipitation (mm/hr), the fourth column is relative humidity (%), the fifth column is global radiation (MJ/$m^2$/day), and the sixth
column is air pressure (kPa).  

.. _fig-meteo-inputs:

.. figure:: /images/meteo_file_format.jpeg
   :alt: meteo file format
   :width: 550px
   :align: center
   
   File format for meteorological data used to simulate surface moisture in AeoLiS where each red dot represents a time series value. 

Default settings
----------------

The AeoLiS model can be configured using a model configuration
file. For any configuration parameters not defined in the model
configuration file, or in case the model configuration file is absent,
the default model configuration is used. The default model
configuration is listed below.

.. literalinclude::   ../../aeolis/constants.py
   :language: python
   :start-after: #: AeoLiS model default configuration
   :end-before: #: Merge initial and model state


Activate/deactivate processes
-----------------------------
After creating the input files that are necessary to run an AeoLiS model, the next step is often to decide which processes and methods to use. Several processes are defined in the configuration file that can be turned on and off. Apart from turning processes on and off, there are also several user-defined thresholds and methods that affect the way in which processes are calculated. For example, there are different sediment transport equations available within *process_transport*. The default is Bagnold, but by defining *method_transport* in the configuration file a different equation can be used. Here, we provide a description of the processes and methods that are defined in configuration file. More detailed descriptions of the processes and their implementation can be found in :ref:model_description

An easy way to look up where these process, threshold and method flags are used is by going to the main page of the AeoLiS github and using the search bar at the top. For instance, searching *process_tide* shows that it is used in :py:mod:`aeolis.threshold.compute`, :py:mod:`aeolis.vegetation.grow`, :py:mod:`aeolis.bed.update`.

process_wind
^^^^^^^^^^^^^
*Process_wind* makes sure the wind file is loaded, and interpolates values of the wind speed and direction to each time step. The model does not work without this flag. Used in :py:mod:`aeolis.wind.interpolate`

process_threshold
^^^^^^^^^^^^^^^^^
*Process_threshold* allows for the alterations of the threshold velocity by processes like grain, moisture and vegetation. This process does not occur if a threshold file is provided as input since this file is used to define the threshold shear velocity. Used in :py:mod:`aeolis.threshold.compute`, documentation of threshold alterations can be found in :ref:model_description

**th_grainsize**: calculates the threshold velocity based on the grain size following Bagnold (:py:mod:`aeolis.threshold.compute_grain_size`)
**th_bedslope**: currently not implemented, but theoretically would include an alteration of the velocity threshold based on the slope of the bed. (:py:mod:`aeolis.threshold.compute_bedslope`)
**th_moisture**: alters the threshold velocity based on the moisture content, many different methods are available (:py:mod:`aeolis.threshold.compute_moisture`). Only works if moisture content is defined, which is calculated when process_moisture is on. 
**th_salt**: alters the wind velocity threshold based on salt content following Nickling and Ecclestone (1981) (:py:mod:`aeolis.threshold.compute_salt`) 
**th_sheltering**: modify the wind velocity threshold based on the presence of roughness elements in the grain size fractions following Raupach (1993) (:py:mod:`aeolis.threshold.compute_sheltering`) 
**th_humidity** and **th_drylayer**: are currently not implemented

process_transport
^^^^^^^^^^^^^^^^^^
*Process_transport* allows the calculation of the equilibrium transport rate based on a user-defined transport method. 

**Method_transport** defines the sediment transport equation used in the calculation of the equilibrium transport rate. Options are: *bagnold, bagnold_gs, kawamura, lettau, dk, sauermann, vanrijn_strypsteen*.

**Method_grainspeed** defines at which speed the sediment transport in the air is occurring. Options are: *duran*/*duran_full*, *windspeed*, and *constant*

Used in :py:mod:`aeolis.transport.equilibrium`

process_bedupdate
^^^^^^^^^^^^^^^^^
Process_bedupdate allows the bed level to change based on calculated erosion/deposition. Used in :py:mod:`aeolis.bed.update`

process_shear
^^^^^^^^^^^^^
Process_shear calls the shear module in shear.py to calculate the shear stress perturbation caused by topography. Used in :py:mod:`aeolis.wind.initialize`

process_tide
^^^^^^^^^^^^
Process_tide changes the threshold velocity to infinity (no aeolian transport) if the bed level is below the water level. This process is not used if th_moist is used. Used in :py:mod:`aeolis.threshold.compute`, :py:mod:`aeolis.vegetation.grow`, :py:mod:`aeolis.bed.update`.

process_wave
^^^^^^^^^^^^
Process_wave allows calculation of the water depth based on the input tide file and interpolates the input wave data to the timesteps of the model run. If the wave file is not available, the wave height and peak period are set to 0. Turning this process flag on also results in the calculation of Hsmix, which is needed for the calculation of the Depth of Disturbance (*process_mixtoplayer*). The initialization/calculation is skipped if external variables are imported from another model.

process_runup
^^^^^^^^^^^^
Process_runup allows calculation of the runup extent based on the wave height, peak period and water level. Process_wave and Process_tide need to be on for this to work. The runup is calculated with the Stockdon equation using a user-defined, static beach slope. The initialization/calculation is skipped if external variables are imported from another model.

process_moist
^^^^^^^^^^^^^
Process_moist allows calculation of the soil moisture content, based on different methods, infiltration or surface_moist

method_moist_process

method_moist_threshold defines the equation used to calculate teh threshold shear veolcity based on the moisture content. Used in :py:mod:`aeolis.threshold.compute_moisture`.

Process_groundwater, Process_seepage_face and Process_scanning are all related to the calculation of the moisture content.

process_mixtoplayer
^^^^^^^^^^^^^^^^^^^
This process flag allows mixing in the layers that are present down to the depth of disturbance. For the calculation of the DoD the process_wave need to be on.

process_wet_bed_reset
^^^^^^^^^^^^^^^^^^^^^
Resets the bed to the original bathymetry if the bed is under water (zs). Used in :py:mod:`aeolis.bed.wet_bed_reset`. The execution of the wet bed reset is dependent on the TWL calculation, which can be turned on process_runup, process_waves and process_tide.

process_meteo
^^^^^^^^^^^^^
This is a place holder and currently has no functionality

process_avalanche
^^^^^^^^^^^^^^^^^^
Simulates the process of avalanching when slopes of the bed become too steep to be realistic (i.e. > a critical static slope).

process_separation
^^^^^^^^^^^^^^^^^^
This enables the calculation of the separation bubble within the shear perturbation module. Before executing the calculation is checks whether steep slopes are present that might lead to a separation bubble. Process_separation will only be used if process_shear is on.

process_vegetation
^^^^^^^^^^^^^^^^^^^
This process flag allows application of shear stress reduction due to vegetation based on Raupach or Okin. It also allows for germination and lateral growth of vegetation if those values are set to larger than 0. This process is actively being developed.

process_fences
^^^^^^^^^^^^^^
This process enables alteration of the shear velocity if fence characteristics are provided as user input. Calculations happen in 1D or 2D depending on grid size following the Okin model. 

process_dune_erosion
^^^^^^^^^^^^^^^^^^^^
This flag turns on dune erosion calculation (:py:mod:`aeolis.erosion`.) based on the Palmsten and Holman (2012) method. After calculating the erosion, the avalanching routine is run in :py:mod:`aeolis.model.update`. This is needed because these modules only get called for aeolian transport in case of winds above threshold.


Model state/output
------------------

The AeoLiS model state is described by a collection of spatial grid
variables with at least one value per horizontal grid cell.  Specific
model state variables can also be subdivided over bed composition
layers and/or grain size fractions.  All model state variables can be
part of the model netCDF4 output. The current model state variables
are listed below.

.. literalinclude:: ../../aeolis/constants.py
   :language: python
   :start-after: #: Aeolis model state variables
   :end-before: #: AeoLiS model default configuration


Solver selection
----------------
