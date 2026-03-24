Model setup
=================
Setting up an AeoLiS model involves configuring various parameters and input files. In this section, we will cover:

- **Model Input**: The input files and their formats.
- **Model Output**: How to analyse the output.
- **Default settings**: Overview of the model parameters, incl. default values.
- **Activate/deactivate processes**: Process flags.
- **Model state/output**: Overview of model state variables.
- **Guidance on solver use**

A general tip for setting up an AeoLiS model is to start simple and build up the complexity. As a starting point, use the default values and then slowly deviate from those and turn on processes. Start with a relatively coarse grid to allow for faster simulation testing. It is also highly recommended to start with constant conditions, such as a one-directional constant wind, making it easier to interpret if the model output is logical. 

While testing, go through the sequential steps AeoLiS takes, plot the variables, and critically assess whether the results are expected. Important model state variables to check include: wind speed and direction (``uws``, ``uwn``), shear velocity (``ustars``, ``ustarn``), velocity threshold (``uth``), sediment concentrations (``Cu``, ``Ct``), pickup (``pickup``), and bed level change (``zb``, ``dzb``).

In case you run into issues, we encourage users to post questions and case studies on the `AeoLiS Discussion Board`_. We use this public forum so our help and advice is available to everyone.

.. _AeoLiS Discussion Board: https://github.com/openearth/aeolis-python/discussions

Model input
-----------

The computational grid and boundary conditions for AeoLiS are specified through external
input files called by the main configuration file. The computational domain is defined
using specific spatial files (``*.grd``), while boundary conditions for wind, wave, and tides
are provided via ``*.txt`` files. An overview of these files is provided in the table below.

.. list-table:: 
   :widths: 15 15 15 15 40
   :header-rows: 1

   * - Input File
     - Keyword
     - Dimensions
     - Requirement
     - File Description
   * - aeolis.txt
     - N/A
     - N/A
     - Mandatory
     - Main file containing parameter definitions
   * - x.grd
     - ``xgrid_file``
     - (ny, nx)
     - Mandatory
     - File containing cross-shore grid coordinates
   * - y.grd
     - ``ygrid_file``
     - (ny, nx)
     - Mandatory
     - File containing alongshore grid coordinates
   * - z.grd
     - ``bed_file``
     - (ny, nx)
     - Mandatory
     - File containing topography and bathymetry data (bed level)
   * - zne.grd
     - ``ne_layer_file``
     - (ny, nx)
     - Optional
     - File containing the non-erodible layer elevation
   * - veg.grd
     - ``veg_file``
     - (ny, nx)
     - Optional
     - Initial vegetation density (if ``process_vegetation = T``)
   * - hveg.grd
     - ``hveg_file``
     - (ny*nx*nspecies)
     - Optional
     - Vegetation height per species (if ``method_vegetation = grass``)
   * - Nt.grd
     - ``Nt_file``
     - (ny*nx*nspecies)
     - Optional
     - Vegetation density/tillers per species (if ``method_vegetation = grass``)     
   * - mass.txt
     - ``mass_file``
     - (nx*ny, nfractions*nlayers)
     - Optional
     - Sediment mass data (for space-varying grain sizes)
   * - wind.txt
     - ``wind_file``
     - (ntimesteps, 3)
     - Mandatory
     - File containing wind speed and direction data
   * - tide.txt
     - ``tide_file``
     - (ntimesteps, 2)
     - Optional
     - Water elevation data (if ``process_tide = T``)
   * - wave.txt
     - ``wave_file``
     - (ntimesteps, 3)
     - Optional
     - Wave height and period data (if ``process_wave = T``)
   * - meteo.txt
     - ``meteo_file``
     - (ntimesteps, 6)
     - Optional
     - Meteorological data (if ``process_groundwater = T``)
   * - mask_*.grd
     - ``mask_tide``, etc.
     - (ny, nx)
     - Optional
     - Spatial masks for tide, wave, or runup boundary conditions

Main configuration file (e.g., aeolis.txt)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The main configuration file controls all processes and boundary conditions.
Parameters in the file are specified by various keywords; each keyword has a pre-defined
default value (see ``constants.py`` and the Default settings tab) that will be used if it is not directly specified. Among the keywords
are those defining external grid files (``xgrid_file``, ``ygrid_file``,
``bed_file``) and external boundary conditions (``tide_file``, ``wave_file``, ``wind_file``).  
Physical processes in AeoLiS can be toggled by setting process keywords to True (``T``) or False (``F``). 
Example parameter files can be found in the examples folder on the AeoLiS GitHub.

Properties of \*.grd files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
All grid files (``x.grd``, ``y.grd``, ``z.grd``, ``zne.grd``, ``veg.grd``, etc.) must have the exact same dimensions. 
Each value (or element) within these matrices represents a single computational cell. This means that element `[i, j]` 
in the ``x.grd`` file corresponds to the exact same physical cell as element `[i, j]` in the ``z.grd`` file. 

The model can be run in 1D and 2D mode depending on the input dimensions. Most processes are implemented in 2D mode. To run the model in 2D mode, all grid files should contain 2D matrices ``(ny, nx)`` of the same size. To run the model in 1D mode, all grid files should contain 1D vectors ``(nx, 1)``. Because some processes are easier to solve in 2D, the 1D model is internally converted to a quasi-2D model by repeating the vectors three times ``(nx, 3)``, assuming the same resolution spacing as the cross-shore direction. Results are then converted back to 1D by extracting the middle vector. 

x.grd and y.grd (Spatial Grid)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The ``x.grd`` and ``y.grd`` files define the computational grid in meters. It is important to understand that the boundary definitions and cross-shore/longshore directions are determined by the shape and structure of the grid matrices, while the physical location and orientation are determined by the specific coordinate values inside those matrices. To ensure correct model execution and boundary alignment, grid generation should follow this step-by-step approach:

1. **Define dimensions and resolution:** Set your domain length and grid spacing. Because the model does not yet support variable grid sizes, the cross-shore resolution must equal the alongshore resolution (e.g., ``dx = 10; dy = 10``).
2. **Create domain axes:** Generate 1D arrays for the cross-shore and alongshore directions. The cross-shore array dictates the domain boundaries and **must be ascending**; the first x-element is always the onshore boundary, and the last is the offshore boundary (e.g., ``x = np.arange(0, Lx + dx, dx)``).
3. **Generate the grid:** Use the 1D axes to create 2D coordinate matrices (e.g., ``X, Y = np.meshgrid(x, y)``).
4. **Apply physical orientation:** Shift and rotate the coordinate values within these matrices to match the desired real-world location and orientation (e.g., ``Xr = X * np.cos(theta) - Y * np.sin(theta) + x0``).

By handling the physical orientation directly within the grid coordinate values, the model orientation keyword ``alfa`` is no longer required and should be left at its default value of ``0``.

z.grd (bed level)
~~~~~~~~~~~~~~~~~
The ``z.grd`` file provides the surface elevation for every cell defined in the spatial grids. Elevation values should be defined such that positive values are above the vertical datum (up) and negative values are below (down). 

zne.grd (non-erodible Layer)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The ``zne.grd`` file defines the elevation of the hard, non-erodible layer beneath the sand surface. The model will not erode sediment below this specified elevation. It follows the same dimensions and elevation datum as ``z.grd``.

veg.grd (vegetation)
~~~~~~~~~~~~~~~~~~~~
The ``veg.grd`` file is an optional grid providing the initial vegetation coverage (density) at each cell. 

.. _fig-veg-inputs:

.. figure:: /images/vegetation_text_file.jpeg
   :alt: vegetation input format
   :width: 200px
   :align: center
   
   File format for a 1D AeoLis vegetation grid.  Each red dot is the vegetation density at a specific location in the computational grid.

Masks (tide, wave, and runup)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Masks (e.g., ``mask_tide``, ``mask_wave``, ``mask_runup``) can be used to spatially modify boundary conditions. This is particularly useful when landward elevations are lower than the offshore water level but remain dry, or when water bodies are disconnected from the offshore (such as a barrier island system where inland water experiences no tidal or wave action).

Masks utilize complex numbers to apply both a scaling multiplier and a static offset to the boundary condition, following the logic: ``Applied Value = (Real Value * Boundary Value) + Complex Value``. To completely zero out the boundary value in a specific area (e.g., no waves inland), use ``0``. To set a fixed static value regardless of the incoming boundary condition (e.g., an inland lake permanently fixed at +2m water level), use a complex value such as ``0 + 2j``.

Example Python script
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    import numpy as np

    # 1. Define dimensions and resolution
    Lx, Ly = 1000, 500  # Lx is always cross-shore, Ly is longshore
    dx = 10             # dx must equal dy
    
    x0, y0 = 5000, 5000 # coordinates for a grid corner point
    theta = np.radians(30) # grid rotation in radians

    # 2. Create axes
    x = np.arange(0, Lx + dx, dx) # must be ascending (first=onshore, last=offshore)
    y = np.arange(0, Ly + dx, dx) 
    
    # 3. Create orthogonal grid
    X, Y = np.meshgrid(x, y)

    # 4. Rotate and shift domain to real-world coordinates
    Xr = X * np.cos(theta) - Y * np.sin(theta) + x0
    Yr = X * np.sin(theta) + Y * np.cos(theta) + y0
    # *Because orientation is baked into the coordinates, set `alfa = 0` in aeolis.txt*

    # 5. Populate specific grid variables
    z = np.zeros_like(X) # example: Flat bed level at 0.0m
    zne = np.full_like(z, -2.0) # example: Flat non-erodible layer at -2.0m
    veg = np.zeros_like(z) # example: no initial vegetation

    # 6. Create masks (e.g., tide mask for an inland lake)
    tide_mask = np.ones_like(X, dtype=complex) 
    tide_mask[:, :10] = 0.0 + 2.0j # 0.0 multiplier (no tidal variation) + 2.0j elevated

    # Save outputs to text files (to be read as .grd)
    np.savetxt('x.grd', Xr)
    np.savetxt('y.grd', Yr)
    np.savetxt('z.grd', z)
    np.savetxt('zne.grd', zne)
    np.savetxt('veg.grd', veg)
    np.savetxt('mask_tide.grd', tide_mask)

Multi-dimensional Inputs (mass, hveg, Nt)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Some variables require multiple dimensions per spatial grid cell, such as spatially varying grain sizes (multiple fractions and bed layers) or complex vegetation cover (multiple species). AeoLiS requires these multi-dimensional arrays to be flattened into two-dimensional formats before saving as textfiles.

* **mass.txt**: Defines the mass of each sediment fraction per bed layer. If the grain size distribution is uniform across the domain, which is most often the case, this file is not needed (see keywords ``grain_dist`` and ``grain_size`` in the configuration file). The 4D shape of `(ny, nx, nlayers, nfractions)` must be reshaped into a 2D matrix of shape `(nx * ny, nfractions * nlayers)`. The rows represent the flattened spatial coordinates, and the columns are grouped by bed layer (e.g., Layer 1: Fraction 1, Fraction 2; Layer 2: Fraction 1, Fraction 2).
* **hveg.grd and Nt.grd**: When using the grass vegetation method (``method_vegetation = grass``), the model needs vegetation height (``hveg``) and density (``Nt``). Because this method supports multiple interacting species, the 3D shape of `(ny, nx, nspecies)` is flattened into a 1D vector of shape `(ny * nx * nspecies)`. 

.. _fig-mass-inputs-2D:

.. figure:: /images/mass_text_file_2D.jpeg
   :alt: mass file format 2D
   :width: 550px
   :align: center
   
   File format for a 2D AeoLis mass file for spatially variable grain size distributions. Each red dot is the mass for a specific sediment fraction within a specific bed layer at a given spatial location.

Example Python script
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Below is a simple Python script demonstrating how to load spatial grid dimensions and properly flatten multi-dimensional arrays for AeoLiS.

.. code-block:: python

    import numpy as np

    # 1. Load spatial grid to get dimensions
    X = np.loadtxt('x.grd')
    ny, nx = X.shape
    n_cells = ny * nx

    # 2. Create mass.txt (sediment fractions and layers)
    n_layers = 3
    n_fractions = 4
    
    # Example: a uniform sediment distribution across the domain
    mass_per_fraction = np.array([0.4, 0.3, 0.2, 0.1]) # kg per fraction
    
    # Tile the fraction distribution to fill all bed layers and spatial cells
    # Resulting shape: (n_cells, n_fractions * n_layers)
    mass_matrix = np.tile(mass_per_fraction, (n_cells, n_layers)) 
    np.savetxt('mass.txt', mass_matrix)

    # 3. Create hveg.grd and Nt.grd (Multi-species vegetation)
    n_species = 2
    
    # Initialize 3D arrays
    hveg = np.zeros((ny, nx, n_species))
    Nt = np.zeros((ny, nx, n_species))

    # Example: species 0 is 0.5m tall, species 1 is 1.0m tall
    hveg[:, :, 0] = 0.5 
    hveg[:, :, 1] = 1.0
    
    # Example: species 0 has 10 tillers/m^2, species 1 has 5 tillers/m^2
    Nt[:, :, 0] = 10.0 
    Nt[:, :, 1] = 5.0  

    # Flatten the 3D arrays into 1D vectors for AeoLiS
    hveg_flat = hveg.reshape(-1)
    Nt_flat = Nt.reshape(-1)

    np.savetxt('hveg.grd', hveg_flat)
    np.savetxt('Nt.grd', Nt_flat)


Time-series (*.txt)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Environmental forcing in AeoLiS—such as wind, water levels, and wave conditions—is provided through external time-series files. For all of these files, the format follows the same structure: the first column represents time in seconds w.r.t. ``refdate`` in the configruation file. The following columns contain the variables at those given times. The model will automatically interpolate these data points to match the modelling time steps.

wind.txt
~~~~~~~~
The ``wind.txt`` file provides the wind boundary conditions driving aeolian transport. It contains three columns: 1. Time (s), 2. Wind speed (m/s), 3. Wind direction (degrees). Wind directions can be specified in either nautical or cartesian convention, which is set in the ``aeolis.txt`` file using the ``wind_convention`` keyword. 

.. _fig-wind-inputs:

.. figure:: /images/wind_text_file_graphic.jpeg
   :alt: wind input format
   :width: 300px
   :align: center
   
   File format for wind boundary conditions file for AeoLis input.

tide.txt
~~~~~~~~
The ``tide.txt`` file contains the water elevation data for the duration of the simulation. It contains two columns: 1. Time (s), 2. Water elevation (m)

.. _fig-tide-inputs:

.. figure:: /images/tide_text_file.jpeg
   :alt: tide input format
   :width: 300px
   :align: center
   
   File format for the water elevation conditions file for AeoLis input.
   
wave.txt
~~~~~~~~
The ``wave.txt`` file provides the wave data used by AeoLiS to calculate runup. It contains three columns: 1. Time (s), 2. Significant wave height (m), 3. Peak wave period (s)

.. _fig-wave-inputs:

.. figure:: /images/wave_text_file_graphic.jpeg
   :alt: wave input format
   :width: 300px
   :align: center
   
   File format for the wave conditions file for AeoLis input.

meteo.txt
~~~~~~~~~
The ``meteo.txt`` file contains meteorological data and is only required if using the groundwater module by Hallin (2023) to simulate surface moisture. It contains six columns: 1. Time (s), 2. Temperature (°C), 3. Precipitation (mm/hr), 4. Relative humidity (%), 5. Global radiation (MJ/$m^2$/day), 6. Air pressure (kPa)

.. _fig-meteo-inputs:

.. figure:: /images/meteo_file_format.jpeg
   :alt: meteo file format
   :width: 550px
   :align: center
   
   File format for meteorological data used to simulate surface moisture in AeoLiS.

Example Python script
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Below is a simple Python script demonstrating how to generate and save these time-series files using NumPy. In this example, we generate a 10-day simulation with a constant onshore wind, a harmonic tide, and constant wave conditions.

.. code-block:: python

    import numpy as np

    # 1. Define time array
    days = 10
    dt = 3600 # 1-hour intervals
    time_sec = np.arange(0, days * 24 * 3600 + dt, dt) 

    # 2. Wind (constant onshore wind)
    # Assuming 270 degrees is straight onshore
    wind_speed = np.full_like(time_sec, 8.0)  # 8 m/s
    wind_dir = np.full_like(time_sec, 270.0)  # 270 degrees
    wind_data = np.column_stack((time_sec, wind_speed, wind_dir))

    # 3. Tide (harmonic tide)
    tide_amp = 1.0 # 1 meter amplitude
    tide_period = 12 * 3600 # 12-hour period in seconds
    water_level = tide_amp * np.sin(2 * np.pi * time_sec / tide_period)
    tide_data = np.column_stack((time_sec, water_level))

    # 4. Waves (constant wave height and period)
    wave_height = np.full_like(time_sec, 1.5) # 1.5m significant wave height
    wave_period = np.full_like(time_sec, 6.0) # 6s peak period
    wave_data = np.column_stack((time_sec, wave_height, wave_period))

    # 5. Save to text files
    np.savetxt('wind.txt', wind_data, fmt='%.2f')
    np.savetxt('tide.txt', tide_data, fmt='%.2f')
    np.savetxt('wave.txt', wave_data, fmt='%.2f')

Visualizating input settings
^^^^^^^^^^^^^^^^^^^^
When running the model, you can automatically generate diagnostic plots by setting the keyword ``visualization = True`` in your ``aeolis.txt`` file. This will automatically generate figures in your simulation folder to help verify your setup:

* ``figure_grid_initialization.png``: Displays the grid orientation and boundaries.
* ``figure_params_initialization.png``: Shows the most relevant spatial parameters mapped onto the domain.
* ``figure_timeseries_initialization.png``: Visualizes the time series of your environmental boundary conditions.


.. _default-settings:

Model output
------------

Default settings
-----------------

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
-------------------------------
After creating the input files that are necessary to run an AeoLiS model, the next step is often to decide which processes and methods to use. Several processes are defined in the configuration file that can be turned on and off. Apart from turning processes on and off, there are also several user-defined thresholds and methods that affect the way in which processes are calculated. For example, there are different sediment transport equations available within *process_transport*. The default is Bagnold, but by defining *method_transport* in the configuration file a different equation can be used. Here, we provide a description of the processes and methods that are defined in configuration file. More detailed descriptions of the processes and their implementation can be found in :ref:`the model description <model_description>`.

An easy way to look up where these process, threshold and method flags are used is by going to the main page of the AeoLiS github and using the search bar at the top. For instance, searching *process_tide* shows that it is used in :py:mod:`aeolis.threshold.compute`, :py:mod:`aeolis.vegetation.grow`, :py:mod:`aeolis.bed.update`.

process_wind
^^^^^^^^^^^^^
*Process_wind* makes sure the wind file is loaded, and interpolates values of the wind speed and direction to each time step. The model does not work without this flag. Used in :py:mod:`aeolis.wind.interpolate`

process_threshold
^^^^^^^^^^^^^^^^^
*Process_threshold* allows for the alterations of the threshold velocity by processes like grain and moisture. This process does not occur if a threshold file is provided as input since this file is used to define the threshold shear velocity. Used in :py:mod:`aeolis.threshold.compute`, documentation of threshold alterations can be found in :ref:model_description

- **th_grainsize**: calculates the threshold velocity based on the grain size following Bagnold (:py:mod:`aeolis.threshold.compute_grain_size`)
- **th_bedslope**: currently not implemented, but theoretically would include an alteration of the velocity threshold based on the slope of the bed. (:py:mod:`aeolis.threshold.compute_bedslope`)
- **th_moisture**: alters the threshold velocity based on the moisture content, many different methods are available (:py:mod:`aeolis.threshold.compute_moisture`). Only works if moisture content is defined, which is calculated when process_moisture is on. 
- **th_salt**: alters the wind velocity threshold based on salt content following Nickling and Ecclestone (1981) (:py:mod:`aeolis.threshold.compute_salt`) 
- **th_sheltering**: modify the wind velocity threshold based on the presence of roughness elements in the grain size fractions following Raupach (1993) (:py:mod:`aeolis.threshold.compute_sheltering`) 
- **th_humidity** and **th_drylayer**: are currently not implemented

process_transport
^^^^^^^^^^^^^^^^^^
*Process_transport* allows the calculation of the equilibrium transport rate based on a user-defined transport method. 

**Method_transport** defines the sediment transport equation used in the calculation of the equilibrium transport rate. Options are: *bagnold, bagnold_gs, kawamura, lettau, dk, sauermann, vanrijn_strypsteen*.

**Method_grainspeed** defines at which speed the sediment transport in the air is occurring. Options are: *duran*/*duran_full*, *windspeed*, and *constant*

Used in :py:mod:`aeolis.transport.equilibrium`

process_bedupdate
^^^^^^^^^^^^^^^^^^
Process_bedupdate allows the bed level to change based on calculated erosion/deposition. Used in :py:mod:`aeolis.bed.update`

process_shear
^^^^^^^^^^^^^^
Process_shear calls the shear module in shear.py to calculate the shear stress perturbation caused by topography. Used in :py:mod:`aeolis.wind.initialize`

process_tide
^^^^^^^^^^^^^
Process_tide changes the threshold velocity to infinity (no aeolian transport) if the bed level is below the water level. This process is not used if th_moist is used. Used in :py:mod:`aeolis.threshold.compute`, :py:mod:`aeolis.vegetation.grow`, :py:mod:`aeolis.bed.update`.

process_wave
^^^^^^^^^^^^^
Process_wave allows calculation of the water depth based on the input tide file and interpolates the input wave data to the timesteps of the model run. If the wave file is not available, the wave height and peak period are set to 0. Turning this process flag on also results in the calculation of Hsmix, which is needed for the calculation of the Depth of Disturbance (*process_mixtoplayer*). The initialization/calculation is skipped if external variables are imported from another model.

process_runup
^^^^^^^^^^^^^^
Process_runup allows calculation of the runup extent based on the wave height, peak period and water level. Process_wave and Process_tide need to be on for this to work. The runup is calculated with the Stockdon equation using a user-defined, static beach slope. The initialization/calculation is skipped if external variables are imported from another model.

process_moist
^^^^^^^^^^^^^^
Process_moist allows calculation of the soil moisture content, based on different methods, infiltration or surface_moist

method_moist_process

method_moist_threshold defines the equation used to calculate teh threshold shear veolcity based on the moisture content. Used in :py:mod:`aeolis.threshold.compute_moisture`.

Process_groundwater, Process_seepage_face and Process_scanning are all related to the calculation of the moisture content.

process_mixtoplayer
^^^^^^^^^^^^^^^^^^^^
This process flag allows mixing in the layers that are present down to the depth of disturbance. For the calculation of the DoD the process_wave need to be on.

process_wet_bed_reset
^^^^^^^^^^^^^^^^^^^^^
Resets the bed to the original bathymetry if the bed is under water (zs). Used in :py:mod:`aeolis.bed.wet_bed_reset`. The execution of the wet bed reset is dependent on the TWL calculation, which can be turned on process_runup, process_waves and process_tide.

process_meteo
^^^^^^^^^^^^^^
This is a place holder and currently has no functionality

process_avalanche
^^^^^^^^^^^^^^^^^^
Simulates the process of avalanching when slopes of the bed become too steep to be realistic (i.e. > a critical static slope).

process_separation
^^^^^^^^^^^^^^^^^^^
This enables the calculation of the separation bubble within the shear perturbation module. Before executing the calculation is checks whether steep slopes are present that might lead to a separation bubble. Process_separation will only be used if process_shear is on.

process_vegetation
^^^^^^^^^^^^^^^^^^^^
This process flag allows application of shear stress reduction due to vegetation based on Raupach or Okin. It also allows for germination and lateral growth of vegetation if those values are set to larger than 0. This process is actively being developed.

process_fences
^^^^^^^^^^^^^^^
This process enables alteration of the shear velocity if fence characteristics are provided as user input. Calculations happen in 1D or 2D depending on grid size following the Okin model. 

process_dune_erosion
^^^^^^^^^^^^^^^^^^^^^
This flag turns on dune erosion calculation (:py:mod:`aeolis.erosion`.) based on the Palmsten and Holman (2012) method. After calculating the erosion, the avalanching routine is run in :py:mod:`aeolis.model.update`. This is needed because these modules only get called for aeolian transport in case of winds above threshold.


Model state/output
-------------------

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


Guidance on solver use
------------------------

Different numerical solvers are available in the latest AeoLiS version. 
The numerical solvers are used to solve the transport equation numerically.
Other modules such as the shear module, vegetation module, and the moisture module
use other equations and numerical implementations that currently do not have
different options for numerical solvers.

The advection equation is implemented in two-dimensional form
following:

.. math::
   :label: apx-advection
   
   \frac{\partial c}{\partial t} +
   u_{z,\mathrm{x}} \frac{\partial c}{\partial x} + 
   u_{z,\mathrm{y}} \frac{\partial c}{\partial y} = 
   \frac{c_{\mathrm{sat}} - c}{T}

in which :math:`c` [:math:`\mathrm{kg/m^2}`] is the sediment mass per
unit area in the air, :math:`c_{\mathrm{sat}}` [:math:`\mathrm{kg/m^2}`] is the
maximum sediment mass in the air that is reached in case of
saturation, :math:`u_{z,\mathrm{x}}` and :math:`u_{z,\mathrm{y}}` are the x- and
y-component of the wind velocity at height :math:`z` [m], :math:`T` [s] is an
adaptation time scale, :math:`t` [s] denotes time and :math:`x` [m] and :math:`y` [m]
denote cross-shore and alongshore distances respectively.

The formulation is discretized in different ways to allow for different types of simulations balancing accuracy vs. computational resources. 
The conservative method combined with a steady state solution is the current default for most simulations.
Non-conservative methods and explicit/implicit Euler forward/backward schemes are also available.

The available solvers are *steadystate*, *trunk*, and *pieter*. Some details are given below. As of version 3 of the AeoLiS
model, the steadystate solver is the default solver. The steadystate solver is most suitable for practical cases. 
The other solvers are still available for specific applications. 

steadystate (default since v3)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The *steadystate* solver is based on the assumption that dc/dt = 0 and uses a finite difference scheme to solve the transport equation.

.. math::
   :label: ss-advection
   
   u_{z,\mathrm{x}} \frac{\partial c}{\partial x} + 
   u_{z,\mathrm{y}} \frac{\partial c}{\partial y} = 
   \frac{c_{\mathrm{sat}} - c}{T}

When solving for equation :eq:`ss-advection` a sweeping algorithm is used that propagates the boundary conditions through the 
4 possible quadrants of the computational grid. The 4 quadrants depend on the signs of the sediment velocities and the 
remaining grid cells that are not part of a quadrant (because winds diverge or converge in that cell) are solved as well.

The steadystate solver is most suitable for case study simulations with larger timeframes and timesteps. All landform 
simulations in the :cite:t:`VANWESTEN2024106093` publication were done with the steadystate solver.

trunk
^^^^^
The *trunk* solver was the first solver that was implemented in AeoLiS. The trunk solver allows a time-varying solution for
sediment concentration with options for explicit and implicit Euler forward/backward schemes. The 1D simulations by :cite:t:`deVries2014a`
were done with the trunk solver in explicit mode. However, the explicit mode is not recommended for most simulations as 
very strict requirements for stability are needed which results in large calculation times. The implicit mode is more stable 
and allows for larger timesteps. However, the implicit numerical scheme lacks accuracy when larger timesteps are used. The 
2D simulations by :cite:t:`Hoonhout2016` were done with the trunk solver in implicit mode.

See :ref:`trunk_num` for details on the numerical implementation of the trunk solver.

Pieter
^^^^^^
The *pieter* solver was built on the basis of using a conservative numerical scheme. This conservative scheme allowed for a better 
implementation of spatially varying wind(/sediment) velocities. In simple cases (spatially non-varying winds) the solver is 
identical to the trunk solver. The solver was built by Professor Pieter Rauwoens, hence the name. 

See :ref:`pieter_num` for details on the numerical implementation of the Pieter solver.
