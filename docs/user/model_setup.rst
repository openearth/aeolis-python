Model in- & output
=================
Setting up an AeoLiS model involves configuring various parameters and input files. In this section, we will cover:

- **Model Input**: The input files and their formats.
- **Model Output**: How to analyse the output.
- **Default settings**: Overview of the model parameters, incl. default values.
- **Activate/deactivate processes**: Process flags.
- **Model state/output**: Overview of model state variables.
- **Guidance on solver use**

A general tip for setting up an AeoLiS model is to start simple and build up the complexity. As a starting point, use the default values and then slowly deviate from those and turn on processes. Start with a relatively coarse grid to allow for faster simulation testing. It is also highly recommended to start with constant conditions, such as a one-directional constant wind, making it easier to interpret if the model output is logical. 

While testing, go through the sequential steps AeoLiS takes, plot the variables, and critically assess whether the results are expected. Important model state variables to check include: wind speed and direction (``uws``, ``uwn``), shear velocity (``ustars``, ``ustarn``), velocity threshold (``uth``), sediment concentrations (``Cu``, ``Ct``), pickup (``pickup``), and, finally, the resulting bed level change (``zb``, ``dzb``).

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

Configuration file (aeolis.txt)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The main configuration file controls all processes and input files.
Parameters in the file are specified by various keywords; each keyword has a pre-defined
default value (see ``constants.py`` and the Default settings tab) that will be used if it is not directly specified. Among the keywords
are those defining grid files (``xgrid_file``, ``ygrid_file``,
``bed_file``) and boundary conditions (``tide_file``, ``wave_file``, ``wind_file``).  
Physical processes in AeoLiS can be toggled by setting process keywords to True (``T``) or False (``F``). 
Example parameter files can be found in the examples folder on the AeoLiS GitHub.

Grid files (\*.grd)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
All grid files (``x.grd``, ``y.grd``, ``z.grd``, ``zne.grd``, ``veg.grd``, etc.) must have the exact same dimensions. 
Each value (or element) within these matrices represents a single computational cell. This means that element `[i, j]` 
in the ``x.grd`` file corresponds to the exact same physical cell as element `[i, j]` in the ``z.grd`` file. 

The model can be run in 1D and 2D mode depending on the input dimensions. Most processes are implemented in 2D mode. To run the model in 2D mode, all grid files should contain 2D matrices ``(ny, nx)`` of the same size. To run the model in 1D mode, all grid files should contain 1D vectors ``(nx, 1)``. Because some processes are easier to solve in 2D, the 1D model is internally converted to a quasi-2D model by repeating the vectors three times ``(nx, 3)``, assuming the same resolution spacing as the cross-shore direction. Results are then converted back to 1D by extracting the middle vector. 

x.grd and y.grd (coordinates)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The ``x.grd`` and ``y.grd`` files define the coordinates of the computational grid in meters. Important: The structure and order of the elements in these files determines the boundary definitions (i.e., which boundary is on which edge, what is longshore, what is cross-shore), while the specific coordinates determine the location and orientation of the domain. To ensure the model correctly interprets model grids, generation should follow this step-by-step approach (see also Python example further down below):

1. **Define dimensions and resolution:** Set your domain length and grid spacing. At this stage, the x-direction always represents cross-shore and y-direction longshore. The model does not (yet) support variable grid sizes; the cross-shore resolution must equal the alongshore resolution (e.g., ``dx = 10; dy = 10``).
2. **Create domain axes:** Generate 1D arrays for the cross-shore (x) and alongshore (y) directions. Both arrays must be ascending. In case of the x-array, the first element (``x[0]``) is at the onshore boundary, and the last element (``x[-1]``) is the offshore boundary (e.g., ``x = np.arange(0, Lx + dx, dx)``).
3. **Generate the grid:** Use the 1D axes to create 2D coordinate matrices (e.g., ``X, Y = np.meshgrid(x, y)``).
4. **Shift and rotate:** Now, shift and rotate the coordinatet to match the desired real-world location and orientation (e.g., ``Xr = X * np.cos(theta) - Y * np.sin(theta) + x0``). Following this order ensures both the structure of the files and the actual coordinates are correct.

By handling the physical orientation through the grid coordinates, the model orientation keyword ``alfa`` is no longer required and should be left untouched.

z.grd (bed level)
~~~~~~~~~~~~~~~~~
The ``z.grd`` file provides the surface elevation for every cell defined in the spatial grids. Elevation values should be defined such that positive values are above the vertical datum (up) and negative values are below (down). 

zne.grd (non-erodible Layer)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The ``zne.grd`` file defines the elevation of the non-erodible layer beneath the sand surface. The model will not erode sediment below this elevation. It follows the same dimensions and elevation datum as ``z.grd``.

veg.grd (vegetation)
~~~~~~~~~~~~~~~~~~~~
The ``veg.grd`` file is an optional grid providing the initial vegetation coverage (density) at each cell (0-1). 

.. _fig-veg-inputs:

.. figure:: /images/vegetation_text_file.jpeg
   :alt: vegetation input format
   :width: 200px
   :align: center
   
   File format for a 1D AeoLis vegetation grid.  Each red dot is the vegetation density at a specific location in the computational grid.

Masks (tide, wave, and runup)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Masks (e.g., ``mask_tide``, ``mask_wave``, ``mask_runup``) can be used to spatially modify boundary conditions. This can be useful when, for instance, onshore elevations are lower than the offshore water level but remain dry, or when water bodies are disconnected from the offshore (such as a barrier island where the lagoon has no tidal or wave action).

Masks use complex numbers to apply both a scaling multiplier (real) and a static offset (complex) to the boundary condition. For instance, for waveheight (Hs): ``Hs =  $\mathbb{R}$ * Hs + $\mathbb{C}$``. To half the boundary value in a specific area (e.g., sheltered waves inland), use ``0.5``. To set a fixed static value regardless of the incoming boundary condition (e.g., an inland lake permanently fixed at +2m water level), use a complex value such as ``0 + 2j``.

*Example Python script* 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Generating the grid files described above

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

Multi-dimensional input
~~~~~~~~~~~~~~~~~~~~~~~
Some variables require multiple dimensions per spatial grid cell, such as spatially varying grain sizes (multiple fractions and bed layers) or complex vegetation cover (multiple species). AeoLiS requires these multi-dimensional arrays to be flattened into lower-dimensional formats before they can be saved as textfiles.

* **mass.txt**: Defines the mass of each sediment fraction per bed layer. If the grain size distribution is uniform across the domain, which is most often the case, this file is not needed (see keywords ``grain_dist`` and ``grain_size`` in the configuration file). The 4D shape of `(ny, nx, nlayers, nfractions)` must be reshaped into a 2D matrix of shape `(nx * ny, nfractions * nlayers)`. The rows represent the flattened spatial coordinates, and the columns are grouped by bed layer (e.g., Layer 1: Fraction 1, Fraction 2; Layer 2: Fraction 1, Fraction 2).
* **hveg.grd and Nt.grd**: When using the grass vegetation method (``method_vegetation = grass``), the model needs vegetation height (``hveg``) and density (``Nt``). Because this method supports multiple interacting species, the 3D shape of `(ny, nx, nspecies)` is flattened into a 1D vector of shape `(ny * nx * nspecies)`. 

.. _fig-mass-inputs-2D:

.. figure:: /images/mass_text_file_2D.jpeg
   :alt: mass file format 2D
   :width: 550px
   :align: center
   
   File format for a 2D AeoLis mass file for spatially variable grain size distributions. Each red dot is the mass for a specific sediment fraction within a specific bed layer at a given spatial location.

Example Python script
~~~~~~~~~~~~~~~~~~~~~
Loading spatial grid dimensions and generating multi-dimensional input grids.

.. code-block:: python

    import numpy as np

    # 1. Load spatial grid to get dimensions
    X = np.loadtxt('x.grd')
    ny, nx = X.shape
    n_cells = ny * nx
    n_species = 2

    # 2a. Create mass.txt (sediment fractions and layers)
    n_layers = 3
    n_fractions = 4
    mass_per_fraction = np.array([0.4, 0.3, 0.2, 0.1]) # kg per fraction (uniform sediment distribution across the domain)
    
    # Create the mass matrix using tile
    mass_matrix = np.tile(mass_per_fraction, (n_cells, n_layers))  # (n_cells, n_fractions * n_layers)

    # Initialize the vegetation arrays
    hveg = np.zeros((ny, nx, n_species))
    Nt = np.zeros((ny, nx, n_species))

    # Fill the vegetation arrays
    hveg[:, :, 0] = 0.5 # first species is 0.5m tall
    hveg[:, :, 1] = 1.0 # second species is 1.0m tall
    Nt[:, :, 0] = 10.0 # tillers/m2
    Nt[:, :, 1] = 5.0  

    # Flatten the vegetation arrays
    hveg_flat = hveg.reshape(-1)
    Nt_flat = Nt.reshape(-1)

    # Save files
    np.savetxt('mass.txt', mass_matrix)
    np.savetxt('hveg.grd', hveg_flat)
    np.savetxt('Nt.grd', Nt_flat)


Time-series (*.txt)
^^^^^^^^^^^^^^^^^^^
Environmental forcing in AeoLiS (wind, water, waves) is provided through time-series files. For all of these files, the format follows the same structure: the first column represents time in seconds w.r.t. ``refdate`` in the configruation file. The other columns contain the variables at those given times. The model will automatically interpolate these data points to match the modelling time steps.

wind.txt
~~~~~~~~
The ``wind.txt`` file provides the wind boundary conditions driving aeolian transport. It contains three columns: 1. Time (s), 2. Wind speed (m/s), 3. Wind direction (degrees). Wind directions can be specified in either nautical or cartesian convention, which is set using the ``wind_convention`` keyword in the configuration file. 

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
The ``meteo.txt`` file contains meteorological data and is only required if using the groundwater module by Hallin (2023) to simulate surface moisture (``process_groundwater = T``). It contains six columns: 1. Time (s), 2. Temperature (°C), 3. Precipitation (mm/hr), 4. Relative humidity (%), 5. Global radiation (MJ/$m^2$/day), 6. Air pressure (kPa)

.. _fig-meteo-inputs:

.. figure:: /images/meteo_file_format.jpeg
   :alt: meteo file format
   :width: 550px
   :align: center
   
   File format for meteorological data used to simulate surface moisture in AeoLiS.

Example Python script
~~~~~~~~~~~~~~~~~~~~
Generating and saving time-series files. In this example, we generate a 10-day simulation with a constant onshore wind, a harmonic tide, and constant wave conditions.

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
When running the model, you can automatically generate diagnostic plots by setting the keyword ``visualization = T`` in the configuration file. This will automatically generate figures in your simulation folder to help verify your setup:

* ``figure_grid_initialization.png``: Displays the grid orientation and boundary definitions.
* ``figure_params_initialization.png``: Shows the most relevant spatial parameters mapped onto the domain.
* ``figure_timeseries_initialization.png``: Visualizes the time series of your environmental boundary conditions.

.. _default-settings:


All parameters and defaults
---------------------------

All available configuration parameters, and their default values, are collected in ``constants.py``. 
For any configuration parameters not defined in the model configuration file,
the default value is used.

.. literalinclude::   ../../aeolis/constants.py
   :language: python
   :start-after: #: AeoLiS model default configuration
   :end-before: #: Merge initial and model state


Processes, Thresholds, Methods, and Boundary Conditions
------------------
Processes in AeoLiS can be activated (``T``) or deactivated (``F``) via the configuration file. Several thresholds (``th_...``) can be set seperately and methods keywords (``method_...``) determine which methods are used to compute these processes. Seperate keywords set which boundary conditions are applied.

More detailed descriptions of these processes can be found in :ref:`the model description <model_description>`. To see exactly where specific process, threshold, and method flags are utilized in the code, an easy approach is to use the search bar on the main page of the AeoLiS GitHub repository.

List overview
^^^^^^^^^^^^^

* ``process_wind``: Interpolates wind speed and direction to each time step. *The model cannot run without this flag.*
* ``process_threshold``: Computes the threshold shear velocity based on various processes. (Note: This is overridden if a static threshold file is provided).
* ``process_transport``: Computes the sediment transport rates. 
* ``process_bedupdate``: Allows the bed level to change based on calculated erosion and deposition.
* ``process_shear``: Computes the shear stress perturbation caused by topography (needed for landform simulations).
* ``process_separation``: Computes flow separation bubbles over steep slopes (requires ``process_shear``).
* ``process_avalanche``: Computes avalanching when bed slopes exceed a critical static angle.
* ``process_tide``: Computes water level elevations, determining whether cells becomes submerged.
* ``process_wave``: Computes wave heights across the domain, including $Hs_{mix}$, required for the mixing of sediment.
* ``process_runup``: Computes the runup extent using the Stockdon equation (requires ``process_wave`` and ``process_tide``).
* ``process_moist``: Computes soil moisture content via different surface moisture methods. (Related flags include ``process_groundwater``, ``process_seepage_face``, and ``process_scanning``).
* ``process_mixtoplayer``: Mixes sediment fractions over several layers down to the Depth of Disturbance (requires ``process_wave``).
* ``process_wet_bed_reset``: Resets the bed to the original bathymetry if submerged. Execution depends on Total Water Level (TWL), requiring ``process_runup``, ``process_wave``, and ``process_tide``.
* ``process_vegetation``: Computes vegetation growth and shear stress reduction.
* ``process_fences``: Computes shear velocity reduction based on user-provided fence characteristics (Okin model).
* ``process_dune_erosion``: Computes dune erosion and triggers avalanching when water levels impact the dunes.
* ``process_meteo``: *(Placeholder flag; currently has no functionality).*

Thresholds
~~~~~~~~~~
* ``th_grainsize``: Computes threshold velocity based on grain size (base value, ``uth0``).
* ``th_moisture``: Modifies threshold velocity based on moisture content (requires ``process_moist``).
* ``th_sheltering``: Modifies threshold based on sheltering by roughness elements (coarser sediment fractions).
* ``th_nelayer``: Modifies threshold based on the presence of a non-erodible layer.
* *(Note: ``th_bedslope``, ``th_drylayer``, ``th_salt``, and ``th_humidity`` are currently not fully implemented or rarely used).*

Methods
~~~~~~~
*Note: More guidance on the decision for ``method_grainspeed``, ``method_shear``, and ``solver`` is given in the "Guidance on advection, shear and grainspeed solvers" section later on.*

* ``method_transport``: Defines the transport equation (Options: ``bagnold``, ``bagnold_gs``, ``kawamura``, ``lettau``, ``dk``, ``sauermann``, ``vanrijn_strypsteen``).
* ``method_grainspeed``: Defines the speed of sediment in the air (Options: ``windspeed``, ``duran``, ``duran_uniform``, ``duran_full``).
* ``method_shear``: Computes topographic effects on wind shear stress (Options: ``fft``, ``1Dstacks``).
* ``method_roughness``: Computes the roughness height (Options: ``constant``, ``constant_nikuradse``, ``vanrijn_strypsteen``).
* ``method_vegetation``: Defines the vegetation formulation (Options: ``duran``, ``grass``). The newest vegetation implementation is called through ``grass``.
* ``method_moist_process``: Computes soil moisture content (Options: ``infiltration``, ``surface_moisture``).
* ``method_moist_threshold``: Computes wind velocity threshold based on soil moisture (Options: ``belly_johnson``).
* ``solver``: Defines the numerical advection solver (Options: ``steadystate``, ``euler_backward``, ``euler_forward``).

Boundary Conditions
~~~~~~~~~~~~~~~~~~~
* ``boundary_offshore``, ``boundary_onshore``, ``boundary_lateral``: Defines the method for handling sediment concentrations at the domain edges. Options include:
  
  * ``constant``: Applies a zero-gradient boundary condition (incoming concentration equals the adjacent internal cell).
  * ``flux``: Applies a user-defined incoming sediment flux based on the equilibrium concentration.
  * ``circular``: Applies periodic boundaries where sediment leaving one side re-enters the opposite side (Note: offshore and onshore boundaries must both be set to ``circular`` together).

* ``offshore_flux``, ``onshore_flux``, ``lateral_flux``: Defines the incoming sediment concentration as a fraction of the equilibrium concentration (``Cu``) when the respective boundary is set to ``flux`` (e.g., ``1.0`` for fully saturated incoming wind, ``0.0`` for clean air with no sediment).

Typical Configurations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
When building an AeoLiS model, it helps to know which processes and methods are typically required for different environments. The table below provides a general baseline for five common simulation types.

.. list-table:: 
   :widths: 25 15 15 15 15 15
   :header-rows: 1

   * - Parameter / Flag
     - Wind-only [1]_
     - Barchan [2]_
     - Parabolic [3]_
     - Beach [4]_
     - Blowout [5]_
   * - **Processes**
     - 
     - 
     - 
     - 
     - 
   * - ``process_wind``
     - ✅
     - ✅
     - ✅
     - ✅
     - ✅
   * - ``process_threshold``
     - ✅
     - ✅
     - ✅
     - ✅
     - ✅
   * - ``process_transport``
     - ✅
     - ✅
     - ✅
     - ✅
     - ✅
   * - ``process_bedupdate``
     - ✅
     - ✅
     - ✅
     - ✅
     - ✅
   * - ``process_shear``
     - ❌
     - ✅
     - ✅
     - ✅
     - ✅
   * - ``process_separation``
     - ❌
     - ✅
     - ✅
     - ❌
     - ❌
   * - ``process_avalanche``
     - ❌
     - ✅
     - ✅
     - ✅
     - ✅
   * - ``process_vegetation``
     - ❌
     - ❌
     - ✅
     - ✅
     - ✅
   * - ``process_tide``
     - ❌
     - ❌
     - ❌
     - ✅
     - ✅
   * - ``process_wave``
     - ❌
     - ❌
     - ❌
     - ✅
     - ✅
   * - ``process_runup``
     - ❌
     - ❌
     - ❌
     - ✅
     - ✅
   * - ``process_moist``
     - ❌
     - ❌
     - ❌
     - ✅
     - ✅
   * - ``process_mixtoplayer``
     - ❌
     - ❌
     - ❌
     - ✅
     - ✅
   * - ``process_wet_bed_reset``
     - ❌
     - ❌
     - ❌
     - ✅
     - ✅
   * - ``process_dune_erosion``
     - ❌
     - ❌
     - ❌
     - ✅ (maybe)
     - ❌
   * - **Thresholds**
     - 
     - 
     - 
     - 
     - 
   * - ``th_grainsize``
     - ✅
     - ✅
     - ✅
     - ✅
     - ✅
   * - ``th_moisture``
     - ❌
     - ❌
     - ❌
     - ✅
     - ✅
   * - ``th_sheltering``
     - ❌
     - ❌
     - ❌
     - ✅
     - ✅
   * - ``th_nelayer``
     - ❌
     - ✅
     - ✅
     - ❌
     - ❌
   * - **Methods**
     - 
     - 
     - 
     - 
     - 
   * - ``solver``
     - ``steadystate``
     - ``steadystate``
     - ``steadystate``
     - ``steadystate``
     - ``steadystate``
   * - ``method_transport``
     - ``bagnold``
     - ``bagnold``
     - ``bagnold``
     - ``bagnold``
     - ``bagnold``
   * - ``method_grainspeed``
     - ``duran_uniform``
     - ``duran``
     - ``duran``
     - ``duran``
     - ``duran_full``
   * - ``method_shear``
     - -
     - ``fft``
     - ``fft``
     - ``1Dstacks``
     - ``fft``
   * - ``method_vegetation``
     - -
     - -
     - ``duran``
     - ``grass``
     - ``grass``
   * - **Boundary Conditions**
     - 
     - 
     - 
     - 
     - 
   * - ``boundary_lateral``
     - ``circular``
     - ``circular``
     - ``circular``
     - ``circular``
     - ``constant``
   * - ``boundary_offshore``
     - ``circular``
     - ``constant``
     - ``constant``
     - ``flux (0)``
     - ``flux (0)``
   * - ``boundary_onshore``
     - ``circular``
     - ``constant``
     - ``constant``
     - ``flux (0)``
     - ``flux (0)``
.. [1] **Wind-only:** A basic flat-bed simulation without topographic, vegetation or marine influences.
.. [2] **Barchan:** A migrating barchan simulation requiring topographic shear steering, flow separation, and avalanching.
.. [3] **Parabolic:** A vegetated dune simulation combining topographic steering with vegetation processes.
.. [4] **Beach:** A simplistic (semi-1D) beach-dune profile simulation describing foredune growth.
.. [5] **Blowout:** A complex coastal dune simulation involving a combination of most AeoLiS processes.

Model output
------------

AeoLiS writes its simulation results to a NetCDF4 file (by default named ``aeolis.nc``). The outputed variables (``output_vars``) and frequency (``output_times``) of this file are configured via the configuration file. You can request any of the variables listed in the :ref:`Model state/output` section, while ``x`` and ``y`` are automatically included.

While most output variables have 2D spatial dimensions combined with a ``time`` dimension (time, ny, nx), some contain additional dimensions. For example, sediment mass (``mass``) is calculated per grid cell, bed layer, and sediment fraction. Outputting such a five-dimensional variable can result in large file sizes. To reduce this, AeoLiS offers ``masstop`` as output, which provides the sediment distribution only for the active top layer.

By default, the output represents the model state at the exact moment defined by the output interval. If the internal time step (``dt``) is smaller than this interval, intermediate calculations are not saved. To evaluate variable behavior between output intervals, statistical summaries can be requested by appending suffixes directly to the variable names in ``output_vars``. Available suffixes include ``_avg`` (average), ``_sum`` (cumulative sum), ``_var`` (variance), ``_min`` (minimum), and ``_max`` (maximum) over the output interval.

.. note::
   **Known Issue:** Because the model parses the underscore (``_``) to identify these statistical requests, variables that contain underscores in their base names will cause parsing conflicts.

Example Python script 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Below is a Python script showing how to read the NetCDF output data en create a single plot.

.. code-block:: python

    import netCDF4 as nc
    import matplotlib.pyplot as plt

    # 1. Open the NetCDF output file
    ncfile = 'aeolis.nc'
    ds = nc.Dataset(ncfile, 'r')

    # 2. Load spatial coordinates and bed level
    x = ds.variables['x'][:, :]
    y = ds.variables['y'][:, :]
    zb_final = ds.variables['zb'][-1, :, :] # (time, y, x)
    ds.close()

    # 3. Plot the final bed level
    fig, ax = plt.subplots(figsize=(8, 6))
    pc = ax.pcolormesh(x, y, zb_final, cmap='viridis', shading='auto')
    ax.set_aspect('equal') 
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.set_title('Final Bed Level (zb)')
    fig.colorbar(pc, ax=ax, label='Bed level (m)')

    plt.tight_layout()
    plt.show()


Example Python script for animation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Below is a second example demonstrating how to animate time-series output of the cross-shore (``ustars``) and alongshore (``ustarn``) shear velocity components.

.. code-block:: python

    import netCDF4 as nc
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation

    # 1. Open the NetCDF output file and load variables
    ncfile = 'aeolis.nc'
    ds = nc.Dataset(ncfile, 'r')

    time = ds.variables['time'][:]
    x = ds.variables['x'][:, :]
    y = ds.variables['y'][:, :]
    ustar = ds.variables['ustar'][:] # (time, y, x)
    ustars = ds.variables['ustars'][:] # (time, y, x)
    ustarn = ds.variables['ustarn'][:] # (time, y, x)
    ds.close()

    # 2. Set up the figure layout
    fig, ax = plt.subplots(figsize=(8, 6))
    title = ax.set_title(f'Shear Velocity at t = {time[0]:.0f} s')
    ax.set_aspect('equal')
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')

    # Initial background mesh, colorbar and overlaying quiver plot (vectors)
    pc = ax.pcolormesh(x, y, ustar[0, :, :], cmap='YlOrRd', shading='auto', vmin=0, vmax=np.max(ustar))
    fig.colorbar(pc, ax=ax, label='Shear velocity magnitude (m/s)')
    Q = ax.quiver(x, y, ustars[0, :, :], ustarn[0, :, :], color='black')
    
    # 3. Define the update function for the animation
    def update(frame):
        pc.set_array(ustar_mag[frame, :, :].ravel())                 # update basemap
        Q.set_UVC(ustars[frame, :, :], ustarn[frame, :, :])          # update quiver
        title.set_text(f'Shear Velocity at t = {time[frame]:.0f} s') # update title
        return pc, Q, title

    # 4. Create the animation
    ani = animation.FuncAnimation(fig, update, frames=len(time), blit=False)
    ani.save('ustar_animation.mp4', writer='ffmpeg', fps=10)



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


Guidance on schematization: advection, shear and grainspeed
---------------------------------------------------

AeoLiS offers several options for schematizing spatial shear and sediment velocity. Higher complexity yields more realistic physics but increases computational cost. 

To illustrate the impact of these choices, we use a consistent demonstration case: A cone-shaped bedform sits on a non-erodible layer, upwind of a vegetated patch designed to capture all incoming sediment. We expect the landform to migrate downwind and evolve into a crescentic barchan dune. We also track the amount of deposition in the vegetated area to monitor mass balance.

Advection solver
^^^^^^^^^^^^^^^^
The advection equation is the core of sediment transport computations in AeoLiS. It balances the spatial change in sediment transport with the pickup or deposition rate:

.. math::

   \frac{\partial C U_s}{\partial x} = \frac{C_{\mathrm{sat}} - C}{T}

Where :math:`C` is the actual sediment concentration, :math:`U_s` is the sediment velocity, :math:`x` is the distance in the transport direction, :math:`C_{\mathrm{sat}}` is the equilibrium sediment concentration, and :math:`T` is the adaptation time scale. 

You can solve this using three different methods (via the ``solver`` keyword):

* **steadystate (default):** Highly recommended for most practical cases. Assumes concentration does not change over time within a single timestep (:math:`dc/dt = 0`) and solves rapidly using a finite difference sweeping algorithm.
* **euler_backward:** An implicit, time-varying solver (formerly "pieter"). Recommended primarily for legacy purposes or when the steady-state assumption is problematic (e.g., :math:`dt < 10` s).
* **euler_forward:** A simple explicit solver. Extremely slow due to strict CFL stability conditions. Included mainly for testing or educational purposes.

Shear and grainspeed schematization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The magnitude of transport is dictated by the equilibrium concentration :math:`C_{\mathrm{sat}}`, which is a function of the shear velocity :math:`u_*`:

.. math::

   C_{\mathrm{sat}} \propto (u_* - u_{th})^3

Where :math:`u_{th}` is the threshold shear velocity. The calculation of sediment velocity (:math:`u_s`) and the spatial variation of shear velocity (:math:`u_*`) are controlled by the ``method_grainspeed`` and ``process_shear`` keywords.

Case 0: Flat Conditions (Original Method)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* **Configuration:** ``process_shear = F``, ``method_grainspeed = windspeed``

Topographic feedback is disabled (:math:`\nabla \tau = 0`), meaning the wind blows over the cone as if it were flat. The sediment velocity simply equals the wind speed (:math:`u_s = u_w`). 

**Result:** No landform development occurs. Sediment is stripped from the cone and deposited directly into the vegetation. This is the fastest method (20 mins), but should only be used for bulk transport calculations where morphodynamics are irrelevant.

.. video:: images/case_00.mp4
   :autoplay:
   :loop:
   :muted:
   :width: 800

Case 1: Uniform Grain Speed
~~~~~~~~~~~~~~~~~~~~~~~~~~~
* **Configuration:** ``process_shear = F``, ``method_grainspeed = duran_uniform``

Sediment in saltation moves significantly slower than the wind. This case introduces the saltation model by Durán (2007) so that :math:`u_s < u_w`, but keeps the flat-bed assumption for shear velocity. 

**Result:** Deposition patterns become more localized, but the landform still does not migrate because there is no topographic steering. It is slightly slower (33 mins) but physically more realistic than Case 0 for static topographies.

.. video:: images/case_01.mp4
   :autoplay:
   :loop:
   :muted:
   :width: 800

    <video autoplay loop muted playsinline width="100%">
        <source src="../_images/case_01.mp4" type="video/mp4">
    </video>

Case 2: Topographic Steering
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* **Configuration:** ``process_shear = T``, ``method_grainspeed = duran_uniform``

The shear velocity vector is now perturbed by the topography, creating spatial gradients in transport capacity. The sediment velocity :math:`u_s` remains uniform.

**Result:** The variation in shear velocity causes the entire landform to migrate towards the vegetation. However, without flow separation or spatial variation in grain speed, it fails to evolve into a crescentic shape. Computational cost increases significantly (2:01 hrs) due to the secondary rotational grid required to calculate shear.

.. video:: images/case_02.mp4
   :autoplay:
   :loop:
   :muted:
   :width: 800

Case 3 (special): Separation Bubble
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* **Configuration:** ``process_separation = T`` (Requires ``process_shear = T``)

Steep lee-side slopes cause airflow to separate, creating a zone of recirculation and low shear. 

**Result:** Activating the separation bubble prevents sediment transport on the lee side, preserving the steep slip face and allowing the crescentic barchan shape to form (2:10 hrs). Note: In highly complex topographies (like dense vegetation), the bubble may produce undesirable morphodynamics, so use it judiciously.

.. video:: images/case_03.mp4
   :autoplay:
   :loop:
   :muted:
   :width: 800

Case 4: Spatially Varying Grain Speed (Analytical)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* **Configuration:** ``process_shear = T``, ``method_grainspeed = duran``

For a complete description of transport around landforms, the sediment velocity :math:`u_s` must also vary spatially. This method uses Durán's analytical approximation, which incorporates slope terms but assumes slopes are relatively gentle to avoid heavy numerical solving.

**Result:** Produces the expected crescentic morphodynamics efficiently (2:03 hrs). **This is the recommended configuration** for most simulations involving bedform evolution where topographic steering is important.

.. video:: images/case_04.mp4
   :autoplay:
   :loop:
   :muted:
   :width: 800

Case 5: Steep Slopes (Numerical Solution)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* **Configuration:** ``process_shear = T``, ``method_grainspeed = duran_full``

The analytical approximation in Case 4 can overestimate upslope transport on very steep slopes (e.g., :math:`> 33^\circ`). This method solves the full sediment velocity equation numerically to account for strong gravitational effects.

**Result:** While differences are subtle in standard dune simulations, this method is crucial for extreme topography like steep blowout cliffs (2:09 hrs). *Caution: Dynamic avalanching is not fully coupled yet, so schematize steep slopes carefully.*

.. video:: images/case_05.mp4
   :autoplay:
   :loop:
   :muted:
   :width: 800

Summary of Parameter Settings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table:: 
   :widths: 5 25 10 10 15 15 10 10
   :header-rows: 1

   * - Case
     - Description
     - ``shear``
     - ``sep.``
     - ``grainspeed``
     - Application
     - Time
     - Dep. (m³)
   * - **0**
     - Flat / Windspeed
     - ``F``
     - ``F``
     - ``windspeed``
     - Bulk transport
     - 0:20
     - 9961
   * - **1**
     - Flat / Grainspeed
     - ``F``
     - ``F``
     - ``duran_uniform``
     - Static topography
     - 0:33
     - 11095
   * - **2**
     - Topo Steering
     - ``T``
     - ``F``
     - ``duran_uniform``
     - Veg-dominated dunes
     - 2:01
     - 13156
   * - **3**
     - Separation Bubble
     - ``T``
     - ``T``
     - ``duran_uniform``
     - Moderate morphology
     - 2:10
     - 11222
   * - **4**
     - Varying $u_s$ (Approx)
     - ``T``
     - ``T``
     - ``duran``
     - **Landform evolution**
     - 2:03
     - 10313
   * - **5**
     - Varying $u_s$ (Full)
     - ``T``
     - ``T``
     - ``duran_full``
     - Extreme topography
     - 2:09
     - 10406

*Note: Deposition tracked out of a total 13,404 m³.*
