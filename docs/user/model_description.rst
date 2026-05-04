.. _model_description:

Model description
=================

Quick Overview
--------------

This section provides a summary of the main processes, equations, and configuration parameters in AeoLiS. For more information, refer to the detailed sections linked in the text (or scroll further down this page).

For guidance on setting up an AeoLiS model, see the :ref:`model in- & output guide <model-input-output>`. The main configuration file (default: ``aeolis.txt``) is the basis of the model setup and contains all parameter settings and process-flags, and serves as the central reference for other input files. The computational domain is constructed using x- and y-coordinates (``xgrid_file``, ``ygrid_file``) alongside the initial bed elevation (``bed_file``). External environmental forcing is defined through continuous time series of wind (``wind_file``), water levels (``tide_file``), and waves (``wave_file``).

The simulation advances sequentially through time steps, repeating all activated processes and continuously updating the morphological model state. The simulation duration runs from a defined start time (``tstart``) to an end time (``tstop``), both specified in seconds relative to a designated reference date (``refdate``). A typical internal time step (``dt``) is 3600 seconds (1 hour). As the model progresses, it exports user-defined variables (``output_vars``) to a NetCDF file (default: ``aeolis.nc``, defined by ``output_file``) at customized intervals (``output_times``).

.. note:: 
   In the AeoLiS source code, all model parameters are stored in two dictionaries for BMI-compatibility. All parameters defined on the computational domain, i.e., with dimension (``ny``,``nx``), are stored in the ``s``-dictionary (e.g., ``s['x']``), while all single parameters are stored in the ``p``-dictionary (e.g., ``p['dt']``).

Sediment Transport
^^^^^^^^^^^^^^^^^^^
For detailed information, see the :ref:`sediment transport section <aeolian-sediment-transport>`.

Aeolian sediment transport is the core of the AeoLiS model. It is computed using a two-dimensional advection scheme, simplified here for one-dimensional transport of a single sediment fraction:

.. math::
   :label: advection_overview
           
   \frac{\partial c}{\partial t} + u_{\mathrm{sed}} \frac{\partial c}{\partial x} = E - D = \min \left ( \frac{\partial m_{\mathrm{a}}}{\partial t} \quad ; \quad \frac{c_{\mathrm{sat}} - c}{T} \right )

The right-hand side of the advection equation represents the net entrainment (``pickup``); the difference between erosion :math:`E` and deposition :math:`D` [:math:`\mathrm{kg/m^2/s}`]. The saturated sediment concentration :math:`c_{\mathrm{sat}}` (``Cu``) [:math:`\mathrm{kg/m^2}`] defines the transport capacity, while :math:`c` (``Ct``) [:math:`\mathrm{kg/m^2}`] is the instantaneous concentration in the air. Transport is activated in the configuration file using ``process_transport``. This net entrainment is governed by the adaptation timescale :math:`T` (``T``) [:math:`\mathrm{s}`], which determines how quickly the concentration reaches equilibrium. To allow sediment to actually erode from or deposit to the bed, ``process_bedupdate`` must be enabled. 

Solving this advection equation is one of the most computationally expensive parts of the model. You can choose different numerical approaches using the ``solver`` keyword. For detailed guidance on these options, see the :ref:`solver guide <solver-guide>`.

Several methods are available to compute the saturated sediment concentration (``method_transport``). The equation by :cite:`Bagnold1937a` (``bagnold``) is the default:

.. math::
   :label: bagnold_overview

   c_{\mathrm{sat}} = \max \left ( 0 \quad ; \quad C \frac{\rho_{\mathrm{a}}}{g} \sqrt{\frac{d_{n}}{D_{n}}} \frac{\left ( u_* - u_{\mathrm{th}} \right )^3}{u_{\mathrm{sed}}} \right )

The sediment velocity :math:`u_{\mathrm{sed}}` (``u``, ``us``, ``un``) [:math:`\mathrm{m/s}`] is determined by the ``method_grainspeed`` parameter. It can either be set equal to the governing wind speed (``windspeed``) or calculated using a saltation model (e.g., ``duran``). For more information on grain speed computations, see the :ref:`sediment velocity <sediment-velocity>` section.

.. note:: 
   For all vector variables (like ``uw``, ``ustar``,  ``tau``, ``u``,  ``q``), the subscripts ``s`` and ``n`` (e.g., ``uws``, ``uwn``) indicate the cross-shore and longshore directions, respectively, and the name without a subscript represents the overall magnitude. 

AeoLiS supports the inclusion of multiple sediment fractions (``grain_size``, ``grain_dist``) across multiple vertical layers (``nlayers``). This allows for the simulation of sediment sorting, mixing, and armoring. More details are provided in the :ref:`multi-fraction sediment transport <multi-fraction-sediment-transport>` section.

Enabling ``process_bedinteraction`` incorporates a bed interaction parameter into the advection equation. For more information, see the :ref:`bed interaction approach <bed-interaction-approach>` section.

Wind and Shear Velocity
^^^^^^^^^^^^^^^^^^^^^^^^
Detailed section: :ref:`wind-shear-velocity`

The shear velocity :math:`u_*` (``ustar``) [:math:`\mathrm{m/s}`] acts as the primary driver of sediment transport. It is initially computed for a flat bed using the Prandtl-Von Kármán Law of the Wall, based on the wind velocity :math:`u_w` (``uw``) [:math:`\mathrm{m/s}`] at a given elevation. Wind conditions are provided via the ``wind_file`` and the computation is enabled via ``process_wind``.

.. math::
   :label: lawofwall_overview

   u_* = \frac{u_w}{\ln \left( \frac{z}{z_0} \right)}\kappa

Topography can steer the wind, causing perturbations in the shear stress :math:`\tau` (``tau``) [:math:`\mathrm{N/m^2}`], where :math:`\tau = \rho_a u_*^2`:

.. math::
   :label: topo_steering_overview

   \vec{\tau}(x,y) = \vec{\tau}_{0} + |\vec{\tau}_{0}|\delta\vec{\tau}(x,y)

This steering process can be activated using the ``process_shear`` keyword. Different methods are available to compute these shear perturbations, which can be selected through ``method_shear``. More information on these computations is given in the :ref:`topographic steering section <topographic-steering>`.

The presence of vegetation can also reduce the effective shear stress. For more information, see the :ref:`vegetation documentation <vegetation>`. 

Shear Velocity Threshold
^^^^^^^^^^^^^^^^^^^^^^^^^
Detailed section: :ref:`shear-velocity-threshold`

Where shear velocity drives transport, the threshold velocity :math:`u_{\mathrm{th}}` (``uth``) [:math:`\mathrm{m/s}`] serves as a supply-limiter. It acts as a collective parameter for all supply-limiting processes, scaling the base threshold :math:`u_{\mathrm{*th,0}}` (``uth0``) [:math:`\mathrm{m/s}`] by various environmental factors. Threshold calculations are enabled via ``process_threshold``.

.. math::
  :label: threshold_overview
  
  u_{\mathrm{* th}} = u_{\mathrm{* th, 0}} \cdot f_{\mathrm{M}} \cdot f_{\mathrm{R}} \cdot f_{\mathrm{S}}

The base threshold :math:`u_{\mathrm{* th, 0}}` [:math:`\mathrm{m/s}`] is computed based on the local grain size [:math:`\mathrm{m}`] and density [:math:`\mathrm{kg/m^3}`] (activated via ``th_grainsize``). This base value is then scaled by supply-limiting factors depending on the enabled model processes. The influence of surface moisture (:math:`f_{\mathrm{M}}`) [:math:`\mathrm{-}`] is activated via ``th_moisture``, while sheltering by non-erodible roughness elements (:math:`f_{\mathrm{R}}`) [:math:`\mathrm{-}`] is configured via ``th_sheltering``. The restricting effect of a non-erodible layer can be included using ``th_nelayer``.

Vegetation
^^^^^^^^^^^
Detailed section: :ref:`vegetation`

Vegetation plays a key role in the development of dunes. The simulation of vegetation growth and spreading, along with the subsequent reduction in shear stress, can be enabled via ``process_vegetation``. The specific vegetation formulation is selected using ``method_vegetation`` (``duran`` for the original implementation, or ``grass`` for the newer framework by :cite:t:`vanWesten2026`). 

In the original description, the vegetation density :math:`\rho_{\mathrm{veg}}` (``rhoveg``) [:math:`\mathrm{-}`] is computed based on the current vegetation height :math:`h_{\mathrm{veg}}` (``hveg``) [:math:`\mathrm{m}`] relative to the maximum attainable height :math:`H_{\mathrm{max}}` (``hveg_max`` or ``Hveg``) [:math:`\mathrm{m}`]:

.. math::
   :label: rhoveg_overview

   \rho_{\mathrm{veg}} = (\frac{h_{\mathrm{veg}}}{H_{\mathrm{max}}})^2

This density determines the magnitude of the shear stress reduction acting on the sand bed, which relies on a vegetation-related roughness parameter :math:`\Gamma` (``gamma_vegshear``) [:math:`\mathrm{-}`] and the basal cover :math:`\rho_{\mathrm{veg}}` (``rhoveg``) [:math:`\mathrm{-}`]:

.. math::
   :label: raupach_overview

   u_{*}= \frac{u_*}{\sqrt{1 + \Gamma \rho_{\mathrm{veg}}}}

The vertical development of the vegetation over time is described by:

.. math::
   :label: dhveg_overview

   \frac{\partial h_{\mathrm{veg}}}{\partial t} = V_{\mathrm{ver}} \left( 1 - \frac{h_{\mathrm{veg}}}{H_{\mathrm{max}}} \right) - \gamma_{\mathrm{veg}} |\Delta z_{\mathrm{burial}}|

This growth is governed by the intrinsic vertical growth rate :math:`V_{\mathrm{ver}}` (``V_ver``) [:math:`\mathrm{m/s}`] and the plant's sensitivity to sediment burial or erosion :math:`\gamma_{\mathrm{veg}}` (``veg_gamma``) [:math:`\mathrm{-}`]. 


Hydrodynamics and Surface Moisture
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Detailed section: :ref:`hydrodynamics-moisture`

Water levels, wave runup, and groundwater can wet the beach, temporarily increasing the shear velocity threshold and mixing sediment fractions. The model reads input water levels (``tide_file``) and wave heights (``wave_file``), which are activated via ``process_tide`` and ``process_wave``. 

The input water level is first projected onto the domain to establish the Still Water Level (``SWL``) [:math:`\mathrm{m}`]. If wave runup is enabled (``process_runup``), the runup height (``R``) [:math:`\mathrm{m}`] is computed and added to form the Total Water Level (``TWL``) [:math:`\mathrm{m}`], where TWL = SWL + R. The local water surface elevation (``zs``) [:math:`\mathrm{m}`] is then determined as the maximum of the bed level (``zb``) [:math:`\mathrm{m}`] and the TWL, from which the actual water depth (``hw``) [:math:`\mathrm{m}`] is derived. Spatial masks (``tide_mask``, ``wave_mask``, ``runup_mask``) can be applied to restrict or modify where these hydrodynamics act.

Inundation wets the bed, increasing the surface moisture (``moist``) [:math:`\mathrm{-}`]. Once the beach is exposed, infiltration and evaporation gradually dry the surface. This moisture tracking is activated via ``process_moist``. A more advanced description of intertidal groundwater fluctuations by :cite:t:`Hallin2023` can also be enabled through ``process_groundwater``. For more information on these computations, see the :ref:`surface moisture section <surface-moisture>`.

Furthermore, wave impacts can mix multiple sediment fractions across several bed layers down to the depth of disturbance. This mixing process is enabled via ``process_mixtoplayer``. For more details on how waves rework the bed, see the :ref:`hydraulic sediment mixing section <hydraulic-sediment-mixing>`.

Morphological Change
^^^^^^^^^^^^^^^^^^^^^
Detailed section: :ref:`morphological-change`

Gradients in aeolian sediment transport result in net erosion or deposition, causing the bed level :math:`z_B` (``zb``) [:math:`\mathrm{m}`] to change over time; :math:`\Delta z_B` (``dzb``) [:math:`\mathrm{m}`]. This morphological updating is enabled by ``process_bedupdate``:

.. math::
   :label: bedupdate_overview

   \frac{\partial z}{\partial t} = - \frac{1}{\rho_{\mathrm{sed}}(1 - p)} (E - D)

This change is driven directly by the net entrainment :math:`(E - D)` (``pickup``) [:math:`\mathrm{kg/m^2/s}`] computed in the advection equation, scaled by the sediment density :math:`\rho_{\mathrm{sed}}` (``rhog``) [:math:`\mathrm{kg/m^3}`] and the sediment porosity :math:`p` (``porosity``) [:math:`\mathrm{-}`]. Together, the density and porosity represent the bulk density of the bed. For more detailed mechanics on this mass balance, see the :ref:`morphological change section <morphological-change>`.

To redistribute sediment when the local slope becomes too steep, avalanching can be enabled through ``process_avalanche``. This routine triggers when the bed slope exceeds the static angle of repose (``theta_stat``) [:math:`\mathrm{^\circ}`] and relaxes the slope back to the dynamic angle of repose (``theta_dyn``) [:math:`\mathrm{^\circ}`].

.. _fig-aeolis-overview:

.. figure:: /images/aeolis_overview.png
   :width: 900px
   :align: center
   
   Overview of the AeoLiS model structure and simulated processes.


.. _aeolian-sediment-transport:
Aeolian Sediment Transport
--------------------------

Calculating aeolian sediment transport is the core of the AeoLiS model. 
The model is based on the approach of :cite:`deVries2014a` which is extended to compute the
spatiotemporal varying sediment availability through simulation of the
process of beach armoring. For this purpose the bed is discretized in
horizontal grid cells and in vertical bed layers (``nlayers``) [:math:`\mathrm{-}`]. Moreover, the
grain size distribution is discretized into fractions (``grain_size`` [:math:`\mathrm{m}`], ``grain_dist`` [:math:`\mathrm{-}`]). This allows the
grain size distribution to vary both horizontally and vertically. A
bed composition module is used to compute the sediment availability
for each sediment fraction individually. This model approach is a
generalization of existing model concepts, like the shear velocity
threshold and critical fetch, and therefore compatible with these
existing concepts.


.. _advection-equation:
Advection Equation
^^^^^^^^^^^^^^^^^^

A 1D advection scheme is adopted in correspondence with
:cite:`deVries2014a` in which :math:`c` (``Ct``) [:math:`\mathrm{kg/m^2}`] is
the instantaneous sediment mass per unit area in transport:

.. math::
   :label: advection
           
   \frac{\partial c}{\partial t} + u_{\mathrm{sed}} \frac{\partial c}{\partial x} = E - D

Here, :math:`t` (``_time``) [:math:`\mathrm{s}`] denotes time, :math:`x` (``x``) [:math:`\mathrm{m}`] denotes the cross-shore
distance, and :math:`u_{\mathrm{sed}}` (``u``, ``us``, ``un``) [:math:`\mathrm{m/s}`] is the sediment velocity. :math:`E` and :math:`D`
[:math:`\mathrm{kg/m^2/s}`] represent the erosion and deposition terms
and hence combined represent the net entrainment of sediment (``pickup``). 

.. note:: 
   Equation :eq:`advection` differs from Equation 9 in
   :cite:`deVries2014a` as they use the saltation height :math:`h` [:math:`\mathrm{m}`]
   and the volumetric sediment concentration :math:`C_{\mathrm{c}}`
   [:math:`\mathrm{kg/m^3}`]. As :math:`h` is not solved for, the
   presented model computes the sediment mass per unit area :math:`c = h
   C_{\mathrm{c}}` rather than the sediment concentration
   :math:`C_{\mathrm{c}}`. For conciseness we still refer to :math:`c` as
   the *sediment concentration*.

The net entrainment is determined based on a balance between the
equilibrium or saturated sediment concentration
:math:`c_{\mathrm{sat}}` (``Cu``) [:math:`\mathrm{kg/m^2}`] and the
instantaneous sediment transport concentration :math:`c`, and is
maximized by the available sediment in the bed :math:`m_{\mathrm{a}}` (``mass``)
[:math:`\mathrm{kg/m^2}`] according to:

.. math::
   :label: erodep
           
   E - D = \min \left ( \frac{\partial m_{\mathrm{a}}}{\partial t} \quad ; \quad \frac{c_{\mathrm{sat}} - c}{T} \right )

:math:`T` (``T``) [:math:`\mathrm{s}`] represents an adaptation time scale that is assumed
to be equal for both erosion and deposition. A time scale of 1 second
is commonly used :cite:`deVries2014a`.

Solving this advection equation is one of the most computationally expensive processes in the model. The numerical approach can be selected using the ``solver`` parameter, which provides three options: ``steadystate``, ``euler_backward``, and ``euler_forward``. 

.. tip::
   The ``steadystate`` solver is the default and recommended option for most applications. It is the most thoroughly tested and fastest option available. The steady-state assumption is generally valid for aeolian transport because wind and sediment transport adjust to local conditions on a scale of seconds to minutes, which is orders of magnitude faster than the timescale of typical timesteps in the AeoLiS model (hours). For detailed guidance on configuring these numerical approaches, see the :ref:`solver guide <solver-guide>`.


.. _saturated-sediment-transport:

Saturated Sediment Transport
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The equilibrium, or saturated, sediment concentration :math:`c_{\mathrm{sat}}` (``Cu``) [:math:`\mathrm{kg/m^2}`] is computed using an empirical sediment transport formulation selected via the ``method_transport`` parameter. The default formulation is based on Bagnold (:cite:`Bagnold1937a`) via the ``bagnold`` setting:

.. math::
   :label: bagnold_qsat

   q_{\mathrm{sat}} = C_{\mathrm{b}} \frac{\rho_{\mathrm{a}}}{g} \left( u_* - u_{\mathrm{th}} \right)^3

in which :math:`q_{\mathrm{sat}}` [:math:`\mathrm{kg/m/s}`] is the saturated sediment transport rate representing the sediment transport capacity. :math:`u_*` (``ustar``) [:math:`\mathrm{m/s}`] is the shear velocity, and :math:`u_{\mathrm{th}}` (``uth``) [:math:`\mathrm{m/s}`] is the velocity threshold. The properties of the sediment and air are represented by a series of parameters: :math:`C_{\mathrm{b}}` (``Cb``) [:math:`\mathrm{-}`] is an empirical constant, :math:`\rho_{\mathrm{a}}` (``rhoa``) [:math:`\mathrm{kg/m^3}`] is the density of the air, and :math:`g` (``g``) [:math:`\mathrm{m/s^2}`] is the gravitational constant.

The saturated sediment concentration :math:`c_{\mathrm{sat}}` is directly related to the saturated sediment transport rate :math:`q_{\mathrm{sat}}` and the sediment velocity :math:`u_{\mathrm{sed}}` (``u``) [:math:`\mathrm{m/s}`] through the relationship :math:`q_{\mathrm{sat}} = c_{\mathrm{sat}} \cdot u_{\mathrm{sed}}`. Therefore, to obtain the mass per unit area, the transport rate is divided by the sediment velocity:

.. math::
   :label: bagnold_csat

   c_{\mathrm{sat}} = \max \left( 0 \quad ; \quad C_{\mathrm{b}} \frac{\rho_{\mathrm{a}}}{g} \frac{\left( u_* - u_{\mathrm{th}} \right)^3}{u_{\mathrm{sed}}} \right)

Depending on the configuration, several other formulations can be selected to compute the saturated sediment concentration, for example:

``kawamura``:

  .. math::
     :label: kawamura

     c_{\mathrm{sat}} = C_{\mathrm{k}} \frac{\rho_{\mathrm{a}}}{g} \frac{\left( u_* + u_{\mathrm{th}} \right)^2 \left( u_* - u_{\mathrm{th}} \right)}{u_{\mathrm{sed}}}

``lettau``:

  .. math::
     :label: lettau

     c_{\mathrm{sat}} = C_{\mathrm{l}} \frac{\rho_{\mathrm{a}}}{g} \frac{u_*^2 \left( u_* - u_{\mathrm{th}} \right)}{u_{\mathrm{sed}}}

``dk``:

  .. math::
     :label: dk

     c_{\mathrm{sat}} = C_{\mathrm{dk}} \frac{\rho_{\mathrm{a}}}{g} \frac{0.8 u_{\mathrm{th}} \left( u_*^2 - \left( 0.8 u_{\mathrm{th}} \right)^2 \right)}{u_{\mathrm{sed}}}


.. _sediment-velocity:

Sediment Transport Velocity
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The horizontal sediment velocity :math:`u_{\mathrm{sed}}` (``u``) [:math:`\mathrm{m/s}`] determines the advection speed of the sediment concentration. It can be computed using different approaches selected via the ``method_grainspeed`` parameter. For detailed guidance on selecting the appropriate method and its impact on computational time and landform evolution, see :ref:`this guide <solver-guide>`.

The simplest approach (``windspeed``) assumes the horizontal sediment velocity is equal to the wind velocity :math:`u_{\mathrm{w}}` (``uw``) [:math:`\mathrm{m/s}`]:

.. math::
   :label: used_windspeed

   u_{\mathrm{sed}} = u_{\mathrm{w}}

While this method is the fastest and most robust, it overpredicts the horizontal sediment velocity because grains in saltation move slower than the wind. Consequently, it fails to capture deposition patterns, making it suitable only for bulk transport calculations.

Predictions of the actual saltation velocity provide a more realistic description of horizontal sediment movement :cite:`sauermann2001continuum`. The sediment velocity can be determined from a momentum balance :cite:`duran2007thesis` consisting of three terms: the drag force acting on the grains, the loss of momentum during grain-bed interaction (splashing), and the downhill gravity force:

.. math::
   :label: used_momentum
   
   \frac{(\vec{v}_{\mathrm{eff}} - \vec{u}_{\mathrm{sed}})|\vec{v}_{\mathrm{eff}} - \vec{u}_{\mathrm{sed}}|}{u_{\mathrm{f}}^2} - \frac{\vec{u}_{\mathrm{sed}}}{2 \alpha |\vec{u}_{\mathrm{sed}}|} - \vec{\nabla z_{\mathrm{B}}} = 0

where :math:`v_{\mathrm{eff}}` [:math:`\mathrm{m/s}`] is the effective wind velocity driving the grains, which depends on the shear velocity :math:`u_*` and the threshold shear velocity :math:`u_{\mathrm{*th}}`. :math:`u_{\mathrm{f}}` [:math:`\mathrm{m/s}`] is the fluid threshold velocity (or grain settling velocity), :math:`\nabla z_{\mathrm{B}}` [:math:`\mathrm{-}`] is the bed slope, and :math:`\alpha` [:math:`\mathrm{-}`] is an effective restitution coefficient for the grain-bed interaction (e.g., :math:`\alpha = 0.42` for :math:`d = 250` :math:`\mathrm{\mu m}`). 

.. note::
   The computed :math:`u_{\mathrm{sed}}` represents the collective horizontal sediment movement, not the velocity of individual grains. 

AeoLiS provides three options based on this momentum balance:

``duran_full`` solves the full momentum balance equation (Equation :eq:`used_momentum`) numerically for each timestep. This method is computationally heavy but necessary for accurate sediment velocities on steep topography (e.g., steep blowout cliffs).

``duran`` uses an analytical approximation of the full momentum balance. It accounts for spatial variations and slope effects by assuming slopes are relatively gentle to avoid computationally expensive numerical solving:

  .. math::
     :label: used_duran_approx

     u_{\mathrm{sed}} \approx \left( v_{\mathrm{eff}} - \frac{u_{\mathrm{f}}}{\sqrt{2\alpha A}} \right) \hat{e}_{\tau} - \frac{\sqrt{2\alpha} u_{\mathrm{f}}}{A} \nabla z_{\mathrm{B}}

  where :math:`\hat{e}_{\tau}` is the wind direction unit vector and :math:`A \equiv |\hat{e}_{\tau} + 2\alpha \nabla z_{\mathrm{B}}|`. The first term points toward the wind direction, while the second is directed along the surface gradient, accounting for the competing effects of wind and gravity.

``duran_uniform`` assumes a flat bed (:math:`\nabla z_{\mathrm{B}} = 0`). The sediment velocity is spatially uniform, simplifying the equation further to:

  .. math::
     :label: used_duran_uniform

     u_{\mathrm{sed}} = v_{\mathrm{eff}} - \frac{u_{\mathrm{f}}}{\sqrt{2\alpha}}

.. tip::
   **Guidance on selecting a grain speed method (see :ref:`this guide <solver-guide>` for more information)**

   * Use ``duran`` for most simulations involving landform evolution where topographic steering is important.
   * Use ``duran_full`` only for simulations involving topography with very steep gradients (e.g., blowout cliffs).
   * Use ``duran_uniform`` for static topography where you need realistic deposition patterns but no landform migration.
   * Use ``windspeed`` only for basic, bulk transport calculations where morphodynamics are irrelevant.



.. _multi-fraction-sediment-transport:
Multi-fraction Sediment Transport
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The formulation for the equilibrium or saturated sediment
concentration :math:`c_{\mathrm{sat}}` (Equation
:eq:`equilibrium-conc`) is capable of dealing with variations in
grain size through the variables :math:`u_{\mathrm{th}}`,
:math:`d_{\mathrm{n}}` and :math:`C` :cite:`Bagnold1937a`. However,
the transport formulation only describes the saturated sediment
concentration assuming a fixed grain size distribution, but does not
define how multiple fractions coexist in transport. If the saturated
sediment concentration formulation would be applied to each fraction
separately and summed up to a total transport, the total sediment
transport would increase with the number of sediment fractions. Since
this is unrealistic behavior the saturated sediment concentration
:math:`c_{\mathrm{sat}}` for the different fractions should be
weighted in order to obtain a realistic total sediment
transport. Equation :eq:`erodep` therefore is modified to include a
weighting factor :math:`\hat{w}_k` in which :math:`k` represents the
sediment fraction index:

.. math::
   :label: erodep_multi
           
   E_k - D_k = \min \left ( \frac{\partial m_{\mathrm{a},k}}{\partial t} \quad ; \quad \frac{\hat{w}_k \cdot c_{\mathrm{sat},k} - c_k}{T} \right )

It is common to use the grain size distribution in the bed as
weighting factor for the saturated sediment concentration
(e.g. :cite:`Delft3DManual`, section 11.6.4). Using the grain size
distribution at the bed surface as a weighting factor assumes, in case
of erosion, that all sediment at the bed surface is equally exposed to
the wind.

Using the grain size distribution at the bed surface as weighting
factor in case of deposition would lead to the behavior where
deposition becomes dependent on the bed composition. Alternatively, in
case of deposition, the saturated sediment concentration can be
weighted based on the grain size distribution in the air. Due to the
nature of saltation, in which continuous interaction with the bed
forms the saltation cascade, both the grain size distribution in the
bed and in the air are likely to contribute to the interaction between
sediment fractions. The ratio between both contributions in the model
is determined by a bed interaction parameter :math:`\zeta` (``bi``).

.. note::
   The bed interaction parameter (``bi``) described here is currently distinct from the recently introduced bed-interaction factor (``zeta``). While conceptually similar, they serve different purposes:

   * ``bi``: A static variable that weights the contribution of the bed composition (sand in the bed) versus the airborne composition (sand already in transport) when computing the transport rate for multiple sediment fractions.
   * ``zeta``: A dynamically computed factor based on surface properties that decouples the air and bed. It determines the extent to which transport is governed by bed conditions (supply-limited) versus airborne conditions (wind-driven capacity) (see :ref:`this section <bed-interaction-approach>`).

   Future updates should aim to consolidate these two parameters.


The weighting of erosion and deposition of individual fractions is
computed according to:

.. math::
   :label: weigh
   
   \begin{align}
     \hat{w}_k &= \frac{w_k}{ \sum_{k=1}^{n_{\mathrm{k}}}{w_k} } \\
     \mathrm{where} \quad w_k &= (1 - \zeta) \cdot w^{\mathrm{air}}_k + (1 - \hat{S}_k) \cdot w^{\mathrm{bed}}_k
   \end{align}

in which :math:`k` represents the sediment fraction index,
:math:`n_{\mathrm{k}}` the total number of sediment fractions, :math:`w_k` is the
unnormalized weighting factor for fraction :math:`k`, :math:`\hat{w}_k` is its
normalized counterpart, :math:`w^{\mathrm{air}}_k` and :math:`w^{\mathrm{bed}}_k`
are the weighting factors based on the grain size distribution in the
air and bed respectively and :math:`\hat{S}_k` is the effective sediment
saturation of the air. The weighting factors based on the grain size
distribution in the air and the bed are computed using mass ratios:

.. math::
   :label: weights
           
   w^{\mathrm{air}}_k = \frac{c_k}{c_{\mathrm{sat},k}} \quad ; \quad
   w^{\mathrm{bed}}_k = \frac{m_{\mathrm{a},k}}{\sum_{k=1}^{n_{\mathrm{k}}}{m_{\mathrm{a},k}}}

The sum of the ratio :math:`w^{\mathrm{air}}_k` over the fractions
denotes the degree of saturation of the air column for fraction
:math:`k`. The degree of saturation determines if erosion of a fraction may
occur. Also in saturated situations erosion of a sediment fraction can
occur due to an exchange of momentum between sediment fractions, which
is represented by the bed interaction parameter :math:`\zeta` (``bi``). The effective
degree of saturation is therefore also influenced by the bed
interaction parameter and defined as:

.. math::
   :label: saturation
   
   \hat{S}_k = \min \left ( 1 \quad ; \quad (1 - \zeta) \cdot \sum_{k=1}^{n_{\mathrm{k}}} w_k^{\mathrm{air}} \right )

When the effective saturation is greater than or equal to unity the
air is (over)saturated and no erosion will occur. The grain size
distribution in the bed is consequently less relevant and the second
term in Equation :eq:`weights` is thus minimized and zero in case
:math:`\zeta = 0`. In case the effective saturation is less than unity erosion
may occur and the grain size distribution of the bed also contributes
to the weighting over the sediment fractions. The weighting factors
for erosion are then composed from both the grain size distribution in
the air and the grain size distribution at the bed surface. Finally,
the resulting weighting factors are normalized to sum to unity over
all fractions (:math:`\hat{w}_k`).

The composition of weighting factors for erosion is based on the
saturation of the air column. The non-saturated fraction determines
the potential erosion of the bed. Therefore the non-saturated fraction
can be used to scale the grain size distribution in the bed in order
to combine it with the grain size distribution in the air according to
Equation :eq:`weights`. The non-saturated fraction of the air column
that can be used for scaling is therefore :math:`1 - \hat{S}_k`.

For example, if bed interaction is disabled (:math:`\zeta = 0`) and
the air is 70\% saturated, then the grain size distribution in the air
contributes 70\% to the weighting factors for erosion, while the grain
size distribution in the bed contributes the other 30\% (Figure
:numref:`fig-bed-interaction-parameter`, upper left panel). In case of
(over)saturation the grain size distribution in transport contributes
100\% to the weighting factors and the grain size distribution in the
bed is of no influence. Transport progresses in downwind direction
without interaction with the bed.

.. _fig-bed-interaction-parameter:

.. figure:: /images/bed_interaction_parameter.png
   :width: 600px
   :align: center

   Contributions of the grain size distribution in the bed and in the
   air to the weighting factors :math:`\hat{w}_k` for the equilibrium
   sediment concentration in Equation :eq:`erodep_multi` for different
   values of the bed interaction parameter.

To allow for bed interaction in saturated situations in which no net
erosion can occur, the bed interaction parameter :math:`\zeta` is used (Figure
:numref:`fig-bed-interaction-parameter`). The bed interaction parameter
can take values between 0.0 and 1.0 in which the weighting factors for
the equilibrium or saturated sediment concentration in an
(over)saturated situation are fully determined by the grain size
distribution in the bed or in the air respectively. A bed interaction
value of 0.2 represents the situation in which the grain size
distribution at the bed surface contributes 20\% to the weighting of
the saturated sediment concentration over the fractions. In the
example situation where the air is 70\% saturated such value for the
bed interaction parameter would lead to weighting factors that are
constituted for :math:`70\% \cdot (100\% - 20\%) = 56\%` based on the grain
size distribution in the air and for the other 44\% based on the grain
size distribution at the bed surface (Figure
:numref:`fig-bed-interaction-parameter`, upper right panel).

The parameterization of the exchange of momentum between sediment
fractions is an aspect of saltation that is still poorly
understood. Therefore calibration of the bed interaction parameter
:math:`\zeta` is necessary. The model parameters in Equation
:eq:`equilibrium-conc` can be chosen in accordance with the
assumptions underlying multi-fraction sediment transport. :math:`C` should
be set to 1.5 as each individual sediment fraction is well-sorted,
:math:`d_{\mathrm{n}}` should be chosen equal to :math:`D_{\mathrm{n}}` as the
grain size dependency is implemented through
:math:`u_{\mathrm{th}}`. :math:`u_{\mathrm{th}}` typically varies between 1 and 6
m/s for sand.

.. _sediment-sorting:
Sediment Sorting and Beach Armoring
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Since the equilibrium or saturated sediment concentration
:math:`c_{\mathrm{sat},k}` is weighted over multiple sediment fractions in
the extended advection model, also the instantaneous sediment
concentration :math:`c_k` is computed for each sediment fraction
individually. Consequently, grain size distributions may vary over the
model domain and in time. These variations are thereby not limited to
the horizontal, but may also vary over the vertical since fine
sediment may be deposited on top of coarse sediment or, reversely,
fines may be eroded from the bed surface leaving coarse sediment to
reside on top of the original mixed sediment. In order to allow the
model to simulate the processes of sediment sorting and beach armoring
the bed is discretized in horizontal grid cells and vertical bed
layers (2DV; Figure :numref:`fig-bedcomposition`).

The discretization of the bed consists of a minimum of three vertical
bed layers with a constant thickness and an unlimited number of
horizontal grid cells. The top layer is the *bed surface layer* and is
the only layer that interacts with the wind and hence determines the
spatiotemporal varying sediment availability and the contribution of
the grain size distribution in the bed to the weighting of the
saturated sediment concentration. One or more *bed composition layers*
are located underneath the bed surface layer and form the upper part
of the erodible bed. The bottom layer is the *base layer* and contains
an infinite amount of erodible sediment according to the initial grain
size distribution. The base layer cannot be eroded, but can supply
sediment to the other layers.

.. _fig-bedcomposition:

.. figure:: /images/bed_composition.png
   :align: center

   Schematic of bed composition discretisation and advection
   scheme. Horizontal exchange of sediment may occur solely through
   the air that interacts with the *bed surface layer*. The detail
   presents the simulation of sorting and beach armoring where the bed
   surface layer in the upwind grid cell becomes coarser due to
   non-uniform erosion over the sediment fractions, while the bed
   surface layer in the downwind grid cell becomes finer due to
   non-uniform deposition over the sediment fractions. Symbols refer
   to Equations :eq:`advection` and :eq:`erodep`.

Each layer in each grid cell describes a grain size distribution (``grain_size``, ``grain_dist``)) over
a predefined number of sediment fractions (``nfractions``) (Figure
:numref:`fig-bedcomposition`, detail). Sediment may enter or leave a
grid cell only through the bed surface layer. Since the velocity
threshold depends among others on the grain size, erosion from the bed
surface layer will not be uniform over all sediment fractions, but
will tend to erode fines more easily than coarse sediment (Figure
:numref:`fig-bedcomposition`, detail, upper left panel). If sediment
is eroded from the bed surface layer, the layer is repleted by
sediment from the lower bed composition layers. The repleted sediment
has a different grain size distribution than the sediment eroded from
the bed surface layer. If more fines are removed from the bed surface
layer in a grid cell than repleted, the median grain size
increases. If erosion of fines continues the bed surface layer becomes
increasingly coarse. Deposition of fines or erosion of coarse material
may resume the erosion of fines from the bed.

In case of deposition the process is similar. Sediment is deposited in
the bed surface layer that then passes its excess sediment to the
lower bed layers (Figure :numref:`fig-bedcomposition`, detail, upper
right panel). If more fines are deposited than passed to the lower bed
layers the bed surface layer becomes increasingly fine.

.. _hydraulic:
Hydraulic Mixing
~~~~~~~~~~~~~~~~~~~~

As sediment sorting due to aeolian processes can lead to armoring of a
beach surface, mixing of the beach surface or erosion of course
material may undo the effects of armoring. To ensure a proper balance
between processes that limit and enhance sediment availability in the
model both types of processes need to be sufficiently represented when
simulating spatiotemporal varying bed surface properties and sediment
availability.

A typical upwind boundary in coastal environments during onshore winds
is the water line. For aeolian sediment transport the water line is a
zero-transport boundary. In the presence of tides, the intertidal
beach is flooded periodically. Hydraulic processes like wave breaking
mix the bed surface layer of the intertidal beach, break the beach
armoring and thereby influence the availability of sediment. 

In the model the mixing of sediment is simulated by averaging the
sediment distribution over the depth of disturbance
(:math:`\Delta z_{\mathrm{d}}`). The depth of disturbance is linearly
related to the breaker height (e.g. :cite:`King1951`, :cite:`Williams1971`, :cite:`Masselink2007`). :cite:`Masselink2007` proposes an empirical factor
:math:`f_{\Delta z_{\mathrm{d}}}` (``facdod``) [-] that relates the depth of disturbance
directly to the local breaker height according to:

.. math::
   :label: disturbance_depth
   
   \Delta z_{\mathrm{d}} = f_{\Delta z_{\mathrm{d}}} \cdot \min \left ( H \quad ; \quad \gamma \cdot d \right )

in which :math:`\Delta z_{\mathrm{d}}` (``DOD``) [m] is the depth of disturbance, and the offshore wave height :math:`H` (``Hsmix``) [m] is taken as the
local wave height maximized by a maximum wave height over depth ratio
:math:`\gamma` (``gamma``) [-]. :math:`d` [m] is the water depth that is provided to the model
through an input time series of water levels. Typical values for
:math:`f_{\Delta z_{\mathrm{d}}}` are 0.05 to 0.4 and 0.5 for :math:`\gamma`.

More information on the computation of hydrodynamic forcing by the AeoLiS model is described in :ref:`this hydrodynamics section <water-levels-waves-run-up>`.


.. _bed-interaction-approach:

Bed-interaction Approach (:math:`\zeta`)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Most existing aeolian models assume that local bed properties dictate the saturated sediment concentration for the entire transport column, relying on a single value for :math:`c_{\mathrm{sat}}` in the advection equation (Equation :eq:`erodep`). However, this assumes that conditions at or near the bed instantly influence the entire air column, which would cause immediate deposition over obstacles like vegetation.

In reality, transport conditions depend on the vertical position of the sediment within the air column. For instance, over a non-erodible surface or a vegetated area, sediment near (or trapped in) the bed experiences severe drag reduction or a lack of supply, while sediment traveling higher in the air column could still "skim" over the top of the surface, retaining its momentum.

To capture this decoupling and allow sediment to bypass non-erodible or vegetated surfaces despite reduced transport conditions at the bed, a bed-interaction approach was recently introduced by :cite:`vanWesten2026`. This approach (enabled via ``process_bedinteraction``) divides the saturated sediment concentration into two distinct transport modes:

* **Bed-dominated transport** :math:`c_{\mathrm{sat,bed}}` (``CuBed``) [:math:`\mathrm{kg/m^2}`]: Sediment directly affected by conditions at or near the bed. This flux responds to local surface constraints, such as vegetation-induced shear reduction (:math:`u_{*,\mathrm{veg}}`) or supply-limiting factors (e.g., moisture, sheltering) that increase the velocity threshold (:math:`u_{\mathrm{th}}`).
* **Airborne transport** :math:`c_{\mathrm{sat,air}}` (``CuAir``) [:math:`\mathrm{kg/m^2}`]: Sediment that remains elevated above the canopy or internal boundary layer, responding primarily to the free-stream wind forcing (:math:`u_*`). This mode operates independently of local surface limitations.

The combined saturated sediment concentration :math:`c_{\mathrm{sat}}` (``Cu``) [:math:`\mathrm{kg/m^2}`] is calculated as a weighted sum of these two modes:

.. math::
   :label: csat_combined

   c_{\mathrm{sat}} = w_{\mathrm{air}} c_{\mathrm{sat,air}} + w_{\mathrm{bed}} c_{\mathrm{sat,bed}}

The dimensionless weights :math:`w_{\mathrm{air}}` and :math:`w_{\mathrm{bed}}` [:math:`\mathrm{-}`] are controlled by a bed-interaction factor :math:`\zeta` (``zeta``) [:math:`\mathrm{-}`], which defines the fraction of the sediment flux actively interacting with the bed. The airborne weight scales with the actual instantaneous sediment concentration :math:`c` (``Ct``) [:math:`\mathrm{kg/m^2}`] relative to the airborne carrying capacity:

.. math::
   :label: weights_zeta

   w_{\mathrm{air}} = (1 - \zeta) \frac{c}{c_{\mathrm{sat,air}}}, \quad \quad w_{\mathrm{bed}} = 1 - w_{\mathrm{air}}

.. note::
   In the case of fully saturated transport (:math:`c = c_{\mathrm{sat,air}}`), Equation :eq:`csat_combined` simplifies to :math:`c_{\mathrm{sat}} = (1 - \zeta) c_{\mathrm{sat,air}} + \zeta c_{\mathrm{sat,bed}}`. The relative saturation term :math:`\frac{c}{c_{\mathrm{sat,air}}}` basically prevents the model from bypassing sediment over an obstacle when the air is actually empty. If the incoming wind carries no sediment (:math:`c = 0`), the airborne weight drops to zero, and all transport must be initiated entirely by local bed conditions.

This formulation allows the model to dynamically simulate different surface types through the parameter :math:`\zeta`:

* Bare sediment surface (:math:`\zeta = 1`): All transport interacts with the bed. If local shear drops or the threshold increases, deposition occurs directly (scaled only by the adaptation time :math:`T` in Equation :eq:`erodep`). Under these conditions, Equation :eq:`csat_combined` simplifies to :math:`c_{\mathrm{sat}} = c_{\mathrm{sat,bed}}`. 
* Non-erodible surface (:math:`\zeta = 0`): Transport is fully decoupled from conditions at the bed. Even if the local pickup capacity is zero (:math:`c_{\mathrm{sat,bed}} = 0`), sediment entering from upwind passes over the cell without depositing, effectively creating an infinite deposition length.
* Vegetation (:math:`0 < \zeta < 1` **): The vegetation canopy intercepts a portion of the sediment flux, while the remainder skims over the top. Consequently, the effective deposition length scale increases, allowing the depositional area to extend further downwind into the vegetation patch.

The actual dynamic computation of :math:`\zeta` based on vegetation height and the vertical transport profile is described further down in the :ref:`vegetation section <vegetation>`.

.. note::
   The bed interaction parameter (``bi``) described earlier is currently distinct from the recently introduced bed-interaction factor (``zeta``). While conceptually similar, they serve different purposes:

   * ``bi``: A static variable that weights the contribution of the bed composition (sand in the bed) versus the airborne composition (sand already in transport) when computing the transport rate for multiple sediment fractions.
   * ``zeta``: A dynamically computed factor based on surface properties that decouples the air and bed. It determines the extent to which transport is governed by bed conditions (supply-limited) versus air conditions (wind-driven capacity).

   Future updates aim to consolidate these two parameters.

.. _wind-shear-velocity:

Wind and Shear Velocity
-----------------------------------

Wind is the primary driver of aeolian sediment transport (``process_wind``). Time-varying wind conditions are provided to the model through the ``wind_file``, which contains both the wind magnitude [:math:`\mathrm{m/s}`] and direction [:math:`\mathrm{^\circ}`]. This wind velocity is measured at a specific reference height :math:`z` (``z``) [:math:`\mathrm{m}`], which defaults to 10 m. The model automatically interpolates the input time series to the internal computational time steps and decomposes the wind into cross-shore :math:`u_{\mathrm{w,s}}` (``uws``) and longshore :math:`u_{\mathrm{w,n}}` (``uwn``) [:math:`\mathrm{m/s}`] components.

.. _shear-velocity:

Shear Velocity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The shear velocity over a flat bed :math:`u_{*,0}` (``ustar0``) [:math:`\mathrm{m/s}`] is first computed from the wind velocity :math:`u_{\mathrm{w}}` (``uw``) [:math:`\mathrm{m/s}`] using the Prandtl-Von Kármán Law of the Wall:

.. math::
   :label: law_of_the_wall

   u_{*,0} = \frac{u_{\mathrm{w}} \kappa}{\ln\left(\frac{z}{z_0}\right)}

where :math:`\kappa` (``kappa``) [:math:`\mathrm{-}`] is the Von Kármán constant, :math:`z` (``z``) [:math:`\mathrm{m}`] is the reference elevation of the wind measurements, and :math:`z_0` (``z0``) [:math:`\mathrm{m}`] is the aerodynamic roughness length.

The model offers multiple methods to compute the roughness length :math:`z_0`, selected via the ``method_roughness`` parameter:

* ``constant`` (default): Uses the user-defined roughness parameter :math:`k` (``k``) [:math:`\mathrm{m}`] directly as the roughness length (:math:`z_0 = k`). This is implemented to ensure backward compatibility and does not follow the standard Nikuradse definition.
* ``constant_nikuradse``: Follows the definition introduced by Nikuradse, scaling the user-defined bed roughness by 30 (:math:`z_0 = k / 30`).
* ``mean_grainsize_initial``: Computes a static roughness based on the initial mean grain size across the domain (:math:`z_0 = d_{\mathrm{mean}} / 30`). This is most applicable to flat beds with a uniform grain size distribution.
* ``mean_grainsize_adaptive``: Dynamically updates the roughness through time and space based on the evolving local mean grain size.
* ``median_grainsize_adaptive``: Uses the local median grain size :math:`d_{50}` (:math:`z_0 = 2d_{50} / 30`). This approach is based on Sherman and Greenwood (1982) and is appropriate for naturally occurring grain size distributions.
* ``vanrijn_strypsteen``: An advanced dynamic formulation based on van Rijn and Strypsteen (2019) and Strypsteen et al. (2021). It calculates the roughness dynamically using the local :math:`d_{50}` and :math:`d_{90}` to account for the additional roughness generated by the saltation layer and ripple formation phases.

.. tip::
   Although the bed roughness :math:`k` (``k``) is a physical parameter, it can be used in practice to calibrate transport rates. A pragmatic workflow is to select the ``bagnold`` transport method, establish flat-bed transport with representative conditions, and use ``k`` (with ``method_roughness = constant``) to calibrate the transport magnitudes before adding more complexity to the model. This approach could be effective if your primary goal is to simulate representative morphological evolution and you're less interested in the accuracy of the physical representation of the aeolian transport.


.. _topographic-steering:

Topographic Steering
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To simulate the topographic steering effects on dunes, Computational Fluid Dynamics (CFD) methods are often used. However, their high computational expense makes them unsuitable for long-term morphodynamic simulations. To reduce computational costs, the topographic steering of the wind due to smooth gradients is implemented following an analytical perturbation theory for turbulent boundary layer flow :cite:`weng1991air, kroy2002minimal` (``process_shear``). 

This approach builds upon the flat bed shear velocity :math:`u_{*0}` (``ustar0``) [:math:`\mathrm{m/s}`] established in the previous step. This velocity provides a baseline, unperturbed shear stress :math:`\vec{\tau}_{0}` (``tau0``) [:math:`\mathrm{N/m^2}`], where :math:`|\vec{\tau}_{0}| = \rho_{\mathrm{a}} u_{*0}^2`. The method then computes the topographically steered shear stress :math:`\vec{\tau}(x,y)` (``tau``) by applying a spatial perturbation :math:`\delta\vec{\tau}(x,y)` to this baseline:

.. math::
   :label: shear_perturbation_base

   \vec{\tau}(x,y) = \vec{\tau}_{0} + |\vec{\tau}_{0}|\delta\vec{\tau}(x,y)

where :math:`\delta\vec{\tau}(x,y)` (``dtau``) is the shear stress perturbation and :math:`\vec{\tau}_{0}` is the computed shear stress on a flat topography. 

For two-dimensional situations, the shear stress perturbation in the x- and y-directions (:math:`\delta\tau_{x}` and :math:`\delta\tau_{y}`) (``dtaux``, ``dtauy``) is computed in Fourier space according to the following equations:

.. math::
   :label: shear_pert_x

   \delta\tilde{\tau}_{x}(\vec{k})=\frac{2\tilde{z}_{b}(\vec{k})}{U^2(l)}
   \frac{k_{x}^2}{|\vec{k}|}\left\lbrace-1+\left(2\ln\frac{l}{z'_{0}}+\frac{|k|^2}{k_{x}^2}\right)\sigma\frac{K_{1}(2\sigma)}{K_{0}(2\sigma)}\right\rbrace

.. math::
   :label: shear_pert_y

   \delta\tilde{\tau}_{y}(\vec{k})=\frac{2\tilde{z}_{b}(\vec{k})}{U^2(l)}
   \frac{k_{x}k_{y}}{|\vec{k}|}2\sqrt{2}\sigma K_{1} (2\sqrt{2}\sigma)

.. math::
   :label: shear_sigma
   
   \sigma=\sqrt{iLk_{x}z'_{0}/l}

where :math:`\tilde{}` indicates the Fourier-transformed components of the parameters, :math:`k_x` and :math:`k_y` are the components of the wave vector :math:`\vec{k}` in Fourier space, and :math:`K_0` and :math:`K_1` are modified Bessel functions. As illustrated in Figure :numref:`fig-concept-topo-steering`, the depth of the inner layer of flow :math:`l` [m], the dimensionless vertical velocity profile :math:`U(l)` [-] at height :math:`l`, and the height of the middle layer of flow :math:`z_m` [m] are defined as:

.. math::
   :label: flow_layer_parameters

   l = \frac{2 \kappa^2 L}{\ln \left( \frac{l}{z'_{0}} \right)} \qquad \qquad U(l) \equiv \frac{\ln \left( \frac{l}{z'_{0}} \right)}{\ln \left( \frac{z_{m}}{z'_{0}} \right)} \qquad \qquad z_{m} = \sqrt{\frac{L^2}{\ln \left( \frac{z_{m}}{z'_{0}} \right)}}

where :math:`L` (``L``) [m] is the typical length scale of the hill.

.. tip:: 
   Tuning the length scale :math:`L` (``L``): The typical length scale of the hill (:math:`L`) influences the strength of the shear perturbation. A higher :math:`L` will result in a stronger shear stress perturbation (higher wind speed-up over the crest or reduction at the toe) and thus in a more outspoken morphodynamic shape. A practical approach is to vary this parameter between 10 and 1000 m depending on your desired morphological outcome.

For one-dimensional situations, a simplified solution of the shear perturbation approach is implemented. By ignoring some minor terms, it provides a less computationally expensive approach :cite:`kroy2002minimal`:

.. math::
   :label: shear_pert_1d

   \delta \tau =\alpha \int_{-\infty}^{\infty}d\xi\frac{\frac{\delta z_b}{\delta x}(x-\xi)}{\pi \xi}+\beta\frac{\delta z_b}{\delta x}(x)

where :math:`\alpha` [-] and :math:`\beta` [-] both depend on :math:`L/z_0`, but are user-defined fixed variables rather than computed in the model. :math:`\xi` [-] is the normalized cross-shore distance :math:`x/L`.

Finally, the computed shear stresses are converted back into shear velocities. The topographically steered shear velocity :math:`u_*` (``ustar``) [:math:`\mathrm{m/s}`], which now includes the computed perturbation, is derived from the total shear stress :math:`\vec{\tau}(x,y)`:

.. math::
   :label: ustar_from_tau

   u_* = \sqrt{\frac{|\vec{\tau}(x,y)|}{\rho_{\mathrm{a}}}}


.. _flow-separation:
Flow separation
^^^^^^^^^^^^^^^^^

The implementation of the shear perturbation theory by :cite:`weng1991air` is only valid in situations with relatively smooth surfaces. The occurrence of steep slopes limits the validity of the approach. To address this, a description of flow separation is used following the Coastal Dune Model (CDM) :cite:`sauermann2001continuum, kroy2002minimal, DuranMoore2013` (``process_separation``). 

A smooth envelope is created, which separates the main flow when a sharp edge is detected in the windward direction, i.e. when the bed slope is steeper than a certain user-defined angle (``mu_b``). This smooth envelope is called a separation bubble, :math:`z_{sep}` (``zsep``) [m] (Figure :numref:`fig-concept-topo-steering`). This separation bubble represents the surface that divides the region of flow reversal from the main flow stream along the smooth hill. Subsequently, in all cells for which the bed level is lower than the separation bubble (:math:`z_b < z_{sep}`), the shear velocity :math:`u_{*}` is set to 0 m/s. This assumes that eventual flow reversal velocities are not significant enough to initiate aeolian transport.

The separation bubble surface :math:`z_{sep}` is modelled by a third-order polynomial. The height of the brinkline, or the location where the separation bubble starts to detach from the bed, is defined by :math:`z_b(x_{\mathrm{brink}}) \equiv z_{\mathrm{brink}}`. Assuming a maximum slope :math:`c` (``c_b``) [deg] for the separation surface that determines the shape of the bubble, the reattachment length :math:`l_r` is obtained by:

.. math::
   :label: reattachment_length

   l_r \approx \frac{3 z_{\mathrm{brink}}}{2c}\left(1+\frac{z_{\mathrm{brink}}}{4c}+2\left(\frac{z_{\mathrm{brink}}}{4c}\right)^2\right)

The separation bubble profile :math:`z_{sep}` is then calculated as:

.. math::
   :label: separation_bubble_poly

   z_{sep}(x)=a_3(x-x_{\mathrm{brink}})^3+a_2(x-x_{\mathrm{brink}})^2+z_{\mathrm{brink}}'(x-x_{\mathrm{brink}})+z_{\mathrm{brink}}

where the polynomial coefficients are:

.. math::
   :label: poly_coeffs

   a_2=-\frac{3 z_{\mathrm{brink}} + 2 z_{\mathrm{brink}}' l_r}{l_r^2} \qquad \qquad a_3=\frac{2 z_{\mathrm{brink}} +  z_{\mathrm{brink}}' l_r}{l_r^3}

.. tip:: 
   Enabling the separation bubble (``process_separation = T``) is recommended for solely aerodynamically dominated landforms (e.g., barchan dunes). However, it can produce undesirable morphodynamics in complex or irregular topographies, such as densely vegetated environments, so use it judiciously.

.. _fig-concept-topo-steering:

.. figure:: /images/concept_topo_steering.jpg
   :alt: concept topographic steering
   :width: 800px
   :align: center

   Schematic overview of the shear perturbation and flow separation approach. Based on :cite:`weng1991air` and :cite:`kroy2002minimal`.

The computed shear stress as a result of the combined influence of the implemented shear stress perturbations and flow separation is shown in Figure :numref:`fig-compare-topo-steering`. These results show the decrease on the windward and lee sides of both Gaussian- and barchan-shaped landforms and an increase over the crest. Additionally, a shear velocity of zero is shown below the separation bubble. 

.. _fig-compare-topo-steering:

.. figure:: /images/compare_topo_steering.jpg
   :alt: topographic steering
   :width: 800px
   :align: center

   Spatial variation in shear stress due to topographic steering of the wind field. The upper panels show the bed level :math:`z_b` [m] and shear stress velocity perturbation :math:`\delta u_*` [m/s] over a uniform Gaussian hill. The lower panels show topographic steering over a barchan dune, including the influence of flow separation. The right panels compare outcomes of the one- and two-dimensional approaches.

Directional winds
~~~~~~~~~~~~~~~~~~~

The underlying implementation of the perturbation theory and separation bubble originally allows only for wind conditions that are perpendicular to the grid. An overlaying computational grid is introduced in AeoLiS, which rotates with the changing wind direction per time step. By doing this, the shear stresses are always estimated in the positive x-direction of the computational grid. The following steps are executed for each time step:

1. Create a computational ('Rotational') grid aligned with the wind direction (``set_computational_grid``).
2. Add and fill a buffer around the original ('Primary') grid.
3. Populate the computational grid by rotating it to the current wind direction and interpolate the original topography onto it. 
4. Compute the morphology-wind induced shear stress by using the perturbation theory.
5. Add the wind-induced shear stresses to the computational grid.
6. Rotate both the grids and the total shear stress results in the opposite direction.
7. Interpolate the total shear stress results from the computational grid to the original grid.
8. Rotate the wind shear stress results and the original grid back to the original orientation.

.. tip:: 
   Generating the rotational grid is computationally expensive. Consider the following optimizations to improve performance:

   * **Resolution:** The secondary rotational grid uses its own resolution parameters (``dx`` and ``dy``) defined in the configuration file. Ensure these are not coarser than your primary grid resolution to prevent information loss.
   * **Buffer zones (FFT):** The FFT method assumes periodic boundaries. If the bed elevations at opposite edges of the domain differ, the FFT will artificially include this steep gradient, causing numerical wiggles. To mitigate this, define a buffer zone (``buffer``) that smoothly connects the edges outside your domain of interest. As a rule of thumb, set the buffer width to at least 5 times the maximum height difference between the domain edges. Check the boundaries for wiggles after an initial run, and reduce the buffer size if possible to save computational time.
   * **Grid shape and 1D alternative:** Because the secondary grid acts as a bounding box, its size increases drastically for highly non-rectangular domains under oblique winds (e.g., a 1x100 grid under a 45° wind requires a secondary grid of approximately 71x71 cells). For (semi-)1D simulations, avoid this overhead by selecting ``method_shear = 1Dstacks``. This uses a direct analytical solution, bypassing the need for FFTs and rotational grids.


.. _vid-rotating-shear:

.. video:: /images/rotating_shear.mp4
   :autoplay:
   :loop:
   :muted:
   :width: 100%

   Animation demonstrating the rotational computational grid aligning with the changing wind direction to solve for topographic steering at each time step.


.. _shear-velocity-threshold:

Shear Velocity Threshold
---------------------------

The shear velocity threshold represents the influence of bed surface
properties in the saturated sediment transport equation (``process_threshold``). The shear
velocity threshold is computed for each grid cell and sediment
fraction separately based on local bed surface properties, like
moisture, roughness elements and salt content. For each bed surface
property supported by the model a factor is computed to increase the
initial shear velocity threshold:

.. math::
  :label: apx-shearvelocity
  
  u_{\mathrm{* th}} = 
  f_{u_{\mathrm{* th}}, \mathrm{M}} \cdot 
  f_{u_{\mathrm{* th}}, \mathrm{R}} \cdot 
  f_{u_{\mathrm{* th}}, \mathrm{NE}} \cdot 
  u_{\mathrm{* th, 0}}

.. _base-threshold-grainsize:

Grainsize (base)
^^^^^^^^^^^^^^^^^^^^^^^^^^

The base shear velocity threshold :math:`u_{\mathrm{* th, 0}}` (``uth0``) [m/s] is
computed based on the grain size following :cite:`Bagnold1937b`:

.. math::
   :label: shear

   u_{\mathrm{* th, 0}} = A \sqrt{ \frac{\rho_{\mathrm{p}} - \rho_{\mathrm{a}}}{\rho_{\mathrm{a}}} \cdot g \cdot d_{\mathrm{n}}}

where :math:`A` (``Aa``) [-] is an empirical constant, :math:`\rho_{\mathrm{p}}` (``rhog``)
[:math:`\mathrm{kg/m^3}`] is the grain density, :math:`\rho_{\mathrm{a}}` (``rhoa``)
[:math:`\mathrm{kg/m^3}`] is the air density, :math:`g` (``g``) [:math:`\mathrm{m/s^2}`] is the
gravitational constant and :math:`d_{\mathrm{n}}` (``grain_size``) [m] is the nominal grain
size of the sediment fraction.


.. _threshold-moisture-content:

Moisture content
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The shear velocity threshold (``th_moisture``) is updated based on moisture content
following :cite:`Belly1964`:

.. math::
  :label: apx-moist
   
  f_{u_{\mathrm{* th}}, \mathrm{M}} = \max(1 \quad ; \quad 1.8 + 0.6 \cdot \log(p_{\mathrm{g}}))

where :math:`f_{u_{\mathrm{* th},M}}` [-] is a factor in Equation :eq:`apx-shearvelocity`, :math:`p_{\mathrm{g}}` [-] is the geotechnical
mass content of water, which is the percentage of water compared to
the dry mass. The geotechnical mass content relates to the volumetric
water content :math:`p_{\mathrm{V}}` [-] according to:

.. math::
  :label: vol-water

  p_{\mathrm{g}} = \frac{p_{\mathrm{V}} \cdot \rho_{\mathrm{w}}}{\rho_{\mathrm{p}} \cdot (1 - p)}

where :math:`\rho_{\mathrm{w}}` [:math:`\mathrm{kg/m^3}`] and
:math:`\rho_{\mathrm{p}}` [:math:`\mathrm{kg/m^3}`] are the water and particle
density respectively and :math:`p` [-] is the porosity. Values for
:math:`p_{\mathrm{g}}` smaller than 0.005 do not affect the shear velocity
threshold :cite:`Pye1990`. Values larger than 0.064 (or 10\%
volumetric content) cease transport :cite:`DelgadoFernandez2010`,
which is implemented as an infinite shear velocity threshold.

.. _threshold-roughness-elements:

Roughness elements
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Sediment sorting may lead to the emergence of non-erodible elements
from the bed. Non-erodible roughness elements may shelter the erodible
bed from wind erosion due to shear partitioning, resulting in a
reduced sediment availability :cite:`Raupach1993` (``th_sheltering``). Therefore the
equation of :cite:`Raupach1993` is implemented according to:

.. math::
   :label: raupach
           
   u_{\mathrm{* th, R}} = u_{\mathrm{* th}} \cdot \sqrt{ \left( 1 - m \cdot \sum_{k=k_0}^{n_{\mathrm{k}}}{w_k^{\mathrm{bed}}} \right) \left( 1 + \frac{m \beta}{\sigma} \cdot \sum_{k=k_0}^{n_{\mathrm{k}}}{w_k^{\mathrm{bed}}} \right) }

in which :math:`\sigma` (``sigma``) is the ratio between the frontal area and the
basal area of the roughness elements and :math:`\beta` (``beta``) is the ratio
between the drag coefficients of the roughness elements and the bed
without roughness elements. :math:`m` (``m``) is a factor to account for the
difference between the mean and maximum shear stress and is usually
chosen 1.0 in wind tunnel experiments and may be lowered to 0.5 for
field applications. The roughness density :math:`\lambda` in the
original equation of :cite:`Raupach1993` is obtained from the mass
fraction in the bed surface layer :math:`w_k^{\mathrm{bed}}` according
to:

.. math::
   :label: rough
   
   \lambda = \frac{\sum_{k=k_0}^{n_{\mathrm{k}}}{w_k^{\mathrm{bed}}}}{\sigma}

in which :math:`k_0` is the index of the smallest non-erodible
sediment fraction in current conditions and :math:`n_{\mathrm{k}}` is the
total number of sediment fractions. It is assumed that the sediment
fractions are ordered by increasing size. Whether a fraction is
erodible depends on the sediment transport capacity.

.. _threshold-non-erodible-layer:

Non-erodible layer
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The model allows for the definition of a spatially varying non-erodible layer beneath the active sand surface (``ne_file``). When the bed elevation :math:`z_{\mathrm{b}}` (``zb``) [:math:`\mathrm{m}`] erodes down to the elevation of this non-erodible layer :math:`z_{\mathrm{ne}}` (``zne``) [:math:`\mathrm{m}`], the sediment supply is completely cut off. 

This process (``th_nelayer``) implements the restriction by forcing the non-erodible scaling factor :math:`f_{u_{\mathrm{* th}}, \mathrm{NE}}` to approach infinity:

.. math::
   :label: nelayer

   f_{u_{\mathrm{* th}}, \mathrm{NE}} = 
   \begin{cases} 
      \infty & \text{if } z_{\mathrm{b}} \leq z_{\mathrm{ne}} \\ 
      1 & \text{if } z_{\mathrm{b}} > z_{\mathrm{ne}} 
   \end{cases}

This effectively raises the velocity threshold to infinity, instantly ceasing any further entrainment from that grid cell. 

.. tip:: 
   AeoLiS does not include groundwater processes on the upper beach or in the dunes. If you suspect a wet layer restricts erosion on the upper beach or within dune slacks, the non-erodible layer can function as a proxy to keep the profile stable. It is also an essential processes for simulating bedforms migrating over hard surfaces, such as barchan dunes.




.. _vegetation:

Vegetation 
----------

The vegetation method is selected using ``method_vegetation``. The model currently supports two approaches:

1. **Original Method** (``duran``): A standard approach typical of established aeolian sediment transport models, relying on a fixed geometric relationship between plant height and cover.
2. **New Ecomorphodynamic Framework** (``grass``): A newly implemented framework designed for dune grasses (e.g., European and American marram grass). It decouples vertical growth from horizontal expansion and introduces concepts like vegetation bending, seed dispersal, competition, wake recovery, and two-layer sediment transport (:ref:`fig-vegetation-overview`).

.. warning::
   The ``grass`` framework is a new addition based on recent research (still in review). It is currently in a more experimental state compared to the more extensively applied ``duran`` method.

.. _fig-vegetation-overview:

.. figure:: /images/vegetation_overview.png
   :width: 800px
   :align: center

   Conceptual overview of the new vegetation framework illustrating its four core components: vegetation metrics, development, shear reduction, and the two-layer sediment transport approach.


.. _vegetation-metrics:

Vegetation Metrics
^^^^^^^^^^^^^^^^^^^

(**Method:** ``duran``) The basal vegetation density :math:`\rho_{\text{veg}}` (``rhoveg``) [:math:`\mathrm{-}`] can vary in space and time. It is determined by the ratio of the actual vegetation height :math:`h_{\text{veg}}` (``hveg``) [:math:`\mathrm{m}`] to the maximum attainable vegetation height :math:`H_{\text{veg}}` (``Hveg``) [:math:`\mathrm{m}`], varying between 0 and 1 :cite:`DuranHerrmann2006`:

.. math::
   :label: Vegetation_density_duran

   \rho_{\text{veg}} = \left( \frac{h_{\text{veg}}}{H_{\text{veg}}} \right)^2

This assumption is based on the idea that burying vegetation reduces its height, which indicates a simultaneous decrease in actual cover. The change in vegetation density per grid cell is directly linked to the alteration in vegetation height within that specific cell. 

(**Method:** ``grass``) The new framework decouples plant structure into two independent state variables: tiller density :math:`N_t` (``N_t``) [:math:`\mathrm{tillers/m^2}`] and tiller height :math:`h_{\text{veg}}` (``hveg``) [:math:`\mathrm{m}`]. This separation enables the explicit representation of distinct morphological states, such as sparse/tall canopies or dense/short canopies. 

To capture fine-scale spatial dynamics, such as clonal expansion, these vegetation metrics and their subsequent developmental processes are resolved on an automatically generated higher-resolution sub-grid (:math:`\Delta x \leq 1` m).

For dune grasses (e.g., European beachgrass), typical tiller densities range from 400 to 1110 tillers/m², heights reach 0.65 to 0.85 m, and tiller diameters :math:`d_t` are 0.004–0.008 m. From these dimensions, the frontal area index :math:`\lambda_{\text{veg}}` (``lambdaveg``) [:math:`\mathrm{m^2/m^2}`] and the basal cover fraction :math:`\rho_{\text{veg}}` (``rhoveg``) [:math:`\mathrm{m^2/m^2}`] are derived explicitly:

.. math::
   :label: Vegetation_density_grass

   \lambda_{\text{veg}} = N_t\,h'_{\text{veg}}\,d_t \quad , \quad \rho_{\text{veg}} = N_t\,\pi\,\left(\frac{d_t}{2}\right)^2

Because grass tillers bend under wind forcing, their effective height interacting with the wind decreases dynamically. The effective tiller height :math:`h'_{\text{veg}}` [:math:`\mathrm{m}`] is computed as a function of the unbent height :math:`h_{\text{veg}}`, wind speed :math:`u_w` [:math:`\mathrm{m/s}`], and tiller density :math:`N_t`:

.. math::
   :label: effective_height_bending

   h'_{\text{veg}} = h_{\text{veg}} \left[ r_{\text{stem}} + (1 - r_{\text{stem}}) \left( \alpha_u u_w + \alpha_N N_t + \alpha_0 \right) \right]

Here, :math:`r_{\text{stem}}` [:math:`\mathrm{-}`] specifies the fraction of the stem that remains rigid, while the empirical constants :math:`\alpha_u` [:math:`\mathrm{s/m}`], :math:`\alpha_N` [:math:`\mathrm{m^2}`], and :math:`\alpha_0` [:math:`\mathrm{-}`] control the sensitivity to wind stress and the structural support provided by neighboring tillers.

.. _vegetation-development:

Vegetation Development
^^^^^^^^^^^^^^^^^^^^^^^

**Method:** ``duran``

Vegetation growth and decay follow the model proposed by :cite:`DuranHerrmann2006`, modified to include an optimal burial rate :math:`\Delta z_{\text{b,opt}}` [:math:`\mathrm{m/yr}`] that shifts the peak of optimal growth:

.. math::
   :label: changes_vegetation_height_duran

   \frac{\partial h_{\text{veg}}}{\partial t} = V_{\text{ver}} \left(1 - \frac{h_{\text{veg}}}{H_{\text{veg}}}\right) - \gamma_{\text{veg}} \left| \Delta z_{\text{burial}} - \Delta z_{\text{b,opt}} \right|

Here, :math:`\gamma_{\text{veg}}` (``veg_gamma``, default = 1) [:math:`\mathrm{-}`] accounts for the impact of sediment burial. :math:`V_{\text{ver}}` (``V_ver``) is the maximum vertical growth rate [:math:`\mathrm{m/yr}`], while the sediment burial rate :math:`\Delta z_{\text{burial}}` [:math:`\mathrm{m/yr}`] is determined as the bed level change averaged over a trailing time window (default is one day) to prevent vegetation from overreacting to instantaneous bed level fluctuations. 

The optimal burial rate for maximum vegetation growth for marram grass is around 0.31 m/year with a burying tolerance of 0.78 to 0.96 m burial/year :cite:`Nolet2018`. Vegetation can begin to grow through lateral propagation or random germination handled on a cell-by-cell basis using a probabilistic approach similar to :cite:`Keijsers2016`.

.. tip:: 
   Meaning of these variables: An intrinsic vertical growth rate of :math:`V_{\text{ver}} = 4` m/year does not mean the vegetation will be 4 meters high after 1 year, as growth follows a logistic curve that slows as it reaches :math:`H_{\text{veg}}`.

**Method:** ``grass``

Vegetation development is simulated through two completely decoupled processes: vertical tiller growth and horizontal tiller establishment (:ref:`fig-vegetation-development`).

.. _fig-vegetation-development:

.. figure:: /images/vegetation_growth.mp4
   :width: 900px
   :align: center

   Simulated spatial and temporal evolution of tiller density and height demonstrating local growth, clonal expansion, seedling dispersal, and inter-species competition.


**1. Vertical Tiller Growth:**
Vertical growth utilizes a generalized logistic growth equation, driven by the intrinsic growth rate :math:`G_h` [:math:`\mathrm{m/yr}`] and an exponent :math:`\phi_h` [:math:`\mathrm{-}`] that provides greater control over the growth trajectory:

.. math::
   :label: growth_height_grass

   \frac{\partial h_{\text{veg}}}{\partial t} = G_h \left(1 - \frac{h_{\text{veg}}}{H_{\text{veg}}}\right)^{\phi_h} + B_h

The response to sediment burial and erosion :math:`B_h` relies on the sensitivity parameter :math:`\gamma_h` [:math:`\mathrm{-}`] and an optimal burial rate :math:`\Delta z_{\text{opt},h}` [:math:`\mathrm{m/yr}`] (e.g., 0.2–1.0 m/yr for marram grass), producing an asymmetric response to burial versus erosion.

**2. Horizontal Tiller Establishment (Density):**
Tiller density evolves through the recruitment of new tillers via seedling germination (:math:`s`) and clonal expansion (:math:`c`). Tiller production occurs in a source cell (:math:`j`) and is distributed to a target cell (:math:`i`) based on dispersal weights :math:`w_{ij}`:

.. math::
   :label: density_growth_grass

   \frac{\partial N_{t,i}}{\partial t} = \sum_j w^{(s)}_{ij} S_{s,j} \;+\; \sum_j w^{(c)}_{ij} S_{c,j} \left(1 - \sum_n \alpha_{AB} \frac{\bar{N}_{t,i}^{(B)}}{N_{t,\text{max}}^{(B)}}\right)

To account for inter-specific competition in multi-species simulations, the logistic saturation term utilizes a Lotka-Volterra approach, where :math:`\alpha_{AB}` [:math:`\mathrm{-}`] defines the competitive effect of species :math:`A` on target species :math:`B`. Production in the source cell is calculated via:

.. math::
   :label: tiller_production

   S_{x,j} = G_x N_{t,j} B_{x,i} \left( \frac{h_{\text{veg},j}}{H_{\text{veg}}} \right)

The spatial dispersal mechanisms differ fundamentally:

* **Clonal Expansion:** Modeled as a short-range, Lévy-like spreading strategy using a truncated Pareto distribution governed by a shape parameter :math:`\mu_c` [:math:`\mathrm{-}`]. Dispersal weights :math:`w^{(c)}` are discretized into a spatial kernel and sampled via a Poisson distribution.
* **Seedling Dispersal:** Stochastically sampled using a two-dimensional Student's *t*-distribution (2Dt), capturing heavy-tailed long-distance transport via a scale parameter :math:`a_s` [:math:`\mathrm{m^2}`] and shape parameter :math:`\nu_s` [:math:`\mathrm{-}`] to define dispersal weights :math:`w^{(s)}`.

.. _vegetation-induced-shear-reduction:

Vegetation-induced Shear Reduction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Method:** ``duran``

Inspired by the Coastal Dune Model (CDM), AeoLiS incorporates vegetation-wind interaction using the simplified expression:

.. math::
   :label: shear_reduction_vegetation_duran

   \frac{u_{*,\text{veg}}}{u_*} = \frac{1}{\sqrt{1 + \Gamma \rho_{\text{veg}}}}

The ratio of shear velocity in the presence of vegetation (:math:`u_{*,\text{veg}}`) [:math:`\mathrm{m/s}`] to the unobstructed shear velocity (:math:`u_*`) [:math:`\mathrm{m/s}`] is driven by the basal vegetation cover :math:`\rho_{\text{veg}}` and a fixed vegetation-related roughness parameter :math:`\Gamma` (``gamma_vegshear``, default = 16) [:math:`\mathrm{-}`].

**Method:** ``grass``

Rather than relying on the basal cover assumption, the updated framework calculates local shear velocity reduction by explicitly using the frontal area index :math:`\lambda_{\text{veg}}`:

.. math::
   :label: shear_reduction_vegetation_grass

   R_{\mathrm{0,veg}} = \frac{u_{*,\text{veg}}}{u_*} = \frac{1}{\sqrt{1 + m \beta_{\text{veg}} \lambda_{\text{veg}}}}

Here, :math:`\beta_{\text{veg}}` (``beta_veg``) [:math:`\mathrm{-}`] represents the drag efficiency of the vegetation elements relative to the bare surface, and :math:`m` [:math:`\mathrm{-}`] accounts for spatial non-uniformity in the surface shear stress distribution. 

Beyond local drag reduction, the framework captures non-local shear stress recovery in the sheltered wake downwind of the plant (:ref:`fig-vegetation-shear-params`). Using a decay function proposed by :cite:`Okin2008`, the spatial reduction factor :math:`R_{\text{veg}}(x)` [:math:`\mathrm{-}`] gradually recovers to free-stream conditions:

.. math::
   :label: wake_recovery

   R_{\text{veg}}(x) = 1 - (1 - R_{\mathrm{0,veg}}) e^{-x c_1 / h'_{\text{veg}}}

Here, :math:`c_1` [:math:`\mathrm{-}`] is a dimensionless calibration constant controlling the wake length.

.. _fig-vegetation-shear-params:

.. figure:: /images/rveg_shear_reduction.png
   :width: 900px
   :align: center

   Influence of varying vegetation metrics and calibration parameters on the spatial distribution of shear reduction and corresponding bed level changes.


.. _vegetation-computing-zeta:

Computing bed-interaction (zeta) over vegetation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Method:** ``duran``

In the standard advection scheme, the model implicitly assumes that local bed properties dictate the saturation concentration for the entire transport column. Therefore, the bed-interaction factor :math:`\zeta` [:math:`\mathrm{-}`] is essentially assumed to be 1, meaning any reduction in shear stress due to vegetation immediately forces the entire sediment flux to deposit.

**Method:** ``grass``

To capture realistic "skimming" flows over dense grass canopies, the new framework divides the saturation concentration :math:`c_{\text{sat}}` into two distinct modes (:ref:`fig-vegetation-sediment-transport`):

1. **Bed-affected transport (**:math:`c_{\text{sat,bed}}`**):** Sediment directly interacting with the canopy and restricted by local drag reduction.
2. **Airborne transport (**:math:`c_{\text{sat,air}}`**):** Sediment elevated above the canopy, responding to the free-stream wind and bypassing the vegetation.

.. _fig-vegetation-sediment-transport:

.. figure:: /images/vegetation_sediment_aeolis.png
   :width: 900px
   :align: center

   Vertical transport distribution over bare sand, non-erodible layers, and varying vegetation canopies, illustrating the computation of the bed-interaction parameter :math:`\zeta`.


The combined saturation is computed via weighted sum: :math:`c_{\text{sat}} = w_{\text{air}} c_{\text{sat,air}} + w_{\text{bed}} c_{\text{sat,bed}}`. These dimensionless weights are controlled by the bed-interaction factor :math:`\zeta` (``zeta``) [:math:`\mathrm{-}`], which dictates the fraction of the flux actively interacting with the bed. 

To estimate :math:`\zeta` over vegetation, the model computes an uplifted vertical transport profile using a Weibull distribution:

.. math::
   :label: weibull_transport

   f(h) = \frac{k}{h_{\text{scale}}} \left(\frac{h}{h_{\text{scale}}}\right)^{k-1} \exp\left[-\left(\frac{h}{h_{\text{scale}}}\right)^k\right]

The extent of the vegetation-induced lift is determined by a lifting coefficient :math:`\alpha_{\text{lift}}` [:math:`\mathrm{-}`] that scales the effective tiller height to define the physical lift height :math:`h_{\text{lift}} = L_h + \alpha_{\text{lift}} h'_{\text{veg}}`. By integrating this profile up to the effective tiller height :math:`h'_{\text{veg}}`, the model calculates the raw trapped fraction :math:`\zeta_0`. Because sparse vegetation does not trigger full skimming, this is adjusted by relative tiller density (using parameter :math:`\theta_\zeta`):

.. math::
   :label: zeta_density

   \zeta_0 = 1 - \left(\frac{N_t}{N_{t,\text{max}}}\right)^{\theta_\zeta} (1 - \zeta_0)

The final bed-interaction factor accounts for airborne sediment bouncing through the canopy using a numerical bounce factor :math:`b` [:math:`\mathrm{-}`]: :math:`\zeta = \zeta_0 (1 - b)`. 

.. _vegetation-mortality:

Vegetation Mortality
^^^^^^^^^^^^^^^^^^^^^^^^

**Method:** ``duran``

Vegetation is subject to destruction caused by hydrodynamic processes. In the event of cell inundation by high water levels, the vegetation density :math:`\rho_{\text{veg}}` in the affected grid cells is instantaneously or proportionally reduced to mimic storm-induced erosion of the canopy.

**Method:** ``grass``

Because tiller height and tiller density are fundamentally coupled, changes in one necessitate updates in the other. Mortality and structural changes occur through three primary drivers:

**1. Burial and Erosion Dieback:** 
When local burial/erosion rates exceed the plant's tolerance range, vertical growth becomes negative (:math:`\partial h_{\text{veg}}/\partial t < 0`). To represent the physical thinning of the dying patch, the tiller density decays proportionally to the relative loss in height, using a thinning factor :math:`\gamma_{N_t}` [:math:`\mathrm{-}`]:

.. math::
   :label: density_dieback

   \frac{\partial N_t}{\partial t} = \gamma_{N_t} \frac{N_t}{h_{\text{veg}}} \frac{\partial h_{\text{veg}}}{\partial t} 

**2. Recruitment Averaging:**
When new tillers establish via seed or clonally (:math:`\partial N_t/\partial t > 0`), they emerge at the bed surface with an initial height of zero. Therefore, the cell-averaged vegetation height mathematically decreases proportionally:

.. math::
   :label: height_correction

   h_{\text{veg}}^{n+1} = h_{\text{veg}}^{n} \left( \frac{N_t}{N_t + \Delta N_t} \right)

**3. Inundation-driven Decay:** 
Tiller density :math:`N_t` decays continuously during inundation events at a rate proportional to the relative water depth :math:`h_w` [:math:`\mathrm{m}`] (computed as Total Water Level minus bed elevation):

.. math:: 
   :label: flood_mortality_grass

   \frac{\partial N_t}{\partial t} = -N_t \frac{h_w}{h_{\text{veg}} T_{\text{flood}}} 

Here, :math:`T_{\text{flood}}` (``T_flood``) [:math:`\mathrm{s}`] is the prescribed inundation decay timescale. While this random thinning removes tillers, it does not alter the average height of the surviving canopy.


.. _hydrodynamics-moisture:
Hydrodynamics and Surface Moisture
-------------------------------------

Aeolis computes .... for ... and ... and ...


.. _water-levels-waves-run-up:

Water levels, Waves and Run-up
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

PLACEHOLDER: TWL, SWL, zs, hw, eta, R, ...

The runup height and wave setup are computed using the Stockdon formulas :cite:`Stockdon2006`. 
Their parameterization differs depending on the dynamic beach steepness expressed through the Irribaren number:

.. math::
   :label: irribaren
   
   \xi  = \tan \beta /\sqrt {{H_0}/{L_0}}

where :math:`{H_0}` is the significant offshore wave height, :math:`{L_0}` is the deepwater wavelength, and :math:`{\tan \beta}` is the foreshore slope.

For dissipative conditions, :math:`{\xi}` < 0.3, the runup, :math:`{R_2}`, is parameterized as,

.. math::
   :label: runup_dissipative
   
   {R_2} = 0.043\sqrt {{H_0}{L_0}}
   
and wave setup:

.. math::
   :label: setup_dissipative
   
   < \eta  >  = 0.02\sqrt {{H_0}{L_0}}

For :math:`{\xi}` > 0.3, runup is paramterized as,

.. math::
   :label: runup
   
   {R_2} = 1.1\left( {0.35\beta \sqrt {{H_0}{L_0}}  + \frac{{\sqrt {{H_0}{L_0}\left( {0.563{\beta ^2} + 0.004} \right)} }}{2}} \right)

and wave setup:

.. math::
   :label: setup
   
   < \eta  >  = 0.35\xi


.. _surface-moisture:
Surface Moisture
^^^^^^^^^^^^^^^^^^^^^

PLACEHOLDER


.. _hydraulic-sediment-mixing:
Hydraulic Sediment Mixing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

PLACEHOLDER

.. _marine-driven-bed-level-change:
Marine-driven Bed Level Change
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Bed-level change due to marine-driven dune erosion occurs when the total water level exceeds the base of the dune, or the dune toe elevation (``dune_toe_elevation``). See :ref:`dune-erosion` for more information.

.. _groundwater-module:
Groundwater Module (Hallin, 2023)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


Wave runup, capillary rise from the beach groundwater, and precipitation periodically wet the intertidal beach
temporally increasing the shear velocity threshold (
:numref:`fig-moisture-processes`). Infiltration and
evaporation subsequently dry the beach.

.. _fig-moisture-processes:

.. figure:: /images/moisture_processes.jpg
   :align: center

   Illustration of processes influencing the volumetric moisture content :math:`\theta` at the beach surface.

The structure of the surface moisture module and included processes are schematized in :numref:`fig-moisture-scheme`. 
The resulting surface moisture is obtained by selecting the largest of the moisture contents computed 
with the water balance approach (right column) and due to capillary rise from the groundwater table (left column). 
The method is based on the assumption that the flow of soil water is small compared to the flow of groundwater 
and that the beach groundwater dynamics primarily is controlled by the water level and wave action at 
the seaward boundary :cite:`Raubenheimer1999`, :cite:`Schmutz2014`. Thus, there is no feedback between the processes 
in the right column of :numref:`fig-moisture-scheme` and the groundwater dynamics described in the left column.

.. _fig-moisture-scheme:

.. figure:: /images/moisture_scheme.jpg
   :width: 600px
   :align: center

   Implementation of surface moisture processes in the AeoLiS.

Groundwater under sandy beaches can be considered as shallow aquifers, with only horizontal groundwater
flow so that the pressure distribution is hydrostatic :cite:`Baird1998,Brakenhoff2019,Nielsen1990,Raubenheimer1999`.
The cross-shore flow dominates temporal variations of groundwater levels. Alongshore, groundwater table variations are typically small :cite:`Schmutz2014`.
Although the surface moisture model can be extended over a two-dimensional grid, the groundwater simulations are performed for 1D transects cross-shore
to avoid numerical instabilities at the seaward boundary and reduce computational time.

The beach aquifers is schematised as a sandy body, with saturated hydraulic conductivity, :math:`K`, and effective porosity, :math:`{n_e}`.
The aquifer is assumed to rest on an impermeable surface, where :math:`D` is the aquifer depth. 
The groundwater elevation relative to the mean sea level (MSL) is denoted :math:`\eta`, and the shore-perpendicular x-axis is positive landwards,
with an arbitrary starting point. The sand is assumed to be homogenous and isotropic. In this context, isotropy implies that hydraulic conductivity
is independent of flow direction.

The horizontal groundwater discharge per unit area, :math:`u`, is then governed by Darcy’s law,

.. math::
   :label: gw-discharge
   
   u =  - K\frac{{\partial \eta }}{{\partial x}}

and the continuity equation (see e.g., :cite:`Nielsen2009`), 

.. math::
   :label: gw-continuity

   \frac{{\partial \eta }}{{\partial t}} =  - \frac{1}{{{n_e}}}\frac{\partial }{{\partial x}}((D + \eta )u)

where :math:`t` is time. 

The groundwater overheight due to runup, :math:`{U_l}`, is computed by :cite:`Kang1994,Nielsen1988`,

.. math::
   :label: gw-runup

   {U_l} = \left\{ \begin{gathered}{C_l}Kf(x)\,\,\,\,{\text{if }}{x_S} \leqslant x \leqslant {x_R} \hfill \\0,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,{\text{if }}x > {x_R} \hfill \\\end{gathered}  \right.

where :math:`{C_l}` is an infiltration coefficient (-), and :math:`f(x)` is a function of :math:`x` ranging from 0 to 1. :math:`{x_S}` is 
the horizontal location of the sum of the still water level and wave setup, and :math:`{x_R}` is the horizontal location of the runup limit:

.. math::
   :label: gw-runup-distribution

   f(x) = \left\{ \begin{gathered}
   \frac{{x - {x_s}}}{{\frac{2}{3}\left( {{x_{ru}} - {x_s}} \right)}}\,\,\,\,\,\,\,\,\,\,\,\,\,if\,{x_s} < x \leqslant {x_s} + \frac{2}{3}\left( {{x_{ru}} - {x_s}} \right)\, \hfill \\
   3 - \frac{{x - {x_s}}}{{\frac{1}{3}\left( {{x_{ru}} - {x_s}} \right)}}\,\,\,\,\,if\,{x_s} + \frac{2}{3}\left( {{x_{ru}} - {x_s}} \right)\, < x < {x_{ru}} \hfill \\ 
   \end{gathered}  \right.

Substitution of :math:`u` (Equation :eq:`gw-discharge`) in the continuity equation (Equation :eq:`gw-continuity`) with the addition of :math:`{U_l}/{n_e}` gives the nonlinear Boussinesq equation:

.. math::
   :label: boussinesq

   \frac{{\partial \eta }}{{\partial t}} = \frac{K}{{{n_e}}}\frac{\partial }{{\partial x}}\left( {(D + \eta )\frac{{\partial \eta }}{{\partial x}}} \right) + \frac{{{U_l}}}{{{n_e}}}

Numerical solution of the Boussinesq groundwater equation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Boussinesq equation is solved numerically with a central finite difference 
method in space and a fourth-order Runge-Kutta integration technique in time:

.. math::
  :label: solve-boussinesq

       f(\eta ) = \frac{K}{{{n_e}}}\left[ {D\underbrace {\frac{{{\partial ^2}\eta }}{{\partial {x^2}}}}_a + \underbrace {\frac{\partial }{{\partial x}}\underbrace {\left\{ {\eta \frac{{\partial \eta }}{{\partial x}}} \right\}}_b}_c} \right]

The Runge-Kutta time-stepping, where :math:`\Delta t` is the length of the timestep, is defined as,

.. math::
  :label: runge-kutta
  
  \begin{gathered}
  \eta _i^{t + 1} = \eta _i^t + \frac{{\Delta t}}{6}\left( {{f_1} + 2{f_2} + 2{f_3} + {f_4}} \right) \hfill \\
  {f_1} = f(\eta _i^t) \hfill \\
  {f_2} = f\left( {\eta _i^t + \frac{{\Delta t}}{2}{f_1}} \right) \hfill \\
  {f_3} = f\left( {\eta _i^t + \frac{{\Delta t}}{2}{f_2}} \right) \hfill \\
  {f_4} = f\left( {\eta _i^t + \Delta t{f_3}} \right) \hfill \\ 
  \end{gathered}

where, :math:`i` is the grid cell in x-direction and :math:`t` is the timestep. The central difference solution to :math:`f(\eta)` is obtained through discretisation of the Boussinesq equation,

.. math::
  :label: a-solve
  
   {a_i} = \frac{{\eta _{i + 1}^{} - 2\eta _i^{} + \eta _{i - 1}^{}}}{{{{(\Delta x)}^2}}}

.. math::
      {b_i} = \frac{{\eta _i^{}\left( {\eta _{i + 1}^{} - \eta _{i - 1}^{}} \right)}}{{\Delta x}}

.. math::
      {c_i} = \frac{{\left( {b_{i + 1}^{} - b_{i - 1}^{}} \right)}}{{\Delta x}}

The seaward boundary condition is defined as the still water level plus the wave setup . 
If the groundwater elevation is larger than the bed elevation, there is a seepage face, 
and the groundwater elevation is set equal to the bed elevation. On the landward boundary, 
a no-flow condition, :math:`\frac{{\partial \eta }}{{\partial t}} = 0` (Neumann condition), or constant head, :math:`\eta = constant` (Dirichlet condition), is prescribed.

Capillary rise
~~~~~~~~~~~~~~~~
Soil water retention (SWR) functions describe the surface moisture due to capillary transport 
of water from the groundwater table :cite:`VanGenuchten1980`:

.. math::
   :label: vangenuchten

   \theta (h) = {\theta _r} + \frac{{{\theta _s} - {\theta _r}}}{{{{\left[ {1 + {{\left| {\alpha h} \right|}^n}} \right]}^m}}}


where :math:`h` is the groundwater table depth, :math:`\alpha` and :math:`n` are fitting parameters 
related to the air entry suction and the pore size distribution. The parameter :math:`m` is commonly 
parameterised as :math:`m = 1 - 1/n`.  

The resulting surface moisture is computed for both drying and 
wetting conditions, i.e., including the 
effect of hysteresis. The moisture contents computed with drying and wetting SWR functions are denoted :math:`{\theta ^d}(h)` and :math:`{\theta ^w}(h)`, respectively. 
When moving between wetting and drying conditions, the soil moisture content follows an intermediate 
retention curve called a scanning curve. The drying scanning curves are scaled from the main 
drying curve and wetting scanning curves from the main wetting curve. The drying scanning curve is then obtained from :cite:`Mualem1974`:

.. math::
   :label: mualem-drying

   {\theta ^d}({h_\Delta },h) = {\theta ^w}(h) + \frac{{\left[ {{\theta ^w}({h_\Delta }) - {\theta ^w}(h)} \right]}}{{\left[ {{\theta _s} - {\theta ^w}(h)} \right]}}\left[ {{\theta ^d}(h) - {\theta ^w}(h)} \right]

where :math:`{h_\Delta}` is the groundwater table depth at the reversal on the wetting curve. 

The wetting scanning curve is obtained from :cite:`Mualem1974`:

.. math::
   :label: mualem-wetting
   
   {\theta ^w}({h_\Delta },h) = {\theta ^w}(h) + \frac{{\left[ {{\theta _s} - {\theta ^w}(h)} \right]}}{{\left[ {{\theta _s} - {\theta ^w}({h_\Delta })} \right]}}\left[ {{\theta ^d}({h_\Delta }) - {\theta ^w}({h_\Delta })} \right]

where :math:`{h_\Delta}` is the groundwater table depth at the reversal on the drying curve.

Infiltration
~~~~~~~~~~~~~~
Infiltration is accounted for by assuming that excess water infiltrates until the moisture content reaches 
field capacity, :math:`{\theta_fc}`. The moisture content at field capacity is the maximum amount of water 
that the unsaturated zone of soil can hold against the pull of gravity. For sandy soils, 
the matric potential at this soil moisture condition is around - 1/10 bar. In equilibrium, 
this potential would be exerted on the soil capillaries at the soil surface when the water 
table is about 100 cm below the soil surface, :math:`{\theta _{fc}} = {\theta ^d}(100)`.

Infiltration is represented by an
exponential decay function that is governed by a drying time scale
:math:`T_{\mathrm{dry}}`. Exploratory model runs of the unsaturated soil with the HYDRUS1D
:cite:`Simunek1998` hydrology model show that the increase of the
volumetric water content to saturation is almost instantaneous with
rising tide. The drying of the beach surface through infiltration
shows an exponential decay. In order to capture this behavior the
volumetric water content is implemented according to:

.. math::
   :label: infiltration

   \frac{{d\theta }}{{dt}} = \left( {\theta  - {\theta _{fc}}} \right)\left( {{e^{ - \ln (2)\frac{{dt}}{{{T_{dry}}}}}}} \right)

An alternative formulation is used for simulations that does not account for ground water and SWR processes,

.. math::
  :label: apx-drying
   
  p_{\mathrm{V}}^{n+1} = \left\{
    \begin{array}{ll}
      p & \mathrm{if} ~ \eta > z_{\mathrm{b}} \\
      p_{\mathrm{V}}^n \cdot e^{\frac{\log \left( 0.5 \right)}{T_{\mathrm{dry}}} \cdot \Delta t^n} - E_{\mathrm{v}} \cdot \frac{\Delta t^n}{\Delta z} & \mathrm{if} ~ \eta \leq z_{\mathrm{b}} \\
    \end{array}
  \right.

where :math:`\eta` [m+MSL] is the instantaneous water level,
:math:`z_{\mathrm{b}}` [m+MSL] is the local bed elevation,
:math:`p_{\mathrm{V}}^n` [-] is the volumetric water content in time step
:math:`n`, :math:`\Delta t^n` [s] is the model time step and :math:`\Delta z` is the bed
composition layer thickness. :math:`T_{\mathrm{dry}}` [s] is the beach
drying time scale, defined as the time in which the beach moisture
content halves.

Precipitation and evaporation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A water balance approach accounts for the effect of precipitation and evaporation,

.. math::
   :label: precipitation

   \frac{{d\theta }}{{dt}} = \frac{{\left( {P - E} \right)\,}}{{\Delta z}}\,

where :math:`P` is the precipitation, :math:`E` is the evaporation, and :math:`\Delta z` is the thickness of the surface layer.

Evaporation is simulated using an adapted version
of the Penman-Monteith equation :cite:`Shuttleworth1993` that is
governed by meteorological time series of solar radiation, temperature
and humidity.

:math:`E_{\mathrm{v}}` [m/s] is the evaporation rate that is
implemented through an adapted version of the Penman equation
:cite:`Shuttleworth1993`:

.. math::
  :label: apx-penman
   
  E_{\mathrm{v}} = \frac{m_{\mathrm{v}} \cdot R_{\mathrm{n}} + 6.43 \cdot \gamma_{\mathrm{v}} \cdot (1 + 0.536 \cdot u_2) \cdot \delta e}
  {\lambda_{\mathrm{v}} \cdot (m_{\mathrm{v}} + \gamma_{\mathrm{v}})} \cdot 9 \cdot 10^7

where :math:`m_{\mathrm{v}}` [kPa/K] is the slope of the
saturation vapor pressure curve, :math:`R_{\mathrm{n}}`
[:math:`\mathrm{MJ/m^2/day}`] is the net radiance, :math:`\gamma_{\mathrm{v}}`
[kPa/K] is the psychrometric constant, :math:`u_2` [m/s] is the wind speed
at 2 m above the bed, :math:`\delta e` [kPa] is the vapor pressure deficit
(related to the relative humidity) and :math:`\lambda_{\mathrm{v}}` [MJ/kg]
is the latent heat vaporization. To obtain an evaporation rate in
[m/s], the original formulation is multiplied by :math:`9 \cdot 10^7`.



.. _morphological-change:

Morphological change
---------------------

PLACEHOLDER

.. _bed-level-update:

Bed Level Update
^^^^^^^^^^^^^^^^^

PLACEHOLDER

.. _avalanching:
Avalanching
^^^^^^^^^^^^^

PLACEHOLDER

.. _dune-erosion

Dune Erosion
^^^^^^^^^^^^

Wave-driven dune erosion occurs when the TWL exceeds the dune toe elevation (``dune_toe_elevation``) [m]. The amount of sediment eroded from the dune is dependent on the frequency of collisions with the dune and the exceedenace of the TWL over the dune toe elevation. The volume of eroded sediment is calculated following the Palmsten and Holman (2012) dune erosion formula: 

.. math::

   V = 4 C_s (TWL - z_{\mathrm{toe}})^2 N_c

where :math:`V` \left[ \frac{m^3}{m} \right] is volume eroded, :math:`C_s` is the dune erodibility coefficient, and :math:`N_c` is the number of bore collisions. The volume of sediment is removed landward of the dune toe elevation contour and the avalanching process prevents formation of vertial scarps. 

Beach Evolution
^^^^^^^^^^^^^^^

Beach shape and size contributes to the overall sediment supply available for aeolian sediment transport. AeoLiS includes numerous approaches to represent temporal and spatial variability in sediment supply related to beach evolution (Figure 1a-d). These approaches do not explicitly currently account for wave-driven processes and their role on beach shape and volume changes, however the available methods are meant to mimic realistic expected behaviors and avoid the need to couple model interfaces with external tools. For the purposes of this beach sediment supply function, the shoreline is defined as the seaward boundary of the beach profile (default :math:'xshoreline' and :math:'zshoreline' are 0 m; e.g., Figure 1e-f)) and shoreline change rate (:math:'shoreline_change_rate') is the rate of change at the 0 m contour. Four specific methods, specified in the input file as ``method_wet_supply``, are implemented, as follows below:

``wet_bed_reset``maintains stability of the bed by assuming any beach volume loss from aeolian sediment transport below the maximum wave runup level is replenished by marine processes. This assumption means the inundated beach profile is continuously reset to its initial morphology (Figure 1a). Input dune toe elevation and beach slope are not used in this method. 

.. note::
   ``wet_bed_reset`` is the default ``method_wet_supply``, while the other methods may only be necessary if modeling a long-term dune evolution case study with    high shoreline change rates.

``vertical_beach_growth`` converts a user input horizontal ``shoreline_change_rate`` (default 0 m contour) into a vertical beach accretion rate with the following equation: 

.. math::
   v_{\mathrm{rate}} = \mathrm{SCR} \cos\left(\frac{\pi}{2} - \tan^{-1}(\mathrm{beach\slope})\right)

The ``vertical_beach_growth`` method models a linear beach with increased elevation every time step (Figure 1c). The beach maintains a fixed input ``beach_slope`` throughout the simulation and has an upper bound of the input ``dune_toe_elevation`` (Figure 1g-h). Though this method simulates sediment supply to the dune, it is important to note that over longer simulation time, the beach width is not maintained. 

``constant_SCR_constant_tanB`` simulates horizontal and vertical shoreline change using a fixed input ``beach_slope`` and fixed input ``dune_toe_elevation``. Like method vertical_beach_growth, ``constant_SCR_constant_tanB`` generates a new linear beach profile at each timestep, based on the input ``beach_slope`` and upper bounded by the input ``dune_toe_elevation``. However, this method also extends the beach seaward to maintain the beach width throughout the simulation (Figure 1b), resulting in the evolution of both shoreline elevation and position (Figure 1e-f). To use ``constant_SCR_constant_tanB``, it is important to note that the initial input grids must have extra seaward x-domain added to allow the beach to prograde.  

``constant_SCR_variable_tanB`` simulates horizontal shoreline change while allowing the ``beach_slope`` to evolve with the shoreline position and dune toe elevation to vary slightly with the upper bound of the input ``dune_toe_elevation`` (Figure 1d-h). A new beach slope is calculated and used to generate a linear beach each timestep. Similar to the ``constant_SCR_constant_tanB``, the x-domain of the input files must be extended seaward to implement this method. This method simulates more natural beach and dune evolution than the other methods available. 

.. _fig-method_wet_supply:

.. figure:: /images/aeolis_wet_supply.jpg
   :alt: wet_supply

   
