.. _model_description:

Model description
=================

Quick Overview
--------------

This section provides a summary of the main processes, equations, and configuration parameters in AeoLiS. For more information, refer to the detailed sections linked in the text (or scroll further down this page).

For guidance on setting up an AeoLiS model, see the :ref:`model input and output guide <model-input-output>`. The main configuration file (default: ``aeolis.txt``) is the basis of the model setup and contains all parameter settings and process-flags, and serves as the central reference for other input files. The computational domain is constructed using x- and y-coordinates (``xgrid_file``, ``ygrid_file``) alongside the initial bed elevation (``bed_file``). External environmental forcing is defined through continuous time series of wind (``wind_file``), water levels (``tide_file``), and waves (``wave_file``).

The simulation advances sequentially through time steps, repeating all activated processes and continuously updating the morphological model state. The simulation duration runs from a defined start time (``tstart``) to an end time (``tstop``), both specified in seconds relative to a designated reference date (``refdate``). A typical internal time step (``dt``) is 3600 seconds (1 hour). As the model progresses, it exports user-defined variables (``output_vars``) to a NetCDF file (default: ``aeolis.nc``, defined by ``output_file``) at customized intervals (``output_times``).

Aeolian Sediment Transport
~~~~~~~~~~~~~~~~~~
Detailed section: :ref:`aeolian-sediment-transport`.

Aeolian sediment transport is the core of the AeoLiS model. It is computed using a two-dimensional advection scheme, simplified here for one-dimensional transport of a single sediment fraction:

.. math::
   :label: advection_overview
           
   \frac{\partial c}{\partial t} + u_{\mathrm{sed}} \frac{\partial c}{\partial x} = \min \left ( \frac{\partial m_{\mathrm{a}}}{\partial t} \quad ; \quad \frac{c_{\mathrm{sat}} - c}{T} \right )

The saturated sediment concentration :math:`c_{\mathrm{sat}}` (``Cu``) defines the transport capacity, while :math:`c` (``Ct``) is the instantaneous concentration in the air. Transport is activated in the configuration file using ``process_transport``. The right-hand side of the advection equation is governed by the adaptation timescale :math:`T` (``T``), which determines how quickly the concentration reaches equilibrium. To allow sediment to actually erode from or deposit to the bed, ``process_bedupdate`` must be enabled. 

Solving this advection equation is one of the most computationally expensive parts of the model. You can choose different numerical approaches using the ``solver`` keyword. For detailed guidance on these options, see the :ref:`solver guide <solver-guide>`.

Several methods are available to compute the saturated sediment concentration (``method_transport``), as explained in the :ref:`saturated sediment transport <saturated-sediment-transport>` section. The equation by :cite:`Bagnold1937a` (``bagnold``) is the default:

.. math::
   :label: bagnold_overview

   c_{\mathrm{sat}} = \max \left ( 0 \quad ; \quad \alpha C \frac{\rho_{\mathrm{a}}}{g} \sqrt{\frac{d_{n}}{D_{n}}} \frac{\left ( u_* - u_{\mathrm{th}} \right )^3}{u_{\mathrm{sed}}} \right )

The sediment velocity :math:`u_{\mathrm{sed}}` (``u``, ``us``, ``un``) is determined by the ``method_grainspeed`` parameter. It can either be set equal to the governing wind speed (``windspeed``) or calculated using a saltation model (e.g., ``duran``). For more information on grain speed computations, see the :ref:`sediment velocity <sediment-velocity>` section.

.. note:: 
   For all vector variables (like ``uw``, ``ustar``,  ``tau``, ``u``,  ``q``), the subscripts ``s`` and ``n`` (e.g., ``uws``, ``uwn``) indicate the cross-shore and longshore directions, respectively, and the name without a subscript represents the overall magnitude. 

AeoLiS supports the inclusion of multiple sediment fractions (``grain_size``, ``grain_dist``) across multiple vertical layers (``nlayers``). This allows for the simulation of sediment sorting, mixing, and armoring. More details are provided in the :ref:`multi-fraction sediment transport <multi-fraction-sediment-transport>` section.

Enabling ``process_bedinteraction`` incorporates a bed interaction parameter into the advection equation. For more information, see the :ref:`bed interaction approach <bed-interaction-approach>` section.

Wind and Shear Velocity
~~~~~~~~~~~~~~~~~~
Detailed section: :ref:`wind-shear-velocity`

The shear velocity :math:`u_*` (``ustar``) acts as the primary driver of transport. It is initially computed for a flat bed using the Prandtl-Von Kármán Law of the Wall, based on wind velocity :math:`u_w` (``uw``) at a given height (provided via the ``wind_file``). This core wind process is required for all simulations and is enabled via ``process_wind``.

.. math::
   :label: lawofwall_overview

   u_* = \frac{u_w}{\ln \left( \frac{z}{z_0} \right)}\kappa

Topography steers the wind, causing perturbations in the shear stress :math:`\tau`, where :math:`\tau={\rho_a}{u_*}^2` (``tau``):

.. math::
   :label: topo_steering_overview

   \vec{\tau}(x,y)=\vec{\tau}_{0}+|\vec{\tau}_{0}|\delta\vec{\tau}(x,y)

This process can be activated using the ``process_shear`` keyword. Different methods are available to compute these shear perturbations, which can be selected through ``method_shear``. More information on topographic steering is given in the :ref:`wind-shear` section.

The presence of vegetation reduces the effective shear stress acting on the bed. This drag reduction is parameterized using the Raupach formulation, which relies on the vegetation-related roughness parameter :math:`\Gamma` (``gamma_veg``) and the basal cover :math:`\rho_{\mathrm{veg}}` (``rhoveg``):

Equation

For more information, see the :ref:`Vegetation <vegetation>` section.

Velocity threshold
~~~~~~~~~~~~~~~~~~

While shear velocity drives transport, the threshold velocity :math:`u_{\mathrm{th}}` (``uth``) serves as a limiter. It acts as a collective parameter for multiple supply-limiting processes, scaling the base threshold :math:`u_{\mathrm{*th,0}}` (``uth0``) by various environmental factors. Threshold calculations are enabled via ``process_threshold``.

.. math::
  :label: threshold_overview
  
  u_{\mathrm{* th}} = u_{\mathrm{* th, 0}} \cdot f_{\mathrm{M}} \cdot f_{\mathrm{R}} \cdot f_{\mathrm{S}}

The base threshold :math:`u_{\mathrm{* th, 0}}` is computed based on the local grain size and density (activated via ``th_grainsize``). This base value is then scaled by supply-limiting factors depending on the enabled model processes. The influence of surface moisture (:math:`f_{\mathrm{M}}`) is activated via ``th_moisture``, while sheltering by non-erodible roughness elements (:math:`f_{\mathrm{R}}`) is configured via ``th_sheltering``. The restricting effect of a non-erodible layer can be included using ``th_nelayer``. For more information, see the :ref:`Velocity threshold <shear-threshold>` section.

Vegetation
~~~~~~~~~~~~~~~~~~

AeoLiS simulates the dynamic growth and spreading of vegetation, enabled via ``process_vegetation``. The specific vegetation formulation is selected using ``method_vegetation`` (``duran`` for the original implementation, ``grass`` for the newest implementation (REF van Westen 2026)). In the original description, vegetation density rhoveg (``rhoveg``) is computed through the vegetation height hveg (``hveg``) w.r.t. the maximum height (``Hveg``):

Equation.

This density determines how much shear reduction (see equation ...). Vegetation growth is described by:

Equation dhveg

Growth is governed by specific parameters for intrinsic vertical growth :math:`V_{\mathrm{ver}}` (``V_ver``) and sensitivity to sediment burial :math:`\gamma_{\mathrm{veg}}` (``veg_gamma``). For more information, see the :ref:`Vegetation <vegetation>` section.

Hydrodynamics and Surface Moisture
~~~~~~~~~~~~~~~~~~

Water levels, wave runup, and groundwater can wet the beach, temporally increasing the shear velocity threshold, and mix sediment fractions. The model reads input water levels (``tide_file``) and wave heights (``wave_file``), activated via ``process_tide`` and ``process_wave``. The water level is first projected to the domain, resulting in the Still Water Level (``SWL``). After computing the wave runup (``R``), enabled via ``process_runup``), the Total Water Level (``TWL`` = ``SWL`` + ``R``) can be computed. The waterlevel ``zs`` maximum of bed level (``zb``) and TWL. Water depth is ``hw``. Applying masks (``tide_mask``, ``wave_mask``, ``runup_mask``) can be used to spatially modify the acting hydrodynamics. For more information on the hydrodynamics, see the ... section.

Inundation wets the bed, increasing the surface moisture (``moist``) and once the beach is exposed, infiltration and evaporation dry the surface. This moisture tracking is activated via ``process_moist``. A more advanced description of intertidal groundwater fluctuations by Hallin (2023) REF can be enabled through ``process_groundwater``. For more informuation on these computations see the :ref:`Surface moisture <surface-moisture>` section.

Wave impact can mix multiple sediment fractions over several bed layers down to the depth of disturbance , the ... section. (``process_mixtoplayer``)

Morphological Change
~~~~~~~~~~~~~~~~~~

Gradients in aeolian sediment transport can cause the bed level (``zb``) to change (``dzb``), enabled by ``process_bedupdate``:

Equation..

Computed from pickup rates (``pickup``), which is ...  More information in the :ref:`Morphological change <morphological-change>` section.

To redistributes sediment when the local slope exceeds the maximum angle of repose (``theta_dyn`` and ``theta_stat``), avalanching can be enabled through ``process_avalanche``. 

Submerged cells are subject to distinct bed level assumptions and marine erosion, managed by configurations like ``process_wet_bed_reset``. 


.. _fig-aeolis-overview:

.. figure:: /images/aeolis_overview.png
   :width: 900px
   :align: center

   Overview of the AeoLiS model


.. _aeolian-sediment-transport:

Aeolian Sediment Transport
-------------------

Calculating aeolian sediment transport is the core of the AeoLiS model. 
The model is based on the approach of :cite:`deVries2014a` which is extended to compute the
spatiotemporal varying sediment availability through simulation of the
process of beach armoring. For this purpose the bed is discretized in
horizontal grid cells and in vertical bed layers (2DV). Moreover, the
grain size distribution is discretized into fractions. This allows the
grain size distribition to vary both horizontally and vertically. A
bed composition module is used to compute the sediment availability
for each sediment fraction individually. This model approach is a
generalization of existing model concepts, like the shear velocity
threshold and critical fetch, and therefore compatible with these
existing concepts.

.. _advection-equation:

Advection Equation
^^^^^^^^^^^^^^^^^^^^^^^^

A 1D advection scheme is adopted in correspondence with
:cite:`deVries2014a` in which :math:`c` [:math:`\mathrm{kg/m^2}`] is
the instantaneous sediment mass per unit area in transport:

.. math::
   :label: advection
           
   \frac{\partial c}{\partial t} + u_z \frac{\partial c}{\partial x} = E - D

:math:`t` [s] denotes time and :math:`x` [m] denotes the cross-shore
distance from a zero-transport boundary. :math:`E` and :math:`D`
[:math:`\mathrm{kg/m^2/s}`] represent the erosion and deposition terms
and hence combined represent the net entrainment of sediment. Note
that Equation :eq:`advection` differs from Equation 9 in
:cite:`deVries2014a` as they use the saltation height :math:`h` [m]
and the sediment concentration :math:`C_{\mathrm{c}}`
[:math:`\mathrm{kg/m^3}`]. As :math:`h` is not solved for, the
presented model computes the sediment mass per unit area :math:`c = h
C_{\mathrm{c}}` rather than the sediment concentration
:math:`C_{\mathrm{c}}`. For conciseness we still refer to :math:`c` as
the *sediment concentration*.

The net entrainment is determined based on a balance between the
equilibrium or saturated sediment concentration
:math:`c_{\mathrm{sat}}` [:math:`\mathrm{kg/m^2}`] and the
instantaneous sediment transport concentration :math:`c` and is
maximized by the available sediment in the bed :math:`m_{\mathrm{a}}`
[:math:`\mathrm{kg/m^2}`] according to:

.. math::
   :label: erodep
           
   E - D = \min \left ( \frac{\partial m_{\mathrm{a}}}{\partial t} \quad ; \quad \frac{c_{\mathrm{sat}} - c}{T} \right )

:math:`T` [s] represents an adaptation time scale that is assumed
to be equal for both erosion and deposition. A time scale of 1 second
is commonly used :cite:`deVries2014a`.

.. _saturated-sediment-transport:

Saturated Sediment Transport
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The equilibrium, or saturated, sediment concentration :math:`c_{\mathrm{sat}}` is computed using an
empirical sediment transport formulation (e.g. :cite:`Bagnold1937a`):

.. math::
   :label: equilibrum-transport
          
   q_{\mathrm{sat}} = \alpha C \frac{\rho_{\mathrm{a}}}{g} \sqrt{\frac{d_{\mathrm{n}}}{D_{\mathrm{n}}}} \left ( u_z - u_{\mathrm{th}} \right )^3

in which :math:`q_{\mathrm{sat}}` [kg/m/s] is the equilibrium or
saturated sediment transport rate and represents the sediment
transport capacity. :math:`u_z` [m/s] is the wind velocity at height :math:`z` [m]
and :math:`u_{\mathrm{th}}` the velocity threshold [m/s]. The properties of
the sediment in transport are represented by a series of parameters:
:math:`C` [--] is a parameter to account for the grain size distribution
width, :math:`\rho_{\mathrm{a}}` [:math:`\mathrm{kg/m^3}`] is the density of the
air, :math:`g` [:math:`\mathrm{m/s^2}`] is the gravitational constant,
:math:`d_{\mathrm{n}}` [m] is the nominal grain size and :math:`D_{\mathrm{n}}`
[m] is a reference grain size. :math:`\alpha` is a constant to account for
the conversion of the measured wind velocity to the near-bed shear
velocity following Prandtl-Von Kármán's Law of the Wall:
:math:`\left(\frac{\kappa}{\ln z / z'} \right)^3` in which :math:`z'` [m] is the
height at which the idealized velocity profile reaches zero and
:math:`\kappa` [-] is the Von Kármán constant.

The equilibrium sediment transport rate :math:`q_{\mathrm{sat}}` is
divided by the wind velocity :math:`u_z` to obtain a mass per unit
area (per unit width):

.. math::
   :label: equilibrium-conc
   
   c_{\mathrm{sat}} = \max \left ( 0 \quad ; \quad \alpha C \frac{\rho_{\mathrm{a}}}{g} \sqrt{\frac{d_{n}}{D_{n}}} \frac{\left ( u_z - u_{\mathrm{th}} \right )^3}{u_z} \right )

in which :math:`C` [--] is an empirical constant to account for
the grain size distribution width, :math:`\rho_{\mathrm{a}}`
[:math:`\mathrm{kg/m^3}`] is the air density, :math:`g` [:math:`\mathrm{m/s^2}`] is the
gravitational constant, :math:`d_{\mathrm{n}}` [m] is the nominal grain
size, :math:`D_{\mathrm{n}}` [m] is a reference grain size, :math:`u_z` [m/s] is
the wind velocity at height :math:`z` [m] and :math:`\alpha` [--] is a constant to
convert from measured wind velocity to shear velocity.

Note that at this stage the spatial variations in wind velocity are
not solved for and hence no morphological feedback is included in the
simulation. The model is initially intended to provide accurate
sediment fluxes from the beach to the dunes rather than to simulate
subsequent dune formation.


.. _sediment-velocity:
Sediment Transport Velocity
^^^^^^^^^^^^^^^^^^^^^^^^^^^

PLACEHOLDER TEXT


.. _multi-fraction-sediment-transport:

Multi-fraction Sediment Transport
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
is determined by a bed interaction parameter :math:`\zeta`.

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
is represented by the bed interaction parameter :math:`\zeta`. The effective
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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

Each layer in each grid cell describes a grain size distribution over
a predefined number of sediment fractions (Figure
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

.. _hydraulic-mixing:

Hydraulic Mixing
~~~~~~~~~~~~~~

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
:math:`f_{\Delta z_{\mathrm{d}}}` [-] that relates the depth of disturbance
directly to the local breaker height according to:

.. math::
   :label: disturbance_depth
   
   \Delta z_{\mathrm{d}} = f_{\Delta z_{\mathrm{d}}} \cdot \min \left ( H \quad ; \quad \gamma \cdot d \right )

in which the offshore wave height :math:`H` [m] is taken as the
local wave height maximized by a maximum wave height over depth ratio
:math:`\gamma` [-]. :math:`d` [m] is the water depth that is provided to the model
through an input time series of water levels. Typical values for
:math:`f_{\Delta z_{\mathrm{d}}}` are 0.05 to 0.4 and 0.5 for :math:`\gamma`.


.. _bed-interaction-approach:

Bed-interaction Approach (zeta)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

PLACEHOLDER



.. _wind-shear-velocity:

Wind and Shear Velocity
-----------------------------------

PLACEHOLDER: WIND READ FROM ..... MAIN PROCESS...

.. _shear-velocity:

Shear Velocity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

PLACEHOLDER: LAW OF THE WALL

.. _topographic-steering:

Topographic Steering
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To simulate the topographic steering effects on dunes, often Computational Fluid Dynamics (CFD) methods are used. However, these methods are computationally expensive, making them not ideal for long-term morphodynamic simulations. To reduce computational costs, the topographic steering of the wind due to smooth gradients is implemented following an analytical perturbation theory for turbulent boundary layer flow :cite:`weng1991air, kroy2002minimal`. This method describes the topographic impact through perturbations in the shear stress :math:`\tau` [:math:`\mathrm{N/m^2}`] (where :math:`\tau={\rho_a}{u_*}^2`):

.. math::
   :label: shear_perturbation_base

   \vec{\tau}(x,y)=\vec{\tau}_{0}+|\vec{\tau}_{0}|\delta\vec{\tau}(x,y)

where :math:`\delta\vec{\tau}(x,y)` is the shear stress perturbation and :math:`\vec{\tau}_{0}` is the computed shear stress on a flat topography. 

For two-dimensional situations, the shear stress perturbation in the x- and y-directions (:math:`\delta\tau_{x}` and :math:`\delta\tau_{y}`) is computed in Fourier space according to the following equations:

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

where :math:`L` [m] is the typical length scale of the hill.

.. tip:: 
   **Tuning the length scale (L):** The typical length scale of the hill (L) significantly influences the shear perturbation. A smaller L will result in...

.. admonition:: Modelling Advice
   When defining the typical length scale of the hill (L), keep in mind that it acts as a smoothing parameter for the topography. If your grid resolution is very fine, setting L too low might...

For one-dimensional situations, a simplified solution of the shear perturbation approach is implemented. By ignoring some minor terms, it provides a less computationally expensive approach :cite:`kroy2002minimal`:

.. math::
   :label: shear_pert_1d

   \delta \tau =\alpha \int_{-\infty}^{\infty}d\xi\frac{\frac{\delta z_b}{\delta x}(x-\xi)}{\pi \xi}+\beta\frac{\delta z_b}{\delta x}(x)

where :math:`\alpha` [-] and :math:`\beta` [-] both depend on :math:`L/z_0`, but are user-defined fixed variables rather than computed in the model. :math:`\xi` [-] is the normalized cross-shore distance :math:`x/L`.

.. _flow-separation:

Flow separation
^^^^^^^^^^^^^^^

The implementation of the shear perturbation theory by :cite:`weng1991air` is only valid in situations with relatively smooth surfaces. The occurrence of steep slopes limits the validity of the approach. To address this, a description of flow separation is used following the Coastal Dune Model (CDM) :cite:`sauermann2001continuum, kroy2002minimal, DuranMoore2013`. A smooth envelope is created, which separates the main flow when a sharp edge is detected in the windward direction. This smooth envelope is called a separation bubble, :math:`z_{sep}` [m] (Figure :numref:`fig-concept-topo-steering`). This separation bubble represents the surface that divides the region of flow reversal from the main flow stream along the smooth hill. Subsequently, in all cells for which the bed level is lower than the separation bubble (:math:`z_b < z_{sep}`), the shear velocity :math:`u_{*}` is set to 0 m/s. This assumes that eventual flow reversal velocities are not significant enough to initiate aeolian transport.

The separation bubble surface :math:`z_{sep}` is modelled by a third-order polynomial. The height of the brinkline, or the location where the separation bubble starts to detach from the bed, is defined by :math:`z_b(x_{\mathrm{brink}}) \equiv z_{\mathrm{brink}}`. Assuming a maximum slope :math:`c` [deg] for the separation surface that determines the shape of the bubble, the reattachment length :math:`l_r` is obtained by:

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
~~~~~~~~~~~~~~~~~

The underlying implementation of the perturbation theory and separation bubble originally allows only for wind conditions that are perpendicular to the grid. An overlaying computational grid is introduced in AeoLiS, which rotates with the changing wind direction per time step. By doing this, the shear stresses are always estimated in the positive x-direction of the computational grid. The following steps are executed for each time step:

1. Create a computational ('Rotational') grid aligned with the wind direction (``set_computational_grid``).
2. Add and fill a buffer around the original ('Primary') grid.
3. Populate the computational grid by rotating it to the current wind direction and interpolate the original topography onto it. 
4. Compute the morphology-wind induced shear stress by using the perturbation theory.
5. Add the wind-induced shear stresses to the computational grid.
6. Rotate both the grids and the total shear stress results in the opposite direction.
7. Interpolate the total shear stress results from the computational grid to the original grid.
8. Rotate the wind shear stress results and the original grid back to the original orientation.

.. note:: 
   The extra rotations in the last two steps are necessary as a simplified, but faster in terms of computational time, interpolation method is used.

.. _vid-rotating-shear:

.. video:: /images/rotating_shear.mp4
   :autoplay:
   :loop:
   :muted:
   :width: 100%

   Animation demonstrating the rotational computational grid aligning with the changing wind direction to solve for topographic steering at each time step.



.. _shear-velocity-threshold:

Shear Velocity Threshold
------------------------

The shear velocity threshold represents the influence of bed surface
properties in the saturated sediment transport equation. The shear
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
  f_{u_{\mathrm{* th}}, \mathrm{S}} \cdot 
  u_{\mathrm{* th, 0}}

.. _base-threshold-grainsize:

Base Threshold: Grainsize
^^^^^^^^^^^^^^^

The initial shear velocity threshold :math:`u_{\mathrm{* th, 0}}` [m/s] is
computed based on the grain size following :cite:`Bagnold1937b`:

.. math::
   :label: shear

   u_{\mathrm{* th, 0}} = A \sqrt{ \frac{\rho_{\mathrm{p}} - \rho_{\mathrm{a}}}{\rho_{\mathrm{a}}} \cdot g \cdot d_{\mathrm{n}}}

where :math:`A` [-] is an empirical constant, :math:`\rho_{\mathrm{p}}`
[:math:`\mathrm{kg/m^3}`] is the grain density, :math:`\rho_{\mathrm{a}}`
[:math:`\mathrm{kg/m^3}`] is the air density, :math:`g` [:math:`\mathrm{m/s^2}`] is the
gravitational constant and :math:`d_{\mathrm{n}}` [m] is the nominal grain
size of the sediment fraction.

.. _threshold-moisture-content:

Threshold: Moisture content
^^^^^^^^^^^^^^^^^

The shear velocity threshold is updated based on moisture content
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

Threshold: Roughness elements
^^^^^^^^^^^^^^^^^^^

Sediment sorting may lead to the emergence of non-erodible elements
from the bed. Non-erodible roughness elements may shelter the erodible
bed from wind erosion due to shear partitioning, resulting in a
reduced sediment availability :cite:`Raupach1993`. Therefore the
equation of :cite:`Raupach1993` is implemented according to:

.. math::
   :label: raupach
           
   u_{\mathrm{* th, R}} = u_{\mathrm{* th}} \cdot \sqrt{ \left( 1 - m \cdot \sum_{k=k_0}^{n_{\mathrm{k}}}{w_k^{\mathrm{bed}}} \right) \left( 1 + \frac{m \beta}{\sigma} \cdot \sum_{k=k_0}^{n_{\mathrm{k}}}{w_k^{\mathrm{bed}}} \right) }

in which :math:`\sigma` is the ratio between the frontal area and the
basal area of the roughness elements and :math:`\beta` is the ratio
between the drag coefficients of the roughness elements and the bed
without roughness elements. :math:`m` is a factor to account for the
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

Threshold: Non-erodible layer
^^^^^^^^^^^^




.. _vegetation:

Vegetation 
------------

The description of the implementation of vegetation dynamics in AeoLiS is based on :cite:`Strypsteen2024`.

In AeoLiS, the influence of vegetation on dune evolution is comprehensively 
addressed. This includes modelling the intrinsic growth of vegetation, 
accounting for factors such as growth and decay due to burial :cite:`DuranMoore2013`, 
lateral expansion and establishment :cite:`Keijsers2016`, as well as simulating 
the destruction of vegetation caused by hydrodynamic processes.In the event of 
cell inundation, vegetation density is reduced as a result. 

.. _vegetation-metrics:

Vegetation Metrics
^^^^^^^^^^^^

The vegetation density :math:`\rho_{\text{veg}}` can vary in space and time and is determined by the ratio of the actual vegetation height (hveg) 
to the maximum vegetation height (:math:`h_{\text{veg,max}}`), and can vary between 0 and 1 (:cite:`DuranHerrmann2006`):

.. math::
   :label: Vegetation_density

   \rho_{\text{veg}} = \left( \frac{h_{\text{veg}}}{h_{\text{veg,max}}} \right)^2


This assumption is based on the idea that burying vegetation reduces its height, which indicates a decrease in 
actual cover. The change in vegetation density per grid cell is directly linked to the alteration in vegetation 
height within that specific cell. This height variation is influenced by both the growth rate of the vegetation and 
the rate of sediment burial. If the vegetation density remains constant over time, it suggests either no 
sedimentation or a growth rate equal to the rate of sediment burial within the cell. 


.. _vegetation-development:

Vegetation Development
^^^^^^^^^^^^

Vegetation growth and decay follow the model proposed by :cite:`DuranHerrmann2006`, modified to include :math:`\delta z_{\text{b,opt}}`
(m/year), representing sediment burial for optimal growth that shifts the peak of optimal growth:

.. math::
   :label: changes_vegetation_height

   \frac{\delta h_{\text{veg}}}{\delta t} = V_{\text{ver}} \left(1 - \frac{h_{\text{veg}}}{h_{\text{veg,max}}}\right) - \gamma_{\text{veg}} \left| \frac{\delta z_{\text{b,veg}}}{\delta t} - \delta z_{\text{b,opt}} \right|

Here, :math:`\gamma_{\text{veg}}` (default = 1) is a sediment burial factor that accounts for the impact of 
sediment burial on vegetation. The height of the vegetation (:math:`h_{\text{veg}}` in m) cannot be less than zero. 
Vver represents the maximum vertical growth rate of vegetation given in m/year, while the sediment burial rate 
:math:`\delta z_{\text{b,veg}}` [m] is determined as the bed level change per time step. By simply converting this 
value to a bed level change per year multiple errors are induced, as the time scale over which the bed level change 
actually occurs is much shorter than this one year. To compare the bed level change per time step with the 
vegetation growth rate per year, an average bed level change is estimated over a specified time (default is one 
day). This average is then extrapolated to an annual rate. This method ensures that sudden changes in the bed level 
change over one time step are not used as an estimate of the total bed level change in one year, which would be far 
too high.

The optimal growth rate for certain vegetation species in dune environments is depending upon sediment burial :cite:`Maun1998`. 
The optimal burial rate for maximum vegetation growth for marram grass for the neighbouring Dutch 
coast is around 0.31 m/year with a burying tolerance of 0.78 to 0.96 m burial/year :cite:`Nolet2018`. This 
optimal value is used in the model. :math:`V_{\text{ver}}` contains information of meteorological and local 
conditions that enhance or inhibit vegetation growth process :cite:`Danin1991`, :cite:`Hesp1991`. 

.. _fig-Veg_growth:

.. figure:: /images/Veg_growth.png
   :width: 600px
   :align: center

   A) The vegetation growth response varies with different vertical growth rates (example for Vver = 1 and 2 m/year). Optimal vegetation growth is determined by a burial rate of 0.31 
   m/year, with a maximum vegetation height set at 1 m and a plant height of 0.5 m. Additionally, the growth response for varying burial factors is depicted (:math:`\lambda_{\text{veg}}`
   = 1 and 2). B) Shear stress reduction for two different vegetation-related roughness parameters and vegetation densities (:math:`\Gamma` = 16 and 32).

Vegetation can begin to grow through lateral propagation or random germination. Once established, it can continue 
to grow and spread laterally. The uncertainties associated with random germination are handled on a cell-by-cell 
basis using a probabilistic approach, similar to the cellular automata method described by Keijsers et al. (2016). 
AeoLiS incorporates a germination probability, denoted as ρ_{ger}, for each grid cell. This probability is constant 
across the domain, except in eroding grid cells (where bed elevation decreases), where ρ_{ger} is set to 0. Lateral 
propagation is determined by identifying the boundaries between vegetated and non-vegetated cells, with the 
parameter ρ_{lat} adjusting the likelihood of lateral propagation at these boundaries.



.. _vegetation-induced-shear-reduction:

Vegetation-induced Shear Reduction
^^^^^^^^^^^^^^^^^^^^^^^^
Inspired by the Coastal Dune Model (CDM) proposed by :cite:`DuranMoore2013`, AeoLiS 
incorporates vegetation-wind interaction using the expression established by :cite:`DuranHerrmann2006`:

.. math::
   :label: shear_reduction_vegetation

   \frac{u_{\text{veg}}}{u_*} = \frac{1}{\sqrt{1 + \Gamma \rho_{\text{veg}}}}


where the ratio of shear velocity in the presence of vegetation (:math:`u_{*,\text{veg}}`) to the unobstructed 
shear velocity (:math:`u_*`) is determined by a vegetation-related roughness parameter (:math:`\Gamma`) and the 
vegetation density within a unit area of the grid cell (:math:`\rho_{\text{veg}}`). In the model, :math:`\Gamma` = 16 is derived 
from plant form drag and geometry values documented for creosote communities :cite:`DuranHerrmann2006`. 
This implementation calculates the expression on each model grid cell, with higher vegetation density 
(expressed by :math:`\rho_{\text{veg}}`) leading to a more substantial reduction in shear velocity compared to sparse 
vegetation. By integrating these physical and ecological processes, AeoLiS simulates spatial patterns and temporal 
variations in sediment transport and morphological changes resulting from aeolian processes in coastal 
environments.


.. _vegetation-computing-zeta:

Computing bed-interaction (zeta) over vegetation
^^^^^^^^^^^^^^^^^^^^^^^^

PLACEHOLDER

.. _vegetation-mortality:

Vegetation Mortality
^^^^^^^^^^^^^^^^^^^^^^^^

PLACEHOLDER



.. _hydrodynamics-moisture:

Hydrodynamics and Surface Moisture
-----------------

Aeolis computes .... for ... and ... and ...


.. _water-levels-waves-run-up:

Water levels, Waves and Run-up
^^^^^^^^^^^^^^^

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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

PLACEHOLDER


.. _hydraulic-sediment-mixing:

Hydraulic Sediment Mixing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

PLACEHOLDER

.. _marine-driven-bed-level-change:

Marine-driven Bed Level Change
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

PLACEHOLDER



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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~
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
~~~~~~~~~~~~
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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

.. _bed-level-update:

Bed Level Update
^^^^^^^^^^^^

PLACEHOLDER

.. _avalanching:
Avalanching
^^^^^^^^^^^^

PLACEHOLDER

Marine-driven Morphodynamics
^^^^^^^^^^^^^^^^

PLACEHOLDER: See section...

