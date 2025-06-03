.. _model:

Model description
=================

The model approach of :cite:`deVries2014a` is extended to compute the
spatiotemporal varying sediment availability through simulation of the
process of beach armoring. For this purpose the bed is discretized in
horizontal grid cells and in vertical bed layers (2DV). Moreover, the
grain size distribution is discretized into fractions. This allows the
grain size distribition to vary both horizontally and vertically. A
bed composition module is used to compute the sediment availability
for each sediment fraction individually. This model approach is a
generalization of existing model concepts, like the shear velocity
threshold and critical fetch, and therefore compatible with these
existing concepts..

Advection Scheme
----------------

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

The saturated sediment concentration :math:`c_{\mathrm{sat}}` is computed using an
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

Multi-fraction Erosion and Deposition
-------------------------------------

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

Simulation of Sediment Sorting and Beach Armoring
-------------------------------------------------

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

Simulation of the Emergence of Non-erodible Roughness Elements
--------------------------------------------------------------

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

Simulation of the Hydraulic Mixing
----------------------------------

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

Simulation of surface moisture
------------------------------

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


Runup and wave setup
^^^^^^^^^^^^^^^^^^^^
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


Tide- and wave-induced groundwater variations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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


Capillary rise
^^^^^^^^^^^^^^
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
^^^^^^^^^^^^
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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
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


Shear velocity threshold
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

Moisture content
^^^^^^^^^^^^^^^^

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


Roughness elements
^^^^^^^^^^^^^^^^^^

The shear velocity threshold is updated based on the presence of
roughness elements following :cite:`Raupach1993`:

.. math::
   :label: shear-rough

  f_{u_{\mathrm{* th},R}} = \sqrt{(1 - m \cdot \sum_{k=k_0}^{n_k}{\hat{w}_k^{\mathrm{bed}}})
    (1 + \frac{m \beta}{\sigma} \cdot \sum_{k=k_0}^{n_k}{\hat{w}_k^{\mathrm{bed}}})}

by assuming:

.. math::
   :label: lambda-rough
   
  \lambda = \frac{\sum_{k=k_0}^{n_k}{\hat{w}_k^{\mathrm{bed}}}}{\sigma}

where :math:`f_{u_{\mathrm{* th},R}}` [-] is a factor in Equation
:eq:`apx-shearvelocity`, :math:`k_0` is the sediment fraction index of
the smallest non-erodible fraction in current conditions and :math:`n_k` is
the number of sediment fractions defined. The implementation is
discussed in detail in section \ref{sec:roughness}.

Salt content
^^^^^^^^^^^^

The shear velocity threshold is updated based on salt content
following :cite:`Nickling1981`:

.. math::
   :label: salt-rough
   
   f_{u_{\mathrm{* th}},S} = 1.03 \cdot \exp(0.1027 \cdot p_{\mathrm{s}})

where :math:`f_{u_{\mathrm{* th},S}}` [-] is a factor in Equation
:eq:`apx-shearvelocity` and :math:`p_{\mathrm{s}}` [-] is the salt
content [mg/g]. Currently, no model is implemented that predicts the
instantaneous salt content. The spatial varying salt content needs to
be specified by the user, for example through the BMI interface.

