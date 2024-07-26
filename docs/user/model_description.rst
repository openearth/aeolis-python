.. _model_description:

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

Numerical implementation of the advection equation
--------------------------------------------------

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

The formulation is discretized in different ways to allow for different types of simulations balancing accuracy vs. computational resources. The conservative method combined with an euler backward scheme (written by Prof. Rauwoens) is the current default for most simulations. Non-conservative methods end explicit Euler forward schemes are also available. 

Default scheme -- Conservative Euler Backward Implicit
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The default numerical method assumes the advection scheme in a conservative form in combination with an euler backward scheme. This scheme is prepared to use a TVD method but this is not implemented yet (add footnote{Total Variance Diminishing, this is explained in the lecture notes by Zijlema p94}) 

The fluxes at the interface of the cells are defined used in the advection terms:

.. math::
   :label: conservative1
   
   \frac{c^{n+1}_{i,j,k} - c^n_{i,j,k}}{\Delta t} + \\
   \frac{u_{\text{x},i+1/2,j} \cdot c^{n+1}_{i+1/2,j,k} - u_{\text{x},i-1/2,j} \cdot c^{n+1}_{i-1/2,j,k}}{\Delta x} + \\
   \frac{u_{\text{y},i,j+1/2} \cdot c^{n+1}_{i,j+1/2,k} - u_{\text{y},i,j-1/2} \cdot c^{n+1}_{i,j-1/2,k}}{\Delta y} +
   \\ = \\
   \frac{\min(\hat{w}^{n+1}_{i,j,k} \cdot c^{n+1}_{\mathrm{sat},i,j,k},m_{i.j.k} + c^{n+1}_{i,j,k} ) - c^{n+1}_{i,j,k}}{T}
   
In which :math:`n` is the time step index, :math:`i` and :math:`j` are the cross-shore and alongshore spatial grid cell indices and :math:`k` is the grain size fraction index. :math:`w` [-] is the weighting factor used for the weighted addition of the saturated sediment concentrations over all grain size fractions. Note that u is spatially varying but has no temporal index. This is because u is a result of a separate wind solver and considered temporally invariant in the advection solver. 

Now we use a correction algorithm where:

.. math::
   :label: conservative2
   
   c^{n+1}_{i,j,k} = c^{n+1 *}_{i,j,k} + \delta c_{i,j,k}
   
where :math:`\delta c_{i,j,k}` is solved for and :math:`*` denotes the previous iteration. 

When now assuming an upwind scheme in space, we can derive 4 concentrations at the cell faces which are dependent on the velocity at the cell faces.

We assume in x direction:

.. math::

    c^{n+1}_{i+1/2,j,k} =
    \begin{cases}
        c^{n+1 *}_{i,j,k} + \delta c_{i,j,k} & \text{if $u_{\text{x},i+1/2,j} > 0$,}\\
        c^{n+1 *}_{i+1,j,k} + \delta c_{i+1,j,k} & \text{if $u_{\text{x},i+1/2,j} < 0$.}
    \end{cases}

.. math::

    c^{n+1}_{i-1/2,j,k} =
    \begin{cases}
        c^{n+1 *}_{i-1,j,k} + \delta c_{i-1,j,k} & \text{if $u_{\text{x},i-1/2,j} > 0$,}\\
        c^{n+1 *}_{i,j,k} + \delta c_{i,j,k} & \text{if $u_{\text{x},i-1/2,j} < 0$.}
    \end{cases}

and in y-direction:

.. math::

    c^{n+1}_{i,j+1/2,k} =
    \begin{cases}
        c^{n+1 *}_{i,j,k} + \delta c_{i,j,k} & \text{if $u_{\text{y},i,j+1/2} > 0$,}\\
        c^{n+1 *}_{i,j+1,k} + \delta c_{i,j+1,k} & \text{if $u_{\text{y},i,j+1/2} < 0$.}
    \end{cases}

.. math::

    c^{n+1}_{i,j-1/2,k} =
    \begin{cases}
        c^{n+1 *}_{i,j-1,k} + \delta c_{i,j-1,k} & \text{if $u_{\text{y},i,j-1/2} > 0$,}\\
        c^{n+1 *}_{i,j,k} + \delta c_{i,j,k} & \text{if $u_{\text{y},i,j-1/2} < 0$.}
    \end{cases}

Now we assume:

- :math:`\Gamma_x = 1` if :math:`u_{\text{x},i+1/2,j,k} > 0` and :math:`\Gamma_x = 0` if :math:`u_{\text{x},i+1/2,j,k} \leq 0`
- :math:`\Gamma_y = 1` if :math:`u_{\text{y},i,j+1/2,k} > 0` and :math:`\Gamma_x = 0` if :math:`u_{\text{y},i,j+1/2,k} \leq 0`
   
(We did not test if this works well with diverging and converging flows. We may need another term that describes the conditions at the negative cell faces if they are of opposite direction than the positive cell faces and vice versa)

Let's continue for the moment so that

.. math::

   \begin{gathered}
   \frac{c^{n+1 *}_{i,j,k} + \delta c_{i,j,k} - c^n_{i,j,k}}{\Delta t} + \\
   \Gamma_x \cdot \frac{u_{\text{x},i+1/2,j} \cdot (c^{n+1 *}_{i,j,k} + \delta c_{i,j,k}) - u_{\text{x},i-1/2,j} \cdot (c^{n+1 *}_{i-1,j,k} + \delta c_{i-1,j,k})}{\Delta x} + \\
   (1-\Gamma_x) \cdot \frac{u_{\text{x},i+1/2,j} \cdot (c^{n+1 *}_{i+1,j,k} + \delta c_{i+1,j,k}) - u_{\text{x},i-1/2,j} \cdot (c^{n+1 *}_{i,j,k} + \delta c_{i,j,k})}{\Delta x} + \\
   \Gamma_y \cdot \frac{u_{\text{y},i,j+1/2} \cdot (c^{n+1 *}_{i,j,k} + \delta c_{i,j,k}) - u_{\text{y},i,j-1/2} \cdot (c^{n+1 *}_{i,j-1,k} + \delta c_{i,j-1,k})}{\Delta y} + \\
   (1-\Gamma_y) \cdot \frac{u_{\text{y},i,j+1/2} \cdot (c^{n+1 *}_{i,j+1,k} + \delta c_{i,j+1,k}) - u_{\text{y},i,j-1/2} \cdot (c^{n+1 *}_{i,j,k} + \delta c_{i,j,k})}{\Delta y} +
   \\ = \\
   \frac{\min(\hat{w}^{n+1}_{i,j,k} \cdot c^{n+1}_{\mathrm{sat},i,j,k},m_{i,j,k}+c^{n+1 *}_{i,j,k} + \xcancel{\delta c_{i,j,k}}) - c^{n+1 *}_{i,j,k} + \delta c_{i,j,k}}{T}
   \end{gathered}

(note that the above does not take converging and diverging flows into account, also :math:`\delta c_{i,j,k}` at the right hand side in the "min" brackets is difficult to solve for. In the code, this term is neglected which may cause some inaccuracy when calculating pickup. Although mass continuity is corrected for in the implicit scheme when calculating pickup using equation ???)

Now we simplify:

.. math::

   \begin{gathered}
    (\frac{\Delta x \Delta y}{\Delta t} + \Gamma_x\Delta y \cdot u_{\text{x},i+1/2,j} - (1-\Gamma_x)\Delta y \cdot u_{\text{x},i-1/2,j} + \Gamma_y\Delta x \cdot u_{\text{y},i,j+1/2}\\ - (1-\Gamma_y)\Delta x \cdot u_{\text{y},i,j-1/2}+\frac{\Delta x \Delta y}{T_s})\cdot \delta c_{i,j,k}\\
    -(\Gamma_x\Delta y \cdot u_{\text{x},i-1/2,j})\cdot \delta c_{i-1,j,k}\\
    +((1-\Gamma_x)\Delta y \cdot u_{\text{x},i+1/2,j})\cdot \delta c_{i+1,j,k}\\ 
    -(\Gamma_y\Delta x \cdot u_{\text{y},i,j-1/2})\cdot \delta c_{i,j-1,k}\\ 
    + ((1-\Gamma_y)\Delta x \cdot u_{\text{y},i,j+1/2})\cdot \delta c_{i,j+1,k}\\ 
   \end{gathered}

or

.. math::

   \begin{gathered}
    A0 \cdot \delta c_{i,j,k}
    + A\text{m1} \cdot \delta c_{i-1,j,k}
    + A\text{p1} \cdot \delta c_{i+1,j,k}\\ 
    + A\text{mx} \cdot \delta c_{i,j-1,k} 
    + A\text{px} \cdot \delta c_{i,j+1,k} = y_{i,j,k}
    \label{eq:lin}
   \end{gathered}

or the linear system of equations in general form:

.. math::
   :label: conservative_general

    A \cdot \delta c_{i,j,k} = y_{i,j,k}

Where :math:`A` is a 3-dimensional sparse matrix that is compiled using the matrix diagonals (:math:`A0, Am1, Ap1, Amx, Apx`) which are defined as: 

.. math:: 

   \begin{aligned}
   A0 = & +\frac{\Delta x \Delta y}{\Delta t} \\
    & +\frac{\Delta x \Delta y}{T_s} \\
    & - (1-\Gamma_x)\Delta y \cdot u_{\text{x},i-1/2,j} \\
    & + \Gamma_x\Delta y \cdot u_{\text{x},i+1/2,j} \\
    & - (1-\Gamma_y)\Delta x \cdot u_{\text{y},i,j-1/2} \\
    & + \Gamma_y\Delta x \cdot u_{\text{y},i,j+1/2} \\
   \end{aligned}

and

.. math::

    A\text{m}1 = -\Gamma_x\Delta y \cdot u_{\text{x},i-1/2,j}

and

.. math::

    A\text{p}1 = (1-\Gamma_x)\Delta y \cdot u_{\text{x},i+1/2,j} 

and

.. math::
    A\text{mx} = -\Gamma_y\Delta x \cdot u_{\text{y},i,j-1/2}

and

.. math::

    A\text{px} = (1-\Gamma_y)\Delta x \cdot u_{\text{y},i,j+1/2}

Let's go towards the RHS

.. math::

   \begin{aligned}  
       y_{i,j,k} = & - \frac{\Delta x \Delta y}{\Delta t}(c^{n+1 *}_{i,j,k}-c^{n}_{i,j,k}) \\
       & + \frac{\Delta x \Delta y}{T_s}(\min(\hat{w}^{n+1}_{i,j,k} \cdot c^{n+1}_{\mathrm{sat},i,j,k},m_{i,j,k} + c^{n+1 *}_{i,j,k}) - c^{n+1 *}_{i,j,k}) \\
       & + \Delta y \cdot u_{\text{x},i-1/2,j} \cdot 
       (\Gamma_x \cdot c^{n+1 *}_{i-1,j,k} + (1-\Gamma_x) c^{n+1 *}_{i,j,k})\\
       & - \Delta y \cdot u_{\text{x},i+1/2,j} \cdot (\Gamma_x \cdot c^{n+1 *}_{i,j,k} + (1-\Gamma_x) c^{n+1 *}_{i+1,j,k})\\
       & + \Delta x \cdot u_{\text{y},i,j-1/2} \cdot (\Gamma_y \cdot c^{n+1 *}_{i,j-1,k} + (1-\Gamma_y) c^{n+1 *}_{i,j,k})\\
       & - \Delta x \cdot u_{\text{y},i,j+1/2} \cdot (\Gamma_y \cdot c^{n+1 *}_{i,j,k} + (1-\Gamma_y) c^{n+1 *}_{i,j+1,k})\\ 
   \end{aligned}

in the python code some intermediate variable is defined to make it easier to shift indexes

.. math::

    \text{Ctxfs} =(\Gamma_x \cdot c^{n+1 *}_{i,j,k} + (1-\Gamma_x) c^{n+1 *}_{i+1,j,k})

and

.. math::

    \text{Ctxfn} =(\Gamma_y \cdot c^{n+1 *}_{i,j,k} + (1-\Gamma_y) c^{n+1 *}_{i,j+1,k})

also Erosion and deposition are defined using seperate variables.

.. math::
    D_{i,j,k}= \frac{\Delta x \Delta y}{T_s}c^{n+1 *}_{i,j,k}

and 

.. math::

    A_{i,j,k}= \frac{\Delta x \Delta y}{T_s}m_{i,j,k} + D_{i,j,k}

and

.. math::
    U_{i,j,k}= \frac{\Delta x \Delta y}{T_s}\hat{w}^{n+1}_{i,j,k} \cdot c^{n+1}_{\mathrm{sat},i,j,k}

and 

.. math::
    E_{i,j,k} = \min(U_{i,j,k},A_{i,j,k})

After solving equation :math:`\delta c_{i,j,k}` using :eq:`conservative_general`, :math:`c^{n+1}_{i,j,k}` can be calculated using equation :eq:`conservative2`.

Also, the pickup per grid cell can be calculated using:

.. math::

   \text{pickup} = \frac{\hat{w}^{n+1}_{i,j,k} \cdot c^{n+1}_{\mathrm{sat},i,j,k}- c^{n+1}_{i,j,k}}{T_s}\Delta t
    \label{eq:pickup}

note that this is only valid when using an Euler backward scheme. 

Solving the Linear System of Equations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The linear system of equations can be elaborated :

.. math::
  :label: apx-system
  
  \left[
    \begin{array}{cccccc}
      A^0_1      & A^{1}_1    & \textbf{0} & \cdots       & \textbf{0}    & A^{n_{\mathrm{y}}+1}_1 \\
      A^{-1}_2   & A^0_2      & \ddots     & \ddots       &               & \textbf{0} \\
      \textbf{0} & \ddots     & \ddots     & \ddots       & \ddots        & \vdots     \\
      \vdots     & \ddots     & \ddots     & \ddots       & \ddots        & \textbf{0} \\
      \textbf{0} &            & \ddots     & \ddots       & A^0_{n_{\mathrm{y}}}      & A^1_{n_{\mathrm{y}}}   \\
      A^{-n_{\mathrm{y}}-1}_{n_{\mathrm{y}}+1} & \textbf{0} & \cdots     & \textbf{0}   & A^{-1}_{n_{\mathrm{y}}+1} & A^0_{n_{\mathrm{y}}+1} \\
    \end{array}
  \right] \left[
    \begin{array}{c}
      \vec{\delta c}_1 \\ \vec{\delta c}_2 \\ \vdots \\ \vdots \\ \vec{\delta c}_{n_{\mathrm{y}}} \\ \vec{\delta c}_{n_{\mathrm{y}}+1} \\
    \end{array} 
  \right] = \left[ 
    \begin{array}{c}
      \vec{y}_1 \\ \vec{y}_2 \\ \vdots \\ \vdots \\ \vec{y}_{n_{\mathrm{y}}} \\ \vec{y}_{n_{\mathrm{y}}+1} \\
    \end{array} 
  \right]
    
where each item in the matrix is again a matrix :math:`A^l_j` and
each item in the vectors is again a vector :math:`\vec{\delta c}_j` and
:math:`\vec{y}_j` respectively. The form of the matrix :math:`A^l_j` depends on
the diagonal index :math:`l` and reads:

.. math::
  :label: apx-diagonal
   
  A^0_j = 
  \left[
    \begin{array}{ccccccc}
      0              & 0               & 0                & 0
      & \cdots           & \cdots           & 0                 \\
      a^{0,-1}_{2,j} & a^{0,0}_{2,j}    & a^{0,1}_{2,j}    & \ddots
      &                  &                  & \vdots            \\
      0              & a^{0,-1}_{3,j}   & a^{0,0}_{3,j}    & a^{0,1}_{3,j}
      & \ddots           &                  & \vdots            \\
      \vdots         & \ddots           & \ddots           & \ddots
      & \ddots           & \ddots           & \vdots            \\
      \vdots         &                  & \ddots           & a^{0,-1}_{n_{\mathrm{x}}-1,j}
      & a^{0,0}_{n_{\mathrm{x}}-1,j} & a^{0,1}_{n_{\mathrm{x}}-1,j} & 0                 \\
      \vdots         &                  &                  & 0
      & a^{0,-1}_{n_{\mathrm{x}},j}  & a^{0,0}_{n_{\mathrm{x}},j}   & a^{0,1}_{n_{\mathrm{x}},j}    \\
      0              & \cdots           & \cdots           & 0
      & 1                & -2               & 1                 \\
    \end{array}
  \right]

for :math:`l = 0` and 

.. math::
  :label: apx-offdiagonal
   
  A^l_j = 
  \left[
    \begin{array}{ccccccc}
      1               & 0                & \cdots           & \cdots
      & \cdots           & \cdots           & 0                 \\
      0               & a^{l,0}_{2,j}    & \ddots           &
      &                  &                  & \vdots            \\
      \vdots          & \ddots           & a^{l,0}_{3,j}    & \ddots
      &                  &                  & \vdots            \\
      \vdots          &                  & \ddots           & \ddots
      & \ddots           &                  & \vdots            \\
      \vdots          &                  &                  & \ddots
      & a^{l,0}_{n_{\mathrm{x}}-1,j} & \ddots           & \vdots            \\
      \vdots          &                  &                  &
      & \ddots           & a^{l,0}_{n_{\mathrm{x}},j}   & 0                 \\
      0               & \cdots           & \cdots           & \cdots  
      & \cdots           & 0                & 1                 \\
    \end{array}
  \right]

for :math:`l \neq 0`. The vectors :math:`\vec{\delta c}_{j,k}` and :math:`\vec{y}_{j,k}`
read:

.. math::
  :label: c-array

  \begin{array}{rclrcl}
    \vec{\delta c}_{j,k} &=& \left[ 
      \begin{array}{c}
        \delta c^{n+1}_{1,j,k} \\
        \delta c^{n+1}_{2,j,k} \\
        \delta c^{n+1}_{3,j,k} \\
        \vdots \\
        \delta c^{n+1}_{n_{\mathrm{x}}-1,j,k} \\
        \delta c^{n+1}_{n_{\mathrm{x}},j,k} \\
        \delta c^{n+1}_{n_{\mathrm{x}}+1,j,k} \\
    \end{array}
    \right] & ~ \mathrm{and} ~
    \vec{y}_{j,k} &=& \left[ 
      \begin{array}{c}
        0 \\
        y^n_{2,j,k} \\
        y^n_{3,j,k} \\
        \vdots \\
        y^n_{n_{\mathrm{x}}-1,j,k} \\
        y^n_{n_{\mathrm{x}},j,k} \\
        0 \\
      \end{array}
    \right] \\
    \end{array}

:math:`n_{\mathrm{x}}` and :math:`n_{\mathrm{y}}` denote the number of
spatial grid cells in x- and y-direction.

Iterations to solve for multiple fractions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The linear system defined in Equation :eq:`apx-system` is solved by a
sparse matrix solver for each sediment fraction separately in
ascending order of grain size. Initially, the weights
:math:`\hat{w}^{n+1}_{i,j,k}` are chosen according to the grain size
distribution in the bed and the air. The sediment availability
constraint is checked after each solve:

.. math::
  :label: solve

     m_{\mathrm{a}} \geq \frac{\hat{w}^{n+1}_{i,j,k} c^{n+1}_{\mathrm{sat},i,j,k} - c^{n+1}_{i,j,k}}{T} \Delta t^n

If the constraint if violated, a new estimate for the weights
is back-calculated following:

.. math::
  :label: solve-weights

  \hat{w}^{n+1}_{i,j,k} = \frac{ c^{n+1}_{i,j,k} + m_{\mathrm{a}} \frac{T}{\Delta t^n} }{c^{n+1}_{\mathrm{sat},i,j,k}}

The system is solved again using the new weights. This
procedure is repeated until a weight is found that does not violate
the sediment availability constraint. If the time step is not too
large, the procedure typically converges in only a few
iterations. Finally, the weights of the larger grains are increased
proportionally as to ensure that the sum of all weights remains
unity. If no larger grains are defined, not enough sediment is
available for transport and the grid cell is truly
availability-limited. This situation should only occur occasionally as
the weights in the next time step are computed based on the new bed
composition and thus will be skewed towards the large fractions. If
the situation occurs regularly, the time step is chosen too large
compared to the rate of armoring.

Euler Schemes in non-conservative form
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Early model results relied on Euler schemes in a non conservative form. This allowed for a relatively easy implementation but did not guarantee mass conservation. In version 2 of AEOLIS the conservative form became the default. However, some users still use the older scheme.

The formulation is discretized following a first order upwind scheme
assuming that the wind velocity :math:`u_z` is positive in both
x-direction and y-direction:

.. math::
  :label: apx-explicit
  
  \frac{c^{n+1}_{i,j,k} - c^n_{i,j,k}}{\Delta t^n} + 
  u^n_{z,\mathrm{x}} \frac{c^n_{i+1,j,k} - c^n_{i,j,k}}{\Delta x_{i,j}} + 
  u^n_{z,\mathrm{y}} \frac{c^n_{i,j+1,k} - c^n_{i,j,k}}{\Delta y_{i,j}} \\ = 
  \frac{\hat{w}^n_{i,j,k} \cdot c^n_{\mathrm{sat},i,j,k} - c^n_{i,j,k}}{T}

in which :math:`n` is the time step index, :math:`i` and :math:`j` are
the cross-shore and alongshore spatial grid cell indices and :math:`k`
is the grain size fraction index. :math:`w` [-] is the weighting
factor used for the weighted addition of the saturated sediment
concentrations over all grain size fractions.

The discretization can be generalized for any wind direction as:

.. math::
  :label: apx-explicit-generalized
   
  \frac{c^{n+1}_{i,j,k} - c^n_{i,j,k}}{\Delta t^n} + 
  u^n_{z,\mathrm{x+}} c^n_{i,j,k,\mathrm{x+}} + 
  u^n_{z,\mathrm{y+}} c^n_{i,j,k,\mathrm{y+}} \\ + 
  u^n_{z,\mathrm{x-}} c^n_{i,j,k,\mathrm{x-}} + 
  u^n_{z,\mathrm{y-}} c^n_{i,j,k,\mathrm{y-}} =
  \frac{\hat{w}^n_{i,j,k} \cdot c^n_{\mathrm{sat},i,j,k} - c^n_{i,j,k}}{T}

in which:

.. math::
  :label: apx-upwind1
  
  \begin{array}{rclcrcl}
    u^n_{z,\mathrm{x+}} &=& \max( 0, u^n_{z,\mathrm{x}} ) &;& u^n_{z,\mathrm{y+}} &=& \max( 0, u^n_{z,\mathrm{y}} ) \\
    u^n_{z,\mathrm{x-}} &=& \min( 0, u^n_{z,\mathrm{x}} ) &;& u^n_{z,\mathrm{y-}} &=& \min( 0, u^n_{z,\mathrm{y}} ) \\
  \end{array}

and 

.. math::
  :label: apx-upwind2
  
  \begin{array}{rclcrcl}
    c^n_{i,j,k,\mathrm{x+}} &=& \frac{c^n_{i+1,j,k} - c^n_{i,j,k}}{\Delta x} &;&
        c^n_{i,j,k,\mathrm{y+}} &=& \frac{c^n_{i,j+1,k} - c^n_{i,j,k}}{\Delta y} \\
    c^n_{i,j,k,\mathrm{x-}} &=& \frac{c^n_{i,j,k} - c^n_{i-1,j,k}}{\Delta x} &;&
        c^n_{i,j,k,\mathrm{y-}} &=& \frac{c^n_{i,j,k} - c^n_{i,j-1,k}}{\Delta y} \\
  \end{array}

Equation :eq:`apx-explicit-generalized` is explicit in
time and adheres to the Courant-Friedrich-Lewis (CFL) condition for
numerical stability. Alternatively, the advection equation can be
discretized implicitly in time for unconditional stability:

.. math::
  :label: apx-implicit-generalized
   
  \frac{c^{n+1}_{i,j,k} - c^n_{i,j,k}}{\Delta t^n} + 
  u^{n+1}_{z,\mathrm{x+}} c^{n+1}_{i,j,k,\mathrm{x+}} + 
  u^{n+1}_{z,\mathrm{y+}} c^{n+1}_{i,j,k,\mathrm{y+}} \\ + 
  u^{n+1}_{z,\mathrm{x-}} c^{n+1}_{i,j,k,\mathrm{x-}} + 
  u^{n+1}_{z,\mathrm{y-}} c^{n+1}_{i,j,k,\mathrm{y-}} =
  \frac{\hat{w}^{n+1}_{i,j,k} \cdot c^{n+1}_{\mathrm{sat},i,j,k} - c^{n+1}_{i,j,k}}{T}

Equation :eq:`apx-explicit-generalized` and
:eq:apx-implicit-generalized` can be rewritten as:

.. math::
  :label: apx-explicit-rewritten
   
  c^{n+1}_{i,j,k} = c^n_{i,j,k} - \Delta t^n \left[ 
  u^n_{z,\mathrm{x+}} c^n_{i,j,k,\mathrm{x+}} + 
  u^n_{z,\mathrm{y+}} c^n_{i,j,k,\mathrm{y+}} \phantom{\frac{c^n_{i,j,k}}{T}} \right. \\ + \left.
  u^n_{z,\mathrm{x-}} c^n_{i,j,k,\mathrm{x-}} + 
  u^n_{z,\mathrm{y-}} c^n_{i,j,k,\mathrm{y-}} +
  \frac{\hat{w}^n_{i,j,k} \cdot c^n_{\mathrm{sat},i,j,k} - c^n_{i,j,k}}{T} \right]

and

.. math::
  :label: apx-implicit-rewritten
   
  c^{n+1}_{i,j,k} + \Delta t^n \left[ 
  u^{n+1}_{z,\mathrm{x+}} c^{n+1}_{i,j,k,\mathrm{x+}} + 
  u^{n+1}_{z,\mathrm{y+}} c^{n+1}_{i,j,k,\mathrm{y+}} \phantom{\frac{c^{n+1}_{i,j,k}}{T}} \right. \\ + \left.
  u^{n+1}_{z,\mathrm{x-}} c^{n+1}_{i,j,k,\mathrm{x-}} + 
  u^{n+1}_{z,\mathrm{y-}} c^{n+1}_{i,j,k,\mathrm{y-}} +
  \frac{\hat{w}^{n+1}_{i,j,k} \cdot c^{n+1}_{\mathrm{sat},i,j,k} - c^{n+1}_{i,j,k}}{T} \right] = c^n_{i,j,k}

and combined using a weighted average:

.. math::
  :label: apx-combined
   
  c^{n+1}_{i,j,k} + \Gamma \Delta t^n \left[ 
  u^{n+1}_{z,\mathrm{x+}} c^{n+1}_{i,j,k,\mathrm{x+}} + 
  u^{n+1}_{z,\mathrm{y+}} c^{n+1}_{i,j,k,\mathrm{y+}} \phantom{\frac{c^{n+1}_{i,j,k}}{T}} \right. \\ + \left.
  u^{n+1}_{z,\mathrm{x-}} c^{n+1}_{i,j,k,\mathrm{x-}} + 
  u^{n+1}_{z,\mathrm{y-}} c^{n+1}_{i,j,k,\mathrm{y-}} +
  \frac{\hat{w}^{n+1}_{i,j,k} \cdot c^{n+1}_{\mathrm{sat},i,j,k} - c^{n+1}_{i,j,k}}{T} \right] \\ =
  c^n_{i,j,k} - (1 - \Gamma) \Delta t^n \left[ 
  u^n_{z,\mathrm{x+}} c^n_{i,j,k,\mathrm{x+}} + 
  u^n_{z,\mathrm{y+}} c^n_{i,j,k,\mathrm{y+}} \phantom{\frac{c^n_{i,j,k}}{T}} \right. \\ + \left.
  u^n_{z,\mathrm{x-}} c^n_{i,j,k,\mathrm{x-}} + 
  u^n_{z,\mathrm{y-}} c^n_{i,j,k,\mathrm{y-}} +
  \frac{\hat{w}^n_{i,j,k} \cdot c^n_{\mathrm{sat},i,j,k} - c^n_{i,j,k}}{T} \right]

in which :math:`\Gamma` is a weight that ranges from 0 -- 1 and
determines the implicitness of the scheme. The scheme is implicit with
:math:`\Gamma = 0`, explicit with :math:`\Gamma = 1` and semi-implicit
otherwise. :math:`\Gamma = 0.5` results in the semi-implicit Crank-Nicolson
scheme.

Equation :eq:`apx-upwind2` is back-substituted in Equation
:eq:`apx-combined`:

.. math::
  :label: apx-combined-substituted
   
  c^{n+1}_{i,j,k} + \Gamma \Delta t^n \left[ 
  u^{n+1}_{z,\mathrm{x+}} \frac{c^{n+1}_{i+1,j,k} - c^{n+1}_{i,j,k}}{\Delta x} + 
  u^{n+1}_{z,\mathrm{y+}} \frac{c^{n+1}_{i,j+1,k} - c^{n+1}_{i,j,k}}{\Delta y} \right. \\ + \left.
  u^{n+1}_{z,\mathrm{x-}} \frac{c^{n+1}_{i,j,k} - c^{n+1}_{i-1,j,k}}{\Delta x} + 
  u^{n+1}_{z,\mathrm{y-}} \frac{c^{n+1}_{i,j,k} - c^{n+1}_{i,j-1,k}}{\Delta y} +
  \frac{\hat{w}^{n+1}_{i,j,k} \cdot c^{n+1}_{\mathrm{sat},i,j,k} - c^{n+1}_{i,j,k}}{T} \right] \\ =
  c^n_{i,j,k} - (1 - \Gamma) \Delta t^n \left[ 
  u^n_{z,\mathrm{x+}} \frac{c^n_{i+1,j,k} - c^n_{i,j,k}}{\Delta x} + 
  u^n_{z,\mathrm{y+}} \frac{c^n_{i,j+1,k} - c^n_{i,j,k}}{\Delta y} \right. \\ + \left.
  u^n_{z,\mathrm{x-}} \frac{c^n_{i,j,k} - c^n_{i-1,j,k}}{\Delta x} + 
  u^n_{z,\mathrm{y-}} \frac{c^n_{i,j,k} - c^n_{i,j-1,k}}{\Delta y} +
  \frac{\hat{w}^n_{i,j,k} \cdot c^n_{\mathrm{sat},i,j,k} - c^n_{i,j,k}}{T} \right]

and rewritten:

.. math::
  :label: apx-combined-rewritten
   
  \left[ 1 - \Gamma \left( 
      u^{n+1}_{z,\mathrm{x+}} \frac{\Delta t^n}{\Delta x} + 
      u^{n+1}_{z,\mathrm{y+}} \frac{\Delta t^n}{\Delta y} - 
      u^{n+1}_{z,\mathrm{x-}} \frac{\Delta t^n}{\Delta x} - 
      u^{n+1}_{z,\mathrm{y-}} \frac{\Delta t^n}{\Delta y} +
      \frac{\Delta t^n}{T}
    \right)
  \right] c^{n+1}_{i,j,k} \\ +
  \Gamma \left(
    u^{n+1}_{z,\mathrm{x+}} \frac{\Delta t^n}{\Delta x} c^{n+1}_{i+1,j,k} + 
    u^{n+1}_{z,\mathrm{y+}} \frac{\Delta t^n}{\Delta y} c^{n+1}_{i,j+1,k} - %\right. \\ - \left.
    u^{n+1}_{z,\mathrm{x-}} \frac{\Delta t^n}{\Delta x} c^{n+1}_{i-1,j,k} - 
    u^{n+1}_{z,\mathrm{y-}} \frac{\Delta t^n}{\Delta y} c^{n+1}_{i,j-1,k}
  \right) \\ =
  \left[ 1 + (1 - \Gamma) \left( 
      u^n_{z,\mathrm{x+}} \frac{\Delta t^n}{\Delta x} + 
      u^n_{z,\mathrm{y+}} \frac{\Delta t^n}{\Delta y} - 
      u^n_{z,\mathrm{x-}} \frac{\Delta t^n}{\Delta x} - 
      u^n_{z,\mathrm{y-}} \frac{\Delta t^n}{\Delta y} +
      \frac{\Delta t^n}{T}
    \right)
  \right] c^n_{i,j,k} \\ +
  (1 - \Gamma) \left(
    u^n_{z,\mathrm{x+}} \frac{\Delta t^n}{\Delta x} c^n_{i+1,j,k} + 
    u^n_{z,\mathrm{y+}} \frac{\Delta t^n}{\Delta y} c^n_{i,j+1,k} - %\right. \\ - \left.
    u^n_{z,\mathrm{x-}} \frac{\Delta t^n}{\Delta x} c^n_{i-1,j,k} - 
    u^n_{z,\mathrm{y-}} \frac{\Delta t^n}{\Delta y} c^n_{i,j-1,k}
  \right) \\ - 
  \Gamma \hat{w}^{n+1}_{i,j,k} \cdot c^{n+1}_{\mathrm{sat},i,j,k} \frac{\Delta t^n}{T} -
  (1 - \Gamma) \hat{w}^n_{i,j,k} \cdot c^n_{\mathrm{sat},i,j,k} \frac{\Delta t^n}{T}

and simplified:

.. math::
  :label: apx-combined-simplified
   
  a^{0,0}_{i,j} c^{n+1}_{i,j,k} +
  a^{1,0}_{i,j} c^{n+1}_{i+1,j,k} + 
  a^{0,1}_{i,j} c^{n+1}_{i,j+1,k} -
  a^{-1,0}_{i,j} c^{n+1}_{i-1,j,k} - 
  a^{0,-1}_{i,j} c^{n+1}_{i,j-1,k} = y_{i,j,k}

where the implicit coefficients are defined as:

.. math::
  :label: apx-implicitcoef
  
  \begin{array}{rclcrcl}
    a^{0,0}_{i,j} &=& \left[1 - \Gamma \left( 
      u^{n+1}_{z,\mathrm{x+}} \frac{\Delta t^n}{\Delta x} + 
      u^{n+1}_{z,\mathrm{y+}} \frac{\Delta t^n}{\Delta y} - 
      u^{n+1}_{z,\mathrm{x-}} \frac{\Delta t^n}{\Delta x} - 
      u^{n+1}_{z,\mathrm{y-}} \frac{\Delta t^n}{\Delta y} +
      \frac{\Delta t^n}{T}
    \right) \right] \\
    a^{1,0}_{i,j} &=& \Gamma u^{n+1}_{z,\mathrm{x+}} \frac{\Delta t^n}{\Delta x} \\
    a^{0,1}_{i,j} &=& \Gamma u^{n+1}_{z,\mathrm{y+}} \frac{\Delta t^n}{\Delta y} \\
    a^{-1,0}_{i,j} &=& \Gamma u^{n+1}_{z,\mathrm{x-}} \frac{\Delta t^n}{\Delta x} \\
    a^{0,-1}_{i,j} &=& \Gamma u^{n+1}_{z,\mathrm{y-}} \frac{\Delta t^n}{\Delta y} \\
  \end{array}

and the explicit right-hand side as:

.. math::
  :label: apx-explicitrhs
   
  y^n_{i,j,k} = \left[ 1 + (1 - \Gamma) \left( 
      u^n_{z,\mathrm{x+}} \frac{\Delta t^n}{\Delta x} + 
      u^n_{z,\mathrm{y+}} \frac{\Delta t^n}{\Delta y} - 
      u^n_{z,\mathrm{x-}} \frac{\Delta t^n}{\Delta x} - 
      u^n_{z,\mathrm{y-}} \frac{\Delta t^n}{\Delta y} +
      \frac{\Delta t^n}{T}
    \right)
  \right] c^n_{i,j,k} \\ +
  (1 - \Gamma) \left(
    u^n_{z,\mathrm{x+}} \frac{\Delta t^n}{\Delta x} c^n_{i+1,j,k} + 
    u^n_{z,\mathrm{y+}} \frac{\Delta t^n}{\Delta y} c^n_{i,j+1,k} -
    u^n_{z,\mathrm{x-}} \frac{\Delta t^n}{\Delta x} c^n_{i-1,j,k} - 
    u^n_{z,\mathrm{y-}} \frac{\Delta t^n}{\Delta y} c^n_{i,j-1,k}
  \right) \\ - 
  \Gamma \hat{w}^{n+1}_{i,j,k} \cdot c^{n+1}_{\mathrm{sat},i,j,k} \frac{\Delta t^n}{T} -
  (1 - \Gamma) \hat{w}^n_{i,j,k} \cdot c^n_{\mathrm{sat},i,j,k} \frac{\Delta t^n}{T}

The offshore boundary is defined to be zero-flux, the
onshore boundary has a constant transport gradient and the lateral
boundaries are circular:

.. math::
  :label: apx-boundaryconditions
  
  \begin{array}{rclcrcl}
    c^{n+1}_{1,j,k} &=& 0 \\
    c^{n+1}_{n_{\mathrm{x}}+1,j,k} &=& 2 c^{n+1}_{n_{\mathrm{x}},j,k} - c^{n+1}_{n_{\mathrm{x}}-1,j,k} \\
    c^{n+1}_{i,1,k} &=& c^{n+1}_{i,n_{\mathrm{y}}+1,k} \\
    c^{n+1}_{i,n_{\mathrm{y}}+1,k} &=& c^{n+1}_{i,1,k} \\
  \end{array}

Shear stress perturbation for non-perpendicular wind directions
---------------------------------------------------------------

The shear stress perturbation 𝛿𝜏 is estimated following the analytical description of the influence of alow and smooth hill in the wind profile by Weng et al. (1991). The perturbation is given by the Fouriertransformed components of the shear stress perturbation in the unperturbed wind direction which are the functions 𝛿𝜏𝑥(𝑘) and 𝛿𝜏𝑦(𝑘). The x-direction is defined by the direction of the wind velocity 𝑣0 on a flat bed, while the y direction is then the transverse.

As a result, the perturbation theory can only estimate the shear stress induced by the morphology-wind interaction in parallel direction of wind. Therefore, model simulations were, up to now, limited to input wind directions parallel to the cross­shore axis of the grid.

To overcome this limitation and to allow for modelling directional winds, an overlaying computational grid is introduced in AeoLiS, which rotates with the changing wind direction per time step. By doing this, the shear stresses are always estimated in the positive x-direction of the computational grid. The following steps are executed for each time step:

1. Create a computational grid alligned with the wind direction (set_computational_grid)
2. Add and fill buffer around the original grid
3. Populate computation grid by rotating it to the current wind direction and interpolate the original topography on it. Additionally, edges around 
4. Compute the morphology-wind induced shear stress by using the perturbation theory
5. Add the only wind induced wind shear stresses to the computational grid
6. Rotate both the grids and the total shear stress results in opposite direction
7. Interpolate the total shear stress results from the computational grid to the original grid
8. Rotate the wind shear stress results and the original grid back to the original orientation

.. note:: 
   The extra rotations in the last two steps are necessary as a simplified, but faster in terms of computational time, interpolation method is used.


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

Vegetation (from :cite:`Strypsteen2024`)
-------------------------------------

In AeoLiS, the influence of vegetation on dune evolution is comprehensively 
addressed. This includes modelling the intrinsic growth of vegetation, 
accounting for factors such as growth and decay due to burial :cite:`DuranMoore2013`, 
lateral expansion and establishment :cite:`Keijsers20106`, as well as simulating 
the destruction of vegetation caused by hydrodynamic processes. In the event of 
cell inundation, vegetation density is reduced as a result. Inspired by the 
Coastal Dune Model (CDM) proposed by :cite:`DuranMoore2013`, AeoLiS incorporates 
vegetation-wind interaction using the expression established by :cite:`DuranHerrmann2006`:

.. math::
   :label: shear_reduction_vegetation

   \frac{u_{\text{veg}}}{u_*} = \frac{1}{\sqrt{1 + \Gamma \rho_{\text{veg}}}}

where the ratio of shear velocity in the presence of vegetation (:math:`u_{*,\text{veg}}`) to the unobstructed 
shear velocity (u∗) is determined by a vegetation-related roughness parameter (Γ) and the 
vegetation density within a unit area of the grid cell (\rho_{\text{veg}}). In the model, Γ = 16 is derived 
from plant form drag and geometry values documented for creosote communities :cite:`DuranHerrmann2006`. 
This implementation calculates the expression on each model grid cell, with higher vegetation density 
(expressed by \rho_{\text{veg}}) leading to a more substantial reduction in shear velocity compared to sparse 
vegetation. By integrating these physical and ecological processes, AeoLiS simulates spatial patterns and temporal 
variations in sediment transport and morphological changes resulting from aeolian processes in coastal 
environments.

\rho_{\text{veg}} can vary in space and time and is determined by the ratio of the actual vegetation height (hveg) 
to the maximum vegetation height (hveg,max), and can vary between 0 and 1 (:cite:`DuranHerrmann2006`):

:math:`\rho_{\text{veg}} = \left( \frac{h_{\text{veg}}}{h_{\text{veg,max}}} \right)^2`

This assumption is based on the idea that burying vegetation reduces its height, which indicates a decrease in 
actual cover. The change in vegetation density per grid cell is directly linked to the alteration in vegetation 
height within that specific cell. This height variation is influenced by both the growth rate of the vegetation and 
the rate of sediment burial. If the vegetation density remains constant over time, it suggests either no 
sedimentation or a growth rate equal to the rate of sediment burial within the cell. Vegetation growth and decay 
follow the model proposed by :cite:`DuranHerrmann2006`, modified to include :math:`\delta z_{\text{b,opt}}`
(m/year), representing sediment burial for optimal growth that shifts the peak of optimal growth:

:math:`\frac{\delta h_{\text{veg}}}{\delta t} = V_{\text{ver}} \left(1 - \frac{h_{\text{veg}}}{h_{\text{veg,max}}}\right) - \gamma_{\text{veg}} \left| \frac{\delta z_{\text{b,veg}}}{\delta t} - \delta z_{\text{b,opt}} \right|`

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

Vegetation can begin to grow through lateral propagation or random germination. Once established, it can continue 
to grow and spread laterally. The uncertainties associated with random germination are handled on a cell-by-cell 
basis using a probabilistic approach, similar to the cellular automata method described by Keijsers et al. (2016). 
AeoLiS incorporates a germination probability, denoted as ρ_{ger}, for each grid cell. This probability is constant 
across the domain, except in eroding grid cells (where bed elevation decreases), where ρ_{ger} is set to 0. Lateral 
propagation is determined by identifying the boundaries between vegetated and non-vegetated cells, with the 
parameter ρ_{lat} adjusting the likelihood of lateral propagation at these boundaries.

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

Numerical solution of the Boussinesq groundwater equation
---------------------------------------------------------
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

