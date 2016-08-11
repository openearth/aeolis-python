Numerical implementation
========================

The numerical implementation of the equations presented in
:ref:`model` is explained here.  The implementation is available as
Python package through the OpenEarth GitHub repository at:
http://www.github.com/openearth/aeolis-python/

Advection equation
------------------

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

These boundary conditions can be combined with Equation
:eq:`apx-combined-simplified`, :eq:`apx-implicitcoef` and
:eq:`apx-explicitrhs` into a linear system of equations:

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
      \vec{c}_1 \\ \vec{c}_2 \\ \vdots \\ \vdots \\ \vec{c}_{n_{\mathrm{y}}} \\ \vec{c}_{n_{\mathrm{y}}+1} \\
    \end{array} 
  \right] = \left[ 
    \begin{array}{c}
      \vec{y}_1 \\ \vec{y}_2 \\ \vdots \\ \vdots \\ \vec{y}_{n_{\mathrm{y}}} \\ \vec{y}_{n_{\mathrm{y}}+1} \\
    \end{array} 
  \right]
    
where each item in the matrix is again a matrix :math:`A^l_j` and
each item in the vectors is again a vector :math:`\vec{c}_j` and
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

for :math:`l \neq 0`. The vectors :math:`\vec{c}_{j,k}` and :math:`\vec{y}_{j,k}`
read:

.. math::

  \begin{array}{rclrcl}
    \vec{c}_{j,k} &=& \left[ 
      \begin{array}{c}
        c^{n+1}_{1,j,k} \\
        c^{n+1}_{2,j,k} \\
        c^{n+1}_{3,j,k} \\
        \vdots \\
        c^{n+1}_{n_{\mathrm{x}}-1,j,k} \\
        c^{n+1}_{n_{\mathrm{x}},j,k} \\
        c^{n+1}_{n_{\mathrm{x}}+1,j,k} \\
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

Implicit solver
---------------

The linear system defined in Equation :eq:`apx-system` is solved by a
sparse matrix solver for each sediment fraction separately in
ascending order of grain size. Initially, the weights
:math:`\hat{w}^{n+1}_{i,j,k}` are chosen according to the grain size
distribution in the bed and the air. The sediment availability
constraint is checked after each solve:

.. math::
     m_{\mathrm{a}} \geq \frac{\hat{w}^{n+1}_{i,j,k} c^{n+1}_{\mathrm{sat},i,j,k} - c^{n+1}_{i,j,k}}{T} \Delta t^n

If the constraint if violated, a new estimate for the weights
is back-calculated following:

.. math::
  \hat{w}^{n+1}_{i,j,k} = \frac{ c^{n+1}_{i,j,k} + m_{\mathrm{a}} \frac{\Delta t^n}{T} }{c^{n+1}_{\mathrm{sat},i,j,k}}

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
  p_{\mathrm{g}} = \frac{p_{\mathrm{V}} \cdot \rho_{\mathrm{w}}}{\rho_{\mathrm{p}} \cdot (1 - p)}

where :math:`\rho_{\mathrm{w}}` [:math:`\mathrm{kg/m^3}`] and
:math:`\rho_{\mathrm{p}}` [:math:`\mathrm{kg/m^3}`] are the water and particle
density respectively and :math:`p` [-] is the porosity. Values for
:math:`p_{\mathrm{g}}` smaller than 0.005 do not affect the shear velocity
threshold (:cite:`Pye1990`). Values larger than 0.064 (or 10\%
volumetric content) cease transport (:cite:`DelgadoFernandez2010`),
which is implemented as an infinite shear velocity threshold.

Exploratory model runs of the unsaturated soil with the HYDRUS1D
(:cite:`Simunek1998`) hydrology model show that the increase of the
volumetric water content to saturation is almost instantaneous with
rising tide. The drying of the beach surface through infiltration
shows an exponential decay. In order to capture this behavior the
volumetric water content is implemented according to:

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
content halves. :math:`E_{\mathrm{v}}` [m/s] is the evaporation rate that is
implemented through an adapted version of the Penman equation
(:cite:`Shuttleworth1993`):

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

Roughness elements
^^^^^^^^^^^^^^^^^^

The shear velocity threshold is updated based on the presence of
roughness elements following :cite:`Raupach1993`:

.. math::
  f_{u_{\mathrm{* th},R}} = \sqrt{(1 - m \cdot \sum_{k=k_0}^{n_k}{\hat{w}_k^{\mathrm{bed}}})
    (1 + \frac{m \beta}{\sigma} \cdot \sum_{k=k_0}^{n_k}{\hat{w}_k^{\mathrm{bed}}})}

by assuming:

.. math::
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
  f_{u_{\mathrm{* th}},S} = 1.03 \cdot \exp(0.1027 \cdot p_{\mathrm{s}})

where :math:`f_{u_{\mathrm{* th},S}}` [-] is a factor in Equation
:eq:`apx-shearvelocity` and :math:`p_{\mathrm{s}}` [-] is the salt
content [mg/g]. Currently, no model is implemented that predicts the
instantaneous salt content. The spatial varying salt content needs to
be specified by the user, for example through the BMI interface.

Basic Model Interface (BMI)
---------------------------

A Basic Model Interface (BMI, :cite:`Peckham2013`) is implemented
that allows interaction with the model during run time. The model can
be implemented as a library within a larger framework as the interface
exposes the initialization, finalization and time stepping
routines. As a convenience functionality the current implementation
supports the specification of a callback function. The callback
function is called at the start of each time step and can be used to
exchange data with the model, e.g. update the topography from
measurements.

An example of a callback function, that is referenced in the model
input file or through the model command-line options as
``callback.py:update``, is:

.. code::

   import numpy as np

   def update(model):
     val = model.get_var('zb')
     val_new = val.copy()
     val_new[:,:] = np.loadtxt('measured_topography.txt')
     model.set_var('zb', val_new)
