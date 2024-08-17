Appendix 
=========

Numerical details of solving the advection equation
-----------------------------------------------------


.. _pieter_num:

Conservative Euler Backward Implicit (pieter solver)
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


.. _trunk_num:

Euler Schemes in non-conservative form (trunk solver)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

