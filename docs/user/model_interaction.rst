Model interaction
=================

ADD INTRODUCTION TO MODEL INTRODUCTION - WHY, HOW

Callback function
-----------------

As a convenience functionality the current implementation
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



Hotstart
--------


Basic Model Interface (BMI)
---------------------------

A Basic Model Interface (BMI, :cite:`Peckham2013`) is implemented
that allows interaction with the model during run time. The model can
be implemented as a library within a larger framework as the interface
exposes the initialization, finalization and time stepping
routines. 


