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
routines. Because of that, we can use this set of functions to control 
and interact with AeoLiS from a Python environment.


.. code::

	import numpy as np
	
	# Timing settigns
	start_time = 0.
	end_time = 157680000.
	dt = 3600. # 1 day (output frequency of dfm output)
	
	# Find file in the same directory as this script that contains aeolis.txt
	configfile = 'aeolis.txt' 
	os.chdir(os.path.dirname(configfile))
	
	# Create AeoLiS BMI Wrapper
	aeolis_wrapper = AeoLiSRunner(configfile)
	
	# Initialize the wrapper
	aeolis_wrapper.initialize()
	
	# Loop over all timesteps to run the model
	for t in np.arange(start_time, end_time, dt):
	
	    # Update AeoLiS
	    aeolis_wrapper.update(dt)
	    aeolis_wrapper.output_write()
	    
		# Make modifications to the variables in the model via data or another model, e.g.:
	    # x_aeolis = aeolis_wrapper.get_var('x')
	    # y_aeolis = aeolis_wrapper.get_var('y’)
	
	    # aeolis_wrapper.set_var('zb', zb_aeolis)
	
	# Finalize the wrapper
	aeolis_wrapper.finalize()


An overview of recent BMI-AeoLiS applications:
   - van Westen, B., Luijendijk, A. P., de Vries, S., Cohn, N., Leijnse, T. W., & de Schipper, M. A. (2024). Predicting marine and aeolian contributions to the Sand Engine’s evolution using coupled modelling. Coastal Engineering, 188, 104444. (https://www.sciencedirect.com/science/article/pii/S0378383923001680)
   - van Westen, B., Leijnse, T., de Schipper, M., Cohn, N., & Luijendijk, A. (2023). Integrated modelling of coastal landforms. In Coastal Sediments 2023: The Proceedings of the Coastal Sediments 2023 (pp. 760-771). (https://pure.tudelft.nl/ws/portalfiles/portal/152631531/Integrated_modelling_Bvw_CS23.pdf)






