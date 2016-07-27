Source code documentation
=========================

Model classes
-------------

The AeoLiS model is based on two main model classes:
:class:`~model.AeoLiS` and
:class:`~model.AeoLiSRunner`. The former is the actual,
low-level, BMI-compatible class that implements the basic model
functions and numerical schemes. The latter is a convenience class
that implements a time loop, netCDF4 output, a progress indicator and
a callback function that allows the used to interact with the model
during runtime.

The additional :class:`~model.WindGenerator` class to generate
realistic wind time series is available from the same module.

AeoLiS
^^^^^^

.. autoclass:: model.AeoLiS
               :special-members:
               :members:

AeoLiSRunner
^^^^^^^^^^^^^

.. autoclass:: model.AeoLiSRunner
               :special-members:
               :members:

WindGenerator
^^^^^^^^^^^^^

.. autoclass:: model.WindGenerator
               :special-members:
               :members:

Physics modules
---------------

Bathymetry and bed composition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: bed
                :members:

Wind velocity and direction
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: wind
                :members:

.. automodule:: shear
                :members:

Wind velocity threshold
^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: threshold
                :members:

Tides, meteorology and soil moisture content
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: hydro
                :members:

Sediment transport
^^^^^^^^^^^^^^^^^^

.. automodule:: transport
                :members:

Helper modules
--------------

Input/Output
^^^^^^^^^^^^

.. automodule:: io
                :members:

netCDF4 output
^^^^^^^^^^^^^^

.. automodule:: netcdf
                :members:

Logging
^^^^^^^

.. automodule:: log
                :members:

Plotting
^^^^^^^^

.. automodule:: plot
                :members:

Command-line tools
^^^^^^^^^^^^^^^^^^

.. automodule:: console
                :members:

Miscellaneous
^^^^^^^^^^^^^

.. automodule:: utils
                :members:
                   
