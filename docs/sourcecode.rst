Source code documentation
=========================

Model classes
-------------

The AeoLiS model is based on two main model classes:
:class:`~aeolis.model.AeoLiS` and
:class:`~aeolis.model.AeoLiSRunner`. The former is the actual,
low-level, BMI-compatible class that implements the basic model
functions and numerical schemes. The latter is a convenience class
that implements a time loop, netCDF4 output, a progress indicator and
a callback function that allows the used to interact with the model
during runtime.

The additional :class:`~aeolis.model.WindGenerator` class to generate
realistic wind time series is available from the same module.

AeoLiS
^^^^^^

.. autoclass:: aeolis.model.AeoLiS
               :special-members:
               :members:

AeoLiSRunner
^^^^^^^^^^^^^

.. autoclass:: aeolis.model.AeoLiSRunner
               :special-members:
               :members:

WindGenerator
^^^^^^^^^^^^^

.. autoclass:: aeolis.model.WindGenerator
               :special-members:
               :members:

Physics modules
---------------

Bathymetry and bed composition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: aeolis.bed
                :members:

Wind velocity and direction
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: aeolis.wind
                :members:

Wind velocity threshold
^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: aeolis.threshold
                :members:

Tides, meteorology and soil moisture content
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: aeolis.hydro
                :members:

Sediment transport
^^^^^^^^^^^^^^^^^^

.. automodule:: aeolis.transport
                :members:

Helper modules
--------------

Input/Output
^^^^^^^^^^^^

.. automodule:: aeolis.io
                :members:

netCDF4 output
^^^^^^^^^^^^^^

.. automodule:: aeolis.netcdf
                :members:

Logging
^^^^^^^

.. automodule:: aeolis.log
                :members:

Plotting
^^^^^^^^

.. automodule:: aeolis.plot
                :members:

Command-line tools
^^^^^^^^^^^^^^^^^^

.. automodule:: aeolis.cmd
                :members:

Miscellaneous
^^^^^^^^^^^^^

.. automodule:: aeolis.utils
                :members:
                   
