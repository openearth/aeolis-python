Source code 
=========================

Here you can find the documentation with direct links to the actual AeoLiS code. You can click on the green [source] button next to the classes and modules below to access the specific source code. You can use ctr-f to look for a specific functionality or variable. It still may be a bit difficult to browse through, in addition you can find an overview of all module code `here <https://aeolis.readthedocs.io/en/latest/_modules/index.html>`_ 


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

Avalanching
^^^^^^^^^^^

.. automodule:: avalanching
                :members:

Vegetation
^^^^^^^^^^^

.. automodule:: vegetation
                :members:

Fences
^^^^^^

.. automodule:: fences
                :members:

Marine Erosion
^^^^^^^^^^^^^^

.. automodule:: erosion
                :members:

Helper modules
--------------

Input/Output
^^^^^^^^^^^^

.. automodule:: inout
                :members:

netCDF4 output
^^^^^^^^^^^^^^

.. automodule:: netcdf
                :members:

Plotting
^^^^^^^^

.. automodule:: plot
                :members:

Command-line tools
^^^^^^^^^^^^^^^^^^

.. automodule:: console
                :members:

Debugging
^^^^^^^^^^

.. automodule:: console_debug
                :members:

.. automodule:: run_console
                :members:

Miscellaneous
^^^^^^^^^^^^^

.. automodule:: utils
                :members:

.. automodule:: gridparams
                :members:


.. note::

    Sierd's favorite function is:
    :py:mod:`aeolis.bed.prevent_tiny_negatives`
