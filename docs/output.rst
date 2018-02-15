Model state/output
==================

The AeoLiS model state is described by a collection of spatial grid
variables with at least one value per horizontal grid cell.  Specific
model state variables can also be subdivided over bed composition
layers and/or grain size fractions.  All model state variables can be
part of the model netCDF4 output. The current model state variables
are listed below.

.. literalinclude:: ../aeolis/constants.py
   :language: python
   :start-after: #: Aeolis model state variables
   :end-before: #: AeoLiS model default configuration
