What's New
==========

v1.1.3 (unreleased)
-------------------

Breaking changes
^^^^^^^^^^^^^^^^

* Removed support for statistical variable names with dot-notation
  (e.g. `.avg` and `.sum`).

Improvements
^^^^^^^^^^^^

None.

New functions/methods
^^^^^^^^^^^^^^^^^^^^^

None.

Bug fixes
^^^^^^^^^

None.

Tests
-----

None.

v1.1.2 (21 December 2017)
-------------------------

Breaking changes
^^^^^^^^^^^^^^^^

* Changed name of statistics variables that describe the average,
  minimum, maximum, cumulative values, or variance of a model state
  variable. The variables names that used to end with `.avg`, `.sum`,
  etc. now end with `_avg`, `_sum`, etc. The new naming convention was
  already adopted in the netCDF output in order to be compatible with
  the CF-1.6 convention, but is now also adopted in, for example, the
  Basic Model Interface (BMI). Old notation is deprecated but still
  supported.

Improvements
^^^^^^^^^^^^

* Prepared for continuous integration through CircleCI.

New functions/methods
^^^^^^^^^^^^^^^^^^^^^

None.

Bug fixes
^^^^^^^^^

* Use percentages (0-100) rather than fractions (0-1) in the
  formulation of Belly and Johnson that describes the effect of soil
  moisture on the shear velocity threshold. Thanks to Dano Roelvink
  and Susana Costas (b3d992b).

Tests
-----

* Reduced required accuracy for mass conservation tests from
  0.00000000000001% to 1%.

v1.1.1 (15 November 2017)
-------------------------

Breaking changes
^^^^^^^^^^^^^^^^

None.

Improvements
^^^^^^^^^^^^

* Made code compatible with Python 3.x.
* Prepared and uploaded package to PyPI.
* Switch back to original working directory after finishing
  simulation.
* Removed double definition of model state. Not only defined in
  `constants.MODEL_STATE`.
* Also write initial model state to output.
* Made netCDF output compatible with CF-1.6 convention.

New functions/methods
^^^^^^^^^^^^^^^^^^^^^

* Added support to run a default model for testing purposes by setting
  the configuration file as "DEFAULT".
* Added generic framework for reading and applying spatial
  masks. Implemented support for wave, tide and threshold masks
  specifically.
* Added option to include a reference date in netCDF output.
* Added experimental option for constant boundary conditions.
* Added support for reading and writing hotstart files to load a
  (partial) model state upon initialisation.
* Added preliminary wind shear perturbation module. Untested.
* Added support to switch on or off specific processes.
* Added support for immutable model state variables. This
  functionality can be combined with BMI or hotstart files to prevent
  external process results to be overwritten by the model.
* Added option to specify wind direction convention (nautical or
  cartesian).

Bug fixes
^^^^^^^^^

* Fixed conversion from volume to mass using porosity and density (fe9aa52).
* Update water level with bed updates to prevent loss of water due to
  bed level change (fe9aa52).
* Fixed mass bug in base layer that drained sediment from bottom
  layers, resulting in empty layers (f612760).
* Made removal of negative concentrations mass conserving by scraping
  the concentrations from all other grid cells (03de813).

Tests
-----

* Added tests to check mass conservation in bed mixing routines.
* Added integration tests.

v1.1.0 (27 July 2016)
---------------------

Initial release
