Installation
============

Requirements
------------

Python packages
^^^^^^^^^^^^^^^

* bmi-python: http://github.com/openearth/bmi-python
* numpy
* scipy
* netCDF4
* docopt

External libraries (Windows)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These libraries are needed on Windows if the Python package netCDF4 is installed manually.

* Microsoft Visual C++ Compiler for Python 2.7: http://aka.ms/vcpython27
* msinttypes for stdint.h: https://code.google.com/archive/p/msinttypes/
* HDF5 headers: https://www.hdfgroup.org/HDF5/release/obtain5.html
* netCDF4 headers: https://github.com/Unidata/netcdf-c/releases

* Set environment variables HDF5_DIR and NETCDF_DIR to the respective installation paths