.. AeoLiS documentation master file, created by
   sphinx-quickstart on Tue Dec 29 12:40:13 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to AeoLiS's documentation!
==================================

AeoLiS is a process-based model for simulating aeolian sediment
transport in situations where supply-limiting factors are important,
like in coastal environments. Supply-limitations currently supported
are soil moisture contents, sediment sorting and armouring, bed slope
effects, air humidity and roughness elements.

This documentation described the Python implementation of the AeoLiS
model. Originally, the model was developed in Fortran and the Python
version only functioned as sandbox for new feature
development. However, the latest implementation of the numerical
solver in Python is as fast as the Fortran implementation. Therefore,
porting the latest developments to Fortran is no priority anymore and
the two versions differ.

AeoLiS is developed and maintained by `Bas Hoonhout
<b.m.hoonhout@tudelft.nl>`_. The source code of the Python
implementation can be found at
`<https://github.com/openearth/aeolis-python>`_. The Fortran
implementation can be found `<https://github.com/openearth/aeolis>`_.

Contents
--------

.. toctree::
   :maxdepth: 3

   sourcecode
   defaults
   

Acknowledgements
================

AeoLiS is developed at `Delft University of Technology
<http://www.tudelft.nl>`_ with support from the ERC-Advanced Grant
291206 Nearshore Monitoring and Modeling (`NEMO
<http://nemo.citg.tudelft.nl>`_) and `Deltares
<http://www.deltares.nl>`_.
   
Indices and tables
==================

* :ref:`genindex`
* :ref:`search`

