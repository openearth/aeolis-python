.. AeoLiS documentation master file, created by
   sphinx-quickstart on Tue Dec 29 12:40:13 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to AeoLiS's Documentation
==================================

.. raw:: html

   <div style="overflow: hidden;">
      <iframe src="https://www.youtube.com/embed/dHr2NlGkSE4?&autoplay=1"  width="100%" height="315" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; web-share"></iframe>
   </div>

-------

AeoLiS is a process-based model for simulating aeolian sediment transport in situations where supply-limiting factors are important,
like in coastal environments. Supply-limitations currently supported are soil moisture contents, sediment sorting and armouring, bed slope effects, air humidity and roughness elements.

This documentation describes the Python implementation of the AeoLiS
model. The source code of the Python implementation can be found in `GitHub <https://github.com/openearth/aeolis-python>`_.

Contents
--------

.. toctree::
   :caption: User Documentation
   :maxdepth: 2

   user/installation
   user/whatsnew
   user/model
   user/implementation
   user/defaults
   user/inputfiles
   user/output
   user/sourcecode
   user/bibliography

.. toctree::
   :caption: Tutorials
   :maxdepth: 2

   tutorials/sandmotor
   tutorials/2D-parabolic

.. toctree::
   :caption: Developer Documentation
   :maxdepth: 2

   developer/quickstart
   developer/testing-introduction
   developer/unit-testing

.. toctree::
   :caption: Current Development
   :maxdepth: 2

   developer/modularity
   developer/domain-decomposition


Acknowledgements
================

* AeoLiS was initially developed at Delft University of Technology with support from the ERC-Advanced Grant 291206 Nearshore Monitoring and Modeling (`NEMO <http://nemo.citg.tudelft.nl>`_) and `Deltares <http://www.deltares.nl>`_.
* AeoLiS is currently maintained by `Bart van Westen <Bart.vanWesten@deltares.nl>`_ at Deltares, `Nick Cohn <nick.cohn@usace.army.mil>`_ at U.S. Army Engineer Research and Development Center (ERDC) and `Sierd de Vries <Sierd.deVries@tudelft.nl>`_ at Delft University of Technology.
* Further developement of AeoLiS is supported by the `Digital Competence Centre <https://dcc.tudelft.nl>`_, Delft University of Technology.


Indices and tables
==================

* :ref:`genindex`
* :ref:`search`

