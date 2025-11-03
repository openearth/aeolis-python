.. _installation:

Installation
============

AeoLiS is a Python package that can be installed from PyPI or from source. For the source code, see the `GitHub repository <https://github.com/openearth/aeolis-python>`_.

Requirements
------------

- Python 3.9 or newer 
- pip 22.0 or newer
- netCDF4

Dependencies
""""""""""""

- docopt
- typer
- bmi-python
- netCDF4
- scipy
- numpy
- matplotlib
- numba

Installing from PyPI
---------------------

On the comand line of your working environment (Bash/Shell, Conda, Mamba, or similar), run the following: 

.. code:: shell

   pip install aeolis

.. attention:: 

   For Windows users, the recommend way to install AeoLiS is to use `Anaconda <https://docs.anaconda.com/free/anaconda/install/windows/>`_.


Installing from source
-----------------------


1. Clone the repository using Git, or download the source code.

2. Go to the `aeolis-python` directory and install using pip
   
   .. code:: shell

    cd aeolis-python/
    pip install .
   

Running AeoLiS
----------------

Example from command line:


.. code:: shell

   aeolis run params.txt

.. note::

   Model parameters and other configuration is passed in a `params.txt`. See the :ref:`default settings` for more details.  
