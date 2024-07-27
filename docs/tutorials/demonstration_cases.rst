Demonstration cases
===================
The `example folder <https://github.com/openearth/aeolis-python/aeolis/examples>`_ contains a collection of input files needed to run example models. The different examples are explained in the `readme-file <https://github.com/openearth/aeolis-python/aeolis/examples/readme.txt>`. They are intended to serve as an example of the model's capabilities and instruction on how to set up and organize models. The model input files can be modified using the pre-processing [ADD LINK] tools, and the model output can be visualized with the post-processing tools [ADD LINK].

To run the examples:

1. Install aeolis

.. code-block:: bash

pip install aeolis

2. Install the examples

To install the examples do the following:

.. code-block:: bash

   aeolis examples .

This will install the examples in the current directory. If you want to install them in a different directory, replace the "." with the path to the directory you want to install the examples in.

3. Run the simulation

Run the simulation with (example parabolic dune model):

.. code-block:: bash

    aeolis run aeolis_examples/parabolic_dune_model/aeolis.txt


4. View/plot the results (example parabolic dune model):

.. code-block:: bash

    [ADD DESCRIPTION OF HOW TO PLOT NETCDF]
