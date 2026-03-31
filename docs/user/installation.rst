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
------------

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

2. AeoLiS users may install the package with only the required dependencies. Go to the `aeolis-python` directory and install using pip
   
   .. code:: shell

    cd aeolis-python/
    pip install .
   
3. AeoLiS users who intend to modify the sourcecode can install additional dependencies for test and documentation as follows. Go to root directory `aeolis-python/` and:
   
   .. code:: shell
   
      pip install -e .[dev]


Running AeoLiS
----------------

Example from command line:


.. code:: shell

   aeolis run params.txt

.. note::

   Model parameters and other configuration is passed in a `params.txt`. See the :ref:`default-settings` for more details.  




Installation with uv + Visual Studio Code
===========================================================

AeoLiS can also be installed using a combination of uv and VS Code:

- **uv** instead of pip/conda — A fast Python package manager (10–100x faster than pip) that handles packages, environments, and Python versions in one tool
- **Visual Studio Code (VS Code)** — All-in-one editor for code, terminal, and Git operations
- **Integrated Git** — Manage version control directly within VS Code (no GitHub Desktop)

With this method, everything is in VS Code (code, terminal, Git)


Step 1: Install uv
------------------

**Note:** This is the ONLY step that uses standalone PowerShell. All subsequent steps use the VS Code terminal.

1.1 Open PowerShell
~~~~~~~~~~~~~~~~~~~

Press ``Windows Key``, type **PowerShell**, and click **Windows PowerShell**.

1.2 Install uv
~~~~~~~~~~~~~~

Run this command:

.. code-block:: powershell

   powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"

1.3 Verify Installation
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: powershell

   uv --version

You should see a version number (e.g., ``uv 0.5.x``).

1.4 Close PowerShell
~~~~~~~~~~~~~~~~~~~~

Close the PowerShell window. You won't need it again — everything else will be done in VS Code.

----

Step 2: Install Visual Studio Code
----------------------------------

1. Download from https://code.visualstudio.com/
2. Run the installer with default settings
3. Launch VS Code

----

Step 3: Install VS Code Extensions
----------------------------------

3.1 Open Extensions Panel
~~~~~~~~~~~~~~~~~~~~~~~~~

- Click the Extensions icon in the left sidebar (four squares icon)
- Or press ``Ctrl+Shift+X``

3.2 Install Python Extension
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Search for **Python**
2. Install the extension by **Microsoft**
3. Wait for installation to complete

----

Step 4: Install Git
-------------------

1. Download Git from https://git-scm.com/
2. Run the installer with **default settings**
3. Restart VS Code after installation

----

Step 5: Clone AeoLiS Repository Using VS Code
---------------------------------------------

**Note:** We'll use VS Code's graphical interface to clone the repository — no terminal commands needed for this step.

5.1 Open Command Palette
~~~~~~~~~~~~~~~~~~~~~~~~

- Click **View** → **Command Palette**
- Or press ``Ctrl+Shift+P``

5.2 Start Git Clone
~~~~~~~~~~~~~~~~~~~

1. Type: **Git: Clone**
2. Press **Enter**

5.3 Enter Repository URL
~~~~~~~~~~~~~~~~~~~~~~~~

Paste:

::

   https://github.com/openearth/aeolis-python.git

Press **Enter**.

5.4 Choose Location
~~~~~~~~~~~~~~~~~~~

1. Select a folder (e.g., ``C:\Users\YourName\Github\``)
2. Click **Select as Repository Destination**

5.5 Open the Project
~~~~~~~~~~~~~~~~~~~~

Click **Open** when prompted.

----

Step 6: Open VS Code Terminal
-----------------------------

All commands from here use the VS Code terminal.

- Click **View** → **Terminal**
- Or press ``Ctrl+` ``

----

Step 7: Create Python Environment with uv
-----------------------------------------

7.1 Create Virtual Environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   uv venv --python 3.12 .venv_aeolis

7.2 Activate Environment
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: powershell

   .venv_aeolis\Scripts\activate

If you get a "scripts disabled" error:

.. code-block:: powershell

   Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser

Then retry activation.

----

Step 8: Install AeoLiS
----------------------

8.1 Install in Editable Mode
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   uv pip install -e .

**Editable mode (``-e``):**

- Changes take effect immediately
- No reinstall needed

8.2 Verify Installation
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   aeolis --help

----

Step 9: Configure VS Code Python Interpreter
--------------------------------------------

1. Press ``Ctrl+Shift+P``
2. Select **Python: Select Interpreter**
3. Choose the one with ``.venv_aeolis``

----

Usage Guide
-------------------------

Option 1: Run via VS Code Terminal
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: powershell

   .venv_aeolis\Scripts\activate
   aeolis run path/to/aeolis.txt

Option 2: Run via Editor (Debug Mode)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Create ``run_console.py``:

.. code-block:: python

   import aeolis

   aeolis.run('path/to/aeolis.txt')

Run with ``F5``.

----

.. tip ::

   If you get:

   ``ModuleNotFoundError: No module named 'pkg_resources'``

   Run:

   .. code-block:: bash

      uv pip install setuptools==70.*

----
