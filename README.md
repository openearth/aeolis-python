![AeoLiS Banner](https://github.com/openearth/aeolis-shortcourse/blob/main/Sandmotor/notebooks/logo.png)

[![ReadTheDocs](http://readthedocs.org/projects/aeolis/badge/?version=latest)](http://aeolis.readthedocs.io/en/latest/)
[![PyPI](https://img.shields.io/pypi/v/aeolis.svg)](https://pypi.python.org/pypi/aeolis)
[![PyPI_versions](https://img.shields.io/pypi/pyversions/aeolis.svg)](https://pypi.python.org/pypi/aeolis)
[![PyPI_status](https://img.shields.io/pypi/status/aeolis.svg)](https://pypi.python.org/pypi/aeolis)
[![PyPI_format](https://img.shields.io/pypi/format/aeolis.svg)](https://pypi.python.org/pypi/aeolis)
[![License](https://img.shields.io/pypi/l/aeolis.svg)](https://pypi.python.org/pypi/aeolis)
[![DOI](https://zenodo.org/badge/7830/openearth/aeolis-python.svg)](https://zenodo.org/badge/latestdoi/7830/openearth/aeolis-python)

# AeoLiS
AeoLiS is a process-based model for simulating aeolian sediment transport in situations where supply-limiting factors are important,
like in coastal environments. Supply-limitations currently supported
are soil moisture contents, sediment sorting and armouring, bed slope
effects, air humidity and roughness elements.

https://github.com/openearth/aeolis-python/assets/14054272/128684d6-73ac-4a5f-a186-51559679bd66

## Installation

**Requirements:**

- Python 3.9 or newer
- pip 22.0 or newer
- netCDF4

### Installing from PyPI

On the comand line of your working environment (Bash/Shell, Conda, Mamba, or similar), run the following:

```shell
pip install aeolis
```

> For Windows users, the recommend way to install AeoLiS is to use [Anaconda](https://docs.anaconda.com/free/anaconda/install/windows/).


### Installing from source

1. Clone the repository using Git, or download the source code.

2. AeoLiS users may install the package with only the required dependencies. Go to `aeolis-python` directory and install using pip
   ```shell
   cd aeolis-python/
   pip install .
   ```

3. AeoLiS users who intend to modify the sourcecode can install additional dependencies for test and documentation as follows. Go to root directory `aeolis-python/` and:

   ```shell
   pip install -e .[dev]
   ```

### Running AeoLiS

Examples from command line:

```shell
aeolis run <path/to/aeolis.txt/>
# or wind module
aeolis wind <path/to/wind.txt> --mean=6 --duration=3600
```

## Documentation
Detailed documentation can be found at [AeoLiS ReadTheDocs](http://aeolis.readthedocs.io/)


## AeoLiS Developer Team
The maintenance and development is done by a group of very enthusiastic people.

**Get Involved:**
Read our [Contribution Guidelines](CONTRIBUTING.md) to know how you can help to develop AeoLiS.

**Current Members:**

- [Bart van Westen](mailto:Bart.vanWesten@deltares.nl) at Deltares
- [Nick Cohn](mailto:nick.cohn@usace.army.mil) at U.S. Army Engineer Research and Development Center (ERDC)
- [Sierd de Vries](mailto:Sierd.deVries@tudelft.nl) (founder) at Delft University of Technology
- [Christa van IJzendoorn](mailto:C.O.vanIJzendoorn@tudelft.nl) at Delft University of Technology
- [Caroline Hallin](mailto:E.C.Hallin@tudelft.nl) at Delft University of Technology
- [Glenn Strypsteen](mailto:glenn.strypsteen@kuleuven.be) at Katholieke Universiteit Leuven
- [Janelle Skaden](mailto:Janelle.E.Skaden@usace.army.mil) at U.S. Army Engineer Research and Development Center (ERDC)

**Previous Members & Contributors:**
- [Bas Hoonhout](mailto:bas@hoonhout.com) (founder)
- Tom Pak
- Pieter Rauwoens
- Lisa Meijer

## Citation

Please, cite this software as follows:

*de Vries, S., Hallin, C., van IJzendoorn, C., van Westen, B., Cohn, N., Strypsteen, G., Skaden, J., Agrawal, N., & Garcia Alvarez, M. (2023). AeoLiS (Version 3.0.0.rc2) [Computer software]. https://github.com/openearth/aeolis-python*

## Acknowlegdements

- AeoLiS is supported by the [Digital Competence Centre](https://dcc.tudelft.nl), Delft University of Technology.
- The contributing guidelines for AeoLiS are derived from the [NLeSC/python-template](https://github.com/NLeSC/python-template) and [numpy contributing guide](https://numpy.org/devdocs/dev/index.html#development-process-summary)

&copy; (2023) AeoLiS Development Team, Delft, The Netherlands.
