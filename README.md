![AeoLiS Banner](https://github.com/openearth/aeolis-shortcourse/blob/main/Sandmotor/notebooks/logo.png)

[![Codecov](https://codecov.io/gh/openearth/aeolis-python/branch/master/graph/badge.svg)](https://codecov.io/gh/openearth/aeolis-python)
[![ReadTheDocs](http://readthedocs.org/projects/aeolis/badge/?version=latest)](http://aeolis.readthedocs.io/en/latest/)
[![PyPI](https://img.shields.io/pypi/v/AeoLiS.svg)](https://pypi.python.org/pypi/AeoLiS)
[![PyPI_versions](https://img.shields.io/pypi/pyversions/AeoLiS.svg)](https://pypi.python.org/pypi/AeoLiS)
[![PyPI_status](https://img.shields.io/pypi/status/AeoLiS.svg)](https://pypi.python.org/pypi/AeoLiS)
[![PyPI_format](https://img.shields.io/pypi/format/AeoLiS.svg)](https://pypi.python.org/pypi/AeoLiS)
[![License](https://img.shields.io/pypi/l/AeoLiS.svg)](https://pypi.python.org/pypi/AeoLiS)
[![DOI](https://zenodo.org/badge/7830/openearth/aeolis-python.svg)](https://zenodo.org/badge/latestdoi/7830/openearth/aeolis-python)

# AeoLiS
AeoLiS is a process-based model for simulating aeolian sediment transport in situations where supply-limiting factors are important,
like in coastal environments. Supply-limitations currently supported
are soil moisture contents, sediment sorting and armouring, bed slope
effects, air humidity and roughness elements.

## Installation

**Requirements:**

- Python 3.8 or higher 
- Setuptools
- pip 

### Installing from PyPI

On the comand line of your working environment (Bash/Shell, Conda, Mamba, or similar), run the following: 

```shell
pip install AeoLiS
```

### Installing from source

1. Clone the repository (AEOLIS_V2 branch) using git, or download the source code.

2. Go to the `aeolis-python` directory and install using pip
   ```shell
   cd aeolis-python/
   pip install .
   ```

### Running AeoLiS

Example from command line:

```shell
aeolis params.txt
aeolis-wind wind.txt --mean=6 --duration=3600
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

*AeoLiS Development Team. AeoLiS (Version 2.1.0) [Computer software].* 

## Acknowlegdements

- AeoLiS is supported by the [Digital Competence Centre](https://dcc.tudelft.nl), Delft University of Technology.
- The contributing guidelines for AeoLiS are derived from the [NLeSC/python-template](https://github.com/NLeSC/python-template) and [numpy contributing guide](https://numpy.org/devdocs/dev/index.html#development-process-summary)

&copy; (2023) AeoLiS Development Team, Delft, The Netherlands. 
