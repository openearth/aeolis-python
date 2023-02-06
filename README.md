[![CircleCI](https://circleci.com/gh/openearth/aeolis-python.svg?style=svg)](https://circleci.com/gh/openearth/aeolis-python)
[![Codecov](https://codecov.io/gh/openearth/aeolis-python/branch/master/graph/badge.svg)](https://codecov.io/gh/openearth/aeolis-python)
[![ReadTheDocs](http://readthedocs.org/projects/aeolis/badge/?version=latest)](http://aeolis.readthedocs.io/en/latest/)

[![PyPI](https://img.shields.io/pypi/v/AeoLiS.svg)](https://pypi.python.org/pypi/AeoLiS)
[![PyPI_versions](https://img.shields.io/pypi/pyversions/AeoLiS.svg)](https://pypi.python.org/pypi/AeoLiS)
[![PyPI_status](https://img.shields.io/pypi/status/AeoLiS.svg)](https://pypi.python.org/pypi/AeoLiS)
[![PyPI_format](https://img.shields.io/pypi/format/AeoLiS.svg)](https://pypi.python.org/pypi/AeoLiS)

[![License](https://img.shields.io/pypi/l/AeoLiS.svg)](https://pypi.python.org/pypi/AeoLiS)
[![DOI](https://zenodo.org/badge/7830/openearth/aeolis-python.svg)](https://zenodo.org/badge/latestdoi/7830/openearth/aeolis-python)

# AeoLiS
AeoLiS is a process-based model for simulating aeolian sediment
transport in situations where supply-limiting factors are important,
like in coastal environments. Supply-limitations currently supported
are soil moisture contents, sediment sorting and armouring, bed slope
effects, air humidity and roughness elements.

## Installation

### Installing from source

Requirements:

- Python 3.x
- setuptools
- pip 


1. Clone the repository (AEOLIS_V2 branch)
2. Go to the root directory and install using pip
   ```shell
   pip install .
   ```

### Running AEOLIS

Example:

```shell
aeolis params.txt
aeolis-wind wind.txt --mean=6 --duration=3600
```

## Documentation
Detailed documentation can be found at http://aeolis.readthedocs.io/


## AEOLIS Developer Team
The maintenance and development is done by the AEOLIS developer team:

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

*de Name, A., & Other, A. AeoLiS (Version 2.0.0) [Computer software].* 

## Acknowlegdements

The AEOLIS project is supported by the [Digital Competence Centre](ddc.tudelft.nl), Delft University of Technology.

&copy; (YEAR) [Author(s)], Delft, The Netherlands. 

