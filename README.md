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

Documentation can be found at
http://aeolis.readthedocs.io/.

The maintenance and development is done by the AEOLIS developer team:

Current members of the AEOLIS developer team are:
- [Bart van Westen](mailto:Bart.vanWesten@deltares.nl) at Deltares
- [Nick Cohn](mailto:nick.cohn@usace.army.mil) at U.S. Army Engineer Research and Development Center (ERDC) 
- [Sierd de Vries](mailto:Sierd.deVries@tudelft.nl) (founder) at Delft University of Technology
- [Christa van IJzendoorn](mailto:C.O.vanIJzendoorn@tudelft.nl) at Delft University of Technology
- [Caroline Hallin](mailto:E.C.Hallin@tudelft.nl) at Delft University of Technology
- [Glenn Strypsteen](mailto:glenn.strypsteen@kuleuven.be) at Katholieke Universiteit Leuven
- [Janelle Skaden](mailto:Janelle.E.Skaden@usace.army.mil) at U.S. Army Engineer Research and Development Center (ERDC)

Previous members and contributors are:
- [Bas Hoonhout](mailto:bas@hoonhout.com) (founder) 
- Tom Pak 
- Pieter Rauwoens
- Lisa Meijer

## Examples

```
aeolis params.txt
aeolis-wind wind.txt --mean=6 --duration=3600
```

## Contributing

If you want to contribute to the development of AeoLiS, have a look at the [contributing guidelines](CONTRIBUTING.md).

## Acknowledgements

The contributing guidelines for AeoLiS are derived from the [NLeSC/python-template](https://github.com/NLeSC/python-template) and [numpy contributing guide](https://numpy.org/devdocs/dev/index.html#development-process-summary)