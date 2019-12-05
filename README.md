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
are soil moisture contents, sediment sorting and armouring and
roughness elements.

Documentation can be found at
http://aeolis.readthedocs.io/.

AeoLiS is initially developed at Delft University of Technology with support from the ERC-Advanced
Grant 291206 Nearshore Monitoring and Modeling
([NEMO](http://nemo.citg.tudelft.nl>)) and
[Deltares](http://www.deltares.nl>). AeoLiS is currently maintained by
[Bart van Westen](mailto:Bart.vanWesten@deltares.nl) at Deltares, 
[Nick Cohn](mailto:nick.cohn@usace.army.mil)
[Sierd de Vries](mailto:Sierd.deVries@tudelft.nl) at Delft University of Technology.

## Examples

```
aeolis params.txt
aeolis-wind wind.txt --mean=6 --duration=3600
```
