# AeoLiS
AeoLiS is a process-based model for simulating aeolian sediment
transport in situations where supply-limiting factors are important,
like in coastal environments. Supply-limitations currently supported
are soil moisture contents, sediment sorting and armouring and
roughness elements.

Documentation can be found at
http://aeolis.readthedocs.io/.

AeoLiS is initially developed by [Bas Hoonhout](b.m.hoonhout@tudelft.nl)
at Delft University of Technology with support from the ERC-Advanced
Grant 291206 Nearshore Monitoring and Modeling
([NEMO](http://nemo.citg.tudelft.nl>)) and
[Deltares](http://www.deltares.nl>). AeoLiS is currently maintained by
[Bas Hoonhout](bas.hoonhout@deltares.nl) at Deltares and
[Sierd de Vries](Sierd.deVries@tudelft.nl) at Delft University of Technology.

## Examples

```
aeolis params.txt
aeolis-wind wind.txt --mean=6 --duration=3600
```
