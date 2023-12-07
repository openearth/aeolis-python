from aeolis.console_debug import aeolis_debug
import cProfile

# configfile = r'c:\Users\weste_bt\aeolis\Tests\RotatingWind\Barchan_Grid270\aeolis.txt'
configfile = r'C:\Users\svries\Documents\GitHub\OE_aeolis-python\examples\2D\Barchan_dune\aeolis_sweep_cons.txt'
# cProfile.run('aeolis_debug(configfile)',sort='tottime')
aeolis_debug(configfile)