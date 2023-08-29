from aeolis.console_debug import aeolis_debug
import cProfile

# configfile = r'c:\Users\weste_bt\aeolis\Tests\RotatingWind\Barchan_Grid270\aeolis.txt'
configfile = r'C:\Users\svries\Documents\GitHub\OE_aeolis-python\examples\deVries2023\Run 2 - Fetch effects and wind directionality\Run2b_play2.txt'
cProfile.run('aeolis_debug(configfile)',sort='tottime')