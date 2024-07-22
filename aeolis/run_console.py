
from aeolis.console_debug import aeolis_debug
import cProfile

def main()-> None:
    '''Runs AeoLiS model in debugging mode.'''

    # configfile = r'c:\Users\weste_bt\aeolis\Tests\RotatingWind\Barchan_Grid270\aeolis.txt'
    configfile = r'C:\Users\svries\Documents\GitHub\OE_aeolis-python\aeolis\examples\grainsizevariations\aeolis_horizontalgradient1.txt'
    aeolis_debug(configfile)

if __name__ == '__main__':
    main()
