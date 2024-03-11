
from aeolis.console_debug import aeolis_debug
import cProfile

def main()-> None:
    '''Runs AeoLiS model in debugging mode.'''

    # configfile = r'c:\Users\weste_bt\aeolis\Tests\RotatingWind\Barchan_Grid270\aeolis.txt'
    configfile = r'C:\Data\OneDrive - Lund University\DUNEFORCE\Aeolis\simple_model\testcase\4_test\aeolis.txt'
    aeolis_debug(configfile)

if __name__ == '__main__':
    main()
