
from aeolis.console_debug import aeolis_debug
import cProfile

def main()-> None:
    '''Runs AeoLiS model in debugging mode. Run this script to start AeoLiS with debugging features enabled, 
    such as step-by-step execution and detailed logging. Useful for development and troubleshooting.
    '''

    # configfile = r'c:\Users\weste_bt\GitHub\AeoLiS\aeolis-python\aeolis\examples\vanWesten2026\01_growth_a\aeolis.txt' # Path to the configuration file
    configfile = r'c:\Users\weste_bt\GitHub\AeoLiS\aeolis-python\voor_floris\01_model\aeolis_base.txt'
    aeolis_debug(configfile)

if __name__ == '__main__':
    main()
