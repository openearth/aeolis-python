import docopt
import numpy as np
from model import AeoLiSRunner, WindGenerator
import log


def aeolis():
    '''aeolis : a process-based model for simulating supply-limited aeolian sediment transport

    Usage:
        aeolis <config> [--callback=FUNC] [--verbose=LEVEL]

    Positional arguments:
        config             configuration file

    Options:
        -h, --help         show this help message and exit
        --callback=FUNC    reference to callback function (e.g. example/callback.py:callback)
        --verbose=LEVEL    write logging messages

    '''
    
    arguments = docopt.docopt(aeolis.__doc__)

    # initialize logger
    if arguments['--verbose'] is not None:
        log.logging.root.setLevel(int(arguments['--verbose']))
    else:
        log.logging.root.setLevel(log.logging.NOTSET)

    # start model
    model = AeoLiSRunner(configfile=arguments['<config>'])
    model.run(callback=arguments['--callback'])


def wind():
    '''aeolis-wind : a wind time series generation tool for the aeolis model

    Usage:
        aeolis-wind <file> [--mean=MEAN] [--max=MAX] [--duration=DURATION] [--timestep=TIMESTEP]

    Positional arguments:
        file               output file

    Options:
        -h, --help           show this help message and exit
        --mean=MEAN          mean wind speed [default: 10]
        --max=MAX            maximum wind speed [default: 30]
        --duration=DURATION  duration of time series [default: 3600]
        --timestep=TIMESTEP  timestep of time series [default: 60]

    '''
    
    arguments = docopt.docopt(wind.__doc__)

    # create random wind time series
    generator = WindGenerator(mean_speed=float(arguments['--mean']),
                              max_speed=float(arguments['--max']),
                              dt=float(arguments['--timestep']))
    generator.generate(duration=float(arguments['--duration']))
    generator.write_time_series(arguments['<file>'])

    u = generator.get_time_series()[1]
    
    fmt = '%-4s : %6.3f m/s'
    print fmt % ('min', np.min(u))
    print fmt % ('mean', np.mean(u))
    print fmt % ('max', np.max(u))


if __name__ == '__main__':
    aeolis()
    
