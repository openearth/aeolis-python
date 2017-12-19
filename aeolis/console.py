'''This file is part of AeoLiS.
   
AeoLiS is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
   
AeoLiS is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
   
You should have received a copy of the GNU General Public License
along with AeoLiS.  If not, see <http://www.gnu.org/licenses/>.
   
AeoLiS  Copyright (C) 2015 Bas Hoonhout

bas.hoonhout@deltares.nl         b.m.hoonhout@tudelft.nl
Deltares                         Delft University of Technology
Unit of Hydraulic Engineering    Faculty of Civil Engineering and Geosciences
Boussinesqweg 1                  Stevinweg 1
2629 HVDelft                     2628CN Delft
The Netherlands                  The Netherlands

'''


from __future__ import absolute_import, division

import os
import docopt
import logging
import numpy as np
from aeolis.model import AeoLiSRunner, WindGenerator


def aeolis():
    '''aeolis : a process-based model for simulating supply-limited aeolian sediment transport

    Usage:
        aeolis <config> [--callback=FUNC] [--restart=FILE] [--verbose=LEVEL]

    Positional arguments:
        config             configuration file

    Options:
        -h, --help         show this help message and exit
        --callback=FUNC    reference to callback function (e.g. example/callback.py:callback)
        --restart=FILE     model restart file
        --verbose=LEVEL    write logging messages [default: 20]

    '''

    print_license()
    
    arguments = docopt.docopt(aeolis.__doc__)

    # initialize file logger
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)-15s %(name)-8s %(levelname)-8s %(message)s',
                        filename='%s.log' % os.path.splitext(arguments['<config>'])[0])

    # initialize console logger
    console = logging.StreamHandler()
    console.setLevel(int(arguments['--verbose']))
    console.setFormatter(logging.Formatter('%(levelname)-8s %(message)s'))
    logging.getLogger('').addHandler(console)

    # start model
    model = AeoLiSRunner(configfile=arguments['<config>'])
    model.run(callback=arguments['--callback'],
              restartfile=arguments['--restart'])


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
    
    print_license()
    
    arguments = docopt.docopt(wind.__doc__)

    # create random wind time series
    generator = WindGenerator(mean_speed=float(arguments['--mean']),
                              max_speed=float(arguments['--max']),
                              dt=float(arguments['--timestep']))
    generator.generate(duration=float(arguments['--duration']))
    generator.write_time_series(arguments['<file>'])

    u = generator.get_time_series()[1]
    
    fmt = '%-4s : %6.3f m/s'
    print(fmt % ('min', np.min(u)))
    print(fmt % ('mean', np.mean(u)))
    print(fmt % ('max', np.max(u)))


def print_license():
    print('AeoLiS  Copyright (C) 2015  Bas Hoonhout')
    print('This program comes with ABSOLUTELY NO WARRANTY.')
    print('This is free software, and you are welcome to redistribute it')
    print('under certain conditions; See LICENSE.txt for details.')
    print('')

    
if __name__ == '__main__':
    aeolis()
    
