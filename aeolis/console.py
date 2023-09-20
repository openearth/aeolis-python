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

AeoLiS  Copyright (C) 2023 AeoLiS Development Team

Delft University of Technology
Faculty of Civil Engineering and Geosciences
Stevinweg 1
2628CN Delft
The Netherlands

'''

from __future__ import absolute_import, division

import os
import logging
import numpy as np
from aeolis.model import AeoLiSRunner, WindGenerator

import typer
from typing import Annotated, Optional

app = typer.Typer(add_completion=False, pretty_exceptions_enable=True)


@app.command(help='A process-based model for simulating supply-limited aeolian sediment transport.')
def aeolis(
    config: Annotated[str, typer.Argument(
        help= 'configuration file'
    )], 
    callback: Annotated[Optional[str], typer.Option(
        help='reference to callback function (e.g. example/callback.py:callback)'
    )]=None,

    restart: Annotated[Optional[str], typer.Option(
         help="model restart file"
    )]=None,

    verbose: Annotated[Optional[int], typer.Option(
         help="logging verbosity"
    )]=20,
    debug: Annotated[bool, typer.Option("--debug",
         help="write debug logs",
    )]=False
):

    print_license()

    logger = logging.getLogger('aeolis')

    # start model
    model = AeoLiSRunner(configfile=config)
    model.run(callback=callback,
              restartfile=restart)


@app.command(name="aeolis-wind",
             help="A wind time series generation tool for the aeolis model.")
def wind(
    file: Annotated[str, typer.Argument(
        help="output file"
    )],
    mean: Annotated[float, typer.Option(
        help="mean wind speed"
    )]=10,
    max: Annotated[float, typer.Option(
        help="maximum wind speed"
    )]=30,
    duration: Annotated[int, typer.Option(
        help="duration of time series"
    )]=3600,
    timestep: Annotated[int, typer.Option(
        help="timestep of time series"
    )]=60
):

    print_license()

    # create random wind time series
    generator = WindGenerator(mean_speed=mean,
                              max_speed=max,
                              dt=timestep)
    generator.generate(duration=duration)
    generator.write_time_series(file)

    u = generator.get_time_series()[1]

    fmt = '%-4s : %6.3f m/s'
    print(fmt % ('min', np.min(u)))
    print(fmt % ('mean', np.mean(u)))
    print(fmt % ('max', np.max(u)))


@app.command(name="aeolis-examples",
             help="Sets up data for some examples of aeolis model.")
def examples(
    dir: Annotated[str, typer.Option(
        help="path to directory for the examples' data. Defaults to current directory."
    )]='./'
):

    print_license()

    current_path = os.getcwd()

    if dir == './':
        destination_path = current_path
    else:    
        destination_path = os.path.join(current_path, dir)

    try:
        os.makedirs(destination_path, exist_ok=False)
    except OSError:
        print("Directory %s already exists." % destination_path)

    with open(os.path.join(destination_path, '.txt'), 'w') as f:
        f.write('it works')

    # TODO: copy example files to destination_path


def print_license():
    print('AeoLiS  Copyright (c) 2023  AeoLiS Development Team')
    print('This program comes with ABSOLUTELY NO WARRANTY.')
    print('This is free software, and you are welcome to redistribute it')
    print('under certain conditions; See LICENSE.txt for details.')
    print('')
if __name__ == '__main__':
    
    app()