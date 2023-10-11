"""This file is part of AeoLiS.

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

"""

from __future__ import absolute_import, division

import os
import sys
import logging
import shutil

import numpy as np
import typer
from typing import Annotated, Optional

import aeolis as aeolis_pkg
from aeolis.model import AeoLiSRunner, WindGenerator


aeolis_app = typer.Typer(
    add_completion=False,
    pretty_exceptions_enable=True,
    help=(
        "AeoLiS, a process-based model for simulating supply-limited aeolian"
        " sediment transport"
    ),
)


@aeolis_app.command(
    help="runs model for simulating supply-limited aeolian sediment transport."
)
def run(
    config: Annotated[str, typer.Argument(help="configuration file")],
    callback: Annotated[
        Optional[str],
        typer.Option(
            help=(
                "reference to callback function (e.g."
                " example/callback.py:callback)"
            )
        ),
    ] = None,
    restart: Annotated[
        Optional[str], typer.Option(help="model restart file")
    ] = None,
    verbose: Annotated[
        Optional[int], typer.Option(help="logging verbosity")
    ] = 20,
    debug: Annotated[
        bool,
        typer.Option(
            "--debug",
            help="write debug logs",
        ),
    ] = False,
):
    print_license()

    logger = logging.getLogger("aeolis")

    # start model
    model = AeoLiSRunner(configfile=config)
    model.run(callback=callback, restartfile=restart)


@aeolis_app.command(
    name="wind", help="A wind time series generation tool for the aeolis model."
)
def wind(
    file: Annotated[str, typer.Argument(help="output file")],
    mean: Annotated[float, typer.Option(help="mean wind speed")] = 10,
    max: Annotated[float, typer.Option(help="maximum wind speed")] = 30,
    duration: Annotated[
        int, typer.Option(help="duration of time series")
    ] = 3600,
    timestep: Annotated[int, typer.Option(help="timestep of time series")] = 60,
):
    print_license()

    # create random wind time series
    generator = WindGenerator(mean_speed=mean, max_speed=max, dt=timestep)
    generator.generate(duration=duration)
    generator.write_time_series(file)

    u = generator.get_time_series()[1]

    fmt = "%-4s : %6.3f m/s"
    print(fmt % ("min", np.min(u)))
    print(fmt % ("mean", np.mean(u)))
    print(fmt % ("max", np.max(u)))


@aeolis_app.command(name="examples", help="sets up some examples for AeoLiS.")
def examples(
    dir: Annotated[
        str,
        typer.Argument(
            help=(
                "directory for the examples' data. It will be created if it"
                " does not exist."
            )
        ),
    ]
):
    print_license()

    current_path = os.getcwd()
    destination_path = os.path.join(current_path, dir)

    destination_path = os.path.join(dir, "aeolis-examples")
    # Path to examples in installed aeolis package directory
    source_path = os.path.join(os.path.dirname(aeolis_pkg.__file__), "examples")

    # Examples that will be included in destination path
    # Key: path where the example is located in the installed aeolis package directory
    # Value: path where the example will be copied to
    _examples = {
        "2D/Parabolic_dune": "Parabolic_dune",
        "sandengine_small_grids": "sandengine_small_grids",
    }
    # Copy examples to target path. Creates destination path if it does not exist.
    for source, target in _examples.items():
        example_path = os.path.join(source_path, source)
        target_path = os.path.join(destination_path, target)
        shutil.copytree(example_path, target_path)

    print("Done. Examples copied to: \n %s" % destination_path)


def print_license():
    print("AeoLiS  Copyright (c) 2023  AeoLiS Development Team")
    print("This program comes with ABSOLUTELY NO WARRANTY.")
    print("This is free software, and you are welcome to redistribute it")
    print("under certain conditions; See LICENSE.txt for details.")
    print("")


def start_aeolis_app():
    if len(sys.argv) == 2 and os.path.isfile(sys.argv[1]):
        print(
            "\nUsage of the command line syntax `aeolis <path_to_aeolis.txt>`"
            " has been deprecated from v3.0.0 onwards.\n\nTo run a model, use"
            " the syntax `aeolis run <path_to_aeolis.txt>`\n"
        )
        sys.exit(1)
    else:
        aeolis_app()


def start_aeolis_wind_app():
    """
    This function serves the purpose of catching the deprecated command line syntx `aeolis-wind <path_to_wind.txt>` and pritning a message on the console instructing the user to use `aeolis wind <path_to_wind.txt>` instead.
    """
    print(
        "\nUsage of the command line syntax `aeolis-wind <path_to_wind.txt>`"
        " has been deprecated from v3.0.0 onwards.\n\nTo run the wind module,"
        " use the syntax `aeolis wind <path_to_wind.txt>`\n"
    )
    sys.exit(1)


if __name__ == "__main__":
    start_aeolis_app()
    start_aeolis_wind_app()
