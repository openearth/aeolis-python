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

import logging
import numpy as np
from aeolis.model import AeoLiSRunner


def aeolis_debug(configfile):
    '''aeolis : a process-based model for simulating supply-limited aeolian sediment transport
    '''

    print_license()
    logger = logging.getLogger('aeolis')

    # start model
    model = AeoLiSRunner(configfile=configfile)
    model.run()


def print_license():
    print('AeoLiS  Copyright (C) 2015  Bas Hoonhout')
    print('This program comes with ABSOLUTELY NO WARRANTY.')
    print('This is free software, and you are welcome to redistribute it')
    print('under certain conditions; See LICENSE.txt for details.')
    print('')


# if __name__ == '__main__':
    # aeolis()
