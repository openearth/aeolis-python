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
from scipy import ndimage, misc
from scipy.stats import norm, mode
import numpy as np
import math
#import matplotlib.pyplot as plt

# package modules
import aeolis.wind

from aeolis.utils import *


# from aeolis.utils import *

# initialize logger
logger = logging.getLogger(__name__)


def calc_runup(Ho, Tp, beta, nsigma, Kd):
    """
    Calculate runup according to Palmsten and Holman, 2011.
    """
    Lo = 9.81 * Tp * Tp /(2 * np.pi)
    eta = 0.35 * beta * np.sqrt(Ho * Lo)
    sigma_s = np.sqrt(Ho * Lo * (0.563 * (beta * beta) + 0.0004)) * nsigma / 2 / 2
    R = 1.1 * (eta + sigma_s)

    return eta, sigma_s, R


def run_ph12(s, p, t):

    Ho = np.interp(t, p['wave_file'][:, 0], p['wave_file'][:, 1])
    Tp = np.interp(t, p['wave_file'][:, 0], p['wave_file'][:, 2])
    wl = np.interp(t, p['tide_file'][:, 0], p['tide_file'][:, 1])

    zToe = p['dune_toe_elevation']
    beta = p['beach_slope']
    dt = p['dt']

    # wave runup calcs
    nsigma = 2  # Note: nsigma=1 for R16 and nsigma=2 for R2.
    Kd = 1.26  # Coefficient to account for higher runup on dune
    eta, sigma_s, R = calc_runup(Ho, Tp, beta, nsigma, Kd)
    twl = R * Kd + wl


    if twl > zToe:
        ny = p['ny']

        for iy in range(ny+1):

            x = s['x'][iy, :]
            zb = s['zb'][iy, :]

            # parameter set up
            dx = np.abs(x[1] - x[0])
            Bt = beta * 0.54  # dune toe slope trajectory
            Cs = 2.5e-3
            dVResidual_prev = 0  # add this parameter in to be consistent with other codes

            # find dune base location
            st0 = np.nanargmin(np.abs(zb - zToe))
            xToe = x[st0]

            # dune toe trajectory
            zbT = np.ones(len(x)) * np.nan
            zbT[st0:] = Bt * (x[st0:] - x[st0]) + zToe

            # correct toe trajectory that exceeds actual bed elevation
            ix = zbT > zb
            zbT[ix] = zb[ix]

            # initial volume calcs
            Vc = np.cumsum(dx * (zb[st0:] - zbT[st0:]))
            Vc = Vc - Vc[0]

            # collision calcs
            p_collision = 1 - norm.cdf(zToe, eta + wl, sigma_s)
            Nc = p_collision * dt / Tp

            # volume change calcs
            dV = 4 * Cs * (np.max(twl - zToe, 0)) ** 2 * Nc
            dVT = dV - dVResidual_prev

            if dVT < 0:
                ds = 0
            else:
                ds = np.nanargmin(np.abs(Vc - dVT))

            st = st0 + ds
            # x_increment = x[st] - xToe
            dVResidual = Vc[ds] - dVT

            # lets add this residual back to the dune toe so have mass conservation
            dz = -dVResidual / dx
            numcells = np.size(np.arange(st0, st))

            # update morphology
            zb_new = zb
            zb_new[st0:st] = zbT[st0:st]

            #approach to redistribute residual sediment to the lower portion of the dune. needs to be tested
            #if numcells <= 1:
            #    zb_new[st] = zb_new[st] + dz
            #elif numcells < 3:
            #    zb_new[st - 1] = zb_new[st - 1] + dz
            #else:
            #    zb_new[(st0 + 1):st] = zb_new[(st0 + 1):st] + dz / [numcells - 1]

            s['zb'][iy,:] = zb_new
            #s['dVT'] = dVT

    return s