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

import logging
import numpy as np
import scipy.special
import scipy.interpolate
from scipy import ndimage, misc
#import matplotlib
#import matplotlib.pyplot as plt
from builtins import range, int
import math
from collections import namedtuple
from copy import copy
from pprint import pprint as pp
import sys
import os
from aeolis.wind import velocity_stress

# package modules
from aeolis.utils import *

# initialize logger
logger = logging.getLogger(__name__)

def initialize(s,p):
    if p['process_fences']:
        s['fence_height'][:,:] = p['fence_file']
        s['fence_base'] = copy(s['zb'])  # initial fence base is the bed elevation
        s['fence_top'] = s['fence_base'] + s['fence_height']
        s['fence_height_init'] = s['fence_height']
        s['zf'] = s['fence_height']
    return s

def update_fences(s,p):
    s = update_fence_height(s, p)
    if p['ny'] > 0:
        s = fence_shear2d(s, p)
    else:
        s = fence_shear1d(s, p)

    s = velocity_stress(s,p)
    return s

def update_fence_height(s, p):
    s['fence_height'] = s['fence_top']-s['zb']
    ix = s['fence_height_init'] < 0.1
    s['fence_height'][ix] = 0
    ix = s['fence_top'] < 0.1
    s['fence_height'][ix] = 0
    ix = s['fence_height'] < 0.1
    s['fence_height'][ix] = 0

    #if exceeds 1.5 m then assume the fence has eroded out
    ix = s['fence_height'] > 1.5
    s['fence_height'][ix] = 0
    s['fence_height_init'][ix] = 0

    return s


def fence_shear2d(s, p):
    x = s['x']
    y = s['y']
    zf = s['fence_height']
    ustarx = s['ustars']
    ustary = s['ustarn']

    #dx = p['dx']/10
    #dy = p['dx']/10

    dx = np.maximum(p['dx']/2, 0.25)
    dy = dx


    udir = s['udir'][0, 0]
    if udir < 0:
        udir = udir + 360
    udir = udir - np.floor(udir / 360) * 360

    if udir == 0 or udir == 360 or udir == -360 or udir == -180 or udir == 180:
        udir += 0.00001

    igrid, cgrid, x0, y0 = initialize_computational_grid(x, y, zf, ustarx, ustary, dx, dy, buffer_width=100)
    igrid, cgrid = calc_fence_shear(igrid, cgrid, udir, x0, y0, p)

    s['ustars'] = igrid['ustarx']
    s['ustarn'] = igrid['ustary']
    s['ustar'] = np.sqrt(s['ustars']**2 + s['ustarn']**2)

    return s


def initialize_computational_grid(x, y, z, ustarx, ustary, dx, dy, buffer_width=100., buffer_relaxation=None):
    if buffer_relaxation is None:
        buffer_relaxation = buffer_width / 4.

    mult_all = np.ones(x.shape)

    igrid = dict(x=x,
                 y=y,
                 z=z,
                 ustarx=ustarx,
                 ustary=ustary,
                 mult_all=mult_all)

    cgrid = dict(dx=dx,
                 dy=dy)

    x0, y0, cgrid = set_computational_grid(igrid, cgrid, buffer_width)

    return igrid, cgrid, x0, y0


def set_computational_grid(igrid, cgrid, buffer_width):
    '''Define computational grid

    The computational grid is square with dimensions equal to the
    diagonal of the bounding box of the input grid, plus twice the
    buffer width.

    '''

    gi = igrid
    gc = cgrid

    # grid center
    x0, y0 = np.mean(gi['x']), np.mean(gi['y'])

    # grid size
    D = np.sqrt((gi['x'].max() - gi['x'].min()) ** 2 +
                (gi['y'].max() - gi['y'].min()) ** 2) + 2 * buffer_width

    # determine equidistant, square grid
    xc, yc = get_exact_grid(x0 - D / 2., x0 + D / 2.,
                            y0 - D / 2., y0 + D / 2.,
                            gc['dx'], gc['dy'])

    gc['xi'] = xc
    gc['yi'] = yc

    return x0, y0, gc


def calc_fence_shear(igrid, cgrid, udir, x0, y0, p):
    '''Compute wind shear for given wind speed and direction

    Parameters
    ----------
    u0 : float
        Free-flow wind speed
    udir : float
        Wind direction in degrees
    process_separattion :

    '''
    gc = cgrid  # computational grid
    gi = igrid  # initial grid

    # Populate computational grid (rotate to wind direction + interpolate input topography)
    populate_computational_grid(igrid, cgrid, udir + 90., x0, y0)

    # Compute wind shear stresses on computational grid
    gc = compute_fenceshear(gi, gc, udir, p)

    ustarx_init = gc['ustarx'][0,0]
    gc['ustarx'] = gc['ustarx'] * gc['mult_all']
    gc['ustary'] = np.zeros(gc['x'].shape)

    #ensure bad data doesnt make it through
    #ix = gc['mindist'] > 20
    #gc['ustarx'][ix] = ustarx_init


    gc['ustarx'], gc['ustary'] = rotate(gc['ustarx'], gc['ustary'], udir + 90)

    # Rotate both (i&c) grids + results in opposite dir.
    gi['x'], gi['y'] = rotate(gi['x'], gi['y'], -(udir + 90.), origin=(x0, y0))

    gc['x'], gc['y'] = rotate(gc['x'], gc['y'], -(udir + 90.), origin=(x0, y0))

    gc['ustary'], gc['ustarx'] = rotate(gc['ustarx'], gc['ustary'], -(udir + 90))

    # Interpolate wind shear results to real grid
    gi['ustarx'] = interpolate(gc['x'], gc['y'], gc['ustarx'],
                               gi['x'], gi['y'])
    gi['ustary'] = interpolate(gc['x'], gc['y'], gc['ustary'],
                               gi['x'], gi['y'])

    # Rotate real grid and wind shear results back to orignal orientation
    gc['x'], gc['y'] = rotate(gc['x'], gc['y'], udir + 90., origin=(x0, y0))
    gi['x'], gi['y'] = rotate(gi['x'], gi['y'], +(udir + 90.), origin=(x0, y0))

    gi['ustarx'], gi['ustary'] = rotate(gi['ustarx'], gi['ustary'], +(udir + 90))

    #avoid any boundary effects
    #gi['ustarx'][1,:] = gi['ustarx'][2,:]
    #gi['ustarx'][:,1] = gi['ustarx'][:,2]
    #gi['ustarx'][-2,:] = gi['ustarx'][-3,:]
    #gi['ustarx'][:,-2] = gi['ustarx'][:,-3]

    #gi['ustary'][1,:] = gi['ustary'][1,:]
    #gi['ustary'][:,1] = gi['ustary'][:,1]
    #gi['ustary'][-2,:] = gi['ustary'][-2,:]
    #gi['ustary'][:,-2] = gi['ustary'][:,-2]

    gi['ustarx'][0,:] = gi['ustarx'][1,:]
    gi['ustarx'][:,0] = gi['ustarx'][:,1]
    gi['ustarx'][-1,:] = gi['ustarx'][-2,:]
    gi['ustarx'][:,-1] = gi['ustarx'][:,-2]

    gi['ustary'][0,:] = gi['ustary'][1,:]
    gi['ustary'][:,0] = gi['ustary'][:,1]
    gi['ustary'][-1,:] = gi['ustary'][-2,:]
    gi['ustary'][:,-1] = gi['ustary'][:,-2]

    return gi, gc

    # Input functions for __call()


def populate_computational_grid(igrid, cgrid, alpha, x0, y0):
    '''Interpolate input topography to computational grid

    Adds and fills buffer zone around the initial grid and
    rotates the computational grid to current wind direction.
    The computational grid is filled by interpolating the input
    topography and initial wind induced shear stresses to it.

    Parameters
    ----------
    alpha : float
        Rotation angle in degrees

    '''
    gi = igrid
    gc = cgrid
    x = gi['x']
    shp = x.shape
    try:
        ny = shp[2]
    except:
        ny = 0

    # Add buffer zone around grid                                           # buffer is based on version bart, sigmoid function is no longer required
    if ny <= 0:
        dxi = gi['x'][0, 0]
        dyi = gi['y'][0, 0]
    else:
        dxi = gi['x'][1, 1] - gi['x'][0, 0]
        dyi = gi['y'][1, 1] - gi['y'][0, 0]

    buf = 100  # amount of cells

    xi, yi = np.meshgrid(
        np.linspace(gi['x'][0, 0] - buf * dxi, gi['x'][-1, -1] + buf * dxi, gi['x'].shape[1] + 2 * buf),
        np.linspace(gi['y'][0, 0] - buf * dyi, gi['y'][-1, -1] + buf * dyi, gi['y'].shape[0] + 2 * buf))

    zi = np.zeros((xi.shape))
    zi[buf:-buf, buf:-buf] = gi['z']

    # Filling buffer zone edges
    zi[buf:-buf, :buf] = np.repeat(zi[buf:-buf, buf + 1][:, np.newaxis], buf, axis=1)
    zi[buf:-buf, -buf:] = np.repeat(zi[buf:-buf, -buf - 1][:, np.newaxis], buf, axis=1)

    zi[:buf, buf:-buf] = np.repeat(zi[buf + 1, buf:-buf][np.newaxis], buf, axis=0)
    zi[-buf:, buf:-buf] = np.repeat(zi[-buf - 1, buf:-buf][np.newaxis], buf, axis=0)

    # Filling buffer zone corners
    zi[:buf, :buf] = zi[buf + 1, buf + 1]
    zi[-buf:, :buf] = zi[-buf - 1, buf + 1]
    zi[:buf, -buf:] = zi[buf + 1, -buf - 1]
    zi[-buf:, -buf:] = zi[-buf - 1, -buf - 1]

    # Rotate computational grid to the current wind direction
    xc, yc = rotate(gc['xi'], gc['yi'], alpha, origin=(x0, y0))

    # Interpolate input topography to computational grid
    zfc = interpolate(gi['x'], gi['y'], gi['z'], xc, yc)
    ustarxc = interpolate(gi['x'], gi['y'], gi['ustarx'], xc, yc)
    ustaryc = interpolate(gi['x'], gi['y'], gi['ustary'], xc, yc)

    # Interpolate input wind - shear
    # tauxc = interpolate(gi['x'], gi['y'], gi['taux'], xc, yc)
    # tauyc = interpolate(gi['x'], gi['y'], gi['tauy'], xc, yc)

    gc['x'] = xc
    gc['y'] = yc
    gc['z'] = zfc

    # gc['taux'] = tauxc
    # gc['tauy'] = tauyc
    gc['zfc'] = zfc
    gc['ustarx'] = ustarxc
    gc['ustary'] = ustaryc

    return gc


def compute_fenceshear(igrid, cgrid, udir, p):
    '''Compute wind shear perturbation for given free-flow wind
    speed on computational grid

    Parameters
    ----------
    u0 : float
        Free-flow wind speed
    nfilter : 2-tuple
        Wavenumber range used for logistic sigmoid filter. See
        :func:`filter_highfrequencies`

    '''
    gc = cgrid
    zf = gc['z']
    ny, nx = gc['z'].shape
    mult_all = np.ones(zf.shape)
    #mindist = np.ones(zf.shape) * 1000

    for iy in range(ny):

        # intialize other grid parameters
        x = gc['x'][iy, :]
        zp = gc['z'][iy, :]
        red_all = np.zeros(x.shape)

        nx2 = x.size
        c1 = p['okin_c1_fence']
        intercept = p['okin_initialred_fence']

        for igrid in range(nx2):

            # only look at cells with a roughness element
            if zp[igrid] > 0:

                # print(zp[igrid])
                # local parameters
                if udir >= 180 and udir <= 360:
                    xrel = x - x[igrid]
                else:
                    xrel = -(x - x[igrid])

                red = np.zeros(x.shape)
                mult = np.ones(x.shape)
                h = zp[igrid]

                for igrid2 in range(nx2):
                    if xrel[igrid2] >= 0 and xrel[igrid2] / h < 20:
                        # apply okin model

                        # apply okin model
                        mult[igrid2] = intercept + (1 - intercept) * (1 - math.exp(-xrel[igrid2] * c1 / h))

                        #ifind = xrel > 0
                        #if np.size(ifind) > 0:
                        #    mindist[iy, igrid2] = np.minimum(np.min(xrel[ifind]), mindist[iy, igrid2])

                red = 1 - mult

                # fix potential issues for summation
                ix = red < 0.00001
                red[ix] = 0
                ix = red > 1
                red[ix] = 1
                ix = xrel < 0
                red[ix] = 0

                # combine all reductions between plants
                red_all = red_all + red

        # cant have more than 100% reduction
        ix = red_all > 1
        red_all[ix] = 1

        # convert to a multiple
        mult_all[iy, :] = 1 - red_all
        #mult_all[iy, igrid2] = mult_all[iy, igrid2] - mult_temp

    #avoid any boundary effects
    mult_all[0,:] = 1
    mult_all[:,0] = 1
    mult_all[-1,:] = 1
    mult_all[:,-1] = 1

    cgrid['mult_all'] = mult_all
    #cgrid['mindist'] = mindist

    return cgrid


def get_exact_grid(xmin, xmax, ymin, ymax, dx, dy):
    '''Returns a grid with given gridsizes approximately within given bounding box'''

    x = np.arange(np.floor(xmin / dx) * dx,
                  np.ceil(xmax / dx) * dx, dx)
    y = np.arange(np.floor(ymin / dy) * dy,
                  np.ceil(ymax / dy) * dy, dy)
    x, y = np.meshgrid(x, y)

    return x, y


def rotate(x, y, alpha, origin=(0, 0)):
    '''Rotate a matrix over given angle around given origin'''

    xr = x - origin[0]
    yr = y - origin[1]

    a = alpha / 180. * np.pi

    R = np.asmatrix([[np.cos(a), -np.sin(a)],
                     [np.sin(a), np.cos(a)]])

    xy = np.concatenate((xr.reshape((-1, 1)),
                         yr.reshape((-1, 1))), axis=1) * R

    return (np.asarray(xy[:, 0].reshape(x.shape) + origin[0]),
            np.asarray(xy[:, 1].reshape(y.shape) + origin[1]))


def interpolate(x, y, z, xi, yi):
    '''Interpolate one grid to an other'''

    xy = np.concatenate((y.reshape((-1, 1)),
                         x.reshape((-1, 1))), axis=1)

    xyi = np.concatenate((yi.reshape((-1, 1)),
                          xi.reshape((-1, 1))), axis=1)

    # version Bart
    inter = scipy.interpolate.RegularGridInterpolator((y[:, 0], x[0, :]), z, bounds_error=False, fill_value=0.)
    zi = inter(xyi).reshape(xi.shape)

    return zi

def fence_shear1d(s, p):


    ustar = s['ustar'].copy()
    ustars = s['ustars'].copy()
    ustarn = s['ustarn'].copy()
    ets = np.zeros(s['zb'].shape)
    etn = np.zeros(s['zb'].shape)
    ix = ustar != 0
    ets[ix] = ustars[ix] / ustar[ix]
    etn[ix] = ustarn[ix] / ustar[ix]
    udir = s['udir'][0,0]+180

    x = s['x'][0,:]
    zp = s['fence_height'][0,:]
    red = np.zeros(x.shape)
    red_all = np.zeros(x.shape)
    nx = x.size
    c1 = p['okin_c1_fence']
    intercept = p['okin_initialred_fence']

    if udir < 360:
        udir = udir + 360

    if udir > 360:
        udir = udir - 360


    #Calculate shear reduction by looking through all cells that have plants present and looking downwind of those features
    for igrid in range(nx):

        if zp[igrid] > 0:         # only look at cells with a roughness element
            mult = np.ones(x.shape)
            h = zp[igrid] #vegetation height at the appropriate cell

            if udir >= 180 and udir <= 360:
                xrel = -(x - x[igrid])
            else:
                xrel = x - x[igrid]

            for igrid2 in range(nx):

                if xrel[igrid2] >= 0 and xrel[igrid2]/h < 20:

                    # apply okin model
                    mult[igrid2] = intercept + (1 - intercept) * (1 - math.exp(-xrel[igrid2] * c1 / h))

            red = 1 - mult

            # fix potential issues for summation
            ix = red < 0.00001
            red[ix] = 0
            ix = red > 1
            red[ix] = 1
            ix = xrel < 0
            red[ix] = 0

            # combine all reductions between plants
            red_all = red_all + red

    # cant have more than 100% reduction
    ix = red_all > 1
    red_all[ix] = 1

    # convert to a multiple
    mult_all = 1 - red_all
    ustarfence = s['ustar'][0,:] * mult_all
    ix = ustarfence < 0.01
    ustarfence[ix] = 0.01 #some small number so transport code doesnt crash

    s['ustar'][0,:] = ustarfence
    s['ustars'][0,:] = s['ustar'][0,:] * ets[0,:]
    s['ustarn'][0,:] = s['ustar'][0,:] * etn[0,:]

    return s