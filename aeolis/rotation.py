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
import aeolis.vegetation
from aeolis.utils import *

# initialize logger
logger = logging.getLogger(__name__)


class rotationClass:
    '''XXXX.
    '''

    igrid = {}
    cgrid = {}
    istransect = False

    def __init__(self, x, y, z, dx, dy, buffer_width):
        '''Class initialization

        Parameters
        ----------
        x : numpy.ndarray
            2D array with x-coordinates of input grid
        y : numpy.ndarray
            2D array with y-coordinates of input grid
        z : numpy.ndarray
            2D array with topography of input grid
        dx : float, optional
            Grid spacing in x dimension of computational grid
            (default: 1)
        dy : float, optional
            Grid spacing of y dimension of computational grid
            (default: 1)
        buffer_width : float, optional
            Width of buffer distance between input grid boundary and
            computational grid boundary (default: 100)
        '''

        if z.shape[0] == 1:
            self.istransect = True

        # Assigning values to original (i) and computational (c) grid
        self.cgrid = dict(dx=dx, dy=dy)

        # Setting buffer settings
        self.buffer_width = buffer_width

    def __call__(self, x, y, z, udir, ustar, tau, hveg, type_flag, params):
        # Reload x and y because of horizontalized input-grid
        self.igrid = dict(x=x, y=y, z=z, ustar=ustar, hveg=hveg, tau=tau)

        # Convert to cartesian to perform all the rotations
        u_angle = 270. - udir  # wind angle

        # Creating the computational grid
        self.set_computational_grid(udir)

        # Storing computational (c) and original (i) grids
        gi = self.igrid  # initial grid
        gc = self.cgrid  # computational grid

        # Rotate computational (c) grid to the current wind direction
        gc['x'], gc['y'] = self.rotate(gc['xi'], gc['yi'], -u_angle, origin=(self.x0, self.y0))
        # Interpolate bed levels and shear to the computational grid
        gc['z'] = self.interpolate(gi['x'], gi['y'], gi['z'], gc['x'], gc['y'], 0)
        gc['x'], gc['y'] = self.rotate(gc['x'], gc['y'], u_angle, origin=(self.x0, self.y0))
        gi['x'], gi['y'] = self.rotate(gi['x'], gi['y'], u_angle, origin=(self.x0, self.y0))

        #Do relevant calculations here
        if type_flag == 1: #okin calcs
            gc['ustar'] = self.interpolate(gi['x'], gi['y'], gi['ustar'], gc['x'], gc['y'], 0)
            gc['hveg'] = self.interpolate(gi['x'], gi['y'], gi['hveg'], gc['x'], gc['y'], 0)
            gc['ustar'] = aeolis.vegetation.compute_okin_shear(gc['x'], gc['ustar'], gc['hveg'], params)
            gi['ustar'] = self.interpolate(gc['x'], gc['y'], gc['ustar'], gi['x'], gi['y'], 0)
        if type_flag == 2: #sandfence calcs
            gc['ustar'] = self.interpolate(gi['x'], gi['y'], gi['ustar'], gc['x'], gc['y'], 0)
            gc['fence_height'] = self.interpolate(gi['x'], gi['y'], gi['fence_height'], gc['x'], gc['y'], 0)
            gc['ustar'] = aeolis.fences.compute_fence_shear(gc['x'], gc['ustar'], gc['fence_height'], params)
            gi['ustar'] = self.interpolate(gc['x'], gc['y'], gc['ustar'], gi['x'], gi['y'], 0)
        if type_flag == 3: #shear perturbation calcs
            gc['tau'] = self.interpolate(gi['x'], gi['y'], gi['tau'], gc['x'], gc['y'], 0)
            if params['shear_type'] == 'duna2d':
                gc['tau'] = aeolis.wind.compute_shear_perturbation_1D(gc['x'], gc['z'], gc['tau'], params)
            if params['shear_type'] == 'quasi2d':
                gc['tau'] = aeolis.wind.compute_shear_perturbation_kroy1D(gc['x'], gc['y'], gc['z'], gc['tau'], params)
            gc['zsep'], gc['hsep'], gc['tau'] = aeolis.wind.separation_quasi2d(gc['x'], gc['y'], gc['z'], gc['tau'], gc['dx'], gc['dy'], params)
            gi['zsep'] = self.interpolate(gc['x'], gc['y'], gc['zsep'], gi['x'], gi['y'], 0)
            gi['hsep'] = self.interpolate(gc['x'], gc['y'], gc['hsep'], gi['x'], gi['y'], 0)
            gi['tau'] = self.interpolate(gc['x'], gc['y'], gc['tau'], gi['x'], gi['y'], 0)

        self.igrid = gi # initial grid
        self.cgrid = gc  # computational grid

        return self

    # Input functions for __call()
    def set_computational_grid(self, udir):
        '''Define computational grid

        The computational grid is square with dimensions equal to the
        diagonal of the bounding box of the input grid, plus twice the
        buffer width.

        '''

        # Copying the original (i) and computational (c) grid
        gi = self.igrid
        gc = self.cgrid

        # Compute grid center, same for both original (i) and computational (c) grid
        x0, y0 = np.mean(gi['x']), np.mean(gi['y'])

        # Initialization
        b_W = np.zeros(4)
        b_L = np.zeros(4)
        xcorner = np.zeros(4)
        ycorner = np.zeros(4)

        # Computing the corner-points of the grid
        xcorner[0] = gi['x'][0, 0]
        ycorner[0] = gi['y'][0, 0]
        xcorner[1] = gi['x'][-1, 0]
        ycorner[1] = gi['y'][-1, 0]
        xcorner[2] = gi['x'][0, -1]
        ycorner[2] = gi['y'][0, -1]
        xcorner[3] = gi['x'][-1, -1]
        ycorner[3] = gi['y'][-1, -1]

        # Preventing vertical lines
        udir_verticals = np.arange(-1080, 1080, 90)
        udir_vertical_bool = False
        for udir_vertical in udir_verticals:
            if (abs(udir - udir_vertical) <= 0.001):
                udir_vertical_bool = True
        if udir_vertical_bool:
            udir -= 0.1

        # Compute slope (m) and intercept (b) from parallel lines along all (4) grids corners
        for i in range(4):
            # Parallel boundaries
            m_W, b_W[i] = np.polyfit([xcorner[i], xcorner[i] - np.sin(np.deg2rad(udir))],
                                     [ycorner[i], ycorner[i] - np.cos(np.deg2rad(udir))], 1)
            # Perpendicular boundaries
            m_L, b_L[i] = np.polyfit([xcorner[i], xcorner[i] - np.sin(np.deg2rad(udir - 90.))],
                                     [ycorner[i], ycorner[i] - np.cos(np.deg2rad(udir - 90.))], 1)

        # Determine the most outer boundaries (for parallel and perpendicular)
        db_W = self.maxDiff(b_W)
        db_L = self.maxDiff(b_L)

        # Compute the distance between the outer boundaries to determine the width (W) and length (L) of the grid
        self.Width = abs(db_W) / np.sqrt((m_W ** 2.) + 1) + self.buffer_width * 2.
        self.Length = abs(db_L) / np.sqrt((m_L ** 2.) + 1) + self.buffer_width * 2.

        # Create the grid
        xc, yc = self.get_exact_grid(x0 - self.Length / 2., x0 + self.Length / 2.,
                                     y0 - self.Width / 2., y0 + self.Width / 2.,
                                     gc['dx'], gc['dy'])

        # Storing grid parameters
        self.x0 = x0
        self.y0 = y0
        gc['xi'] = xc
        gc['yi'] = yc

        return self

    @staticmethod
    def get_exact_grid(xmin, xmax, ymin, ymax, dx, dy):
        '''Returns a grid with given gridsizes approximately within given bounding box'''

        x = np.arange(np.floor(xmin / dx) * dx,
                      np.ceil(xmax / dx) * dx, dx)
        y = np.arange(np.floor(ymin / dy) * dy,
                      np.ceil(ymax / dy) * dy, dy)
        x, y = np.meshgrid(x, y)

        return x, y

    @staticmethod
    def get_borders(x):
        '''Returns borders of a grid as one-dimensional array'''

        return np.concatenate((x[0, :].T,
                               x[1:-1, -1],
                               x[-1, ::-1].T,
                               x[-1:1:-1, 0],
                               x[0, :1]), axis=0)

    @staticmethod
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

    def get_ustar(self):
        '''Returns ustar
        '''

        ustar = self.igrid['ustar']
        return ustar

    def get_tau(self):
        '''Returns ustar
        '''

        tau = self.igrid['tau']
        return tau

    def maxDiff(self, arr):

        result = 0
        n = len(arr)

        # Iterate through all pairs.
        for i in range(0, n):
            for j in range(i, n):

                if (abs(arr[i] - arr[j]) + abs(i - j) > result):
                    result = abs(arr[i] - arr[j]) + abs(i - j)

        return result

    def interpolate(self, x, y, z, xi, yi, z0):
        '''Interpolate one grid to an other'''

        # First compute angle with horizontal
        dx = x[0, 1] - x[0, 0]
        dy = y[0, 1] - y[0, 0]

        angle = np.rad2deg(np.arctan(dy / dx))

        if dx <= 0 and dy <= 0:
            angle += 180.

        # Rotate grids to allign with horizontal
        x, y = self.rotate(x, y, angle, origin=(self.x0, self.y0))
        xi, yi = self.rotate(xi, yi, angle, origin=(self.x0, self.y0))

        # Rotate 180 deg if necessary
        if not np.all(sorted(y[:, 0]) == y[:, 0]) and not np.all(sorted(x[0, :]) == x[0, :]):
            x, y = self.rotate(x, y, 180, origin=(self.x0, self.y0))
            xi, yi = self.rotate(xi, yi, 180, origin=(self.x0, self.y0))

        # Concatenate
        xy = np.concatenate((y.reshape((-1, 1)),
                             x.reshape((-1, 1))), axis=1)

        xyi = np.concatenate((yi.reshape((-1, 1)),
                              xi.reshape((-1, 1))), axis=1)

        # Interpolate
        pad_w = np.maximum(np.shape(x)[0], np.shape(x)[1])
        x_pad = np.pad(x, ((pad_w, pad_w), (pad_w, pad_w)), 'reflect', reflect_type='odd')
        y_pad = np.pad(y, ((pad_w, pad_w), (pad_w, pad_w)), 'reflect', reflect_type='odd')
        z_pad = np.pad(z, ((pad_w, pad_w), (pad_w, pad_w)), 'edge')

        if self.istransect:
            zi = np.interp(xi.flatten(), x_pad.flatten(), z_pad.flatten()).reshape(xi.shape)
        else:
            # in the scipy 1.10 version the regular grid interpolator does not work with non c-contigous arrays.
            # Here we make a copy as a dirty solution feeding the interpolator with ordered copies
            inter = scipy.interpolate.RegularGridInterpolator(
                (y_pad[:, 0].copy(order='C'), x_pad[0, :].copy(order='C')), z_pad, bounds_error=False, fill_value=z0)
            zi = inter(xyi).reshape(xi.shape)

        return zi