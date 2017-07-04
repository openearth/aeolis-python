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


import numpy as np
import scipy.special
import scipy.interpolate


class WindShear:
    '''Class for computation of 2DH wind shear perturbations over a topography.
        
    The class implements a 2D FFT solution to the wind shear
    perturbation on curvilinear grids.  As the FFT solution is only
    defined on an equidistant rectilinear grid with circular boundary
    conditions that is aligned with the wind direction, a rotating
    computational grid is automatically defined for the computation.
    The computational grid is extended in all directions using a
    logistic sigmoid function as to ensure full coverage of the input
    grid for all wind directions, circular boundaries and preservation
    of the alongshore uniformity.  An extra buffer distance can be
    used as to minimize the disturbence from the borders in the input
    grid.  The results are interpolated back to the input grid when
    necessary.

    Frequencies related to wave lengths smaller than a computational
    grid cell are filtered from the 2D spectrum of the topography
    using a logistic sigmoid tapering. The filtering aims to minimize
    the disturbance as a result of discontinuities in the topography
    that may physically exists, but cannot be solved for in the
    computational grid used.

    Example
    -------
    >>> w = WindShear(x, y, z)
    >>> w(u0=10., udir=30.).add_shear(taux, tauy)

    Notes
    -----
    To do:

    * Actual resulting values are still to be compared with the results
       from Kroy et al. (2002)
    * Grid interpolation can still be optimized
    * Separation bubble is still to be implemented
    * Avalanching is still to be implemented

    '''

    
    igrid = {}
    cgrid = {}
    istransect = False
    
    
    def __init__(self, x, y, z, dx=1., dy=1.,
                 buffer_width=100., buffer_relaxation=None,
                 L=100., z0=.001, l=10.):
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
        buffer_relaxation : float, optional
            Relaxation of topography in buffer from input grid
            boundary to computational grid boundary (default:
            buffer_width / 4)
        L : float, optional
            Length scale of topographic features (default: 100)
        z0 : float, optional
            Aerodynamic roughness (default: .001)
        l : float, optional
            Height of inner layer (default: 10)

        '''
        
        if buffer_relaxation is None:
            buffer_relaxation = buffer_width / 4.

        if z.shape[0] == 1:
            self.istransect = True
        
        self.igrid = dict(x = x,
                          y = y,
                          z = z)
            
        self.cgrid = dict(dx = dx,
                          dy = dy)
                          
        self.buffer_width = buffer_width
        self.buffer_relaxation = buffer_relaxation
                          
        self.L = L
        self.z0 = z0
        self.l = l
                          
        self.set_computational_grid()


    def __call__(self, u0, udir):
        '''Compute wind shear for given wind speed and direction
        
        Parameters
        ----------
        u0 : float
            Free-flow wind speed
        udir : float
            Wind direction in degrees
        
        '''
            
        self.populate_computational_grid(udir)
        self.compute_shear(u0)
                    
        gc = self.cgrid
        gi = self.igrid
                            
        dtaux, dtauy = self.rotate(gc['dtaux'], gc['dtauy'], udir)
                                
        self.cgrid['dtaux'] = dtaux
        self.cgrid['dtauy'] = dtauy
                                        
        self.igrid['dtaux'] = self.interpolate(gc['x'], gc['y'], dtaux,
                                               gi['x'], gi['y'])
        self.igrid['dtauy'] = self.interpolate(gc['x'], gc['y'], dtauy,
                                               gi['x'], gi['y'])
        
        return self


    def get_shear(self):
        '''Returns wind shear perturbation field
        
        Returns
        -------
        dtaux : numpy.ndarray
            Wind shear perturbation in x-direction
        dtauy : numpy.ndarray
            Wind shear perturbation in y-direction
        
        '''

        dtaux = self.igrid['dtaux']
        dtauy = self.igrid['dtauy']
            
        return dtaux, dtauy
        
        
    def add_shear(self, taux, tauy):
        '''Add wind shear perturbations to a given wind shear field
        
        Parameters
        ----------
        taux : numpy.ndarray
            Wind shear in x-direction
        tauy : numpy.ndarray
            Wind shear in y-direction

        Returns
        -------
        taux : numpy.ndarray
            Wind shear including perturbations in x-direction
        tauy : numpy.ndarray
            Wind shear including perturbations in y-direction
        
        '''

        tau = np.sqrt(taux**2 + tauy**2)
        ix = tau != 0.

        dtaux = self.igrid['dtaux']
        dtauy = self.igrid['dtauy']
        
        taux[ix] = tau[ix] * (taux[ix] / tau[ix] + dtaux[ix])
        tauy[ix] = tau[ix] * (tauy[ix] / tau[ix] + dtauy[ix])

        return taux, tauy


    def set_topo(self, z):
        '''Update topography

        Parameters
        ----------
        z : numpy.ndarray
            2D array with topography of input grid

        '''

        self.igrid['z'] = z

        return self
        

    def populate_computational_grid(self, alpha):
        '''Interpolate input topography to computational grid
            
        Rotates computational grid to current wind direction and
        interpolates the input topography to the rotated grid. Any
        grid cells that are not covered by the input grid are filled
        using a sigmoid function.
            
        Parameters
        ----------
        alpha : float
            Rotation angle in degrees

        '''
        
        gc = self.cgrid
        gi = self.igrid
        
        xc, yc = self.rotate(gc['xi'], gc['yi'], alpha, origin=(self.x0, self.y0))
        zc = self.interpolate(gi['x'], gi['y'], gi['z'], xc, yc)
        self.cgrid['z'] = zc
        self.cgrid['x'] = xc
        self.cgrid['y'] = yc
        
        bx = self.get_borders(gi['x'])
        by = self.get_borders(gi['y'])
        bz = self.get_borders(gi['z'])

        ix = np.isnan(zc)
        if np.any(ix):
            
            d = np.zeros((np.sum(ix),))
            z = np.zeros((np.sum(ix),))
            
            for i, (xn, yn) in enumerate(zip(xc[ix], yc[ix])):
                
                distances = np.hypot(bx - xn, by - yn)
                idx = np.argmin(distances)
                d[i] = np.min(distances)
                z[i] = bz[idx]
                
                for j in range(2):
                    i1 = idx+j-1
                    i2 = idx+j
                    
                    k = self.interpolate_projected_point((bx[i1], by[i1], bz[i1]),
                                                         (bx[i2], by[i2], bz[i2]),
                                                         (xn, yn))
                    
                    if k:
                        d[i] = k[0]
                        z[i] = k[1]
                        break
                    
            self.cgrid['z'][ix] = z * self.get_sigmoid(d)

        
    def compute_shear(self, u0, nfilter=(1.,2.)):
        '''Compute wind shear perturbation for given free-flow wind speed on computational grid
        
        Parameters
        ----------
        u0 : float
            Free-flow wind speed
        nfilter : 2-tuple
            Wavenumber range used for logistic sigmoid filter. See
            :func:`filter_highfrequencies`

        '''
            
        g = self.cgrid
                
        if u0 == 0.:
            self.cgrid['dtaux'] = np.zeros(g['z'].shape)
            self.cgrid['dtauy'] = np.zeros(g['z'].shape)
            return
                                
        ny, nx = g['z'].shape
        kx, ky = np.meshgrid(2. * np.pi * np.fft.fftfreq(nx+1, g['dx'])[1:],
                             2. * np.pi * np.fft.fftfreq(ny+1, g['dy'])[1:])
        hs = np.fft.fft2(g['z'])
        hs = self.filter_highfrequenies(kx, ky, hs, nfilter)

        k = np.sqrt(kx**2 + ky**2)
        sigma = np.sqrt(1j * self.L / 4. * kx * self.z0 / self.l)
        
        dtaux_t = hs * kx**2 / k * 2 / u0**2 * \
                  (-1 + (2 * np.log(self.l/self.z0) + k**2/kx**2) * sigma * \
                   scipy.special.kv(1, 2 * sigma) / scipy.special.kv(0, 2 * sigma))
        dtauy_t = hs * kx * ky / k * 2 / u0**2 * \
                  2 * np.sqrt(2) * sigma * scipy.special.kv(1, 2 * np.sqrt(2) * sigma)
        
        self.cgrid['dtaux'] = np.real(np.fft.ifft2(dtaux_t))
        self.cgrid['dtauy'] = np.real(np.fft.ifft2(dtauy_t))
        
        
    def set_computational_grid(self):
        '''Define computational grid
        
        The computational grid is square with dimensions equal to the
        diagonal of the bounding box of the input grid, plus twice the
        buffer width.

        '''
            
        g = self.igrid
                
        # grid center
        x0, y0 = np.mean(g['x']), np.mean(g['y'])
                    
        # grid size
        self.D = np.sqrt((g['x'].max() - g['x'].min())**2 +
                         (g['y'].max() - g['y'].min())**2) + 2 * self.buffer_width
                        
        # determine equidistant, square grid
        xc, yc = self.get_exact_grid(x0 - self.D/2., x0 + self.D/2.,
                                     y0 - self.D/2., y0 + self.D/2.,
                                     self.cgrid['dx'], self.cgrid['dy'])
        
        self.x0 = x0
        self.y0 = y0
        self.cgrid['xi'] = xc
        self.cgrid['yi'] = yc
        
        
    def get_sigmoid(self, x):
        '''Get sigmoid function value
        
        Get bed level multiplication factor in buffer area based on
        buffer specificationa and distance to input grid boundary.
        
        Parameters
        ----------
        x : float or numpy.ndarray
            Distance(s) to input grid boundary
        
        Returns
        -------
        float or numpy.ndarray
            Bed level multiplication factor (z = factor * z_boundary)

        '''
            
        return 1. / (1. + np.exp(-(self.buffer_width-x) / self.buffer_relaxation))
        

    def filter_highfrequenies(self, kx, ky, hs, nfilter=(1, 2), p=.01):
        '''Filter high frequencies from a 2D spectrum

        A logistic sigmoid filter is used to taper higher frequencies
        from the 2D spectrum. The range over which the sigmoid runs
        from 0 to 1 with a precision ``p`` is given by the 2-tuple
        ``nfilter``. The range is defined as wavenumbers in terms of
        gridcells, i.e. a value 1 corresponds to a wave with length
        ``dx``.

        Parameters
        ----------
        kx : numpy.ndarray
            Wavenumbers in x-direction
        ky : numpy.ndarray
            Wavenumbers in y-direction
        hs : numpy.ndarray
            2D spectrum
        nfilter : 2-tuple
            Wavenumber range used for logistic sigmoid filter
        p : float
            Precision of sigmoid range definition

        Returns
        -------
        hs : numpy.ndarray
            Filtered 2D spectrum

        '''

        if nfilter is not None:
            n1 = np.min(nfilter)
            n2 = np.max(nfilter)
            px = 2 * np.pi / self.cgrid['dx'] / np.abs(kx)
            py = 2 * np.pi / self.cgrid['dy'] / np.abs(ky)
            s1 =  n1 / np.log(1. / .01 - 1.)
            s2 = -n2 / np.log(1. / .99 - 1.)
            f1 = 1. / (1. + np.exp(-(px + n1 - n2) / s1))
            f2 = 1. / (1. + np.exp(-(py + n1 - n2) / s2))
            hs *= f1 * f2

        return hs
    
                                                                                                                                
    def plot(self, ax=None, cmap='Reds', stride=10, computational_grid=False, **kwargs):
        '''Plot wind shear perturbation
            
        Parameters
        ----------
        ax : matplotlib.pyplot.Axes, optional
            Axes to plot onto
        cmap : matplotlib.cm.Colormap or string, optional
            Colormap for topography (default: Reds)
        stride : int, optional
            Stride to apply to wind shear vectors (default: 10)
        computational_grid : bool, optional
            Plot on computational grid rather than input grid
            (default: False)
        kwargs : dict
            Additional arguments to :func:`matplotlib.pyplot.quiver`
            
        Returns
        -------
        ax : matplotlib.pyplot.Axes
            Axes used for plotting

        '''
        
        d = stride
        
        if ax is None:
            fig, ax = subplots()
        
        if computational_grid:
            g = self.cgrid
        else:
            g = self.igrid
        
        ax.pcolormesh(g['x'], g['y'], g['z'], cmap=cmap)
        ax.quiver(g['x'][::d,::d], g['y'][::d,::d], 
                  g['dtaux'][::d,::d], g['dtauy'][::d,::d], **kwargs)
                  
        if computational_grid:
            ax.plot(self.get_borders(self.igrid['x']),
                    self.get_borders(self.igrid['y']), '-k')
                  
        return ax


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
        
        return np.concatenate((x[0,:].T, 
                               x[1:-1,-1], 
                               x[-1,::-1].T, 
                               x[-1:1:-1,0],
                               x[0,:1]), axis=0)
    
    
    @staticmethod
    def rotate(x, y, alpha, origin=(0,0)):
        '''Rotate a matrix over given angle around given origin'''
        
        xr = x - origin[0]
        yr = y - origin[1]
        
        a = alpha / 180. * np.pi
        
        R = np.asmatrix([[np.cos(a), -np.sin(a)],
                         [np.sin(a),  np.cos(a)]])
        
        xy = np.concatenate((xr.reshape((-1,1)), 
                             yr.reshape((-1,1))), axis=1) * R
                         
        return (np.asarray(xy[:,0].reshape(x.shape) + origin[0]),
                np.asarray(xy[:,1].reshape(y.shape) + origin[1]))
    
    
    def interpolate(self, x, y, z, xi, yi):
        '''Interpolate a grid onto another grid'''
        
        xy = np.concatenate((x.reshape((-1,1)),
                             y.reshape((-1,1))), axis=1)
        xyi = np.concatenate((xi.reshape((-1,1)),
                              yi.reshape((-1,1))), axis=1)

        if self.istransect:
            z = np.interp(xi.flatten(), x.flatten(), z.flatten()).reshape(xi.shape)
        else:
            z = scipy.interpolate.griddata(xy, z.reshape((-1,1)), xyi, method='cubic').reshape(xi.shape)
                             
        return z


    @staticmethod
    def interpolate_projected_point(a, b, p):
        '''Project point to line segment and return distance and interpolated value
        
        Parameters
        ----------
        a : iterable
            Start vector for line segment
        b : iterable
            End vector for line segment
        p : iterable
            Point vector to be projected
            
        Returns
        -------
        d : numpy.ndarray
            Distance from point p to projected point q
        z : float
            Interpolated value at projected point q
            
        '''
        
        a = np.asarray(a)
        b = np.asarray(b)
        p = np.asarray(p)
        
        ab = b[:-1]-a[:-1]                     # line segment
        ab2 = np.dot(ab, ab)                   # length of line segment squared
        
        if ab2 > 0.:
            ap = p[:-1]-a[:-1]
            t = np.dot(ap, ab) / ab2           # location of projected point along line segment as fraction of its length
            if t >= 0. and t <= 1.:
                q = a[:-1] + t * ab            # projected point
                pq = p-q                       # vector from original to projected point
                d = np.sqrt(np.dot(pq, pq))    # distance from original to projected point
                z = a[-1] * (1.-t) + b[-1] * t # linearly interpolated height of projected point
                return d, z
            
        return
