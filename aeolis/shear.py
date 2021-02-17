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
import matplotlib
import matplotlib.pyplot as plt
#import scipy.interpolate as spint
#import scipy.spatial.qhull as qhull
#import time

# package modules
from aeolis.utils import *


# initialize logger
logger = logging.getLogger(__name__)


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
    * Separation bubble is still to be improved                                                                 

    '''

    
    igrid = {}
    cgrid = {}
    istransect = False
    
    
    def __init__(self, x, y, z, dx, dy, L, l, z0,
                 buffer_width=100., buffer_relaxation=None):
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
        l : float, optional
            Height of inner layer (default: 10)
        z0 : float, optional
            Aerodynamic roughness (default: .001)

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
        self.l = l
        self.z0 = z0
                          
        self.set_computational_grid()
        

    def __call__(self, u0, udir, process_separation, c, mu_b):
        '''Compute wind shear for given wind speed and direction
        
        Parameters
        ----------
        u0 : float
            Free-flow wind speed
        udir : float
            Wind direction in degrees
        process_separattion : 
        
        '''
        gc = self.cgrid # computational grid
        gi = self.igrid # initial grid
        
        
        # Populate computational grid (rotate to wind direction + interpolate input topography)
        self.populate_computational_grid(udir+90.)
        
        # Compute separation bubble
        if process_separation:
            zsep = self.separation(c, mu_b)
            z_origin = gc['z'].copy()
            gc['z'] = np.maximum(gc['z'], zsep)
                    
        # Compute wind shear stresses on computational grid 
        
        self.compute_shear(u0)
                
        gc['dtaux'], gc['dtauy'] = self.rotate(gc['dtaux'], gc['dtauy'], udir+90)
        
        # Add shear and apply reduction factor for shear in sep. bubble
        self.add_shear()
        
        if process_separation:
            gc['hsep'] = gc['z'] - z_origin
            self.separation_shear(gc['hsep'])
        
        # Rotate both (i&c) grids + results in opposite dir.
        gi['x'], gi['y'] = self.rotate(gi['x'], gi['y'], -(udir+90.), origin=(self.x0, self.y0))
        
        gc['x'], gc['y'] = self.rotate(gc['x'], gc['y'], -(udir+90.), origin=(self.x0, self.y0))
        
        gc['taux'], gc['tauy'] = self.rotate(gc['taux'], gc['tauy'], -(udir+90))
        
    
        # Interpolate wind shear results to real grid
        gi['taux'] = self.interpolate(gc['x'], gc['y'], gc['taux'],
                                              gi['x'], gi['y'])
        gi['tauy'] = self.interpolate(gc['x'], gc['y'], gc['tauy'],
                                              gi['x'], gi['y'])
                
        if process_separation:
            gi['hsep'] = self.interpolate(gc['x'], gc['y'], gc['hsep'], 
                                          gi['x'], gi['y'] )
            
        # Rotate real grid and wind shear results back to orignal orientation
        gc['x'], gc['y'] = self.rotate(gc['x'], gc['y'], udir+90., origin=(self.x0, self.y0))
        gi['x'], gi['y'] = self.rotate(gi['x'], gi['y'], +(udir+90.), origin=(self.x0, self.y0))
                       
        gi['taux'], gi['tauy'] = self.rotate(gi['taux'], gi['tauy'], +(udir+90))

        
        # PLOTTING! -------------------------------------------------------
        # d = 10
        # plt.figure(figsize=(8,4))
        # plt.pcolormesh(gi['x'], gi['y'], gi['z'], cmap='copper_r')  
        # bar = plt.colorbar()
        # bar.set_label('z [m]')                        
        #plt.quiver(gi['x'][::d, ::d], gi['y'][::d, ::d],
        #           gi['taux'][::d, ::d], gi['tauy'][::d, ::d], color='white')
        # plt.xlabel('x [m]')
        # plt.ylabel('y [m]')
        # plt.title('Shear stresses i.grid including separation bubble')
        # plt.show()
        # -------------------------------------------------        
        
        return self
 
    
    # Input functions for wind.py
    def set_topo(self, z):
        '''Update topography

        Parameters
        ----------
        z : numpy.ndarray
            2D array with topography of input grid

        '''

        self.igrid['z'] = z

        return self
    
    def set_shear(self, taus, taun):
        '''Update shear

        Parameters
        ----------
        tau : numpy.ndarray
            array with wind shear stresses of input grid

        '''
        self.igrid['taux'] = taus
        self.igrid['tauy'] = taun
        
        return self
    
    def get_shear(self):
        '''Returns wind shear perturbation field
        
        Returns
        -------
        taux : numpy.ndarray
            Wind shear perturbation in x-direction
        tauy : numpy.ndarray
            Wind shear perturbation in y-direction
        
        '''

        taux = self.igrid['taux']
        tauy = self.igrid['tauy']
            
        return taux, tauy
        
        
    def add_shear(self):
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
        taux = self.cgrid['taux']
        tauy = self.cgrid['tauy']
        
        tau = np.sqrt(taux**2 + tauy**2)
        ix = tau != 0.

        dtaux = self.cgrid['dtaux']
        dtauy = self.cgrid['dtauy']
        
        self.cgrid['taux'][ix] = tau[ix] * (taux[ix] / tau[ix] + dtaux[ix])
        self.cgrid['tauy'][ix] = tau[ix] * (tauy[ix] / tau[ix] + dtauy[ix])
        
        return self

    
    def get_separation(self):
        '''Returns difference in height between z-coordinate of 
        the separation polynomial and of the bed level 
        
        Returns
        -------
        hsep : numpy.ndarray
            Height of seperation bubble
            
        '''  
        hsep = self.igrid['hsep']
        
        return hsep
    
    
    # Input functions for __call()
    def set_computational_grid(self):
        '''Define computational grid
        
        The computational grid is square with dimensions equal to the
        diagonal of the bounding box of the input grid, plus twice the
        buffer width. 

        '''
            
        gi = self.igrid
        gc = self.cgrid
                
        # grid center
        x0, y0 = np.mean(gi['x']), np.mean(gi['y'])
                    
        # grid size
        self.D = np.sqrt((gi['x'].max() - gi['x'].min())**2 +
                         (gi['y'].max() - gi['y'].min())**2) + 2 * self.buffer_width
                        
        # determine equidistant, square grid
        xc, yc = self.get_exact_grid(x0 - self.D/2., x0 + self.D/2.,
                                     y0 - self.D/2., y0 + self.D/2.,
                                     gc['dx'], gc['dy'])
        
        self.x0 = x0
        self.y0 = y0
        gc['xi'] = xc
        gc['yi'] = yc
        
        return self
    
        
    def populate_computational_grid(self, alpha):                               
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
        gi = self.igrid
        gc = self.cgrid
        x = gi['x']
        shp = x.shape
        try:
            ny = shp[2]
        except:
            ny = 0

        # Add buffer zone around grid                                           # buffer is based on version bart, sigmoid function is no longer required
        if ny <= 0:
            dxi = gi['x'][0,0]
            dyi = gi['y'][0,0]
        else:
            dxi = gi['x'][1,1] - gi['x'][0,0]
            dyi = gi['y'][1,1] - gi['y'][0,0]

        buf = 200                                                                # amount of cells


        xi, yi = np.meshgrid(np.linspace(gi['x'][0,0]-buf*dxi, gi['x'][-1,-1]+buf*dxi, gi['x'].shape[1]+2*buf),
                            np.linspace(gi['y'][0,0]-buf*dyi, gi['y'][-1,-1]+buf*dyi, gi['y'].shape[0]+2*buf))
        
        
        zi = np.zeros((xi.shape))
        zi[buf:-buf, buf:-buf] = gi['z']
        
        # Filling buffer zone edges
        zi[buf:-buf,:buf] = np.repeat(zi[buf:-buf,buf+1][:,np.newaxis], buf, axis = 1)
        zi[buf:-buf,-buf:] = np.repeat(zi[buf:-buf,-buf-1][:,np.newaxis], buf, axis = 1)

        zi[:buf,buf:-buf] = np.repeat(zi[buf+1,buf:-buf][np.newaxis], buf, axis = 0)
        zi[-buf:,buf:-buf] = np.repeat(zi[-buf-1,buf:-buf][np.newaxis], buf, axis = 0)
        
        # Filling buffer zone corners
        zi[:buf,:buf] = zi[buf+1,buf+1]
        zi[-buf:,:buf] = zi[-buf-1,buf+1]
        zi[:buf,-buf:] = zi[buf+1,-buf-1]
        zi[-buf:,-buf:] = zi[-buf-1,-buf-1]
        
        # Rotate computational grid to the current wind direction
        xc, yc = self.rotate(gc['xi'], gc['yi'], alpha, origin=(self.x0, self.y0))
        
        # Interpolate input topography to computational grid
        zc = self.interpolate(xi, yi, zi, xc, yc)
        
        # Interpolate input wind - shear
        tauxc = self.interpolate(gi['x'], gi['y'], gi['taux'], xc, yc)
        tauyc = self.interpolate(gi['x'], gi['y'], gi['tauy'], xc, yc)
        
        gc['x'] = xc
        gc['y'] = yc
        gc['z'] = zc
        
        gc['taux'] = tauxc
        gc['tauy'] = tauyc

        return self
    
    def separation(self, c, mu_b):
        
        # Initialize grid and bed dimensions
        
        gc = self.cgrid
         
        x = gc['x']
        y = gc['y']
        z = gc['z']
                
        nx = len(gc['z'][1])
        ny = len(gc['z'][0])
        dx = gc['dx']
        dy = gc['dy']
    
        # Initialize arrays

        dzx = np.zeros(gc['z'].shape)  

        dzdx0 = np.zeros(gc['z'].shape)
        dzdx1 = np.zeros(gc['z'].shape)
                
        stall = np.zeros(gc['z'].shape)
        bubble = np.zeros(gc['z'].shape)
        
        k = np.array(range(0, nx))
        
        zsep =  z.copy()                                                        # total separation bubble      
        
        zsep0 = np.zeros(z.shape)                                               # zero-order separation bubble surface      
        zsep1 = np.zeros(z.shape)                                               # first-oder separation bubble surface
                
        zfft = np.zeros((ny,nx), dtype=np.complex)

        # Compute bed slope angle in x-dir
        dzx[:,:-1] = np.rad2deg(np.arctan((z[:,1:]-z[:,:-1])/dx))
        dzx[:,0] = dzx[:,1]
        dzx[:,-1] = dzx[:,-2]
              
        # Determine location of separation bubbles
        '''Separation bubble exist if bed slope angle (lee side) 
        is larger than max angle that wind stream lines can 
        follow behind an obstacle (mu_b = ..)'''

        stall += np.logical_and(abs(dzx) > mu_b, dzx < 0.) 
        
        
        stall[:,1:-1] += np.logical_and(stall[:,1:-1]==0, stall[:,:-2]>0., stall[:,2:]>0.)
        
        # Define separation bubble
        bubble[:,:-1] = np.logical_and(stall[:,:-1] == 0., stall[:,1:] > 0.) 
        
        # Shift bubble back to x0: start of separation bubble 
        p = 2
        bubble[:,:-p] = bubble[:,p:]
        bubble[:,:p] = 0
        
        bubble = bubble.astype(int)
        
        # Count separation bubbles
        n = np.sum(bubble)
        #print ('number of sep bubbles', n)
        bubble_n = np.asarray(np.where(bubble == True)).T

        
        # Walk through all separation bubbles and determine polynoms
        
        for k in range(0, n):
            
            i = bubble_n[k,1]
            j = bubble_n[k,0]       

            ix_neg = (dzx[j, i+5:] >= 0)                                         # i + 5??
                                    
            if np.sum(ix_neg) == 0:
                zbrink = z[j,i]                                                 # z level of brink at z(x0) 
            else:
                zbrink = z[j,i] - z[j,i+5+np.where(ix_neg)[0][0]]

            # Zero order polynom
            dzdx0 = (z[j,i-1] - z[j,i-2])/dx
            
            #if dzdx0 > 0.1:
            #    dzdx0 = 0.1
            
            a = dzdx0 / c
        
            ls = np.minimum(np.maximum((3.*zbrink/(2.*c) * (1. + a/4. + a**2/8.)), 0.1), 200.)
            # print ('ls:', ls)
            
            a2 = -3 * zbrink/ls**2 - 2 * dzdx0 / ls
            # print('a2:', a2)
            a3 =  2 * zbrink/ls**3 +     dzdx0 / ls**2
            # print('a3:', a3)
          
            i_max = min(i+int(ls/dx),int(nx-1))
            # print('imax:', i_max)

            xs = x[j,i:i_max] - x[j,i]
            
            zsep0[j,i:i_max] = (a3*xs**3 + a2*xs**2 + dzdx0*xs + z[j,i])

            # Zero order filter
            Cut = 1.5
            dk = 2.0 * np.pi / (np.max(x))
            zfft[j,:] = np.fft.fft(zsep0[j,:])
            zfft[j,:] *= np.exp(-(dk*k*dx)**2/(2.*Cut**2))
            zsep0[j,:] = np.real(np.fft.ifft(zfft[j,:]))
            
            # First order polynom
            dzdx1 = (zsep0[j,i-1] - zsep0[j,i-2])/dx
                
            a = dzdx1 / c
        
            ls = np.minimum(np.maximum((3.*z[j,i]/(2.*c) * (1. + a/4. + a**2/8.)), 0.1), 200.)
            # print ('ls:', ls)
            
            a2 = -3 * z[j,i]/ls**2 - 2 * dzdx1 / ls
            # print('a2:', a2)
            a3 =  2 * z[j,i]/ls**3 +     dzdx1 / ls**2
            # print('a3:', a3)
          
            i_max1 = min(i+int(ls/dx),int(nx-1))
            # print('imax:', i_max)

            xs1 = x[j,i:i_max1] - x[j,i]
            
            zsep1[j,i:i_max1] = (a3*xs1**3 + a2*xs1**2 + dzdx1*xs1 + z[j,i])
                                                          
            #plt.figure()
            #plt.plot(x[j,i:i_max], zsep0[j,i:i_max], label='zsep0')
            #plt.plot(x[j,i:i_max], zsep1[j,i:i_max], label='zsep1')
            #plt.plot(x[j,:], z[j,:], 'grey', label='zb')
            #plt.xlabel('x [m]')
            #plt.ylabel('z [m]')
            #plt.title('Separation bubble')
            #plt.legend()
            #plt.show()     
            
            zsep[j,i:i_max] = np.maximum(zsep0[j,i:i_max], z[j,i:i_max])
                          
        
        # Smooth surface of separation bubbles over y direction
        zsep = ndimage.gaussian_filter1d(zsep, sigma=0.2, axis=0)
        
        #plt.pcolormesh(x, y, np.maximum(zsep,z))#, vmin=0, vmax=0.00001)
        #bar = plt.colorbar()
        #bar.set_label('z [m]')
        #plt.title('Zsep gaussian filter in y')
        #plt.show()
        
            
        return zsep
                
    
    def compute_shear(self, u0, nfilter=(1.5,6.)):                               
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
        gc = self.cgrid
        
        if u0 == 0.:
            self.cgrid['dtaux'] = np.zeros(gc['z'].shape)
            self.cgrid['dtauy'] = np.zeros(gc['z'].shape)
            return
                                
        ny, nx = gc['z'].shape
        kx, ky = np.meshgrid(2. * np.pi * np.fft.fftfreq(nx+1, gc['dx'])[1:],
                             2. * np.pi * np.fft.fftfreq(ny+1, gc['dy'])[1:])
        
        hs = np.fft.fft2(gc['z'])
        hs = self.filter_highfrequenies(kx, ky, hs, nfilter, p=0.001)
        
        z0 = self.z0            # roughness length which takes into account saltation
        L  = self.L /4.         # typical length scale of the hill (=1/kx) ??
        
        # Inner layer height
        l  = self.l         
        
        for i in range(5):
            l = 2 * 0.41**2 * L /np.log(l/z0)
        
        # Middle layer height
        hm = 1.0
        for i in range(5):
            hm = L / np.sqrt(np.log(hm/z0))
            
        # Non-dimensional velocity    
        ul = np.log(l/z0) / np.log(hm/z0)
        
        # Arrays in Fourier 
        k = np.sqrt(kx**2 + ky**2)
        sigma = np.sqrt(1j * L * kx * z0 /l)
        
        # Shear stress perturbation
        dtaux_t = hs * kx**2 / k * 2 / ul**2 * \
                  (-1. + (2. * np.log(l/z0) + k**2/kx**2) * sigma * \
                   scipy.special.kv(1., 2. * sigma) / scipy.special.kv(0., 2. * sigma))
        
        dtauy_t = hs * kx * ky / k * 2 / ul**2 * \
                  2. * np.sqrt(2.) * sigma * scipy.special.kv(1., 2. * np.sqrt(2.) * sigma)
        
        gc['dtaux'] = np.real(np.fft.ifft2(dtaux_t))
        gc['dtauy'] = np.real(np.fft.ifft2(dtauy_t))
        
        
    def separation_shear(self, hsep):
        '''Reduces the computed wind shear perturbation below the 
        separation surface to mimic the turbulence effects in the 
        separation bubble
        
        Parameters
        ----------
        hsep : numpy.ndarray
            Height of seperation bubble (in x direction)

        '''
        
        tau_sep = 0.5 
        slope = 0.2                                                             # according to Durán 2010 (Sauermann 2001: c = 0.25 for 14 degrees)
        delta = 1./(slope*tau_sep)
        
        zsepdelta = np.minimum(np.maximum(1. - delta * hsep, 0.), 1.)
        
        self.cgrid['taux'] *= zsepdelta
        self.cgrid['tauy'] *= zsepdelta
        
        
    
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
            fig, ax = plt.subplots()
        
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
        '''Interpolate one grid to an other'''
        
        xy = np.concatenate((y.reshape((-1,1)),
                             x.reshape((-1,1))), axis=1)

        xyi = np.concatenate((yi.reshape((-1,1)),
                              xi.reshape((-1,1))), axis=1)

        if self.istransect:
            zi = np.interp(xi.flatten(), x.flatten(), z.flatten()).reshape(xi.shape)
        else:
            # Trunk interpolation
            # zi = scipy.interpolate.griddata(xy, z.reshape((-1,1)), xyi, method='cubic').reshape(xi.shape)
            
            # version Bart
            inter = scipy.interpolate.RegularGridInterpolator((y[:,0], x[0,:]), z, bounds_error = False, fill_value = 0.)
            zi = inter(xyi).reshape(xi.shape)
            
        return zi
    