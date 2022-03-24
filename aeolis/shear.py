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
import matplotlib.pyplot as plt
#import scipy.interpolate as spint
#import scipy.spatial.qhull as qhull
import time

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
        
        # Assigning values to original (i) and computational (c) grid
        self.igrid = dict(x = x, y = y, z = z)
        self.cgrid = dict(dx = dx, dy = dy)
        
        # Setting buffer settings                  
        self.buffer_width = buffer_width
        self.buffer_relaxation = buffer_relaxation
        
        # Setting shear perturbation settings                  
        self.L = L
        self.l = l
        self.z0 = z0
       

    def __call__(self, u0, udir, process_separation, c, mu_b, taus0, taun0, zero_order_filter, plot=True):
        '''Compute wind shear for given wind speed and direction
        
        Parameters
        ----------
        u0 : float
            Free-flow wind speed
        udir : float
            Wind direction in degrees
        process_separattion : 
        
        '''
        
        # Convert to cartesian to perform all the rotations
        u_angle = 270. - udir # wind angle
        
        if plot:
            fig, axs = plt.subplots(2, 3)
            self.plot(ax=axs[0,0], cmap='Reds', stride=10, computational_grid=False)
            axs[0,0].set_title('Original input grid')

        # =====================================================================
        # Creating, rotating and filling computational grid
        # =====================================================================
         
        # Creating the computational grid         
        self.set_computational_grid(udir)
        
        # Storing computational (c) and original (i) grids
        gi = self.igrid # initial grid
        gc = self.cgrid # computational grid
           
        # Rotate computational (c) grid to the current wind direction
        gc['x'], gc['y'] = self.rotate(gc['xi'], gc['yi'], -u_angle, 
                                       origin=(self.x0, self.y0))
        
        # =====================================================================
        # Filling the computational grid with bedlevel and shear stress
        # =====================================================================
        
        # For now turned off because caused problems.
        # Just normal extrapolation applied
        # xi_buff, yi_buff, zi_buff = self.buffer_original_grid()
        # gc['z'] = self.interpolate(xi_buff, yi_buff, zi_buff, gc['x'], gc['y'], 0)
        
        # Interpolate bed levels and shear to the computational grid
        gc['z'] = self.interpolate(gi['x'], gi['y'], gi['z'], gc['x'], gc['y'], None)
        
        # Project the taus0 and taun0 on the computational grid
        gc['taux'] = np.full(np.shape(gc['x']), taus0)
        gc['tauy'] = np.full(np.shape(gc['x']), taun0)

        if plot:
            self.plot(ax=axs[0,1], cmap='Reds', stride=10, computational_grid=True)
            axs[0,1].set_title('Interpolated values on computational grid')

        # =====================================================================
        # Computing  bubble and add it the bedlevel for  shear perturbation.
        # Afterwards, computing the change in shear stress (dtaux and dtauy), 
        # rotate to horizontal computational grid and add to tau0 
        # =====================================================================

        # Compute separation bubble
        if process_separation:
            zsep = self.separation(c, mu_b, zero_order_filter)
            z_origin = gc['z'].copy()
            gc['z'] = np.maximum(gc['z'], zsep)
                    
        # Compute wind shear stresses on computational grid 
        self.compute_shear(u0)

        # Rotate to the horizontal rotational grid        
        gc['dtaux'], gc['dtauy'] = self.rotate(gc['dtaux'], gc['dtauy'], -u_angle)
        
        # Add shear and apply reduction factor for shear in sep. bubble
        self.add_shear()

        # Compute the influence of the separation on the shear stress
        if process_separation:
            gc['hsep'] = gc['z'] - z_origin
            self.separation_shear(gc['hsep'])
            
        if plot:
            self.plot(ax=axs[0,2], cmap='Reds', stride=10, computational_grid=True)
            axs[0,2].set_title('Computed shear stress and bubble effect')

        if plot:
            pc = axs[1,0].pcolormesh(gc['x'], gc['y'], gc['taux'])
            plt.colorbar(pc, ax=axs[1,0])
            axs[1,0].set_title('Rotate grids, such that computational is horizontal')
        
        # =====================================================================
        # Interpolation from the computational grid back to the original
        # =====================================================================
    
        # Interpolate wind shear results to real grid
        gi['taux'] = self.interpolate(gc['x'], gc['y'], gc['taux'], gi['x'], gi['y'], taus0)
        gi['tauy'] = self.interpolate(gc['x'], gc['y'], gc['tauy'], gi['x'], gi['y'], taun0)
  
        if process_separation:
            gi['hsep'] = self.interpolate(gc['x'], gc['y'], gc['hsep'], gi['x'], gi['y'], 0. )
            
        # Final plots and lay-out    
        if plot:
            pc = axs[1,1].pcolormesh(gi['x'], gi['y'], gi['taux'])
            plt.colorbar(pc, ax=axs[1,1])
            axs[1,1].set_title('Interpolate back onto original grid')

            self.plot(ax=axs[1,2], cmap='Reds', stride=10, computational_grid=False)
            axs[1,2].set_title('Rotate original grid back')

            for axr in axs:
                for ax in axr:
                    ax.set_xlim([-400, 400])
                    ax.set_ylim([-400, 400])
                    ax.set_aspect('equal')
        
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
        b_W =       np.zeros(4)
        b_L =       np.zeros(4)
        xcorner =   np.zeros(4)
        ycorner =   np.zeros(4)
        
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
        udir_verticals = np.arange(-360, 360, 90)  
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
            m_L, b_L[i] = np.polyfit([xcorner[i], xcorner[i] - np.sin(np.deg2rad(udir-90.))],
                                     [ycorner[i], ycorner[i] - np.cos(np.deg2rad(udir-90.))], 1)
 
        # Determine the most outer boundaries (for parallel and perpendicular)
        db_W = self.maxDiff(b_W)
        db_L = self.maxDiff(b_L)
        
        # Compute the distance between the outer boundaries to determine the width (W) and length (L) of the grid
        self.W = abs(db_W) / np.sqrt((m_W**2.) + 1) + self.buffer_width * 2.
        self.L = abs(db_L) / np.sqrt((m_L**2.) + 1) + self.buffer_width * 2.
        
        # Create the grid
        xc, yc = self.get_exact_grid(x0 - self.L/2., x0 + self.L/2.,
                                     y0 - self.W/2., y0 + self.W/2.,
                                     gc['dx'], gc['dy'])
        
        # Storing grid parameters
        self.x0 = x0
        self.y0 = y0
        gc['xi'] = xc
        gc['yi'] = yc
        
        return self
    
        
    
    def separation(self, c, mu_b, zero_order_filter):
        
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
        bubble[:,:-1] = (stall[:,:-1] == 0.) * (stall[:,1:] > 0.) 
        
        # Shift bubble back to x0: start of separation bubble 
        p = 1
        bubble[:,p:] = bubble[:,:-p]
        bubble[:,-p:] = 0
        
        bubble = bubble.astype(int)
        
        # Count separation bubbles
        n = np.sum(bubble)
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
            dzdx0 = (z[j,i] - z[j,i-1]) / dx
            # dzdx0 = 1.
            
            a = dzdx0 / c
        
            ls = np.minimum(np.maximum((3.*zbrink/(2.*c) * (1. + a/4. + a**2/8.)), 0.1), 200.)
            
            a2 = -3 * zbrink/ls**2 - 2 * dzdx0 / ls
            a3 =  2 * zbrink/ls**3 +     dzdx0 / ls**2
          
            i_max = min(i+int(ls/dx),int(nx-1))

            xs = x[j,i:i_max] - x[j,i]
            
            zsep0[j,i:i_max] = (a3*xs**3 + a2*xs**2 + dzdx0*xs + z[j,i])

            # Zero order filter
            if zero_order_filter:
                
                Cut = 1.5
                dk = 2.0 * np.pi / (np.max(x))
                zfft[j,:] = np.fft.fft(zsep0[j,:])
                zfft[j,:] *= np.exp(-(dk*k*dx)**2/(2.*Cut**2))
                zsep0[j,:] = np.real(np.fft.ifft(zfft[j,:]))
                
                # First order polynom
                dzdx1 = (zsep0[j,i] - zsep0[j,i-1])/dx
                   
                a = dzdx1 / c
            
                ls = np.minimum(np.maximum((3.*z[j,i]/(2.*c) * (1. + a/4. + a**2/8.)), 0.1), 200.)
                # print ('ls:', ls)
                
                a2 = -3 * z[j,i]/ls**2 - 2 * dzdx1 / ls
                a3 =  2 * z[j,i]/ls**3 +     dzdx1 / ls**2
              
                i_max1 = min(i+int(ls/dx),int(nx-1))
    
                xs1 = x[j,i:i_max1] - x[j,i]
                
                zsep1[j,i:i_max1] = (a3*xs1**3 + a2*xs1**2 + dzdx1*xs1 + z[j,i])
            

            # Pick the maximum seperation bubble hieght at all locations
            if zero_order_filter:
                zsep[j,i:i_max] = np.maximum(zsep1[j,i:i_max], zsep[j,i:i_max])
            else:
                zsep[j,i:i_max] = np.maximum(zsep0[j,i:i_max], zsep[j,i:i_max])

            
        # Smooth surface of separation bubbles over y direction
        zsep = ndimage.gaussian_filter1d(zsep, sigma=0.2, axis=0)

        #Correct for any seperation bubbles that are below the bed surface following smoothing
        ilow = zsep < z
        zsep[ilow] = z[ilow]

            
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
        
        
        time_start_perturbation = time.time()
        
        # Shear stress perturbation
        
        dtaux_t = hs * kx**2 / k * 2 / ul**2 * \
                  (-1. + (2. * np.log(l/z0) + k**2/kx**2) * sigma * \
                    sc_kv(1., 2. * sigma) / sc_kv(0., 2. * sigma))

        
        dtauy_t = hs * kx * ky / k * 2 / ul**2 * \
                    2. * np.sqrt(2.) * sigma * sc_kv(1., 2. * np.sqrt(2.) * sigma)

        
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
        
        
        
    def maxDiff(self, arr):

        result = 0
        n = len(arr)
     
        # Iterate through all pairs.
        for i in range(0,n):
            for j in range(i, n):
     
                if (abs(arr[i] - arr[j]) + abs(i - j) > result):
                    result = abs(arr[i] - arr[j]) + abs(i - j)
             
        return result
        
        
    
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
                  g['taux'][::d,::d], g['tauy'][::d,::d], **kwargs)
                  
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
        
    
    def interpolate(self, x, y, z, xi, yi, z0):
        '''Interpolate one grid to an other'''

        
        # First compute angle with horizontal
        dx = x[0,1] - x[0,0]
        dy = y[0,1] - y[0,0]
        
        if dx != 0:
            dydx = dy/dx
        elif dy >= 0:
            dydx = np.inf
        elif dy < 0:
            dydx = -np.inf

        angle = np.rad2deg(np.arctan(dy/dx))
        
        if dx <= 0 and dy<=0:
            angle += 180.
            
        # Rotate grids to allign with horizontal
        x, y = self.rotate(x, y, angle, origin=(self.x0, self.y0))
        xi, yi = self.rotate(xi, yi, angle, origin=(self.x0, self.y0))
        
        # Rotate 180 deg if necessary
        if not np.all(sorted(y[:,0]) == y[:,0]) and not np.all(sorted(x[0,:]) == x[0,:]):
            x, y = self.rotate(x, y, 180, origin=(self.x0, self.y0))
            xi, yi = self.rotate(xi, yi, 180, origin=(self.x0, self.y0))

        # Concatenate
        xy = np.concatenate((y.reshape((-1,1)),
                             x.reshape((-1,1))), axis=1)

        xyi = np.concatenate((yi.reshape((-1,1)),
                              xi.reshape((-1,1))), axis=1)

        if self.istransect:
            zi = np.interp(xi.flatten(), x.flatten(), z.flatten()).reshape(xi.shape)
        else:
            inter = scipy.interpolate.RegularGridInterpolator((y[:,0], x[0,:]), z, bounds_error = False, fill_value = z0)
            zi = inter(xyi).reshape(xi.shape)
            
        return zi
    

        
        