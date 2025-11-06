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
import importlib.machinery
import importlib.metadata
import time
import glob
import logging
import warnings
import operator
import numpy as np
import scipy.sparse
import pickle
import scipy.sparse.linalg
import matplotlib.pyplot as plt
from datetime import timedelta
from bmi.api import IBmi
from functools import reduce
from numba import njit


# package modules
import aeolis.inout
import aeolis.bed
import aeolis.avalanching
import aeolis.wind
import aeolis.threshold
import aeolis.transport
import aeolis.hydro
import aeolis.netcdf
import aeolis.constants
import aeolis.erosion
import aeolis.vegetation
import aeolis.fences
import aeolis.gridparams

# type hints
from typing import Any, Union, Tuple

from aeolis.utils import *


class StreamFormatter(logging.Formatter):
    """A formater for log messages"""

    def format(self, record) -> str:
        """Formats log messages for console standard output"""

        if record.levelname == 'INFO':
            return record.getMessage()
        else:
            return '%s: %s' % (record.levelname, record.getMessage())


# initialize logger
logger = logging.getLogger(__name__)

__version__ = importlib.metadata.version('aeolis')

if not __version__ :
    logger.warning('WARNING: Unknown model version.')


class ModelState(dict):
    """Dictionary-like object to store model state.
    Model state variables are mutable by default, but can be set
    to immutable. In the latter case any actions that set the immutable
    model state variable are ignored.
    """

    def __init__(self, *args, **kwargs) -> None:
        """Initialize dictionary for model states"""
        self.ismutable = set()
        super(ModelState, self).__init__(*args, **kwargs)


    def __setitem__(self, k:str, v:Any) -> None:
        """sets new and existing mutable model variables and thier values 
        to the model state
        """
        
        if k not in self.keys() or k in self.ismutable:
            super(ModelState, self).__setitem__(k, v)
            self.set_mutable(k)


    def set_mutable(self, k:str) -> None:
        """Adds mutable variable to the model state"""
        self.ismutable.add(k)


    def set_immutable(self, k:str) -> None:
        """Removes existing model varialble from the set of mutable variables"""
        if k in self.ismutable:
            self.ismutable.remove(k)


class AeoLiS(IBmi):
    """AeoLiS model class

    AeoLiS is a process-based model for simulating supply-limited
    aeolian sediment transport. This model class is compatible with
    the Basic Model Interface (BMI) and provides basic model
    operations, like initialization, time stepping, finalization and
    data exchange. For higher level operations, like a progress
    indicator and netCDF4 output is refered to the AeoLiS model
    runner class, see :class:`~model.AeoLiSRunner`.

    Examples
    --------
    >>> with AeoLiS(configfile='aeolis.txt') as model:
    >>>     while model.get_current_time() <= model.get_end_time():
    >>>         model.update()

    >>> model = AeoLiS(configfile='aeolis.txt')
    >>> model.initialize()
    >>> zb = model.get_var('zb')
    >>> model.set_var('zb', zb + 1)
    >>> for i in range(10):
    >>>     model.update(60.) # step 60 seconds forward
    >>> model.finalize()
    """

    def __init__(self, configfile: str) -> None:
        '''Initialize class

        Parameters
        ----------
        configfile:
            Path to model configuration file. See :func:`~inout.read_configfile()`.

        '''

        self.t: float = 0
        self.dt: float = 0
        self.configfile: str = ''

        self.l: dict = {} # previous spatial grids
        self.s: ModelState = ModelState() # spatial grids
        self.p: dict = {} # parameters
        self.c: dict = {} # counters

        self.configfile = configfile

    def __enter__(self):
        """Initialize AeoliS class"""
        self.initialize()
        return self

    def __exit__(self, *args) -> None:
        """Calls to terminate AeoLiS class"""
        self.finalize()

    def initialize(self)-> None:
        '''
        Read model configuration file, initialize parameters and
        spatial grids dictionary, and load bathymetry and bed
        composition.
        '''

        # read configuration file
        self.p = aeolis.inout.read_configfile(self.configfile)
        aeolis.inout.check_configuration(self.p)

        # set nx, ny and nfractions and make 1D into quasi 2D
        if self.p['xgrid_file'].ndim == 2:
            self.p['ny'], self.p['nx'] = self.p['xgrid_file'].shape
            
            # change from number of points to number of cells
            self.p['nx'] -= 1  
            self.p['ny'] -= 1
            
        else:
            if 0:
                #this is the old 1D stuff
                self.p['nx'] = len(self.p['xgrid_file'])
                self.p['nx'] -= 1 
                self.p['ny'] = 0

            if 1:
                # this is the new quasi 2D stuff
                # this is where we make the 1D grid into a 2D grid to ensure process compatibility
                
                # define size of 2D grid 3 is minumum, variable defined for debugging purposes
                qnr = 3
                
                self.p['xgrid_file'] = np.transpose(np.stack([self.p['xgrid_file'] for i in range (qnr)],axis=1))
                
                # redefine shape
                self.p['ny'], self.p['nx'] = self.p['xgrid_file'].shape

                # repeat the above for ygrid_file assume dy is equal to dx
                dy = self.p['xgrid_file'][1,2]-self.p['xgrid_file'][1,1]
                self.p['ygrid_file'] = np.transpose(np.stack([self.p['ygrid_file'] + i * dy for i in range(qnr)], axis=1))
                
                # repeat the above for bed_file
                self.p['bed_file'] = np.transpose(np.stack([self.p['bed_file'] for i in range (qnr)],axis=1))

                # change from number of points to number of cells
                self.p['nx'] -= 1  
                self.p['ny'] -= 1

        #self.p['nfractions'] = len(self.p['grain_dist'])
        self.p['nfractions'] = len(self.p['grain_size'])
		
        # initialize time
        self.t = self.p['tstart']

        # get model dimensions
        nx = self.p['nx']
        ny = self.p['ny']
        nl = self.p['nlayers']
        nf = self.p['nfractions']

        # initialize spatial grids
        for var, dims in self.dimensions().items():
            self.s[var] = np.zeros(self._dims2shape(dims))
            self.l[var] = self.s[var].copy()
            
        # initialize grid parameters
        self.s, self.p = aeolis.gridparams.initialize(self.s, self.p)

        # initialize bed composition
        self.s = aeolis.bed.initialize(self.s, self.p)

        # initialize wind model
        self.s = aeolis.wind.initialize(self.s, self.p)
         
        #initialize vegetation model
        self.s = aeolis.vegetation.initialize(self.s, self.p)                  

        #initialize fence model
        self.s = aeolis.fences.initialize(self.s, self.p)

        # Create interpretation information
        if self.p['visualization']:
            aeolis.inout.visualize_grid(self.s, self.p)
            aeolis.inout.visualize_timeseries(self.p, self.t)


    def update(self, dt:float=-1) -> None:
        '''Time stepping function

        Takes a single step in time. Interpolates wind and
        hydrodynamic time series to the current time, updates the soil
        moisture, mixes the bed due to wave action, computes wind
        velocity threshold and the equilibrium sediment transport
        concentration. Subsequently runs one of the available
        numerical schemes to compute the instantaneous sediment
        concentration and pickup for the next time step and updates
        the bed accordingly.

        For explicit schemes the time step is maximized by the
        Courant-Friedrichs-Lewy (CFL) condition. See
        :func:`~model.AeoLiS.set_timestep()`.

        Parameters
        ----------
        dt : optional
            Time step in seconds. The time step specified in the model
            configuration file is used in case dt is smaller than
            zero. For explicit numerical schemes the time step is
            maximized by the CFL confition.

        '''

        self.p['_time'] = self.t 
  

        
        # store previous state
        self.l = self.s.copy()
        self.l['zb'] = self.s['zb'].copy()
        self.l['dzbavg'] = self.s['dzbavg'].copy()
        
        # interpolate wind time series
        self.s = aeolis.wind.interpolate(self.s, self.p, self.t)  
    
        # Rotate gridparams, such that the grids is alligned horizontally
        self.s = self.grid_rotate(self.p['alpha'])
      
        if np.sum(self.s['uw']) != 0:
            self.s = aeolis.wind.shear(self.s, self.p)

        #compute sand fence shear
        if self.p['process_fences']:
            self.s = aeolis.fences.update_fences(self.s, self.p)

        # compute vegetation shear
        if self.p['process_vegetation']: 
            self.s = aeolis.vegetation.vegshear(self.s, self.p)
        
        # determine optimal time step
        self.dt_prev = self.dt
        if not self.set_timestep(dt):
            return

        # interpolate hydrodynamic time series
        self.s = aeolis.hydro.interpolate(self.s, self.p, self.t)
        self.s = aeolis.hydro.update(self.s, self.p, self.dt, self.t)

        # mix top layer
        self.s = aeolis.bed.mixtoplayer(self.s, self.p)
        
        # compute threshold
        # if np.sum(self.s['uw']) != 0:
        self.s = aeolis.threshold.compute(self.s, self.p)

        # compute saltation velocity and equilibrium transport
        self.s = aeolis.transport.equilibrium(self.s, self.p)

        # compute instantaneous transport
        if self.p['scheme'] == 'euler_forward':
            self.s.update(self.euler_forward())
        elif self.p['scheme'] == 'euler_backward':
            self.s.update(self.euler_backward())
        elif self.p['scheme'] == 'crank_nicolson':
            self.s.update(self.crank_nicolson())
        else:
            logger.log_and_raise('Unknown scheme [%s]' % self.p['scheme'], exc=ValueError)

        # update bed
        self.s = aeolis.bed.update(self.s, self.p)

        # avalanching
        self.s = aeolis.avalanching.angele_of_repose(self.s, self.p)
        self.s = aeolis.avalanching.avalanche(self.s, self.p)
        
        # reset original bed in marine zone (wet)
        self.s = aeolis.bed.wet_bed_reset(self.s, self.p)

        # calculate average bed level change over time
        self.s = aeolis.bed.average_change(self.l, self.s, self.p)

        # compute dune erosion
        if self.p['process_dune_erosion']:
            self.s = aeolis.erosion.run_ph12(self.s, self.p, self.t)
            self.s = aeolis.avalanching.angele_of_repose(self.s, self.p) #Since the aeolian module is 
            # only run for winds above threshold, also run avalanching routine here
            self.s = aeolis.avalanching.avalanche(self.s, self.p)

        # grow vegetation
        if self.p['process_vegetation']:
            self.s = aeolis.vegetation.germinate(self.s, self.p)
            self.s = aeolis.vegetation.grow(self.s, self.p)

        # increment time
        self.t += self.dt * self.p['accfac']
        self._count('time')
        
        # Rotate gridparams back to original grid orientation
        self.s = self.grid_rotate(-self.p['alpha'])

        # Visualization of the model results after the first time step as a check for interpretation
        if self.c['time'] == 1 and self.p['visualization']:
            aeolis.inout.visualize_spatial(self.s, self.p)
            

    def finalize(self):
        '''Finalize model'''
        pass


    def get_current_time(self) -> float:
        '''
        Returns
        -------
        float
            Current simulation time
        '''

        return self.t

    def get_end_time(self) -> float:
        '''
        Returns
        -------
            Final simulation time
        '''

        return self.p['tstop']

    def get_start_time(self) -> float:
        '''
        Returns
        -------
            Initial simulation time
        '''

        return self.p['tstart']

    def get_var(self, var:str) -> Union[np.ndarray, int, float, str, list]:
        '''Returns spatial grid or model configuration parameter

        If the given variable name matches with a spatial grid, the
        spatial grid is returned. If not, the given variable name is
        matched with a model configuration parameter. If a match is
        found, the parameter value is returned. Otherwise, nothing is
        returned.

        Parameters
        ----------
        var:
            Name of spatial grid or model configuration parameter

        Returns
        -------
        np.ndarray or int, float, str or list
            Spatial grid or model configuration parameter

        Examples
        --------
        >>> # returns bathymetry grid
        ... model.get_var('zb')

        >>> # returns simulation duration
        ... model.get_var('tstop')

        See Also
        --------
        model.AeoLiS.set_var
        '''

        if var in self.s:
            if var in ['Ct', 'Cu']:
                return self.s[var] / self.p['accfac']
            else:
                return self.s[var]
        elif var in self.p:
            return self.p[var]
        else:
            return None


    def get_var_count(self) -> int:
        '''
        Returns
        -------
            Number of spatial grids
        '''

        return len(self.s)


    def get_var_name(self, i: int) -> Union[str, int]:
        '''Returns name of spatial grid by index (in alphabetical order)

        Parameters
        ----------
        i : 
            Index of spatial grid

        Returns
        -------
        str or -1
            Name of spatial grid or -1 in case index exceeds the number of grids
        '''

        if len(self.s) > i:
            return sorted(self.s.keys())[i]
        else:
            return -1

    def get_var_rank(self, var: str) -> int:
        '''Returns rank of spatial grid

        Parameters
        ----------
        var :
            Name of spatial grid

        Returns
        -------
            Rank of spatial grid or -1 if not found
        '''

        if var in self.s:
            return len(self.s[var].shape)
        else:
            return -1


    def get_var_shape(self, var:str) -> Union[Tuple, str]:
        '''Returns shape of spatial grid

        Parameters
        ----------
        var : 
            Name of spatial grid

        Returns
        -------
            Dimensions of spatial grid or -1 if not found
        '''

        if var in self.s:
            return self.s[var].shape
        else:
            return -1


    def get_var_type(self, var:str) -> Union[str, int]:
        '''Returns data type of variable in spatial grid

        Parameters
        ----------
        var :
            Name of spatial grid

        Returns
        -------
            Variable type of spatial grid or -1 if not found
        '''

        if var in self.s:
            return 'double'
        else:
            return -1

    def inq_compound(self) -> None:
        logger.log_and_raise('Method not yet implemented [inq_compound]', exc=NotImplementedError)


    def inq_compound_field(self) -> None:
        logger.log_and_raise('Method not yet implemented [inq_compound_field]', exc=NotImplementedError)


    def set_var(self, var:str, val:Union[np.ndarray, int, float, str, list]) -> None:
        '''Sets spatial grid or model configuration parameter

        If the given variable name matches with a spatial grid, the
        spatial grid is set. If not, the given variable name is
        matched with a model configuration parameter. If a match is
        found, the parameter value is set. Otherwise, nothing is set.

        Parameters
        ----------
        var : 
            Name of spatial grid or model configuration parameter
        val : 
            Spatial grid or model configuration parameter

        Examples
        --------
        >>> # set bathymetry grid
        ... model.set_var('zb', np.array([[0.,0., ... ,0.]]))

        >>> # set simulation duration
        ... model.set_var('tstop', 3600.)

        See Also
        --------
        model.AeoLiS.get_var
        '''

        if var in self.s:
            self.s[var] = val
        elif var in self.p:
            self.p[var] = val


    def set_var_index(self, i:int, val:np.ndarray) -> None:
        '''Set spatial grid by index (in alphabetical order)

        Parameters
        ----------
        i :
            Index of spatial grid
        val : 
            Spatial grid
        '''

        var = self.get_var_name(i)
        self.set_var(var, val)


    def set_var_slice(self) -> None:
        logger.log_and_raise('Method not yet implemented [set_var_slice]', exc=NotImplementedError)


    def set_timestep(self, dt:float=-1.0) -> bool:
        '''Determine optimal time step

        If no time step is given the optimal time step is
        determined. For explicit numerical schemes the time step is
        based in the Courant-Frierichs-Lewy (CFL) condition. For
        implicit numerical schemes the time step specified in the
        model configuration file is used. Alternatively, a preferred
        time step is given that is maximized by the CFL condition in
        case of an explicit numerical scheme.

        Returns True except when:

        1. No time step could be determined, for example when there is
        no wind and the numerical scheme is explicit. In this case the
        time step is set arbitrarily to one second.

        2. Or when the time step is smaller than -1. In this case the
        time is updated with the absolute value of the time step, but
        no model execution is performed. This funcionality can be used
        to skip fast-forward in time.

        Parameters
        ----------
        df :
            Preferred time step

        Returns
        -------
            False if determination of time step was unsuccessful, True otherwise

        '''

        if dt > 0.:
            self.dt = dt
        elif dt < -1:
            self.dt = dt
            self.t += np.abs(dt)
            return False
        else:
            self.dt = self.p['dt']

        if self.p['scheme'] == 'euler_forward':
            if self.p['CFL'] > 0.:
                dtref = np.max(np.abs(self.s['uws']) / self.s['ds']) + \
                        np.max(np.abs(self.s['uwn']) / self.s['dn'])
                if dtref > 0.:
                    self.dt = np.minimum(self.dt, self.p['CFL'] / dtref)
                # else:
                #     print(self.dt)
                    #self.t += self.dt
                    #self.dt = np.minimum(self.dt, 1.)
                    #return False
           
        if self.p['max_bedlevel_change'] != 999. and np.max(self.s['dzb']) != 0. and self.dt_prev != 0.:
            
            dt_zb = self.dt_prev * self.p['max_bedlevel_change'] / np.max(self.s['dzb'])
            self.dt = np.minimum(self.dt, dt_zb)
            
        self.p['dt_opt'] = self.dt

        return True

    def grid_rotate(self, angle:float) -> ModelState:
        
        s = self.s
        p = self.p
        
        s['x'], s['y'] = rotate(s['x'], s['y'], angle, origin=(np.mean(s['x']), np.mean(s['y'])))
        
        s['taus'], s['taun'] = rotate(s['taus'], s['taun'], angle, origin=(0, 0))
        s['taus0'], s['taun0'] = rotate(s['taus0'], s['taun0'], angle, origin=(0, 0))
        
        s['ustars'], s['ustarn'] = rotate(s['ustars'], s['ustarn'], angle, origin=(0, 0))
        s['ustars0'], s['ustarn0'] = rotate(s['ustars0'], s['ustarn0'], angle, origin=(0, 0))
        
        s['uws'], s['uwn'] = rotate(s['uws'], s['uwn'], angle, origin=(0, 0))
    
        for i in range(p['nfractions']):
            s['qs'][:,:,i], s['qn'][:,:,i] = rotate(s['qs'][:,:,i], s['qn'][:,:,i], angle, origin=(0, 0))
            s['us'][:,:,i], s['un'][:,:,i] = rotate(s['us'][:,:,i], s['un'][:,:,i], angle, origin=(0, 0))
        
        self.s['udir'] += angle
        
        return s
    
    def euler_forward(self) -> Any:
        '''Convenience function for explicit solver based on Euler forward scheme

        See Also
        --------
        model.AeoLiS.solve

        '''
        
        if self.p['solver'].lower() == 'trunk':
            solve = self.solve_EF(alpha=0., beta=1)
            #solve = self.solve(alpha=0., beta=1)
        elif self.p['solver'].lower() == 'pieter': 
            solve = self.solve_pieter(alpha=0., beta=1)
        elif self.p['solver'].lower() == 'steadystate':
            solve = self.solve_SS()
        elif self.p['solver'].lower() == 'steadystatepieter':
            solve = self.solve_steadystatepieter()

        return solve


    def euler_backward(self) -> Any:
        '''Convenience function for implicit solver based on Euler backward scheme

        See Also
        --------
        model.AeoLiS.solve

        '''
        
        if self.p['solver'].lower() == 'trunk':
            solve = self.solve(alpha=1., beta=1)
        elif self.p['solver'].lower() == 'pieter': 
            solve = self.solve_pieter(alpha=1., beta=1)
        elif self.p['solver'].lower() == 'steadystate':
            solve = self.solve_SS()
        elif self.p['solver'].lower() == 'steadystatepieter':
            solve = self.solve_steadystatepieter()
            
        return solve

    def crank_nicolson(self) -> Any:
        '''Convenience function for semi-implicit solver based on Crank-Nicolson scheme

        See Also
        --------
        model.AeoLiS.solve
        '''

        if self.p['solver'].lower() == 'trunk':
            solve = self.solve(alpha=.5, beta=1)
        elif self.p['solver'].lower() == 'pieter': 
            solve = self.solve_pieter(alpha=.5, beta=1)
        elif self.p['solver'].lower() == 'steadystate':
            solve = self.solve_steadystate()
        elif self.p['solver'].lower() == 'steadystatepieter':
            solve = self.solve_steadystatepieter()

        return solve


    def solve_steadystate(self) -> dict:
        '''Implements the steady state solution
        '''
        # upwind scheme:
        beta = 1. 
        
        l = self.l
        s = self.s
        p = self.p

        Ct = s['Ct'].copy()
        pickup = s['pickup'].copy()

        # compute transport weights for all sediment fractions
        w_init, w_air, w_bed = aeolis.transport.compute_weights(s, p)

        if self.t == 0.:
            # use initial guess for first time step
            if p['grain_dist'] != None:
                w = p['grain_dist'].reshape((1,1,-1))
                w = w.repeat(p['ny']+1, axis=0)
                w = w.repeat(p['nx']+1, axis=1)
            else:
                w = w_init.copy()
        else:
            w = w_init.copy()

        # set model state properties that are added to warnings and errors
        logprops = dict(minwind=s['uw'].min(),
                        maxdrop=(l['uw']-s['uw']).max(),
                        time=self.t,
                        dt=self.dt)
            
        nf = p['nfractions']     
                
        us = np.zeros((p['ny']+1,p['nx']+1))
        un = np.zeros((p['ny']+1,p['nx']+1))
        
        us_plus = np.zeros((p['ny']+1,p['nx']+1))
        un_plus = np.zeros((p['ny']+1,p['nx']+1))
        
        us_min = np.zeros((p['ny']+1,p['nx']+1))
        un_min = np.zeros((p['ny']+1,p['nx']+1))

        Cs = np.zeros(us.shape)
        Cn = np.zeros(un.shape)
        
        Cs_plus = np.zeros(us.shape)
        Cn_plus = np.zeros(un.shape)
        
        Cs_min = np.zeros(us.shape)
        Cn_min = np.zeros(un.shape)
        
        for i in range(nf):
            us[:,:] = s['us'][:,:,i] 
            un[:,:] = s['un'][:,:,i] 
            
            us_plus[:,1:] = s['us'][:,:-1,i] 
            un_plus[1:,:] = s['un'][:-1,:,i] 
            
            us_min[:,:-1] = s['us'][:,1:,i]
            un_min[:-1,:] = s['un'][1:,:,i]
        
            #boundary values
            us[:,0]  = s['us'][:,0,i]
            un[0,:]  = s['un'][0,:,i]
            
            us_plus[:,0]  = s['us'][:,0,i]
            un_plus[0,:]  = s['un'][0,:,i]
            
            us_min[:,-1]  = s['us'][:,-1,i]
            un_min[-1,:]  = s['un'][-1,:,i]
            
            
            # define matrix coefficients to solve linear system of equations        
            Cs = s['dn'] * s['dsdni'] * us[:,:]  
            Cn = s['ds'] * s['dsdni'] * un[:,:] 

            Cs_plus = s['dn'] * s['dsdni'] * us_plus[:,:]  
            Cn_plus = s['ds'] * s['dsdni'] * un_plus[:,:]
            
            Cs_min = s['dn'] * s['dsdni'] * us_min[:,:]  
            Cn_min = s['ds'] * s['dsdni'] * un_min[:,:]
            
            
            Ti = 1 / p['T']
            
            beta = abs(beta)
            if beta >= 1.:
                # define upwind direction
                ixs = np.asarray(us[:,:] >= 0., dtype=float)
                ixn = np.asarray(un[:,:] >= 0., dtype=float)
                sgs = 2. * ixs - 1.
                sgn = 2. * ixn - 1.
            
            else:
                # or centralizing weights
                ixs = beta + np.zeros(us)
                ixn = beta + np.zeros(un)
                sgs = np.zeros(us)
                sgn = np.zeros(un)

            # initialize matrix diagonals
            A0 = np.zeros(s['zb'].shape)
            Apx = np.zeros(s['zb'].shape)
            Ap1 = np.zeros(s['zb'].shape)
            Ap2 = np.zeros(s['zb'].shape)
            Amx = np.zeros(s['zb'].shape)
            Am1 = np.zeros(s['zb'].shape)
            Am2 = np.zeros(s['zb'].shape)

            # populate matrix diagonals
            A0  = sgs * Cs + sgn * Cn + Ti
            Apx = Cn_min * (1. - ixn)
            Ap1 = Cs_min * (1. - ixs)
            Amx = -Cn_plus * ixn
            Am1 = -Cs_plus * ixs    

            # add boundaries
            A0[:,0] = 1.
            Apx[:,0] = 0.
            Amx[:,0] = 0.
            Am2[:,0] = 0.
            Am1[:,0] = 0.

            A0[:,-1] = 1.
            Apx[:,-1] = 0.
            Ap1[:,-1] = 0.
            Ap2[:,-1] = 0.
            Amx[:,-1] = 0.

            if p['boundary_offshore'] == 'flux':
                Ap2[:,0] = 0.
                Ap1[:,0] = 0.
            elif p['boundary_offshore'] == 'constant':
                Ap2[:,0] = 0.
                Ap1[:,0] = 0.
            elif p['boundary_offshore'] == 'uniform':
                Ap2[:,0] = 0.
                Ap1[:,0] = -1.
            elif p['boundary_offshore'] == 'gradient':
                Ap2[:,0] = s['ds'][:,1] / s['ds'][:,2]
                Ap1[:,0] = -1. - s['ds'][:,1] / s['ds'][:,2]
            elif p['boundary_offshore'] == 'circular':
                logger.log_and_raise('Cross-shore cricular boundary condition not yet implemented', exc=NotImplementedError)
            else:
                logger.log_and_raise('Unknown offshore boundary condition [%s]' % self.p['boundary_offshore'], exc=ValueError)

            if p['boundary_onshore'] == 'flux':                              
                Am2[:,-1] = 0.
                Am1[:,-1] = 0.            
            elif p['boundary_onshore'] == 'constant':                              
                Am2[:,-1] = 0.
                Am1[:,-1] = 0.
            elif p['boundary_onshore'] == 'uniform':
                Am2[:,-1] = 0.
                Am1[:,-1] = -1.
            elif p['boundary_onshore'] == 'gradient':
                Am2[:,-1] = s['ds'][:,-2] / s['ds'][:,-3]
                Am1[:,-1] = -1. - s['ds'][:,-2] / s['ds'][:,-3]
            elif p['boundary_offshore'] == 'circular':
                logger.log_and_raise('Cross-shore cricular boundary condition not yet implemented', exc=NotImplementedError)
            else:
                logger.log_and_raise('Unknown onshore boundary condition [%s]' % self.p['boundary_onshore'], exc=ValueError)

            if p['boundary_lateral'] == 'constant':
                A0[0,:] = 1.
                Apx[0,:] = 0.
                Ap1[0,:] = 0.
                Amx[0,:] = 0.
                Am1[0,:] = 0.
                
                A0[-1,:] = 1.
                Apx[-1,:] = 0.
                Ap1[-1,:] = 0.
                Amx[-1,:] = 0.
                Am1[-1,:] = 0.
            
                #logger.log_and_raise('Lateral constant boundary condition not yet implemented', exc=NotImplementedError)
            elif p['boundary_lateral'] == 'uniform':
                logger.log_and_raise('Lateral uniform boundary condition not yet implemented', exc=NotImplementedError)
            elif p['boundary_lateral'] == 'gradient':
                logger.log_and_raise('Lateral gradient boundary condition not yet implemented', exc=NotImplementedError)
            elif p['boundary_lateral'] == 'circular':
                pass
            else:
                logger.log_and_raise('Unknown lateral boundary condition [%s]' % self.p['boundary_lateral'], exc=ValueError)

            # construct sparse matrix
            if p['ny'] > 0:
                j = p['nx']+1
                A = scipy.sparse.diags((Apx.flatten()[:j],
                                        Amx.flatten()[j:],
                                        Am2.flatten()[2:],
                                        Am1.flatten()[1:],
                                        A0.flatten(),
                                        Ap1.flatten()[:-1],
                                        Ap2.flatten()[:-2],
                                        Apx.flatten()[j:],
                                        Amx.flatten()[:j]),
                                       (-j*p['ny'],-j,-2,-1,0,1,2,j,j*p['ny']), format='csr')
            else:
                A = scipy.sparse.diags((Am2.flatten()[2:],
                                        Am1.flatten()[1:],
                                        A0.flatten(),
                                        Ap1.flatten()[:-1],
                                        Ap2.flatten()[:-2]),
                                       (-2,-1,0,1,2), format='csr')

            # solve transport for each fraction separately using latest
            # available weights

            # renormalize weights for all fractions equal or larger
            # than the current one such that the sum of all weights is
            # unity
            w = aeolis.transport.renormalize_weights(w, i)

            # iteratively find a solution of the linear system that
            # does not violate the availability of sediment in the bed
            for n in range(p['max_iter']):
                self._count('matrixsolve')

                # compute saturation levels
                ix = s['Cu'] > 0.
                S_i = np.zeros(s['Cu'].shape)
                S_i[ix] = s['Ct'][ix] / s['Cu'][ix]
                s['S'] = S_i.sum(axis=-1)

                # create the right hand side of the linear system
                y_i = np.zeros(s['zb'].shape)
                
                y_i[:,1:-1] = (
                    (w[:,1:-1,i] * s['Cuf'][:,1:-1,i] * Ti) * (1. - s['S'][:,1:-1]) +
                    (w[:,1:-1,i] * s['Cu'][:,1:-1,i] * Ti) * s['S'][:,1:-1]
                    )

                # add boundaries
                if p['boundary_offshore'] == 'flux':
                    y_i[:,0] = p['offshore_flux'] * s['Cu0'][:,0,i] 
                if p['boundary_onshore'] == 'flux':
                    y_i[:,-1] = p['onshore_flux'] * s['Cu0'][:,-1,i] 
                    
                if p['boundary_offshore'] == 'constant':
                    y_i[:,0] = p['constant_offshore_flux'] / s['u'][:,0,i] 
                if p['boundary_onshore'] == 'constant':
                    y_i[:,-1] = p['constant_onshore_flux'] / s['u'][:,-1,i]

                # solve system with current weights
                Ct_i = scipy.sparse.linalg.spsolve(A, y_i.flatten())
                Ct_i = prevent_tiny_negatives(Ct_i, p['max_error'])
                   
                # check for negative values
                if Ct_i.min() < 0.:
                    ix = Ct_i < 0.

                    logger.warning(format_log('Removing negative concentrations',
                                              nrcells=np.sum(ix),
                                              fraction=i,
                                              iteration=n,
                                              minvalue=Ct_i.min(),
                                              coords=np.argwhere(ix.reshape(y_i.shape)),
                                              **logprops))

                    Ct_i[~ix] *= 1. + Ct_i[ix].sum() / Ct_i[~ix].sum()
                    Ct_i[ix] = 0.

                # determine pickup and deficit for current fraction
                Cu_i = s['Cu'][:,:,i].flatten()
                mass_i = s['mass'][:,:,0,i].flatten()
                w_i = w[:,:,i].flatten()
                pickup_i = (w_i * Cu_i - Ct_i) / p['T'] * self.dt
                deficit_i = pickup_i - mass_i
                ix = (deficit_i > p['max_error']) \
                     & (w_i * Cu_i > 0.)

                # quit the iteration if there is no deficit, otherwise
                # back-compute the maximum weight allowed to get zero
                # deficit for the current fraction and progress to
                # the next iteration step
                if not np.any(ix):
                    logger.debug(format_log('Iteration converged',
                                            steps=n,
                                            fraction=i,
                                            **logprops))
                    pickup_i = np.minimum(pickup_i, mass_i)
                    break
                else:
                    w_i[ix] = (mass_i[ix] * p['T'] / self.dt \
                               + Ct_i[ix]) / Cu_i[ix]
                    w[:,:,i] = w_i.reshape(y_i.shape)

            # throw warning if the maximum number of iterations was reached
            if np.any(ix):
                logger.warning(format_log('Iteration not converged',
                                          nrcells=np.sum(ix),
                                          fraction=i,
                                          **logprops))

            # check for unexpected negative values
            if Ct_i.min() < 0:
                logger.warning(format_log('Negative concentrations',
                                          nrcells=np.sum(Ct_i<0.),
                                          fraction=i,
                                          minvalue=Ct_i.min(),
                                          **logprops))
            if w_i.min() < 0:
                logger.warning(format_log('Negative weights',
                                          nrcells=np.sum(w_i<0),
                                          fraction=i,
                                          minvalue=w_i.min(),
                                          **logprops))

            Ct[:,:,i] = Ct_i.reshape(y_i.shape)
            pickup[:,:,i] = pickup_i.reshape(y_i.shape)

        # check if there are any cells where the sum of all weights is
        # smaller than unity. these cells are supply-limited for all
        # fractions. Log these events.
        ix = 1. - np.sum(w, axis=2) > p['max_error']
        if np.any(ix):
            self._count('supplylim')
            logger.warning(format_log('Ran out of sediment',
                                      nrcells=np.sum(ix),
                                      minweight=np.sum(w, axis=-1).min(),
                                      **logprops))
           
        
        qs = Ct * s['us'] 
        qn = Ct * s['un'] 
        q = np.hypot(qs, qn)


        return dict(Ct=Ct,
                    qs=qs,
                    qn=qn,
                    pickup=pickup,
                    w=w,
                    w_init=w_init,
                    w_air=w_air,
                    w_bed=w_bed,
                    q=q)
        
        
    def solve(self, alpha:float=.5, beta:float=1.) -> dict:
        '''Implements the explicit Euler forward, implicit Euler backward and semi-implicit Crank-Nicolson numerical schemes

        Determines weights of sediment fractions, sediment pickup and
        instantaneous sediment concentration. Returns a partial
        spatial grid dictionary that can be used to update the global
        spatial grid dictionary.

        Parameters
        ----------
        alpha :
            Implicitness coefficient (0.0 for Euler forward, 1.0 for Euler backward or 0.5 for Crank-Nicolson, default=0.5)
        beta : 
            Centralization coefficient (1.0 for upwind or 0.5 for centralized, default=1.0)

        Returns
        -------
            Partial spatial grid dictionary

        Examples
        --------
        >>> model.s.update(model.solve(alpha=1., beta=1.) # euler backward

        >>> model.s.update(model.solve(alpha=.5, beta=1.) # crank-nicolson

        See Also
        --------
        model.AeoLiS.euler_forward
        model.AeoLiS.euler_backward
        model.AeoLiS.crank_nicolson
        transport.compute_weights
        transport.renormalize_weights

        '''

        l = self.l
        s = self.s
        p = self.p

        Ct = s['Ct'].copy()
        pickup = s['pickup'].copy()

        # compute transport weights for all sediment fractions
        w_init, w_air, w_bed = aeolis.transport.compute_weights(s, p)

        if self.t == 0.:
            if type(p['bedcomp_file']) == np.ndarray:
                w = w_init.copy()
            else:
                # use initial guess for first time step
                w = p['grain_dist'].reshape((1,1,-1))
                w = w.repeat(p['ny']+1, axis=0)
                w = w.repeat(p['nx']+1, axis=1)
        else:
            w = w_init.copy()

        # set model state properties that are added to warnings and errors
        logprops = dict(minwind=s['uw'].min(),
                        maxdrop=(l['uw']-s['uw']).max(),
                        time=self.t,
                        dt=self.dt)
            
        nf = p['nfractions']
        
        us = np.zeros((p['ny']+1,p['nx']+1))
        un = np.zeros((p['ny']+1,p['nx']+1))
        
        us_plus = np.zeros((p['ny']+1,p['nx']+1))
        un_plus = np.zeros((p['ny']+1,p['nx']+1))
        
        us_min = np.zeros((p['ny']+1,p['nx']+1))
        un_min = np.zeros((p['ny']+1,p['nx']+1))

        Cs = np.zeros(us.shape)
        Cn = np.zeros(un.shape)
        
        Cs_plus = np.zeros(us.shape)
        Cn_plus = np.zeros(un.shape)
        
        Cs_min = np.zeros(us.shape)
        Cn_min = np.zeros(un.shape)
        
        
        for i in range(nf):
            
            us[:,:] = s['us'][:,:,i] 
            un[:,:] = s['un'][:,:,i] 
            
            us_plus[:,1:] = s['us'][:,:-1,i] 
            un_plus[1:,:] = s['un'][:-1,:,i] 
            
            us_min[:,:-1] = s['us'][:,1:,i]
            un_min[:-1,:] = s['un'][1:,:,i]
        
            #boundary values            
            us_plus[:,0]  = s['us'][:,0,i]
            un_plus[0,:]  = s['un'][0,:,i]
            
            us_min[:,-1]  = s['us'][:,-1,i]
            un_min[-1,:]  = s['un'][-1,:,i]
            
            
            # define matrix coefficients to solve linear system of equations        
            Cs = self.dt * s['dn'] * s['dsdni'] * us[:,:]  
            Cn = self.dt * s['ds'] * s['dsdni'] * un[:,:] 

            Cs_plus = self.dt * s['dn'] * s['dsdni'] * us_plus[:,:]  
            Cn_plus = self.dt * s['ds'] * s['dsdni'] * un_plus[:,:]
            
            Cs_min = self.dt * s['dn'] * s['dsdni'] * us_min[:,:]  
            Cn_min = self.dt * s['ds'] * s['dsdni'] * un_min[:,:]
            
            Ti = self.dt / p['T']          

            
            beta = abs(beta)
            if beta >= 1.:
                # define upwind direction
                ixs = np.asarray(s['us'][:,:,i] >= 0., dtype=float)
                ixn = np.asarray(s['un'][:,:,i] >= 0., dtype=float)
                sgs = 2. * ixs - 1.
                sgn = 2. * ixn - 1.
            
            else:
                # or centralizing weights
                ixs = beta + np.zeros(Cs.shape)
                ixn = beta + np.zeros(Cn.shape)
                sgs = np.zeros(Cs.shape)
                sgn = np.zeros(Cn.shape)
                
            # initialize matrix diagonals
            A0 = np.zeros(s['zb'].shape)
            Apx = np.zeros(s['zb'].shape)
            Ap1 = np.zeros(s['zb'].shape)
            Ap2 = np.zeros(s['zb'].shape)
            Amx = np.zeros(s['zb'].shape)
            Am1 = np.zeros(s['zb'].shape)
            Am2 = np.zeros(s['zb'].shape)

            # populate matrix diagonals
            A0  = 1. + (sgs * Cs + sgn * Cn + Ti) * alpha
            Apx = Cn_min * alpha * (1. - ixn)
            Ap1 = Cs_min * alpha * (1. - ixs)
            Amx = -Cn_plus * alpha * ixn
            Am1 = -Cs_plus * alpha * ixs    

            # add boundaries
            A0[:,0] = 1.
            Apx[:,0] = 0.
            Amx[:,0] = 0.
            Am2[:,0] = 0.
            Am1[:,0] = 0.

            A0[:,-1] = 1.
            Apx[:,-1] = 0.
            Ap1[:,-1] = 0.
            Ap2[:,-1] = 0.
            Amx[:,-1] = 0.

            if (p['boundary_offshore'] == 'flux') | (p['boundary_offshore'] == 'noflux'):
                Ap2[:,0] = 0.
                Ap1[:,0] = 0.
            elif p['boundary_offshore'] == 'constant':
                Ap2[:,0] = 0.
                Ap1[:,0] = 0.
            elif p['boundary_offshore'] == 'uniform':
                Ap2[:,0] = 0.
                Ap1[:,0] = -1.
            elif p['boundary_offshore'] == 'gradient':
                Ap2[:,0] = s['ds'][:,1] / s['ds'][:,2]
                Ap1[:,0] = -1. - s['ds'][:,1] / s['ds'][:,2]
            elif p['boundary_offshore'] == 'circular':
                logger.log_and_raise('Cross-shore cricular boundary condition not yet implemented', exc=NotImplementedError)
            else:
                logger.log_and_raise('Unknown offshore boundary condition [%s]' % self.p['boundary_offshore'], exc=ValueError)

            if (p['boundary_onshore'] == 'flux') | (p['boundary_offshore'] == 'noflux'):                              
                Am2[:,-1] = 0.
                Am1[:,-1] = 0.            
            elif p['boundary_onshore'] == 'constant':                              
                Am2[:,-1] = 0.
                Am1[:,-1] = 0.
            elif p['boundary_onshore'] == 'uniform':
                Am2[:,-1] = 0.
                Am1[:,-1] = -1.
            elif p['boundary_onshore'] == 'gradient':
                Am2[:,-1] = s['ds'][:,-2] / s['ds'][:,-3]
                Am1[:,-1] = -1. - s['ds'][:,-2] / s['ds'][:,-3]
            elif p['boundary_offshore'] == 'circular':
                logger.log_and_raise('Cross-shore cricular boundary condition not yet implemented', exc=NotImplementedError)
            else:
                logger.log_and_raise('Unknown onshore boundary condition [%s]' % self.p['boundary_onshore'], exc=ValueError)

            if p['boundary_lateral'] == 'constant':
                A0[0,:] = 1.
                Apx[0,:] = 0.
                Ap1[0,:] = 0.
                Amx[0,:] = 0.
                Am1[0,:] = 0.
                
                A0[-1,:] = 1.
                Apx[-1,:] = 0.
                Ap1[-1,:] = 0.
                Amx[-1,:] = 0.
                Am1[-1,:] = 0.
            
                #logger.log_and_raise('Lateral constant boundary condition not yet implemented', exc=NotImplementedError)
            elif p['boundary_lateral'] == 'uniform':
                logger.log_and_raise('Lateral uniform boundary condition not yet implemented', exc=NotImplementedError)
            elif p['boundary_lateral'] == 'gradient':
                logger.log_and_raise('Lateral gradient boundary condition not yet implemented', exc=NotImplementedError)
            elif p['boundary_lateral'] == 'circular':
                pass
            else:
                logger.log_and_raise('Unknown lateral boundary condition [%s]' % self.p['boundary_lateral'], exc=ValueError)

            # construct sparse matrix
            if p['ny'] > 0:
                j = p['nx']+1
                A = scipy.sparse.diags((Apx.flatten()[:j],
                                        Amx.flatten()[j:],
                                        Am2.flatten()[2:],
                                        Am1.flatten()[1:],
                                        A0.flatten(),
                                        Ap1.flatten()[:-1],
                                        Ap2.flatten()[:-2],
                                        Apx.flatten()[j:],
                                        Amx.flatten()[:j]),
                                       (-j*p['ny'],-j,-2,-1,0,1,2,j,j*p['ny']), format='csr')
            else:
                A = scipy.sparse.diags((Am2.flatten()[2:],
                                        Am1.flatten()[1:],
                                        A0.flatten(),
                                        Ap1.flatten()[:-1],
                                        Ap2.flatten()[:-2]),
                                       (-2,-1,0,1,2), format='csr')

            # solve transport for each fraction separately using latest
            # available weights

            # renormalize weights for all fractions equal or larger
            # than the current one such that the sum of all weights is
            # unity
            # Christa: seems to have no significant effect on weights, 
            # numerical check to prevent any deviation from unity
            w = aeolis.transport.renormalize_weights(w, i)            

            # iteratively find a solution of the linear system that
            # does not violate the availability of sediment in the bed
            for n in range(p['max_iter']):
                self._count('matrixsolve')

                # compute saturation levels
                ix = s['Cu'] > 0.
                S_i = np.zeros(s['Cu'].shape)
                S_i[ix] = s['Ct'][ix] / s['Cu'][ix]
                s['S'] = S_i.sum(axis=-1)

                # create the right hand side of the linear system
                y_i = np.zeros(s['zb'].shape)
                y_im = np.zeros(s['zb'].shape)  # implicit terms
                y_ex = np.zeros(s['zb'].shape)  # explicit terms
                
                y_im[:,1:-1] = (
                    (w[:,1:-1,i] * s['Cuf'][:,1:-1,i] * Ti) * (1. - s['S'][:,1:-1]) +
                    (w[:,1:-1,i] * s['Cu'][:,1:-1,i] * Ti) * s['S'][:,1:-1]
                    )
                
                y_ex[:,1:-1] = (
                    (l['w'][:,1:-1,i] * l['Cuf'][:,1:-1,i] * Ti) * (1. - s['S'][:,1:-1]) \
                    + (l['w'][:,1:-1,i] * l['Cu'][:,1:-1,i] * Ti) * s['S'][:,1:-1] \
                    - (
                        sgs[:,1:-1] * Cs[:,1:-1] +\
                        sgn[:,1:-1] * Cn[:,1:-1] + Ti
                    ) * l['Ct'][:,1:-1,i] \
                    + ixs[:,1:-1] * Cs_plus[:,1:-1] * l['Ct'][:,:-2,i] \
                    - (1. - ixs[:,1:-1]) * Cs_min[:,1:-1] * l['Ct'][:,2:,i] \
                    + ixn[:,1:-1] * Cn_plus[:,1:-1] * np.roll(l['Ct'][:,1:-1,i], 1, axis=0) \
                    - (1. - ixn[:,1:-1]) * Cn_min[:,1:-1] * np.roll(l['Ct'][:,1:-1,i], -1, axis=0) \
                    )
                
                y_i[:,1:-1] = l['Ct'][:,1:-1,i] + alpha * y_im[:,1:-1] + (1. - alpha) * y_ex[:,1:-1]

                # add boundaries
                if p['boundary_offshore'] == 'flux':
                    y_i[:,0] = p['offshore_flux'] * s['Cu0'][:,0,i] 
                if p['boundary_onshore'] == 'flux':
                    y_i[:,-1] = p['onshore_flux'] * s['Cu0'][:,-1,i] 
                    
                if p['boundary_offshore'] == 'constant':
                    y_i[:,0] = p['constant_offshore_flux'] / s['u'][:,0,i] 
                if p['boundary_onshore'] == 'constant':
                    y_i[:,-1] = p['constant_onshore_flux'] / s['u'][:,-1,i]

                # solve system with current weights
                Ct_i = scipy.sparse.linalg.spsolve(A, y_i.flatten())
                Ct_i = prevent_tiny_negatives(Ct_i, p['max_error'])
                
                # check for negative values
                if Ct_i.min() < 0.:
                    ix = Ct_i < 0.

                    logger.warning(format_log('Removing negative concentrations',
                                              nrcells=np.sum(ix),
                                              fraction=i,
                                              iteration=n,
                                              minvalue=Ct_i.min(),
                                              coords=np.argwhere(ix.reshape(y_i.shape)),
                                              **logprops))

                    if Ct_i[~ix].sum() != 0:
                        Ct_i[~ix] *= 1. + Ct_i[ix].sum() / Ct_i[~ix].sum()
                    else:
                        Ct_i[~ix] = 0

                    #Ct_i[~ix] *= 1. + Ct_i[ix].sum() / Ct_i[~ix].sum()
                    Ct_i[ix] = 0.

                # determine pickup and deficit for current fraction
                Cu_i = s['Cu'][:,:,i].flatten()
                mass_i = s['mass'][:,:,0,i].flatten()
                w_i = w[:,:,i].flatten()
                pickup_i = (w_i * Cu_i - Ct_i) / p['T'] * self.dt
                deficit_i = pickup_i - mass_i
                ix = (deficit_i > p['max_error']) \
                     & (w_i * Cu_i > 0.)

                # quit the iteration if there is no deficit, otherwise
                # back-compute the maximum weight allowed to get zero
                # deficit for the current fraction and progress to
                # the next iteration step
                if not np.any(ix):
                    logger.debug(format_log('Iteration converged',
                                            steps=n,
                                            fraction=i,
                                            **logprops))
                    pickup_i = np.minimum(pickup_i, mass_i)
                    break
                else:
                    w_i[ix] = (mass_i[ix] * p['T'] / self.dt \
                               + Ct_i[ix]) / Cu_i[ix]
                    w[:,:,i] = w_i.reshape(y_i.shape)

            # throw warning if the maximum number of iterations was reached
            if np.any(ix):
                logger.warning(format_log('Iteration not converged',
                                          nrcells=np.sum(ix),
                                          fraction=i,
                                          **logprops))

            # check for unexpected negative values
            if Ct_i.min() < 0:
                logger.warning(format_log('Negative concentrations',
                                          nrcells=np.sum(Ct_i<0.),
                                          fraction=i,
                                          minvalue=Ct_i.min(),
                                          **logprops))
            if w_i.min() < 0:
                logger.warning(format_log('Negative weights',
                                          nrcells=np.sum(w_i<0),
                                          fraction=i,
                                          minvalue=w_i.min(),
                                          **logprops))

            Ct[:,:,i] = Ct_i.reshape(y_i.shape)
            pickup[:,:,i] = pickup_i.reshape(y_i.shape)

        # check if there are any cells where the sum of all weights is
        # smaller than unity. these cells are supply-limited for all
        # fractions. Log these events.
        ix = 1. - np.sum(w, axis=2) > p['max_error']
        if np.any(ix):
            self._count('supplylim')
            # logger.warning(format_log('Ran out of sediment',
            #                           nrcells=np.sum(ix),
            #                           minweight=np.sum(w, axis=-1).min(),
            #                           **logprops))
           
        qs = Ct * s['us'] 
        qn = Ct * s['un'] 

        return dict(Ct=Ct,
                    qs=qs,
                    qn=qn,
                    pickup=pickup,
                    w=w,
                    w_init=w_init,
                    w_air=w_air,
                    w_bed=w_bed)
        
    #@njit
    def solve_EF(self, alpha:float=0., beta:float=1.) -> dict:
        '''Implements the explicit Euler forward, implicit Euler backward and semi-implicit Crank-Nicolson numerical schemes

        Determines weights of sediment fractions, sediment pickup and
        instantaneous sediment concentration. Returns a partial
        spatial grid dictionary that can be used to update the global
        spatial grid dictionary.

        Parameters
        ----------
        alpha :
            Implicitness coefficient (0.0 for Euler forward, 1.0 for Euler backward or 0.5 for Crank-Nicolson, default=0.5)
        beta : 
            Centralization coefficient (1.0 for upwind or 0.5 for centralized, default=1.0)

        Returns
        -------
            Partial spatial grid dictionary

        Examples
        --------
        >>> model.s.update(model.solve(alpha=1., beta=1.) # euler backward

        >>> model.s.update(model.solve(alpha=.5, beta=1.) # crank-nicolson

        See Also
        --------
        model.AeoLiS.euler_forward
        model.AeoLiS.euler_backward
        model.AeoLiS.crank_nicolson
        transport.compute_weights
        transport.renormalize_weights

        '''

        l = self.l
        s = self.s
        p = self.p

        Ct = s['Ct'].copy()
        pickup = s['pickup'].copy()
        Ts = p['T']

        # compute transport weights for all sediment fractions
        w_init, w_air, w_bed = aeolis.transport.compute_weights(s, p)

        if self.t == 0.:
            if type(p['bedcomp_file']) == np.ndarray:
                w = w_init.copy()
            else:
                # use initial guess for first time step
                w = p['grain_dist'].reshape((1,1,-1))
                w = w.repeat(p['ny']+1, axis=0)
                w = w.repeat(p['nx']+1, axis=1)
        else:
            w = w_init.copy()

        # set model state properties that are added to warnings and errors
        logprops = dict(minwind=s['uw'].min(),
                        maxdrop=(l['uw']-s['uw']).max(),
                        time=self.t,
                        dt=self.dt)
            
        nf = p['nfractions']
            
        
        for i in range(nf):

            if 1:
                #define 4 quadrants based on wind directions
                ix1 = ((s['us'][:,:,0]>=0) & (s['un'][:,:,0]>=0))
                ix2 = ((s['us'][:,:,0]<0) & (s['un'][:,:,0]>=0))
                ix3 = ((s['us'][:,:,0]<0) & (s['un'][:,:,0]<0))
                ix4 = ((s['us'][:,:,0]>0) & (s['un'][:,:,0]<0))
                
                # initiate solution matrix including ghost cells to accomodate boundaries
                Ct_s = np.zeros((Ct.shape[0]+2,Ct.shape[1]+2))
                # populate solution matrix with previous concentration results
                Ct_s[1:-1,1:-1] = Ct[:,:,i]
                
                #set upwind boundary condition
                Ct_s[:,0:2]=0
                #circular boundary condition in lateral directions
                Ct_s[0,:]=Ct_s[-2,:]
                Ct_s[-1,:]=Ct_s[1,:]
                # using the Euler forward scheme we can calculate pickup first based on the previous timestep
                # there is no need for iteration
                pickup[:,:,i] = self.dt*(np.minimum(s['Cu'][:,:,i],s['mass'][:,:,0,i]+Ct[:,:,i])-Ct[:,:,i])/Ts
                
                #solve for all 4 quadrants in one step using logical indexing
                Ct_s[1:-1,1:-1] = Ct_s[1:-1,1:-1] + \
                    ix1*(-self.dt*s['us'][:,:,i]*(Ct_s[1:-1,1:-1]-Ct_s[1:-1,:-2])/s['ds'] \
                         -self.dt*s['un'][:,:,i]*(Ct_s[1:-1,1:-1]-Ct_s[:-2,1:-1])/s['dn']) +\
                    ix2*(+self.dt*s['us'][:,:,i]*(Ct_s[1:-1,1:-1]-Ct_s[1:-1,2:])/s['ds'] \
                         -self.dt*s['un'][:,:,i]*(Ct_s[1:-1,1:-1]-Ct_s[:-2,1:-1])/s['dn']) +\
                    ix3*(+self.dt*s['us'][:,:,i]*(Ct_s[1:-1,1:-1]-Ct_s[1:-1,2:])/s['ds'] \
                         +self.dt*s['un'][:,:,i]*(Ct_s[1:-1,1:-1]-Ct_s[2:,1:-1])/s['dn']) +\
                    ix4*(-self.dt*s['us'][:,:,i]*(Ct_s[1:-1,1:-1]-Ct_s[1:-1,:-2])/s['ds'] \
                         +self.dt*s['un'][:,:,i]*(Ct_s[1:-1,1:-1]-Ct_s[2:,1:-1])/s['dn']) \
                    + pickup[:,:,i]
                
                # define Ct as a subset of Ct_s (eliminating the boundaries)
                Ct[:,:,i] = Ct_s[1:-1,1:-1] 
         
        qs = Ct * s['us'] 
        qn = Ct * s['un'] 

        return dict(Ct=Ct,
                    qs=qs,
                    qn=qn,
                    pickup=pickup,
                    w=w,
                    w_init=w_init,
                    w_air=w_air,
                    w_bed=w_bed)
        
    #@njit
    def solve_SS(self, alpha:float=0., beta:float=1.) -> dict:
        '''Implements the explicit Euler forward, implicit Euler backward and semi-implicit Crank-Nicolson numerical schemes

        Determines weights of sediment fractions, sediment pickup and
        instantaneous sediment concentration. Returns a partial
        spatial grid dictionary that can be used to update the global
        spatial grid dictionary.

        Parameters
        ----------
        alpha :
            Implicitness coefficient (0.0 for Euler forward, 1.0 for Euler backward or 0.5 for Crank-Nicolson, default=0.5)
        beta : 
            Centralization coefficient (1.0 for upwind or 0.5 for centralized, default=1.0)

        Returns
        -------
            Partial spatial grid dictionary

        Examples
        --------
        >>> model.s.update(model.solve(alpha=1., beta=1.) # euler backward

        >>> model.s.update(model.solve(alpha=.5, beta=1.) # crank-nicolson

        See Also
        --------
        model.AeoLiS.euler_forward
        model.AeoLiS.euler_backward
        model.AeoLiS.crank_nicolson
        transport.compute_weights
        transport.renormalize_weights

        '''

        l = self.l
        s = self.s
        p = self.p

        Ct = s['Ct'].copy()
        pickup = s['pickup'].copy()
        Ts = p['T']

        # compute transport weights for all sediment fractions
        w_init, w_air, w_bed = aeolis.transport.compute_weights(s, p)

        if self.t == 0.:
            if type(p['bedcomp_file']) == np.ndarray:
                w = w_init.copy()
            else:
                # use initial guess for first time step
                # when p['grain_dist'] has 2 dimensions take the first row otherwise take the only row
                if len(p['grain_dist'].shape) == 2:
                    w = p['grain_dist'][0,:].reshape((1,1,-1))
                else:
                    w = p['grain_dist'].reshape((1,1,-1))
                    
                w = w.repeat(p['ny']+1, axis=0)
                w = w.repeat(p['nx']+1, axis=1)
        else:
            w = w_init.copy()

        # set model state properties that are added to warnings and errors
        logprops = dict(minwind=s['uw'].min(),
                        maxdrop=(l['uw']-s['uw']).max(),
                        time=self.t,
                        dt=self.dt)
            
        nf = p['nfractions']
            
        
        for i in range(nf):
            

            if 1:
                #print('sweep')

                # initiate emmpty solution matrix, this will effectively kill time dependence and create steady state.
                Ct = np.zeros(Ct.shape)
                
                if p['boundary_offshore'] == 'flux':
                    Ct[:,0,0] =  s['Cu0'][:,0,0] 

                if p['boundary_onshore'] == 'flux':
                    Ct[:,-1,0] =  s['Cu0'][:,-1,0] 

                if p['boundary_offshore'] == 'circular':
                    Ct[:,0,0] =  -1                
                    Ct[:,-1,0] =  -1                  
                
                if p['boundary_offshore'] == 're_circular':
                    Ct[:,0,0] =  -2                
                    Ct[:,-1,0] =  -2         
                
                if p['boundary_lateral'] == 'circular':
                    Ct[0,:,0] =  -1                
                    Ct[-1,:,0] =  -1
                
                if p['boundary_lateral'] == 're_circular':
                    Ct[0,:,0] =  -2                
                    Ct[-1,:,0] =  -2

                Ct, pickup = sweep(Ct, s['Cu'].copy(), s['mass'].copy(), self.dt, p['T'], s['ds'], s['dn'], s['us'], s['un'],w)

        qs = Ct * s['us'] 
        qn = Ct * s['un'] 
        q = np.hypot(qs, qn)


        return dict(Ct=Ct,
                    qs=qs,
                    qn=qn,
                    pickup=pickup,
                    w=w,
                    w_init=w_init,
                    w_air=w_air,
                    w_bed=w_bed,
                    q=q)
        
        
    def solve_steadystatepieter(self) -> dict:
        
        beta = 1. 
        
        l = self.l
        s = self.s
        p = self.p

        Ct = s['Ct'].copy()
        qs = s['qs'].copy()
        qn = s['qn'].copy()
        pickup = s['pickup'].copy()
        
        Ts = p['T']
        
        # compute transport weights for all sediment fractions
        w_init, w_air, w_bed = aeolis.transport.compute_weights(s, p)

        if self.t == 0.:
            # use initial guess for first time step
            w = p['grain_dist'].reshape((1,1,-1))
            w = w.repeat(p['ny']+1, axis=0)
            w = w.repeat(p['nx']+1, axis=1)
            return dict(w=w)
        else:
            w = w_init.copy()

        # set model state properties that are added to warnings and errors
        logprops = dict(minwind=s['uw'].min(),
                        maxdrop=(l['uw']-s['uw']).max(),
                        time=self.t,
                        dt=self.dt)
        
        nf = p['nfractions']
        
        ufs = np.zeros((p['ny']+1,p['nx']+2))
        ufn = np.zeros((p['ny']+2,p['nx']+1))    
        
        for i in range(nf): #loop over fractions
        
            #define velocity fluxes
            
            ufs[:,1:-1] = 0.5*s['us'][:,:-1,i] + 0.5*s['us'][:,1:,i]
            ufn[1:-1,:] = 0.5*s['un'][:-1,:,i] + 0.5*s['un'][1:,:,i]
            
            #boundary values
            ufs[:,0]  = s['us'][:,0,i]
            ufs[:,-1] = s['us'][:,-1,i]
            
            if p['boundary_lateral'] == 'circular':
                ufn[0,:] = 0.5*s['un'][0,:,i] + 0.5*s['un'][-1,:,i]
                ufn[-1,:] = ufn[0,:]
            else:
                ufn[0,:]  = s['un'][0,:,i]
                ufn[-1,:] = s['un'][-1,:,i]
        
            beta = abs(beta)
            if beta >= 1.:
                # define upwind direction
                ixfs = np.asarray(ufs >= 0., dtype=float)
                ixfn = np.asarray(ufn >= 0., dtype=float)
            else:
                # or centralizing weights
                ixfs = beta + np.zeros(ufs)
                ixfn = beta + np.zeros(ufn)

            # initialize matrix diagonals
            A0 = np.zeros(s['zb'].shape)
            Apx = np.zeros(s['zb'].shape)
            Ap1 = np.zeros(s['zb'].shape)
            Amx = np.zeros(s['zb'].shape)
            Am1 = np.zeros(s['zb'].shape)

            # populate matrix diagonals
            #A0         += s['dsdn'] / self.dt                                        #time derivative
            A0         += s['dsdn'] / Ts                                        #source term
            A0[:,1:]   -= s['dn'][:,1:]  * ufs[:,1:-1] * (1. - ixfs[:,1:-1])    #lower x-face
            Am1[:,1:]  -= s['dn'][:,1:]  * ufs[:,1:-1] *       ixfs[:,1:-1]     #lower x-face
            A0[:,:-1]  += s['dn'][:,:-1] * ufs[:,1:-1] *       ixfs[:,1:-1]     #upper x-face
            Ap1[:,:-1] += s['dn'][:,:-1] * ufs[:,1:-1] * (1. - ixfs[:,1:-1])    #upper x-face
            A0[1:,:]   -= s['ds'][1:,:]  * ufn[1:-1,:] * (1. - ixfn[1:-1,:])    #lower y-face
            Amx[1:,:]  -= s['ds'][1:,:]  * ufn[1:-1,:] *       ixfn[1:-1,:]     #lower y-face
            A0[:-1,:]  += s['ds'][:-1,:] * ufn[1:-1,:] *       ixfn[1:-1,:]     #upper y-face
            Apx[:-1,:] += s['ds'][:-1,:] * ufn[1:-1,:] * (1. - ixfn[1:-1,:])    #upper y-face
        
            # add boundaries
            # offshore boundary (i=0)

            if p['boundary_offshore'] == 'flux':
                #nothing to be done
                pass
            elif p['boundary_offshore'] == 'constant':
                #constant sediment concentration (Ct) in the air
                A0[:,0] = 1.
                Apx[:,0] = 0.
                Amx[:,0] = 0.
                Ap1[:,0] = 0.
                Am1[:,0] = 0.
            elif p['boundary_offshore'] == 'gradient':
                #remove the flux at the inner face of the cell
                A0[:,0]  -= s['dn'][:,0] * ufs[:,1] *       ixfs[:,1]           #upper x-face
                Ap1[:,0] -= s['dn'][:,0] * ufs[:,1] * (1. - ixfs[:,1])          #upper x-face
            elif p['boundary_offshore'] == 'circular':
                raise NotImplementedError('Cross-shore cricular boundary condition not yet implemented')
            else:
                raise ValueError('Unknown offshore boundary condition [%s]' % self.p['boundary_offshore'])

            #onshore boundary (i=nx)

            if p['boundary_onshore'] == 'flux':
                #nothing to be done
                pass
            elif p['boundary_onshore'] == 'constant':
                #constant sediment concentration (hC) in the air
                A0[:,-1] = 1.
                Apx[:,-1] = 0.
                Amx[:,-1] = 0.
                Ap1[:,-1] = 0.
                Am1[:,-1] = 0.
            elif p['boundary_onshore'] == 'gradient':
                #remove the flux at the inner face of the cell
                A0[:,-1]   += s['dn'][:,-1]  * ufs[:,-2]   * (1. - ixfs[:,-2])      #lower x-face
                Am1[:,-1]  += s['dn'][:,-1]  * ufs[:,-2]   *       ixfs[:,-2]       #lower x-face
            elif p['boundary_onshore'] == 'circular':
                raise NotImplementedError('Cross-shore cricular boundary condition not yet implemented')
            else:
                raise ValueError('Unknown offshore boundary condition [%s]' % self.p['boundary_onshore'])
        
            #lateral boundaries (j=0; j=ny)    

            if p['boundary_lateral'] == 'flux':
                #nothing to be done
                pass
            elif p['boundary_lateral'] == 'constant':
                #constant sediment concentration (hC) in the air
                A0[0,:] = 1.
                Apx[0,:] = 0.
                Amx[0,:] = 0.
                Ap1[0,:] = 0.
                Am1[0,:] = 0.
                A0[-1,:] = 1.
                Apx[-1,:] = 0.
                Amx[-1,:] = 0.
                Ap1[-1,:] = 0.
                Am1[-1,:] = 0.
            elif p['boundary_lateral'] == 'gradient':
                #remove the flux at the inner face of the cell
                A0[0,:]   -= s['ds'][0,:] * ufn[1,:]   *       ixfn[1,:]        #upper y-face
                Apx[0,:]  -= s['ds'][0,:] * ufn[1,:]   * (1. - ixfn[1,:])       #upper y-face
                A0[-1,:]  += s['ds'][-1,:] * ufn[-2,:] * (1. - ixfn[-2,:])      #lower y-face
                Amx[-1,:] += s['ds'][-1,:] * ufn[-2,:] *       ixfn[-2,:]       #lower y-face
            elif p['boundary_lateral'] == 'circular':   
                A0[0,:]   -= s['ds'][0,:]  * ufn[0,:]  * (1. - ixfn[0,:])       #lower y-face
                Amx[0,:]  -= s['ds'][0,:]  * ufn[0,:]  *       ixfn[0,:]        #lower y-face
                A0[-1,:]  += s['ds'][-1,:] * ufn[-1,:] *       ixfn[-1,:]       #upper y-face
                Apx[-1,:] += s['ds'][-1,:] * ufn[-1,:] * (1. - ixfn[-1,:])      #upper y-face
            else:
                raise ValueError('Unknown lateral boundary condition [%s]' % self.p['boundary_lateral'])
         
            # construct sparse matrix
            if p['ny'] > 0:
                j = p['nx']+1
                A = scipy.sparse.diags((Apx.flatten()[:j],
                                        Amx.flatten()[j:],
                                        Am1.flatten()[1:],
                                        A0.flatten(),
                                        Ap1.flatten()[:-1],
                                        Apx.flatten()[j:],
                                        Amx.flatten()[:j]),
                                       (-j*p['ny'],-j,-1,0,1,j,j*p['ny']), format='csr')
            else:
                j = p['nx']+1
                ny = 0
                A = scipy.sparse.diags((Am1.flatten()[1:],
                                        A0.flatten(),
                                        Ap1.flatten()[:-1]),
                                       (-1, 0, 1), format='csr')

            # solve transport for each fraction separately using latest
            # available weights


            # renormalize weights for all fractions equal or larger
            # than the current one such that the sum of all weights is
            # unity
            w = aeolis.transport.renormalize_weights(w, i)

            # iteratively find a solution of the linear system that
            # does not violate the availability of sediment in the bed
            for n in range(p['max_iter']):
                self._count('matrixsolve')
                
                # define upwind face value
                # sediment concentration
                Ctxfs_i = np.zeros(ufs.shape)
                Ctxfn_i = np.zeros(ufn.shape)
                
                Ctxfs_i[:,1:-1] = ixfs[:,1:-1] * Ct[:,:-1,i] \
                                    + (1. - ixfs[:,1:-1]) * Ct[:,1:,i] 
                Ctxfn_i[1:-1,:] = ixfn[1:-1,:] * Ct[:-1,:,i] \
                                    + (1. - ixfn[1:-1,:]) * Ct[1:,:,i] 

                if p['boundary_lateral'] == 'circular':
                    Ctxfn_i[0,:] = ixfn[0,:] * Ct[-1,:,i] \
                                    + (1. - ixfn[0,:]) *  Ct[0,:,i] 
                
                # calculate pickup
                D_i = s['dsdn'] / Ts * Ct[:,:,i]                                          
                A_i = s['dsdn'] / Ts * s['mass'][:,:,0,i] + D_i # Availability
                U_i = s['dsdn'] / Ts *  w[:,:,i] *  s['Cu'][:,:,i] 
                                            
                #deficit_i = E_i - A_i
                E_i= np.minimum(U_i, A_i)
                #pickup_i = E_i - D_i

                # create the right hand side of the linear system
                # sediment concentration
                yCt_i = np.zeros(s['zb'].shape)
                                
                yCt_i         += E_i - D_i                                      #source term
                yCt_i[:,1:]   += s['dn'][:,1:]  * ufs[:,1:-1] * Ctxfs_i[:,1:-1] #lower x-face
                yCt_i[:,:-1]  -= s['dn'][:,:-1] * ufs[:,1:-1] * Ctxfs_i[:,1:-1] #upper x-face
                yCt_i[1:,:]   += s['ds'][1:,:]  * ufn[1:-1,:] * Ctxfn_i[1:-1,:] #lower y-face
                yCt_i[:-1,:]  -= s['ds'][:-1,:] * ufn[1:-1,:] * Ctxfn_i[1:-1,:] #upper y-face
                
                # boundary conditions
                # offshore boundary (i=0)

                if p['boundary_offshore'] == 'flux':
                    yCt_i[:,0]  += s['dn'][:,0] * ufs[:,0] * s['Cu0'][:,0,i] * p['offshore_flux'] 
                elif p['boundary_offshore'] == 'constant':
                    #constant sediment concentration (Ct) in the air 
                    yCt_i[:,0]  = p['constant_offshore_flux']

                elif p['boundary_offshore'] == 'gradient':
                    #remove the flux at the inner face of the cell
                    yCt_i[:,0]  += s['dn'][:,1] * ufs[:,1] * Ctxfs_i[:,1] 

                elif p['boundary_offshore'] == 'circular':
                    raise NotImplementedError('Cross-shore cricular boundary condition not yet implemented')
                else:
                    raise ValueError('Unknown offshore boundary condition [%s]' % self.p['boundary_offshore'])
                    
                # onshore boundary (i=nx)

                if p['boundary_onshore'] == 'flux':
                    yCt_i[:,-1]  += s['dn'][:,-1]  * ufs[:,-1] * s['Cu0'][:,-1,i] * p['onshore_flux']

                elif p['boundary_onshore'] == 'constant':
                    #constant sediment concentration (Ct) in the air 
                    yCt_i[:,-1]  = p['constant_onshore_flux']

                elif p['boundary_onshore'] == 'gradient':
                    #remove the flux at the inner face of the cell
                    yCt_i[:,-1]  -= s['dn'][:,-2] * ufs[:,-2] * Ctxfs_i[:,-2] 

                elif p['boundary_onshore'] == 'circular':
                    raise NotImplementedError('Cross-shore cricular boundary condition not yet implemented')
                else:
                    raise ValueError('Unknown onshore boundary condition [%s]' % self.p['boundary_onshore'])
                    
                #lateral boundaries (j=0; j=ny)    

                if p['boundary_lateral'] == 'flux':
                    
                    yCt_i[0,:]   += s['ds'][0,:] * ufn[0,:]  * s['Cu0'][0,:,i] * p['lateral_flux'] #lower y-face
                    yCt_i[-1,:]  -= s['ds'][-1,:] * ufn[-1,:] * s['Cu0'][-1,:,i] * p['lateral_flux'] #upper y-face                    
                elif p['boundary_lateral'] == 'constant':
                    #constant sediment concentration (hC) in the air
                    yCt_i[0,:]  = 0.
                    yCt_i[-1,:] = 0.
                elif p['boundary_lateral'] == 'gradient':
                    #remove the flux at the inner face of the cell
                    yCt_i[-1,:] -= s['ds'][-2,:] * ufn[-2,:] * Ctxfn_i[-2,:] #lower y-face
                    yCt_i[0,:]  += s['ds'][1,:]  * ufn[1,:]  * Ctxfn_i[1,:]  #upper y-face
                elif p['boundary_lateral'] == 'circular':
                    yCt_i[0,:]  += s['ds'][0,:]  * ufn[0,:]  * Ctxfn_i[0,:]  #lower y-face
                    yCt_i[-1,:] -= s['ds'][-1,:] * ufn[-1,:] * Ctxfn_i[-1,:] #upper y-face
                else:
                    raise ValueError('Unknown lateral boundary condition [%s]' % self.p['boundary_lateral'])
                
                # print("ugs = %.*g" % (3,s['ugs'][10,10]))
                # print("ugn = %.*g" % (3,s['ugn'][10,10]))
                # print("%.*g" % (3,np.amax(np.absolute(y_i))))
                
                # solve system with current weights
                Ct_i = Ct[:,:,i].flatten()
                Ct_i += scipy.sparse.linalg.spsolve(A, yCt_i.flatten())
                Ct_i = prevent_tiny_negatives(Ct_i, p['max_error'])
                
                # check for negative values
                if Ct_i.min() < 0.:
                    ix = Ct_i < 0.
                    
#                    logger.warn(format_log('Removing negative concentrations',
#                                           nrcells=np.sum(ix),
#                                           fraction=i,
#                                           iteration=n,
#                                           minvalue=Ct_i.min(),
#                                           **logprops))

                    Ct_i[~ix] *= 1. + Ct_i[ix].sum() / Ct_i[~ix].sum()
                    Ct_i[ix] = 0.

                # determine pickup and deficit for current fraction
                Cu_i = s['Cu'][:,:,i].flatten()
                mass_i = s['mass'][:,:,0,i].flatten()
                w_i = w[:,:,i].flatten()
                Ts_i = Ts
                
                pickup_i = (w_i * Cu_i - Ct_i) / Ts_i * self.dt # Dit klopt niet! enkel geldig bij backward euler
                deficit_i = pickup_i - mass_i
                ix = (deficit_i > p['max_error']) \
                     & (w_i * Cu_i > 0.)

                pickup[:,:,i] = pickup_i.reshape(yCt_i.shape)
                Ct[:,:,i] = Ct_i.reshape(yCt_i.shape)
                     
                # quit the iteration if there is no deficit, otherwise
                # back-compute the maximum weight allowed to get zero
                # deficit for the current fraction and progress to
                # the next iteration step
                if not np.any(ix):
                    logger.debug(format_log('Iteration converged',
                                            steps=n,
                                            fraction=i,
                                            **logprops))
                    pickup_i = np.minimum(pickup_i, mass_i)
                    break
                else:
                    w_i[ix] = (mass_i[ix] * Ts_i / self.dt \
                               + Ct_i[ix]) / Cu_i[ix]
                    w[:,:,i] = w_i.reshape(yCt_i.shape)

            # throw warning if the maximum number of iterations was
            # reached
            if np.any(ix):
                logger.warn(format_log('Iteration not converged',
                                       nrcells=np.sum(ix),
                                       fraction=i,
                                       **logprops))
            
            # check for unexpected negative values
            if Ct_i.min() < 0:
                logger.warn(format_log('Negative concentrations',
                                       nrcells=np.sum(Ct_i<0.),
                                       fraction=i,
                                       minvalue=Ct_i.min(),
                                       **logprops))
            if w_i.min() < 0:
                logger.warn(format_log('Negative weights',
                                       nrcells=np.sum(w_i<0),
                                       fraction=i,
                                       minvalue=w_i.min(),
                                       **logprops))
        # end loop over frations

        # check if there are any cells where the sum of all weights is
        # smaller than unity. these cells are supply-limited for all
        # fractions. Log these events.
        ix = 1. - np.sum(w, axis=2) > p['max_error']
        if np.any(ix):
            self._count('supplylim')
#            logger.warn(format_log('Ran out of sediment',
#                                   nrcells=np.sum(ix),
#                                   minweight=np.sum(w, axis=-1).min(),
#                                   **logprops))
        qs = Ct * s['us'] 
        qn = Ct * s['un']
                    
        return dict(Ct=Ct,
                    qs=qs,
                    qn=qn,
                    pickup=pickup,
                    w=w,
                    w_init=w_init,
                    w_air=w_air,
                    w_bed=w_bed)
    
    
    def solve_pieter(self, alpha:float=.5, beta:float=1.) -> dict:
        '''Implements the explicit Euler forward, implicit Euler backward and semi-implicit Crank-Nicolson numerical schemes

        Determines weights of sediment fractions, sediment pickup and
        instantaneous sediment concentration. Returns a partial
        spatial grid dictionary that can be used to update the global
        spatial grid dictionary.

        Parameters
        ----------
        alpha : 
            Implicitness coefficient (0.0 for Euler forward, 1.0 for Euler backward or 0.5 for Crank-Nicolson, default=0.5)
        beta : float, optional
            Centralization coefficient (1.0 for upwind or 0.5 for centralized, default=1.0)

        Returns
        -------
            Partial spatial grid dictionary

        Examples
        --------
        >>> model.s.update(model.solve(alpha=1., beta=1.) # euler backward

        >>> model.s.update(model.solve(alpha=.5, beta=1.) # crank-nicolson

        See Also
        --------
        model.AeoLiS.euler_forward
        model.AeoLiS.euler_backward
        model.AeoLiS.crank_nicolson
        transport.compute_weights
        transport.renormalize_weights
        '''

        l = self.l
        s = self.s
        p = self.p

        Ct = s['Ct'].copy()
        qs = s['qs'].copy()
        qn = s['qn'].copy()
        pickup = s['pickup'].copy()
        
        Ts = p['T']
        
        # compute transport weights for all sediment fractions
        w_init, w_air, w_bed = aeolis.transport.compute_weights(s, p)

        if self.t == 0.:
            # use initial guess for first time step
            w = p['grain_dist'].reshape((1,1,-1))
            w = w.repeat(p['ny']+1, axis=0)
            w = w.repeat(p['nx']+1, axis=1)
            return dict(w=w)
        else:
            w = w_init.copy()

        # set model state properties that are added to warnings and errors
        logprops = dict(minwind=s['uw'].min(),
                        maxdrop=(l['uw']-s['uw']).max(),
                        time=self.t,
                        dt=self.dt)
        
        nf = p['nfractions']

        ufs = np.zeros((p['ny']+1,p['nx']+2))
        ufn = np.zeros((p['ny']+2,p['nx']+1))               
        
        for i in range(nf): #loop over fractions
        
            #define velocity fluxes
            ufs[:,1:-1] = 0.5*s['us'][:,:-1,i] + 0.5*s['us'][:,1:,i]
            ufn[1:-1,:] = 0.5*s['un'][:-1,:,i] + 0.5*s['un'][1:,:,i]
            
            #boundary values
            ufs[:,0]  = s['us'][:,0,i]
            ufs[:,-1] = s['us'][:,-1,i]
            
            if p['boundary_lateral'] == 'circular':
                ufn[0,:] = 0.5*s['un'][0,:,i] + 0.5*s['un'][-1,:,i]
                ufn[-1,:] = ufn[0,:]
            else:
                ufn[0,:]  = s['un'][0,:,i]
                ufn[-1,:] = s['un'][-1,:,i]    
        
            beta = abs(beta)
            if beta >= 1.:
                # define upwind direction
                ixfs = np.asarray(ufs >= 0., dtype=float)
                ixfn = np.asarray(ufn >= 0., dtype=float)
            else:
                # or centralizing weights
                ixfs = beta + np.zeros(ufs)
                ixfn = beta + np.zeros(ufn)

            # initialize matrix diagonals
            A0 = np.zeros(s['zb'].shape)
            Apx = np.zeros(s['zb'].shape)
            Ap1 = np.zeros(s['zb'].shape)
            Amx = np.zeros(s['zb'].shape)
            Am1 = np.zeros(s['zb'].shape)

            # populate matrix diagonals
            A0         += s['dsdn'] / self.dt                                        #time derivative
            A0         += s['dsdn'] / Ts                                     * alpha #source term
            A0[:,1:]   -= s['dn'][:,1:]  * ufs[:,1:-1] * (1. - ixfs[:,1:-1]) * alpha #lower x-face
            Am1[:,1:]  -= s['dn'][:,1:]  * ufs[:,1:-1] *       ixfs[:,1:-1]  * alpha #lower x-face
            A0[:,:-1]  += s['dn'][:,:-1] * ufs[:,1:-1] *       ixfs[:,1:-1]  * alpha #upper x-face
            Ap1[:,:-1] += s['dn'][:,:-1] * ufs[:,1:-1] * (1. - ixfs[:,1:-1]) * alpha #upper x-face
            A0[1:,:]   -= s['ds'][1:,:]  * ufn[1:-1,:] * (1. - ixfn[1:-1,:]) * alpha #lower y-face
            Amx[1:,:]  -= s['ds'][1:,:]  * ufn[1:-1,:] *       ixfn[1:-1,:]  * alpha #lower y-face
            A0[:-1,:]  += s['ds'][:-1,:] * ufn[1:-1,:] *       ixfn[1:-1,:]  * alpha #upper y-face
            Apx[:-1,:] += s['ds'][:-1,:] * ufn[1:-1,:] * (1. - ixfn[1:-1,:]) * alpha #upper y-face
        
            # add boundaries
            # offshore boundary (i=0)

            if p['boundary_offshore'] == 'flux':
                #nothing to be done
                pass
            elif p['boundary_offshore'] == 'constant':
                #constant sediment concentration (Ct) in the air
                A0[:,0] = 1.
                Apx[:,0] = 0.
                Amx[:,0] = 0.
                Ap1[:,0] = 0.
                Am1[:,0] = 0.
            elif p['boundary_offshore'] == 'gradient':
                #remove the flux at the inner face of the cell
                A0[:,0]  -= s['dn'][:,0] * ufs[:,1] *       ixfs[:,1]  * alpha #upper x-face
                Ap1[:,0] -= s['dn'][:,0] * ufs[:,1] * (1. - ixfs[:,1]) * alpha #upper x-face
            elif p['boundary_offshore'] == 'circular':
                raise NotImplementedError('Cross-shore cricular boundary condition not yet implemented')
            else:
                raise ValueError('Unknown offshore boundary condition [%s]' % self.p['boundary_offshore'])

            #onshore boundary (i=nx)

            if p['boundary_onshore'] == 'flux':
                #nothing to be done
                pass
            elif p['boundary_onshore'] == 'constant':
                #constant sediment concentration (hC) in the air
                A0[:,-1] = 1.
                Apx[:,-1] = 0.
                Amx[:,-1] = 0.
                Ap1[:,-1] = 0.
                Am1[:,-1] = 0.
            elif p['boundary_onshore'] == 'gradient':
                #remove the flux at the inner face of the cell
                A0[:,-1]   += s['dn'][:,-1]  * ufs[:,-2]   * (1. - ixfs[:,-2]) * alpha #lower x-face
                Am1[:,-1]  += s['dn'][:,-1]  * ufs[:,-2]   *       ixfs[:,-2]  * alpha #lower x-face
            elif p['boundary_onshore'] == 'circular':
                raise NotImplementedError('Cross-shore cricular boundary condition not yet implemented')
            else:
                raise ValueError('Unknown offshore boundary condition [%s]' % self.p['boundary_onshore'])
        
            #lateral boundaries (j=0; j=ny)    

            if p['boundary_lateral'] == 'flux':
                #nothing to be done
                pass
            elif p['boundary_lateral'] == 'constant':
                #constant sediment concentration (hC) in the air
                A0[0,:] = 1.
                Apx[0,:] = 0.
                Amx[0,:] = 0.
                Ap1[0,:] = 0.
                Am1[0,:] = 0.
                A0[-1,:] = 1.
                Apx[-1,:] = 0.
                Amx[-1,:] = 0.
                Ap1[-1,:] = 0.
                Am1[-1,:] = 0.
            elif p['boundary_lateral'] == 'gradient':
                #remove the flux at the inner face of the cell
                A0[0,:]   -= s['ds'][0,:] * ufn[1,:]   *       ixfn[1,:]   * alpha #upper y-face
                Apx[0,:]  -= s['ds'][0,:] * ufn[1,:]   * (1. - ixfn[1,:])  * alpha #upper y-face
                A0[-1,:]  += s['ds'][-1,:] * ufn[-2,:] * (1. - ixfn[-2,:]) * alpha #lower y-face
                Amx[-1,:] += s['ds'][-1,:] * ufn[-2,:] *       ixfn[-2,:]  * alpha #lower y-face
            elif p['boundary_lateral'] == 'circular':
                A0[0,:]   -= s['ds'][0,:]  * ufn[0,:]  * (1. - ixfn[0,:])  * alpha #lower y-face
                Amx[0,:]  -= s['ds'][0,:]  * ufn[0,:]  *       ixfn[0,:]   * alpha #lower y-face
                A0[-1,:]  += s['ds'][-1,:] * ufn[-1,:] *       ixfn[-1,:]  * alpha #upper y-face
                Apx[-1,:] += s['ds'][-1,:] * ufn[-1,:] * (1. - ixfn[-1,:]) * alpha #upper y-face
            else:
                raise ValueError('Unknown lateral boundary condition [%s]' % self.p['boundary_lateral'])
         
            # construct sparse matrix
            if p['ny'] > 0:
                j = p['nx']+1
                A = scipy.sparse.diags((Apx.flatten()[:j],
                                        Amx.flatten()[j:],
                                        Am1.flatten()[1:],
                                        A0.flatten(),
                                        Ap1.flatten()[:-1],
                                        Apx.flatten()[j:],
                                        Amx.flatten()[:j]),
                                       (-j*p['ny'],-j,-1,0,1,j,j*p['ny']), format='csr')
            else:
                A = scipy.sparse.diags((Am1.flatten()[1:],
                                        A0.flatten(),
                                        Ap1.flatten()[:-1]),
                                       (-1,0,1), format='csr')

            # solve transport for each fraction separately using latest
            # available weights
        
            # renormalize weights for all fractions equal or larger
            # than the current one such that the sum of all weights is
            # unity
            w = aeolis.transport.renormalize_weights(w, i)

            # iteratively find a solution of the linear system that
            # does not violate the availability of sediment in the bed
            for n in range(p['max_iter']):
                self._count('matrixsolve')
#                print("iteration nr = %d" % n)
                # define upwind face value
                # sediment concentration
                Ctxfs_i = np.zeros(ufs.shape)
                Ctxfn_i = np.zeros(ufn.shape)
                
                Ctxfs_i[:,1:-1] = ixfs[:,1:-1] * ( alpha * Ct[:,:-1,i] \
                                                  + (1. - alpha ) * l['Ct'][:,:-1,i] ) \
                    + (1. - ixfs[:,1:-1]) * ( alpha * Ct[:,1:,i] \
                                             + (1. - alpha ) * l['Ct'][:,1:,i] )
                Ctxfn_i[1:-1,:] = ixfn[1:-1,:] * (alpha * Ct[:-1,:,i] \
                                                  + (1. - alpha ) * l['Ct'][:-1,:,i] ) \
                    + (1. - ixfn[1:-1,:]) * ( alpha * Ct[1:,:,i] \
                                             + (1. - alpha ) * l['Ct'][1:,:,i] )
                    
                if p['boundary_lateral'] == 'circular':
                    Ctxfn_i[0,:] = ixfn[0,:] * (alpha * Ct[-1,:,i] \
                                                + (1. - alpha ) * l['Ct'][-1,:,i] ) \
                        + (1. - ixfn[0,:]) * ( alpha * Ct[0,:,i] \
                                               + (1. - alpha ) * l['Ct'][0,:,i] )
                    Ctxfn_i[-1,:] = Ctxfn_i[0,:]                   
                
                # calculate pickup
                D_i = s['dsdn'] / Ts * ( alpha * Ct[:,:,i]  \
                                            + (1. - alpha ) * l['Ct'][:,:,i] )
                A_i = s['dsdn'] / Ts * s['mass'][:,:,0,i] + D_i # Availability
                U_i = s['dsdn'] / Ts * ( w[:,:,i] * alpha * s['Cu'][:,:,i] \
                                            + (1. - alpha ) * l['w'][:,:,i] * l['Cu'][:,:,i] )
                #deficit_i = E_i - A_i
                E_i= np.minimum(U_i, A_i)
                #pickup_i = E_i - D_i

                # create the right hand side of the linear system
                # sediment concentration
                yCt_i = np.zeros(s['zb'].shape)
                yCt_i         -= s['dsdn'] / self.dt * ( Ct[:,:,i] \
                                                        - l['Ct'][:,:,i] )      #time derivative
                yCt_i         += E_i - D_i                                      #source term
                yCt_i[:,1:]   += s['dn'][:,1:]  * ufs[:,1:-1] * Ctxfs_i[:,1:-1] #lower x-face
                yCt_i[:,:-1]  -= s['dn'][:,:-1] * ufs[:,1:-1] * Ctxfs_i[:,1:-1] #upper x-face
                yCt_i[1:,:]   += s['ds'][1:,:]  * ufn[1:-1,:] * Ctxfn_i[1:-1,:] #lower y-face
                yCt_i[:-1,:]  -= s['ds'][:-1,:] * ufn[1:-1,:] * Ctxfn_i[1:-1,:] #upper y-face
                  
                # boundary conditions
                # offshore boundary (i=0)

                if p['boundary_offshore'] == 'flux':
                    yCt_i[:,0]  += s['dn'][:,0] * ufs[:,0] * s['Cu0'][:,0,i] * p['offshore_flux'] 

                elif p['boundary_offshore'] == 'constant':
                    #constant sediment concentration (Ct) in the air (for now = 0)
                    yCt_i[:,0]  = 0.

                elif p['boundary_offshore'] == 'gradient':
                    #remove the flux at the inner face of the cell
                    yCt_i[:,0]  += s['dn'][:,1] * ufs[:,1] * Ctxfs_i[:,1] #upper x-face

                elif p['boundary_offshore'] == 'circular':
                    raise NotImplementedError('Cross-shore cricular boundary condition not yet implemented')
                else:
                    raise ValueError('Unknown offshore boundary condition [%s]' % self.p['boundary_offshore'])
                    
                # onshore boundary (i=nx)

                if p['boundary_onshore'] == 'flux':
                    yCt_i[:,-1]  += s['dn'][:,-1]  * ufs[:,-1] * s['Cu0'][:,-1,i] * p['onshore_flux'] 

                elif p['boundary_onshore'] == 'constant':
                    #constant sediment concentration (Ct) in the air (for now = 0)
                    yCt_i[:,-1]  = 0.

                elif p['boundary_onshore'] == 'gradient':
                    #remove the flux at the inner face of the cell
                    yCt_i[:,-1]  -= s['dn'][:,-2] * ufs[:,-2] * Ctxfs_i[:,-2] #lower x-face

                elif p['boundary_onshore'] == 'circular':
                    raise NotImplementedError('Cross-shore cricular boundary condition not yet implemented')
                else:
                    raise ValueError('Unknown onshore boundary condition [%s]' % self.p['boundary_onshore'])
                    
                #lateral boundaries (j=0; j=ny)    

                if p['boundary_lateral'] == 'flux':
                    
                    yCt_i[0,:]  += s['ds'][0,:] * ufn[0,:]  * s['Cu0'][0,:,i] * p['lateral_flux'] #lower y-face
                    yCt_i[-1,:] -= s['ds'][-1,:] * ufn[-1,:] * s['Cu0'][-1,:,i] * p['lateral_flux'] #upper y-face
                    
                elif p['boundary_lateral'] == 'constant':
                    #constant sediment concentration (hC) in the air
                    yCt_i[0,:]  = 0.
                    yCt_i[-1,:] = 0.
                elif p['boundary_lateral'] == 'gradient':
                    #remove the flux at the inner face of the cell
                    yCt_i[-1,:] -= s['ds'][-2,:] * ufn[-2,:] * Ctxfn_i[-2,:] #lower y-face
                    yCt_i[0,:]  += s['ds'][1,:]  * ufn[1,:]  * Ctxfn_i[1,:]  #upper y-face
                elif p['boundary_lateral'] == 'circular':
                    yCt_i[0,:]  += s['ds'][0,:]  * ufn[0,:]  * Ctxfn_i[0,:]  #lower y-face
                    yCt_i[-1,:] -= s['ds'][-1,:] * ufn[-1,:] * Ctxfn_i[-1,:] #upper y-face
                else:
                    raise ValueError('Unknown lateral boundary condition [%s]' % self.p['boundary_lateral'])
                
                # print("ugs = %.*g" % (3,s['ugs'][10,10]))
                # print("ugn = %.*g" % (3,s['ugn'][10,10]))
                # print("%.*g" % (3,np.amax(np.absolute(y_i))))
                
                # solve system with current weights
                Ct_i = Ct[:,:,i].flatten()
                Ct_i += scipy.sparse.linalg.spsolve(A, yCt_i.flatten())
                Ct_i = prevent_tiny_negatives(Ct_i, p['max_error'])
                
                # check for negative values
                if Ct_i.min() < 0.:
                    ix = Ct_i < 0.
                    
#                    logger.warn(format_log('Removing negative concentrations',
#                                           nrcells=np.sum(ix),
#                                           fraction=i,
#                                           iteration=n,
#                                           minvalue=Ct_i.min(),
#                                           **logprops))
                    
                    if 0: #Ct_i[~ix].sum()>0.:
                        # compensate the negative concentrations by distributing them over the positives.
                        # I guess the idea is to conserve mass but it is not sure if this is needed, 
                        # mass continuity in the system is guaranteed by exchange with bed.
                        Ct_i[~ix] *= 1. + Ct_i[ix].sum() / Ct_i[~ix].sum()
                    Ct_i[ix] = 0.

                # determine pickup and deficit for current fraction
                Cu_i = s['Cu'][:,:,i].flatten()
                mass_i = s['mass'][:,:,0,i].flatten()
                w_i = w[:,:,i].flatten()
                Ts_i = Ts
                
                pickup_i = (w_i * Cu_i - Ct_i) / Ts_i * self.dt # Dit klopt niet! enkel geldig bij backward euler
                deficit_i = pickup_i - mass_i
                ix = (deficit_i > p['max_error']) \
                     & (w_i * Cu_i > 0.)

                pickup[:,:,i] = pickup_i.reshape(yCt_i.shape)
                Ct[:,:,i] = Ct_i.reshape(yCt_i.shape)
                     
                # quit the iteration if there is no deficit, otherwise
                # back-compute the maximum weight allowed to get zero
                # deficit for the current fraction and progress to
                # the next iteration step
                if not np.any(ix):
                    logger.debug(format_log('Iteration converged',
                                            steps=n,
                                            fraction=i,
                                            **logprops))
                    pickup_i = np.minimum(pickup_i, mass_i)
                    break
                else:
                    w_i[ix] = (mass_i[ix] * Ts_i / self.dt \
                               + Ct_i[ix]) / Cu_i[ix]
                    w[:,:,i] = w_i.reshape(yCt_i.shape)

            # throw warning if the maximum number of iterations was
            # reached
            
            if np.any(ix):
                logger.warn(format_log('Iteration not converged',
                                       nrcells=np.sum(ix),
                                       fraction=i,
                                       **logprops))
            
            if 0: #let's disable these warnings
                # check for unexpected negative values
                if Ct_i.min() < 0:
                    logger.warn(format_log('Negative concentrations',
                                        nrcells=np.sum(Ct_i<0.),
                                        fraction=i,
                                        minvalue=Ct_i.min(),
                                        **logprops))
                if w_i.min() < 0:
                    logger.warn(format_log('Negative weights',
                                        nrcells=np.sum(w_i<0),
                                        fraction=i,
                                        minvalue=w_i.min(),
                                        **logprops))
        # end loop over frations

        # check if there are any cells where the sum of all weights is
        # smaller than unity. these cells are supply-limited for all
        # fractions. Log these events.
        ix = 1. - np.sum(w, axis=2) > p['max_error']
        if np.any(ix):
            self._count('supplylim')
#            logger.warn(format_log('Ran out of sediment',
#                                   nrcells=np.sum(ix),
#                                   minweight=np.sum(w, axis=-1).min(),
#                                   **logprops))

        qs = Ct * s['us'] 
        qn = Ct * s['un']
        qs = Ct * s['us'] 
        qn = Ct * s['un'] 
        q = np.hypot(qs, qn)
        
                    
        return dict(Ct=Ct,
                    qs=qs,
                    qn=qn,
                    pickup=pickup,
                    w=w,
                    w_init=w_init,
                    w_air=w_air,
                    w_bed=w_bed,
                    q=q)        

    def get_count(self, name:str) -> int:
        '''Get counter value

        Parameters
        ----------
        name : str
            Name of counter

        '''

        if name in self.c:
            return self.c[name]
        else:
            return 0


    def _count(self, name:str, n:int=1) -> None:
        '''Increase counter

        Parameters
        ----------
        name : 
            Name of counter
        n : optional
            Increment of counter (default: 1)
        '''

        if name not in self.c:
            self.c[name] = 0
        self.c[name] += n


    def _dims2shape(self, dims) -> Tuple:
        '''Converts named dimensions to numbered shape

        Supports only dimension names that can be found in the model
        parameters dictionary. The dimensions ``nx`` and ``ny`` are
        increased by one, so they match the size of the spatial grids
        rather than the number of spatial cells in the model.

        Parameters
        ----------
        dims : iterable
            Iterable with strings specifying dimension names

        Returns
        -------
            Shape of spatial grid

        '''

        shape = []
        for dim in dims:
            shape.append(self.p[dim])
            if dim in ['nx', 'ny']:
                shape[-1] += 1
        return tuple(shape)


    @staticmethod
    def dimensions(var:str=None) -> Union[Tuple, dict]:
        '''Static method that returns named dimensions of all spatial grids

        Parameters
        ----------
        var : optional
            Name of spatial grid

        Returns
        -------
        tuple or dict
            Tuple with named dimensions of requested spatial grid or
            dictionary with all named dimensions of all spatial
            grids. Returns nothing if requested spatial grid is not
            defined.
        '''

        dims = {s:d
                for d, states in aeolis.constants.MODEL_STATE.items()
                for s in states}

        if var is not None:
            if var in dims:
                return dims[var]
            else:
                return None
        else:
            return dims


class AeoLiSRunner(AeoLiS):
    '''AeoLiS model runner class

    This runner class is a convenience class for the BMI-compatible
    AeoLiS model class (:class:`~model.AeoLiS()`). It implements a
    time loop, a progress indicator and netCDF4 output. It also
    provides the definition of a callback function that can be used to
    interact with the AeoLiS model during runtime.

    The command-line function ``aeolis`` is available that uses this
    class to start an AeoLiS model run.

    Examples
    --------
    >>> # run with default settings
    ... AeoLiSRunner().run()

    >>> AeoLiSRunner(configfile='aeolis.txt').run()

    >>> model = AeoLiSRunner(configfile='aeolis.txt')
    >>> model.run(callback=lambda model: model.set_var('zb', zb))

    >>> model.run(callback='bar.py:add_bar')

    See Also
    --------
    console.aeolis
    '''

    def __init__(self, configfile:str='aeolis.txt') -> None:
        '''Initialize class

        Reads model configuration file without parsing all referenced
        files for the progress indicator and netCDF output. If no
        configuration file is given, the default settings are used.

        Parameters
        ----------
        configfile : str, optional
            Model configuration file. See :func:`~inout.read_configfile()`.

        '''

        super(AeoLiSRunner, self).__init__(configfile=configfile)

        self.t0 = None
        self.tout = 0.
        self.tlog = 0.
        self.plog = -1.
        self.trestart = 0.

        self.n = 0 # time step counter
        self.o = {} # output stats

        self.changed = False
        self.cwd = None

        self.set_configfile(configfile)
        if os.path.exists(self.configfile):
            self.p = aeolis.inout.read_configfile(self.configfile, parse_files=False)
            self.changed = False
        elif self.configfile.upper() == 'DEFAULT':
            self.changed = True
            self.configfile = os.path.abspath('aeolis.txt')
            self.p = aeolis.constants.DEFAULT_CONFIG

            # add default profile and time series
            self.p.update(dict(nx         = 99,
                               ny         = 0,
                               xgrid_file = np.arange(0.,100.,1.),
                               bed_file   = np.linspace(-5.,5.,100.),
                               wind_file  = np.asarray([[0.,10.,0.],
                                                        [3601.,10.,0.]])))
        else:
            logger.log_and_raise('Configuration file not found [%s]' % self.configfile, exc=IOError)


    def run(self, callback=None, restartfile:str=None) -> None:
        '''Start model time loop

        Changes current working directory to the model directory,
        prints model configuration parameters and progress indicator
        to the screen, writes netCDF4 output and calls a callback
        function upon request.

        Parameters
        ----------
        callback : str or function
            The callback function is called at the start of every
            single time step and takes the AeoLiS model object as
            input. The callback function can be used to interact with
            the model during simulation (e.g. update the bed with new
            measurements). See for syntax
            :func:`~model.AeoLiSRunner.parse_callback()`.
        restartfile : str
            Path to previously written restartfile. The model state is
            loaded from this file after initialization of the model.

        See Also
        --------
        model.AeoLiSRunner.parse_callback

        '''

        # http://www.patorjk.com/software/taag/
        # font: Colossal

        if (logger.hasHandlers()):
            logger.handlers.clear()
        logger.setLevel(logging.DEBUG)

        # initialize file logger
        filehandler = logging.FileHandler('%s.log' % os.path.splitext(self.configfile)[0], mode='w')
        filehandler.setLevel(logging.INFO)
        filehandler.setFormatter(logging.Formatter('%(asctime)-15s %(name)-8s %(levelname)-8s %(message)s'))
        logger.addHandler(filehandler)

        # initialize console logger
        streamhandler = logging.StreamHandler()
        streamhandler.setLevel(20)
        streamhandler.setFormatter(StreamFormatter())
        logger.addHandler(streamhandler)


        logger.info('**********************************************************')
        logger.info('                                                          ')
        logger.info('         d8888                   888      d8b  .d8888b.   ')
        logger.info('        d88888                   888      Y8P d88P  Y88b  ')
        logger.info('       d88P888                   888          Y88b.       ')
        logger.info('      d88P 888  .d88b.   .d88b.  888      888  "Y888b.    ')
        logger.info('     d88P  888 d8P  Y8b d88""88b 888      888     "Y88b.  ')
        logger.info('    d88P   888 88888888 888  888 888      888       "888  ')
        logger.info('   d8888888888 Y8b.     Y88..88P 888      888 Y88b  d88P  ')
        logger.info('  d88P     888  "Y8888   "Y88P"  88888888 888  "Y8888P"   ')
        logger.info('                                                          ')
        logger.info('  Version:  %-45s' % __version__)
        # logger.info('  Git hash: %-45s' % __gitversion__) # commenting for now until we implement pyproject.toml
        logger.info('                                                          ')

        # set working directory
        fpath, fname = os.path.split(self.configfile)
        if fpath != os.getcwd():
            self.cwd = os.getcwd()
            os.chdir(fpath)
            logger.info('Changed working directory to: %s\n', fpath)

        # print settings
        self.print_params()

        # write settings
        self.write_params()

        # parse callback
        if callback is not None:
            callback = self.parse_callback(callback)
        else:
            callback = self.parse_callback(self.p['callback'])
        if callback is not None:
            logger.info('Applying callback function: %s()\n', callback.__name__)

        # initialize model
        self.initialize()
        self.load_hotstartfiles()

        # load restartfile
        if self.load_restartfile(restartfile):
            logger.info('Loaded model state from restart file: %s\n', restartfile)

        # start model loop
        self.t0 = time.time()
        self.output_write()
        while self.t <= self.p['tstop']:
            if callback is not None:
                callback(self)
            self.update()
            self.output_write()
            self.print_progress()

        # finalize model
        self.finalize()

        self.print_stats()

        if self.cwd is not None:
            os.chdir(self.cwd)

        logging.shutdown()

        logging.shutdown()


    def set_configfile(self, configfile:str) -> None:
        '''Set model configuration file name'''

        self.changed = False
        if os.path.exists(configfile):
            self.configfile = os.path.abspath(configfile)
        else:
            self.configfile = configfile


    def set_params(self, **kwargs) -> None:
        '''Set model configuration parameters'''

        if len(kwargs) > 0:
            self.changed = True
            self.p.update(kwargs)


    def get_statistic(self, var:str, stat:str='avg') -> Union[np.ndarray, None]:
        '''Return statistic of spatial grid

        Parameters
        ----------
        var : str
            Name of spatial grid
        stat : str
            Name of statistic (avg, sum, var, min or max)

        Returns
        -------
            Statistic of spatial grid
        '''

        if stat in ['min', 'max']:
            return self.o[var][stat]
        elif stat == 'sum':
            return self.o[var][stat]*self.dt
        elif stat == 'avg':
            if self.n > 0:
                return self.o[var]['sum'] / self.n
            else:
                return np.zeros(self.o[var]['sum'].shape)
        elif stat == 'var':
            if self.n > 1:
                return (self.o[var]['var'] - self.o[var]['sum']**2 / self.n) \
                    / (self.n - 1)
            else:
                return np.zeros(self.o[var]['var'].shape)
        else:
            return None


    def get_var(self, var:str, clear:bool=True) -> Union[np.ndarray, int, float, str, list]:
        '''Returns spatial grid, statistic or model configuration parameter

        Overloads the :func:`~model.AeoLiS.get_var()` function and
        extends it with the functionality to return statistics on
        spatial grids by adding a postfix to the variable name
        (e.g. Ct_avg). Supported statistics are avg, sum, var, min and
        max.

        Parameters
        ----------
        var : str
            Name of spatial grid or model configuration
            parameter. Spatial grid name can be extended with a
            postfix to request a statistic (_avg, _sum, _var, _min or
            _max).
        clear : bool
            Clear output statistics afterwards.

        Returns
        -------
        np.ndarray or int, float, str or list
            Spatial grid, statistic or model configuration parameter

        Examples
        --------
        >>> # returns average sediment concentration
        ... model.get_var('Ct_avg')

        >>> # returns variance in wave height
        ... model.get_var('Hs_var')

        See Also
        --------
        model.AeoLiS.get_var
        '''

        self.clear = clear

        if '_' in var:
            var, stat = var.split('_')
            if var in self.o:
                return self.get_statistic(var, stat)

        # TODO: delete in future releases
        if '.' in var:
            warnings.warn('The use of "%s" is deprecated, use '
                             '"%s" instead.' % (var, var.replace('.','_')), DeprecationWarning)
            var, stat = var.split('.')
            if var in self.o:
                return self.get_statistic(var, stat)

        return super(AeoLiSRunner, self).get_var(var)


    def initialize(self) -> None:
        '''Initialize model

        Overloads the :func:`~model.AeoLiS.initialize()` function, but
        also initializes output statistics.

        '''

        super(AeoLiSRunner, self).initialize()
        self.output_init()

    def update(self, dt:float=-1) -> None:
        '''Time stepping function

        Overloads the :func:`~model.AeoLiS.update()` function,
        but also updates output statistics and clears output
        statistics upon request.

        Parameters
        ----------
        dt : float, optional
            Time step in seconds.

        '''

        if self.clear or self.dt < -1:
            self.output_clear()
            self.clear = False

        super(AeoLiSRunner, self).update(dt=dt)
        self.output_update()


    def write_params(self) -> None:
        '''Write updated model configuration to configuration file

        Creates a backup in case the model configration file already
        exists.

        See Also
        --------
        inout.backup

        '''

        if self.changed:
            aeolis.inout.backup(self.configfile)
            aeolis.inout.write_configfile(self.configfile, self.p)
            self.changed = False


    def output_init(self) -> None:

        '''Initialize netCDF4 output file and output statistics dictionary'''

        self.p['output_vars'] = makeiterable(self.p['output_vars'])
        self.p['output_types'] = makeiterable(self.p['output_types'])

        # determine unique combinations of variables and types
        self.p['_output_vars'] = {}
        for var in self.p['output_vars']:
            if '_' in var:
                var0, ext = var.split('_')
            # TODO: delete in future release
            elif '.' in var:
                warnings.warn('The use of "%s" is deprecated, use '
                                 '"%s" instead.' % (var, var.replace('.','_')), DeprecationWarning)
                var0, ext = var.split('.')
            else:
                var0, ext = var, None
            if var0 not in self.p['_output_vars']:
                self.p['_output_vars'][var0] = []
            if ext not in self.p['_output_vars'][var0]:
                self.p['_output_vars'][var0].append(ext)
            for ext in self.p['output_types']:
                if ext not in self.p['_output_vars'][var0]:
                    self.p['_output_vars'][var0].append(ext)

        aeolis.netcdf.initialize(self.p['output_file'],
                                 self.p['_output_vars'],
                                 self.s,
                                 self.p,
                                 self.dimensions())

        self.output_clear()


    def output_clear(self) -> None:
        '''Clears output statistics dictionary

        Creates a matrix for minimum, maximum, variance and summed
        values for each output variable and sets the time step counter
        to zero.

        '''

        for k in self.p['_output_vars'].keys():
            s = self.get_var_shape(k)
            self.o[k] = dict(min=np.zeros(s) + np.inf,
                             max=np.zeros(s) - np.inf,
                             var=np.zeros(s),
                             sum=np.zeros(s))

        self.n = 0


    def output_update(self) -> None:
        '''Updates output statistics dictionary

        Updates matrices with minimum, maximum, variance and summed
        values for each output variable with current spatial grid
        values and increases time step counter with one.

        '''

        for k, exts in self.p['_output_vars'].items():
            v = self.get_var(k, clear=False).copy()
            if 'min' in exts:
                self.o[k]['min'] = np.minimum(self.o[k]['min'], v)
            if 'max' in exts:
                self.o[k]['max'] = np.maximum(self.o[k]['max'], v)
            if 'sum' in exts or 'avg' in exts or 'var' in exts:
                self.o[k]['sum'] = self.o[k]['sum'] + v
            if 'var' in exts:
                self.o[k]['var'] = self.o[k]['var'] + v**2

        #also update the q variable here
        # self.s['q'] = np.hypot(self.s['qs'], self.s['qn'])

        self.n += 1


    def output_write(self) -> None:
        '''Appends output to netCDF4 output file

        If the time since the last output is equal or larger than the
        set output interval, append current output to the netCDF4
        output file. Computes the average and variance values based on
        available output statistics and clear output statistics
        dictionary.
        '''

        if self.t - self.tout >= self.p['output_times'] or self.t == 0.:

            variables = {}
            variables['time'] = self.t
            for k, exts in self.p['_output_vars'].items():
                for ext in exts:
                    if ext is None:
                        variables[k] = self.get_var(k, clear=False).copy()
                    else:
                        variables['%s_%s' % (k, ext)] = self.get_statistic(k, ext)

            aeolis.netcdf.append(self.p['output_file'], variables)

            self.output_clear()
            self.tout = self.t

        if self.p['restart'] and self.t - self.trestart >= self.p['restart']:
            self.dump_restartfile()
            self.trestart = self.t


    def load_hotstartfiles(self) -> None:
        '''Load model state from hotstart files

        Hotstart files are plain text representations of model state
        variables that can be used to hotstart the (partial) model
        state. Hotstart files should have the name of the model state
        variable it contains and have the extension
        `.hotstart`. Hotstart files differ from restart files in that
        restart files contain entire model states and are pickled
        Python objects.

        See Also
        --------
        model.AeoLiSRunner.load_restartfile

        '''

        for fname in glob.glob('*.hotstart'):
            var = os.path.splitext(fname)[0]
            if var in self.s.keys():
                shp = self.s[var].shape
                self.s[var] = np.loadtxt(fname).reshape(shp)
                self.s.set_immutable(var)
                logger.info('Loaded "%s" from hotstart file.' % var)
            else:
                logger.warning('Unrecognized hotstart file [%s]' % fname)


    def load_restartfile(self, restartfile:str) -> bool:
        '''Load model state from restart file

        Parameters
        ----------
        restartfile : str
            Path to previously written restartfile.
        
        Returns
        -------
            True if model state from restartfile is loaded successfully
        '''

        if restartfile:
            if os.path.exists(restartfile):
                # Load model state
                with open(restartfile, 'r') as fp:
                    state = pickle.load(fp)

                    self.t = state['t']
                    self.p = state['p']
                    self.s = state['s']
                    self.l = state['l']
                    self.c = state['c']

                    self.trestart = self.t

                    return True
            else:
                logger.log_and_raise('Restart file not found [%s]' % restartfile, exc=IOError)
                return False


    def dump_restartfile(self) -> None:
        '''Dumps model state to restart file'''

        restartfile = '%s.r%010d' % (os.path.splitext(self.p['output_file'])[0], int(self.t))
        with open(restartfile, 'w') as fp:
            pickle.dump({'t':self.t,
                         'p':self.p,
                         's':self.s,
                         'l':self.l,
                         'c':self.c}, fp)

        logger.info('Written restart file [%s]' % restartfile)


    def parse_callback(self, callback):
        '''Parses callback definition and returns function

        The callback function can be specified in two formats:

        - As a native Python function
        - As a string refering to a Python script and function,
          separated by a colon (e.g. ``example/callback.py:function``)

        Parameters
        ----------
        callback : str or function
            Callback definition

        Returns
        -------
        function
            Python callback function
        '''

        if isinstance(callback, str):
            if ':' in callback:
                fname, func = callback.split(':')
                if os.path.exists(fname):
                    mod = importlib.machinery.SourceFileLoader('callback', fname).load_module()
                    if hasattr(mod, func):
                        return getattr(mod, func)
                else:
                    logger.error('Invalid callback definition [%s]', callback)
                    raise IOError('Check definition in input file [%s]' % callback)
                #     logger.error(f"Callback function not found: {fname}. Check the reference to the callback in the input file")
        elif hasattr(callback, '__call__'):
            return callback
        elif callback is None:
            return callback
        else:
            return None


    def print_progress(self, fraction:float=.01, min_interval:float=1., max_interval:float=60.) -> None:
        '''Print progress to screen

        Parameters
        ----------
        fraction : optional
            Fraction of simulation at which to print progress (default: 1%)
        min_interval :  optional
            Minimum time in seconds between subsequent progress prints (default: 1s)
        max_interval : optional
            Maximum time in seconds between subsequent progress prints (default: 60s)
        '''

        p = (self.t-self.p['tstart']) / (self.p['tstop']-self.p['tstart'])
        pr = np.ceil(p/fraction)*fraction
        t = time.time()
        interval = t - self.tlog
                
        if self.get_count('time') == 1:
            logger.info('        Time elapsed / Total time / Time remaining / Average Timestep')
            self.dt_array = []
            
        

        self.dt_array.append(self.dt)

        if (np.mod(p, fraction) < .01 and self.plog != pr) or interval > max_interval:
            t1 = timedelta(0, round(t-self.t0))
            t2 = timedelta(0, round((t-self.t0)/p))
            t3 = timedelta(0, round((t-self.t0)*(1.-p)/p))
            dt_avg = np.average(self.dt_array)
            logger.info('%05.1f%%  %12s / %10s / %14s / %0.1f' % (p * 100., t1, t2, t3, dt_avg))
            self.tlog = time.time()
            self.plog = pr
            self.dt_array = []
            


    def print_params(self) -> None:
        '''Print model configuration parameters to screen'''

        maxl = np.max([len(par) for par in self.p.keys()])
        fmt1 = '  %%-%ds = %%s' % maxl
        fmt2 = '  %%-%ds   %%s' % maxl

        logger.info('**********************************************************')
        logger.info('PARAMETER SETTINGS                                        ')
        logger.info('**********************************************************')

        for par, val in sorted(self.p.items()):
            if isiterable(val):
                if par.endswith('_file'):
                    logger.info(fmt1 % (par, '%s.txt' % par.replace('_file', '')))
                elif len(val) > 0:
                    logger.info(fmt1 % (par, aeolis.inout.print_value(val[0])))
                    for v in val[1:]:
                        logger.info(fmt2 % ('', aeolis.inout.print_value(v)))
                else:
                    logger.info(fmt1 % (par, ''))
            else:
                logger.info(fmt1 % (par, aeolis.inout.print_value(val)))

        logger.info('**********************************************************')
        logger.info('')


    def print_stats(self) -> None:
        '''Print model run statistics to screen'''

        n_time = self.get_count('time')
        n_matrixsolve = self.get_count('matrixsolve')
        n_supplylim = self.get_count('supplylim')

        logger.info('')
        logger.info('**********************************************************')

        fmt = '%-20s : %s'
        logger.info(fmt % ('# time steps', aeolis.inout.print_value(n_time)))
        logger.info(fmt % ('# matrix solves', aeolis.inout.print_value(n_matrixsolve)))
        logger.info(fmt % ('# supply lim', aeolis.inout.print_value(n_supplylim)))
        logger.info(fmt % ('avg. solves per step',
                           aeolis.inout.print_value(float(n_matrixsolve) / n_time)))
        logger.info(fmt % ('avg. time step',
                           aeolis.inout.print_value(float(self.p['tstop']) / n_time)))

        logger.info('**********************************************************')
        logger.info('')


class WindGenerator():
    '''Wind velocity time series generator

    Generates a random wind velocity time series with given mean and
    maximum wind speed, duration and time resolution. The wind
    velocity time series is generated using a Markov Chain Monte Carlo
    (MCMC) approach based on a Weibull distribution. The wind time
    series can be written to an AeoLiS-compatible wind input file
    assuming a constant wind direction of zero degrees.

    The command-line function ``aeolis-wind`` is available that uses
    this class to generate AeoLiS wind input files.

    Examples
    --------
    >>> wind = WindGenerator(mean_speed=10.).generate(duration=24*3600.)
    >>> wind.write_time_series('wind.txt')
    >>> wind.plot()
    >>> wind.hist()

    See Also
    --------
    console.wind

    '''

    # source:
    # http://www.lutralutra.co.uk/2012/07/02/simulating-a-wind-speed-time-series-in-python/

    
    def __init__(self,
                 mean_speed:float=9.0,
                 max_speed:float=30.0,
                 dt:float=60., # size of time step
                 n_states:int=30,
                 shape:float=2.,
                 scale:float=2.) -> None:

        self.mean_speed=mean_speed
        self.max_speed=max_speed
        self.n_states=n_states

        self.t=0.
        self.dt=dt

        # setup matrix
        n_rows = n_columns = n_states
        self.bin_size = float(max_speed)/n_states

        # weibull parameters
        weib_shape=shape
        weib_scale=scale*float(mean_speed)/np.sqrt(np.pi);

        # wind speed bins
        self.bins = np.arange(self.bin_size/2.0,
                              float(max_speed) + self.bin_size/2.0,
                              self.bin_size)

        # distribution of probabilities, normalised
        fdpWind = self.weibullpdf(self.bins, weib_scale, weib_shape)
        fdpWind = fdpWind / sum(fdpWind)

        # decreasing function
        G = np.empty((n_rows, n_columns,))
        for x in range(n_rows):
            for y in range(n_columns):
                G[x][y] = 2.0**float(-abs(x-y))

        # initial value of the P matrix
        P0 = np.diag(fdpWind)

        # initital value of the p vector
        p0 = fdpWind

        P, p = P0, p0
        rmse = np.inf
        while rmse > 1e-10:
            pp = p
            r = self.matmult4(P,self.matmult4(G,p))
            r = r/sum(r)
            p = p+0.5*(p0-r)
            P = np.diag(p)

            rmse = np.sqrt(np.mean((p - pp)**2))

        N=np.diag([1.0/i for i in self.matmult4(G,p)])
        MTM=self.matmult4(N,self.matmult4(G,P))
        self.MTMcum = np.cumsum(MTM,1)


    def __getitem__(self, s:int)-> np.ndarray:
        """
        Retrieves item from list of wind speeds as an array 
        
        Parameters:
        -----------
        s: 
            index in a list
        """

        return np.asarray(self.wind_speeds[s])


    def generate(self, duration:float=3600.):
        """Generates wind velocities

        Parameters:
        -----------
            duration: duration of time series in seconds
        """

        # initialise series
        self.state = 0
        self.states = []
        self.wind_speeds = []
        self.randoms1 = []
        self.randoms2 = []

        self.update()
        self.t = 0.

        while self.t < duration:
            self.update()

        return self


    def update(self) -> None:
        """
        Generates wind speeds base on ramdom uniform distribution 
        """

        r1 = np.random.uniform(0,1)
        r2 = np.random.uniform(0,1)

        self.randoms1.append(r1)
        self.randoms2.append(r2)

        self.state = next(j for j,v in enumerate(self.MTMcum[self.state]) if v > r1)
        self.states.append(self.state)

        u = np.maximum(0., self.bins[self.state] - 0.5 + r2 * self.bin_size)
        self.wind_speeds.append(u)

        self.t += self.dt


    def get_time_series(self) -> Any:
        """
        returns time series
        """
        u = np.asarray(self.wind_speeds)
        t = np.arange(len(u)) * self.dt

        return t, u


    def write_time_series(self, fname:str) -> None:
        """
        Saves time series to file

        Parameters:
        -----------
            fname: file name to write time series to
        """
        t, u = self.get_time_series()
        M = np.concatenate((np.asmatrix(t),
                            np.asmatrix(u),
                            np.zeros((1, len(t)))), axis=0).T

        np.savetxt(fname, M)


    def plot(self):
        """
        Creates plot of wind speed vs. time
        
        Returns:

        """
        t, u = self.get_time_series()

        fig, axs = plt.subplots(figsize=(10,4))
        axs.plot(t, u

    , '-k')
        axs.set_ylabel('wind speed [m/s]')
        axs.set_xlabel('time [s]')
        axs.set_xlim((0, np.max(t)))
        axs.grid()

        return fig, axs


    def hist(self):
        """
        Creates histogram for wind speeds

        Returns:

        """
        fig, axs = plt.subplots(figsize=(10,4))
        axs.hist(self.wind_speeds, bins=self.bins, normed=True, color='k')
        axs.set_xlabel('wind speed [m/s]')
        axs.set_ylabel('occurence [-]')
        axs.grid()

        return fig, axs


    @staticmethod
    def weibullpdf(data, scale, shape):
        # TODO: add description
        return [(shape/scale)
                * ((x/scale)**(shape-1))
                * np.exp(-1*(x/scale)**shape)
                for x in data]


    @staticmethod
    def matmult4(m, v):
        # TODO: add description
        return [reduce(operator.add, map(operator.mul,r,v)) for r in m]
