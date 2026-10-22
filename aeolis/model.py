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

import warnings
warnings.filterwarnings(
    "ignore",
    message="pkg_resources is deprecated as an API",
    category=UserWarning
)

import os
import importlib.machinery
import importlib.metadata
import time
import glob
import logging
import warnings
import operator
import numpy as np
import pickle
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
import aeolis.grass
import aeolis.fences
import aeolis.gridparams
import aeolis.zeta
import aeolis.onerow

from aeolis.advection import (
    solve_EF,
    solve_SS,
    solve_EB,
)

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


# Solvers for the advection-pickup equation and the legacy names that map onto them
SOLVERS = ('steadystate', 'euler_forward', 'euler_backward', 'crank_nicolson')
LEGACY_SOLVERS = ('trunk', 'pieter', 'steadystatepieter')


def resolve_solver(p:dict) -> str:
    '''Return the solver name to use for a model configuration.

    The current solver names are listed in :data:`SOLVERS`. Older
    configurations combine ``solver = trunk | pieter`` with a separate
    ``scheme = euler_forward | euler_backward | crank_nicolson`` setting,
    or use ``solver = steadystatepieter``. These are mapped onto the
    current solvers with a deprecation warning, so that existing input
    files keep working.

    Parameters
    ----------
    p : dict
        Model configuration parameters

    Returns
    -------
    str
        One of :data:`SOLVERS`

    '''

    solver = str(p['solver']).lower()
    if solver in SOLVERS:
        return solver

    if solver == 'steadystatepieter':
        resolved = 'steadystate'
    elif solver in ('trunk', 'pieter'):
        resolved = str(p.get('scheme') or 'euler_backward').lower()
        if resolved not in SOLVERS:
            logger.log_and_raise('Unknown scheme [%s] for solver [%s]' % (resolved, solver), exc=ValueError)
    else:
        logger.log_and_raise('Unknown solver [%s], use one of: %s' % (solver, ', '.join(SOLVERS)),
                             exc=ValueError)

    msg = ('Solver "%s" (scheme "%s") is deprecated and is mapped to solver = %s. '
           'Please update the configuration file.' % (solver, p.get('scheme'), resolved))
    warnings.warn(msg, DeprecationWarning)
    logger.warning(msg)
    return resolved


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
        
        # (called some 70 times per step: plain dict and set operations, same rule)
        if k in self.ismutable or not dict.__contains__(self, k):
            dict.__setitem__(self, k, v)
            self.ismutable.add(k)


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

    # Bind imported functions as methods
    solve_EB = solve_EB
    solve_EF = solve_EF
    solve_SS = solve_SS
    
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
        self.tveg: float = 0. # time since last vegetation update (grass)
        self._one_row = None # transect on one row (aeolis.onerow): undecided, False or OneRow

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
        self._one_row = None
        self.p = aeolis.inout.read_configfile(self.configfile)
        aeolis.inout.check_configuration(self.p)
        self.solver = resolve_solver(self.p)

        # the time series are read one column at a time every step; column-major
        # storage makes each column contiguous, so numpy.interp does not copy it
        for key in ('wind_file', 'tide_file', 'wave_file', 'meteo_file'):
            if isinstance(self.p.get(key), np.ndarray) and self.p[key].ndim == 2:
                self.p[key] = np.asfortranarray(self.p[key])

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

        # Determine number of grain size fractions
        self.p['nfractions'] = len(self.p['grain_size'])

        # Initilize number of species and size of subgrid for grass model
        if self.p['method_vegetation'] == 'grass':

            # Determine number of grass species
            if isinstance(self.p['species_names'], np.ndarray):
                self.p['nspecies'] = len(self.p['species_names'])  # number of grass species
            else:
                self.p['nspecies'] = 1 

            # Determine size of vegetation subgrid
            self.p['nx_vsub'] = self.p['nx'] * self.p['veg_res_factor']
            self.p['ny_vsub'] = self.p['ny'] * self.p['veg_res_factor']

        else:
            self.p['nspecies'] = 1  # default to one species

            # Not used if "grass" is not used, but set to avoid errors
            self.p['nx_vsub'] = self.p['nx'] * 1
            self.p['ny_vsub'] = self.p['ny'] * 1

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
         
        # initialize vegetation model
        if self.p['method_vegetation'] == 'duran':
            self.s = aeolis.vegetation.initialize(self.s, self.p)
        elif self.p['method_vegetation'] == 'grass':
            self.s, self.p = aeolis.grass.initialize(self.s, self.p)  
        else:
            logger.log_and_raise('Unknown vegetation method [%s]' % self.p['method_vegetation'], exc=ValueError)                

        # initialize fence model
        self.s = aeolis.fences.initialize(self.s, self.p)

        # initilize bed-interaction parameter zeta
        self.s = aeolis.zeta.initialize(self.s, self.p)

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

        # a transect of identical rows is computed on one row (transect_1d)
        if self._one_row is None:
            self._one_row = aeolis.onerow.start(self)
        
        # store previous state
        self.l = self.s.copy()
        self.l['zb'] = self.s['zb'].copy()
        self.l['dzbavg'] = self.s['dzbavg'].copy()
        
        # interpolate wind time series
        self.s = aeolis.wind.interpolate(self.s, self.p, self.t)  
    
        # Rotate gridparams, such that the grids is alligned horizontally
        self.s = self.grid_rotate(self.p['alpha'])

        #compute sand fence shear
        if self.p['process_fences']:
            self.s = aeolis.fences.update_fences(self.s, self.p)
        
        # topographic steering (including Okin on rotating grid)
        if np.sum(self.s['uw']) != 0:
            self.s = aeolis.wind.shear(self.s, self.p)

        # apply vegetation shear
        if self.p['process_vegetation']:
            if self.p['method_vegetation'] == 'duran':
                self.s = aeolis.vegetation.vegshear(self.s, self.p)
            if self.p['method_vegetation'] == 'grass':
                self.s = aeolis.grass.compute_shear_reduction(self.s, self.p)
                self.s = aeolis.grass.apply_shear_reduction(self.s, self.p)

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

        # compute bed-interaction parameter zeta
        if self.p['process_bedinteraction']:
            self.s = aeolis.zeta.compute_zeta(self.s, self.p)

        # compute saltation velocity and equilibrium transport
        self.s = aeolis.transport.equilibrium(self.s, self.p)

        # compute instantaneous transport
        if self.solver == 'steadystate':
            solve_dict = self.solve_SS()
        elif self.solver == 'euler_backward':
            solve_dict = self.solve_EB(alpha=1.)
        elif self.solver == 'crank_nicolson':
            solve_dict = self.solve_EB(alpha=.5)
        elif self.solver == 'euler_forward':
            solve_dict = self.solve_EF()
        self.s.update(solve_dict)

        # update bed
        self.s = aeolis.bed.update(self.s, self.p)
        if self.p['dt_limiter'] != 'legacy':
            # the aeolian bed change of this step, for the time step limiter
            # (average_change overwrites dzb with the total change)
            self.dzb_aeolian = self.s['dzb'].copy()

        # avalanching
        self.s = aeolis.avalanching.angele_of_repose(self.s, self.p)
        self.s = aeolis.avalanching.avalanche(self.s, self.p)
        
        # reset original bed in marine zone (wet)
        self.s = aeolis.bed.wet_bed_reset(self.s, self.p)

        # calculate average bed level change over time (old method)
        if self.p['process_vegetation'] and self.p['method_vegetation'] == 'duran':
            self.s = aeolis.bed.average_change(self.l, self.s, self.p)

        # compute dune erosion
        if self.p['process_dune_erosion']:
            self.s = aeolis.erosion.run_ph12(self.s, self.p, self.t)
            self.s = aeolis.avalanching.angele_of_repose(self.s, self.p) #Since the aeolian module is 
            # only run for winds above threshold, also run avalanching routine here
            self.s = aeolis.avalanching.avalanche(self.s, self.p)

        # grow vegetation
        if self.p['process_vegetation']:
            if self.p['method_vegetation'] == 'duran':
                self.s = aeolis.vegetation.germinate(self.s, self.p)
                self.s = aeolis.vegetation.grow(self.s, self.p)
            # update grass when dt_veg has passed
            elif self.p['method_vegetation'] == 'grass':
                self.tveg += self.dt * self.p['accfac']
                if self.tveg >= self.p['dt_veg']:
                    self.s = aeolis.grass.update(self.s, self.p)
                    self.tveg -= self.p['dt_veg']

        # increment time
        self.t += self.dt * self.p['accfac']
        if getattr(self, '_t_target', None) is not None:
            # the step was cut to end on an output time or tstop: land on it exactly
            self.t = self._t_target
            self._t_target = None
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
            v = self._one_row.get(self, var) if self._one_row else self.s[var]
            if var in ['Ct', 'Cu']:
                return v / self.p['accfac']
            else:
                return v
        elif var in self.p:
            if self._one_row and var in ('ny', 'ny_vsub'):
                return (self._one_row.n - 1) * (self._one_row.f if var == 'ny_vsub' else 1)
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
            if self._one_row:
                return self._one_row.shape(self, var)
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
            if self._one_row:
                self._one_row.put(self, {var: val}, assigned={var})
            else:
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

        # dt is kept a float: an integer ``dt`` from the configuration file
        # would make the compiled kernels compile a second variant and switch
        # to the float one as soon as the time step limiter first applies
        if dt > 0.:
            self.dt = float(dt)
        elif dt < -1:
            self.dt = float(dt)
            self.t += np.abs(dt)
            return False
        else:
            self.dt = float(self.p['dt'])

        if self.solver == 'euler_forward':
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
           
        if self.p['dt_limiter'] != 'legacy':
            self._set_timestep_aeolian()
            self.p['dt_opt'] = self.dt
            return True

        if self.p['max_bedlevel_change'] != 999. and np.max(self.s['dzb']) != 0. and self.dt_prev != 0.:

            # # limit dt based on pickup rate (if non-limited pickup is larger than available mass in toplayer)
            # # this does not account for individual fractions with small availability!
            # mass_layer = np.sum(self.s['mass'][:,:,0,:], axis=2)
            # ix_zero = mass_layer == 0.
            # r_pickup = np.ones_like(mass_layer)
            # r_pickup[~ix_zero] = self.s['pickup'][:,:,0][~ix_zero] / mass_layer[~ix_zero]
            # dt_pickup = self.dt_prev * 0.7 / np.max(r_pickup) # 0.7 factor to be conservative
            
            # limit dt based pickup rate to max_bedlevel_change
            dt_zb = self.dt_prev * self.p['max_bedlevel_change'] / np.max(np.abs(self.s['dzb']))
            self.dt = np.min([self.dt, dt_zb])#, dt_pickup])
            
        self.p['dt_opt'] = self.dt

        return True

    def _set_timestep_aeolian(self) -> None:
        '''Time step limiter ``dt_limiter = aeolian`` (see constants).

        Limits the aeolian bed level change of the next step to
        ``max_bedlevel_change`` and, with ``dt_layer_limit``, to the sediment
        available in the top layer, using the magnitude of the aeolian bed
        change of the previous step (avalanching and the wet-bed reset do
        not scale with the time step and are left out). The step may grow
        by at most ``dt_growth`` per step (if set) and ends exactly on the next output
        time and on ``tstop``.
        '''
        p = self.p
        limit = p['max_bedlevel_change']
        if p['dt_layer_limit']:
            limit = min(limit, p.get('_layer_availability', limit))
        dzb = getattr(self, 'dzb_aeolian', None)
        if dzb is None:
            dzb = self.s['dzb']
        if limit != 999. and self.dt_prev > 0.:
            dzmax = np.max(np.abs(dzb))
            if dzmax > 0.:
                self.dt = min(self.dt, self.dt_prev * limit / dzmax)
            if p['dt_growth'] > 0.:
                self.dt = min(self.dt, self.dt_prev * p['dt_growth'])

        # end the step on the next output time or on tstop
        acc = p['accfac']
        out = p['output_times']
        t_next = p['tstop']
        if out and out > 0.:
            n_out = np.floor((self.t - p['tstart']) / out + 1e-9) + 1
            t_next = min(t_next, p['tstart'] + n_out * out)
        remaining = (t_next - self.t) / acc
        self._t_target = None
        if 0. < remaining < self.dt * (1. + 1e-9):
            self.dt = remaining
            self._t_target = t_next

    def grid_rotate(self, angle:float) -> ModelState:

        s = self.s
        p = self.p

        if angle == 0.:
            return s    # identity (avoids rounding the coordinates through (x - mean) + mean)

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
        self.tveg = 0.
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


    def init_logger(self):
        """
        Initialize file and console logging for AeoLiS.

        This method is safe to call once before initialize()
        and is reused by run().
        """

        # IMPORTANT: this is the module-level logger used everywhere
        logger = logging.getLogger(__name__)

        # Clear existing handlers (interactive / repeated runs)
        if logger.hasHandlers():
            logger.handlers.clear()

        logger.setLevel(logging.DEBUG)

        # ---- file logger ----
        logfile = '%s.log' % os.path.splitext(self.configfile)[0]
        filehandler = logging.FileHandler(logfile, mode='w')
        filehandler.setLevel(logging.INFO)
        filehandler.setFormatter(
            logging.Formatter(
                '%(asctime)-15s %(name)-8s %(levelname)-8s %(message)s'
            )
        )
        logger.addHandler(filehandler)

        # ---- console logger ----
        streamhandler = logging.StreamHandler()
        streamhandler.setLevel(logging.INFO)
        streamhandler.setFormatter(StreamFormatter())
        logger.addHandler(streamhandler)

        # Ensure submodule loggers propagate here
        logger.propagate = True

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
        logger.info('  Version:  %-45s', __version__)
        logger.info('                                                          ')
        logger.info('**********************************************************\n')


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

        self.init_logger()

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
        while (self.t <= self.p['tstop']) if self.p['tstop_inclusive'] else (self.t < self.p['tstop']):
            if callback is not None:
                if self._one_row:
                    self._one_row.call(self, callback)      # the callback sees all rows
                else:
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

        o = self.o[var]
        if not o:
            # statistics were not requested in the configuration: allocate
            # them now (they stay zero/inf, as before the lazy allocation)
            o = self.o[var] = self._statistic_arrays(var)

        if stat in ['min', 'max']:
            return o[stat]
        elif stat == 'sum':
            return o[stat]*self.dt
        elif stat == 'avg':
            if self.n > 0:
                return o['sum'] / self.n
            else:
                return np.zeros(o['sum'].shape)
        elif stat == 'var':
            if self.n > 1:
                return (o['var'] - o['sum']**2 / self.n) \
                    / (self.n - 1)
            else:
                return np.zeros(o['var'].shape)
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

        self.clear = True

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

        for k, exts in self.p['_output_vars'].items():
            if any(ext in self._STATISTICS for ext in exts):
                self.o[k] = self._statistic_arrays(k)
            else:
                # no statistics requested: nothing to accumulate (allocated
                # on demand by get_statistic if a user asks anyway)
                self.o[k] = {}

        self.n = 0


    _STATISTICS = ('min', 'max', 'sum', 'avg', 'var')

    def _statistic_arrays(self, k:str) -> dict:
        s = self.get_var_shape(k)
        return dict(min=np.zeros(s) + np.inf,
                    max=np.zeros(s) - np.inf,
                    var=np.zeros(s),
                    sum=np.zeros(s))


    def output_update(self) -> None:
        '''Updates output statistics dictionary

        Updates matrices with minimum, maximum, variance and summed
        values for each output variable with current spatial grid
        values and increases time step counter with one.

        '''

        for k, exts in self.p['_output_vars'].items():
            o = self.o[k]
            if not o:
                continue            # no statistics requested for this variable
            v = self.get_var(k, clear=False)        # read only, no copy needed
            if 'min' in exts:
                np.minimum(o['min'], v, out=o['min'])
            if 'max' in exts:
                np.maximum(o['max'], v, out=o['max'])
            if 'sum' in exts or 'avg' in exts or 'var' in exts:
                o['sum'] += v
            if 'var' in exts:
                o['var'] += v**2

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
                        variables[k] = self.get_var(k, clear=False)  # read only
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
                with open(restartfile, 'rb') as fp:
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
        s, l, p = self._one_row.full_state(self) if self._one_row else (self.s, self.l, self.p)
        with open(restartfile, 'wb') as fp:
            pickle.dump({'t':self.t,
                         'p':p,
                         's':s,
                         'l':l,
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
        if n_time > 0:
            logger.info(fmt % ('avg. solves per step',
                               aeolis.inout.print_value(float(n_matrixsolve) / n_time)))
            logger.info(fmt % ('avg. time step',
                               aeolis.inout.print_value(float(self.p['tstop'] - self.p['tstart']) / n_time)))

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
        import matplotlib.pyplot as plt
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
        import matplotlib.pyplot as plt
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
