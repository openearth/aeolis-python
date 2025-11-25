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
import re
import time
import shutil
import logging
from webbrowser import UnixBrowser
import numpy as np
from matplotlib import pyplot as plt
from scipy.io import savemat

# package modules
from aeolis.utils import *
from aeolis.constants import *
# from regex import S

# initialize logger
logger = logging.getLogger(__name__)


def read_configfile(configfile, parse_files=True, load_defaults=True):
    '''Read model configuration file

    Updates default model configuration based on a model configuration
    file. The model configuration file should be a text file with one
    parameter on each line. The parameter name and value are seperated
    by an equal sign (=). Any lines that start with a percent sign (%)
    or do not contain an equal sign are omitted.

    Parameters are casted into the best matching variable type. If the
    variable type is ``str`` it is optionally interpreted as a
    filename. If the corresponding file is found it is parsed using
    the ``numpy.loadtxt`` function.

    Parameters
    ----------
    configfile : str
        Model configuration file
    parse_files : bool
        If True, files referred to by string parameters are parsed
    load_defaults : bool
        If True, default settings are loaded and overwritten by the
        settings from the configuration file

    Returns
    -------
    dict
        Dictionary with casted and optionally parsed model
        configuration parameters

    See Also
    --------
    write_configfile
    check_configuration

    '''

    if load_defaults:
        p = DEFAULT_CONFIG.copy()
    else:
        p = {}
    
    if os.path.exists(configfile):
        with open(configfile, 'r') as fp:
            for line in fp:
                if '=' in line and not line.strip().startswith('%'):
                    key, val = line.split('=')[:2]
                    p[key.strip()] = parse_value(val, parse_files=parse_files)
    else:
        logger.log_and_raise('File not found [%s]' % configfile, exc=IOError)
       
    # normalize grain size distribution
    if 'grain_dist' in p:
        #p['grain_dist'] = normalize(p['grain_dist']) # commented to allow distribution for multiple layers.
        p['grain_dist'] = makeiterable(p['grain_dist'])
        p['grain_size'] = makeiterable(p['grain_size'])

    # set default output file, if not given
    if 'output_file' in p and not p['output_file']:
        p['output_file'] = '%s.nc' % os.path.splitext(configfile)[0]

    # set default value for h, if not given
    if 'h' in p and not p['h']:
        p['h'] = p['z']

    if 'process_fences' in p and not p['process_fences']:
        p['process_fences'] = False

    # set default for nsavetimes, if not given
    if 'nsavetimes' in p and not p['nsavetimes']:
        p['nsavetimes'] = int(p['dzb_interval']/p['dt'])   
    
    return p


def write_configfile(configfile, p=None):
    '''Write model configuration file

    Writes model configuration to file. If no model configuration is
    given, the default configuration is written to file. Any
    parameters with a name ending with `_file` and holding a matrix
    are treated as separate files. The matrix is then written to an
    ASCII file using the ``numpy.savetxt`` function and the parameter
    value is replaced by the name of the ASCII file.

    Parameters
    ----------
    configfile : str
        Model configuration file
    p : dict, optional
        Dictionary with model configuration parameters

    Returns
    -------
    dict
        Dictionary with casted and optionally parsed model
        configuration parameters

    See Also
    --------
    read_configfile

    '''

    if p is None:
        p = DEFAULT_CONFIG.copy()

    fmt = '%%%ds = %%s\n' % np.max([len(k) for k in p.keys()])
        
    with open(configfile, 'w') as fp:

        fp.write('%s\n' % ('%' * 70))
        fp.write('%%%% %-64s %%%%\n' % 'AeoLiS model configuration')
        fp.write('%%%% Date: %-58s %%%%\n' % time.strftime('%Y-%m-%d %H:%M:%S'))
        fp.write('%s\n' % ('%' * 70))
        fp.write('\n')
        
        for k, v in sorted(p.items()):
            if k.endswith('_file') and isiterable(v):
                fname = '%s.txt' % k.replace('_file', '')
                backup(fname)
                np.savetxt(fname, v)
                fp.write(fmt % (k, fname))
            else:
                fp.write(fmt % (k, print_value(v, fill='')))


def check_configuration(p):
    '''Check model configuration validity

    Checks if required parameters are set and if the references files
    for bathymetry, wind, tide and meteorological input are
    valid. Throws an error if one or more requirements are not met.

    Parameters
    ----------
    p : dict
        Model configuration dictionary with parsed files

    See Also
    --------
    read_configfile

    '''

    # check validity of configuration
    if not isarray(p['xgrid_file']) or \
       not isarray(p['bed_file']) or (not isarray(p['ygrid_file']) and p['ny'] > 0):
        logger.log_and_raise('Incomplete bathymetry definition', exc=ValueError)

    if isarray(p['wind_file']):
        if p['wind_file'].ndim != 2 or p['wind_file'].shape[1] < 3:
            logger.log_and_raise('Invalid wind definition file', exc=ValueError)

    if isarray(p['tide_file']):
        if p['tide_file'].ndim != 2 or p['tide_file'].shape[1] < 2:
            logger.log_and_raise('Invalid tide definition file', exc=ValueError)
            
    if isarray(p['meteo_file']):
        if p['meteo_file'].ndim != 2 or p['meteo_file'].shape[1] < 6:
            logger.log_and_raise('Invalid meteo definition file', exc=ValueError)
            
    if p['th_humidity']:
        logger.warning('Wind velocity threshold based on air humidity following Arens (1996) '
                       'is implemented for testing only. Use with care.')

    if p['th_salt']:
        logger.warning('Wind velocity threshold based on salt content following Nickling and '
                       'Ecclestone (1981) is implemented for testing only. Use with care.')
        
    if p['method_roughness'] == 'constant':        
        logger.warning('Warning: the used roughness method (constant) defines the z0 as '
                       'k (z0 = k), this was implemented to ensure backward compatibility '
                       'and does not follow the definition of Nikuradse (z0 = k / 30).')
    
    # check if steadystate solver is used with multiple sediment fractions
    if p['solver'].lower() in ['steadystate', 'steadystatepieter']:
        if len(p['grain_size']) > 1:
            logger.log_and_raise('The steadystate solver is not compatible with multiple sediment fractions. '
                                 'Please use a single sediment fraction or switch to a different solver (e.g., trunk or pieter).', 
                                 exc=ValueError)

        
def parse_value(val, parse_files=True, force_list=False):
    '''Casts a string to the most appropriate variable type

    Parameters
    ----------
    val : str
        String representation of value
    parse_files : bool
        If True, files referred to by string parameters are parsed by
        ``numpy.loadtxt``
    force_list:
        If True, interpret the value as a list, even if it consists of
        a single value

    Returns
    -------
    misc
        Casted value

    Examples
    --------
    >>> type(parse_value('T'))
        bool
    >>> type(parse_value('F'))
        bool
    >>> type(parse_value('123'))
        int
    >>> type(parse_value('123.2'))
        float
    >>> type(parse_value('euler_forward'))
        str
    >>> type(parse_value(''))
        NoneType
    >>> type(parse_value('zb zs Ct'))
        numpy.ndarray
    >>> type(parse_value('zb', force_list=True))
        numpy.ndarray
    >>> type(parse_value('0.1 0.2 0.3')[0])
        float
    >>> type(parse_value('wind.txt'), parse_files=True)
        numpy.ndarray
    >>> type(parse_value('wind.txt'), parse_files=False)
        str

    '''

    # New initial steps to filter out commented lines (with %)
    if '%' in val:
        val = val.split('%')[0]
        
    val = val.strip()
    
    if ' ' in val or force_list:
        return np.asarray([parse_value(x) for x in val.split(' ')])
    elif re.match('^[TF]$', val):
        return val == 'T'
    elif re.match(r'^-?\d+$', val):
        return int(val)
    elif re.match(r'^-?[\d.]+$', val):
        return float(val)
    elif re.match('None', val):
        return None
    elif os.path.isfile(val) and parse_files:
        for dtype in [float, complex]:
            try:
                val = np.loadtxt(val, dtype=dtype)
                break
            except:
                pass
        return val
    elif val == '':
        return None
    else:
        return val



def backup(fname):
    '''Creates a backup file of the provided file, if it exists'''
    
    if os.path.exists(fname):
        backupfile = get_backupfilename(fname)
        shutil.copyfile(fname, backupfile)

        
def get_backupfilename(fname):
    '''Returns a non-existing backup filename'''
    
    for n in range(1, 1000):
        backupfile = '%s~%03d' % (fname, n)
        if not os.path.exists(backupfile):
            break

    if os.path.exists(backupfile):
        logger.log_and_raise('Too many backup files in use! Please clean up...', exc=ValueError)
    
    return backupfile
            

        
def visualize_grid(s, p):
    '''Create figures and tables for the user to check whether the grid-input is correctly interpreted'''
    
    # Read the x,y,z dimensions
    x = s['x']
    y = s['y']
    zb = s['zb']

    # Read the angle of the rotated grid and avoid negative values
    alpha = p['alpha']
    if alpha < 0.:
        alpha += 360.

    # Determine the maximum dimensions in x- and y-direction
    xlen = np.max(x)-np.min(x)
    ylen = np.max(y)-np.min(y)
    xylen = np.maximum(xlen, ylen)

    # Compute the coordinates for the arc of the angle alpha
    arc_angles = np.linspace(270., 270. + alpha, 1+int(alpha))
    radius = np.minimum(xlen, ylen) * 0.05
    arc_x = x[0,0] + radius * np.cos(np.deg2rad(arc_angles))
    arc_y = y[0,0] + radius * np.sin(np.deg2rad(arc_angles))
    
    # Compute coordinates of labels to indicate boundary types
    x_offshore = np.mean([x[0,0], x[-1,0]])
    y_offshore = np.mean([y[0,0], y[-1,0]])
    x_onshore = np.mean([x[0,-1], x[-1,-1]])
    y_onshore = np.mean([y[0,-1], y[-1,-1]])
    x_lateralA = np.mean([x[0,0], x[0,-1]])
    y_lateralA = np.mean([y[0,0], y[0,-1]])
    x_lateralB = np.mean([x[-1,0], x[-1,-1]])
    y_lateralB = np.mean([y[-1,0], y[-1,-1]])

    # Create plots
    fig, ax = plt.subplots()

    if p['ny'] > 0:
        pc = ax.pcolormesh(x,y,zb)
    else:
        pc = ax.scatter(x,y,c=zb)

    # Plot all the texts
    plottxts = []
    plottxts.append(ax.text(x_offshore, y_offshore, 'Offshore: ' + p['boundary_offshore'], rotation=alpha + 90, ha = 'center', va='center'))
    plottxts.append(ax.text(x_onshore, y_onshore, 'Onshore: ' + p['boundary_onshore'], rotation=alpha + 270, ha = 'center', va='center'))
    plottxts.append(ax.text(x_lateralA, y_lateralA, 'Lateral: ' + p['boundary_lateral'], rotation=alpha + 0, ha = 'center', va='center'))
    plottxts.append(ax.text(x_lateralB, y_lateralB, 'Lateral: ' + p['boundary_lateral'], rotation=alpha + 180, ha = 'center', va='center'))
    plottxts.append(ax.text(x[0,0], y[0,0], '(0,0)', ha = 'right', va='top'))
    plottxts.append(ax.text(x[0,-1], y[0,-1], '(0,' + str(len(x[0,:])-1) + ')', ha = 'right', va='top'))
    plottxts.append(ax.text(x[-1,0], y[-1,0], '(' + str(len(x[:,0])-1) + ',0)', ha = 'right', va='top'))
    plottxts.append(ax.text(x[-1,-1], y[-1,-1], '(' + str(len(x[:,0])-1) + ',' + str(len(x[0,:])-1) + ')', ha = 'right', va='top'))
    plottxts.append(ax.text(x[0,0], y[0,0]-0.1*ylen, r'$\alpha$ :' + str(int(alpha)) + r'$^\circ$', ha='center', va='center'))

    # Set boxes around the texts
    for txt in plottxts:
        txt.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='black'))

    # Plot dots to indicate the corner-points
    ax.plot(x[0,0], y[0,0], 'ro')
    ax.plot(x[0,-1], y[0,-1], 'ro')
    ax.plot(x[-1,0], y[-1,0], 'ro')
    ax.plot(x[-1,-1], y[-1,-1], 'ro')
    
    # Plot the arc to indicate angle
    ax.plot(arc_x, arc_y, color = 'red')
    ax.plot([x[0,0], x[0,0]], [y[0,0],y[0,0]-0.08*ylen], '--', color = 'red', linewidth=3)
    ax.plot([x[0,0], arc_x[-1]], [y[0,0], arc_y[-1]], color = 'red', linewidth=3)
    
    # Figure lay-out settings
    fig.colorbar(pc, ax=ax)
    ax.axis('equal')
    ax.set_xlim([np.min(x) - 0.15*xylen, np.max(x) + 0.15*xylen])
    ax.set_ylim([np.min(y) - 0.15*xylen, np.max(y) + 0.15*xylen])
    height = 8.26772 # A4 width
    width = 11.6929 # A4 height
    fig.set_size_inches(width, height)
    plt.tight_layout()

    # Saving and plotting figure
    fig.savefig('figure_grid_initialization.png', dpi=200)
    plt.close()

    return 

def visualize_timeseries(p, t):
    '''Create figures and tables for the user to check whether the timeseries-input is correctly interpreted'''
   
    # Initiate plots
    fig, axs = plt.subplots(5, 1)

    # Start and stop times
    tstart = p['tstart']
    tstop = p['tstop']
        
    # Read the user input (wind)
    uw_t = p['wind_file'][:,0]
    uw_s = p['wind_file'][:,1]
    uw_d = p['wind_file'][:,2]
    axs[0].plot(uw_t, uw_s, 'k')
    axs[1].plot(uw_t, uw_d, 'k')
    axs[0].set_title('Wind velocity at height z, uw (m/s)')
    axs[1].set_title('Wind direction, udir (deg)')

    # Read the user input (waves)
    
    if np.shape(p['wave_file'])[1]== 3:
        w_t = p['wave_file'][:,0]
        w_Hs = p['wave_file'][:,1]
        w_Tp = p['wave_file'][:,2]
        axs[2].plot(w_t, w_Hs, 'k')
        axs[3].plot(w_t, w_Tp, 'k')
        axs[2].set_title('Wave height, Hs (m)')
        axs[3].set_title('Wave period, Tp (sec)')

    # Read the user input (tide)
    if np.shape(p['tide_file'])[1]==2:
        T_t = p['tide_file'][:,0]
        T_zs = p['tide_file'][:,1]
        axs[4].plot(T_t, T_zs, 'k')
        axs[4].set_title('Water level, zs (m)')




    # Plotting
   

   

    # Assiging titles
    
    

    for ax in axs:
        ax.set_xlim([tstart, tstop])
        ax.set_xlabel('Time since refdate (s) from tstart (=' + str(tstart) + ') to tstop (=' + str(tstop) + ')')

    width = 8.26772 # A4 width
    height = 11.6929 # A4 height
    fig.set_size_inches(width, height)
    plt.tight_layout()

    # Saving and plotting figure
    fig.savefig('figure_timeseries_initialization.png', dpi=200)
    plt.close()


def visualize_spatial(s, p):
    '''Create figures and tables for the user to check whether the input is correctly interpreted'''
    
    # Read the x,y dimensions
    x = s['x']
    y = s['y']
    
    # Reading masks and if constant, fill 2D-array
    uth_mask_multi = np.ones(np.shape(x)) * np.real(s['threshold_mask'])
    tide_mask_multi = np.ones(np.shape(x)) * np.real(s['tide_mask'])
    wave_mask_multi = np.ones(np.shape(x)) * np.real(s['wave_mask'])

    uth_mask_add = np.ones(np.shape(x)) * np.imag(s['threshold_mask'])
    tide_mask_add = np.ones(np.shape(x)) * np.imag(s['tide_mask'])
    wave_mask_add = np.ones(np.shape(x)) * np.imag(s['wave_mask'])

    # Determine the maximum dimensions in x- and y-direction
    xlen = np.max(x)-np.min(x)
    ylen = np.max(y)-np.min(y)

    # Creating values
    fig, axs = plt.subplots(5, 3)
    pcs = [[None for _ in range(3)] for _ in range(5)]

    # In the plotting below, prevent the UserWarning: The input coordinates to pcolormesh are interpreted as cell centers, but are not monotonically increasing or decreasing (...)
    import warnings 
    warnings.filterwarnings("ignore", category=UserWarning)

    # Plotting colormeshes
    if p['ny'] > 0:
        pcs[0][0] = axs[0,0].pcolormesh(x, y, s['zb'], cmap='viridis')
        pcs[0][1] = axs[0,1].pcolormesh(x, y, s['zne'], cmap='viridis')
        pcs[0][2] = axs[0,2].pcolormesh(x, y, s['rhoveg'], cmap='Greens', clim= [0, 1])
        pcs[1][0] = axs[1,0].pcolormesh(x, y, s['uw'], cmap='plasma')
        pcs[1][1] = axs[1,1].pcolormesh(x, y, s['ustar'], cmap='plasma')
        pcs[1][2] = axs[1,2].pcolormesh(x, y, s['u'][:, :, 0], cmap='plasma')
        pcs[2][0] = axs[2,0].pcolormesh(x, y, s['moist'], cmap='Blues', clim= [0, 0.4])
        pcs[2][1] = axs[2,1].pcolormesh(x, y, s['gw'], cmap='viridis')
        pcs[2][2] = axs[2,2].pcolormesh(x, y, s['uth'][:,:,0], cmap='plasma')
        pcs[3][0] = axs[3,0].pcolormesh(x, y, uth_mask_multi, cmap='binary', clim= [0, 1])
        pcs[3][1] = axs[3,1].pcolormesh(x, y, tide_mask_multi, cmap='binary', clim= [0, 1])
        pcs[3][2] = axs[3,2].pcolormesh(x, y, wave_mask_multi, cmap='binary', clim= [0, 1])
        pcs[4][0] = axs[4,0].pcolormesh(x, y, uth_mask_add, cmap='binary', clim= [0, 1])
        pcs[4][1] = axs[4,1].pcolormesh(x, y, tide_mask_add, cmap='binary', clim= [0, 1])
        pcs[4][2] = axs[4,2].pcolormesh(x, y, wave_mask_add, cmap='binary', clim= [0, 1])
    else:
        pcs[0][0] = axs[0,0].scatter(x, y, c=s['zb'], cmap='viridis')
        pcs[0][1] = axs[0,1].scatter(x, y, c=s['zne'], cmap='viridis')
        pcs[0][2] = axs[0,2].scatter(x, y, c=s['rhoveg'], cmap='Greens', clim= [0, 1])
        pcs[1][0] = axs[1,0].scatter(x, y, c=s['uw'], cmap='plasma')
        pcs[1][1] = axs[1,1].scatter(x, y, c=s['ustar'], cmap='plasma')
        pcs[1][2] = axs[1,2].scatter(x, y, c=s['tau'], cmap='plasma')
        pcs[2][0] = axs[2,0].scatter(x, y, c=s['moist'], cmap='Blues', clim= [0, 0.4])
        pcs[2][1] = axs[2,1].scatter(x, y, c=s['gw'], cmap='viridis')
        pcs[2][2] = axs[2,2].scatter(x, y, c=s['uth'][:,:,0], cmap='plasma')
        pcs[3][0] = axs[3,0].scatter(x, y, c=uth_mask_multi, cmap='binary', clim= [0, 1])
        pcs[3][1] = axs[3,1].scatter(x, y, c=tide_mask_multi, cmap='binary', clim= [0, 1])
        pcs[3][2] = axs[3,2].scatter(x, y, c=wave_mask_multi, cmap='binary', clim= [0, 1])
        pcs[4][0] = axs[4,0].scatter(x, y, c=uth_mask_add, cmap='binary', clim= [0, 1])
        pcs[4][1] = axs[4,1].scatter(x, y, c=tide_mask_add, cmap='binary', clim= [0, 1])
        pcs[4][2] = axs[4,2].scatter(x, y, c=wave_mask_add, cmap='binary', clim= [0, 1])

    # Re-allow the UserWarning
    warnings.filterwarnings("default", category=UserWarning)

    # Quiver for vectors
    skip = 10
    axs[1,0].quiver(x[::skip, ::skip], y[::skip, ::skip], s['uws'][::skip, ::skip], s['uwn'][::skip, ::skip])
    axs[1,1].quiver(x[::skip, ::skip], y[::skip, ::skip], s['ustars'][::skip, ::skip], s['ustarn'][::skip, ::skip])
    axs[1,2].quiver(x[::skip, ::skip], y[::skip, ::skip], s['us'][::skip, ::skip, 0], s['un'][::skip, ::skip, 0])

    # Adding titles to the plots
    axs[0,0].set_title('Bed level, zb (m)')
    axs[0,1].set_title('Non-erodible layer, zne (m)')
    axs[0,2].set_title('Vegetation density, rhoveg (-)')
    axs[1,0].set_title('Wind velocity, uw (m/s)')
    axs[1,1].set_title('Shear velocity, ustar (m/s)')
    axs[1,2].set_title('Grain velocity, u (m/s)')
    axs[2,0].set_title('Soil moisture content, (-)')
    axs[2,1].set_title('Ground water level, gw (m)')
    axs[2,2].set_title('Velocity threshold (0th fraction), uth (m/s)')
    axs[3,0].set_title('Threshold multiplication mask (-)')
    axs[3,1].set_title('Tide multiplication mask (-)')
    axs[3,2].set_title('Wave multiplication mask (-)')
    axs[4,0].set_title('Threshold addition mask (-)')
    axs[4,1].set_title('Tide addition mask (-)')
    axs[4,2].set_title('Wave addition mask (-)')

    # Formatting the plot
    for irow, ax_rows in enumerate(axs):
        for icol, ax in enumerate(ax_rows):
        # Figure lay-out settings
            fig.colorbar(pcs[irow][icol], ax=ax)
            ax.axis('equal')
            xylen = np.maximum(xlen, ylen)
            ax.set_xlim([np.min(x) - 0.15*xylen, np.max(x) + 0.15*xylen])
            ax.set_ylim([np.min(y) - 0.15*xylen, np.max(y) + 0.15*xylen])
            ax.axes.xaxis.set_visible(False)
            ax.axes.yaxis.set_visible(False)
            width = 8.26772*2 # A4 width
            height = 11.6929*2 # A4 height
            fig.set_size_inches(width, height)
            plt.tight_layout()
    
    # Saving and plotting figure
    fig.savefig('figure_params_initialization.png', dpi=300)
    plt.close()

    return 

def output_sedtrails(s, p):
    '''Create additional output for SedTRAILS and save as mat-files.
    Chosen for seperate files, such that only relevant (Ct > 0) cells 
    are exported for memory and speed efficiency''' 

    nf = p['nfraction_sedtrails']

    # Speed and concetration: Only for the first fraction now
    x = s['x'].flatten()
    y = s['y'].flatten()
    us = s['usST'][:,:,nf].flatten()
    un = s['unST'][:,:,nf].flatten()
    pickup = s['pickup'][:,:,nf].flatten()
    dzb = s['dzb'].flatten() # Store the bed level change (AEOLIAN ONLY) for every timestep

    os.makedirs('sedtrails_output', exist_ok=True) 
    
    time = p['_time']
    if time == 0: # Save the x and y coordinates only once to save memory
        mdic = {'x': x, 'y': y, 'us': us, 'un': un, 'dzb': dzb, 'pickup': pickup}
        savemat(os.path.join('sedtrails_output', str(int(time)).zfill(12) + '.mat'), mdic)
    else:
        mdic = {'us': us, 'un': un, 'dzb': dzb, 'pickup': pickup}
        savemat(os.path.join('sedtrails_output', str(int(time)).zfill(12) + '.mat'), mdic)
