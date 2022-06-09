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
import numpy as np
from matplotlib import pyplot as plt

# package modules
from aeolis.utils import *
from aeolis.constants import *
from regex import S

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

    fmt = '%%%ds = %%s\n' % np.max([len(k) for k in p.iterkeys()])
        
    with open(configfile, 'w') as fp:

        fp.write('%s\n' % ('%' * 70))
        fp.write('%%%% %-64s %%%%\n' % 'AeoLiS model configuration')
        fp.write('%%%% Date: %-58s %%%%\n' % time.strftime('%Y-%m-%d %H:%M:%S'))
        fp.write('%s\n' % ('%' * 70))
        fp.write('\n')
        
        for k, v in sorted(p.iteritems()):
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
    elif re.match('^-?\d+$', val):
        return int(val)
    elif re.match('^-?[\d\.]+$', val):
        return float(val)
    elif re.match('None', val):
        return None
    elif os.path.isfile(val) and parse_files:
        for dtype in [np.float, np.complex]:
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
            

        
def interpretation_visualization(s, p):
    '''Create figures and tables for the user to check whether the input is correctly interpreted'''
    
    # Grid definition + boundaries + wind direction at t=0 + index / sizes + alpha
    x = s['x']
    y = s['y']
    zb = s['zb']

    uws = s['uws']
    print(uws[0])

    alpha = p['alpha']
    if alpha < 0.:
        alpha += 360.

    xlen = np.max(x)-np.min(x)
    ylen = np.max(y)-np.min(y)

    arc_angles = np.linspace(270., 270. + alpha, int(alpha))
    radius = np.minimum(xlen, ylen) * 0.05
    arc_x = x[0,0] + radius * np.cos(np.deg2rad(arc_angles))
    arc_y = y[0,0] + radius * np.sin(np.deg2rad(arc_angles))
    
    x_offshore = np.mean([x[0,0], x[-1,0]])
    y_offshore = np.mean([y[0,0], y[-1,0]])
    x_onshore = np.mean([x[0,-1], x[-1,-1]])
    y_onshore = np.mean([y[0,-1], y[-1,-1]])
    x_lateralA = np.mean([x[0,0], x[0,-1]])
    y_lateralA = np.mean([y[0,0], y[0,-1]])
    x_lateralB = np.mean([x[-1,0], x[-1,-1]])
    y_lateralB = np.mean([y[-1,0], y[-1,-1]])

    fig, ax = plt.subplots()
    pc = ax.pcolormesh(x,y,zb)

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

    for txt in plottxts:
        txt.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='black'))

    ax.plot(x[0,0], y[0,0], 'ro')
    ax.plot(x[0,-1], y[0,-1], 'ro')
    ax.plot(x[-1,0], y[-1,0], 'ro')
    ax.plot(x[-1,-1], y[-1,-1], 'ro')
    
    ax.plot(arc_x, arc_y, color = 'red')
    ax.plot([x[0,0], x[0,0]], [y[0,0],y[0,0]-0.08*ylen], '--', color = 'red', linewidth=3)
    ax.plot([x[0,0], arc_x[-1]], [y[0,0], arc_y[-1]], color = 'red', linewidth=3)
    
    fig.colorbar(pc, ax=ax)
    ax.axis('equal')
    ax.set_xlim([np.min(x) - 0.15*xlen, np.max(x) + 0.15*xlen])
    ax.set_ylim([np.min(y) - 0.15*ylen, np.max(y) + 0.15*ylen])

    # Draw arc for angle


    plt.show()
    print('test')


    # Spatial values and masks
    # Bed level
    # Ne level
    # Vegetation
    # Threshold mask
    # Tide mask
    # Wave mask

    # Timeseries
    # Wind
    # Tide
    # Waves

    # Grainsizes

    return 
