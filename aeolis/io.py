import os
import re
import log
import numpy as np

# package modules
from utils import *

#: AeoLiS model default configuration
DEFAULT_CONFIG = dict(th_grainsize=True,
                      th_bedslope=False,
                      th_moisture=True,
                      th_humidity=False,
                      th_roughness=False,
                      mixtoplayer=True,
                      xgrid_file=None,
                      ygrid_file=None,
                      bed_file=None,
                      wind_file=None,
                      tide_file=None,
                      meteo_file=None,
                      bedcomp_file=None,
                      tstart=0.,
                      tstop=3600.,
                      outputtimes=60.,
                      outputfile='aeolis.nc',
                      outputvars=['zb', 'zs', 'Ct', 'Cu', 'uw', 'uth', 'mass', 'pickup'],
                      outputtypes=[],
                      grain_size=[225e-6],
                      grain_dist=[1.],
                      nfractions=1,
                      nlayers=3,
                      layer_thickness=.01,
                      g=9.81,
                      rhoa=1.25,
                      rhop=2650.,
                      rhow=1025.,
                      porosity=.4,
                      A=100.,
                      z0=1.,
                      k=0.01,
                      Cb=1.5,
                      bi=1.,
                      T=1.,
                      F=1e-4,
                      eps=1e-3,
                      facDOD=.1,
                      scheme='euler_backward',
                      method_moist='belly_johnson',
                      max_error=1e-6,
                      max_iter=1000,
                      accfac=1.)

#: Required model configuration paramaters
REQUIRED_CONFIG = ['nx', 'ny']


def read_configfile(configfile, parse_files=True):
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

    Returns
    -------
    dict
        Dictionary with casted and optionally parsed model
        configuration parameters

    See Also
    --------
    aeolis.io.check_configuration
    '''

    p = DEFAULT_CONFIG.copy()
    
    if os.path.exists(configfile):
        with open(configfile, 'r') as fp:
            for line in fp:
                if '=' in line and not line.strip().startswith('%'):
                    key, val = line.split('=')[:2]
                    p[key.strip()] = parse_value(val, parse_files=parse_files)
    else:
        err = 'File not found [%s]' % configfile
        log.logging.error(err)
        raise IOError(err)

    return p


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
    aeolis.io.read_configfile

    '''
    
    logger = log.Logger()
    
    # check for missing parameters
    missing = [k for k in REQUIRED_CONFIG if not p.has_key(k)]
    if len(missing) > 0:
        log.error('Missing required parameters [%s]' % ', '.join(missing))

    # check validity of configuration
    if not isarray(p['xgrid_file']) or \
       not isarray(p['bed_file']) or (isarray(p['ygrid_file']) and p['ny'] > 0):
        logger.error('Incomplete bathymerty definition')

    if not isarray(p['wind_file']) or \
       p['wind_file'].ndim != 2 or p['wind_file'].shape[1] < 3:
        logger.error('Invalid wind definition file')

    if p['wind_file'][-1,0] < p['tstop']:
        logger.error('Wind definition file too short')

    if isarray(p['tide_file']):
        if p['tide_file'].ndim != 2 or p['tide_file'].shape[1] < 2:
            logger.error('Invalid tide definition file')
            
        if p['tide_file'][-1,0] < p['tstop']:
            logger.error('Tide definition file too short')

    if isarray(p['meteo_file']):
        if p['meteo_file'].ndim != 2 or p['meteo_file'].shape[1] < 7:
            logger.error('Invalid meteo definition file')
            
        if p['meteo_file'][-1,0] < p['tstop']:
            logger.error('Meteo definition file too short')

    logger.flush(exception=ValueError('Invalid model configuration'))


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
        a single valye

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

    val = val.strip()
    
    if ' ' in val or force_list:
        return np.asarray([parse_value(x) for x in val.split(' ')])
    elif re.match('^[TF]$', val):
        return val == 'T'
    elif re.match('^-?\d+$', val):
        return int(val)
    elif re.match('^-?[\d\.]+$', val):
        return float(val)
    elif os.path.isfile(val) and parse_files:
        try:
            return np.loadtxt(val)
        except:
            return val
    elif val == '':
        return None
    else:
        return val
