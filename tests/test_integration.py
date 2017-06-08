from nose.tools import *
from .generate import *

import os
import glob
import netCDF4
import logging

from aeolis.model import AeoLiSRunner


TNY = 1e-3


def test_integration():

    generate_models()

    for fname in glob.glob(os.path.join(os.path.split(__file__)[0], 'models', '*', 'aeolis.txt')):
        AeoLiSRunner(configfile=fname).run()

        fname1 = '%s.nc' % os.path.splitext(fname)[0]
        fname2 = '%s.json' % os.path.splitext(fname)[0]

        results = check_continuity(fname1)

        with open(fname2, 'r') as fp:
            data = json.load(fp)
            data.update(dict(results=results))
        with open(fname2, 'w') as fp:
            json.dump(data, fp, indent=4)
                
        yield assert_true, results['fraction'] < TNY


def check_continuity(fname):

    with netCDF4.Dataset(fname, 'r') as ds:

        # dimensions
        dt = np.diff(ds.variables['time'][:])[0]
        dx = np.diff(ds.variables['x'][0,:])[0]
        dy = np.diff(ds.variables['y'][:,0])[0]
        
        # bed level change
        zb = ds.variables['zb'][:,:,:]
        dz = zb[-1,:,:] - zb[0,:,:]
        V1 = dz * dx * dy
        
        # loss over border
        Ct = np.sum(ds.variables['Ct.avg'][:,:,:,:], axis=-1)
        uws = ds.variables['uws.avg'][:,:,:]
        qs = Ct * uws
        V2 = qs[:,:,-1] * dt * dy / (2650. * .6)
        
        # in saltation
        V3 = Ct[-1,:,:] * dx * dy / (2650. * .6)
        
        # metrics
        E = np.abs(np.minimum(0., V1).sum())
        D = np.abs(np.maximum(0., V1).sum())
        dV = np.abs(V1.sum() + V2.sum())# + V3.sum())
        frac = dV/E

        return dict(erosion=float(E),
                    deposition=float(D),
                    loss=float(V1.sum()),
                    loss_over_border=float(V2.sum()),
                    in_saltation=float(V3.sum()),
                    deficit=float(dV),
                    fraction=float(frac))
