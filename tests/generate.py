# 1d/2d
# grid size
# wind speed
# wind direction
# fractions
# T
# dt
# bi

import os
import json
import numpy as np


COMBINATIONS = [
    [dict(grain_dist = 1.0,
          grain_size = 0.000225,
          nfractions = 1),
     dict(grain_dist = [0.95, 0.05],
          grain_size = [0.000225, 0.032000],
          nfractions = 2)],
    
    [dict(dt = 60),
     dict(dt = 600)],
    
    [dict(T = 3.0),
     dict(T = 1.0),
     dict(T = 0.1)],
    
    [dict(bi = 0.05),
     dict(bi = 1.0)],
    
    [dict(tstop = 3600,
          output_vars = ['zb', 'Cu.avg', 'Ct.avg', 'uws.avg', 'pickup.sum', 'w'])]
    
]


def generate_bathy(fpath, dx=10., dy=10., **kwargs):

    x = np.arange(0, 100).reshape((1,-1)) * dx
    y = np.arange(0, 10) * dy

    X, Y = np.meshgrid(x, y)
    Z = np.zeros(X.shape)

    TH = np.ones(Z.shape)
    TH[:,:25] = 100.
    TH[:,-25:] = 100.

    p = dict(nx = Z.shape[1]-1,
             ny = Z.shape[0]-1,
             xgrid_file = 'x.txt',
             ygrid_file = 'y.txt',
             bed_file = 'z.txt',
             threshold_mask = 'th.txt')

    np.savetxt(os.path.join(fpath, p['xgrid_file']), X)
    np.savetxt(os.path.join(fpath, p['ygrid_file']), Y)
    np.savetxt(os.path.join(fpath, p['bed_file']), Z)
    np.savetxt(os.path.join(fpath, p['threshold_mask']), TH)

    return p


def generate_wind(fpath, u=10., udir=0., **kwargs):

    data = np.asarray([[0., u, udir],
                       [3600., u, udir]])

    p = dict(wind_file = 'wind.txt')
    
    np.savetxt(os.path.join(fpath, p['wind_file']), data)

    return p

    
def generate_model(fpath, p, **kwargs):

    if not os.path.exists(fpath):
        os.makedirs(fpath)
    
    p.update(generate_bathy(fpath, **kwargs))
    p.update(generate_wind(fpath, **kwargs))

    with open(os.path.join(fpath, 'aeolis.txt'), 'w') as fp:
        for k, v in p.items():
            fp.write('%-10s = %s\n' % (k, format_value(v)))

    with open(os.path.join(fpath, 'aeolis.json'), 'w') as fp:
        data = dict(settings=kwargs,
                    parameters=p)
        json.dump(data, fp, indent=4)


def generate_models():

    n = 0
    for dx in [1, 10]:
        for dy in [10]:
            for u in [10]:
                for udir in [0, 20, 45, 60]:
                    
                    for i, p in enumerate(iterate_models()):
                        fpath = os.path.join(os.path.split(__file__)[0], 'models', 'model%06d' % (1000*n + i))
                        generate_model(fpath, p, dx=dx, dy=dy, u=u, udir=udir),

                    n += 1


def iterate_models(p={}, combinations=COMBINATIONS):

    for c in combinations[0]:
        p.update(c)
        if len(combinations) > 1:
            for p in iterate_models(p=p, combinations=combinations[1:]):
                yield p
        else:
            yield p


def format_value(v):
    if type(v) is list:
        return ' '.join([format_value(vi) for vi in v])
    elif type(v) is int:
        return '%d' % v
    elif type(v) is float:
        return '%0.6f' % v
    else:
        return v

    
if __name__ == '__main__':
    generate_models()
