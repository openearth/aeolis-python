import os
import numpy as np


#root = '/Users/hoonhout/GitHub/aeolis-python/example/'
root = '.'

x0 = np.loadtxt(os.path.join(root, 'x.txt')).reshape((1,-1))
z0 = np.loadtxt(os.path.join(root, 'z.txt')).reshape((1,-1))


def move_bar(t, zb, zs, a=1., sigma=50.):
    global x0, z0

    mu = np.minimum(450., 200. + 1./3600. * t)
    dz = a * np.exp(-(x0-mu)**2/sigma**2)
    ix = np.where((x0 > mu-2*sigma) & (x0 < mu+2*sigma))[1]
    
    zb_new = (z0 + dz).reshape(zb.shape)
    ix = ix[np.argmax(zb_new[0,ix])]
    zb_new[0,ix:] = np.maximum(zb_new[0,ix], zb_new[0,ix:])
    
    if t > 0:
        zs_new = zs + zb_new - zb
    else:
        zs_new = zs
        
    return zb_new, zs_new
    
    
def add_bar(model):
    zb = model.get_var('zb')
    zs = model.get_var('zs')
    
    zb_new, zs_new = move_bar(model.get_current_time(), zb, zs)
    
    model.set_var('zb', zb_new)
    model.set_var('zs', zs_new)
