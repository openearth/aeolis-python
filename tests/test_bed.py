'''This module tests the function in bed.py. Most importantly
continuity over the fractions needs to be ensured even in extreme
cases where some fractions erode and others accrete. Erosion is
limited in the model by the sediment availability in the top bed
composition layer, but deposition is not. Therefore deposition is
unbounded and the most complex of the two phenomena.

'''

from nose.tools import *
from tools import *

import numpy as np
import copy

import aeolis


# dimensions
NX = 0
NY = 0
NL = 3
NF = 4

# parameters
P = {
    'bedupdate':True,
    'nx':NX,
    'ny':NY,
    'nlayers':NL,
    'nfractions':NF,
    'rhop':2650.,
    'porosity':.4,
    'grain_dist':np.ones((NF,)),
}

# variables
S = {
    'zb':np.zeros((NY+1, NX+1)),
    'pickup':np.zeros((NY+1, NX+1, NF)),
    'mass':np.ones((NY+1, NX+1, NL, NF)),
}


def assert_continuity(s1, s2=S):
    '''Convenience function to test whether sediment mass in bed layers is positive and constant to ensure continuity

    Parameters
    ----------
    s : dict
        Result structure from bed.update()

    '''

    print s1['mass']
    
    assert_true(np.all(s1['mass'] >= 0.),
                msg='Layer mass is negative')

    assert_almost_equal_array(s1['mass'].sum(axis=3),
                              s2['mass'].sum(axis=3),
                              msg='Layer mass not constant')

    
def test_trivial():
    '''Test if zero pickup leads to no changes in bed composition and level'''

    s = copy.deepcopy(S)
    s = aeolis.bed.update(s, P)
    assert_continuity(s)

    assert_equal_array(s['mass'],
                       S['mass'],
                       msg='Bed composition changed')
    
    assert_equal_array(s['zb'],
                       S['zb'],
                       msg='Bed level changed')


def test_erosion_uniform():
    '''Test if uniform erosion on a uniform bed leads to no changes in bed composition and a decrease in bed level'''

    s = copy.deepcopy(S)
    s['pickup'][0,0,:] = .25 / 4.
    s = aeolis.bed.update(s, P)
    assert_continuity(s)

    assert_almost_equal_array(s['mass'],
                              S['mass'],
                              msg='Bed composition changed')

    assert_less_array(s['zb'],
                      S['zb'],
                      msg='Bed level did not decrease')


def test_erosion_singlefraction():
    '''Test if erosion of a single fraction from a uniform bed leaves the other fractions unaffected'''

    s = copy.deepcopy(S)
    s['pickup'][0,0,0] = .25
    s = aeolis.bed.update(s, P)
    assert_continuity(s)

    assert_almost_equal_array(np.abs(np.diff(s['mass'][:,:,:,1:], axis=3)),
                              np.zeros((NY+1, NX+1, NL, NF-2)),
                              msg='Non-erodible fractions changed')
    

def test_erosion_mixed():
    '''Test if continuity is ensured in a net erosion cell with a single accretive fraction'''

    s = copy.deepcopy(S)
    s['pickup'][:,:,:] = [.95, .95, -.95, 0.]
    s = aeolis.bed.update(s, P)
    assert_continuity(s)
    
    
def test_erosion_progressive():
    '''Test if progressive erosion only affects top layer and continiously decrease the bed level'''

    s = copy.deepcopy(S)
    s['pickup'][0,0,:] = .25 * np.asarray([.6, .3, .1, 0.]) # sum: .25

    for i in range(NL):
        s = aeolis.bed.update(s, P)
        assert_continuity(s)

        assert_almost_equal_array(s['mass'][:,:,1:,:],
                                  S['mass'][:,:,1:,:],
                                  msg='Other layers than top layer affected')

        assert_less_array(s['zb'],
                          S['zb'],
                          msg='Bed level did not decrease')


def test_deposition_uniform():
    '''Test if uniform deposition on a uniform bed leads to no changes in bed composition and an increase in bed level'''

    s = copy.deepcopy(S)
    s['pickup'][0,0,:] = -.25 / 4.
    s = aeolis.bed.update(s, P)
    assert_continuity(s)

    assert_almost_equal_array(s['mass'],
                              S['mass'],
                              msg='Bed composition changed')

    assert_greater_array(s['zb'],
                         S['zb'],
                         msg='Bed level did not increase')


def test_deposition_huge():
    '''Test if continuity is ensured if an amount of sediment larger than the total contents of a bed composition layer is deposited'''

    s1 = copy.deepcopy(S)
    s1['mass'][:,:,:,0] -= .5
    s1['mass'][:,:,:,-1] += .5
    s1['pickup'][:,:,:] = -10. * s1['mass'][:,:,0,:].sum(axis=2) / NF
    s2 = aeolis.bed.update(copy.deepcopy(s1), P)
    assert_continuity(s2, s1)


def test_deposition_mixed():
    '''Test if continuity is ensured in a net deposition cell with a single erosive fraction'''

    s = copy.deepcopy(S)
    s['pickup'][:,:,:] = [-.95, -.95, .95, 0.]
    s = aeolis.bed.update(s, P)
    assert_continuity(s)
    
    
def test_deposition_progressive():
    '''Test if progressive deposition only affects an increasing number of top layer and continiously increase the bed level'''

    s = copy.deepcopy(S)
    s['pickup'][0,0,:] = -.25 * np.asarray([.6, .3, .1, 0.]) # sum: -.25

    for i in range(NL):
        s = aeolis.bed.update(s, P)
        assert_continuity(s)
        
        assert_almost_equal_array(s['mass'][:,:,(i+1):,:],
                                  S['mass'][:,:,(i+1):,:],
                                  msg='Other layers than top #%d layers affected' % (i+1))

        assert_greater_array(s['zb'],
                             S['zb'],
                             msg='Bed level did not increase')
