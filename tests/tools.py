import numpy as np
from nose.tools import *


TNY = 1e-10


def assert_equal_array(a, b, **kwargs):
    assert_equal(list(np.asarray(a).flatten()),
                 list(np.asarray(b).flatten()), **kwargs)


def assert_less_array(a, b, **kwargs):
    assert_less(list(np.asarray(a).flatten()),
                list(np.asarray(b).flatten()), **kwargs)


def assert_greater_array(a, b, **kwargs):
    assert_greater(list(np.asarray(a).flatten()),
                   list(np.asarray(b).flatten()), **kwargs)


def assert_almost_equal_array(a, b, **kwargs):
    assert_true(np.all(np.abs(np.asarray(a).flatten() -
                              np.asarray(b).flatten()) < TNY), **kwargs)
