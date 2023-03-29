"""
This file is part of AeoLiS test suit. It uses the Pytest framework.

To run all tests, use:
    pytest

To run only the test in this file, use:
    pytest aeolis/tests/test_utils.py

To run all test in a group (class), use:
    pytest aeolis/tests/test_utils.py::<TestClass>

To run an specific test, use:
    - If in a group:
        pytest aeolis/tests/test_utils.py::<TestClass>::<test-method>
    
    - If not in any group:
        pytest aeolis/tests/test_utils.py::<test-function>
"""

from aeolis import utils
import numpy as np


# Use test classes to group test-cases that belong to the 
# same scenario; unless a single test-case is required

class TestIsIterable:
    """
    Test if iterable variables can be correctly distinguished from non-iterable variables
    """
    
    def test_isiterable(self):
        x_iterable = np.random.rand(3,2)
        """
        Check if an iterable object as input returns True.
        """
        assert utils.isiterable(x_iterable) == True

    def test_not_iterable(self):
        """Check if non-iterable object as input returns False"""
        x_iterable = None
        assert utils.isiterable(x_iterable) == False

class TestMakeIterable:

    """Test if no-iterable variables are correctly returned as iterable variables"""

    def test_make_iterable(self):
        """Check if a non-iterable object as input returns an iterable object"""
        x_iterable = None
        assert utils.isiterable(x_iterable) == False
        assert utils.isiterable(utils.makeiterable(x_iterable)) == True

class TestIsArray:
    """ Test if array variables can be correctly distinguished from non-array variables"""

    def test_is_array(self):
        """Check if an array object as input returns True"""
        x_array = np.random.rand(3,2)
        assert utils.isarray(x_array) == True

    def test_not_array(self):
        """Check if non-array object as input returns False"""
        x_array = "Hello World!"
        assert utils.isarray(x_array) == False

class TestInterpArray:
    """Test if interpolation of multiple time series is correctly performed"""

    def test_interp_array(self):
        """Check if interpolation of multiple time series is correctly performed"""
        x = np.random.rand(3,2)
        y = np.random.rand(3,2)
        x_new = np.random.rand(3,2)
        assert utils.isarray(utils.interp_array(x, y, x_new)) == True

# TODO: Implement other unit-tests.
