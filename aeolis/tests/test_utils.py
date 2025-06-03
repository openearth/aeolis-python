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
import pytest

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

@pytest.fixture
def y():
    """Return a 2D array"""
    return np.random.rand(3,3)

class TestInterpArray:
    """Test if interpolation of multiple time series is correctly performed"""

    def test_interp_array(self, y):
        pass
        #TODO: Fix this test
        # """Check if interpolation of multiple time series is correctly performed"""
        # x = np.random.rand(3,2)
        # y = np.random.rand(3,2)
        
        # assert utils.isarray(utils.interp_array(x, y, x_new)) == True


class TestInterpCircular:
    """Test the interp_circular function"""

    def test_value_error_raised(self):
        """Test if a value error is raised when xp and f have different lengths"""
        
        x = np.random.rand(10)
        xp = np.random.rand(10)
        fp = np.random.rand(9)        

        with pytest.raises(ValueError):
            utils.interp_circular(x, xp, fp)


class TestPreventTinyNegatives:
    """Test the prevent tiny negatives function"""

    def test_prevent_tiny_negatives(self):
        """Test if tiny negative values in an array are replace a the replacement value"""
        
        x = np.array([1, 2, -1e-10, 3, -1e-10, 4, 5])
        x = utils.prevent_tiny_negatives(x, max_error=1e-9, replacement=-0.01)

        assert np.all(x >= -0.01)


# TODO: Implement other unit-tests.
