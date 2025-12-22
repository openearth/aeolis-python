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

# TODO: Implement other unit-tests.
