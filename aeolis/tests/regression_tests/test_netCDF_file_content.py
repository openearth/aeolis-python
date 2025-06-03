"""
This is a blackbox (regression) test for the AeoLiS model. It tests if the
output produced by the model for a given input is identical to the expected
output for the same input across code changes.

The test checks if the dimension, shape, and the array values of all the
variables in the generated output netCDF file for a given input model
configuration file is identical to the expected values in the reference
output netCDF file for the same input.

The test inputs and reference outputs used by the test are stored in the
inputs/ directory. The test inputs are organized in subdirectories based on
the dimension of the model. Each dimension has a list of test cases, each of
which is a directory that contains an input model configuration file and its
corresponding reference output netCDF file.

test_netCDF_file_content() is the main test function that contains the test logic.

path_test_input() is a helper function that does the preparation work needed
by test_netCDF_file_content().

"""


import os
import netCDF4
import numpy as np
import pytest

from aeolis import model


@pytest.fixture
def path_test_input():
    """
    Returns the path to the directory where the test inputs are stored.

    A helper function that does the preparation work needed by
    test_netCDF_file_content().

    """
    path_regression_test_dir = os.path.dirname(os.path.abspath(__file__))
    path_test_input = path_regression_test_dir + "/inputs/"
    return path_test_input


def test_netCDF_file_content(path_test_input):
    """
    Run simulation for a list of input model configurations and check if the
    dimension, shape, and the array values of all the variables in the
    generated output netCDF file is identical to the expected values
    in the reference output netCDF file for the same input.

    This function requires the path to the test input directory as a
    prerequisite. This is accomplished by the function path_test_input()
    that is marked as a fixture for this test. When pytest executes this test,
    it detects path_test_input() mentioned as a fixture, and executes it before
    executing this function. The return value of path_test_input() is passed as
    an argument to this function.

    Parameters
    ----------
    path_test_input : str
        Path to the test input directory.

    """
    for dimension in os.listdir(path_test_input):
        for case in os.listdir(path_test_input + "/" + dimension):
            os.chdir(path_test_input + "/" + dimension + "/" + case)
            tmp_model = initialize_model()
            run_model(tmp_model)
            assert_netCDF_content_consistency(dimension, case)
            delete_output_files()


def assert_netCDF_content_consistency(dimension: str, case: str):
    """
    Compare the dimension, shape, and the array values of all the variables
    in the generated output netCDF file with the expected values in the
    reference output netCDF file.
    """
    with netCDF4.Dataset("aeolis.nc", "r") as ds, netCDF4.Dataset(
        "output_expected.nc", "r"
    ) as ds_expected:
        for variable in ds.variables.values():
            # Collect the array values
            variable_value = variable[:]
            variable_value_expected = ds_expected.variables[variable.name][:]

            assert (
                variable.shape == ds_expected.variables[variable.name].shape
            ), (
                f"Array shape of the paremeter '{variable.name}' is expected to"
                " remain consistent across simulations for the same model"
                f" parameter file for {dimension} {case}"
            )
            assert variable.ndim == ds_expected.variables[variable.name].ndim, (
                f"Dimension of the parameter '{variable.name}' are expected to"
                " remain consistent across simulations for the same model"
                f" parameter file for {dimension} {case}"
            )

            # Unmask the arrays to get the underlying data
            variable_value_data = np.ma.getdata(variable_value)
            variable_value_expected_data = np.ma.getdata(
                variable_value_expected
            )

            # Compare unmasked elements using numpy.allclose
            assert np.allclose(
                variable_value_data, variable_value_expected_data
            ), (
                f"Array values of the parameter '{variable.name}' are expected"
                " to remain consistent across simulations for the same model"
                f" parameter file for {dimension} {case}"
            )


def initialize_model():
    tmp_model = model.AeoLiSRunner()
    return tmp_model


def run_model(tmp_model):
    tmp_model.run()


def delete_output_files():
    for file_name in os.listdir(os.getcwd()):
        if file_name.endswith((".log", "aeolis.nc")):
            os.remove(file_name)
