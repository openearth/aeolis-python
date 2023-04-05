"""
Test if the dimension, shape, and the array values of all the variables in the generated output netCDF file 
mirrors the expected values in the reference output netCDF file.
"""


import os
import netCDF4
import numpy as np
import pytest

from aeolis import model


@pytest.fixture
def path_test_input():
    path_integration_test_dir = os.path.dirname(os.path.abspath(__file__))
    path_test_input = path_integration_test_dir + "/inputs/"
    return path_test_input


def test_netCDF_file_content(path_test_input):
    for dimension in os.listdir(path_test_input):
        for case in os.listdir(path_test_input + "/" + dimension):
            os.chdir(path_test_input + "/" + dimension + "/" + case)
            tmp_model = initialize_model()
            run_model(tmp_model)
            verify_netCDF_file_content(dimension, case)
            delete_output_files()


def verify_netCDF_file_content(dimension: str, case: str):
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

            # Unmask the arrays and collect the underlying data
            variable_value_data = np.ma.getdata(variable_value)
            variable_value_expected_data = np.ma.getdata(
                variable_value_expected
            )

            # Extract the mask
            variable_value_mask = np.ma.getmask(variable_value)
            variable_value_expected_mask = np.ma.getmask(
                variable_value_expected
            )

            # Compare non-masked elements using numpy.allclose

            # Using numpy.allclose to compare two masked arrays element-wise generates the
            # error "TypeError: ufunc 'bitwise_and' not supported for the input types,
            # and the inputs could not be safely coerced to any supported types according
            # to the casting rule safe". Unmasking the arrays converts them to regular
            # numpy arrays, which are compatible with element-wise comparisons of numpy.allclose.

            assert np.allclose(
                variable_value_data[~variable_value_mask],
                variable_value_expected_data[~variable_value_expected_mask],
                rtol=1e-5,
                atol=1e-8,
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
