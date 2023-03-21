"""
Integration test to verify the model's behavior for 1D and 2D cases.

"""

import os
import netCDF4
import numpy as np

from aeolis import model


def test_integration():
    """
    Black-box test to verify if the model produces the same output for the
    same input.
    """
    path_integration_test_dir = os.path.dirname(os.path.abspath(__file__))
    path_inputs_dir = path_integration_test_dir + "/inputs/"

    for dimension in os.listdir(path_inputs_dir):
        for case in os.listdir(path_inputs_dir + "/" + dimension):
            os.chdir(path_inputs_dir + "/" + dimension + "/" + case)
            tmp_model = initialize_model()
            run_model(tmp_model)
            verify_netCDF_file_creation()
            verify_netCDF_file_content()
            delete_output_files()


def initialize_model():
    tmp_model = model.AeoLiSRunner()
    return tmp_model


def run_model(tmp_model):
    tmp_model.run()


def verify_netCDF_file_creation():
    filepath_aeolis_log = os.getcwd() + "/aeolis.log"
    filepath_aeolis_nc = os.getcwd() + "/aeolis.nc"
    assert os.path.isfile(filepath_aeolis_log) == True, (
        "A successful simulation should generate, upon its completion, a log"
        " file with the name 'aeolis.log' in the same directory as the"
        " configuration file"
    )
    assert os.path.isfile(filepath_aeolis_nc) == True, (
        "A successful simulation should generate, upon its completion, a netCDF"
        " file with the name 'aeolis.nc' in the same directory as the"
        " configuration file"
    )


def verify_netCDF_file_content():
    """
    Checks if the dimension, shape, and the array values of all the
    variables in the generated output netCDF file mirrors the expected
    values in the reference output netCDF file.
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
                " parameter file"
            )
            assert variable.ndim == ds_expected.variables[variable.name].ndim, (
                f"Dimension of the parameter '{variable.name}' are expected to"
                " remain consistent across simulations for the same model"
                " parameter file"
            )
            assert (
                np.array_equal(variable_value, variable_value_expected) == True
            ), (
                f"Array values of the parameter '{variable.name}' are expected"
                " to remain consistent across simulations for the same model"
                " parameter file"
            )


def delete_output_files():
    for file_name in os.listdir(os.getcwd()):
        if file_name.endswith((".log", "aeolis.nc")):
            os.remove(file_name)


def main():
    test_integration()


if __name__ == "__main__":
    main()
