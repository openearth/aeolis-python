"""
WIP: Test to verify the model's expected behavior for 1D and 2D base cases.
Add more description...
"""

from aeolis import model
import os
import netCDF4
import numpy as np


def test_model():
    """
    Run simulations using example configurations and run test cases to compare
    the results against expected results.
    Loop through input configurations, run simulation and compare results by
    running testcases.
    """

    # Locate the path for the directory containing the input configurations to test
    path_integration_test_dir = os.path.dirname(os.path.abspath(__file__))
    path_examples_dir = path_integration_test_dir + "/examples/"

    # Loop through all cases in 1D and 2D and run the simulation for each case
    for dimension in os.listdir(path_examples_dir):
        for case in os.listdir(path_examples_dir + "/" + dimension):
            # Simulations must run from the respective case directory
            os.chdir(path_examples_dir + "/" + dimension + "/" + case)

            temp_model = initialize_model()

            run_model(temp_model)

            # Check if output netCDF files are created successfully
            verify_netCDF_file_creation()

            # The contents of the output netCDF file should be consistent with
            # the contents of the reference netCDF file generated using the same input
            verify_netCDF_content()

            # Check if continuity/mass conservation holds
            # check_continuity()

            # Check number of parameters in the output netCDF file
            # check_num_params_netcdf()

            # Generated output files become obsolete after all the tests cases are run
            delete_output_files()


def initialize_model():
    temp_model = model.AeoLiSRunner()
    return temp_model


def run_model(temp_model):
    temp_model.run()


def verify_netCDF_file_creation():
    """Verify the successful generation of the output netCDF file and the log
    file post the completion of the simulation.
    """
    filenames = os.listdir(os.getcwd())
    # check for correct filenames
    assert "aeolis.log" in filenames
    assert "aeolis.nc" in filenames
    filepath_aeolis_log = os.getcwd() + "/aeolis.log"
    filepath_aeolis_nc = os.getcwd() + "/aeolis.nc"
    # check for a valid file
    assert os.path.isfile(filepath_aeolis_log) == True
    assert os.path.isfile(filepath_aeolis_nc) == True


def verify_netCDF_content():
    """Compare the values of all the variables in the output netCDF file with
    the expected values in the reference output file.

    To do: Remove duplicated code for each variable by using a for loop
    """
    with netCDF4.Dataset("aeolis.nc", "r") as ds, netCDF4.Dataset(
        "aeolis_expected.nc", "r"
    ) as ds_expected:
        zb = ds.variables["zb"]
        zb_expected = ds_expected.variables["zb"]
        assert np.array_equal(zb[:], zb_expected[:]) == True

        Ct = ds.variables["Ct"]
        Ct_expected = ds_expected.variables["Ct"]
        assert np.array_equal(Ct[:], Ct_expected[:]) == True

        zs = ds.variables["zs"]
        zs_expected = ds_expected.variables["zs"]
        assert np.array_equal(zs[:], zs_expected[:]) == True

        w = ds.variables["w"]
        w_expected = ds_expected.variables["w"]
        assert np.array_equal(w[:], w_expected[:]) == True

        mass = ds.variables["mass"]
        mass_expected = ds_expected.variables["mass"]
        assert np.array_equal(mass[:], mass_expected[:]) == True

        Hs = ds.variables["Hs"]
        Hs_expected = ds_expected.variables["Hs"]
        assert np.array_equal(Hs[:], Hs_expected[:]) == True

        zs = ds.variables["zs"]
        zs_expected = ds_expected.variables["zs"]
        assert np.array_equal(zs[:], zs_expected[:]) == True


def delete_output_files():
    for file_name in os.listdir(os.getcwd()):
        if file_name.endswith((".log", "aeolis.nc")):
            os.remove(file_name)


# def check_dimensions():
#     pass


# def check_num_parameters():
#     pass


def main():
    test_model()


if __name__ == "__main__":
    main()
