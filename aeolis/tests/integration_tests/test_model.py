"""
WIP: Test to check 2D and 1D base cases.
Add more description...
"""

from aeolis import model
import os
import netCDF4
import numpy as np
import pytest


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

    # Run simulation for all cases in 1D and 2D examples
    # Simulations must run from the respective case directory
    for dimension in os.listdir(path_examples_dir):
        for case in os.listdir(path_examples_dir + "/" + dimension):
            os.chdir(path_examples_dir + "/" + dimension + "/" + case)
            # print(os.getcwd())  # for debugging only
            temp_model = initialize_model()
            run_model(temp_model)
            are_output_files_generated()
            check_continuity()
            # check_num_params_in_output_file()
            delete_output_files()


def initialize_model():
    temp_model = model.AeoLiSRunner()
    return temp_model


def run_model(temp_model):
    temp_model.run()


def are_output_files_generated():
    filenames = os.listdir(os.getcwd())
    # check for correct filenames
    assert "aeolis.log" in filenames
    assert "aeolis.nc" in filenames
    filepath_aeolis_log = os.getcwd() + "/aeolis.log"
    filepath_aeolis_nc = os.getcwd() + "/aeolis.nc"
    # check for a valid file
    assert os.path.isfile(filepath_aeolis_log) == True
    assert os.path.isfile(filepath_aeolis_nc) == True


def check_continuity():
    """Check for mass conservation.

    Total erosion = outbound flux + sediment in transport
    cum_sum(zb) = cum_sum(qs_total(at boundaries)) + sum(Ct(at t_end))
    cm _sum(zb): sum of all elements in z dimension on the 50th row (50th time step)
    cum_sum(qs): sum of 1st and last values in z dimension for all 50 rows

    To do:
        Order of precision to compare the above equation
    """
    with netCDF4.Dataset("aeolis.nc", "r") as ds:
        # get spatial dimensions and bed levels
        zb = ds.variables["zb"][...]
        qs_sum = ds.variables["qs_sum"][:, :]
        zb_initial = zb[0, 0, :]
        zb_end = zb[-1, 0, :]
        delta_zb = np.sum(zb_initial, axis=0) - np.sum(zb_end, axis=0)
        print("delta zb is: ", delta_zb)
        qs_sum_trimmed = qs_sum[:, :, 0 :: ((np.ma.size(qs_sum, axis=2)) - 1)]
        cum_sum_qs = np.sum(qs_sum_trimmed, axis=2).sum(axis=0)
        print("cum_sum_qs is: ", np.ndarray.item(cum_sum_qs))
        # assert (delta_zb - np.ndarray.item(cum_sum_qs)) == pytest.approx(0.01)
        # plt.plot(zb[49, 0, :])
        # print("zb trimmed array is :", zb_trimmed)
        # print("qs_sum trimmed array is :", qs_sum_trimmed)
        # print("cum sum zb is: ", cum_sum_zb)
        # print("cum sum qs is: ", np.ndarray.item(cum_sum_qs) * 3600)


def delete_output_files():
    for file_name in os.listdir(os.getcwd()):
        if file_name.endswith((".log", ".nc")):
            os.remove(file_name)


# def test_output_file_contents():
#     check_dimensions()
#     check_num_parameters()
#     check continuity()


# def check_dimensions():
#     pass


# def check_num_parameters():
#     pass


def main():
    test_model()


if __name__ == "__main__":
    main()


## Examples directory can be manually copied in integration test directory for
## now. Real time copying can be added later if needed. Examples for testing
## can differ from the user examples.
# def prepare_inputs_for_integration_test():
#     """
#     This function copies the example folders in the root of the repository
#     to the integration test/ directory. The example folder would then serve as
#     a collection of inputs to feed to integration test.
#     """

#     # find relative path to the examples directory
#     # Should we do this via command line? by launching a temporary shell


#     # copy example folder
#     shutil.copytree()

# # #     # move to current dir
# # #     # paste the example folder here
