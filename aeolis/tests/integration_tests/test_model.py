"""
WIP: Test to verify the model's expected behavior for 1D and 2D base cases.
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

    # Loop through all cases in 1D and 2D and run the simulation for each case
    for dimension in os.listdir(path_examples_dir):
        for case in os.listdir(path_examples_dir + "/" + dimension):
            # Simulations must run from the respective case directory
            os.chdir(path_examples_dir + "/" + dimension + "/" + case)

            temp_model = initialize_model()

            run_model(temp_model)

            # Check if output netCDF files are created successfully
            verify_netCDF_file_creation()

            # Check if continuity/mass conservation holds
            check_continuity()

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
        Ct = ds.variables["Ct"][:, :]
        pickup_sum = ds.variables["pickup_sum"][:, :]
        print("Ct shape is", Ct.shape)
        zb_initial = zb[0, 0, :]
        zb_end = zb[-1, 0, :]
        delta_zb = np.sum(zb_end, axis=0) - np.sum(zb_initial, axis=0)
        print("delta zb is: ", delta_zb)
        q_in = qs_sum[:, :, 0]
        q_out = qs_sum[:, :, -1]
        q_in_sum = np.sum(q_in, axis=0)
        q_out_sum = np.sum(q_out, axis=0)
        cum_sum_qs = q_in_sum - q_out_sum
        print("cum_sum_qs is: ", cum_sum_qs)
        Ct_final = Ct[-1, 0, :, 0]
        Ct_final_sum = np.sum(Ct_final, axis=0)
        print("Ct final sum is: ", Ct_final_sum)
        # plt.plot(zb[49, 0, :])
        density = 2650  # global constant?
        porosity = 0.4  # global constant?
        delta_mass = (delta_zb * density * 0.5) / (1 - porosity)
        delta_pos = (
            np.sum(np.absolute(zb_end - zb_initial), axis=0) / 2 * density * 0.5
        ) / (1 - porosity)
        print(
            "pickup sum", np.sum(np.sum(pickup_sum[:, 0, :, 0], axis=0), axis=0)
        )
        print("delta pos ", delta_pos)
        print("delta mass", delta_mass)


def delete_output_files():
    for file_name in os.listdir(os.getcwd()):
        if file_name.endswith((".log", ".nc")):
            os.remove(file_name)


# def check_dimensions():
#     pass


# def check_num_parameters():
#     pass


def main():
    test_model()


if __name__ == "__main__":
    main()
