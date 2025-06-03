"""
This is a blackbox (regression) test for the AeoLiS model. It tests whether the
model produces a netCDF file and a log file upon the completion of the
simulation.

The test inputs used by the test are stored in the inputs/ directory. The test
inputs are organized in subdirectories based on the dimension of the model.
Each dimension has a list of test cases, each of which is a directory that
contains an input model configuration file.

test_netCDF_file_creation() is the main test function that contains the test logic.

path_test_input() is a helper function that does the preparation work needed
by test_netCDF_file_creation().

"""


import os
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


def test_output_files_generation(path_test_input):
    """
    Run simulation for a list of input model configurations and check if the
    model produces a netCDF file and a log file upon the completion of the
    simulation.


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
            assert_output_files_existence()
            delete_output_files()


def initialize_model():
    tmp_model = model.AeoLiSRunner()
    return tmp_model


def run_model(tmp_model):
    tmp_model.run()


def assert_output_files_existence():
    """
    Checks whether the model produces a netCDF file and a log file upon the
    completion of the simulation.

    """
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


def delete_output_files():
    for file_name in os.listdir(os.getcwd()):
        if file_name.endswith((".log", "aeolis.nc")):
            os.remove(file_name)
