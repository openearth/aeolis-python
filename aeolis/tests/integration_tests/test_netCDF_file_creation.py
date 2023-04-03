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


def test_netCDF_file_creation(path_test_input):
    for dimension in os.listdir(path_test_input):
        for case in os.listdir(path_test_input + "/" + dimension):
            os.chdir(path_test_input + "/" + dimension + "/" + case)
            tmp_model = initialize_model()
            run_model(tmp_model)
            verify_netCDF_file_creation()
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


def delete_output_files():
    for file_name in os.listdir(os.getcwd()):
        if file_name.endswith((".log", "aeolis.nc")):
            os.remove(file_name)
