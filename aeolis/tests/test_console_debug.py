import os
from aeolis.console import aeolis_app
from typer.testing import CliRunner

runner = CliRunner()


def debug_run_command():
    model_config_file_dir = os.path.join(
        "regression_tests", "inputs", "1D", "case1_small_waves"
    )
    abs_path = os.path.dirname(os.path.abspath(__file__))

    path = abs_path + "/" + model_config_file_dir + "/aeolis.txt"
    print(path)
    # run the model
    result = runner.invoke(
        aeolis_app,
        ["run", abs_path + "/" + model_config_file_dir + "/aeolis.txt"],
    )

    print(result.stdout)

    # # check if aeolis.log and aeolis.nc got created in model config file path
    # self.assertTrue(
    #     os.path.isfile(os.path.join(model_config_file_path, "aeolis.log"))
    # )
    # self.assertTrue(
    #     os.path.isfile(os.path.join(model_config_file_path, "aeolis.nc"))
    # )

    # delete the generate log and nc files
    os.remove(os.path.join(abs_path, model_config_file_dir, "aeolis.log"))
    os.remove(os.path.join(abs_path, model_config_file_dir, "aeolis.nc"))


if __name__ == "__main__":
    debug_run_command()
