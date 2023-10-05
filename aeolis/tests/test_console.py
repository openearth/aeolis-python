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

from aeolis.console import aeolis_app
import unittest
import shutil
from typer.testing import CliRunner
import os

runner = CliRunner()


class TestAeolisCommandLineInterface(unittest.TestCase):
    def test_help_command(self):
        result = runner.invoke(aeolis_app, ["--help"])
        self.assertEqual(result.exit_code, 0)

    def test_run_help_command(self):
        result = runner.invoke(aeolis_app, ["run", "--help"])
        self.assertEqual(result.exit_code, 0)

    def test_run_command(self):
        model_config_file_dir = os.path.join(
            "regression_tests",
            "inputs",
            "1D",
            "case1_small_waves",
        )

        # path to parent directory
        cwd = os.path.dirname(os.path.abspath(__file__))

        result = runner.invoke(
            aeolis_app,
            ["run", cwd + "/" + model_config_file_dir + "/aeolis.txt"],
        )

        self.assertEqual(result.exit_code, 0)

        # check if aeolis.log and aeolis.nc got created in model config file path
        self.assertTrue(
            os.path.isfile(
                os.path.join(cwd, model_config_file_dir, "aeolis.log")
            )
        )
        self.assertTrue(
            os.path.isfile(
                os.path.join(cwd, model_config_file_dir, "aeolis.nc")
            )
        )
        # delete the generate log and nc files
        os.remove(os.path.join(cwd, model_config_file_dir, "aeolis.log"))
        os.remove(os.path.join(cwd, model_config_file_dir, "aeolis.nc"))

    # def test_missing_run_command(self):
    #     model_config_file_dir = os.path.join(
    #         "regression_tests",
    #         "inputs",
    #         "1D",
    #         "case1_small_waves",
    #     )
    #     result = runner.invoke(
    #         aeolis_app, [model_config_file_dir + "/aeolis.txt"]
    #     )

    #     self.assertIn(
    #         "Error: You entered an incorrect command.\n To run a model,",
    #         result.stdout,
    #     )

    def test_examples_command(self):
        # delete existing aeolis-examples directory if it exists
        if os.path.isdir("aeolis-examples"):
            shutil.rmtree("aeolis-examples")

        result = runner.invoke(aeolis_app, ["examples", "."])
        self.assertEqual(result.exit_code, 0)

        # test the creation of aeolis-examples directory in the same directory from where the command was run
        self.assertTrue(os.path.isdir("aeolis-examples"))

        # do cleanup;delete the aeolis-examples directory
        shutil.rmtree("aeolis-examples")


if __name__ == "__main__":
    unittest.main()
