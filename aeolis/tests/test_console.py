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

from aeolis.console import aeolis_app, start_aeolis_app
import unittest
import sys
from io import StringIO
from unittest.mock import patch
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
        """
        Test if executing `aeolis run <path_to_configfile>` on the command line runs the simulaition successfully.

        Use the model configuration file in regression_tests/inputs/1D/case1_small_waves directory as the test input.

        To verify if simulation was completed successfully, check if aeolis.log file and a aeolis.nc file are generated in the same directory as the model configuration file.

        """
        path_model_dir = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "regression_tests",
            "inputs",
            "1D",
            "case1_small_waves",
        )

        path_model_config_file = os.path.join(path_model_dir, "aeolis.txt")

        result = runner.invoke(
            aeolis_app,
            ["run", path_model_config_file],
        )

        self.assertEqual(result.exit_code, 0)

        # Check if aeolis.log and aeolis.nc got generated in the directory containing the model configuration file
        self.assertTrue(os.path.isfile(path_model_dir + "/aeolis.log"))
        self.assertTrue(os.path.isfile(path_model_dir + "/aeolis.nc"))

        # delete the generated aeolis.log and aeolis.nc files
        os.remove(path_model_dir + "/aeolis.log")
        os.remove(path_model_dir + "/aeolis.nc")

    def test_run_command_with_invalid_config_file(self):
        result = runner.invoke(
            aeolis_app,
            ["run", "invalid_config_file.txt"],
        )

        self.assertEqual(result.exit_code, 1)

    def test_examples_command(self):
        # delete existing aeolis-examples directory if it exists
        if os.path.isdir("aeolis-examples"):
            shutil.rmtree("aeolis-examples")

        result = runner.invoke(aeolis_app, ["examples", "."])
        self.assertEqual(result.exit_code, 0)

        # test the creation of aeolis-examples directory in the same directory from where the command was run
        self.assertTrue(os.path.isdir("aeolis-examples"))

        # Do cleanup: delete the aeolis-examples directory
        shutil.rmtree("aeolis-examples")

    def test_missing_run_command(self):
        """
        Test if executing `aeolis <path_to_configfile>` on the command line produces an error message telling the user to use `aeolis run <path_to_configfile>` instead.
        """
        path_model_dir = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "regression_tests",
            "inputs",
            "1D",
            "case1_small_waves",
        )

        path_model_config_file = os.path.join(path_model_dir, "aeolis.txt")
        with patch("sys.argv", [path_model_config_file]):
            with self.assertRaises(SystemExit) as cm:
                result = start_aeolis_app()
                self.assertNotEqual(cm.exception.code, 0)
                self.assertIn(
                    (
                        "Error: You entered an incorrect command.\n To run a"
                        " model, type: aeolis run <path_to_configfile>\n From"
                        " aeolis v2.2.0 onwards, the command run needs to be"
                        " passed to run a model.\n"
                    ),
                    result,
                ),


if __name__ == "__main__":
    unittest.main()
