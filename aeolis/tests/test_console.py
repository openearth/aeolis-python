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

import unittest
from unittest.mock import patch
from io import StringIO
import shutil
import os
import filecmp

from typer.testing import CliRunner

from aeolis.console import aeolis_app, start_aeolis_app, start_aeolis_wind_app

runner = CliRunner()


class TestAeolisCommandLineInterface(unittest.TestCase):
    def test_aeolis_help(self):
        result = runner.invoke(aeolis_app, ["--help"])
        self.assertEqual(result.exit_code, 0)

    def test_aeolis_run_help(self):
        result = runner.invoke(aeolis_app, ["run", "--help"])
        self.assertEqual(result.exit_code, 0)

    def test_aeolis_run(self):
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

    def test_aeolis_examples(self):
        """
        Test if executing `aeolis examples <path_to_directory>` on the command line creates a directory named 'aeolis-examples' in <path_to_directory> and copies aeolis/examples/sandengine_small_grids and aeolis/examples/2D/Parabolic_dune directories into it.
        """

        # path to source directory: aeolis/examples/sandengine_small_grids directory

        path_sandengine_small_grids = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "..",
            "examples",
            "sandengine_small_grids",
        )

        # path to source directory: aeolis/examples/2D/Parabolic_dune directory
        path_parabolic_dune = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "..",
            "examples",
            "2D",
            "Parabolic_dune",
        )

        # delete existing aeolis-examples directory if it exists
        if os.path.isdir("aeolis-examples"):
            shutil.rmtree("aeolis-examples")

        result = runner.invoke(
            aeolis_app,
            ["examples", "."],
        )

        self.assertEqual(result.exit_code, 0)

        # test the creation of aeolis-examples directory in the same directory from where the test was run. The test was run from aeolis/tests directory.
        self.assertTrue(os.path.isdir("aeolis-examples"))

        # test if the 'aeolis-examples' directory contains the sub directories named 'sandengine_small_grids' and 'Parabolic_dune'.
        self.assertTrue(os.path.isdir("aeolis-examples/sandengine_small_grids"))
        self.assertTrue(os.path.isdir("aeolis-examples/Parabolic_dune"))

        # path to destination: aeolis-examples/sandengine_small_grids directory.
        path_sandengine_small_grids_aeolis_examples = os.path.os.path.abspath(
            "aeolis-examples/sandengine_small_grids"
        )

        # path to destination: aeolis-examples/Parabolic_dune directory.
        path_parabolic_dune_aeolis_examples = os.path.os.path.abspath(
            "aeolis-examples/Parabolic_dune"
        )

        # compare the contents of source directory with the destination directory
        self.assertTrue(
            filecmp.dircmp(
                path_sandengine_small_grids,
                path_sandengine_small_grids_aeolis_examples,
            ).same_files
        )

        self.assertTrue(
            filecmp.dircmp(
                path_parabolic_dune, path_parabolic_dune_aeolis_examples
            ).same_files
        )

        # Do cleanup: delete the aeolis-examples directory
        shutil.rmtree("aeolis-examples")

    def test_deprecation_message_aeolis_run(self):
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

        # Testcase: Create two dummy command line arguments with the first argument as 'aeolis' and second as the path to the model configuration file. Invoking the aeolis CLI with these two arguments is expected to fail and produce an error message instructing the user to use `aeolis run <path_to_aeolis.txt>` instead.
        with patch("sys.argv", ["aeolis", path_model_config_file]):
            with patch("sys.exit") as mock_exit:
                with patch("sys.stdout", new_callable=StringIO) as mock_stdout:
                    start_aeolis_app()
                    output_error_text = mock_stdout.getvalue().strip()

        expected_error_text = (
            "Usage of the command line syntax `aeolis <path_to_aeolis.txt>` has"
            " been deprecated from v3.0.0 onwards.\n\n"
            "To run a model, use the syntax `aeolis run <path_to_aeolis.txt>`"
        )

        self.assertIn(expected_error_text, output_error_text)

        mock_exit.assert_called_with(1)

    def test_deprecation_message_aeolis_wind(self):
        """
        Test if executing `aeolis-wind <path_to_wind.txt>` on the command line produces an error message telling the user to use `aeolis wind <path_to_wind.txt>` instead.
        """
        path_model_dir = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "regression_tests",
            "inputs",
            "1D",
            "case1_small_waves",
        )

        path_wind_config_file = os.path.join(path_model_dir, "wind.txt")

        # Testcase: Create two dummy command line arguments with the first argument as 'aeolis-wind' and second as the path to the wind configuration file. Invoking aeolis-wind from the command line should print a message on the console reporting the deprecated usage of the syntax `aeolis-wind` instructing the user to use `aeolis wind <path_to_wind.txt>` instead.
        with patch("sys.argv", ["aeolis-wind", path_wind_config_file]):
            with patch("sys.exit") as mock_exit:
                with patch("sys.stdout", new_callable=StringIO) as mock_stdout:
                    start_aeolis_wind_app()
                    output_error_text = mock_stdout.getvalue().strip()

        expected_error_text = (
            "Usage of the command line syntax `aeolis-wind <path_to_wind.txt>`"
            " has been deprecated from v3.0.0 onwards.\n\nTo run the wind"
            " module, use the syntax `aeolis wind <path_to_wind.txt>`"
        )

        self.assertIn(expected_error_text, output_error_text)

        mock_exit.assert_called_with(1)


if __name__ == "__main__":
    unittest.main()
