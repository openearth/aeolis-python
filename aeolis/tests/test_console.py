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

from aeolis.console import aeolis_app, start_aeolis_app

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
        try:
            os.remove(path_model_dir + "/aeolis.log")
            os.remove(path_model_dir + "/aeolis.nc")
        except Exception as exception:
            print(f"Failed to delete files: {exception}")

    def test_examples_command(self):
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

        # Testcase: Create two dummy command line arguments with the first argument as 'aeolis' and second as the path to the model configuration file. Invoking the aeolis CLI with these two arguments is expected to fail and produce an error message instructing the user to use `aeolis run <path_to_configfile>` instead.
        with patch("sys.argv", ["aeolis", path_model_config_file]):
            with patch("sys.exit") as mock_exit:
                with patch("sys.stdout", new_callable=StringIO) as mock_stdout:
                    start_aeolis_app()
                    output_error_text = mock_stdout.getvalue().strip()

        expected_error_text = (
            "Error: You entered an incorrect command.\n To run a model, type:"
            " `aeolis run <path_to_configfile>`\n From aeolis v2.2.0 onwards,"
            " the command run needs to be passed to run a model."
        )

        self.assertIn(expected_error_text, output_error_text)

        mock_exit.assert_called_with(1)


if __name__ == "__main__":
    unittest.main()
