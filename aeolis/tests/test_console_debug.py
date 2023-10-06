# import os
# from aeolis.console import aeolis_app, start_aeolis_app
# from unittest.mock import patch
# import unittest
# from typer.testing import CliRunner
# from io import StringIO

# runner = CliRunner()


# def debug_run_command():
#     path_model_config_file = os.path.join(
#         os.path.dirname(os.path.abspath(__file__)),
#         "regression_tests",
#         "inputs",
#         "1D",
#         "case1_small_waves",
#         "aeolis.txt",
#     )

#     print(path_model_config_file)

#     # run the model
#     result = runner.invoke(
#         aeolis_app,
#         ["run", path_model_config_file],
#     )

#     print(result.stdout)

#     # check if the model produced a log file and a netCDF file
#     assert os.path.isfile(
#         os.path.join(
#             os.path.dirname(os.path.abspath(__file__)),
#             "regression_tests",
#             "inputs",
#             "1D",
#             "case1_small_waves",
#             "aeolis.log",
#         )
#     )
#     assert os.path.isfile(
#         os.path.join(
#             os.path.dirname(os.path.abspath(__file__)),
#             "regression_tests",
#             "inputs",
#             "1D",
#             "case1_small_waves",
#             "aeolis.nc",
#         )
#     )

#     # delete the generated log and nc files
#     os.remove(
#         os.path.join(
#             os.path.dirname(os.path.abspath(__file__)),
#             "regression_tests",
#             "inputs",
#             "1D",
#             "case1_small_waves",
#             "aeolis.log",
#         )
#     )
#     os.remove(
#         os.path.join(
#             os.path.dirname(os.path.abspath(__file__)),
#             "regression_tests",
#             "inputs",
#             "1D",
#             "case1_small_waves",
#             "aeolis.nc",
#         )
#     )


# class TestConsole(unittest.TestCase):
#     def test_start_aeolis_app(self):
#         """
#             Test case: Execute aeolis/console.py as a script with only two arguments, for example python aeolis/console.py aeolis <path_config file>. 1st argument is aeolis and 2nd argument is a path to a model configuration file. DO not use runner invoke. Use mock patch to provide dummy arguments
#         assert for the print message in stdout: Following message should get printed on console ""Error: You entered an incorrect command.\n To run a model, type:"
#                 " aeolis run <path_to_configfile>\n From aeolis v2.2.0 onwards, the"
#                 " command run needs to be passed to run a model.\n""
#         """
#         path_model_dir = os.path.join(
#             os.path.dirname(os.path.abspath(__file__)),
#             "regression_tests",
#             "inputs",
#             "1D",
#             "case1_small_waves",
#         )

#         path_model_config_file = os.path.join(path_model_dir, "aeolis.txt")

#         with patch("sys.argv", ["aeolis", path_model_config_file]):
#             start_aeolis_app()
#             self.assertNotEqual(start_aeolis_app(), 0)


# if __name__ == "__main__":
#     unittest.main()
