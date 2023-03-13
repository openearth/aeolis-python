"""Sandbox for brainstorming methods to compare the output netCDF files.
This is a teporary file for experimental purposes only. Finished ideas will end 
up in test_model.py module.
"""

# import netCDF4
# import numpy as np


# def verify_netCDF_content_atomic():
#     with netCDF4.Dataset("aeolis.nc", "r") as ds, netCDF4.Dataset(
#         "aeolis_expected.nc", "r"
#     ) as ds_expected:
#         zb = ds.variables["zb"]
#         zb_expected = ds_expected.variables["zb"]
#         assert np.array_equal(zb[:], zb_expected[:]) == True

#         Ct = ds.variables["Ct"]
#         Ct_expected = ds_expected.variables["Ct"]
#         assert np.array_equal(Ct[:], Ct_expected[:]) == True

#         zs = ds.variables["zs"]
#         zs_expected = ds_expected.variables["zs"]
#         assert np.array_equal(zs[:], zs_expected[:]) == True

#         w = ds.variables["w"]
#         w_expected = ds_expected.variables["w"]
#         assert np.array_equal(w[:], w_expected[:]) == True

#         mass = ds.variables["mass"]
#         mass_expected = ds_expected.variables["mass"]
#         assert np.array_equal(mass[:], mass_expected[:]) == True

#         Hs = ds.variables["Hs"]
#         Hs_expected = ds_expected.variables["Hs"]
#         assert np.array_equal(Hs[:], Hs_expected[:]) == True


# # def verify_netCDF_content_atomic_extended():
# #     with netCDF4.Dataset("aeolis.nc", "r") as ds_new, netCDF4.Dataset(
# #         "aeolis_expected.nc", "r"
# #     ) as ds_expected:
# #         # filling of variables
# #         # zb_new = ds_new.variables["zb"]

# #         for variable in ds_new.variables.values():


# #         # zb_expected = ds_expected.variables["zb"]

# #         # assert np.array_equal(zb_new[:], zb_expected[:]) == True


# def verify_netCDF_content():
#     with netCDF4.Dataset("aeolis.nc", "r") as ds:
#         variable_list = [variable for variable in ds.variables]
#         # for variable in ds.variables:
#         #     variable_list.append(variable)
#         # variable_value = ds.variables[variable]
#     with netCDF4.Dataset("aeolis_expected.nc", "r") as ds_expected:
#         for variable_expected in ds_expected.variables:
#             variable_value_expected = ds_expected.variables[variable_expected]
#             # assert (
#             #     np.array_equal(variable_value[:], variable_value_expected[:])
#             #     == True
#             # )
#     print(variable_list)


# def verify_netCDF_content_merged():
#     with netCDF4.Dataset("aeolis.nc", "r") as ds, netCDF4.Dataset(
#         "aeolis_expected.nc", "r"
#     ) as ds_expected:
#         variable_list = [variable for variable in ds.variables]
#         variable_list_expected = [
#             variable for variable in ds_expected.variables
#         ]

#         # check num parameters
#         assert len(variable_list) == len(variable_list_expected)
#         assert variable_list == variable_list_expected

#         # check dimensions
#         # for variable_value in ds.variables.values():
#         #     print(variable_value.ndim)
#         #     print(variable_value.shape)

#         # print(ds.variables["zb"].name)

#         # check array values
#         for variable in ds.variables:
#             variable_value = ds.variables[variable]
#             variable_name = variable_value.name
#             # print(
#             #     f"Variable {variable_name} has the type"
#             #     f" {type(variable_value[:])}"
#             # )

#         variable_values_list = [
#             variable_value[:] for variable_value in ds.variables.values()
#         ]
#         # print(variable_values_list)
#         variable_value_dict = dict(zip(variable_list, variable_values_list))

#         for variable in ds_expected.variables:
#             variable_value = ds_expected.variables[variable]
#             variable_name = variable_value.name
#             # print(type(variable_value[:]))

#         variable_values_list_expeted = [
#             variable_value[:]
#             for variable_value in ds_expected.variables.values()
#         ]
#         # print(variable_values_list)
#         variable_value_dict_expected = dict(
#             zip(variable_list_expected, variable_values_list_expeted)
#         )

#         # Compare array values of all variables

#         # assert variable_value_dict == variable_value_dict_expected

#         print(variable_value_dict)


# if __name__ == "__main__":
#     verify_netCDF_content_atomic()
#     # verify_netCDF_content()
#     # verify_netCDF_content_merged()
#     # verify_netCDF_content_atomic_extended()
