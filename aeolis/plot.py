import netCDF4
import matplotlib.pyplot as plt


def profile(outputfile, var, ix=-1, subplots_kw={'figsize':(10,4)}):
    '''Plots profile
    
    Parameters
    ----------
    outputfile : str
        Name of netCDF4 output file
    var : str
        Name of spatial grid
    ix : int, optional
        Time slice (default: -1)
    subplots_kw : dict
        Keyword options to subplots function

    '''

    with netCDF4.Dataset(outputfile, 'r') as ds:

        x = ds.variables['x'][...]
        y = ds.variables[var][ix,...]

        fig, axs = plt.subplots(**subplots_kw)
        axs.plot(x, y, '-k')
        axs.set_xlabel('x [m]')
        axs.set_title(var)

        return fig, axs

    
