'''This file is part of AeoLiS.
   
AeoLiS is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
   
AeoLiS is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
   
You should have received a copy of the GNU General Public License
along with AeoLiS.  If not, see <http://www.gnu.org/licenses/>.
   
AeoLiS  Copyright (C) 2015 Bas Hoonhout

bas.hoonhout@deltares.nl         b.m.hoonhout@tudelft.nl
Deltares                         Delft University of Technology
Unit of Hydraulic Engineering    Faculty of Civil Engineering and Geosciences
Boussinesqweg 1                  Stevinweg 1
2629 HVDelft                     2628CN Delft
The Netherlands                  The Netherlands

'''


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

    
