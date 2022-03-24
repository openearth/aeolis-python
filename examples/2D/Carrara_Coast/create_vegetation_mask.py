# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 14:28:53 2020

@author: tspak
"""

import numpy as np

"""
create a vegetation mask file for AeoLiS, this mask can be implemented into AeoLiS by defining:
    threshold_mask = vegetation_mask.txt
in the configuration file

visit :
    https://aeolis.readthedocs.io/en/latest/sourcecode.html#utils.apply_mask
to see how the mask works

input parameters:
shape_file (str): name of a spatial grid file which determines the shape of the mask file 

"""

shape_file = 'z.txt'

# open shape_file and create a vegetation mask with the same shape
shape = np.loadtxt(shape_file)
vegetation_mask = (1+0j) * np.ones(shape.shape)

#you can make alterations to the vegetation mask below
vegetation_mask[0,0] = 1+0j

# save vegetation mask
np.savetxt('vegetation_mask.txt', vegetation_mask)