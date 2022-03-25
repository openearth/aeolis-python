# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 11:04:27 2022

@author: cijzendoornvan
"""

import numpy as np
from aeolis.inout import *

def create_4Dgraindist(nx, ny, nl, nf, bedcomp_settings):
    '''Create grain size distribution matrix

    The grain size distribution matrix is created in which the weight of each fraction
	is defined for each grid cell, so for all x- and y-coordinates and layers.

    Parameters
    ----------
    nx : int
        Number of grid cells in x-direction
    ny : int
        Number of grid cells in y-direction
    nl : int
        Number of vertical layers
    nf : int
        Number of grain size fractions
    bedcomp_settings : numpy.ndarray
        Parameters for creation of grain distribution matrix

    Returns
    -------
    grain_dist : numpy.ndarray
        Array or matrix containing the grain distribution

    '''
    #nx = nx +1
    #ny = ny +1
	
    # Prepare grain dist matrix
    grain_dist = np.zeros((ny, nx, nl, nf))
	
	# Create grain distribution matrix with horizontally uniform layers
    if bedcomp_settings['uni_vert'] == False and bedcomp_settings['uni_x'] == True and bedcomp_settings['uni_y'] == True:
        for y in range(ny):
            for x in range(nx):
                for i in range(nl):
                    print(bedcomp_settings['dist1_layer'])
                    if i+1 == bedcomp_settings['dist1_layer']: #i+1 in settings['dist1_layer'] should be used if len(settings['dist1_layer']) > 1
                        grain_dist[y, x, i,:] = bedcomp_settings['dist1']
                    else:
                        grain_dist[y, x, i,:] = bedcomp_settings['base_dist']
						
	# Create grain distribution matrix with cross-shore gradient, fine-coarse-fine
    if bedcomp_settings['uni_vert'] == True and bedcomp_settings['uni_x'] == False and bedcomp_settings['uni_y'] == True and bedcomp_settings['gradient'] == True and bedcomp_settings['mid_gradient'] == True:
        for y in range(ny):
            for x in range(nx):
                for i in range(nl):
                    for f in range(nf):
                        gradient1 = np.linspace(bedcomp_settings['start_grad'][f], bedcomp_settings['mid_grad'][f], int(nx*bedcomp_settings['mid_location'])+1, endpoint=True)
                        gradient2 = np.linspace(bedcomp_settings['mid_grad'][f], bedcomp_settings['end_grad'][f], int(nx*(1-bedcomp_settings['mid_location'])), endpoint=True)
                        gradient = np.concatenate((gradient1, gradient2))
                        grain_dist[y, :, i, f] = gradient

	# Create grain distribution matrix with cross-shore gradient
    if bedcomp_settings['uni_vert'] == True and bedcomp_settings['uni_x'] == False and bedcomp_settings['uni_y'] == True and bedcomp_settings['gradient'] == True and bedcomp_settings['mid_gradient'] == False:
        for y in range(ny):
            for x in range(nx):
                for i in range(nl):
                    for f in range(nf):
                        grain_dist[y, :, i, f] = np.linspace(bedcomp_settings['start_grad'][f], bedcomp_settings['end_grad'][f], nx, endpoint=True)

	# Create grain distribution matrix with vertically uniform patches
    if bedcomp_settings['uni_vert'] == True and bedcomp_settings['uni_x'] == False and bedcomp_settings['uni_y'] == True and bedcomp_settings['patches'] == True:
        patchmask = np.zeros((nx))
        patch_factor = int(bedcomp_settings['patch_perc']*nx)
        patchmask[:patch_factor] = 1
        np.random.shuffle(patchmask)  		
        for y in range(ny):
            for x in range(nx):
                for i in range(nl):
                    if patchmask[x] == 1:
                        grain_dist[y, x, i,:] = bedcomp_settings['patch_dist']
                    else:
                        grain_dist[y, x, i,:] = bedcomp_settings['base_dist']

	# Create grain distribution matrix with vertically varying patches
    if bedcomp_settings['uni_vert'] == False and bedcomp_settings['uni_x'] == False and bedcomp_settings['uni_y'] == True and bedcomp_settings['patches'] == True:
        patchmask = np.zeros((nx, nl))
        patch_factor = int(bedcomp_settings['patch_perc']*nx)
        patchmask[:patch_factor, :] = 1
        for i in range(nl):
            np.random.shuffle(patchmask[:, i])
        for y in range(ny):
            for x in range(nx):
                for i in range(nl):
                    if patchmask[x, i] == 1:
                        grain_dist[y, x, i,:] = bedcomp_settings['patch_dist']
                    else:
                        grain_dist[y, x, i,:] = bedcomp_settings['base_dist']

    return grain_dist


def graindist_to_mass(nx, ny, nl, nf, grain_dist, model_settings):
    '''Convert grain size distribution matrix to mass distribution matrix

    The mass distribution per fraction is created by calculating the mass per 
    x, y, layer and fraction.
	

    Parameters
    ----------
    nx : int
        Number of grid cells in x-direction
    ny : int
        Number of grid cells in y-direction
    nl : int
        Number of vertical layers
    nf : int
        Number of grain size fractions
    grain_dist : numpy.ndarray
        Array or matrix containing the grain distribution

    Returns
    -------
    mass : numpy.ndarray
        Array or matrix containing the mass distribution

    '''

    if ny != grain_dist.shape[0] or nx != grain_dist.shape[1] or nl != grain_dist.shape[2] or nf != grain_dist.shape[3]:
        logger.log_and_raise('Grain size distribution not assigned for each grid cell and fraction.', exc=ValueError)	

    # Prepare mass matrix
    mass = np.zeros((ny, nx, nl, nf))
        
    for y in range(ny):
        for x in range(nx):
            for i in range(nl):
                gs = makeiterable(grain_dist[y, x, i,:])
                gs = gs / np.sum(gs)
                for j in range(nf):
                    mass[y,x,i,j] = model_settings['rhog'] * (1. - model_settings['porosity']) * model_settings['layer_thickness'] * gs[j]

    return mass

def conv4Dto2D(arr_4D, nl, nf):
    '''Create 2D array from 4D array

    Convert a 4D array to a 2D array using numpy's reshape.
    Needed to save the array in a text file. Created to save the grain distribution. 

    Parameters
    ----------
    arr_4D : numpy.ndarray
        4D array

    Returns
    -------
    arr_2D : numpy.ndarray
        2D array
    '''
    
    arr_2D = arr_4D.reshape(-1, nl*nf)
    
    return arr_2D


def conv2Dto4D(arr_2D, ny, nx, nl, nf):
    '''Create 2D array from 4D array

    Convert a 4D array to a 2D array using numpy's reshape.
    Needed to save the array in a text file. Created to save the grain distribution. 

    Parameters
    ----------
    arr_2D : numpy.ndarray
        2D array

    Returns
    -------
    arr_4D : numpy.ndarray
        4D array
    '''
    
    arr_4D = arr_2D.reshape((ny, nx, nl, nf))	
    
    return arr_4D


def save_2Dtotxt(arr_2D, save_dir, filename):
    '''Save 2D array as txt file

    Save data in 2D array in text file in the assigned directory. 

    Parameters
    ----------
    arr_2D : numpy.ndarray
        2D array

    '''
    
    np.savetxt(save_dir + filename, arr_2D, delimiter= ' ')

    return
    

model_dir = "C:/Users/cijzendoornvan/OneDrive - Delft University of Technology/Documents/DuneForce/AEOLIS/aeolis-python/examples/grainsizevariations"
save_dir = model_dir #+ "/../data"

model_code = 'aeolis_horizontalgradient2'

model_settings_file = model_dir + "/" + model_code + ".txt"
model_settings = read_configfile(model_settings_file) # text file with model run settings

bedcomp_settings_file = model_dir + "/data/gs_settings_horzgrad2.txt" #../data/
bedcomp_settings = read_configfile(bedcomp_settings_file) # text file with bed comp creation settings

xgrid_file = np.loadtxt(model_dir + '/' + model_settings['xgrid_file'])
if xgrid_file.ndim == 2:
    ny, nx = xgrid_file.shape         
else:
    nx = len(xgrid_file)
    ny = 1

if '.txt' in model_settings['grain_size'][0]:
    nf = len(np.loadtxt(model_dir + '/' + model_settings['grain_size'][0]))
else:
    nf = len(model_settings['grain_size'])
nl = model_settings['nlayers']

graindist_4D = create_4Dgraindist(nx, ny, nl, nf, bedcomp_settings)
mass_4D = graindist_to_mass(nx, ny, nl, nf, graindist_4D, model_settings)
mass_2D = conv4Dto2D(mass_4D, nl, nf)

mass_file = '/data/bedcomp_' + model_code + '.txt'
save_2Dtotxt(mass_2D, save_dir, mass_file)

