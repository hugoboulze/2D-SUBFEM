#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 17:41:27 2021

@author: boulze, garaud

Create the surface slab as a csv file.
Extrapole the surface of the slab in Patagonia.

"""

# =============================================================================
import os
import sys
sys.path.append('./../..')
import numpy as np
import matplotlib.pyplot as plt
from tools import PtLatiLongiToXYZ
# =============================================================================

DIR = './'

def xyz_to_grid(raw_slab):

    '''
    raw_slab: sam_slab2_dep_02.23.18.xyz
    area : list
        [long_min,long_max,lat_min,lat_max]

    step=0.05 par defaut

    Returns
    -------
    None.

    '''

    raw_slab = np.genfromtxt(raw_slab, delimiter=',', dtype=float)

    depth_grid = raw_slab[:,2]
    depth_grid = np.reshape(depth_grid , (1281,581))

    return depth_grid


def extract_slab(grid_slab, lon_min, lon_max, lat_min, lat_max, max_depth, res, start_with_trench=False, plot=False):
    
    '''
    extract a piece of the slab in the grid according to lon_min/max, lat_min/max,  max_depth and res
    
    grid_slab: = [[lon_i, lat_i, depth_i], ...]
    step = 0.05 by default in Slab2.0

    '''

    lon_0 = 274
    lat_0 = 15


    grid_res = 0.05

    step = np.int64(res/grid_res)

    LAT_SLAB = []
    LON_SLAB = []
    DEPTH_SLAB = []

    #on parcourt la grille en latitudes
    for i in np.arange(0, grid_slab.shape[0], step):

        lat = lat_0 - grid_res*i

        #quand la latitude est bonne on prend la premiere valeur non 'nan' en longitude de maniere a recup la fosse
        if lat_min <= lat <= lat_max:

            for j in np.arange(0, grid_slab.shape[1], 1):

                # print(grid_slab[i,j])
                if not np.isnan(grid_slab[i,j]):
                    # start_lon = lon_0 + j*grid_res
                    start_lon_idx = j
                    break

            for l in range(start_lon_idx, grid_slab.shape[1], step):

                if not np.isnan(grid_slab[i,l]):

                    lon = lon_0 + l*grid_res
                    
                    if grid_slab[i,l] > -max_depth:

                        LAT_SLAB.append(round(lat,4))
                        LON_SLAB.append(round(lon-360,4))
                        DEPTH_SLAB.append(round(grid_slab[i,l],4))

    if plot:
        plt.figure()
        plt.axis('equal')
        plt.xlabel('Longitude [deg]')
        plt.ylabel('Latitude [deg]')
        plt.scatter(LON_SLAB, LAT_SLAB, c=DEPTH_SLAB)
        plt.show()
    

    return np.array(LON_SLAB), np.array(LAT_SLAB), np.array(DEPTH_SLAB)

