#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np

DIR = './'
MESH_NAME = 'mesh_test'

# /!\ In case of encoutering trouble generating the mesh: please check first that the dip angle,
# slab max depth and extent of the mesh are compatible. /!\ 

# -------------------------------
# GEOPHYSICAL REGION PARAMETERS
# -------------------------------

SLAB_MAX_DEPTH           = -700   # max depth [km]

OC_LITHOSPHERE_THICKNESS = -50    # slab thickness [km]
CO_LITHOSPHERE_THICKNESS = -70    # continental lithosphere base depth [km]
H_A1                     = -300   # asthenosphere base depth [km]
H_A2                     = -670   # lower upper-mantle base depth [km]

DIST_TRENCH_COAST        = 100
TRENCH_DEPTH             = -6

# MESH EXTENT (suppose that trench is at (0,0), so oceanic side is negative while continenale side is positive)
MESH_EXTENT_OCEAN        = -1000  # (right mesh extent) [km]
MESH_EXTENT_CONTI        =  3000  # (left mesh extent) [km]
MESH_DEPTH               = -2900  # depth of the mesh bottomkm [km]

# SLAB DIP
CONSTANT_DIP = True

if CONSTANT_DIP:
    
    SLAB_DIP_ANGLE = 22 #deg
    
else:
    
    # example of a slab profile at 31°S between 76°W et 65°W
    SLAB_XYZ_FILE = './sam_slab2_dep_02.23.18.xyz' #come from Slab2.0, Hayes 2018
    
    LAT_MIN = -31
    LAT_MAX = -31
    LON_MIN = -76
    LON_MAX = -65

# =============================================================================
# MESHING PARAMETERS
# =============================================================================

XC, YC, ZC = 101, -33, 0 #reference point (X,Y)[km] for the remeshing function
ES = 5 #size of the base element [km]

MEAN_LAT = -31 #mean latitude of the mesh
DEG_TO_KM = 111.32 * np.cos(np.deg2rad(MEAN_LAT)) #to convert degrees in mesh coordinates (roughly 1° lon ≈ 111 km * cos(lat))

