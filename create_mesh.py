#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import gmsh
import numpy as np
import create_surface_slab as surfSlab
from scipy.interpolate import interp1d
import os
import sys
import copy
import zset
import conf as cf


def offset_curve_constant_thickness(x, z, thickness_km):
    
    """
    Offsets the curve (x,z) by thickness_km perpendicular to the local tangent.
    x, z in meters.
    thickness_km in km.
    Returns (xb, zb) for the bottom curve.
    """
    
    thickness = thickness_km ##* 1e3  # to meters
    xb, zb = [], []
    for i in range(len(x)):
        # Central difference for tangent except at ends
        if i == 0:
            dx = x[1] - x[0]
            dz = z[1] - z[0]
        elif i == len(x) - 1:
            dx = x[-1] - x[-2]
            dz = z[-1] - z[-2]
        else:
            dx = x[i+1] - x[i-1]
            dz = z[i+1] - z[i-1]

        # Normal vector
        length = np.sqrt(dx**2 + dz**2)
        nx = -dz / length
        nz = dx / length

        # Offset point
        xb.append(x[i] + nx * thickness)
        zb.append(z[i] + nz * thickness)

    return np.array(xb), np.array(zb)



# =============================================================================
#  MESHING FUNCTION (can be modified)
# =============================================================================
def meshSizeFunction(dim, tag, x, y, z, lc):
    distance_to_point = np.sqrt((x-cf.XC)**2 + (y-cf.XC)**2 + (z-cf.ZC)**2)
    
    if distance_to_point < 200:
        return cf.ES
    
    if 200 < distance_to_point < 500:
        return 2*cf.ES
    
    if 500 < distance_to_point:
        return 10*cf.ES #100


# =============================================================================
# SLAB SHAPE
# =============================================================================

if cf.CONSTANT_DIP:

    dx = 10
    x_max = (cf.TRENCH_DEPTH - cf.SLAB_MAX_DEPTH) / np.tan(np.radians(cf.SLAB_DIP_ANGLE))
    x = np.arange(0, x_max + dx, dx)
    DEPTH_SLAB = cf.TRENCH_DEPTH - x * np.tan(np.radians(cf.SLAB_DIP_ANGLE))
    
else:
    
    grid_slab = surfSlab.xyz_to_grid(cf.SLAB_XYZ_FILE)
    LON_SLAB, LAT_SLAB, DEPTH_SLAB = surfSlab.extract_slab(grid_slab, cf.LON_MIN, cf.LON_MAX,  cf.LAT_MIN,  cf.LAT_MAX, -cf.SLAB_MAX_DEPTH , 0.1)
    x = (LON_SLAB - LON_SLAB[0]) * cf.DEG_TO_KM
    DEPTH_SLAB[0] = cf.TRENCH_DEPTH
    
    
slab_top = DEPTH_SLAB
slab_bot = (DEPTH_SLAB + cf.OC_LITHOSPHERE_THICKNESS)


# =============================================================================
# INITIATE GMSH MESH
# =============================================================================

reorient = True

gmsh.initialize()
gmsh.model.add("model")


# 1. CREATE SLAB POINTS AND EDGES

slab_top_func = interp1d(x, slab_top, kind='linear', fill_value='extrapolate')
slab_down_x, slab_bot_z = offset_curve_constant_thickness(x, slab_top, cf.OC_LITHOSPHERE_THICKNESS)
slab_down_func = interp1d(slab_down_x, slab_bot_z, kind='linear', fill_value='extrapolate')

step = 10

DEPTH_TOP_SLAB = []
DEPTH_DOWN_SLAB = []

for i, xi in enumerate(np.arange(0, np.ceil(x_max), step)):
    
    d_top = slab_top_func(xi) 
    DEPTH_TOP_SLAB.append(float(d_top))
    d_down = slab_down_func(xi)
    DEPTH_DOWN_SLAB.append(float(d_down))

DEPTH_TOP_SLAB = np.array(DEPTH_TOP_SLAB)
DEPTH_DOWN_SLAB = np.array(DEPTH_DOWN_SLAB)

depth = {'litho_conti': cf.CO_LITHOSPHERE_THICKNESS, 'astheno_1': cf.H_A1, 'astheno_2': cf.H_A2}

#top slab
arg_top_slab = {'litho_conti': 0, 'astheno_1': 0, 'astheno_2': 0}
for k in depth.keys():
    arg_top = np.argmin(np.abs(depth[k] - DEPTH_TOP_SLAB))
    arg_top_slab[k] = arg_top

#down slab
arg_down_slab = {'astheno_1': 0, 'astheno_2': 0}
for k in depth.keys():
    if k != 'litho_conti':
        arg_down = np.argmin(np.abs(depth[k] - DEPTH_DOWN_SLAB))
        arg_down_slab[k] = arg_down
    

LINK_TOP_SLAB =  {'litho_conti': 0, 'astheno_1': 0, 'astheno_2': 0}
LINK_DOWN_SLAB = {'astheno_1': 0, 'astheno_2': 0}

top_slab = []
down_slab = []

line_top_slab = []
line_down_slab = []

line_slab_top = {'litho_conti': [], 'astheno_1': [], 'astheno_2': []}
line_slab_down = {'astheno_1': [], 'astheno_2': []}

last_arg_conti = 0
last_arg_ocean = 0

for i, xi in enumerate(np.arange(0, np.ceil(x_max), step)):

    if slab_top_func(xi) > cf.SLAB_MAX_DEPTH:
            
        zone = ''
        
        depth_i = slab_top_func(xi)
        
        link = 0

        for k in arg_top_slab.keys():
            
            if i == arg_top_slab[k]:
                
                link = 1
                p = gmsh.model.geo.addPoint(xi, depth[k], 0, cf.ES)
                top_slab.append(p)
                LINK_TOP_SLAB[k] = p
                zone = k
                
        if link != 1:

            p = gmsh.model.geo.addPoint(xi, slab_top_func(xi), 0, cf.ES)
            top_slab.append(p)

        
        if i !=0:
                  
            line_top_slab.append(gmsh.model.geo.addLine(top_slab[i-1], top_slab[i]))
            
            ### affcetation des groupes de lignes pour creer les elsets plus tard
            if link == 1:

                line_slab_top[zone] = line_top_slab[last_arg_conti:]
                last_arg_conti += len(line_slab_top[zone])
        
    if slab_down_func(xi) > cf.SLAB_MAX_DEPTH:   
        
        depth_i = slab_down_func(xi)
        
        link = 0
        zone = ''
        
        for k in arg_down_slab.keys():
            
            if i == arg_down_slab[k]:
                link=1
                p = gmsh.model.geo.addPoint(xi, depth[k], 0, cf.ES)
                down_slab.append(p)
                LINK_DOWN_SLAB[k] = p
                zone = k
                
        if link != 1:
    
            p = gmsh.model.geo.addPoint(xi, slab_down_func(xi), 0, cf.ES)
            down_slab.append(p)
  
        if i !=0:
            line_down_slab.append(gmsh.model.geo.addLine(down_slab[i-1], down_slab[i]))
                
            if link == 1:

                line_slab_down[zone] = line_down_slab[last_arg_ocean:]
                last_arg_ocean += len(line_slab_down[zone])
         

line_slab_top['mantle'] = line_top_slab[last_arg_conti:]
line_slab_down['mantle'] = line_down_slab[last_arg_ocean:]          

# =============================================================================
# 2. CREATION OF THE OCEANIC SIDE OF THE MESH
# =============================================================================

l1_ocean = gmsh.model.geo.addPoint(cf.MESH_EXTENT_OCEAN, -4, 0, cf.ES)
l2_ocean = gmsh.model.geo.addPoint(cf.MESH_EXTENT_OCEAN, -4 + cf.OC_LITHOSPHERE_THICKNESS , 0, cf.ES)

slab1 = gmsh.model.geo.addPoint(-50, -4, 0, cf.ES)
slab2 = gmsh.model.geo.addPoint(-50, -4 + cf.OC_LITHOSPHERE_THICKNESS, 0, cf.ES)

a1_ocean = gmsh.model.geo.addPoint(cf.MESH_EXTENT_OCEAN, cf.H_A1, 0, cf.ES)
a2_ocean = gmsh.model.geo.addPoint(cf.MESH_EXTENT_OCEAN, cf.H_A2, 0, cf.ES)
bb1 = gmsh.model.geo.addPoint(cf.MESH_EXTENT_OCEAN, cf.MESH_DEPTH, 0, cf.ES)

litho_ocean = gmsh.model.geo.addLine(l1_ocean, l2_ocean)
astheno_ocean = gmsh.model.geo.addLine(a1_ocean, a2_ocean)
litho_astheno_ocean = gmsh.model.geo.addLine(l2_ocean, a1_ocean)

litho_slab1 = gmsh.model.geo.addLine(l1_ocean, slab1)
litho_slab2 = gmsh.model.geo.addLine(l2_ocean, slab2)

closing_top_slab = gmsh.model.geo.addLine(slab1, slab2) #terminaison haute du slab (rrelier haut et bas)
closing_down_slab = gmsh.model.geo.addLine(top_slab[-1], down_slab[-1]) #terminaison basse du slab (rrelier haut et bas)

# =============================================================================
# 3. CREATION OF THE CONTINENTAL SIDE OF THE MESH
# =============================================================================
l1_conti = gmsh.model.geo.addPoint(cf.MESH_EXTENT_CONTI, 0 , 0, cf.ES)
l2_conti = gmsh.model.geo.addPoint(cf.MESH_EXTENT_CONTI, cf.CO_LITHOSPHERE_THICKNESS, 0, cf.ES)
a1_conti = gmsh.model.geo.addPoint(cf.MESH_EXTENT_CONTI, cf.H_A1, 0, cf.ES)
a2_conti = gmsh.model.geo.addPoint(cf.MESH_EXTENT_CONTI, cf.H_A2, 0, cf.ES)
bb2 = gmsh.model.geo.addPoint(cf.MESH_EXTENT_CONTI, cf.MESH_DEPTH, 0, cf.ES)

litho_conti = gmsh.model.geo.addLine(l1_conti, l2_conti)
astheno_conti = gmsh.model.geo.addLine(a1_conti, a2_conti)
litho_astheno_conti = gmsh.model.geo.addLine(l2_conti, a1_conti)

prism = gmsh.model.geo.addPoint(25, -3, 0, cf.ES)
coast = gmsh.model.geo.addPoint(cf.DIST_TRENCH_COAST, 0 , 0, cf.ES)

prism_coast = gmsh.model.geo.addLine(prism, coast)
coast_east_border = gmsh.model.geo.addLine(coast, l1_conti)

# =============================================================================
# 4. JONCTION SLAB-LITHOS CONTI/OCEAN & ASTHENO + BOUNDING BOX
# =============================================================================

slab_prism = gmsh.model.geo.addLine(prism, top_slab[0])
LOC_slab_top = gmsh.model.geo.addLine(slab1, top_slab[0])
LOC_slab_down = gmsh.model.geo.addLine(slab2, down_slab[0])

LOC_top_slab = gmsh.model.geo.addLine(l2_conti, LINK_TOP_SLAB['litho_conti'])
A1_top_slab = gmsh.model.geo.addLine(a1_conti, LINK_TOP_SLAB['astheno_1'])
A2_top_slab = gmsh.model.geo.addLine(a2_conti, LINK_TOP_SLAB['astheno_2'])

A1_down_slab = gmsh.model.geo.addLine(a1_ocean, LINK_DOWN_SLAB['astheno_1'])
A2_down_slab = gmsh.model.geo.addLine(a2_ocean, LINK_DOWN_SLAB['astheno_2'])

    
### BOTTOM
l_bb1 = gmsh.model.geo.addLine(a2_conti, bb2)
l_bb2 = gmsh.model.geo.addLine(bb1, bb2)
l_bb3 = gmsh.model.geo.addLine(a2_ocean, bb1)

# =============================================================================
# 5. DEFINITION OF GEOPHYSICAL REGIONS
# =============================================================================

gmsh.model.geo.addCurveLoop([LOC_slab_top, LOC_slab_down, closing_top_slab, closing_down_slab] + line_down_slab + line_top_slab,  1, reorient=reorient)
slab = gmsh.model.geo.addPlaneSurface([1], 1)
gmsh.model.addPhysicalGroup(2, [slab], name="SLAB")

gmsh.model.geo.addCurveLoop([litho_ocean, litho_slab1, litho_slab2, closing_top_slab],  2, reorient=reorient)
litho_oc = gmsh.model.geo.addPlaneSurface([2], 2)
gmsh.model.addPhysicalGroup(2, [litho_oc], name="LITHO-OCEAN")

gmsh.model.geo.addCurveLoop([litho_conti, LOC_top_slab, slab_prism, prism_coast, coast_east_border] + line_slab_top['litho_conti'],  3, reorient=reorient)
litho_co = gmsh.model.geo.addPlaneSurface([3], 3)
gmsh.model.addPhysicalGroup(2, [litho_co], name="LITHO-CONTI")

gmsh.model.geo.addCurveLoop([litho_astheno_conti, A1_top_slab, LOC_top_slab ] + line_slab_top['astheno_1'],  4, reorient=reorient)
astheno_1_conti = gmsh.model.geo.addPlaneSurface([4], 4)
gmsh.model.addPhysicalGroup(2, [astheno_1_conti], name="ASTHENO-1-CONTI")

gmsh.model.geo.addCurveLoop([astheno_conti, A1_top_slab, A2_top_slab ] + line_slab_top['astheno_2'],  5, reorient=reorient)
astheno_2_conti  = gmsh.model.geo.addPlaneSurface([5], 5)
gmsh.model.addPhysicalGroup(2, [astheno_2_conti], name="ASTHENO-2-CONTI")

gmsh.model.geo.addCurveLoop([litho_astheno_ocean, LOC_slab_down, A1_down_slab, litho_slab2] + line_slab_down['astheno_1'],  6, reorient=reorient)
astheno_1_ocean  = gmsh.model.geo.addPlaneSurface([6], 6)
gmsh.model.addPhysicalGroup(2, [astheno_1_ocean], name="ASTHENO-1-OCEAN")

gmsh.model.geo.addCurveLoop([astheno_ocean, A2_down_slab, A1_down_slab ] + line_slab_down['astheno_2'],  7, reorient=reorient)
astheno_2_ocean  = gmsh.model.geo.addPlaneSurface([7], 7)
gmsh.model.addPhysicalGroup(2, [astheno_2_ocean], name="ASTHENO-2-OCEAN")

gmsh.model.geo.addCurveLoop([closing_down_slab, A2_top_slab, A2_down_slab, l_bb1, l_bb2, l_bb3] + line_slab_down['mantle'] + line_slab_top['mantle'],  8, reorient=reorient)
mantle  = gmsh.model.geo.addPlaneSurface([8], 8)
gmsh.model.addPhysicalGroup(2, [mantle], name="MANTLE")

gmsh.model.mesh.setSizeCallback(meshSizeFunction)

gmsh.model.geo.synchronize()

gmsh.option.setNumber("Mesh.SaveGroupsOfNodes", 1)
gmsh.option.setNumber("Mesh.SaveGroupsOfElements", 1)

gmsh.model.mesh.generate(dim=2)

gmsh.write("{}.inp".format(cf.MESH_NAME))

# if '-nopopup' not in sys.argv:
#       gmsh.fltk.run()

# gmsh.model.occ.synchronize()

# gmsh.finalize()

# =============================================================================
# using zset
# =============================================================================


mesh = zset.Mesh("{}.inp".format(cf.MESH_NAME), format='abaqus')

mesh.transform("**remove_set \
                  *nsets_start_with LINE \
                  *elsets_start_with SURFACE LINE INTERFACE EAST WEST BOTTOM EARTH_SURFACE ALL_ELEMENT")
                                     
mesh.transform('**cleanup_bsets')

# print('\n ELSET \n')
# print(mesh.elsets)

# print('\n NSET \n')
# print(mesh.nsets)

# print('\n BSET \n')
# print(mesh.bsets)


coords = mesh.nodes_coordinates()

s = mesh.nsets['SLAB'].mask()
cl = mesh.nsets['LITHO-CONTI'].mask()
fault = np.logical_and.reduce((s,cl))
mesh.nsets.add('int_subduction', mask=fault)

mesh.transform('**sort_nset *nset_name int_subduction *criterion (y2>y1);')
mesh.transform('**bset fault_plane *use_nset int_subduction *function 1; *use_dimension 0')



mesh.transform('**nset bord_plan *use_bset fault_plane *function y > -6 - 0.1;' )

mesh.transform('**open_bset *bset fault_plane *surface bord_plan')

mesh.transform('**rename_set *bsets SIDE1 fault_plane_A_Elargi SIDE0 fault_plane_B_Elargi')

mesh.transform('**nset fault_plane_A_Elargi *use_bset fault_plane_A_Elargi *function 1. ;')
mesh.transform('**nset fault_plane_B_Elargi  *use_bset fault_plane_B_Elargi  *function 1. ;')


mesh.transform('**nset fault_plane_A *use_nset fault_plane_A_Elargi *function 1. ;')
mesh.transform('**remove_nodes_from_nset *nset_name fault_plane_A *nsets_to_remove FRONT')
mesh.transform('**nset fault_plane_B  *use_nset fault_plane_B_Elargi  *function 1. ;')
mesh.transform('**remove_nodes_from_nset *nset_name fault_plane_B  *nsets_to_remove FRONT')


tol = 1e-1

west_nodes  = np.abs(coords[:,0] + cf.MESH_EXTENT_OCEAN) < tol
mesh.nsets.add('WEST', mask=west_nodes)

east_nodes  = np.abs(coords[:,0] - cf.MESH_EXTENT_CONTI) < tol
mesh.nsets.add('EAST', mask=east_nodes)

bottom_nodes  = np.abs(coords[:,1] + cf.MESH_DEPTH) < tol
mesh.nsets.add('BOTTOM', mask=bottom_nodes)

mesh.transform("**remove_set \
                  *nsets_start_with FRONT SFRONT SLAB ASTHENO LITHO BRANCHING MANTLE \
                  *elsets_start_with SIDE side")
                           
mesh.transform('**check_orientation')
mesh.transform('**scale 1000.')
mesh.transform('**to_2d')

mesh.save('{}.geof'.format(cf.MESH_NAME))
