#%%
from __future__ import unicode_literals
import os, sys
import numpy as np
import pickle as pk
import meshio as msh
import time
import compaction as cpt

import warnings
warnings.filterwarnings("ignore")

########################################################################
# Auxiliary functions
########################################################################
#%%
def getvtkData(filename):
    """
    This function uploads coordinates, and pressure field 
    """
    data = msh.read(filename)
    
    # Here the nodes coordinates
    # coordinates 3 lines, Nnodes columns
    points = data.points.T
    
    nodes = np.zeros((3,len(points[0][:])))
    for k in range(len(points[0][:])):
        nodes[0][k] = points[1][k]
        nodes[1][k] = points[0][k]
        nodes[2][k] = points[2][k]
    
    print('Nodes shape',nodes.shape)
    
    #x-nodes: nodes[1][:]
    #y-nodes: nodes[0][:]
    #z-nodes: nodes[2][:]
    
    # Upload the connectivity of the mesh
    # prisms Nelem lines, 8 columns (nodes per elem) =data.cells_dict['hexahedron'].shape[1]
    print(data.cells)
    connectivity = data.cells_dict['hexahedron']
    
    print(connectivity[0])
    print(connectivity.shape)
    
    elem_coords = np.zeros((len(connectivity),data.cells_dict['hexahedron'].shape[1]))
    
    for k in range(len(connectivity)):
        prism = connectivity[k]
        xloc , yloc , zloc = [] , [] , []
        for i in prism:
            yloc.append(nodes[0][i])
            xloc.append(nodes[1][i])
            zloc.append(nodes[2][i])
        elem_coords[k][0] = np.amin(yloc)
        elem_coords[k][1] = np.amax(yloc)
        
        elem_coords[k][2] = np.amin(xloc)
        elem_coords[k][3] = np.amax(xloc)

        elem_coords[k][4] = np.amax(zloc)
        elem_coords[k][5] = np.amin(zloc)
    
    print('elems',elem_coords.shape)

    # Upload pressure but change the nan for -1
    # pressure Nelem lines
    ptest = data.cell_data['pressure']
    
    p = np.zeros(len(ptest[0]))
    for k in range(len(p)):
        if(np.isnan(ptest[0][k])):
            p[k] = -1
        else:
            p[k] = ptest[0][k]

    return nodes, elem_coords, p
#%%
if __name__ == "__main__":
    
    time_0 = time.process_time() # Here start count time
    
    # Load the data from the vtk file
    
    # Make a loop for time (pressure = p - p_init)
    coordinates, prisms, pressure = getvtkData('data/vtk_data/data_ts0.vtk')
    
    # Define elastic constants
    poisson = 0.25
    young   = 3300.
    
    # Now we solve the reservoir's displacement
    ux = cpt.displacement_x_component(coordinates, prisms, pressure, poisson, young)
    uy = cpt.displacement_y_component(coordinates, prisms, pressure, poisson, young)
    uz = cpt.displacement_z_component(coordinates, prisms, pressure, poisson, young)
    
    # Now we solve the reservoir's stress state
    sx = cpt.stress_x_component(coordinates, prisms, pressure, poisson, young)
    sy = cpt.stress_y_component(coordinates, prisms, pressure, poisson, young)
    sz = cpt.stress_z_component(coordinates, prisms, pressure, poisson, young)
    
#    # Save the results
#    fields_res = dict()
#    fields_res['x'] = coordinates[1][:]
#    fields_res['y'] = coordinates[0][:]
#    fields_res['z'] = coordinates[2][:]
#    # Displacement field
#    fields_res['u_x'] = ux
#    fields_res['u_y'] = uy
#    fields_res['u_z'] = uz
#    # Stress field
#    fields_res['u_x'] = sx
#    fields_res['u_y'] = sy
#    fields_res['u_z'] = sz
#    
#    # Save the fields
#    file_name = 'data_ts0.pickle'
#    with open(file_name, 'wb') as f:
#        pk.dump(fields_d, f)
    
    time_1 = time.process_time() # Here end counting time
    
    print("Elapsed time to plot: %0.3f mins." % ((time_1-time_0)/60.0))

    

# %%
