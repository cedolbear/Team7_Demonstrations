
"""
Created on Wed Mar 30 09:38:39 2022

@author: inrae
"""

import numpy as np
import solidspy.postprocesor as pos
import solidspy.assemutil as assm
import solidspy.solutil as sol
import solidspy.uelutil as uel


# nodes information
nodes = np.array([
    [0, 10*np.sqrt(3)/2, 0.0, 0, 0],
    [1, 0.0, 5.0, -1, 0],    #Last two elements indicate boundary conditin for each node 
    [2, 0.0, 0.0, -1, -1]])   #0: Unconstrained, -1:Constraained

# material profiles
mats = np.array([
    [1e6, 0.01],    #Elastic Modulus and Cross-sectional area
    [1e6, 0.01], 
    [1e6, 0.01]])              

# element information
eles = np.array([
    [0, 6, 0, 2, 0],     # Second element indicated element type
    [1, 6, 1, 0, 1],     # 6 is element type for truss
    [2, 6, 2, 1, 2]])
    

# load information
loads = np.array([
    [0, 0.0, -200.0],
    [1, 0.0, 0.0]])

# generate stiffness of each element
for k in range(3):
    Kloc = uel.ueltruss2D(nodes[eles[k, 3:], 1:3], *mats[k, :])
    print("\nElement {} Stiffness:\n".format(k), np.round(Kloc))

# generate global stiffness of whole system
DME , IBC , neq = assm.DME(nodes, eles) 
# Count active equations, create boundary conditions array and the assembly operator
# DME = Assembly operator , IBC = boundary condition array, neq = number of active equation
stiff = assm.assembler(eles, mats, nodes, neq, DME, sparse=False)
load_vec = assm.loadasem(loads, IBC, neq)     #Load vector in matrix form
print('\nStiffness Matrix of Structure:\n{}'.format(np.round(stiff)))

# system solution
disp = sol.static_sol(stiff, load_vec)     #Static solver that solves linear equation, Ax = B
disp_complete = pos.complete_disp(IBC, nodes, disp)     #Formats solution in matrix form 
print('\nDisplacement Vectors:\n{}'.format(disp_complete))

stress = pos.stress_truss(nodes, eles, mats, disp_complete)   #Calculate stress for each element
print('\nStresses:\n{}'.format(stress))

# generate displaced nodes
nodes_disp = nodes.copy()
nodes_disp[:,1] = nodes[:,1]+disp_complete[:,0]     #Adding calculated displacement to node 0
nodes_disp[:,2] = nodes[:,2]+disp_complete[:,1]     #Adding calculated displacement to node 1

# plots to compare original and deformed truss
pos.plot_truss(nodes, eles, mats, title = 'Original Truss')   
pos.plot_truss(nodes_disp,eles,mats,stresses = stress, title = 'Deformed Truss')
