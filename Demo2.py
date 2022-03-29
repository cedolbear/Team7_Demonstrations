#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 12:23:47 2022

@author: caradolbear
"""



from __future__ import division, print_function
import numpy as np
import solidspy.postprocesor as pos
import solidspy.assemutil as ass
import solidspy.solutil as sol
import solidspy.uelutil as uel



nodes = np.array([
  #nodal identifier x-coord y-coord x-constr y-constr   
    [0, 10*np.sqrt(3)/2, 0.0, 0, 0],
    [1, 0.0, 5.0, -1, 0], # change 
    [2, 0.0, 0.0, -1, -1]])

mats = np.array([
    [1e6, 0.01],
    [1e6, 0.01],
    [1e6, 0.01]])

eles = np.array([
    [0, 6, 0, 2, 0],
    [1, 6, 1, 0, 1],
    [2, 6, 2, 1, 2]])
# 6 is element type for truss

loads = np.array([
    [0, 0.0, -200.0],
    [1, 0.0, 0.0]])


for cont in range(3):
    Kloc = uel.ueltruss2D(nodes[eles[cont, 3:], 1:3], *mats[cont, :])
    print("Element {}\n".format(cont), np.round(Kloc))

# System assembly
DME , IBC , neq = ass.DME(nodes, eles)
stiff = ass.assembler(eles, mats, nodes, neq, DME, sparse=False)
load_vec = ass.loadasem(loads, IBC, neq) # right-hand-side
print(np.round(stiff))

# print(load_vec)

# System solution
disp = sol.static_sol(stiff, load_vec)
disp_complete = pos.complete_disp(IBC, nodes, disp)
print(disp_complete)

stress = pos.stress_truss(nodes, eles, mats, disp_complete)

print(stress)

nodes_disp = nodes.copy()
nodes_disp[:,1] = nodes[:,1]+disp_complete[:,0]
nodes_disp[:,2] = nodes[:,2]+disp_complete[:,1]

pos.plot_truss(nodes, eles, mats, title = 'Original Truss')
pos.plot_truss(nodes_disp,eles,mats,stresses = stress, title = 'Deformed Truss')
