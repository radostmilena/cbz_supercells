#!/usr/bin/env python
# coding: utf-8

# ### This notebook creates a supercell from a unit cell described by a cif-file
# 
# In this example the input is a cif-file containing the structure of hexagonal ice. A supercell is constructed by replicating the unit cell three times along each unit cell vector. The unit cell and supercell are then written to PDB files (unitcell.pdb and supercell.pdb).

# In[1]:


#import numpy
import numpy as np
from numpy.linalg import norm
import re

#import atomic simulation environment
import ase
from ase.io import read, write, cif
from ase.atoms import Atoms

#modified pdb writer
from read_pdb import write_proteindatabank

#custom functions
from cbz_functions import get_nr_of_atoms, get_connectivity_list 


# In[2]:


#read crystal structure from cif file
structure = read('form1.cif')

pos = structure.get_positions() #positions from cif file

sym = structure.get_chemical_symbols() #element symbols from cif file

#unit cell vectors
a = structure.cell[0]
b = structure.cell[1]
c = structure.cell[2]

chemform, atnr = get_nr_of_atoms(structure)
print(chemform, atnr)


# In[3]:


#repeat unit cell
[na, nb, nc] = [6, 2, 2]
var = structure.repeat(rep=([na, nb, nc]))
varpos = var.get_positions()
symbol = var.get_chemical_symbols()
print(len(varpos))


# In[4]:


#get distance matrix
dist = var.get_all_distances(mic=True)
print(len(dist))


# In[5]:


#print new coordinates

new_index, atom_names = get_connectivity_list(dist, var, varpos, symbol, atnr)

#with open('testfile.txt', 'w') as fh:
new_pos = []
symb = []
atns = []
resnrs = []

for i in range(len(new_index)):
    for j in range(atnr):
        resnrs.append(i+1)
        atns.append(atom_names[j])
        new_pos.append(varpos[new_index[i][j]]) 
        symb.append(symbol[new_index[i][j]])
        
atms = Atoms('O{}'.format(str(len(new_pos))),positions=new_pos) #create atoms object

atms.set_chemical_symbols(symb) #set new chemical symbols


# In[6]:


print(len(atms.get_positions()))


# In[7]:


#write unit cell and super cell to PDB files
#write_proteindatabank('unitcell.pdb',structure,write_arrays=True)
write_proteindatabank('supercell.pdb', atms, atns, resnrs, write_arrays=True)


# In[8]:


#fix periodic boundary conditions in PDB file for supercell

x = '%.3f' % (na*norm(a)) #box x dimension
y = '%.3f' % (nb*norm(b)) #box y dimension
z = '%.3f' % (nc*norm(c)) #box z dimension

alpha = '%.2f' % ( np.rad2deg(np.arccos(np.dot(b,c)/(norm(b)*norm(c)))) ) #angle between b and c vector
beta = '%.2f' % ( np.rad2deg(np.arccos(np.dot(a,c)/(norm(a)*norm(c)))) ) #angle between a and c vector
gamma = '%.2f' % ( np.rad2deg(np.arccos(np.dot(a,b)/(norm(a)*norm(b)))) ) #angle between a and b vector

with open('supercell.pdb', 'r') as original: data = original.read()
with open('supercell.pdb', 'w') as modified: modified.write("CRYST1   {}   {}   {}  {}  {} {} P 1\n".format(x,y,z,alpha,beta,gamma) + data)

