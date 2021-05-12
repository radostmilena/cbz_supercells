#import numpy
import numpy as np
from numpy.linalg import norm
import re

#import atomic simulation environment
import ase
from ase.io import read, write, cif
from ase.atoms import Atoms


#get number of atoms in one residue
def get_nr_of_atoms(structure):
    resnr = structure.get_chemical_formula(mode='hill', empirical=True)
    res = re.findall(r'\d+', resnr)
    all_str = re.findall(r'\w', str(res))
    res_array = list(map(int, res))

    if (len(resnr)-len(all_str)) > len(res_array):
        rest = (len(resnr)-len(all_str))-len(res_array)
        res_sum = sum(res_array)+rest
    else:
        res_sum = sum(res_array)

    return(resnr, res_sum)

#get connectivity in PBZ
def get_connectivity_list(dist, var, varpos, symbol, atnr):
    dist_a = lambda a : [(d) for d in dist[a] if d <= 1.6 and d!= 0]

    index = lambda a, dist0 : [(i) for i, d_i in enumerate(dist[a]) for j, d_j in enumerate(dist0) if d_i == d_j
                            and d_i != 0 and (i in new_sym2) != True]

    atom_names = ['C1', 'H1', 'C2', 'H2', 'C3', 'H3', 'C4', 'H4', 'C5', 'H5', 'C6', 'N1', 'C7', 'O1', 'C8', 'C9', 
              'C10', 'C11', 'C12', 'H6', 'H7', 'H8', 'H9', 'H10', 'H11', 'H12', 'H13', 'H14', 'H15', 'C13', 
              'O2', 'N2', 'C14', 'C15', 'H16', 'C16', 'H17', 'C17', 'H18', 'C18', 'H19', 'C19', 'H20']

    new_index = []
    new_sym2 = []

    for i in range(int(len(var)/atnr)):

        new_sym = np.zeros(43).astype(int)
        o_index = [(o) for o, s_o in enumerate(symbol) if s_o == 'O']
    
        #--------O1--------
        init_list = [(o) for i, o in enumerate(o_index) if (o in new_sym2) != True]
        new_sym[13] = init_list[0]
        init = init_list[0]

        #---------C7----------
        dist0 = dist_a(init)
        c_index = index(init, dist0)
        new_sym[12] = c_index[0]
        new_sym2.append(c_index[0])
        new_sym2.append(init)

        #------------C8--------------------
        dist0 = dist_a(c_index[0])
        c_index = index(c_index[0], dist0)
        sym2 = [symbol[i] for i in c_index]
        c_index2 = [(i) for i, sym in enumerate(sym2) if sym == 'C']
        c8 = c_index[c_index2[0]]
        new_sym[14] = c8
        new_sym2.append(c8)

        #---------------H15----------------
        dist0 = dist_a(c_index[c_index2[0]])
        c_index = index(c_index[c_index2[0]], dist0)
        sym2 = [symbol[i] for i in c_index]

        h_index = [(i) for i, sym in enumerate(sym2) if sym == 'H']
        new_sym[28] = c_index[h_index[0]]
        new_sym2.append(c_index[h_index[0]])

        #---------------C9-----------------------------
        c_index2 = [(i) for i, sym in enumerate(sym2) if sym == 'C']
        for c in c_index2:
    
            dist0 = dist_a(c_index[c])
            c_index3 = index(c_index[c], dist0)
            sym2 = [symbol[i] for i in c_index3]
    
            if any("H" in s for s in sym2):
        
                new_sym[15] = c_index[c]
                c_i = c_index[c]
                new_sym2.append(c_i)
    
            else:
        
                c13 = c_index[c]

        #----------C10, C11, C12-------------------        
        for i in range(16, 19, 1):
    
            dist0 = dist_a(c_i)
            c_index = index(c_i, dist0)
            sym2 = [symbol[i] for i in c_index]
            c_index2 = [(i) for i, sym in enumerate(sym2) if sym == 'C']
            c_i = c_index[c_index2[0]]
            new_sym[i] = c_i
            new_sym2.append(c_i)

        #----------- H6, H7, H8 ---------------
        c_ter = new_sym2[-1]
        c_next = [new_sym2[-2], new_sym2[-3], new_sym2[-4]]
        dist0 = dist_a(c_ter)
        c_index = index(c_ter, dist0)
        sym2 = [symbol[i] for i in c_index]
        h_index = [(i) for i, sym in enumerate(sym2) if sym == 'H'] 

        for h in h_index:
    
            new_sym[h+19] = c_index[h]
            new_sym2.append(c_index[h])
    
        #----------H9, H10, H11, H12, H13, H14---------------
        i = 22
        for c in c_next:

            dist0 = dist_a(c)
            c_index = index(c, dist0)
            sym2 = [symbol[i] for i in c_index]
            h_index = [(i) for i, sym in enumerate(sym2) if sym == 'H'] 
    
            for h in h_index:
        
                new_sym[i+h] = c_index[h]
                new_sym2.append(c_index[h])
        
            i+=2

        #---------------C13-------------------
        dist0 = dist_a(c8)
        c_index = index(c8, dist0)
        new_sym[29] = c_index[0]
        new_sym2.append(c_index[0])

        #---------------O2, N2----------------
        dist0 = dist_a(c_index[0])
        no_index = index(c_index[0], dist0)
        sym2 = [symbol[i] for i in no_index]
        n_index = [(i) for i, sym in enumerate(sym2) if sym == 'N'] 
        o_index = [(i) for i, sym in enumerate(sym2) if sym == 'O'] 
        new_sym[30] = no_index[o_index[0]]
        n2 = no_index[n_index[0]]
        new_sym[31] = n2
        new_sym2.append(no_index[o_index[0]])
        new_sym2.append(n2)

        #---------------second ring (C14-C19)------------------------
        dist0 = dist_a(n2)
        c_index = index(n2, dist0)
        sym2 = [symbol[i] for i in c_index]
        c_i = [(i) for i, sym in enumerate(sym2) if sym == 'C']
        n_i = [(i) for i, sym in enumerate(sym2) if sym == 'N']
        c14 = c_index[c_i[0]]
        n1 = c_index[n_i[0]]
        new_sym[32] = c14
        new_sym2.append(c14)

        dist0 = dist_a(c14)
        c_index = index(c14, dist0)
        new_sym[33] = c_index[0]
        new_sym2.append(c_index[0])
        dist0 = dist_a(c_index[0])
        nb_index = index(c_index[0], dist0)
        h_i = [(j) for i, j in enumerate(nb_index) if symbol[j] == 'H'][0]
        new_sym[34] = h_i
        new_sym2.append(h_i)
        c14 = c_index[0]

        j = 0
        for i in range(8):

            dist0 = dist_a(c14)
            c_index = index(c14, dist0)
    
            for i in range(len(c_index)):

                c_i = c_index[i]
                dist0 = dist_a(c_i)
                nb_index = index(c_i, dist0)
                sym2 = [symbol[i] for i in nb_index]

                if any("H" in s for s in sym2):

                    c14 = c_i
                    h_i = [(j) for i, j in enumerate(nb_index) if symbol[j] == 'H'][0]
                    new_sym[35+j] = c_i
                    new_sym2.append(c14)
                    new_sym[35+j+1] = h_i
                    new_sym2.append(h_i)
                    j+=2
            
                else:

                    if len(c_index) == 1:
   
                        break
                
        #----------------first ring (C1-C6)----------------
        new_sym[11] = n1
        new_sym2.append(n1)
        dist0 = dist_a(n1)
        c_index = index(n1, dist0)
        sym2 = [symbol[i] for i in c_index]
        c_i = [(i) for i, sym in enumerate(sym2) if sym == 'C']
        c6 = c_index[c_i[0]]
        new_sym[10] = c6
        new_sym2.append(c6)

        dist0 = dist_a(c6)
        c_index = index(c6, dist0)
        new_sym[8] = c_index[0]
        new_sym2.append(c_index[0])
        dist0 = dist_a(c_index[0])
        nb_index = index(c_index[0], dist0)
        h_i = [(j) for i, j in enumerate(nb_index) if symbol[j] == 'H'][0]
        new_sym[9] = h_i
        new_sym2.append(h_i)
        c6 = c_index[0]

        j = 2
        for i in range(8):

            dist0 = dist_a(c6)
            c_index = index(c6, dist0)
    
            for i in range(len(c_index)):

                c_i = c_index[i]
                dist0 = dist_a(c_i)
                nb_index = index(c_i, dist0)
                sym2 = [symbol[i] for i in nb_index]

                if any("H" in s for s in sym2):

                    c6 = c_i
                    h_i = [(j) for i, j in enumerate(nb_index) if symbol[j] == 'H'][0]
                    new_sym[8-j] = c_i
                    new_sym2.append(c6)
                    new_sym[(8-j)+1] = h_i
                    new_sym2.append(h_i)
                    j+=2
            
                else:

                    if len(c_index) == 1:

                        break
        
        #------------append new residue--------------
    
        for i in range(len(new_sym)):
        
            new_index.append(new_sym[i])
               
    return(new_index, atom_names)

#create arrays for PDB writing
def create_arrays(new_index, varpos, symbol, atnr, atom_names):

    new_pos = []
    symb = []
    atns = []
    resnrs = []

    for i in range(len(new_index)):

        new_pos.append(varpos[new_index[i]])
        symb.append(symbol[new_index[i]])
        resnrs.append(int(i/atnr)+1)
        atns.append(atom_names[(i%atnr)])

    return(new_pos, symb, resnrs, atns)
