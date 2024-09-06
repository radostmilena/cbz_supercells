#import numpy
import numpy as np
from numpy.linalg import norm
import re

# import atomic simulation environment
import ase
from ase.io import read, write, cif
from ase.atoms import Atoms


# get number of atoms in one residue
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

    return (resnr, res_sum)



# get connectivity in CBZ
def get_connectivity_list(dist, var, varpos, symbol, atnr):
    def dist_a(a): return [(d) for d in dist[a] if d <= 1.5 and d != 0]

    def index1(a, dist0): return [(i) for i, d_i in enumerate(dist[a]) for j, d_j in enumerate(dist0) if d_i == d_j
                                  and d_i != 0 and (i in new_sym) != True]

    def index2(a, dist0): return [(i) for i, d_i in enumerate(dist[a]) for j, d_j in enumerate(dist0) if d_i == d_j
                                  and (i in new_sym) != True]

    def index3(a, dist0): return [(i) for i, d_i in enumerate(dist[a]) for j, d_j in enumerate(dist0) if d_i == d_j
                                  and (i in new_sym) != True and (i in new_sym2) != True]

    o_index = [(o) for o, s_o in enumerate(symbol) if s_o == 'O']

    atom_names = ['C1', 'C2', 'H1', 'C3', 'H2', 'C4', 'H3', 'C5', 'H4', 'C6', 'C7', 'H5', 'C8', 'H6', 'C9', 'C10',
                  'H7', 'C11', 'H8', 'C12', 'H9', 'C13', 'H10', 'C14', 'N1', 'C15', 'N2', 'H11', 'H12', 'O1']

    new_index = []
    h2o_index = []

    #loops over residue number
    for i in range(int(len(var)/atnr)*3):  

        new_sym = np.zeros(6)
        init = o_index[i]
        dist0 = dist_a(init)
        check_index = index1(init, dist0)

        #collect water molecules here
        if len(check_index) > 1:

            tmp_index = []
            tmp_index.append(init)
            tmp_index.append(check_index[0])
            tmp_index.append(check_index[1])

            h2o_index.append(np.array(tmp_index))

        #loop for CBZ
        else:

            ############### functional group ######################

            new_sym[-1] = init

            c_index = index1(init, dist0)
            new_sym[1] = c_index[0]

            c_i = c_index[0]
            dist0 = dist_a(c_i)
            n_index = index2(c_i, dist0)

            for i in range(len(n_index)):

                n_i = n_index[i]
                dist0 = dist_a(n_i)
                h_index = index2(n_i, dist0)

                if symbol[h_index[0]] == 'H' and symbol[h_index[1]] == 'H':
                    new_sym[2] = n_i
                    new_sym[3] = h_index[0]
                    new_sym[4] = h_index[1]

                else:
                    new_sym[0] = n_i

            new_sym = new_sym.astype(int)

            ############### ring system ######################

            init2 = new_sym[0]
            dist0 = dist_a(init2)
            c1_index = index1(init2, dist0)

            c1 = c1_index[0]
            new_sym2 = []
            new_sym2.append(c1)

            for i in range(24):

                dist0 = dist_a(c1)
                index = index3(c1, dist0)
                sym1 = [symbol[i] for i in index]

                for i in range(len(index)):

                    c_i = index[i]
                    dist0 = dist_a(c_i)
                    nb_index = index3(c_i, dist0)
                    sym2 = [symbol[i] for i in nb_index]

                    if any("H" in s for s in sym2):

                        c1 = c_i
                        h_i = [(j) for i, j in enumerate(
                            nb_index) if symbol[j] == 'H'][0]
                        new_sym2.append(c1)
                        new_sym2.append(h_i)

                    else:

                        if len(index) == 1:

                            c1 = index[0]
                            new_sym2.append(c1)

            new_index.append(np.concatenate((new_sym2, new_sym)))

    for i in range(len(h2o_index)):
        new_index.append(h2o_index[i])

    h2o_counter = len(h2o_index)

    return (new_index, atom_names, h2o_counter)



# create arrays for PDB writing
def create_arrays(new_index, varpos, symbol, atnr, atom_names, h2o_counter):

    new_pos = []
    symb = []
    atns = []
    resnrs = []
    cbz_length = int(len(new_index) - h2o_counter)

    for i in range(cbz_length):
        for j in range(atnr-6):
            resnrs.append(i+1)
            atns.append(atom_names[j])
            new_pos.append(varpos[new_index[i][j]])
            symb.append(symbol[new_index[i][j]])

    water_names = ["Ow", "Hw1", "Hw2"]
    for i in range(cbz_length, len(new_index), 1):
        for j in range(3):
            resnrs.append(i+1)
            atns.append(water_names[j])
            new_pos.append(varpos[new_index[i][j]])
            symb.append(symbol[new_index[i][j]])

    return (new_pos, symb, resnrs, atns)
