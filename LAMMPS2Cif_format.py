## This python script is to convert LAMMPS production dump file to CIF file to generate 21 CIF files that can be used to model flexibility using GCMC simulations

import numpy as np
import os
import sys
import subprocess
from subprocess import call
import ase
from ase.io import read, write
from ase import *

starts = 0
stops =0
mov_no = 1

def line_prepend(filename, line):
    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0,0)
        f.write(line.rstrip('\r\n') + '\n' + content)
        
prependThis='3648\n \n'

                

x = np.array([],dtype=float)
y = np.array([],dtype = float)
z = np.array([],dtype = float)
symbol = np.array([],dtype = float)
charge = np.array([],dtype = float)
sym_char = np.array([])
big_charge = np.array([],dtype = float)

with open('dump.production.lammpstrj','r') as f:
    f1 = f.readlines()
    for line in range(len(f1)):
#         if "TIMESTEP" in f1[line]:
#             starts += 1
        if "id" in f1[line]:
            x = np.array([],dtype=float)
            y = np.array([],dtype = float)
            z = np.array([],dtype = float)
            symbol = np.array([],dtype = float)
            charge = np.array([],dtype = float)
            sym_char = np.array([])
            starts = starts + 9
            for n in np.r_[starts:starts + 4544]:  ## 4544 because 7 acetone moelcules are in the lammps dump that I am reading
                i = f1[n].split()
                if i[1]=='1':
                    symbol = np.append(symbol,i[2])
                    x = np.append(x,float(i[7]))
                    y = np.append(y,float(i[8]))
                    z = np.append(z,float(i[9]))
                    charge = np.append(charge,i[6]).astype(float)
            big_charge = np.append(big_charge,charge)[:,np.newaxis]
            for k in symbol:
                if k == '1':
                    k = 'C'
                if k == '2':
                    k = 'C'
                if k == '3':
                    k = 'C'
                if k == '4':
                    k = 'H'
                if k == '5':  ###mu3H
                    k = 'Cu'
                if k == '6':
                    k = 'O'
                if k == '7':  ###mu3O
                    k = 'Ne'
                if k == '8':
                    k = 'O'
                if k == '9':
                    k = 'Zr'
                sym_char = np.append(sym_char,k)
            
            os.system('mkdir Movies-xyz')
            np.savetxt(f'Movies-xyz/fraccoord_MOV-{mov_no}.xyz',np.c_[sym_char,x,y,z], fmt = '%s')
            line_prepend(f'Movies-xyz/fraccoord_MOV-{mov_no}.xyz', prependThis)
            mov_no += 1
            starts += 4544  # Finished the movie


def line_prepend(filename, line):
    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0,0)
        f.write(line.rstrip('\r\n') + '\n' + content)
        
prependThis1 = '# CIF file generated by openbabel 2.4.1, see http://openbabel.sf.net \ndata_I \n_chemical_name_common \'image0\' \n_cell_length_a 41.9568 \n_cell_length_b 41.9568 \n_cell_length_c 41.9568 \n_cell_angle_alpha 90 \n_cell_angle_beta 90 \n_cell_angle_gamma 90 \n_space_group_name_H-M_alt \'P 1\' \n_space_group_name_Hall \'P 1\' \nloop_ \n  \t _symmetry_equiv_pos_as_xyz \n  \t x,y,z \nloop_ \n   \t_atom_site_label \n   \t _atom_site_type_symbol \n   \t _atom_site_fract_x \n   \t _atom_site_fract_y \n   \t _atom_site_fract_z \n   \t _atom_site_charge  \n   \t _atom_site_occupancy \n'

xx = 0
num = 1
bb_num = 0
for i in range(21):
    xx = 0
    with open(f'MOV_{num}.cif', 'r+') as fd:
        with open(f'FinalMovies/movie-nomu3OH-{num}.cif','w+') as f5:
            contents = fd.readlines()

            for line in contents:
                line = line.split()
                if len(line) == 7 and xx<3648 and line[0]!='_chemical_formula_sum':
                    newline = '\t{0} \t{1} \t{2} \t{3} \t{4} \t{5} \t{6} \n'.format(line[1],line[0],line[3],line[4],line[5],bb[bb_num][xx],line[6])
                    f5.write(newline)
                    xx+=1
        line_prepend(f'FinalMovies/movie-nomu3OH-{num}.cif', prependThis1)
        num+=1
        bb_num+=1
