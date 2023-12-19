# Residence time vizualization code.
# Calculates center of mass for each adsorbate molecule, storing the values.
#Then counts the number of adsorbate molecules in each "bin" forming a 3D distribution of the adsorbate locations.
# COM function loads 1 snapshot (1 timestep) from trajectories at a time, storing COM data in a list, before moving on to
# the next snapshot. This avoids proccessing issues which would occur from the same molecule appearing in different snapshots.
# Thoughts - Maybe read in 1 snapshot at a time (read in total number of atoms at a time). This will prevent
# loading too much at once, and prevent confusion when grouping by molecule number (as equal molecule numbers will
# appear once for each snapshot)
import numpy as np
import math as m
import pandas as pd

# Center of Mass function
# Input:
# Input_file - Name of desired input. Should be a trajectory file
# Mass_dict - Dictionary containing element masses in the form {'Element Abbrviation': mass} (i.e. {'C': 12.01, 'O': 15.999})
# Snapshots to read - Number of snapshots to read from file. 1 would just use atom locations in the first set of trajectories
# num_atoms - Number of atoms in the system (also the size of a single snapshot in lines without the header)
# May update num_atoms to read num_atoms from the file, as it appears in the header of lammps trajectory files
# atoms_per_adsorbate - Number of atoms in one adsorbate molecule as it appears in the trajectory file(for example IPA = 5: 3 carbons, 1 hydrogen, 1 oxygen)
# Ouput:
# COMs - list of all center of mass values, each is an (x,y,z) coordinate.
# Err - number of molecules where number of atoms != atoms_per_adsorbate. This should be zero, otherwise check that the read in line is reading
# the correct number of lines in at a time.
def calc_center_of_masses(input_file, mass_dict, snapshots_to_read, num_atoms, atoms_per_adsorbate, MOF_atoms):
    # stores all center of mass values
    COMs = []
    # err counts the number of molecule groupings that do not match the total number of atoms_per_adsorbate.
    err = 0
    # Iterate through every desired snapshot
    for snp in range(snapshots_to_read):
        print('Snapshot Number:', snp + 1)
        # Reads in 1 snapshot at a time
        if snp == 0:
            df = pd.read_csv(input_file, skiprows = 9, nrows = num_atoms, names = ['mol', 'type', 'x', 'y', 'z'], usecols=[1,2,3,4,5], delim_whitespace=True)
        else:
            skips = num_atoms*snp + 9*(snp+1)
            df = pd.read_csv(input_file, skiprows = skips, nrows = num_atoms, error_bad_lines = False, names = ['mol', 'type', 'x', 'y', 'z'], usecols=[1,2,3,4,5], delim_whitespace=True)
        print(df.head())
        # Removes all atoms that belong to the MOF (in UiO-66 potential framework, mol = 444)
        # If Framework has more than one mol number, passing a list of mol values will work
        df = df[~df['mol'].isin(MOF_atoms)].sort_values('mol')
        print(df.head(10))
        print(len(df))
        # Sorting of dataframe so entries are grouped by molecule number (each is 1 adsorbate molecule) and sorted
        # by element symbol within those groups.
        df_grouped = df.groupby(df['mol'])
        df_grouped = df_grouped.apply(lambda x: x.sort_values(['type'], ascending = False))
        df_grouped = df_grouped.reset_index(drop = True)
        df_grouped = df_grouped.groupby('mol')
	# print(df)
        # Group counter

        # Iterate over each group (molecule), calculating Center of Mass for each
        n = 0
        for mol, group in df_grouped:
            if n == 0:
                # Mass sum, the denominator of the COM equation, calculated during the first molecule run through.
                # Becuase the atoms are grouped by mol number and sorted by element, this is identical for all molecules
		# and only needs to be computed once
                mass_sum = 0
                for j in range(atoms_per_adsorbate):
                    temp_atm = df_grouped.get_group(mol).iloc[j]
                    #print(temp_atm)
                    mass_sum += masses[temp_atm['type']]
                    print(mass_sum)

            n += 1
            # Check for groups that have atoms which do not add to the number of atoms in the adsorbate, ideally = 0
            if len(df_grouped.get_group(mol)) != atoms_per_adsorbate:
                err += 1
                continue

	    # List to hold each atom_mass * atm_position value for an individual molecule
            molecule = []
            #Iterate over each atom in the adsorbate, calculating the position * mass. Add to the molecule list
            for k in range(atoms_per_adsorbate):
                atm = df_grouped.get_group(mol).iloc[k]
                #print('\n',atm['element'])
                atom_calc = masses[atm['type']] * atm[['x', 'y', 'z']]
                molecule.append(atom_calc)

            # Calculate the center of mass using the values in the molecule list and the overall mass
            # COM calculation: SUM(Mass of atom n * (xn,yn,zn)) for n = 1 through n = num_adsorbate / mass_sum

        COM_SUM = 0
        for atom in range(atoms_per_adsorbate):
            COM_SUM += molecule[atom]
            COM = COM_SUM / mass_sum
            #print(mol)
            #print(COM)
            # Add center of mass to list of COMs
            COMs.append(list(COM))
        #print(n)
    return COMs, err

# Mass dictionary parameter. Store each unique mass value with its corresponding element, to be called in the COM fucnction
# These are parameters to change based on each system
masses = {10: 1.007, 11: 15.999, 12: 12.01, 13: 15.03, 14: 15.03}
snapshots_to_read = 5
loading = 1
atoms_per_adsorbate = 5
num_atoms = 4768 #3808

MOF_atoms = [1]
input_file = '/bgfs/kjohnson/pbs13/FlexibleMOFs-GCMC/NVT-MD-all/7loadedIPA-UiO-66/dump.production.lammpstrj'
# Call to center of mass calculator
coms, err = calc_center_of_masses(input_file, masses, snapshots_to_read, num_atoms, atoms_per_adsorbate, MOF_atoms)
# Center of mass to 3D histogram data file function.
# Input - center_of_masses - array of center of masses for adsorbate molecules
# cell_dimensions - dimensions of the overall system, given in dictionary format {'x': xdim, 'y': ydim, 'z': zdim}
# box_size - size for each box/bin of the histogram system. Given in the same format as cell_dimensions
# Output - None. Saves a data file of the 3d distribution
print("COM",coms)
print(np.array(coms).shape)
print("length",len(coms))
def COM_histogram(center_of_masses, cell_dimensions, box_size):
    # Calculate number of bins in each direction
    x_bins = int(m.ceil(cell_dimensions['x'] / box_size['x']))
    y_bins = int(m.ceil(cell_dimensions['y'] / box_size['y']))
    z_bins = int(m.ceil(cell_dimensions['z'] / box_size['z']))
    #print(x_bins)
    # Create 3D array of adsorbate distrubution
    hist_array = np.zeros((x_bins, y_bins, z_bins))

    # For each adsorbate, calculate its bin (location in distribution)
    # Increase the value in the corresponding bin by 1
    for com in center_of_masses:
        ind_xbin = int(m.floor(com[0]/box_size['x']))
        ind_ybin = int(m.floor(com[1]/box_size['y']))
        ind_zbin = int(m.floor(com[2]/box_size['z']))
        hist_array[ind_xbin][ind_ybin][ind_zbin] += 1

    #print(hist_array[1])
    print("hist array:",hist_array)
    final = np.array(hist_array)
    # Save distrivution data in a redable format.
    with open('adsorbate_distribution.dat', 'w') as f:
	    f.write('# Array Shape: {0}\n# To load in python\n# dist = numpy.loadtxt(\'adsorbate_distribution.dat\')\n# dist.reshape(xbins, ybins, zbins)\n# New Slice\n'.format(final.shape))

	    for dat_slice in final:
		    np.savetxt(f, dat_slice, fmt = '%-3d')
		    f.write('# New Slice\n')

# Set cell dimensions (assuming 0,0,0 start)
cell_dimensions = {'x': 41.4008, 'y': 41.4008, 'z': 41.4008}
# Select size of bins (cell_dim/del_x should equal an integer)
box_size = {'x': 1.03502, 'y': 1.03502, 'z': 1.03502}

# Call to center of mass histogram function
COM_histogram(coms, cell_dimensions, box_size)
# Save center of mass data. Store in a trajectory style file to allow for easy vizualization in VMD/Ovito
com_df = pd.DataFrame(np.array(coms))
com_df.columns = ['x', 'y', 'z']
#com_df.to_csv('COM.csv')
com_df['id'] = list(range(len(com_df['x'])))
com_df['mol'] = [100 for x in range(len(com_df['x']))]
com_df['type'] = [1 for x in range(len(com_df['x']))]

with open('dump.COMIPALowLoading.lammpstrj', 'w') as f:
    f.write('ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n80000\nITEM: BOX BOUNDS pp pp pp\n0.000e+00 4.14008e+01\n0.000e+00 4.14008e+01\n0.000e+00 4.14008e+01\nITEM: ATOMS ')
com_df.to_csv('dump.COMIPALowLoading.lammpstrj', sep = ' ', columns = ['id', 'mol', 'type', 'x', 'y', 'z'], index = False, mode = 'a')

print('Number of Molecules in Distribution:', len(coms))
print('Number of Molecules with less than the number of atoms in adsorbate molecule and thus not in the distribution:',err)
