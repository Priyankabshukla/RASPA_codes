import numpy as np
import pandas as pd
import math as m

## Function for center of mass formula
def calc_com(pos):
    return sum(pos*mass)/sum(mass)
    
## Mass of atoms in IPA molecule: CH:13.018, CH3: 15.03, H:1.007, O:15.999
# masses = {10: 1.007, 11: 15.999, 12: 12.011, 13: 15.03, 14: 15.03}
mass = np.array([13.018,15.999,1.007,15.03,15.03])

def IPA_COM(filename1,name):
    starts, stops = [], []
    with open(f'{filename1}') as f:
        for i, line in enumerate(f):
            if line.startswith('MODEL'):
                starts += [i]
            if line.startswith('ENDMDL'):
                stops += [i]


    dataframes = [pd.read_csv(f'{filename1}',
                              delimiter='\s+',
                             skiprows=start + 2,
                             nrows=stop - start - 2,
                              index_col=False,
                             names=['ATOM','ANO','SYM1','Mol','X','Y','Z','RES1','RES2','SYM2'])
                      for start, stop in zip(starts, stops)]

    X = []
    Y = []
    Z = []
    with open('isopropanol_%s_3.3.3_COM.xyz'%name,'w+') as f1:

        for i in range(len(dataframes)):
            f1.write(str(int(len(dataframes[i])/5))+'\n'+str('i=')+' '+str(i+1)+'\n')
            for k in range(int(len(dataframes[i])/5)):
                xf = dataframes[i]['X'][k*5:k*5+5]
                yf = dataframes[i]['Y'][k*5:k*5+5]
                zf = dataframes[i]['Z'][k*5:k*5+5]
                xcom = calc_com(xf)
                ycom = calc_com(yf)
                zcom = calc_com(zf)
                f1.write(str('C')+' '+str(xcom)+' '+str(ycom)+' '+str(zcom)+'\n')
                
## Call IPA_COM to calculate com and store in xyz file

IPA_COM('/bgfs/kjohnson/pbs13/FlexibleMOFs-GCMC/IPA/7loading/1/18/Movies/System_0/Movie_finalmovie-18_1.1.1_291.000000_1.000000_component_isopropanol_0.pdb','1Pa')

        
  COM_list_1Snap=np.loadtxt('isopropanol_1Pa_3.3.3_COM.xyz',skiprows=2,usecols=(1,2,3),max_rows=90) ## COM of a single snapshot from COM xyz file
  
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
cell_dimensions = {'x': 41.957, 'y': 41.957, 'z': 41.957}
# Select size of bins (cell_dim/del_x should equal an integer)
box_size = {'x': 1.048925, 'y': 1.048925, 'z': 1.048925}

# Call to center of mass histogram function
COM_histogram(COM_list_1Snap, cell_dimensions, box_size)
