import ase
from ase.io import read, write
import os
import numpy as np

pressure = ['0.064' ]
snap = ['3','5','7','9','12','13','15','17','18','20','21']

for j in pressure:
    j_temp = float(j)*100
    j_temprounded =round(j_temp,1) 
    print('j:',j)
    print('j_temp:',j_temp)
    x=j_temp*10**6
    y = format(x/(10**6),'.6f')
    print('Y:',y)
    print(round(j_temp,1))
    os.chdir(str(j_temprounded))
    for i in snap:
        os.chdir(i)
        os.mkdir('Hydrogenbonding')
        os.chdir('Hydrogenbonding')
        os.system('cp /ihome/kjohnson/pbs13/FlexibleFramework-resultanalysis/acetone-HB/HB_7loading/get-framework-step2.py .')
        os.system('cp /ihome/kjohnson/pbs13/FlexibleFramework-resultanalysis/acetone-HB/HB_7loading/convert_file_type.ipynb .')
        os.system('cp /ihome/kjohnson/pbs13/FlexibleFramework-resultanalysis/acetone-HB/HB_7loading/mainscript-step4.py .')
        os.system(f'cp /bgfs/kjohnson/pbs13/FlexibleMOFs-GCMC/acetone_Rogge_7Loadingacetone_final/acetone_severalsnaps/{str(j_temprounded)}/{i}/Movies/System_0/Framework_0_final_1_1_1_P1.cif .')
        os.system(f'cp /bgfs/kjohnson/pbs13/FlexibleMOFs-GCMC/acetone_Rogge_7Loadingacetone_final/acetone_severalsnaps/{str(j_temprounded)}/{i}/Movies/System_0/Movie_finalmovie-{i}_1.1.1_298.000000_{y}_component_acetone_0.pdb .')
        os.rename(f'Movie_finalmovie-{i}_1.1.1_298.000000_{y}_component_acetone_0.pdb','acetone.pdb')
        
        os.chdir('..')
        os.chdir('..')
    os.chdir('..')
