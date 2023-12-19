import csv
import numpy as np
from random import random
import pandas as pd
import ast
import os
import pandas as pd

pressure = ['0.001','0.064','0.216']
snaps = ['3','5','7','9','12','13','15','17','18','20','21']


for i in pressure:
    i_temp = float(i)*100
    i_temprounded =round(i_temp,1) 
    x=i_temp*10**6
    y = format(x/(10**6),'.6f')
    print(round(i_temp,1))
    os.chdir(str(i_temprounded))
    for j in snaps:
        os.chdir(j)
        os.chdir('Hydrogenbonding')
        infile = open(f'{j}-Framework_0_final_1_1_1_P1.cif')
        newopen = open('Framework.txt', 'w')
        for line in infile :
                newopen.write(line)
        newopen.close()
        inFile = pd.read_csv('Framework.txt', sep='\s+')
        inFile.to_csv('Framework.csv', sep='\t')
        os.remove("Framework.txt")
        
        os.chdir('..')
        os.chdir('..')
    os.chdir('..')
