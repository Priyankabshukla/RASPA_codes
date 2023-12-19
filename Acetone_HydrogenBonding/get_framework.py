import numpy as np
import os


pressure = ['0.001','0.064','0.216' ]
snaps = ['3','5','7','9','12','13','15','17','18','20','21']

def line_prepend(filename, line):
    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0,0)
        f.write(line.rstrip('\r\n') + '\n' + content)
        
prependThis='1\n'

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
        with open('Framework_0_final_1_1_1_P1.cif','r') as fin:
            data = fin.read().splitlines(True)
        with open(f'{j}-Framework_0_final_1_1_1_P1.cif','w') as fout:
            fout.writelines(data[28:])
        line_prepend(f'{j}-Framework_0_final_1_1_1_P1.cif',prependThis)
            
        os.chdir('..')
        os.chdir('..')
    os.chdir('..')
