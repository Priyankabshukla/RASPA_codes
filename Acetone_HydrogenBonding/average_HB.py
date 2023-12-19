import ase
from ase.io import read, write
import os
import numpy as np

pressure = ['0.001','0.064','0.216']

snaps = ['3','5','7','9','12','13','15','17','18','20','21']
largest = np.array([],dtype=float)
big=np.array([],dtype=float)
error=np.array([],dtype=float)
k=1
# error1=np.array([],dtype=float)

for i in pressure:
    i_temp = float(i)*100
    i_temprounded =round(i_temp,1) 
    x=i_temp*10**6
    y = format(x/(10**6),'.6f')
#     print(round(i_temp,1))
    os.chdir(str(i_temprounded))
    k=1
    print(k)
    for j in snaps:
        os.chdir(j)
        os.chdir('Hydrogenbonding')
        data=np.array([],dtype=float)
        with open("Hydrogenbond_acetoneresults.txt") as f:
            d = [line.split() for line in f]
        for l in range(len(d)):
            data = np.append(data,float(d[l][2]))
        avg = np.mean(data)
        print("avg:",avg)
        big=np.append(big,avg)
        if k == 10:
            print(k)
            largest_avg= np.mean(big)
            largest=np.append(largest,largest_avg)
        print("big:",big)
        print("largest:",largest)
        np.savetxt('/ihome/kjohnson/pbs13/FlexibleFramework-resultanalysis/acetone-HB/HB_7loading/Average-7load-allsnaps-acetone-HB.txt', np.c_[big])
        np.savetxt('/ihome/kjohnson/pbs13/FlexibleFramework-resultanalysis/acetone-HB/HB_7loading/average-7loading-allpressure-acetone-HB.txt', np.c_[largest])
        
        k=k+1
        os.chdir('..')
        os.chdir('..')
   
    os.chdir('..')
