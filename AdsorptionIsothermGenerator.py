import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib
import subprocess
import math


directory_to_check= os.getcwd()
averageN = np.array([],dtype=float)
final_averageN = np.array([],dtype=float)
Nerror = np.array([],dtype = float)
Nexcesslist = np.array([],dtype = float)
Nexcess_error = np.array([],dtype = float)

def myfile(filePath1):
    global averageN, Nlist, Nerror, Nexcess, Nexcesslist, Nexcess_error, final_averageN, Nerror_list
    if os.path.isfile(filePath1):
        f = open(filePath1,"r")
        f1 = f.readlines()
        Nlist = np.array([],dtype=float)
        Nexcess = np.array([],dtype=float)
        
        for n in np.r_[21755:len(f1)]:
            if "Average loading absolute" in f1[n] and "mol/kg" in f1[n] and "excess" not in f1[n]:
#                 print('Hi')
                i3 = f1[n].split()
                print(i3)
                Nlist = np.append(Nlist,float(i3[5]))
          
        averageN = np.append(averageN,Nlist)
        Nerror = np.append(Nerror,float(i3[7]))
    
Pressure1 = []
snapshot_number = []
P = np.array([],dtype = float)

def myfunction(directory):
    global Pressure1,P,final_averageN, snapshot_number
    
    
    for filename in os.listdir(directory):
        
        if filename.endswith(".data") and "1.1.1" in filename:
            print(directory+"/"+filename)
            myfile(directory+"/"+filename)
    
    X=directory.split('/')
    Pressure1 = np.append(Pressure1,float(X[7]))
    snapshot_number = np.append(snapshot_number,int(X[8]))
    
            
directories = [os.path.abspath(x[0]) for x in os.walk(directory_to_check)]
                                                                                                                

for i in directories:
        if "Output" in i and "System_0" in i:
                os.chdir(i)         # Change working Directory
                myfunction(i)
                
