import itertools
from math import sqrt, ceil
import math
import csv
from random import random
import pandas as pd
import numpy as np
import os
import ast
import matplotlib.pyplot as plt
import re



def load_positions(filename,newline=''):
    with open(filename) as csvfile:
        output_data = csv.reader(csvfile, delimiter='\t')
        output = []
        for i in output_data:
            row = i
            output.append(row)
        return output
    
def calculate_distance(p1,p2):
    r_vect = p1 - p2
    r_mag = np.linalg.norm(r_vect)
    return r_mag


# change file name here 


file = 'Framework.csv'

output=load_positions(file)
# change x y z based on framework file and change index numbers until line 90 if x y z changed
x = 2
y = 3
z = 4

for i in range(1,len(output)):
    ax=float(output[i][x])
    ay=float(output[i][y])
    az=float(output[i][z])
    output[i][2]=ax
    output[i][3]=ay
    output[i][4]=az
a = len(output)


# if any atom is within 3 A of edges, I will duplicate it to its periodic position
ratio = 3/41.9568

duplicate = []


# I have fractional coordinates in the framework input file

for i in range (1,len(output)):
    if output[i][2]>1-ratio:
        array= [output[i][0], output[i][1], output[i][2]-1, output[i][3], output[i][4], output[i][5]]
        output.append(array)
    if output[i][2]<ratio:
        array= [output[i][0], output[i][1], output[i][2]+1, output[i][3], output[i][4], output[i][5]]
        output.append(array)

for i in range (1,len(output)):
    if output[i][3]>1-ratio:
        array= [output[i][0], output[i][1], output[i][2], output[i][3]-1, output[i][4], output[i][5]]
        output.append(array)
    if output[i][3]<ratio:
        array= [output[i][0], output[i][1], output[i][2], output[i][3]+1, output[i][4], output[i][5]]
        output.append(array)

for i in range (1,len(output)):
    if output[i][4]>1-ratio:
        array= [output[i][0], output[i][1], output[i][2], output[i][3], output[i][4]-1, output[i][5]]
        output.append(array)
    if output[i][4]<ratio:
        array= [output[i][0], output[i][1], output[i][2], output[i][3], output[i][4]+1, output[i][5]]
        output.append(array)





framework_position = []
for i in range (1,len(output)):
    position=np.array([float(item) for item in output[i][2:5]])
    framework_position.append(position)


#convert framework atoms' positions from fractional to cartesian

def get_fractional_to_cartesian_matrix(a, b, c, alpha, beta, gamma, angle_in_degrees=True):
    if angle_in_degrees:
        alpha = np.deg2rad(alpha)
        beta = np.deg2rad(beta)
        gamma = np.deg2rad(gamma)
    cosa = np.cos(alpha)
    sina = np.sin(alpha)
    cosb = np.cos(beta)
    sinb = np.sin(beta)
    cosg = np.cos(gamma)
    sing = np.sin(gamma)
    volume = 1.0 - cosa**2.0 - cosb**2.0 - cosg**2.0 + 2.0 * cosa * cosb * cosg
    volume = np.sqrt(volume)
    r = np.zeros((3, 3))
    r[0, 0] = a
    r[0, 1] = b * cosg
    r[0, 2] = c * cosb
    r[1, 1] = b * sing
    r[1, 2] = c * (cosa - cosb * cosg) / sing
    r[2, 2] = c * volume / sing
    return r

 # get the a b c and degree number from framwork cif file (not necessaty 1 by 1 by 1)
r=get_fractional_to_cartesian_matrix(41.9568,41.9568,41.9568, 90, 90, 90,angle_in_degrees=True)




for i in range (0,len(framework_position)):
    framework_position[i] = np.dot(r, framework_position[i])

#group added O, added H, mu3H, mu3O separately 

positions_Oa = []
for i in range (1,a):
    if 'OH' in output[i]:
        positions_Oa.append(framework_position[i-1])
        
positions_Ha = []
for i in range (1,a):
    if 'HO' in output[i]:
        positions_Ha.append(framework_position[i-1])

positions_mu3H = []
for i in range (1,a):
    if 'mu3H' in output[i]:
        positions_mu3H.append(framework_position[i-1])

positions_mu3O = []
for i in range (1,a):
    if 'mu3O' in output[i]:
        positions_mu3O.append(framework_position[i-1])


replicated_Oa = []
for i in range (a, len(output)):
    if 'OH' in output[i]:
        replicated_Oa.append(framework_position[i-1])
        
replicated_Ha = []
for i in range (a, len(output)):
    if 'HO' in output[i]:
        replicated_Ha.append(framework_position[i-1])

replicated_mu3H = []
for i in range (a, len(output)):
    if 'mu3H' in output[i]:
        replicated_mu3H.append(framework_position[i-1])

replicated_mu3O = []
for i in range (a, len(output)):
    if 'mu3O' in output[i]:
        replicated_mu3O.append(framework_position[i-1])



# get data from IPA file
with open("acetone.pdb", "r+") as f:
    p = re.compile("ENDMDL|CRYST")
    lines = [line for line in f.readlines() if p.search(line) is None]
    f.seek(0)
    f.truncate(0)
    f.writelines(lines)
inFile = pd.read_csv('acetone.pdb', sep='\s+')
inFile.to_csv('acetone.txt', sep='\t')

filename = 'acetone.txt'

#split the whole file into snapshots
with open(filename) as fin:
    a=1
    fout = open("output0.txt","w")
    for i,line in enumerate(fin):
        if 'MODEL' not in line:
            fout.write(line)
        if 'MODEL' in line:
            fout.close()
            fout = open("output%d.txt"%a,"w")
            a=a+1

    fout.close()
os.remove("output0.txt")


IPA_addedOH=[]
IPA_mu3OH=[]
IPA_both=[]
sum_addedOH=0
sum_mu3OH=0
sum_both=0

for snapshot in range (1,a):
    file = 'output%d.txt'%snapshot


    output_IPA=load_positions(file)



# group IPA_O and IPA_H separately
    positions_IPA_O = []
    for i in range (1,len(output_IPA)):
        if 'O' in output_IPA[i]:
            position=np.array([float(item) for item in output_IPA[i][4:7]])
            positions_IPA_O.append(position)

    positions_IPA_H = []
    for i in range (1,len(output_IPA)):
        if 'H' in output_IPA[i]:
            position=np.array([float(item) for item in output_IPA[i][4:7]])
            positions_IPA_H.append(position)


    total_IPA = len(positions_IPA_O)

    

    sigma_mix = 2.0 #setup cutoff value

    count_addedOH = 0
    count_u3OH = 0
    count_addedOH_u3OH =0
    
    #loop every IPA to check if it bonded to H or O

    for j in range (0, len(positions_IPA_O)):

        distance_a = 0
        distance_u = 0
        distance = 0
        list_a = []
        list_u = []
#         p3 = positions_IPA_H[j]
        p4 = positions_IPA_O[j]
        for i in range (0, len(positions_Oa)): 
            p1 = positions_Oa[i] 
            p2 = positions_Ha[i]
#             distance_Oa_IPAH=calculate_distance(p1,p3)
            distance_Ha_IPAO=calculate_distance(p2,p4)
            distance = distance_Ha_IPAO
#             distance = min(distance_Oa_IPAH, distance_Ha_IPAO)
            if distance < sigma_mix:
                list_a.append(distance)
        for i in range (0, len(replicated_Oa)):
            p1 = replicated_Oa[i]
#             distance_Oa_IPAH=calculate_distance(p1,p3)
            if distance<sigma_mix:
                list_a.append(distance)
        for i in range (0, len(replicated_Ha)):
            p2 = replicated_Ha[i]
            distance_Ha_IPAH=calculate_distance(p2,p4)
            if distance<sigma_mix:
                list_a.append(distance)
        if list_a:       
            distance_a = float(min(list_a))
        else: 
            distance_a = 0

        for i in range (0, len(positions_mu3O)):
            p1 = positions_mu3O[i] 
            p2 = positions_mu3H[i]
#             distance_u3O_IPAH=calculate_distance(p1,p3)
            distance_Hu3O_IPAO=calculate_distance(p2,p4)
            distance = distance_Hu3O_IPAO
#             distance = min(distance_u3O_IPAH,  distance_Hu3O_IPAO)
            if distance < sigma_mix:
                list_u.append(distance)
        for i in range (0, len(replicated_mu3O)):
            p1 = replicated_mu3O[i]
#             distance_u3O_IPAH=calculate_distance(p1,p3)
            if distance<sigma_mix:
                list_u.append(distance)
        for i in range (0, len(replicated_mu3H)):
            p2 = replicated_mu3H[i]
            distance_u3O_IPAH=calculate_distance(p2,p4)
            if distance<sigma_mix:
                list_u.append(distance)
        if list_u:
            distance_u = float(min(list_u))
        else: 
            distance_u = 0




        if distance_a == 0 and distance_u != 0: 
            count_u3OH +=1
        if distance_u == 0 and distance_a != 0:
            count_addedOH +=1
        if distance_u != 0 and distance_a != 0:
            count_u3OH +=1
            count_addedOH +=1
            count_addedOH_u3OH += 1



    total = count_u3OH+count_addedOH-count_addedOH_u3OH




    if total_IPA != 0:
        percentage_addedOH = count_addedOH/total_IPA
        percentage_u3OH = count_u3OH/total_IPA
        percentage_hydrogenbond=total/total_IPA
        IPA_addedOH.append(percentage_addedOH)
        IPA_mu3OH.append(percentage_u3OH)
        IPA_both.append(percentage_hydrogenbond)
        sum_addedOH+=percentage_addedOH
        sum_mu3OH+=percentage_u3OH
        sum_both+=percentage_hydrogenbond
    else:
        a=a-1
        


np.savetxt('Hydrogenbond_acetoneresults.txt', np.c_[IPA_addedOH,IPA_mu3OH,IPA_both])






