import numpy as np
import ase
from ase import Atoms
from ase.io import read, write
from collections import Counter
import matplotlib.pyplot as plt

shift = 0
dump=open('/bgfs/kjohnson/pbs13//dump.7LoadIPA_UiO66_Framework_Only.lammpstrj')
# out = open(out_file,'w+')
count = 0
direct_count = 0
start = 0
tag = 0
frame = []
atm = []
for l in dump:
#     l = lines.split()
    if l.split()[0] == 'ITEM:':
        if l.split()[1] == 'TIMESTEP':
            count = 0

    if tag == 1:
        atm.append([int(l.split()[0]),int(l.split()[2]),float(l.split()[7]) + shift,float(l.split()[8]) + shift,float(l.split()[9]) + shift])

        count += 1

    if count == 3648:
        tag = 0
        frame.append(atm)
#         break
        atm = []


    if l.split()[0] == 'ITEM:':
        if l.split()[1] == 'ATOMS':
            direct_count += 1
            tag = 1
            count = 0
frame = np.array(frame)
frame_sort = []
for i in range (len(frame)):
        fr_i = frame[i]
        fr_is = fr_i[fr_i[:,0].argsort()]
        frame_sort.append(fr_is)
frame_sort = np.array(frame_sort)

### Sort rows by column 2 (Type of atoms in dump.7loadedIPA)

frame_sort_AtomTypes = []
for i in range (len(frame_sort)):
        fr_i_AtomTypes = frame_sort[i]
        fr_is_AtomTypes = fr_i_AtomTypes[fr_i_AtomTypes[:,1].argsort()]
        frame_sort_AtomTypes.append(fr_is_AtomTypes)
frame_sort_AtomTypes = np.array(frame_sort_AtomTypes)

## Count number of occurrences of distinct atom types

items = Counter(frame_sort_AtomTypes[:,:,1].tolist()[0]).keys()
print(items)
values = Counter(frame_sort_AtomTypes[:,:,1].tolist()[0]).values()
print(values)

## Once sorted the data file by atom types, Use ase to convert fraction coordinates to cartesian coordinates

# Give atom types distinct identity. #1 (He):C1, #2 (Ca):C2, #3 (Be):Cca,
                                    #4 B:H2, #5 (H): Hoh, #6 P:Oca, #7 O: Ooh, #8: K (Oox), #9: Zr(Zr)
## For a single file
frac_to_cart_frame_sort_AtomTypes = Atoms('He384Ca768Be384B768H128P768O128K128Zr192',pbc=True,cell=[41.9568,41.9568,41.9568])
frac_to_cart_frame_sort_AtomTypes.set_scaled_positions(frame_sort_AtomTypes[:,:,2:][0])

## Scaled positions for snapshot #18 and its distance matrix
big_frame=[]
big_distance=[]
for i in range(18,19):
    frac_to_cart_frame_sort_AtomTypes = Atoms('He384Ca768Be384B768H128P768O128K128Zr192',pbc=True,cell=[41.9568,41.9568,41.9568])
    frac_to_cart_frame_sort_AtomTypes.set_scaled_positions(frame_sort_AtomTypes[:,:,2:][i])
    Distance = frac_to_cart_frame_sort_AtomTypes.get_all_distances(mic=True)
    big_frame+=[frac_to_cart_frame_sort_AtomTypes.get_positions()] # Add positions as elements in the list
    big_distance += [Distance]
    
    
## Find mu3OH H starting position
#He384Ca768Be384B768H128P768O128K128Zr192
# Give atom types distinct identity. #1 (He):C1, #2 (Ca):C2, #3 (Be):Cca,
                                    #4 B:H2, #5 (H): Hoh, #6 P:Oca, #7 O: Ooh, #8: K (Oox), #9: Zr(Zr)
mu3O_start = 384+768+384+768+128+768+128

def Print6Smallest(a,n):
    VALUE = max(a) + 1
    idx_min0 = a.index(min(a))
    print("idx_min0: ",idx_min0)
    min0 = a[idx_min0]
    print("min0: ", min0)
    a[idx_min0] = VALUE
    print(a[idx_min0])

    idx_min1 = a.index(min(a))
    print("idx_min1:",idx_min1)
    min1 = a[idx_min1]
    a[idx_min1] = VALUE

    idx_min2 = a.index(min(a))
    min2 = a[idx_min2]
    a[idx_min2] = VALUE
    
    idx_min3=a.index(min(a))
    min3=a[idx_min3]
    a[idx_min3] = VALUE
    
    idx_min4=a.index(min(a))
    min4 = a[idx_min4]
    a[idx_min4] = VALUE
    
    idx_min5 = a.index(min(a))
    min5=a[idx_min5]
    
    

    a[idx_min0] = min0
    a[idx_min1] = min1
    a[idx_min2] = min2
    a[idx_min3] = min3
    a[idx_min4] = min4
    a[idx_min5] = min5
    
    print([min0,min1,min2,min3,min4,min5],[idx_min0,idx_min1,idx_min2,idx_min3,idx_min4,idx_min5],[frac_to_cart_frame_sort_AtomTypes.get_chemical_symbols()[idx_min0],
                                                        frac_to_cart_frame_sort_AtomTypes.get_chemical_symbols()[idx_min1],frac_to_cart_frame_sort_AtomTypes.get_chemical_symbols()[idx_min2],
                                                                      frac_to_cart_frame_sort_AtomTypes.get_chemical_symbols()[idx_min3],
                                                                          frac_to_cart_frame_sort_AtomTypes.get_chemical_symbols()[idx_min4],frac_to_cart_frame_sort_AtomTypes.get_chemical_symbols()[idx_min5]])
    
    return [min0,min1,min2,min3,min4,min5],[idx_min0,idx_min1,idx_min2,idx_min3,idx_min4,idx_min5],[frac_to_cart_frame_sort_AtomTypes.get_chemical_symbols()[idx_min0],
                                                        frac_to_cart_frame_sort_AtomTypes.get_chemical_symbols()[idx_min1],frac_to_cart_frame_sort_AtomTypes.get_chemical_symbols()[idx_min2],
                                                                      frac_to_cart_frame_sort_AtomTypes.get_chemical_symbols()[idx_min3],
                                                                          frac_to_cart_frame_sort_AtomTypes.get_chemical_symbols()[idx_min4],frac_to_cart_frame_sort_AtomTypes.get_chemical_symbols()[idx_min5]]

# Given atom types distinct identity. #1 (He):C1, #2 (Ca):C2, #3 (Be):Cca,
                                    #4 B:H2, #5 (H): Hoh, #6 P:Oca, #7 O: Ooh, #8: K (Oox), #9: Zr(Zr)
## Print 6 Nearest neighbors for all mu3O H in rigid UiO66
big_NN_Atomdist=[]
big_NN_AtomType=[]

for i in range(mu3O_start,mu3O_start+128):
    NN_Atomdist_mu3O, index_NN, NN_AtomType_mu3O= Print6Smallest(big_distance[0][i].tolist(),len(big_distance[0][i]))
    print('\n')
    big_NN_Atomdist += [NN_Atomdist_mu3O]
    big_NN_AtomType += [NN_AtomType_mu3O]
    
 # Given atom types distinct identity. #1 (He):C1, #2 (Ca):C2, #3 (Be):Cca,
                                    #4 B:H2, #5 (H): Hoh, #6 P:Oca, #7 O: Ooh, #8: K (Oox), #9: Zr(Zr)
for i in range(len(big_NN_AtomType)):
    for j in range(0,6):
        if big_NN_AtomType[i][j]== 'He':
            big_NN_AtomType[i][j]='C1'
        elif big_NN_AtomType[i][j] == 'Ca':
            big_NN_AtomType[i][j] = 'C2'
        elif big_NN_AtomType[i][j] == 'Be':
            big_NN_AtomType[i][j] = 'Cca'
        elif big_NN_AtomType[i][j] == 'B':
            big_NN_AtomType[i][j] = 'H2'
        elif big_NN_AtomType[i][j] == 'H':
            big_NN_AtomType[i][j] = 'Hoh'
        elif big_NN_AtomType[i][j] == 'P':
            big_NN_AtomType[i][j] = 'OCa'
        elif big_NN_AtomType[i][j] == 'O':
            big_NN_AtomType[i][j]= 'Ooh'
        elif big_NN_AtomType[i][j] =='K':
            big_NN_AtomType[i][j] = 'Oox'
    
for i in range(len(big_NN_AtomType)):
    plt.plot(big_NN_AtomType[i],'ro')
plt.ylabel('NN Atom Types (7 loaded IPA snpashot #18)')
plt.xlabel('NN')

plt.figure(figsize=(10,10))
plt.plot(np.array(big_NN_Atomdist)[:,1],'bo',label='1 NN - Zr',alpha=1,markersize=5,markerfacecolor='none')
plt.plot(np.array(big_NN_Atomdist)[:,2],'bo',label='2 NN - Zr',alpha=1,markersize=5,markerfacecolor='none')
plt.plot(np.array(big_NN_Atomdist)[:,3],'bo',label='3 NN - Zr',alpha=1,markersize=5,markerfacecolor='none')
plt.plot(np.array(big_NN_Atomdist)[:,4][np.array(big_NN_AtomType)[:,4]=='Ooh'],'ro',label='4 NN - Ooh',alpha=1,markersize=5,markerfacecolor='none')
plt.plot(np.array(big_NN_Atomdist)[:,4][np.array(big_NN_AtomType)[:,4]=='Oox'],'go',label='4 NN - Oox',alpha=1,markersize=5,markerfacecolor='none')


plt.plot(np.array(big_NN_Atomdist)[:,5][np.array(big_NN_AtomType)[:,5]=='Ooh'],'ro',label='5 NN - Ooh',alpha=1,markersize=5,markerfacecolor='none')
plt.plot(np.array(big_NN_Atomdist)[:,5][np.array(big_NN_AtomType)[:,5]=='Oox'],'go',label='5 NN - Oox',alpha=1,markersize=5,markerfacecolor='none')


plt.ylabel(r'NN Atom Distances ($\rm \AA$) ',size=25)
plt.xlabel('$\mu_3$-O index',size=25)
plt.legend(frameon=False, prop = {'size':23},fontsize='x-large',bbox_to_anchor=(1, 1))
plt.tick_params(direction = 'in', right = False, top = False)
plt.tick_params(axis='both',labelsize = 25,length=10, width=2)




