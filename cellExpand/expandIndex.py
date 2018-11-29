import numpy as np 

index=[-1,0,1]
lattices=np.array([[[]]])
for i in range(3):
    for j in range(3):
        for k in range(3):
            #if not (i==0 and j==0 and k==0):
            for l in range(3):
                lattice=np.ones((3,3))
                lattice[l][0]+=index[i]
                lattice[l][1]+=index[j]
                lattice[l][2]+=index[k]
                lattices=np.append(lattices,lattice)
print(lattices)

            
