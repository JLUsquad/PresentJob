import numpy as np 
def latticeExpand(oriLattice):
    index=[-1,0,1]
    lattices=[]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                #if not (i==0 and j==0 and k==0):
                lattice=np.array(oriLattice)
                # print(lattice)
                for l in range(3):
                    
                    lattice[l][0]+=index[i]
                    lattice[l][1]+=index[j]
                    lattice[l][2]+=index[k]
                
            
                lattices.append(lattice)
    return np.array(lattices)

            
