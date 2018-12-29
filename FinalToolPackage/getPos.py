'''Cell Expanding and Getting Final Atom Positions'''
#Author:William Song

import numpy as np 
def latticeExpand(oriPos):
    print("expanding the origin lattice to 3x3.")
    index=np.array([[x,y,z] for x in range(-1,2) for y in range(-1,2) for z in range(-1,2)])
    expandedPos=[]
    Positions=np.array(oriPos)
    
    
    for posn in Positions:
        for matrix in index:
            posm=posn+matrix
            expandedPos.append(posm)
   
    
    return np.array(expandedPos)
def getPos(allPos,lattice):
    
    newAllPos=latticeExpand(allPos)
    
    posList=np.array([])
    
    for pos in newAllPos:
        newPos = np.dot(pos,lattice)
        if len(posList) == 0:
                
            posList = newPos 
            
        else :    
            
            posList=np.vstack((posList,newPos))
            
    return posList