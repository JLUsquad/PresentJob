import numpy as np 
class Structure(object):
    def __init__(self,atoms,atomNums,posMatrix,lattice):
        self.posMatrix=posMatrix
        self.atoms=atoms
        self.atomNums=atomNums
        self.lattice=lattice
        
    def getStructure(self):
        structure={}
        natoms=len(self.atoms)
        for eleNumber in range(1,natoms+1):
            for atomNumber in range(1,self.atomNums[eleNumber-1]+1):
                atomName=str(self.atoms[eleNumber-1])+str(atomNumber)
                atomNumInAll=0
                for i in range(0,eleNumber-1):
                    atomNumInAll+=self.atomNums[i]
                atomNumInAll+=atomNumber
                structure[atomName]=self.getBonds(self.posMatrix[atomNumInAll-1],self.posMatrix,self.lattice)
        return structure
    def getBonds(self,atomPos,allPos,lattice):
        atomPos=np.dot(np.array(atomPos),lattice)
        allPos=np.array(allPos)
        minLength=None
        allLength=np.array([])
        bonds=np.array([])
        for pos in allPos:
            
            length=np.linalg.norm(np.dot(pos,lattice)-atomPos)
            allLength=np.append(allLength,length)
            if length > 0:
                if minLength==None:
                    minLength=length

                elif length < minLength:
                    minLength=length
        for i in range(0,len(allLength)):
            if allLength[i] == minLength :
                if len(bonds)==0:
                    
                    bonds=np.array([np.dot(allPos[i],lattice)-pos])
                else:
                    
                    bonds=np.vstack((bonds,np.dot(allPos[i],lattice)-pos))
        #print(bonds,'a')
        return bonds


            
            
