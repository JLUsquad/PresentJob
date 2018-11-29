from read import Read
import numpy as np 
class Structure(object):
    def __init__(self,filePath,fileType=None):
        self.filePath=filePath
        self.fileType=fileType
        self.structure = Read(self.filePath,self.fileType).getStructure()
        self.lattice=self.structure['lattice']
        self.positions = self.structure['positions']
        self.elements=self.structure['elements']
        self.atoms=self.structure['numbers']
        
        
    def getStructure(self):
        structure={}
        natoms=len(self.atoms)
        for eleNumber in range(1,natoms+1):
            for atomNumber in range(1,self.atoms[eleNumber-1]+1):
                atomName=str(self.elements[eleNumber-1])+str(atomNumber)
                atomNumInAll=0
                for i in range(0,eleNumber-1):
                    atomNumInAll+=self.atoms[i]
                atomNumInAll+=atomNumber
                structure[atomName]=self.getBonds(self.positions[atomNumInAll-1],self.positions,self.lattice)
        return structure
    def getBonds(self,atomPos,allPos,lattice):
        index=[-1,0,1]
        lattices=[]
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        latticeO=lattice
                        latticeO[l][0]+=index[i]
                        latticeO[l][1]+=index[j]
                        latticeO[l][2]+=index[k]
                        lattices.append(latticeO)
        atomPos=np.dot(np.array(atomPos),lattice)
        allposlist=[np.dot(allPos,lattices[i]) for i in range(len(lattices))]
        allPos=np.stack(allposlist)
        shape=allPos.shape
        allPos=allPos.reshape(shape[0]*shape[1],shape[2])
        minLength=None
        allLength=np.array([])
        bonds=[]
        for pos in allPos:
            
            length=np.linalg.norm(pos-atomPos)
            allLength=np.append(allLength,length)
            if length > 0:
                if minLength==None:
                    minLength=length
                    #print(1)

                elif length < minLength:
                    minLength=length
                    #print(2)
        for i in range(0,len(allLength)):
            if allLength[i] == minLength :
                
                    
                bonds.append(np.dot(allPos[i],lattice)-pos)
                
        #print(len(allLength))
        #print(len(bonds))
        return np.stack(bonds)
s=Structure('/home/william/workspace/boo/1-ICSD-108.cif').getStructure()
print(s)