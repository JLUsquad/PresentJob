'''Distance Distribution Class of A Cif File'''
#Author:William Song

from read import Read
import numpy as np 
import getPos
class distanceDistribution(object):
    def __init__(self,filePath,fileType=None):
        self.filePath=filePath
        self.fileType=fileType
        self.structure = Read(self.filePath,self.fileType).getStructure()
        self.lattice=self.structure['lattice']
        self.positions = self.structure['positions']
        self.elements=self.structure['elements']
        self.atoms=self.structure['numbers']
        print("reading cif file:"+self.filePath)
    def getStructure(self):
        print("reading atoms")
        structure={}
        natoms=len(self.atoms)
        for eleNumber in range(1,natoms+1):
            for atomNumber in range(1,self.atoms[eleNumber-1]+1):
                atomName=str(self.elements[eleNumber-1])+str(atomNumber)
                atomNumInAll=0
                for i in range(0,eleNumber-1):
                    atomNumInAll+=self.atoms[i]
                atomNumInAll+=atomNumber
                print("reading the bonds of "+atomName)
                structure[atomName]=self.getBonds(self.positions[atomNumInAll-1],self.positions,self.lattice)
        return structure
    def getBonds(self,atomPos,allPos,lattice):
        
        
        atomPos=np.dot(atomPos,lattice)
        
        allPos = getPos.getPos(allPos,lattice)
        
        allLength=np.array([])
       
        for pos in allPos:
            
            length=np.linalg.norm(pos-atomPos)
            
            if length > 0.001:
                allLength=np.append(allLength,length)
            

        return allLength
        
        