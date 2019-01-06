'''the normal structure class of a cif file'''
#Author:William Song
from read import Read
import numpy as np 

import getPos
class Structure(object):
    def __init__(self,lattice,positions,elements,atomNums):
        self.lattice=lattice
        self.positions = positions
        self.elements= elements
        self.atoms=atomNums
        
    # def __init__(self,filePath,fileType=None):
    #     self.filePath=filePath
    #     self.fileType=fileType
    #     self.structure = Read(self.filePath,self.fileType).getStructure()
    #     self.lattice=self.structure['lattice']
    #     self.positions = self.structure['positions']
    #     self.elements=self.structure['elements']
    #     self.atoms=self.structure['numbers']
        
        
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
        
        
        atomPos=np.dot(np.array(atomPos),lattice)
        
        allPos=getPos.getPos(allPos,lattice)
        
        minLength=None
        allLength=np.array([])
        bonds=[]
        for pos in allPos:
            
            length=np.linalg.norm(pos-atomPos)
           
            allLength=np.append(allLength,length)
            if length > 0:
                if minLength==None :
                    if length > 0.1:
                        # A bond should be longer than 0.1 A
                        minLength=length
                    
                
                if length < minLength:
                    if length > 0.1:
                        minLength=length
                    
        for i in range(0,len(allLength)):
            if (allLength[i] - minLength )<0.0001:
                bond = allPos[i]-atomPos
                if  np.linalg.norm(bond) > 0.1:   
                    bonds.append(bond)
         
        return np.stack(bonds)
