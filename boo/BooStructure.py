from read import Read
from boo import Boo
from structure import Structure
import numpy as np 
class BooStructure(object):
    def __init__(self,filePath,fileType=None):
        self.filePath=filePath
        self.fileType=fileType
        self.structure = Read(self.filePath,self.fileType).getStructure()
        self.lattice=self.structure['lattice']
        self.positions = self.structure['positions']
        self.elements=self.structure['elements']
        self.atoms=self.structure['numbers']
    def getBooStructure(self):
        natoms=sum(self.atoms)
        strcVecs=Structure(self.elements,self.atoms,self.positions,self.lattice).getStructure()
        keys=strcVecs.keys()
        booStructure={}
        for key in keys:
            vecs=strcVecs[key]
            boos=np.array([])
            for vec in vecs:
                boos=np.append(boos,Boo(vec).getBoo())
            booStructure[key]=boos
        return booStructure
        
        
