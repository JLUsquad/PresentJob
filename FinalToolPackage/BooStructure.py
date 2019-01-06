'''the boo class of a cif file.
this class inits with a path of a cif file.
the getBoostructure return a dic for which the keys are the names of atoms the values are the boo repersentations of relevant atom.  
'''
#Author:William Song
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
        
        bondsDics=Structure(self.lattice,self.positions,self.elements,self.atoms).getStructure()
        
        
        keys=bondsDics.keys()
        #the keys are the atoms.
        booStructure={}
        for key in keys:
            vecs=bondsDics[key]
            boos=[]
           
            for vec in vecs:
                boos.append(Boo(vec).getBoo())
            booStructure[key]=np.array(boos)
        return booStructure
        
        
