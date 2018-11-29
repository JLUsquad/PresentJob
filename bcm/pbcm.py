# coding: utf-8

from read import Read
import fbcm
import numpy as np
import os


class Bcm(object):
    '''

    '''

    def __init__(self,type=None):

        self.type = type


    def bcm_calc(self,filepath):

        structure = Read(filepath,self.type).getStructure()
        nele = len(structure['elements'])        # nele is numble of elements
        natom = sum(structure['numbers'])
        lattice = structure['lattice']
        positions = structure['positions']
        numbers = structure['numbers']
        bcmatrix = fbcm.bcm(nele,natom,positions,lattice,numbers)
        bcmatrix = np.array(bcmatrix)

        return bcmatrix

    def calc_folder(self,folderpath):
        resultstore = []
        k = 0
        for root, dirs, files in os.walk(folderpath):
            for file in files:
                try:
                    test = Read(file)
                except Exception:
                    continue
                k = k + 1
                filepath = os.path.join(root, file)
                print (filepath)
                bcmstore = []
                bcmatrix = self.bcm_calc(filepath)
                bcmstore.append(bcmatrix)
                bcmstore.append(filepath)
                resultstore.append(bcmstore)
        print ('total num = ', k)
        return resultstore




    def distxy(self,x,y):

        if x.shape != y.shape:
            print ('can not compare two matrix')
            exit()
            d = np.linalg.norm(x-y)
            return d

    def dist0(self,x):

        a = np.linalg.norm(x)

        return a


   
