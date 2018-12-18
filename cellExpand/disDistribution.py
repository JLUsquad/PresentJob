from distanceDistribution import distanceDistribution
import numpy as np
from matplotlib import pyplot as plt
DB=distanceDistribution('/home/william/workspace/boo/1-ICSD-108.cif').getStructure()
atoms = DB.keys()

for atom in atoms:
	plt.figure(figsize=(13.66,7.68))
	distances = DB[atom]
	plt.plot(distances,'o-',label=atom)
	plt.xlabel('Atoms')
	plt.ylabel('Distances(Å)')
	path='./images/'+atom+'.png'
	plt.savefig(path)
	plt.close()
