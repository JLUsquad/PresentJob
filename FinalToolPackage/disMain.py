'''an example of using booStructure'''
#Author:William Song
from distanceDistribution import distanceDistribution
import numpy as np
from matplotlib import pyplot as plt
DB=distanceDistribution('./1-ICSD-60650.cif').getStructure()
atoms = DB.keys()

for atom in atoms:
	plt.figure(figsize=(13.66,7.68))
	disDics={}
	distances = DB[atom]
	distances=sorted(distances)
	
	for dist in distances:
		if dist in disDics.keys():
			disDics[dist]+=1
		

		else:
			disDics[dist]=1
	mindis=min(disDics.keys())
	
	
	disRratio=[x/mindis for x in disDics.keys()]
	
	disValues=list(disDics.values())
	x=[]
	y=[]
	index=0
	for i in range(len(disRratio)):
		
		
		if i ==0:
			x.append(disRratio[i])
			y.append(disValues[i])
		else :
			if disRratio[i]-x[index] < 0.001:
				
				y[index]+=disValues[i]
			else :
				x.append(disRratio[i])
				y.append(disValues[i])
				index+=1
	
	plt.plot(x,y,'o')
	
	plt.xlabel('Distances(D/Dmin)')
	plt.ylabel('Atoms')
	path='./images/'+str(atom)+'.png'
	
	plt.savefig(path)
	plt.close()
print("saving images")