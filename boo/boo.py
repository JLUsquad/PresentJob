import numpy as np
import math
from sphereHarmonics import Ylm
class Boo(object):
    def __init__(self,vec=None):
        # self.origin=origin
        # self.pos=pos
        self.vec=vec
        self.sphericalPos=self.turnIntoSphere()


    def turnIntoSphere(self):
        normalvec=np.array([0,0,1])
        # if self.vec ==None :
        #     vec3d=self.pos-self.origin
        # elif self.vec !=None :
        vec3d = np.array(self.vec) 
        # else:
        #     raise ValueError('wrong parameter')
        r=np.linalg.norm(vec3d)
        #print(vec3d)
        cosTheta=np.dot(vec3d,normalvec)/r/1.0
        theta=math.acos(cosTheta)
        vec2d=vec3d[:2]
        cosPhi=vec2d[0]/np.linalg.norm(vec2d)
        phi=math.acos(cosPhi)
        
        return (r,theta,phi)

    def getBoo(self):
        r,theta,phi=self.sphericalPos
        ySeries=[0,2,4,6,8,10]
        qSeries=np.array([])
        for l in ySeries:
            for m in range(-l,l+1):
                Qlm=abs(Ylm(l,m,theta,phi))
                qSeries=np.append(qSeries,math.sqrt(4*math.pi/(2*l+1)*(Qlm)**2))
        return qSeries

#print(Boo(np.array([1,1,1])).getBoo())
