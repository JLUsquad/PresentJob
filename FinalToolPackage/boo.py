'''
turn a vec into boo vec.
make sure you have install pyshtools (which is quite easy with 'sudo pip install pyshtools')

'''
#Author:William Song
import numpy as np
import math
import pyshtools

def step(n):
    if n==1 or n==0:
        return 1
    else:
        return n*step(n-1)
class Boo(object):
    def __init__(self,vec=None):
        # self.origin=origin
        # self.pos=pos
        self.vec=vec
        self.sphericalPos=self.turnIntoSphere()


    def turnIntoSphere(self):
        normalvec=np.array([0,0,1])
        #normalvec is normal vection
        vec3d = np.array(self.vec) 
       
        r=np.linalg.norm(vec3d)
       
        cosTheta=np.dot(vec3d,normalvec)/r/1.0
        theta=math.acos(cosTheta)
        vec2d=vec3d[:2]
        mod=np.linalg.norm(vec2d)
        phi=0
        if mod != 0:
            cosPhi=vec2d[0]/mod
            phi=math.acos(cosPhi)
        #if the subvector in the x-y plane the phi will be 0.
        
        return (r,theta,phi)

    def getBoo(self):
        r,theta,phi=self.sphericalPos
        ySeries=[0,2,4,6,8,10]
        qSeries=np.array([])
        for l in ySeries:
            pSeries=pyshtools.shtools.PlmSchmidt(lmax=l,z=math.cos(theta))
            #pSeries is a set of values with specific l and cos theta.
            ql=0
            for m in range(-l,l+1):
                index=int(l*(l+1)/2+m)
                #index is the location of schmidt legendre with specific m. 
                plm=pSeries[index]
                yModule=(-1)**abs(m)*math.sqrt(step(l-abs(m))*(2*l+1)/step(l+abs(m))/4/math.pi)*plm
                #yModule is the module of spherical fuction with given l, m and cos theta.
                if m>=0:
                    Ylm=complex(yModule*math.cos(m*phi),yModule*math.sin(m*phi))
                else:
                    Ylm=(-1)**m*complex(yModule*math.cos((-m)*phi),yModule*math.sin(-m*phi)*(-1))
                ql+=math.sqrt(4*math.pi/(2*l+1)*abs(Ylm))
            qSeries=np.append(qSeries,ql)
        return qSeries
