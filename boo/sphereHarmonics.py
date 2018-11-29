'''
This module include legendre function(legendre),associate legendre function (Plm),
and sphere harmonics(Ylm) which return a complex.
'''
import sympy
import math
import pyshtools
import plm
def step(n):
    if n==1 or n==0:
        return 1
    else:
        return n*step(n-1)
def legendre(l):
    x=sympy.Symbol('x')
    # if n==0:
    #     return 1
    # elif n==1:
    #     return x
    # else:
    
    #     result=(2*n-1)*x*legendre(n-1)/n-(n-1)*legendre(n-2)/n
        
    #     return result 
    origin=(x**2-1)**l
    for i in range(l):
        origin = sympy.diff(origin,x)
    result = origin/2**l/step(l)
    return result 

def assLegendre(l,m):
    # x=sympy.Symbol('x')
    # pl=legendre(l)
    # if m==0:
    #     return pl
    # else:
    #     return sympy.diff(assLegendre(l,m-1),x)
    x=sympy.Symbol('x')
    origin=(x**2-1)**l
    for i in range(l+abs(m)):
        origin = sympy.diff(origin,x)
    result = origin/(2**l)/step(l)*(1-x**2)**(abs(m)/2)
    return result
def Plm(l,m,cos):
    index=int(l*(l+1)/2 +m)
    
    result=pyshtools.shtools.PlmSchmidt(lmax=l,z=cos,csphase=1,cnorm=0)[index]
    return result
def newYlm(l,m,theta,phi):
    result=(-1)**abs(m)*math.sqrt(step(l-abs(m))*(2*l+1)/step(l+abs(m))/4/math.pi)*Plm(l,m,math.cos(theta))
    if m >=0:
        return complex(result*math.cos(m*phi),result*math.sin(m*phi))
    else:
        return (-1)**m*complex(result*math.cos((-m)*phi),result*math.sin(-m*phi)*(-1))

def Ylm(l,m,theta,phi):
    x=sympy.Symbol('x')
    module=(-1)**abs(m)*math.sqrt(step(l-abs(m))*(2*l+1)/step(l+abs(m))/4/math.pi)*assLegendre(l,abs(m))
    
    if m >=0 :
        return complex(module.evalf(subs={x:math.cos(theta)})*math.cos(m*phi),module.evalf(subs={x:math.cos(theta)})*math.sin(m*phi))
    else :
        result = (-1)**m*complex(module.evalf(subs={x:math.cos(theta)})*math.cos((-m)*phi),module.evalf(subs={x:math.cos(theta)})*math.sin(-m*phi)*(-1))
        return result





#print(assLegendre(2,0))
#print(Plm(2,0))
print(Ylm(2,-1,math.pi/4,math.pi/4))
print(newYlm(2,-1,math.pi/4,math.pi/4))
#print(math.sqrt(3/4/math.pi)*0.5)
#print(legendre(3))
#print(Plm(10,5,0.5))