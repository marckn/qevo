"""
This module can be loaded to add more matrix definitions for the 
harmonic manybody wavefunctions. Just change the parameters below and 
import or python2 compute_integrals.py

"""

import basis.harmonic_precomputed as hp
import basis.harmonic_analyticals as ha

mcpoints=1000000
dim=2
eta=5.0

wx = (2/(1+eta**2))**0.5
wy = eta*wx

bpars=(wx,wy)
caps=(8,5)

print "Building integral tables for basis set parameters ", bpars
print "level caps: ", caps



def v1((x,y)):
   return 0.5*(x**2+y**2)

def v1x((x,y)):
   return 0.5*x*x


def v1_xanalytical((l1,l2)):
   return 0.5*ha.powx(2,l1[0],l2[0],bpars[0])
   
def v1_yanalytical((l1,l2)):
   return 0.5*ha.powx(2,l1[1],l2[1],bpars[1])


def v1x_zeros((l1,l2)):
   if abs(l1[0]-l2[0])%2 != 0:
      return True

   if l1[1] != l2[1]:
      return True
      
   
   return False

name="x**2_eta"+str(eta)
print "Building ",name
hp.writeFunc(name,v1x,bpars,caps,mcpoints,
              iszero=v1x_zeros, analytical =v1_xanalytical)


def v1y((x,y)):
   return 0.5*y*y

def v1y_zeros((l1,l2)):
   if abs(l1[1]-l2[1])%2 != 0:
      return True

   if l1[0] != l2[0]:
      return True
      
   
   return False

name="y**2_eta"+str(eta)
print "Building ",name
hp.writeFunc(name,v1y,bpars,caps,mcpoints,
             iszero=v1y_zeros,analytical=v1_yanalytical)


def v2((x1,y1),(x2,y2)):
   dist = (x1-x2)**2 + (y1-y2)**2
   return 1.0/dist**0.5

def v2_zeros((l1,l2,l3,l4)):
   ax=l1[0]+l2[0]
   bx=l3[0]+l4[0]
   ay=l1[0]+l2[0]
   by=l3[0]+l4[0]
   
   if abs(ax+bx)%2 != 0:
      return True
   
   if abs(ay+by)%2 != 0:
      return True
      
   if abs(l1[0]-l2[0]) > 2:
      return True
     
   if abs(l1[1]-l2[1]) > 2:
      return True

   if abs(l3[0]-l4[0]) > 2:
      return True
     
   if abs(l3[1]-l4[1]) > 2:
      return True
      
   return False

def v2_analytical((l1,l2,l3,l4)):
   return ha.c2d(l1,l2,l3,l4,bpars)

name="r12**-1_eta"+str(eta)
print "Building ",name
hp.writeFunc(name,v2,bpars,caps,mcpoints,iszero=v2_zeros, analytical=v2_analytical)
