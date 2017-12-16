import basis.planewaves.manybody as bpw
import quantum.framework as q
import math
from os import mkdir
import numpy as np
from aux.utils import printData as fprint
import sys
import scipy.constants as const


L=[25.0,25.0]
nks=[32,32]
N=1
lamb=0.5
dodens=True
b=bpw.bset(L,nks,N,"Distinguishable",256)

try:
   mkdir("output")
except:
   pass
   

def vdwell((x,y)):
   a=0.5
   b=4.0
   
   c=0.5
   d=4.0
   
   return a*x**4 - b*x**2 + c*y**4 - d*y**2

def vho((x,y)):
   kx=1.0
   ky=3.0
   return 0.5*kx*x**2 + 0.5*ky*y**2


def gauss(X,X0,SIG):
   grt=1.0
   for x,x0,sig in zip(X,X0,SIG):
      grt = grt*math.exp(-(x-x0)**2/(2*sig**2))
      
   return grt

def vper((x,y)):
   def vcom(X):
      g1= gauss(X,(-L[0]/4,L[1]/3),(1.0,3.0))
      g2= gauss(X,(-L[0]/5,-L[1]/4),(3.0,1.0))
      g3= gauss(X,(L[0]/4,L[1]/3),(0.7,2.0))
      g4= gauss(X,(L[0]/5,-L[1]/4),(0.7,0.7))
      g5= gauss(X,(L[0]/3,-L[1]/10),(1.2,1.2))
      g6= gauss(X,(2*L[0]/5,-2*L[1]/5),(0.4,0.9))
      g7= gauss(X,(0,0),(1.2,0.6))
      g8= gauss(X,(-1.6,8.8),(0.3,0.3))
      g9= gauss(X,(1.6,8.8),(0.3,0.3))
      g10= gauss(X,(-1.6,5.0),(0.3,0.3))
      g11= gauss(X,(1.6,5.0),(0.3,0.3))
      
      
      v =      2.5*g1 + 3.0*g2 + 4*g3 + 8*g4 + 2*g5 + 3*g6 + 8*g7
      v = v +  8*g8 + 8*g9 + 8*g10 + 8*g11
      
      return 10*v

   r0=(x,y)
   r1 = (x+L[0],y)
   r2 = (x-L[0],y)
   r3 = (x,y+L[1])
   r4 = (x,y-L[1])
   r5 = (x+L[0],y+L[1])
   r6 = (x+L[0],y-L[1])
   r7 = (x-L[0],y+L[1])
   r8 = (x-L[0],y-L[1])
   
   return (vcom(r0)+vcom(r1)+vcom(r2)+vcom(r3)+vcom(r4)+vcom(r5)+vcom(r6)+vcom(r7)+vcom(r8))/9


Npoints=(200,200)
dr=(L[0]/(Npoints[0]-1),L[1]/(Npoints[1]-1))
Xs = [ (-L[0]/2 + i*dr[0]) for i in xrange(Npoints[0]) ]
Ys = [ (-L[1]/2 + i*dr[1]) for i in xrange(Npoints[1]) ]

Uext = np.zeros(Npoints)
for i,x in enumerate(Xs):
   for j,y in enumerate(Ys):
      Uext[i,j]= vper((x,y))

fprint("output/uext",Uext,(Xs,Ys)) 
print "Uext saved in output"
   


def vext((x,y)):
   return vper((x,y))


def v2((x1,y1),(x2,y2)):
   dist = (x1-x2)**2 + (y1-y2)**2
   return 1.0/dist**0.5


b.prepareMatrices(vext,None)
qmech = q.qfr(b,vext,lamb)


tolerance=1E-2
states=[]
ens=qmech.Hdiag.tolist()
curr = ens[0]
avgen=curr
lidx=[0]
for ii,en in enumerate(ens[1:]):
   if abs(curr-en) < tolerance:
      lidx.append(ii+1)
      avgen=avgen+curr
   else:
      avgen = avgen/len(lidx)
      states.append((avgen,tuple(lidx)))
      avgen=en
      lidx=[ii+1]
   
   curr=en   

if len(lidx) > 0:
   avgen = avgen/len(lidx)
   states.append((avgen,tuple(lidx)))
   
print "List of states:"
for i,s in enumerate(states):
   ln="{0:4d}) En={1:3f} Degeneracy: {2:3d}".format(i,s[0],len(s[1]))
   print ln



# initial state, N=1

def igauss(x,x0,sigma):
   retval=1.0
   norm=1.0
   for cx,cx0,csig in zip(x,x0,sigma):
      norm = norm*csig*(math.pi)**0.5
      retval = retval*math.exp(-(cx-cx0)**2/(2*csig**2))

   return retval/norm**0.5

from basis.planewaves.utils import fft
sigma=(1.0,1.0)
x0=(0,-6.9)
def wrapg(x):
   return igauss(x,x0,sigma)

fpnts=128
ft = fft(wrapg,fpnts,L)
cstate=np.zeros(b.size,dtype=complex)
cstatenorm=0
for i in xrange(b.size):
   ks = b._qnumbers[i][0]
   ks=(ks[0]-nks[0]/2,ks[1]-nks[0]/2)
   
   ks=(ks[0]-fpnts/2,ks[1]-fpnts/2)
   cstate[i]=ft[ks]
   cstatenorm= cstatenorm + cstate[i].real**2 + cstate[i].imag**2

cstatenorm = math.sqrt(cstatenorm)
cstate=cstate/cstatenorm

rho_points=(120,120)
rho_boundaries=[ (-12.5,12.5),
                 (-12.5,12.5) ]


dt=0.002
Nt=500

for i in xrange(Nt):
   print "Processing frame ",i
   ct = i*dt
   fout="snap_"+str(i)
   if i > 122:
      b.print_PSIn(cstate.tolist(),0,fout,"mov",rho_boundaries,rho_points,TOLERANCE=1E-3)
  
   qmech.settime(ct)
   
   cstate = np.dot(qmech.timeop,cstate) 
   cstate=cstate.tolist()
   cstate=cstate[0]
   cstate=np.array(cstate)
   #should really learn how to use ndarrays...
   
   
   
   
