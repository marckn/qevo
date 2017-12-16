import basis.planewaves.manybody as bpw
import quantum.framework as q
import math
from os import mkdir
import numpy as np
from aux.utils import printData as fprint
import sys
import scipy.constants as const


L=[70.0,]
nks=[512,]
N=1
lamb=0.8
dodens=True
b=bpw.bset(L,nks,N,"Distinguishable",2048)


def gauss(X,X0,SIG):
   grt=1.0
   for x,x0,sig in zip(X,X0,SIG):
      grt = grt*math.exp(-(x-x0)**2/(2*sig**2))
      
   return grt

def vext((x,)):
   p_wall=33
   sig_wall= 1.0
   eps = 1E-2
   
   a_barr = 1.0
   sig_barr=4.0
   a_offset=4.0
   sig_offset=2.0
   
   if abs(x) > p_wall:
      vwall = sig_wall**12/eps
   else:
      dwall = (abs(x) - p_wall)**12
      vwall = sig_wall**12/(dwall+eps)
   
   vbar= a_barr*math.exp(-x**2/(2*sig_barr**2))
   voffset = a_offset*(0.5 - math.atan(sig_offset*x)/math.pi)
   
   return vwall + vbar + voffset
   


Npoints=(2000,)
dr=(L[0]/(Npoints[0]-1),)
Xs = [ (-L[0]/2 + i*dr[0]) for i in xrange(Npoints[0]) ]

Uext = np.zeros(Npoints)
for i,x in enumerate(Xs):
      Uext[i,]= vext((x,))

fprint("output/uext",Uext,(Xs,)) 
print "Uext saved in output"
   

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
sigma=(1.0,)
x0=(-20,)
def wrapg(x):
   return igauss(x,x0,sigma)

fpnts=nks[0]*4
ft = fft(wrapg,fpnts,L)
cstate=np.zeros(b.size,dtype=complex)
cstatenorm=0
for i in xrange(b.size):
   ks = b._qnumbers[i][0]
   ks=(ks[0] - nks[0]/2,)

   ks=(fpnts/2-ks[0],)
   cstate[i]=ft[ks]
   cstatenorm= cstatenorm + cstate[i].real**2 + cstate[i].imag**2
   
   
cstatenorm = math.sqrt(cstatenorm)
cstate=cstate/cstatenorm

rho_points=(1000,)
rho_boundaries=[ (-L[0]/2,L[0]/2),]


dt=0.002
Nt=2000
for i in xrange(Nt):
   print "Processing frame ",i
   ct = i*dt
   fout="snap_"+str(i)+".dat"
   b.print_PSIn(cstate.tolist(),0,fout,"mov",rho_boundaries,rho_points,TOLERANCE=1E-3)
   qmech.settime(ct)
   
   cstate = np.dot(qmech.timeop,cstate) 
   cstate=cstate.tolist()
   cstate=cstate[0]
   cstate=np.array(cstate)
   #should really learn how to use ndarrays...
   
   
   
   
