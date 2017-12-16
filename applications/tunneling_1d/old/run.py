import basis.planewaves.manybody as bpw
import quantum.framework as q
import math
from os import mkdir
import numpy as np
from aux.utils import printData as fprint
import sys
import scipy.constants as const


L=[70.0,]
nks=[1024,]
N=1
lamb=0.5
dodens=True
b=bpw.bset(L,nks,N,"Distinguishable",5096)


def gauss(X,X0,SIG):
   grt=1.0
   for x,x0,sig in zip(X,X0,SIG):
      grt = grt*math.exp(-(x-x0)**2/(2*sig**2))
      
   return grt

def vext((x,)):
   a_conf=0.04
   a_barr=150.0
   sig_barr=2.0
   a_offset=20.0
   sig_offset=1.0
   
   vovr = a_conf*x**4
   vbar= a_barr*math.exp(-x**2/(2*sig_barr**2))
   voffset = a_offset*(0.5 - math.atan(sig_offset*x)/math.pi)
   
   return vovr + vbar + voffset
   


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
sigma=(0.5,)
x0=(6,)
def wrapg(x):
   return igauss(x,x0,sigma)

fpnts=5096
ft = fft(wrapg,fpnts,L)
cstate=np.zeros(b.size,dtype=complex)
cstatenorm=0
for i in xrange(b.size):
   ks = b._qnumbers[i][0]
   ks=(ks[0]-fpnts/2,)
   cstate[i]=ft[ks]
   cstatenorm= cstate[i].real**2 + cstate[i].imag**2

cstatenorm = math.sqrt(cstatenorm)
cstate=cstate/cstatenorm

rho_points=(1000,)
rho_boundaries=[ (-L[0]/2,L[0]/2),]


dt=0.0001
Nt=500
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
   
   
   
   

"""
rho_points=(50,50)
rho_boundaries=[ (-L[0]/2,L[0]/2),
                 (-L[1]/2,L[1]/2) ]
		 
dr=[]
for n,cb in zip(rho_points,rho_boundaries):
   v = float(cb[1]-cb[0])/(n-1)
   dr.append(v)

Xs=[]
Ys=[]

for i in range(rho_points[0]):
   cx = rho_boundaries[0][0] + i*dr[0]
   Xs.append(cx)

for i in range(rho_points[1]):
   cy = rho_boundaries[1][0] + i*dr[1]
   Ys.append(cy)
   

dr = tuple(dr)

try:
   mkdir("output")
except:
   pass
   

T=5
beta=1.0/T

f = open("output/states","w+")
f.write("List of states:\n")
Z=0
icut=len(states)+1
for i,s in enumerate(states):
   ln="{0:4d}) En={1:3f} Degeneracy: {2:3d}\n".format(i,s[0],len(s[1]))
   f.write(ln)
   cw = len(s[1])*math.exp(-beta*(s[0]-states[0][0]))
   if cw < 1E-3:
      icut=i
      break
   Z = Z + cw

lw = math.exp(-beta*(states[-1][0]-states[0][0]))
if lw > 1E-2:
   print "Warning: at this temperature, weight of highest state is ",lw
   




print "Computing densities for ", icut," energy states."   
   
rho_T = np.zeros(rho_points)
for ii,s in enumerate(states):
   if ii>=icut:
      break

   print "Processing energy level ",ii
   try:
      mkdir("output/state_"+str(ii))
   except:
      pass
      
   st=s[1]
   Ei=s[0]
   cn = math.exp(-beta*(Ei-states[0][0]))
   rho_i = np.zeros(rho_points)
   
   for idx,ids in enumerate(st):
      fdir = "output/state_"+str(ii)
      fnm="psi_"+str(idx)+".dat"
      b.print_PSIn(qmech.U[:,ids].tolist(),0,fnm,fdir,[(-20,20),(-20,20)],(100,100),1E-2)
      

      if dodens:
         fdnm="rho_"+str(idx)+".dat"
         (cdens,axes) = b.print_density(qmech.U[:,ids].tolist(),fdnm,fdir,rho_boundaries,rho_points,TOLERANCE=0.5E-1)
	 rho_i = rho_i + cdens/len(st)

   
   if dodens:
      rho_T = rho_T + len(st)*cn*rho_i/Z
      fprint("output/state_"+str(ii)+"/rho.dat",rho_i,(Xs,Ys))

if dodens:
   fprint("output/thermal_rho.dat",rho_T,(Xs,Ys))   

"""
