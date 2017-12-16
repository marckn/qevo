import basis.planewaves.manybody as bpw
import quantum.framework as q
import math
from os import mkdir
import numpy as np
from aux.utils import printData as fprint
import sys
import scipy.constants as const


import aux.cells.hcp as crystal

# hcp cell ratios: y/x = sqrt(3), z/s = 1.633
# Berillium hcp cell: x = 4.257a0
# Z = 4

lcell = 4.257
L=[lcell,lcell*3**0.5,lcell*1.633]
nks=[10,14,14]
N=1
lamb=0.5
dodens=False
b=bpw.bset(L,nks,N,"Distinguishable",64)


vper = crystal.factory(L,(4.0,4.0,4.0,4.0),0,(0.5,0.5,0.5))


Npoints=(64,64,64)
dr=(L[0]/(Npoints[0]-1),L[1]/(Npoints[1]-1),L[2]/(Npoints[2]-1))
Xs = [ (-L[0]/2 + i*dr[0]) for i in xrange(Npoints[0]) ]
Ys = [ (-L[1]/2 + i*dr[1]) for i in xrange(Npoints[1]) ]
Zs = [ (-L[2]/2 + i*dr[2]) for i in xrange(Npoints[2]) ]
try:
   mkdir("output")
except:
   pass


Uext = np.zeros(Npoints)
print "Tabulating potential [ ",
sys.stdout.flush()
knt=0
kntmax= reduce(lambda i,j : int(i)*int(j), Npoints, 1)
kntmax=kntmax/10

for i,x in enumerate(Xs):
   for j,y in enumerate(Ys):
      for k,z in enumerate(Zs):
         knt = knt + 1
	 if knt%kntmax == 0:
	    knt=0
	    print "*",
	    sys.stdout.flush()
	    
	    
         Uext[i,j,k]= vper((x,y,z))

print " ]."
fprint("output/uext",Uext,(Xs,Ys,Zs),format="csv") 
print "Uext saved in output"
   


def vext(r):
   (ix,iy,iz) = map(lambda i, j, k : int((i-j/2)/k), r, L, dr)
   return Uext[ix,iy,iz]

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








rho_points=(20,20,20)
rho_boundaries=[ (-L[0]/2,L[0]/2),
                 (-L[1]/2,L[1]/2),
		 (-L[2]/2,L[2]/2) ]
		 
dr=[]
for n,cb in zip(rho_points,rho_boundaries):
   v = float(cb[1]-cb[0])/(n-1)
   dr.append(v)

Xs=[]
Ys=[]
Zs=[]

for i in range(rho_points[0]):
   cx = rho_boundaries[0][0] + i*dr[0]
   Xs.append(cx)

for i in range(rho_points[1]):
   cy = rho_boundaries[1][0] + i*dr[1]
   Ys.append(cy)

for i in range(rho_points[2]):
   cz = rho_boundaries[2][0] + i*dr[2]
   Zs.append(cz)
   

dr = tuple(dr)


   

T=5
beta=1.0/T
minstates=10

f = open("output/states.dat","w+")
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

icut = max(icut,minstates)
icut = min(icut,len(states)+1)

lw = math.exp(-beta*(states[-1][0]-states[0][0]))
if lw > 1E-2:
   print "Warning: at this temperature, weight of highest state is ",lw
   




print "Computing densities for ", icut," energy states."   
   
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
      b.print_PSIn(qmech.U[:,ids].tolist(),0,fnm,fdir,
                   [(-L[0]/2,L[0]/2),(-L[1]/2,L[1]/2),(-L[2]/2,L[2]/2)],
		      (50,50,50),1E-2,format="csv", modquad=True)
      

      if dodens:
         fdnm="rho_"+str(idx)+".dat"
         (cdens,axes) = b.print_density(qmech.U[:,ids].tolist(),fdnm,fdir,rho_boundaries,rho_points,TOLERANCE=0.5E-1,format="csv")
	 rho_i = rho_i + cdens/len(st)

   
   if dodens:
      fprint("output/state_"+str(ii)+"/rho.dat",rho_i,(Xs,Ys,Zs),format="csv")

