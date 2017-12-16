import basis.harmonic_manybody as hm
import quantum.framework as q
import math
import basis.harmonic_precomputed as hp
from os import mkdir
import numpy
from aux.utils import printData as fprint
import sys
import scipy.constants as const
import fit_state as fit

###########################
###########################
######### INPUT ###########

npart=2
symmetry = "Distinguishable"
interaction = "Off"
nlevs=(8,5)
eta = 3.0
Rw = 1.34
gamma_scr = 0.9
epsilon_r = 12.5
effective_mass = 0.07
T=40

rho_points=(50,50)

# in nm
rho_boundaries=((-20,20),
                (-20,20))


###########################
###########################
###########################

hbarw = gamma_scr**2 * const.e**4 * effective_mass * const.m_e
hbarw = hbarw / (4*math.pi*epsilon_r*const.epsilon_0*const.hbar*Rw)**2
length_unit = const.hbar/(effective_mass * const.m_e * hbarw)**0.5

#converting to meV and nm
hbarw = hbarw*1000/const.e
length_unit = length_unit * 1E9

# In unit of hbarw
kb = (1000*const.k/const.e)/hbarw
beta = 1./(kb*T) 

print "hbarw = ", hbarw, " meV"
print "length unit is ", length_unit," nm"
print ">> NOTE: output lengths are in nm"
wx = (2/(1+eta**2))**0.5
wy = eta*wx


bpar=(wx,wy)
bn=str(eta)

def v1((x,y)):
   return 0.5*(wx**2*x*x+wy**2*y*y)

def v1x((x,y)):
   return 0.5*x*x
   
def v1y((x,y)):
   return 0.5*y*y
   

def v2((x1,y1),(x2,y2)):
   dist = (x1-x2)**2 + (y1-y2)**2
   return Rw/dist**0.5
   
namex="x**2_eta"+bn
namey="y**2_eta"+bn
namer="r12**-1_eta"+bn

hp.loadFunc(namex,v1x)
hp.loadFunc(namey,v1y)
hp.compoundFunc("Uext",v1,[namex,namey],[wx**2,wy**2])

if interaction == "On":
   hp.loadFunc(namer,v2,Rw)

basis=hm.bset(bpar,nlevs,npart,[(-10,10),(-10,10)],10000,symmetry)


lamb = 0.5
if interaction == "On":
   basis.prepareMatrices(v1,v2)
   qmech = q.qfr(basis, v1, lamb,v2)
else:
   basis.prepareMatrices(v1,None)
   qmech = q.qfr(basis, v1, lamb)
   


#psi0=qmech.U[:,0].tolist()
#psi1=qmech.U[:,1].tolist()



# from eigenvalues list regroup with degenerate states

# this is an array of tuples. Each tuple is (E_n, states), where 
# states is a tuple of indices corresponding to the relative eigenstates.
# Since E_n are approximated we need to specify a tolerance for the definition of 
# degenerate states (i.e. E_m = E_n). Note that E_n are already ordered.

tolerance=1E-1
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
   ln="{0:4d}) En={1:3f} hbarw  =  {2:3f} meV Degeneracy: {3:3d}".format(i,s[0],s[0]*hbarw,len(s[1]))
   print ln
   
   
dr=[]
for np,b in zip(rho_points,rho_boundaries):
   v = float(b[1]-b[0])/(np-1)
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
   


neig = qmech.U.shape[0]
thstate=numpy.zeros(neig)

cnorm=0
for i,e in enumerate(ens):
   dE = ens[i]-ens[0]
   ce = math.exp(-beta*dE)
   cnorm = cnorm + ce*ce
   thstate = thstate + ce*qmech.U[:,i]
   
thstate = thstate/cnorm

tmax=20.0
ntm = 40
dt = tmax/(ntm-1)
fconv=1.0

times= [ i*dt*fconv for i in xrange(ntm) ]
for it,ct in enumerate(times):
   qmech.settime(ct)
   psit = numpy.dot(qmech.timeop, thstate)
   
   psit = psit.tolist()
   psit = psit[0]
   
   rho_i = numpy.zeros(rho_points)
   for i in range(rho_points[0]):
      for j in range(rho_points[1]):
         cx = Xs[i]/length_unit
	 cy = Ys[j]/length_unit
	 rho_i[i,j] = basis.density((cx,cy),psit)
   
   fprint("output/rho_"+str(it)+".dat",rho_i,(Xs,Ys))
   print it,
   sys.stdout.flush()
