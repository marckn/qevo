import basis.harmonic.manybody as hm
import quantum.framework as q
import math
import basis.precompute as hp
from os import mkdir
import numpy
from aux.utils import printData as fprint
import sys
import scipy.constants as const

###########################
###########################
######### INPUT ###########

npart=2
symmetry = "Fermi"
interaction = "On"
nlevs=(4,4)
eta = 1.38
Rw = 1.3
gamma_scr = 0.9
epsilon_r = 12.5
effective_mass = 0.07
T=11.6045

rho_points=(100,50)

# in nm
rho_boundaries=((-40,40),
                (-20,20))


# energy vs temperature
Tf = 100
Nt = 200
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


###########################
## TEST UNITS
#hbwx=4.23
#hbwy=5.84

#eta = hbwy/hbwx
#hbw0 = ((hbwx**2 + hbwy**2)/2)**0.5

#hbw0is = hbw0 * const.e / 1000
#Rsis = gamma_scr*const.e**2 * (effective_mass*const.m_e/(const.hbar**2*hbw0is))**0.5
#Rsis = Rsis/(4*math.pi*const.epsilon_0*epsilon_r)

#print "ETA ", eta, "hbarw0 = ", hbw0, " Rw = ", Rsis

#print "\n\n\n"
###########################



dt = float(Tf)/(Nt-1)
Temps=[t*dt+dt/2 for t in range(Nt)] 

hbarw=5.1
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

hp.loadFunc(namex,v1x,fdir="../../data/harmonic")
hp.loadFunc(namey,v1y,fdir="../../data/harmonic")
hp.compoundFunc("Uext",v1,[namex,namey],[wx**2,wy**2])
hp.loadFunc(namer,v2,Rw,fdir="../../data/harmonic")

basis=hm.bset(bpar,nlevs,npart,[(-10,10),(-10,10)],10000,symmetry)


lamb = 0.5
if interaction == "On":
   basis.prepareMatrices(v1,v2)
   qmech = q.qfr(basis, v1, lamb,v2)
else:
   basis.prepareMatrices(v1,v2)
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


"""
################################
##### Perturbation theory ######
coulfac=Rw
E0 = ens[0]
def psibrak(ni,nj):
   bret = 0
   for i,ci in enumerate(qmech.U[:,ni].tolist()):
      if abs(ci)<1E-2:
         continue
      
      for j,cj in enumerate(qmech.U[:,nj].tolist()):
         if abs(cj)<1E-2:
	    continue
	 
	 bret = bret + ci*cj*basis.braket(i,j,v2)/Rw
	 

   return bret


F = numpy.zeros(shape=(len(ens),len(ens)))
print "lmax = ", len(ens)
for i in range(len(ens)):
   print i,
   sys.stdout.flush()
   for j in range(i+1):
      v = psibrak(i,j)
      F[i,j]=v
      F[j,i]=v

print "\n"
E1 = F[0,0]
E2=0
E3=0
E3term=0
for k2 in range(1,len(ens)):
   den = E0 - ens[k2]
   num = F[k2,0]
   E2 = E2 + num*num/den
   E3term = E3term + num*num/den**2
   for k3 in range(1, len(ens)):
      den = (ens[0]-ens[k2])*(ens[0]-ens[k3])
      E3 = E3 + F[0,k3]*F[k3,k2]*F[k2,0]/den 
   
E3 = E3 - E3term*F[0,0]   
   
print "First order energy = ", (E0+coulfac*E1)
print "Second order energy = ",(E0+coulfac*E1+coulfac**2 * E2)
print "Third order energy = ",(E0+coulfac*E1+coulfac**2 * E2 + coulfac**3 * E3)


####################################
####################################
"""


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
   
   
# Now computing and printing the density for each of these states and 
# for the thermal state at some temperature T
# For a given eigenvalue the density considering the degeneracy of states is simply the 
# sum of the density of the single states. Just integrate the imaginary phases and see 
# that they cancel out.

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
   
   
f = open("output/states.dat","w+")
f.write("List of states:\n")

fth = open("output/edt.dat", "w+")
fth.write("#T[K]   E[meV]\n")
for cT in Temps:
   cbeta = 1./(kb*cT)
   cZ = 0
   cE = 0
   for i,s in enumerate(states):
      cw = len(s[1])*math.exp(-cbeta*(s[0]-states[0][0]))
      cZ = cZ + cw
      cE = cE + s[0]*cw
      
   cE = cE/cZ
   fth.write("{0:5.2f} {1:10.5f}\n".format(cT,cE*hbarw))

fth.close()   
      

Z=0
icut=len(states)+1
for i,s in enumerate(states):
   ln="{0:4d}) En={1:3f} hbarw  =  {2:3f} meV Degeneracy: {3:3d}\n".format(i,s[0],s[0]*hbarw,len(s[1]))
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
   
rho_T = numpy.zeros(rho_points)
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
   rho_i = numpy.zeros(rho_points)
   knt=0
   kmax = len(st)*rho_points[1]*rho_points[0]
   dec = kmax/10
   print "[",
   sys.stdout.flush()
   
   for ids in st:
      for i in range(rho_points[0]):
         for j in range(rho_points[1]):
            cx=Xs[i]/length_unit
	    cy=Ys[j]/length_unit
	    rho_i[i,j] = rho_i[i,j] + basis.density((cx,cy),qmech.U[:,ids].tolist())/len(st)
	    if (knt+1)%dec == 0:
	       print "*",
	       sys.stdout.flush()
	       
	    knt = knt+1
   
   print "]"	    
   rho_T = rho_T + len(st)*cn*rho_i/Z
   fprint("output/state_"+str(ii)+"/rho.dat",rho_i,(Xs,Ys))

fprint("output/thermal_rho.dat",rho_T,(Xs,Ys))   
   
         
      
