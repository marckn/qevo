import basis.harmonic as bh
import quantum.framework as q
import aux.ft as ft
import math

borders=(-20,20)
gamma = 1.0
size=30
lamb = 0.5  #hbar^2/2m
basis = bh.bset(gamma,size,borders)
print "Basis initialized"


print basis.braket(0,0,lambda x : 1.0)

print basis.braket(1,1,lambda x : 1.0)

print basis.braket(2,2,lambda x : 1.0)

print basis.braket(3,3,lambda x : 1.0)

print basis.braket(4,4,lambda x : 1.0)

print basis.braket(5,5,lambda x : 1.0)

print "##### Second test #####"




qmech = q.qfr(basis, lambda x : 0.5*x*x + 0.1*x*x*x*x, lamb)

print "E0=",qmech.E0

print "PSI0=",qmech.PSI0

print "PSI1=",qmech.U[:,1]

print "#### Third test ####"

#qmech.settime(0.5)

beta=0.1
qmech.settemperature(beta)

fo = open("fdc.dat","w+")

ts=[]
xxtre=[]
xxtim=[]

pptre=[]
pptim=[]

Nt=1000
dt = 2*math.pi/(Nt+1)
for i in range(0,Nt):
   t = int(i)*dt
   qmech.settime(t)
   #x = qmech.getXavg()
   #p = qmech.getPavg()
   xxt = qmech.getXcorr()
   ppt = qmech.getPcorr()
   ln=str(t)+"\t"+str(xxt.real)+"\t"+str(xxt.imag)+"\t"
   ln = ln + str(ppt.real)+"\t"+str(ppt.imag) + "\n"
   fo.write(ln)
   
   ts.append(t)
   xxtre.append(xxt.real)
   xxtim.append(xxt.imag)
   pptre.append(ppt.real)
   pptim.append(ppt.imag)
   

fo.close()


print "### Ok. Fourier transforming...###"

func_xxt = (ts,xxtre,xxtim)
func_ppt = (ts,pptre,pptim)

xxw = ft.ft(0,100,100,beta,func_xxt)
ppw = ft.ft(0,100,100,beta,func_ppt)

fo = open("fdw.dat","w+")

for w,vx,vp in zip(xxw[0],xxw[1],ppw[1]):
   ln=str(w)+"\t"+str(vx)+"\t"+str(vp)+"\n"
   fo.write(ln)
   
fo.close()
