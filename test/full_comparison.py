import basis.harmonic_analyticals as ha
import basis.harmonic_util as hu
def v2((x1,y1),(x2,y2)):
   vl = (x1-x2)**2 + (y1-y2)**2
   return 1.0/vl**0.5

def compare(l1,l2,l3,l4,b,nbl=10,mcsteps=10000):
   van = ha.c2d(l1,l2,l3,l4,b)
   vnum = 0
   vnumq = 0
   for i in range(nbl):
      v=hu.integral_2b(b,l1,l2,l3,l4,v2,None,mcsteps)
      vnum=vnum+v
      vnumq=vnumq+v*v
      
   vnum = vnum/nbl
   vnumq=vnumq/nbl
   err = ((vnumq-vnum**2)/(nbl-1))**0.5
   print van,"  VS  ",vnum," +/- ",err
   
   return (van,vnum,err)

l1=(0,0)
l2=(0,0)
l3=(0,0)
l4=(0,0)
b=(1.0,1.0)

caps=(2,2)
from aux.enum import enumerator
enu = enumerator(caps)
mx = enu.maxidx

for i in range(mx):
   for j in range(mx):
      for k in range(mx):
         for l in range(mx):
	    l1 = enu.tolevels(i)
	    l2 = enu.tolevels(j)
	    l3 = enu.tolevels(k)
	    l4 = enu.tolevels(l)
	    print "LEVELS: ",l1,l2,l3,l4
	    rr=compare(l1,l2,l3,l4,b,nbl=15)
	    res="PASSED"
	    if abs(rr[0]-rr[1]) > 2*rr[2]:
	       res="FAILED"
	    print "   ========> ",res
