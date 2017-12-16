import polynome as p
import math
import numpy as np

trigint_pre={}


def build_uptoc(c2, bx,by, func):
   """
   Build the integral I_{\theta} for all the combinations up to n2 = c2 and n4 = c2.
   These values are stored in a dictionary of numpy arrays, trigint_pre.
   """
   
   #print "Trigint: building table for n2+n4 = ",c2
      
   nentry = np.zeros(shape=(c2+1,c2+2), dtype=float)
   for n2 in range(0,c2+1,2):
      for n4 in range(0,c2+1,2):
         if n2+n4 == c2:   
            nentry[n2,n4] = func(n2,n4,bx,by)
	    
   trigint_pre[str(c2)]=nentry
	    
	       
   

def ptrigint(n2,n4,bx,by, func):
   """
   Calls one of the two variant of the trigint 
   function and uses previously stored values.
   """
   
   if (n2+n4)%2 != 0:
      return 0
   
   c2 = int(n2+n4)
   if not str(c2) in trigint_pre:
      build_uptoc(c2,bx,by, func)
      
   arrvals=trigint_pre[str(c2)]
   return arrvals[n2,n4]



def trigint(n2,n4, bx, by):
   """
   The horrible trigonometric integral... 
   (comment to be written)
   """
   
   c = int(n2+n4)/2
   alpha = by
   beta = bx - by
   
   pfac = 2**(n4+1)*2*math.pi*complex(0,1)/bx**(c+0.5)
   
   
   if beta > 0:
      zre = math.sqrt(beta)/(math.sqrt(alpha+beta))
      zimg = math.sqrt(alpha)/(math.sqrt(alpha+beta))
   
      z1 = complex(zre,zimg)
      z2 = complex(-zre,-zimg)
      z3 = complex(zre,-zimg)
      z4 = complex(-zre,zimg)
   elif beta < 0:
      impm = math.sqrt(alpha) - math.sqrt(abs(beta))
      impp = math.sqrt(alpha) + math.sqrt(abs(beta))
      den = math.sqrt(alpha+beta)
      z1 = complex(0,impm/den)
      z2 = complex(0,impp/den)
      z3 = complex(0,-impm/den)
      z4 = complex(0,-impp/den)
   else:
      if n2 % 2 != 0 or n4 % 2 != 0:
         return 0
	 
      rfac=1
      tn4=1*n4
      tn2=1*n2
      for i in range(n4/2):
         rfac = rfac*(1.0 + float(tn2+1)/(tn4-1))
	 tn4 = tn4 - 2      

      for i in range(n2/2):
         rfac = rfac*(1.0 + float(tn4+1)/(tn2-1))
	 tn2 = tn2 - 2
	 
      return 2*math.pi/rfac      

  
   idef = complex(0,1)
   
   cfac = math.factorial(c)**(-1)
   
   f1 = p.polynome([-1.0,0,1.0])
   f2 = p.polynome([1.0,0.0])
   fz1 =p.polynome([1.0,-z1])
   fz2 =p.polynome([1.0,-z2])
   fz3 =p.polynome([1.0,-z3])
   fz4 =p.polynome([1.0,-z4])
   
   
   if beta >= 0:
      pres_z1 = p.polynome_product([(f1,n2),(f2,n4),(fz2,-c-0.5),(fz3,-c-0.5),(fz4,-c-0.5)])
      pres_z4 = p.polynome_product([(f1,n2),(f2,n4),(fz2,-c-0.5),(fz3,-c-0.5),(fz1,-c-0.5)])
   
      derz1 = pres_z1.nderive(c)
      derz4 = pres_z4.nderive(c)
      
      res_z1 = cfac*derz1.eval(z1)
      res_z4 = cfac*derz4.eval(z4)
      
      
      integ = pfac*(res_z1 + res_z4)
   else:
      pres_z1 = p.polynome_product([(f1,n2),(f2,n4),(fz2,-c-0.5),(fz3,-c-0.5),(fz4,-c-0.5)])
      pres_z2 = p.polynome_product([(f1,n2),(f2,n4),(fz4,-c-0.5),(fz3,-c-0.5),(fz1,-c-0.5)])
   
      derz1 = pres_z1.nderive(c)
      derz2 = pres_z2.nderive(c)
   
      res_z1 = cfac*derz1.eval(z1)
      res_z2 = cfac*derz2.eval(z2)
   
      integ = pfac*(res_z1 + res_z2)
      
   
   if abs(integ.imag) < 1E-8:
      return integ.real
   else:
      return integ
   


def numint(n2,n4,bx,by, npoints=5000):
   dx = math.pi/(npoints-1)
   
   if n2 % 2!= 0 or n4 % 2 != 0:
      return 0
   
   cm = (n2+n4)/2
   integ=0
   for i in range(npoints):
      ct = dx*(i+0.5)
      c = math.cos(ct)
      s = math.sin(ct)
      
      num = math.pow(c,n2)*math.pow(s,n4)
      den = (by + (bx-by)*c*c)**(cm+0.5)
      
      integ = integ + dx*num/den
      
   return 2*integ
      
     



#print trigint(4,0,1,1)
#print numint(4,0,1,1,10000)
