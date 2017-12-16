"""
This module contains the analytical formula for some 
<phi_i|f|phi_j> that can be carried out without numerical 
integration.  

The formulas so far are for:
<i|x^n|j>

[ij|r_12|kl]   (chemist notation)

Feel free to add more if you know how to compute them.
As a first step when this module is loaded the coefficients of the hermite 
polynomials up to order 20 are computed via an iterative function. 
If more are needed the analytical formulas should realize and iterate
up to the required level.
"""

import math
import coulomb.trigint as hc
from scipy.special import binom

hermite_coeff=[(1,),
               (0,2), 
	       (-2,0,4)]
	       
def add_hermite_coeff():
   new_ord = len(hermite_coeff)
   ncoeff = [-hermite_coeff[-1][1]]
   for k in range(1,new_ord+1):
      ncoeff.append(2*hermite_coeff[-1][k-1])
      
   for k in range(1,new_ord-1):
      ncoeff[k]=ncoeff[k] - (k+1)*hermite_coeff[-1][k+1]
      
   hermite_coeff.append(tuple(ncoeff))
   

for i in range(17):
   add_hermite_coeff()
   

powx_prod_data=[1,1]
def add_powx_prod_data():
   shalf = int(len(powx_prod_data))
   ns = shalf*2
   val = ns-1
   for n in range(1, shalf-1):
      val = val*int(ns-1-2*n)
      
   powx_prod_data.append(val)
   
   
   
def powx_prod(s):
   """
   Productory for powx, computes 
   2^{-s/2} * \prod_{n=0}^{s/2-1} (s-1-2n)
   When computed for an s, stores the value in powx_prod_data 
   and calls it again if needed. So we skip this heavy operation.
   
   powx_prod_data contains values of the sum for s/2 = 0, 1, ...
   """
   if s%2 != 0:
      print "ERROR: PRODUCTORY ONLY DEFINED FOR EVEN ARGUMENTS"
      raise ValueError
   
   s = int(s)
   sh = s/2
   
   if sh < len(powx_prod_data):
      return powx_prod_data[sh]*(2**(-sh))
   else:
      toadd = sh - len(powx_prod_data)+1
      for n in range(toadd):
         add_powx_prod_data()
	 
      return powx_prod_data[-1]*(2**(-sh))
   
   
   	       
def powx(m, l, k, b):
   """
   This is <l|x^m|k>, with obvious naming of variables.
   b is the basis set parameter
   """
   ordmax = len(hermite_coeff)-1
   idmax = max(k,l)
   
   if idmax > ordmax:
      for i in range(idmax-ordmax):
         add_hermite_coeff()
	 
   
   if m % 2 != 0:
      if abs(l-k) % 2 == 0:
         return 0
   else:
      if abs(l-k) % 2 != 0:
         return 0
	 
   prefac = math.pow(2,0.5*(l+k))*math.pow(b,0.5*m)
   fct1 = math.sqrt(math.factorial(l))
   fct2 = math.sqrt(math.factorial(k))
   
   prefac = prefac*fct1*fct2
   
   koff=0
   loff=0
   
   if k % 2 != 0:
      koff = 1
   
   if l % 2 != 0:
      loff = 1
   
   
   val = 0
   for i in range(koff,k+1,2):
      for j in range(loff,l+1,2):
         s = i + j + m
	 pd = hermite_coeff[k][i]*hermite_coeff[l][j]*powx_prod(s)
	 val = val + pd
	 
   return float(val)/prefac   	 
   

#### coulomb 2d ####

def gaussint(n,alpha):
   """
   Gaussian integral \int_{-inf}^{+inf} dx x^n exp(-alfa*x*x)
   """
   if n%2 != 0:
      return 0
   
   prod = powx_prod(n)*math.sqrt(math.pi)/math.pow(alpha,(n+1.0)/2)
   return prod

def c2d(i,j,k,l,b):
   """
   Return the Coulomb integral in the chemist notation, whose matrix element is
   [ij|r12|kl] = \int d1d2 psi_i(1)psi_j(1) r12^-1 psi_k(2)psi_l(2), where 
   psi_n is a 2d harmonic oscillator eigenstate.
   
   i,j,k,l are (x,y) tuples for the levels.
   """
   idlist = list(i)+list(j)+list(k)+list(l)
   idsum = reduce(lambda ii, jj: ii+jj, idlist, 0)
   prefac=math.pow(2,float(idsum)/2)*math.pi**2/(b[0]*b[1])

   for ci in idlist:
      prefac = prefac*math.sqrt(math.factorial(ci))
      
   idmax = max(i[0],i[1],j[0],j[1],k[0],k[1],l[0],l[1])
   ordmax=len(hermite_coeff)-1
   if idmax > ordmax:
      for ii in range(idxmax-ordmax):
         add_hermite_coeff()
	
   integ = 0
   
   
   for i1 in range(i[0]+1):
      for i2 in range(i[1]+1):
         for j1 in range(j[0]+1):
	    for j2 in range(j[1]+1):
	       for k1 in range(k[0]+1):
	          for k2 in range(k[1]+1):
		     for l1 in range(l[0]+1):
		        for l2 in range(l[1]+1):
			   hfac =       hermite_coeff[i[0]][i1]*hermite_coeff[i[1]][i2]
			   hfac = hfac * hermite_coeff[j[0]][j1]*hermite_coeff[j[1]][j2]
			   hfac = hfac * hermite_coeff[k[0]][k1]*hermite_coeff[k[1]][k2]
			   hfac = hfac * hermite_coeff[l[0]][l1]*hermite_coeff[l[1]][l2]
			   bcoeff = getProduct(i1+j1,k1+l1,i2+j2,k2+l2)
			   
			   for trm in bcoeff:
			      cc = trm[0]*math.pow(b[0],float(i1+j1+k1+l1)/2)*math.pow(b[1],float(i2+j2+k2+l2)/2)
			      n1 = trm[1][0]
			      n2 = trm[1][1]
			      n3 = trm[1][2]
			      n4 = trm[1][3]
			      
			      cvl = n2+n4
			      cc2 = powx_prod(cvl)*math.sqrt(math.pi)/math.pow(2,float(cvl+1.0)/2)
			      
			      fac = hfac*cc2*gaussint(n1,2*b[0])*gaussint(n3,2*b[1])*hc.ptrigint(n2,n4,b[0],b[1], func = hc.numint)
			      integ = integ + cc*fac
			    
   return integ/prefac
   
   
   
products_list={}
def getProduct(a,b,c,d):
   """
   From x1^a * x2^b * y1^a * y2^b compute 
   (tx+sx)^a * (tx-sx)^b * (ty+sy)^c * (ty-sy)^d 
   Return a list of monomes, each element consisting of a tuple 
   (a_i , [n1,n2,n3,n4] )  that is intended to mean 
   a_i * tx^n1 * sx ^ n2 * ty^n3 * sy^n4 
   
   where tx = (x1+x2)/2 and sx = (x1-x2)/2
   NOTE however that here the terms with n1 n2, n3 or n4 odd won't be 
   passed in the list because they would give zero contribution to the 
   coulomb integral.
   The parity check on n2 and n4 is already built in the respective for loops 
   by letting them span only over even terms. On n1 and n3 a condition on 
   bt and dt is that they are even: otherwise n1 and n3 would always be odd in the loops.
   """
   
   key=str(a)+str(b)+str(c)+str(d)
   if key in products_list:
      return products_list[key]
      
   
   if a >= b:
      at = b
      bt = a-b
   else:
      at=a
      bt = b-a
      
   if c >= d:
      ct = d
      dt = c-d
   else:
      ct=c
      dt = d-c
   
   # this way (t+s)(t-s) are converted to (t^2-s^2) 
   # as far as possible. Should save on iterations.
   
   clist=[]
   
   if bt%2 != 0 or dt%2 != 0: #means that n1 or n3 are always odd
      return clist
   
   for k1 in range(at+1):
      for k2 in range(0,bt+1,2):
         for k3 in range(ct+1):
	    for k4 in range(0,dt+1,2):
	       coeff = (-1)**(k1+k2+k3+k4)
	       coeff = coeff*binom(at,k1)*binom(bt,k2)
	       coeff = coeff*binom(ct,k3)*binom(dt,k4)
	       n1 = 2*(at-k1)+bt-k2
	       n2 = 2*k1+k2
	       n3 = 2*(ct-k3)+dt-k4
	       n4 = 2*k3+k4
	       if n1%2 != 0 or n3%2 != 0:
	          continue
		  
	       lel = [n1,n2,n3,n4]
	       repeat=False
	       for i,pe in enumerate(clist):
	          if lel == pe[1]:
		     repeat=True
		     clist[i][0] = clist[i][0] + coeff
		     break
		     
               
	       if repeat is False:
	          clist.append([coeff, lel])
		  
		  
		  
   clist = map(lambda i : tuple(i), clist)
   products_list[key]=clist
   return clist
	       
	       
	       
