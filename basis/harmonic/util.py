"""
"""

import math
import numpy as np
import basis.precompute as hp

def integral(par,n1,n2,func,borders, nint=1000):
   """
   IN:  par = the parameter for this basis set
        n = order of the basis function
        func = the function appearing in the integrand
	borders = (x0,x1) tuple for the box borders
	nint = number of points for the numerical integral
	
   Computes <phi_n|f|phi_n> and returns it
   
   """

   dx = float(borders[1]-borders[0])/(nint-1)
   ipoints = map(lambda i : borders[0] + i*dx, range(0,nint))
   
   integral=0.0
   
   for xlo in ipoints:
      cx = xlo+dx/2
      fval = func(cx)
      phi1= base_phi(par,n1,cx)
      if n1 != n2:
         phi2= base_phi(par,n2,cx)
      else:
         phi2=phi1

      integral = integral + fval*phi1*phi2
      
   return integral*dx

def base_phi(par,n,cx, nogauss=False):
   """
   IN:   par     = the basis set parameter
         n       = the number of the function
	 cx      = the coordinate
	 nogauss = do not add the gaussian factor
	 
   evaluate the basis function at a determinate coord
   """
   
   if n > 1:
      t=int(n)
      nfac=int(n)**0.5
      for v in reversed(range(2,n)):
         nfac = nfac*math.sqrt(v)
   else:
      nfac=1
   
   dn=2**(0.5*n)
   
   fac1= 1.0/(dn*nfac)
   fac2=1
   if nogauss is False:
      fac1 = (par/math.pi)**0.25/(dn*nfac)
      fac2= math.exp(-0.5*par*cx*cx)
      
   
   
   hx = cx*math.sqrt(par)  # hermite coordinate

   h0 = 1
   h1 = 2*hx
   h2 = 4*hx*hx-2
   
   Hsec=[h0, h1, h2]
   
   if n < 3:
      return fac1*fac2*Hsec[n]
   else:
      for nn in range(3,n+1):
         Hsec[0]=1*Hsec[1]
	 Hsec[1]=1*Hsec[2]
	 Hsec[2] = 2*hx*Hsec[1] - 2*(nn-1)*Hsec[0]
	 
      return fac1*fac2*Hsec[2]


def integral_1b(bpars, l1, l2, func, borders, nint):
   """
   Manybody integral of the form:
   \int d1 \phi_1(1) V1(1) \phi_2(1). 
   
   IN: bpars = tuple with basis set parameters
       l1,l2 = quantum numbers
       V1 = one body function expecting a tuple the same size of bpars,l1,l2
       borders = borders of integration. Since we use Monte Carlo 
                 with gaussian sampling this wont be used
       nint = number of integration points
       
   
   """
   idx=hp.find_precomputed_1b(bpars,l1,l2,func)
   if not idx is False:
      return hp.precomputed_integral_1b(l1,l2,idx)
      
   
   val=float(0)
   for i in range(0,nint):
      point=[]
      phiprod=1
      for par,ll1,ll2 in zip(bpars,l1,l2):
         # for each dimension two factors phi1phi2
	 # have the same exp(-b*x*x/2) brought outside 
	 # as sampling distirbution
         s = np.random.normal(0,1.0/math.sqrt(2*par))
	 phi1 = base_phi(par,ll1,s,nogauss=True)
	 phi2 = base_phi(par,ll2,s,nogauss=True)
	 phiprod = phiprod*phi1*phi2
	 point.append(s)
	 
      veval = func(tuple(point))
      val = val+ veval*phiprod
      
   return float(val)/nint   


def integral_2b(bpars, l1,l2,l3,l4, func, borders, nint):
   """
   Same as integral_1b, only difference is that now V2 expects 
   two tuples because it is a 2d integral.
   
   The notation used here is the "chemist notation", so an element 
   marked l1,l2,l3,l4 is the integral
   \int d1d2 phi_1(1) phi_2(1) V(12) phi_3(2) phi_4(2)
		  
   """
   idx=hp.find_precomputed_2b(bpars,l1,l2,l3,l4,func)
   if not idx is False:
      return hp.precomputed_integral_2b(l1,l2,l3,l4,idx)


   val=float(0)
   for i in range(0,nint):
      point1=[]
      point2=[]
      
      phiprod=1
      for par,ll1,ll2,ll3,ll4 in zip(bpars,l1,l2,l3,l4):
         s1 = np.random.normal(0,1.0/math.sqrt(2*par))
	 s2 = np.random.normal(0,1.0/math.sqrt(2*par))
	 
	 phi1 = base_phi(par,ll1,s1,nogauss=True)
	 phi2 = base_phi(par,ll2,s1,nogauss=True)
	 phi3 = base_phi(par,ll3,s2,nogauss=True)
	 phi4 = base_phi(par,ll4,s2,nogauss=True)
	 
	 
	 phiprod = phiprod*phi1*phi2*phi3*phi4
	 point1.append(s1)
	 point2.append(s2)
	 
      veval = func(tuple(point1),tuple(point2))
      val = val+ veval*phiprod
      
   return float(val)/nint   
