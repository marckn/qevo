import math

def factory(cell,charges, cellsum=0, center=(0,0,0), coulfac=1.0,ltrunc=1E-2):
   """
   Return a function f((x,y,z)) that returns the lattice potential energy
   INPUT:
      cell     = cell size as a tuple (Lx,Ly,Lz)
      charges  = the charge of the ions. It is a tuple containing the 
                 charge of each site. Ordering is: 
		 (0,0,0) 
		 (0.5,0.5,0)
		 (0.5,1/6,0.5)
		 (0,2/3,0.5)
		 
      cellsum  = how many periodic images to consider in the summation for vpot
      
      center   = the center of the cell, in fractional coordinates.
      
      coulfac  = the coulomb interaction strength (ie 1/4pi*eps_0)
      
      ltrunc   = the short-range cutoff
		 
   """

   if len(charges) != 4:
      raise ValueError("HCP lattice has only four atoms per cell")
      
   (Lx,Ly,Lz)=map(lambda i : float(i), cell)
   (q1,q2,q3,q4)=map(lambda i : float(i),charges)
   cmin=-cellsum
   cmax=cellsum
   Ltot=((2*cellsum+1)*Lx,(2*cellsum+1)*Ly,(2*cellsum+1)*Lz)   
   def vpot(extr):
      pot=0
      r = map(lambda i,j,k: i + j*k, extr, center, cell)
      	 
      for ccx in xrange(cmin,cmax+1):
         for ccy in xrange(cmin,cmax+1):
	    for ccz in xrange(cmin,cmax+1):
	       zx = Lx*ccx
	       zy = Ly*ccy
	       zz = Lz*ccz
	       
	       v1 = (zx,zy,zz)
	       v2 = (Lx/2+zx,Ly/2+zy,zz)
	       v3 = (Lx/2+zx,Ly/6+zy,Lz/2+zz)
	       v4 = (zx,2*Ly/2+zy,Lz/2+zz)
	       
	       
	       vdist1=[]
	       vdist2=[]
	       vdist3=[]
	       vdist4=[]
	       
	       for ri,vi,lti in zip(r,v1,Ltot):
	          pdst = abs(ri-vi)
		  if pdst > lti/2:
		     pdst = pdst - lti
		     
		  vdist1.append(pdst**2)
               
	       for ri,vi,lti in zip(r,v2,Ltot):
	          pdst = abs(ri-vi)
		  if pdst > lti/2:
		     pdst = pdst - lti
		     
		  vdist2.append(pdst**2)
				  
	       
	       for ri,vi,lti in zip(r,v3,Ltot):
	          pdst = abs(ri-vi)
		  if pdst > lti/2:
		     pdst = pdst - lti
		     
		  vdist3.append(pdst**2)
		
	       
	       for ri,vi,lti in zip(r,v4,Ltot):
	          pdst = abs(ri-vi)
		  if pdst > lti/2:
		     pdst = pdst - lti
		     
		  vdist4.append(pdst**2)
	       
	       
	       d1 = reduce(lambda i,j : i+j, vdist1,0)
	       d2 = reduce(lambda i,j : i+j, vdist2,0)
	       d3 = reduce(lambda i,j : i+j, vdist3,0)
	       d4 = reduce(lambda i,j : i+j, vdist4,0)
	       
	       
	       d1 = math.sqrt(d1)+ltrunc
	       d2 = math.sqrt(d2)+ltrunc
	       d3 = math.sqrt(d3)+ltrunc
	       d4 = math.sqrt(d4)+ltrunc
	       
	       
	       pot = pot + q1/d1 + q2/d2 + q3/d3 + q4/d4
	       
      
      return -coulfac*pot
      
   return vpot
   
   
