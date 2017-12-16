import math

def factory(cell,charges, cellsum=0, center=(0,0,0), coulfac=0.5, ltrunc=1E-2):
   """
   Return a function f((x,y,z)) that returns the lattice potential energy
   INPUT:
      cell     = cell size as a tuple (Lx,Ly,Lz)
      charges  = the charge of the ions. It is a tuple containing the 
                 charge of each site. Ordering is: (0,0,0) and (Lx/2,Ly/2,Lz/2)
		 
      cellsum  = how many periodic images to consider in the summation for vpot

      center   = the center of the cell, in fractional coordinates.
      
      coulfac  = the coulomb interaction strength (ie 1/4pi*eps_0)
      
      ltrunc   = the short-range cutoff		 
   """

   if len(charges) != 2:
      raise ValueError("BCC lattice has only two atoms per cell")
      
   (Lx,Ly,Lz)=map(lambda i : float(i), cell)
   (q1,q2)=map(lambda i : float(i),charges)
   cmin=-cellsum
   cmax=cellsum   
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
	       v2 = (Lx/2+zx,Ly/2+zy,Lz/2+zz)
	       
	       vdist1=map(lambda i,j: (i-j)**2, r , v1)
	       vdist2=map(lambda i,j: (i-j)**2, r , v2)
	       
	       d1 = reduce(lambda i,j : i+j, vdist1,0)
	       d2 = reduce(lambda i,j : i+j, vdist2,0)
	       
	       d1 = math.sqrt(d1)+ltrunc
	       d2 = math.sqrt(d2)+ltrunc
	       
	       pot = pot + q1/d1 + q2/d2
	       
      
      return -coulfac*pot
      
   return vpot
   
   
