import math

def factory(cell,charges, cellsum=0, center=(0,0,0), coulfac=1.0, cutoff=1E-2):
   """
   Return a function f((x,y,z)) that returns the lattice potential energy
   INPUT:
      cell     = cell size as a tuple (Lx,Ly,Lz)
      charges  = the charge of the ions. It is a tuple containing the 
                 charge of each site. Ordering is: (0,0,0)
		 
      cellsum  = how many periodic images to consider in the summation for vpot
      
      center   = the center of the cell, in fractional coordinates.
      
      coulfac  = the coulomb interaction strength (ie 1/4pi*eps_0)
      
      cutoff   = the short-range cutoff
		 
   """

   if len(charges) != 1:
      raise ValueError("Simple cubic lattice has only one atom per cell")
      
   (Lx,Ly,Lz)=map(lambda i : float(i), cell)
   q1 = float(charges[0])
   
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
	       vdist1=map(lambda i,j: (i-j)**2, r , v1)
	       d1 = reduce(lambda i,j : i+j, vdist1,0) + 1E-10
	       d1 = math.sqrt(d1) + cutoff
	       
	       pot = pot + q1/d1
	       
      
      return -coulfac*pot
      
   return vpot
   
   
