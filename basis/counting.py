def define_quantumnumbers(sizes,size,N,symm):
   """
   Build a list containing the correspondence between 
   the quantum many-body state index and the values of 
   all the quantum numbers. By definition, quantum numbers 
   start from 0, eventual translations should be considered 
   when invoking the wavefunction.
   
   IN:  sizes = tuple of max quantum number in each direction
        N = number of particles
	symm = either Bose, Fermi or None
	
   
   OUT: list of of list of tuples: [q1,q2,q3]
        where q is the list of quantum numbers for state 1...etc etc
	and q1 = [(nx,ny..),(nx,ny,...),...].
   """

   lo=[]
   d = len(sizes)
   nset=[0]*N*d
   
   def nset_to_tuples(ns):
      """
      Reformat a nset list to a list of tuples to be used as output format
      """
      
      lou=[]
      for i in range(0,N):
         vls=[]
	 for dd in range(0,d):
	    vls.append(ns[d*i+dd])
	 
	 lou.append(tuple(vls))
   
      return lou
      
   def nset_increment(ns):
      """
      Increment the quantum numbers by 1 respecting each maximum value. 
      Do this starting from the left.
      """
      
      nns=1*ns
      nns[0]=nns[0]+1
      
      for i in range(0,N):
         for dd in range(0,d):
	    idx = i*d + dd
	    curmax=sizes[dd]
	    if nns[idx] >= curmax and idx < len(nns)-1:
	       nns[idx]=0
	       nns[idx+1]=nns[idx+1]+1
	    else:
	       return nns
	    
   def nset_overlapping(ns):
      """
      Scan through the quantum number in search of identical values. 
      Return true if found
      """
      
      for i in range(0,N):
         for j in range(0,i):
	    different=False
	    for dd in range(0,d):
	       if ns[d*i+dd] != ns[d*j+dd]:
	          different=True
		  break
	
	    if different is False:
	       return True
	      
      return False
   
   def is_permutation(l1,l2):
      """
      Check whether l1 can become l2 with permutations
      """
      tl1=list(l1)
      tl2=list(l2)
      
      s = len(tl1)
      for el in l1:
         if el in tl2:
	    tl1.remove(el)
	    tl2.remove(el)
	 else:
	    return False
      
      if len(tl1) == 0:
         return True
	 
      return False
      
      

   
   # now we have everything to make the enumeration of states
   for i in range(0,size):
      if symm == "Fermi":
         while nset_overlapping(nset):
	    nset = nset_increment(nset)

      ent = nset_to_tuples(nset)
      if symm in ("Fermi", "Bose"):
         samevec=False
	 for oldp in lo:
	    samevec = is_permutation(ent,oldp)
	    if samevec is True:
	       break
	       
	 if samevec is False: 
            lo.append(ent)
	    
      
      else:
         lo.append(ent)
      
      nset = nset_increment(nset)      

   return lo
