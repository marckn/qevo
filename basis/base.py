"""
"""
from aux.enum import enumerator as enum
import scipy.sparse as sparse
import numpy as np
from inspect import getargspec as nargs
import aux.utils as utils
import sys
import math

class bset(object):
   """
   Base class that implements a basis-set of N-particle *orthogonal* states.
   
   A N-particle state can be identified by N tuples of quantum numbers 
   dimensioned as the dimensionality of the basis. 
   PSI(r1,r2,...,rn)=PHI_i1(r1)PHI_i2(r2)...PHI_iN(rN)
   Then there is the symmetrization.
   For distinguishable and bosonic particles there are sizes^N different basis.
   For fermions we have to build them with no repetitions in 
   the quantum numbers, so will be defined 
   _sizesing that is the product of the max quantum number for each direction.
   
   enu is an enumerator object: stores a mapping between 
   harmonic oscillator levels (i1,i2,...id) and an 
   enumeration of the levels from 0 to max. 
   
   
   This class contains the universal operations of symmetrization, output and 
   bra-ket management. The particular choice of single particle wavefunction 
   have to be provided in a derived class. The derived class should override the following functions:
   
   * integrator_1b
   * integrator_2b
   * first_der_distinguishable
   * sec_der_distinguishable
   * base_phi
   
   See the comment in the dummy functions below for details.
      
   """
   def ___init___(self):
      """
      Basic structure for a basis set object with 
      all the dummy functions required.
      """
   
      self.size=0
      self.c=0
      self.borders=0
      self.nint=0
   
      self.secder_matrix=[]
      self.firstder_matrix=[]
      self.v2_matrix=[]
      self._enu = enum((1,1,1)) # dummy object
      self._symmetry="Distinguishable"
      self._qnumbers = [1,1] # dummy
      self.dim=2
      self._npart=1

   def integrator_1b(self,l1,l2,func):
      """
      One-body integral \int dr phi*_l1(r) phi_l2(r) f(r) where 
      f(r) is a function that takes a tuple.
      """
      return 1.0

   def integrator_2b(self,l1,l2,l3,l4,func):
      """
      Two-body integral in the chemist notation:
       \int dr1dr2 phi*_l1(r1)phi_l2(r1) phi*_l3(r2) phi_l4(r2) func(|r1-r2|) 
      where func is a function that takes two tuples and depends only on their distance
      """
      return 1.0

   def first_der_distinguishable(self, qn1, qn2, locked=False):
      """
      Given two sets of quantum numbers (translated to start from 0), 
      iqn1, iqn2 each of which characterizes 
      a basis element compute <qn1| sum_{i=1}^N d/dx_i |qn2>. 
      Return this element for each space direction in a list.
      
      if locked is true the sum has to be truncated after the first term, so that 
      basically the return value is [ <d/dx_1>, <d/dy_1>, ...]. This is 
      used when the particles are indistinguishable.
      
      """
      return 0
      
   def sec_der_distinguishable(self,qn1,qn2,locked=False):
      """
      Same as first_der_distinguishable(...) but for the second derivative. 
      Note that the return format is still a list divided by space direction.
      """
      return 0
      
   def base_phi(self,qnum,coord,d):
      """
      Return the single-particle wave function value at position coord, 
      with quantum number qnum for the space direction d:
       phi^([x,y,z])_qnum (coord).
      
      """
      return 1
      
   def densm(self, x1, x2, ket): # maybe will do later
      """
      The off-diagonal part of the density. It's very similar to the density 
      in way of computing but has yet to be implemented.
      """
      return 1      
   
   def braket(self,i,j,func,dtype=float):
      """
      Wrapper function that calls the distinguishable product directly or 
      in a (anti)symmetrized way. The cat doesn't want me to write this comment.
      """
      
      arf = nargs(func)
      bodies = int(len(arf.args))
      
      found=False
      for tf in self._regmatrices:
         if tf[1] == func:
	    found=True
	    Fmatrix=tf[0]
	    break

      if not found:
         print "+++Basis set+++ Registering new function ",func
	 self.registerFunction(func,bodies)
	 return self.braket(i,j,func)

      
      qn1 = self._qnumbers[i]
      qn2 = self._qnumbers[j]
      
      
      
      if self._symmetry == None or self._symmetry == "Distinguishable":
         return self.braket_distinguishable(qn1,qn2,Fmatrix,dtype)
	 
      
      pqn1 = utils.signed_permutations(qn1)
      pqn2 = utils.signed_permutations(qn2)
      nperm1 = len(pqn1)
      nperm2 = len(pqn2)
      
      
      if bodies == 1:
         fac = float(self._npart)/(nperm1*nperm2)**0.5
      else:
         fac = float(self._npart*(self._npart-1))/(2*(nperm1*nperm2)**0.5)
      
      
      val=dtype(0)
      
      for el1 in pqn1:
         q1 = el1[1]
	 sign1=1
	 
	 for el2 in pqn2:
	    q2 = el2[1]
	    sign2 = 1
	    if self._symmetry == "Fermi":
	       sign1 = el1[0]
	       sign2 = el2[0]
	       
	       
	    val = val + sign1*sign2*self.braket_distinguishable(q1,q2,Fmatrix,dtype,locked=True)
            
      return val*fac
      
      


   def braket_distinguishable(self, qn1, qn2, Fmatrix,dtype, locked=False):
      """
      
      Wrapper to the numerical integral function.
      Computes <i|func|j> where now i> and j> are 
      manybody wavefunctions basis whose quantum numbers are qn1 and qn2
      
      
      Scan _regmatrices to find the tuple containing the same 
      function passed here. If not found build the elements and 
      register them in _regmatrices.
      
      If it's a one-body product then this function computes 
      <i|sum^N f_n | j>. If it is a two-body the expectation value 
      becomes <i| sum_{l<n} f_{ln} | j> where as usual 
      |i> and |j> are basis manybody vectors
      
      
      locked evaluates the function only on r1 (or r1,r2 for 2-body), 
      this is meant to be used by the braket wrapper. 
      When evaluating in "locked" mode two things change:
      First, different quantum numbers for n > 1 or n> 2 (for 2-body) 
      are always zero because the operator couples only atom 
      variables for n = 1 or n = 2 (for 2-body).
      Second, there is no sum over atoms.
      """
      
      
      bodies = int(len(Fmatrix.shape))/2
      diffcoord=[]
      
      for n in range(0,self._npart):
         t1=qn1[n]
	 t2=qn2[n]
	 # quantum numbers for particle n
	 
	 
	 if t1 != t2:
	    if locked and n > bodies-1:
	       return 0
	       
	    diffcoord.append(n)
	    if len(diffcoord)>bodies:
	       return 0
	

      retval= dtype(0)
      
      if bodies == 1:
         if len(diffcoord) == 0:
	    for t1,t2 in zip(qn1,qn2):
	       id1 = self._enu.toidx(t1)
	       id2 = self._enu.toidx(t2)
	       retval = retval + Fmatrix[id1,id2]
	       if locked:
	          break
		       
	 else:
	    cn = diffcoord[0]
	    id1 = self._enu.toidx(qn1[cn])
	    id2 = self._enu.toidx(qn2[cn])
	    retval = Fmatrix[id1,id2]
	    
      else:       
         if len(diffcoord) == 0:
	    for n1 in range(1,self._npart):
	       idn1_i = self._enu.toidx(qn1[n1])
	       idn1_j = self._enu.toidx(qn2[n1])
	       for n2 in range(0,n1):
	          idn2_i = self._enu.toidx(qn1[n2])
	          idn2_j = self._enu.toidx(qn2[n2])
	          retval = retval + Fmatrix[idn1_i,idn1_j,idn2_i,idn2_j]
		  if locked:
		     return retval
		     
         elif len(diffcoord) == 1:
	    cn = diffcoord[0]
	    nnz = range(0,self._npart)
	    nnz = nnz[0:cn] + nnz[cn+1:]
	    idn1_i = self._enu.toidx(qn1[cn])
	    idn1_j = self._enu.toidx(qn2[cn])
	    for n2 in nnz:
	       idn2_i = self._enu.toidx(qn1[n2])
	       idn2_j = self._enu.toidx(qn2[n2])
	       retval = retval + Fmatrix[idn1_i,idn1_j,idn2_i,idn2_j]
	       if locked:
	          break
	 
	 else:
	    cn1 = diffcoord[0]
	    cn2 = diffcoord[1]
	    idn1_i = self._enu.toidx(qn1[cn1])
	    idn1_j = self._enu.toidx(qn2[cn1])
	    idn2_i = self._enu.toidx(qn1[cn2])
	    idn2_j = self._enu.toidx(qn2[cn2])
	    retval = Fmatrix[idn1_i,idn1_j,idn2_i,idn2_j]
	       
         
      
      
      return retval
      

   def registerFunction(self, func, bodies, dtype=complex):
      """
      Add function func to the list of registered functions. 
      That means that self_regmatrices is updated with the tuple (M,func) 
      where M are the matrix elements of the function.
      This is similar to preparematrices function, and the 
      format of M is defined in the same way of V1 and V2. 
       
      
      """
      if bodies > 2:
         print "No more than 2-body operators"
	 raise ValueError
	 
      sz = self._enu.maxidx
      
      if bodies == 1:
	 idt = int(sz*(sz+1)/20)
	 knt=0
	 M = sparse.lil_matrix((sz,sz), dtype=dtype)
	 print "[",
	 sys.stdout.flush()
	 for i in xrange(0,self._enu.maxidx):
            l1 = self._enu.tolevels(i)
            for j in xrange(0,i+1):
	       l2 = self._enu.tolevels(j)
	       Mele = self.integrator_1b(l1,l2,func)
	       knt = knt + 1
	       if knt % idt == 0:
	          knt = 0
		  print "*",
		  sys.stdout.flush()
		  
	       if abs(Mele) > 1E-5:
	          M[i,j]=Mele
		  Mele2=Mele
		  if isinstance(Mele,complex):
		     Mele2 = complex(Mele.real, -Mele.imag)
		     
		  M[j,i]=Mele2
	
         self._regmatrices.append((M,func))
         print "]."
      else:
         M2 = np.zeros(shape=(sz,sz,sz,sz), dtype=dtype)
         for i in xrange(0,self._enu.maxidx):
            l1 = self._enu.tolevels(i)
            for j in xrange(0,i+1):
	       l2 = self._enu.tolevels(j)
	       for k in xrange(0,i+1):
	          l3 = self._enu.tolevels(k)
	          for l in xrange(0,k+1):
	             l4 = self._enu.tolevels(l)
	          
	             # chemist notation: 
		     # \int d1d2 phi_1*(1) phi_2(1) V(12) phi_3*(2) phi_4(2)
		     M2ele = self.integrator_2b(l1,l2,l3,l4,func)
	             if abs(M2ele) < 1E-5:
		        continue
			
		     M2ele_1 = M2ele
		     M2ele_2 = M2ele
		     if isinstance(M2ele,complex):
		        M2ele_1 = self.integrator_2b(l2,l1,l3,l4,func)
			M2ele_2 = self.integrator_2b(l1,l2,l4,l3,func)
		     
		     
		     M2ele=complex(M2ele)
		     M2eleconj = complex(M2ele.real,-M2ele.imag)
		     
		     
		     M2[i,j,k,l]=M2ele
		     M2[j,i,k,l]=M2ele_1
		     M2[j,i,l,k]=M2eleconj
		     M2[i,j,l,k]=M2ele_2
		     M2[k,l,i,j]=M2ele
		     M2[k,l,j,i]=M2ele_1
		     M2[l,k,j,i]=M2eleconj
		     M2[l,k,i,j]=M2ele_2

         self._regmatrices.append((M2,func))


   def _deriv(self,i,j, targ):
      """
      Wrapper to the respective distinguishable function that takes into account 
      quantum symmetry.
      """
      
      qn1 = self._qnumbers[i]
      qn2 = self._qnumbers[j]
      
      
      if self._symmetry == None or self._symmetry == "Distinguishable":
         return targ(qn1,qn2)
	 
      
      pqn1 = utils.signed_permutations(qn1)
      pqn2 = utils.signed_permutations(qn2)
      nperm = (len(pqn1)**0.5)*len(pqn2)**0.5
      
      fac = float(self._npart)/nperm
      val=[0]*len(self._sizes)
      
      for el1 in pqn1:
         q1 = el1[1]
	 sign1=1
	 
	 for el2 in pqn2:
	    q2 = el2[1]
	    sign2 = 1
	    if self._symmetry == "Fermi":
	       sign1 = el1[0]
	       sign2 = el2[0]
	       
	    nval= targ(q1,q2,locked=True)
	    val=map(lambda i,j : i + sign1*sign2*j, val,nval)
	    
	    
      return map(lambda i : i*fac, val)


   def first_der(self, i, j):
      """
      First derivative <i|susm_n^N d/dx_n |j>
      """
      return self._deriv(i,j,self.first_der_distinguishable)
      
   def sec_der(self, i, j):
      """
      Second derivative <i|susm_n^N d^2/dx_n^2 |j>
      """
      return self._deriv(i,j,self.sec_der_distinguishable)



   def basis_element(self, coord, n):
      """
      Evaluate the basis element Psi_n (coord), where coord is a tuple of d-dim coordinates.
      """
      
      retval=complex(0)
      
      plevels = [(1,self._qnumbers[n])]
      if not self._symmetry in (None,"Distinguishable"):
         plevels = utils.signed_permutations(self._qnumbers[n])
      
      for tp in plevels:
         levels = tp[1]
	 sign=1
	 if self._symmetry == "Fermi":
	    sign=tp[0]
	    
	 val = 1
         for xn, ln in zip(coord, levels):
	    for d in xrange(self.dim):
	       x,l=xn[d],ln[d]
	       val = val*self.base_phi(l,x,d)
	 
	 retval = retval + sign*val   
	    
      return retval
      

   def PSI(self, coord, ket):
      """
      Given a vector of coefficients and a N-body coordinate compute PSI(coord) = sum_n c_n Psi_n(coord)
      
      IN: ket   = a list of coefficients defining a ket in this vector space
          coord = tuple of manybody coordinates, each consisting of d-dimensional tuples.
	  
      OUT: complex valued number
      """
      
      retval=complex(0)
      for n,cn in enumerate(ket):
         retval = retval + cn*self.basis_element(coord,n)
	 
      return retval

   def PSIn(self,x,ket,natm,TOLERANCE):
      """
      Return the natm-particle wavefunction, ignoring everything else in the system.
      """
      retval=complex(0)
      for n,cn in enumerate(ket):
         if abs(cn) < TOLERANCE:
	    continue
	       
	 qn=self._qnumbers[n]
	 qni = qn[natm]
	    
	 psin=complex(1,0)
	 for d,lev in enumerate(qni):
	    psin = psin*self.base_phi(lev,x[d],d)
	    
         retval = retval + cn*psin
         
      return retval


   def density(self, x, ket, TOLERANCE=1E-1):
      """
      Probability density defined as rho(x)= int dx2...dxN |PSI(x,x2,...,xN)|^2
      Making use of orthogonality of single particle wavefunctions we are in luck as this is not
      even an integral. It's "just" a combinatorics exercise.
      It's easier to separate the integral with diagonal and offdiagonal contributions:
      rho(x) = \int sum_n c_n^2 psi_n^2 + 2\int sum_n1<n2 c_n1*c_n2*PSI_n1*PSI_n2
      
      For distinguishable particles: the only non-zero elements in the second term are those that 
      differe by only n_1 as that wavefunction is not integrated.
      
      For indistinguishable particles: the first term is nonzero only when the same permutation applies both 
      on the two PSI. This gives N terms with sign +1 regardless of the quantum symmetry. The second term applies only to 
      n1 and n2 differing PRECISELY by ONE single-atom number. 
      
      Diagonal and offdiagonal contributions are defined in functions below.
      """
      retval=complex(0)

      for n,cn in enumerate(ket):
         if abs(cn) < TOLERANCE:
	    continue
	 
	 cn=complex(cn)
	 cnc=complex(cn.real,-cn.imag)
	    
	 
	 
	 retval = retval + cn*cnc*self._diagonal_rho(x,n)
	    
	 for n2,cn2 in enumerate(ket[0:n]):
	    if abs(cn2)<TOLERANCE:
	       continue
	    
	    cn2=complex(cn2)
	    cn2=complex(cn2.real,-cn2.imag)
	       
	    t1 = cn*cn2*self._offdiagonal_rho(x,n,n2)
	    t1conj = complex(t1.real,-t1.imag)
	    retval = retval + t1 + t1conj
       
      if abs(retval.imag) > TOLERANCE:
         print "WARNING: non-neglgible (>",TOLERANCE,") imaginary part in density"
	 
      return retval.real  




   def _diagonal_rho(self,x,n):
      """
      Private function: computes the diagonal terms in the probability density.
      """
      levels = self._qnumbers[n]
      
      val=1
	 
      if self._symmetry in (None,"Distinguishable"):
	 for d in xrange(self.dim):
	    cx,cl=x[d],levels[0][d]
	    phi=self.base_phi(cl,cx,d)
	    phi=complex(phi)
	    val = val*(phi.real**2 + phi.imag**2)
	    
	 return val
	 
      sterms=0
      for lev in levels:
         val=1
	 for d in xrange(self.dim):
	    cx,cl=x[d],lev[d] 
            phi = self.base_phi(cl,cx,d)
	    phi=complex(phi)
	    val = val*(phi.real**2 + phi.imag**2)
	    
	 sterms = sterms + val

      	 
      plev = utils.signed_permutations(levels)
      norm = len(plev)  # this is laziness
      return float(sterms)/norm

   def _offdiagonal_rho(self,x,n,n2):
      """
      Private function: computes the offdiagonal terms in the probability density.
      """
      l1 = self._qnumbers[n]
      l2 = self._qnumbers[n2]
   
      if self._symmetry in (None,"Distinguishable"):
         for el1,el2 in zip(l1[1:],l2[1:]):
            if el1 != el2:
	       return 0
   
         val=complex(1,0)
         for d in xrange(self.dim):
	    cx,cl1,cl2=x[d],l1[0][d],l2[0][d]
	    phi1 = self.base_phi(cl1,cx,d)
	    phi2 = self.base_phi(cl2,cx,d)
	    phi1=complex(phi1)
	    phi2=complex(phi2)
	    phi2=complex(phi2.real,-phi2.imag)
	    val = val*phi1*phi2
	    
	 
         return val
	 
   
      cl1 = 1*l1
      cl2 = 1*l2
      
      for el in l1:
         if el in cl2:
	    cl2.remove(el)
	    cl1.remove(el)
	    
      if len(cl1) > 1:
         return 0
      elif len(cl1) > 0:
         diff1=[l1.index(cl1[0])]
         diff2=[l2.index(cl2[0])]	
      else:
         diff1=[]
	 diff2=[]
	 	
      p1 = utils.signed_permutations(l1)
      p2 = utils.signed_permutations(l2)
      norm1 = len(p1)**0.5
      norm2 = len(p2)**0.5
   
      # I feel really ashamed of myself
      
      if len(diff1) == 0:
         #they are just exchanged quantum numbers, proceed as the diagonal case
         # btw: this should not happen if they are properly symmetrized.
	 raise AssertionError("Shouldnt be here if basis set is properly symmetrized")
	 
   
      
      idx1=diff1[0]
      idx2=diff2[0]
      
      tmpl1 = 1*l1
      tmpl2 = 1*l2
      
      tmpl1[0],tmpl1[idx1] = tmpl1[idx1],tmpl1[0]
      tmpl2[0],tmpl2[idx2] = tmpl2[idx2],tmpl2[0]
      
      # tmpl1(2) are now with the different element at the beginning
      sign=1
      if self._symmetry == "Fermi":
         sign = int((0.5-idx1%2)*2)*int((0.5-idx2%2)*2)
         sign = sign*utils.sign_overlap(tmpl1[1:], tmpl2[1:])
      
      val=complex(1,0)
      for d in xrange(self.dim):
         cx,cl1,cl2=x[d],tmpl1[0][d],tmpl2[0][d]
         phi1 = self.base_phi(cl1,cx,d)
         phi2 = self.base_phi(cl2,cx,d)
	 phi1=complex(phi1)
	 phi2=complex(phi2)
	 phi2=complex(phi2.real,-phi2.imag)
	 val = val*phi1*phi2
	 
      return sign*val/(norm1*norm2)
   
   def print_density(self, ket, fname, path=".", borders = None, npoints = None,TOLERANCE=1E-1,format="dat",modquad=False):
      """
      Wrapper for print_func, this prints the diagonal density matrix.
      """
      return self.print_func(ket,fname,path,borders,npoints,self.density,TOLERANCE,dtype=float,format=format)

   def print_PSIn(self, ket, natm, fname, path="", borders = None, npoints = None,TOLERANCE=1E-1,format="dat",modquad=True):
      """
      Print the single-particle wavefunction. Note that this is not always making sense.
      """
      
      def fpsin(x,ket,TOLERANCE):
         return self.PSIn(x,ket,natm,TOLERANCE)
      
      return self.print_func(ket,fname,path,borders,npoints,fpsin,TOLERANCE,dtype=complex,format=format,modquad=modquad)
   

   def print_func(self, ket, fname, path, borders, npoints, func,TOLERANCE,dtype=complex,format="dat",modquad=False):
      """
      Print a given function (now it's either density or single particle wavefunction) 
      of the given ket to a file.   
      IN:   ket     =   the list of coefficients of the vector
            borders =   list of touble containing the table borders
	    npoints =   for each direction the points of the grid
	 
      OUT: prints to file in X Y ... Z   VALUE   format.
      """  
      
      print "Printing plot [",
      sys.stdout.flush()
   
      if borders is None:
         borders = self._borders
      
      if npoints is None:
         npoints = tuple([100]*len(self._borders))
      
      dx = []
      for cb, cnp in zip(borders,npoints):
         cdx = float(cb[1]-cb[0])/(cnp-1)
         dx.append(cdx)
      
      dx = tuple(dx)
   
   
        
      def increment(_knt):
         knt = 1*_knt
         knt[0] = knt[0]+1
      
         for i in range(0,len(npoints)-1):
            if knt[i] >= npoints[i]:
	       knt[i] = 0
	       knt[i+1]= knt[i+1] + 1
	    else:
	       return knt
	       
	 return knt

   
      f = open(path+"/"+fname+"."+format, "w+")
   
      nlines = reduce(lambda i,j : i*j, npoints, 1)
      knt = [0]*len(dx)
     
      ires = nlines/10
      if ires == 0:
         ires = 1
	 
      densvals= np.zeros(npoints,dtype)
      Xarr = [ [ borders[d][0] + dx[d]*i for i in xrange(nd) ] for d,nd in enumerate(npoints) ]
      for i in range(0,nlines):
      
         fln=""
	 X = [borders[j][0]+knt[j]*dx[j] for j in range(0,len(dx))]
      
         val = func(X,ket,TOLERANCE)
	 densvals[tuple(knt)]=val
	 
	 if format == "dat":
	    entry = reduce(lambda ii,jj : ii+" {0:10f} ".format(jj), X, "")
	    if dtype == complex:
	       entry = entry +"{0:15f}  {1:15f}".format(val.real,val.imag)
	       if modquad:
	          entry = entry +"  {0:15f}".format(math.sqrt(val.real**2+val.imag**2))
	    else:
               entry = entry + "{0:15f}".format(val)
         
         elif format == "csv":
	    entry = reduce(lambda ii,jj : ii+"{0:.4f},".format(jj), X, "")
	    if dtype == complex:
	       entry = entry +"{0:.8f},{1:.8f}".format(val.real,val.imag)
	       if modquad:
	          entry = entry +",{0:.8f}".format(val.real**2+val.imag**2)
		     
	    else:
               entry = entry + "{0:.8f}".format(val)
         
	 
	 
	 f.write(entry+"\n")
      
	 if knt[0] == npoints[0]-1:
            f.write("\n")
	    
	    
	    
	 knt=increment(knt)
	 
	 if i%ires == 0:
	    print "*",
	    sys.stdout.flush()
	 
      
      print "] Done."
      f.close()
      
      return (densvals,Xarr)
