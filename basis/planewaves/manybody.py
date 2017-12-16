"""
"""
import math
import cmath
from scipy.special import binom
import scipy.sparse as sparse
from aux.enum import enumerator as enum
from inspect import getargspec as nargs
import numpy as np

import basis.base as base
import utils as hu
import aux.utils as utils
from basis.counting import define_quantumnumbers

class bset(base.bset):
   def __init__(self, Lbox,ksizes, npart, symmetry=None, fftpoints=128):
      """
      Initialize a set of N-particle basis functions made of planewaves.
      IN: Lbox      = tuple containing the size of the box for each direction
          ksizes    = tuple for the number of single-particle wavevectors for each dimension. Should be an even quantity.
          npart     = number of particles
          symmetry  = can be "Bose", "Fermi" or None
	  fftpoints = points for each direction of the fourier transforms.  

      
      Single particle in pbc wavefunction phi_n (r) = exp(-i k_n * r)/sqrt(V). For n=-ksize/2,...,ksize/2: this 
      is considered in the definition of wavefunctions, where the levels are translated accordingly.
      Box is considered symmetric and centered so that the edges are at +/- L/2 for each direction.
      """

      # if ksizes not even translate accordingly.
      ksizes=map(lambda i : i + i%2 , ksizes)
      self._ksizes = ksizes
      self._enu = enum(self._ksizes)
      self._sizesing = reduce(lambda i,j : int(i)*int(j),self._ksizes,1) 
      self._kbase = map(lambda i : 2*math.pi/i, Lbox)
      self.Lbox = Lbox
      self.Vbox = 1
      for L in Lbox:
         self.Vbox = self.Vbox*L
      
      self._npart = int(npart)   
      self._symmetry = symmetry
      self.dim = len(self._ksizes)
      
      if fftpoints % 2 != 0:
         print "++++ ODD fftpoints given. Making it even"
	 fftpoints = fftpoints + 1
	 
	 
      if fftpoints/2 < 2*max(ksizes):
         print "++++ WARNING: fftpoints is too small, making it the minimum allowed size."
	 fftpoints = 4*max(ksizes)
	 
      self.fftpoints=fftpoints
      
      
      if not symmetry in ("Bose","Fermi","Distinguishable",None):
         print "Unrecognized symmetry=",symmetry
	 raise ValueError
      
      
      maxcount = int(self._sizesing**self._npart)
      
      if self._symmetry is "Fermi":
         if self._sizesing < npart:
            raise ValueError
	 
	 
	 
         maxcount = int(binom(self._sizesing,self._npart)*math.factorial(self._npart))
      
      
      self._qnumbers=define_quantumnumbers(self._ksizes,maxcount,npart,symmetry)
      self.size = len(self._qnumbers)
      self.firstder_matrix=None
      self.secder_matrix=None
      self.Umatrix=None
      self.V2matrix=None
      
      # anytime braket is called with a different function
      # you register and compute the elements for a new 
      # matrix. Elements should be (func,nbodies)
      self._regmatrices=[]
      
      # integrator functions will register fourier transforms the first time 
      # <func> is asked for. The elements of these arrays are tuples (func, ndarray)
      self._regft1b=[]
      self._regft2b=[]
      

   def integrator_1b(self,l1,l2,func):
      """
      
      """
      
      lidx = map(lambda i,j : i-j + self.fftpoints/2, l2,l1)   
      lidx=tuple(lidx)
      for (f,ft) in self._regft1b:
         if f == func:
	    return ft[lidx]   
	    
      nfft = hu.fft(func,self.fftpoints,self.Lbox)
      self._regft1b.append((func,nfft))
      
      return nfft[lidx]
       
      
      
   def integrator_2b(self,l1,l2,l3,l4,func):
      """
      This integral can be decoupled and becomes 
      2^d delta(l2+l4 -l1-l3)* (1/V)\int_{L/2}^{L/2} dt f(2|t|) exp(-i (k2+k3 - k1-k4)*t)
      """
      
      ldiff=[ l2[d]+l4[d]-l1[d]-l3[d] for d in xrange(self.dim) ]
      if max(ldiff) != 0:
         return complex(0)
      
      cnk = [ l2[d]+l3[d]-l1[d]-l4[d]+self.fftpoints/2 for d in xrange(self.dim) ]
      cnk=tuple(cnk)
      
      def fwrap(t):
         """
	 From f(x,y) wrap to f(2|t|) where t is given as a tuple
	 """
         p0 = [0]*self.dim
	 targ = map(lambda i : i*2,t)
         return func(p0,t)
	 
      
      for (f,ft) in self._regft2b:
         if f == func:
	    return 2**self.dim*ft[cnk]
	    
      nfft = hu.fft(fwrap,fftpoints,self.Lbox)
      self._regft2b.append((func,nfft))
      return 2**self.dim*nnft[cnk]
      
      
   def braket(self,i,j,func):
      return super(bset,self).braket(i,j,func,complex)


   def base_phi(self,l,x,d):
      return cmath.exp(complex(0,-self._kbase[d]*(l-self._ksizes[d]/2)*x))/math.sqrt(self.Lbox[d])


   def prepareMatrices(self,V1,V2):
      """
      Brakets are more easily computed using pre-computed basic integrals of the kind 
      <PHI_i U PHI_j> for each dimension, where i,j are two single particle 
      quantum numbers. 
      
      For derivatives: 
      first derivative, second derivative for each direction
      Format: [d1,d2,...] where d1 is a non-zero elements matrix.
      Ok: planewaves are trivial but to keep the format d1 is still a 
      matrix.
      
      V1:
      Two entry matrix with elements <i V1 j>, where 
      i,j come from the mapping between quantum numbers and index.
      
      For 2 body:
      V2 function applied to <ij V2 kl>.
      Format: matrix with four component i,j,k,l. These components are 
      from the 1-1 mapping between quantum numbers and index.
      """
      
      self.firstder_matrix=[]
      self.secder_matrix=[]

      for d,size in enumerate(self._ksizes):
         d1m=[]
	 d2m=[]
	 
	 
	 for i in xrange(0,size):
            first_nnz = [complex(0,-self._kbase[d]*(i-size/2))]
	    sec_nnz = [complex(-(self._kbase[d]*(i-size/2))**2,0)]
	    # i+1 because k_n start from n=1
	 
	    d1m.append(first_nnz)
	    d2m.append(sec_nnz)

         self.firstder_matrix.append(d1m)
	 self.secder_matrix.append(d2m)
	 
      
      
      
      print "Registering external potential..."
      self.registerFunction(V1,1,dtype=complex)
      self.Umatrix = self._regmatrices[0][0]

      if not V2 is None:
         print "Registering interaction..."
         self.registerFunction(V2,2,dtype=complex)
         self.V2matrix = self._regmatrices[1][0]

		  
      print "Base matrix elements prepared."


   def first_der_distinguishable(self, qn1, qn2, locked = False):
      """
      IN: the base numbers i,j
      
      Returns [<i|sum_N d/dx|j>,...] for each dimension. 
      
      The result for each space direction is a sum over atoms and 
      there are no offdiagonal elements.
      
      As in braket, locked will lock the evaluation of the derivatives 
      only on the first particle
      """
      
      
      diffcoord=[]
      rval= [complex(0)]*self.dim
      
      if qn1 != qn2:
         return rval
      

      for t1 in qn1:
         for d,l1 in enumerate(t1):
	    rval[d] = rval[d]+self.firstder_matrix[d][l1][0]
	    
	 if locked:
	    break	 
	
      return rval
   

   def sec_der_distinguishable(self, qn1, qn2, locked = False):
      """
      IN: the base numbers i,j
      
      Returns <i|d^2/dx^2|j>. Same as before.
      
      locked is defined as in the first derivative function
      """
      
      
      diffcoord=[]
      rval= [complex(0)]*self.dim
      
      if qn1 != qn2:
         return rval
      
      for t1 in qn1:
         for d,l1 in enumerate(t1):
	    rval[d] = rval[d]+self.secder_matrix[d][l1][0]
	    
	 if locked:
	    break	 
	
      return rval

