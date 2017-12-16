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
import util as hu
import aux.utils as utils
from basis.counting import define_quantumnumbers

class bset(base.bset):
   def __init__(self, bpar,sizes, npart, borders, nint=1000, symmetry=None):
      """
      Initialize a set of N-particle basis functions.
      IN: bpar    = tuple containing the harmonic oscillator basis set 
                    parameter (m*w/hbar) for each space direction.
          sizes   = tuple for the number of single-particle states for each dimension
          npart   = number of particles
          borders = list of integration borders for each direction:
                    [(x0,x1),(y0,y1)...]
          nint    = number of points for numerical integrals
          symmetry= can be "Bose", "Fermi" or None
       
      
      """

      self._sizes = sizes
      self._enu = enum(sizes)
      self._sizesing = reduce(lambda i,j : int(i)*int(j),sizes,1) 
      self._npart = int(npart)   
      self._symmetry = symmetry
      self._bpar = bpar
      self._borders = borders
      self._nint = nint
      self.dim = len(self._sizes)
      
      
      
      if not symmetry in ("Bose","Fermi","Distinguishable",None):
         print "Unrecognized symmetry=",symmetry
	 raise ValueError
      
      
      maxcount = int(self._sizesing**self._npart)
      
      if self._symmetry is "Fermi":
         if self._sizesing < npart:
            raise ValueError
	 
	 
	 
         maxcount = int(binom(self._sizesing,self._npart)*math.factorial(self._npart))
      
      
      self._qnumbers=define_quantumnumbers(sizes,maxcount,npart,symmetry)
      self.size = len(self._qnumbers)
      self.firstder_matrix=None
      self.secder_matrix=None
      self.Umatrix=None
      self.V2matrix=None
      
      # anytime braket is called with a different function
      # you register and compute the elements for a new 
      # matrix. Elements should be (func,nbodies)
      self._regmatrices=[]

   def integrator_1b(self,l1,l2,func):
      return hu.integral_1b(self._bpar,l1,l2,func,self._borders,self._nint)
      
   def integrator_2b(self,l1,l2,l3,l4,func):
      return hu.integral_2b(self._bpar,l1,l2,l3,l4,func,self._borders,self._nint)
      
      
   def braket(self,i,j,func):
      return super(bset,self).braket(i,j,func,float)


   def base_phi(self,l,x,d):
      return hu.base_phi(self._bpar[d],l,x,nogauss=False)


   def prepareMatrices(self,V1,V2):
      """
      Brakets are more easily computed using pre-computed basic integrals of the kind 
      <PHI_i U PHI_j> for each dimension, where i,j are two single harmonic oscillator 
      quantum numbers. 
      
      For derivatives: 
      first derivative, second derivative for each direction
      Format: [d1,d2,...] where d1 is a non-zero elements matrix.
      Non zero elements only for i<=j are computed (see harmonic.py)
      
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

      for size,bp,bord in zip(self._sizes,self._bpar,self._borders):
         d1m=[]
	 d2m=[]
	 
	 def bket(i,j,func):
	    return hu.integral(bp,i,j,func,bord,self._nint)
	    
	 
	 for i in xrange(0,size):
            xintegral =bket(i,i,lambda x : x)
	    xintegral_iip = bket(i,i+1,lambda x : x)
	    xintegral_iim = bket(i,i-1,lambda x : x)
	 
	    xxintegral = bket(i,i,lambda x : x*x)
	    xxintegral_iipp = bket(i,i+2,lambda x : x*x)
	 
	    d1_ii = -bp*xintegral
	    d1_iip = math.sqrt(2*bp*(i+1)) - bp*xintegral_iip
	 
	    d2_ii = bp*(bp*xxintegral -2*math.sqrt(2*bp*i)*xintegral_iim -1.0)
	    d2_iip = -2*bp*math.sqrt(2*bp*(i+1))*xintegral
	    d2_iipp = 2*bp*math.sqrt((i+1)*(i+2))
	    d2_iipp = d2_iipp + bp*bp*xxintegral_iipp
	    d2_iipp = d2_iipp - 2*bp*math.sqrt(2*bp*(i+2))*xintegral_iip
	 
	 
	    first_nnz = [d1_ii, d1_iip]
	    sec_nnz = [d2_ii, d2_iip, d2_iipp]
	 
	    d1m.append(first_nnz)
	    d2m.append(sec_nnz)

         self.firstder_matrix.append(d1m)
	 self.secder_matrix.append(d2m)
	 
      
      
      
      self.registerFunction(V1,1,dtype=float)
      self.Umatrix = self._regmatrices[0][0]

      if not V2 is None:
         self.registerFunction(V2,2,dtype=float)
         self.V2matrix = self._regmatrices[1][0]

		  
      print "Base matrix elements prepared."


   def first_der_distinguishable(self, qn1, qn2, locked = False):
      """
      IN: the base numbers i,j
      
      Returns [<i|sum_N d/dx|j>,...] for each dimension. 
      
      See harmonic.py for discussion of the properties of this quantity that 
      applies to each space direction. 
      
      The result for each space direction is a sum over atoms, 
      
      
      
      Given a space direction, the result is nonzero only 
      if the quantum number are all the same or at most only 
      one differs by +/- 1 in the same direction. One can play with this.
       
      
      As in braket, locked will lock the evaluation of the derivatives 
      only on the first particle
      """
      
      
      diffcoord=[]
      rval= [0]*len(self._sizes)
      
      for n in range(0,self._npart):
         t1=qn1[n]
	 t2=qn2[n]
	 # quantum numbers for particle n
	 
	 for d in range(0,len(t1)):
	    q1=t1[d]
	    q2=t2[d]
	    
	    if abs(q1-q2)>1:
	       return rval
	       
	    if abs(q1-q2)==1:
	       if locked and n>0:
	          return rval
		  
	       coords=(n,d)
	       diffcoord.append(coords)
	       if len(diffcoord)>1:
	          return rval
		  
      if len(diffcoord) == 1:
         cn = diffcoord[0][0]
	 cdim = diffcoord[0][1]
	 
	 l1= qn1[cn][cdim]
	 l2= qn2[cn][cdim]
	 
	 sign=1
	 if l1 > l2:
	    l1,l2 = l2,l1
	    sign=-1
	    
	 diff = 1
	 rval[cdim]=self.firstder_matrix[cdim][l1][diff]*sign
	 return rval 

      # if here then we are doing <i|grad|i>
      
      for t1 in qn1:
         for d,l1 in enumerate(t1):
	    rval[d] = rval[d]+self.firstder_matrix[d][l1][0]
	    
	 if locked:
	    break	 
	
      return rval
   

   def sec_der_distinguishable(self, qn1, qn2, locked = False):
      """
      IN: the base numbers i,j
      
      Returns <i|d^2/dx^2|j>. Same as before but now the matrix element 
      is nonzero only for |l-m|<=2, and being proportional to 
      the kinetic energy it is also hermitian (and real valued 
      because of the basis-set), so <i,j> = <j,i>.
      
      We go as for first_der but now we have more offdiagonal couplings.
      
      Note: It is a bit redundant and probably a utility function shared with 
      first_der can be defined.
      
      locked is defined as in the first derivative function
      """
      
      
      diffcoord=[]
      rval= [0]*len(self._sizes)
      
      for n in range(0,self._npart):
         t1=qn1[n]
	 t2=qn2[n]
	 # quantum numbers for particle n
	 
	 for d in range(0,len(t1)):
	    q1=t1[d]
	    q2=t2[d]
	    
	    if abs(q1-q2)>2:
	       return rval
	       
	    if abs(q1-q2)>0:
	       if locked and n > 0:
	          return rval
	
	       coords=(n,d,abs(q1-q2))
	       diffcoord.append(coords)
	       if len(diffcoord)>1:
	          return rval
      
      if len(diffcoord) == 1:
         cn = diffcoord[0][0]
	 cdim = diffcoord[0][1]
	 
	 l1 = qn1[cn][cdim]
	 l2 = qn2[cn][cdim]
	 if l1 > l2:
	    l1,l2 = l2,l1
	    
	 diff = diffcoord[0][2]
	 rval[cdim]=self.secder_matrix[cdim][l1][diff]
	 return rval 

      # if here then we are doing <i|lapl|i>
      
      for t1 in qn1:
         for d,l1 in enumerate(t1):
	    rval[d] = rval[d]+self.secder_matrix[d][l1][0]
	    
	 if locked:
	    break	 
	
      return rval

