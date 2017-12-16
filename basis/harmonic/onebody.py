"""
"""

import util as hu
import basis.base as base
import math

class bset(base.bset):
   def __init__(self, bpar,size, borders, nint=1000):
      """
      IN: bpar = the only parameter for this basis set: mw/hbar
          size = the number of functions to use
	  borders = (x0,x1) tuple for the box borders
          nint = number of points in the numerical integrals
	  
      The matrix elements of the derivative and second derivative 
      operators are evaluated during this init function. 
      Considered the sparsity of the elements, the matrices are stored 
      in the format [i][j] where j spans over the non-zero elements.
      d^2/dx^2 is real symmetric and d/dx is antisymmetric (ihbar*d/dx is 
      the momentum operator and has to be hermitian) and thus only 
      the i<=j parts are computed.
      	  
      """
      self.c=bpar
      self.size=size
      self.borders=borders
      self.nint = nint
      
      self.secder_matrix = []
      self.firstder_matrix = []
      
      for i in range(0,size):
         xintegral = self.braket(i,i,lambda x : x)
	 xintegral_iip = self.braket(i,i+1,lambda x : x)
	 xintegral_iim = self.braket(i,i-1,lambda x : x)
	 
	 xxintegral = self.braket(i,i,lambda x : x*x)
	 xxintegral_iipp = self.braket(i,i+2,lambda x : x*x)
	 
	 d1_ii = -bpar*xintegral
	 d1_iip = math.sqrt(2*bpar*(i+1)) - bpar*xintegral_iip
	 
	 d2_ii = bpar*(bpar*xxintegral -2*math.sqrt(2*bpar*i)*xintegral_iim -1.0)
	 d2_iip = -2*bpar*math.sqrt(2*bpar*(i+1))*xintegral
	 d2_iipp = 2*bpar*math.sqrt((i+1)*(i+2))
	 d2_iipp = d2_iipp + bpar*bpar*xxintegral_iipp
	 d2_iipp = d2_iipp - 2*bpar*math.sqrt(2*bpar*(i+2))*xintegral_iip
	 
	 
	 first_nnz = [d1_ii, d1_iip]
	 sec_nnz = [d2_ii, d2_iip, d2_iipp]
	 
	 self.firstder_matrix.append(first_nnz)
	 self.secder_matrix.append(sec_nnz)

   
   def braket(self, i, j, func):
      """
      
      Wrapper to the numerical integral function.
      Computes <i|func|j>
      
      """
      
      
      return hu.integral(self.c,i,j,func,self.borders,self.nint)
      
   
   def first_der(self, i, j):
      """
      IN: the base numbers i,j
      
      Returns <i|d/dx|j>. The matrix element 
      is nonzero only for |i-j|<=1, and being proportional to 
      the momentum it has to be antisymmetric (and real valued 
      because of the basis-set), so <i,j> = -<j,i>. 
      That means that eigenvalues are purely imaginary, as expected 
      because p is an observable (thus with real eigenvalues) 
      and is ihbar*d/dx
      
      The elements are computed during initialization and stored 
      in a matrix whose structure is self-evident from the return 
      statement.
      
      """
      
      sign=1;
      if i > j:
         i,j = j,i
         sign=-1
	 
      diff = j - i
      
      if diff > 1:
         return 0
      
     
      return self.firstder_matrix[i][diff]*sign  
   
   
   def sec_der(self, i, j):
      """
      IN: the base numbers i,j
      
      Returns <i|d^2/dx^2|j>. The matrix element 
      is nonzero only for |i-j|<=2, and being proportional to 
      the kinetic energy it is also hermitian (and real valued 
      because of the basis-set), so <i,j> = <j,i>.
      The elements are computed during initialization and stored 
      in a matrix whose structure is self-evident from the return 
      statement.
      
      """
      
      if i > j:
         i,j = j,i

      diff = j - i
      
      if diff > 2:
         return 0
      
     
      return self.secder_matrix[i][diff]     
