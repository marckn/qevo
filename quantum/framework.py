"""
"""

import math
import cmath
import numpy as np
from numpy import linalg as la
import scipy.sparse as sparse
import scipy.sparse.linalg as las
import utils
import sys


class qfr:
   def __init__(self, bset, Vfunc, lamb, Vtb = None,dtype=complex, sparse_maxstates=200):
      """
      IN:  
          bset = the basis set object
          Vfunc = potential energy function suitable for the basis set.
	  Vtb = two-body potential function
	  lamb = hbar^2/2m in the relevant units
	  sparse_maxstates = if doing sparse eigendecomposition recover this number of eigenvalues.
      
      This init function computes the H, p and x operators, 
      diagonalizes the hamiltonian and stores the unitary 
      transformations. The ground state energy and ket are also 
      computed and saved for ease of use.
      
      Note that the matrix P is given in hbar units.
      """
      
      self.bset = bset
      self.Vfunc = Vfunc
      self.lamb = lamb
      
      #note: add support for sparse linear algebra
      # with a flag in basis sets indicating the nnz elements
      
      pp=[]
      xx=[]
      
      self.H = np.zeros((bset.size,bset.size),dtype=dtype)
      
      print "Preparing matrix: [",
      sys.stdout.flush()
      knt=int(0)
      dec = bset.size*(bset.size+1)/(20)
      sparsity=float(0)
      for i in range(0,bset.size):
         for j in range(0,i+1):
	    if (knt+1) % dec == 0:
	       print "*",
	       sys.stdout.flush()
	    
	    knt = knt + 1
	    
	    
	    d2 = bset.sec_der(i,j)
	    T=dtype(0)
	    if isinstance(d2,list):
	       for v in d2:
	          T = T - lamb*v
	    else:
	       T = -lamb*d2
	       
	    V = bset.braket(i,j,Vfunc)
	    V2=dtype(0)
	    if not Vtb is None:
	       V2 = bset.braket(i,j,Vtb)
	    
	    E=T+V+V2
	    if abs(E) > 1E-12:
	       E=complex(E)
	       Econj=complex(E.real,-E.imag)
	       self.H[i,j]=E
	       self.H[j,i]=Econj
	       
	       if i == j:
	          sparsity = sparsity + 1
	       else:
	          sparsity = sparsity + 2
	       
	    #X = bset.braket(i,j,lambda x: x)
	    #pim = bset.first_der(i,j)
	    #P = complex(0,pim)
	    
	    
      #self.P=np.matrix(pp)
      #self.X = np.matrix(xx)
      
      sparsity = sparsity/(bset.size**2)
      sparsity = 1.0 - sparsity	    
      print "]. sparsity = ",sparsity
            
      if sparsity > 0.8 and bset.size > 1000:
         print "Using sparse eigendec"
         hsparse=sparse.csr_matrix(self.H)
         self.Hdiag,self.U = las.eigsh(hsparse, sparse_maxstates,which="SM")
	 
	 
	 nidxs = self.Hdiag.argsort()
	 nhdiag = np.zeros(self.Hdiag.shape)
	 nu = np.zeros(self.U.shape)
	 for i,idx in enumerate(nidxs):
	    nhdiag[i] = self.Hdiag[idx]
	    nu[:,i]= self.U[:,idx]
	 
	 self.Hdiag=nhdiag
	 self.U=nu	 
      else:	 
	 self.Hdiag, self.U = la.eigh(self.H)
      
      self.E0 = self.Hdiag[0]
      self.PSI0 = self.U[:,0]   
      
      self.timeop = np.matrix(np.identity(bset.size)) 
      self.densm = None
      self.densm_trace=1
      self.ctime = 0
      
   def settime(self,t):
      """
      IN:
         t = the new time
	 
      This function sets the new time at which operators are computed.
      Basically, redefines the time evolution operator with a new time 
      t, to be given in 1/hbar units.
      """
      
      delements=[]
      for En in self.Hdiag:
         v = complex(0,-t*(En-self.E0))
	 delements.append(cmath.exp(v))
      
      self.timeop = np.matrix(np.diag(delements))
      self.timeop = np.dot(self.U,self.timeop)
      self.timeop = np.dot(self.timeop,self.U.conjugate().transpose())
      
      
      
   def settemperature(self,beta=None, thres=1E-4):
      """
      IN:
         beta = the new temperature given as beta=1/k_bT, can be passed None to go back to T=0
	 thres = the threshold for the last eigenvalue to be considered as zero. If above then 
	         you'll need a greater basis set
	 
      This function sets the new temperature at which operators are computed.
      If we are at T=0 then the variable densm becomes None.
      """
      
      if beta is None:
         self.densm = None
	 self.densm_trace=1
	 return

      delements=[]
      for En in self.Hdiag:
         v = -beta*(En-self.E0)
	 delements.append(math.exp(v))
      
      if delements[-1] > thres:
         print "Warning: last eigenvalue in density matrix is ", delements[-1]
	 
	 
      self.densm = np.matrix(np.diag(delements))
      self.densm=np.dot(self.U,self.densm)
      self.densm=np.dot(self.densm,self.U.conjugate().transpose())
      
      self.densm_trace = utils.format(self.densm.trace())
      
   def heisenbergO(self, O):
      """
      IN:
         the operator to be evolved at time t
	 
      returns a np.matrix with the operator O(t) in 
      the Heisenberg picture.
      
      """
      Oheis = self.timeop.conjugate().transpose()*O*self.timeop
      
      return Oheis
   
   def getOavg(self,O,istate=0):
      """
      IN:
         O is the np.matrix operator at t=0
         istate is the eigenstate index. If "thermal" is passed as value, then 
	 the thermal expval is computed.
	 
      return <istate(t)|O|istate(t)>	 
      """
      
      if istate is None:
         if self.densm is None:
	    istate=0
	 else:
	    istate="thermal"
      
      Oheis = self.heisenbergO(O)
      
      if istate != "thermal":
         ket = np.matrix(self.U[:,istate])
         bra = ket.conjugate().transpose()
         expval = bra*Oheis*ket
      else:
         if self.densm is None:
	    raise ValueError
	 
         prod = self.densm*Oheis
         expval = prod.trace()/self.densm_trace
      
      return utils.format(expval)
      
   def getOcorr(self, O, istate):
      """
      IN: 
         O is the operator at t=0
	 istate is the eigenstate index. If "thermal" is passed as value, then 
	 the thermal expval is computed.
	 
      Compute <O(0)O(t)> = <PSI|O*U^dag(t)*O*U(t)|PSI>
      """   
      
      if istate is None:
         if self.densm is None:
	    istate=0
	 else:
	    istate="thermal"
      
      Oheis = self.heisenbergO(O)
      if istate != "thermal":
         ket = np.matrix(self.U[:,istate])
         bra = ket.conjugate().transpose()
         expval = bra*O*Oheis*ket
      else:
         if self.densm is None:
	    raise ValueError
	    
	 prod = self.densm*O*Oheis
	 expval = prod.trace()/self.densm_trace
      
      return utils.format(expval)
      
   def getXavg(self, istate=None):
      """
      Wrapper for getOavg.
      """
      
      return self.getOavg(self.X,istate)
      

   def getPavg(self, istate=None):
      """
      Wrapper for getOavg.
      """
      
      return self.getOavg(self.P,istate)
   

   def getXcorr(self, istate=None):
      """
      Wrapper for getOcorr.
      """
      
      return self.getOcorr(self.X,istate)
      

   def getPcorr(self, istate=None):
      """
      Wrapper for getOcorr.
      """
      
      return self.getOcorr(self.P,istate)
