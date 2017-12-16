"""
"""

import math
import numpy as np
from itertools import permutations
from enum import enumerator as enum

def ftw(w, beta, func):
   """
   IN: 
      beta = 1/k_bT
      func = a tuple containing t, Re[func], Img[func]. 
             Note that this also fixes the number of integration points.
	     
      It is expected that Re[func(x)]=Re[func(-x)] and Img[func(x)]=-Img[func(-x)]
      so that ftw returns a real number. (units of 1/hbar)
   
   """

   t = func[0]
   fre = func[1]
   fim= func[2]
   
   dt = t[1] - t[0]
   
   rval=0
   
   for t,re,im in zip(t,fre,fim):
      rval = rval + re*math.cos(w*t) - im*math.sin(w*t)
   
   return 2*rval*dt/beta   
   
   

def signed_permutations(lst):
   """
   Return a list of permutations of the given list of quantum numbers.
   Repetitions are removed, so that the correct counting is for both bosons 
   and fermions.
   The format is: [(sign,(qnumbers))]
   """
   
   ltmp = list(permutations(lst))
   ltmp = list(set(ltmp))
   
   lout=[]
   
   for i in range(0, len(ltmp)):
      sign=sign_overlap(lst,ltmp[i])
      lout.append((sign,ltmp[i]))
   
   return lout	 
	 
def sign_overlap(targ1, targ2):
   """
   Return the sign of the permutation that 
   puts targ2 on targ1
   """	 
   cper=list(targ2)
   sign=0
   swapped=True
   while swapped:
      swapped=False
      for ctpl,ccorr in zip(cper,targ1):
         if ctpl == ccorr:
	    continue
	 
	 idc = cper.index(ctpl)
	 idcorr = targ1.index(ctpl)
	 
	 sign = sign + abs(idc-idcorr)
	 cper[idc],cper[idcorr]=cper[idcorr],cper[idc]
	    
	 swapped = True
	 break
      
      
   sign = sign % 2
   sign = int(-2*(sign-0.5))
   return sign



def printData(fname,data, lcoords, format="dat"):
   """
   Given a data file with coordinates list print a gnuplot .dat file
   IN:    fname   =   the output file name
          data    =   numpy ndarray containing data in the format value=V[i,j,...]
	  lcoords =   tuple of lists containing the coordinates (x1,x2,x3):
	              this are the values assigned by the enum function
		      
          format  =   the output format, supports: .dat (gnuplot-like), csv (comma separated values)
   """
   
   f = open(fname+"."+format,"w+")
   if len(data.shape) != len(lcoords):
      print "ERROR: data and coordinates are inconsistent"
      raise ValueError
      
   caps=[]
   for l in lcoords:
      caps.append(len(l))
      
   enu = enum(caps)
   
   for i in range(enu.maxidx):
      levels = enu.tolevels(i)
      
      lout=""
      for ci, cc in zip(levels,lcoords):
         if format == "dat":
	    lout= lout + "{0:4f}    ".format(float(cc[ci]))
	 elif format == "csv":
	    lout= lout + "{0:.4f},".format(float(cc[ci]))
	 
	 
      value = data[levels]
      if format == "dat":
         lout = lout + "{0:8f}\n".format(float(value))
         if levels[-1] == caps[-1]-1:
            lout = lout + "\n"
	
      elif format == "csv":
         lout = lout + "{0:.8f}\n".format(float(value))
      
       
      f.write(lout)
         

def bfromrweta(Rw, eta, hbar, coulfac, effective_mass, me):
   """
   Harmonic basis set parameter from Rw and eta parameters 
   of a quantum dot.
   """
   A= coulfac*math.sqrt(effective_mass)/math.pow(hbar,1.5)
   
   wy = A**2 / (Rw**2 * math.sqrt(1.0 + eta**2))
   wx = eta*wy
   print wx,wy
   return (wx*math.sqrt(effective_mass)*me/hbar,wy*math.sqrt(effective_mass)*me/hbar)
   
