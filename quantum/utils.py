"""
"""

import numpy as np

def format(value,CUTOFF=1E-10):
   """
   Change type of value removing infinitesimal imaginary parts or 
   1x1 matrix structure.
   """
   
   if isinstance(value,np.matrix):
      if value.size == 1:
         return format(value[0,0],CUTOFF)

   if isinstance(value,complex):
      if abs(value.imag) < CUTOFF:
         value=float(value.real)
      
      return value	 
   
   return value
	       
