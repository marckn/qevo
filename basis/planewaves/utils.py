"""
"""

import numpy as np
from aux.enum import enumerator as enum
import math
import cmath

def fft(func, fftpoints, L):
   """
   Does the integral (1/L**d) \int_{-L/2}^{L/2} dr exp(-ikr)f(r) 
   in d dimensions, where k = (2*pi/L)n. 
   With some algebra one can see that this is equivalent to 
   exp(i * pi * n) \int_0^1 ds exp(-i 2*pi n*s) f(L(s-0.5)).
   This notation is intended to be vectorial in the obvious way. In this form 
   the numpy fft can be used.
   
   Returns a ndarray whose elements are wavenumbers translated by 
   nk[d]/2, so that a[n/2,n/2,..] is the k=0 component.
   """
   dim = len(L)
   pcaps = [fftpoints] * dim
   enu = enum(pcaps)
   fdata = np.zeros(pcaps)
   for samp in xrange(enu.maxidx):
      lsamp=enu.tolevels(samp)
      s = map(lambda i,j : float(i)/j, lsamp,pcaps)
      
      xsamp = [ (L[d]*float(lsamp[d])/pcaps[d] - float(L[d])/2) for d in xrange(dim) ]
      
      fdata[lsamp] = func(xsamp)
      
   ft = np.fft.fftn(fdata)
   ft = np.fft.fftshift(ft)
   
   for nks in xrange(enu.maxidx):
      ks = enu.tolevels(nks)
      nrm= (fftpoints) ** dim
      for ck in ks:
         ck = ck - fftpoints/2
	 if ck%2 != 0:
	    nrm = -nrm
	 
      ft[ks] = ft[ks]/nrm
   
   
   return ft
