"""
"""

import utils

def ft(wmin,wmax,nw,beta, func):
   """
   IN: 
      wmin = the initial w
      wmax = the last w
      nw = the number of points in frequency domain
      beta = 1/k_bT
      func = a tuple containing t, Re[func], Img[func]. 
             Note that this also fixes the number of integration points
      
   The fourier transform of func in unit of hbar for the case 
   of Re [func] symmetric and Img[func] antisymmetric, so that 
   there's only a real part in output. 
   
   OUT:
      The fourier transform as a tuple of lists (w, ft[func])
   """

   
   ws = range(0,nw)
   dw = float(wmax-wmin)/(nw-1)
   ws = map(lambda x : wmin + dw*x, ws) 
   ft = []
  
   for w in ws:
      cft = utils.ftw(w,beta,func)
      ft.append(cft)
  
   return (ws,ft)
   
