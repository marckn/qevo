"""
"""

class enumerator:
   def __init__(self,caps,translation=0):
      """
      Mapping between enumeration (i1,i2,i3,...) 
      with maximum values given by the tuple caps (N1,N2,...) to 
      the enumeration j = 0... N1*N2*N3-1.
      
      By convention every index starts from 0. The argument translation will 
      change this.
      """
      
      self.caps = caps
      self.translation=translation
      self._factors = []
      
      
      ftot = reduce(lambda i,j : int(i)*int(j), caps, 1)
      
      self.maxidx = 1*ftot
      
      for v in caps:
         ftot = int(ftot)/int(v)
	 self._factors.append(ftot)

   def toidx(self,levels):
      """
      From a tuple of quantum numbers return the corresponding index.
      """
      idx = int(0)
      for l,f in zip(levels,self._factors):
         idx = idx + int(l-self.translation)*int(f)
	 
      return idx
      

   def tolevels(self,idx):
      """
      From an index return the tuple of quantum numbers
      """	
      
      lev=[]
      cidx = 1*idx
      
      for f in self._factors:
         cl = int(cidx)/int(f)
	 cidx = cidx - cl*f
	 lev.append(cl+self.translation)
	 
      return tuple(lev) 
         
