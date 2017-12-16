import math
import cmath
import sys

class polynome:
   def __init__(self,lcoeff=(0,)):
      """
      The polynome class, defined by a passed argument of the kind:
      (a_n,a_n-1,...,a_0), so that the polynome is a_nx^n....a_0
      """
      
      self.coeff=lcoeff
      self.order = len(lcoeff)-1
      
      self._trim()
      
      
   def printout(self):
      """
      Print the polynome formatted as a_nx^n + ... + a0
      """
      
      for i,cn in enumerate(self.coeff):
         if abs(cn)<1E-12:
	    continue
	 
	 if isinstance(cn,complex):
	    sou = " + "
	 elif cn > 0:
	    sou = " + "
	 else:
	    sou = " "
	    
         if isinstance(cn,complex):
	    sou = sou + "({0:8f})".format(cn)
	 else:
	    sou = sou + "{0:8f}".format(cn)
	 
	 if self.order-i > 1:
	    sou=sou + "x^({0:2d})".format(int(self.order-i))
	 elif self.order-i > 0:
	    sou=sou + "x"
	 
	 print sou,
	 
      sys.stdout.flush()
	 
   
   def sum(self,pol2, isign=1):
      """
      Sum self with pol2 and return a new polynome.
      If isign -1 then make the difference self - pol2
      """
      
      if pol2.order > self.order:
         ncoeff = isign*list(pol2.coeff)
	 for i,cv in enumerate(reversed(self.coeff)):
	    ncoeff[-i-1] = ncoeff[-i-1] + cv

         return polynome(tuple(ncoeff))	     
      else:
         ncoeff = list(self.coeff)
	 for i,cv in enumerate(reversed(pol2.coeff)):
	    print i, cv
	    ncoeff[-i-1] = ncoeff[-i-1] + isign*cv

         return polynome(tuple(ncoeff))	     
       
       
   def sub(self,pol2):
      """
      self - pol2
      """       
      return self.sum(pol2,-1)
      
   def mult(self,pol2):
      """
      Returns a new polynome containing the product
      """
      
      ord2 = pol2.order
      ordmax = self.order + ord2
      
      ncoeff=[0]*(ordmax+1)
      
      for i in range(self.order+1):
         for j in range(pol2.order+1):
	    c_i = self.coeff[-i-1]
	    c_j = pol2.coeff[-j-1]
	    cord = i+j
	    ncoeff[-cord-1] = ncoeff[-cord-1] + c_i*c_j
	    
      return polynome(tuple(ncoeff))
      
   def pow(self,n):
      """
      Return the polynome p^n, with n>=0
      """
      
      if n<0:
         print "ERROR: cannot take negative powers of polynomes"
	 raise ValueError
	 
      if n == 0:
         return polynome((1,))
	 
      cpol = self.mult(self)
      for i in range(1,n-1):
         cpol = cpol.mult(self)
	 
      return cpol
      
   def eval(self, z):
      """
      Evaluate p(z), where z is [generally] a complex number
      """
      
      zret = complex(self.coeff[-1])
      for i, vls in enumerate(reversed(self.coeff[:-1])):
         exp = i+1
	 zret = zret + complex(vls*z**exp)      
      
      if abs(zret.imag) < 1E-12 and not z is complex:
         zret = float(zret.real)
      
      return zret
      
   def derive(self):
      """
      Return a new polynome with the derivative of this.
      """   
      
      ncoeff=[]
      for i,el in enumerate(self.coeff[:-1]):
         deg = self.order - i
	 der = el*deg
	 ncoeff.append(der)
	 
      return polynome(tuple(ncoeff)) 
      

   def _trim(self):
      """
      If after operations higher order coefficients become 0 you can 
      trim the polynome.
      """      
      
      ncoeff=[]
      canstart=False
      
      for e in self.coeff:
         if abs(e)>1E-16:
	    canstart=True
	    
	 if canstart:
	    ncoeff.append(e)
	    
      
      self.coeff = tuple(ncoeff)
      self.order = len(self.coeff)-1
      
      if len(self.coeff) == 0:
         self.coeff=(0,)
	 self.order=0

   def null(self):
      """
      Determine whether the polynome is the null element or not.
      """
      self._trim()
      if self.order == 0 and abs(self.coeff[0]) < 1E-16:
         return True

      else:
         return False
      
   def equal(self,pol2,TOLERANCE=1E-10):
      """
      Determine if two polynomes are equal or not. 
      If they are not equal return 0, otherwise
      return the factor so that 
      self = fac * pol2
      """
      
      if self.null() and pol2.null():
         return 1
      elif self.null() or pol2.null():
         return 0
      
      if self.order != pol2.order:
         return 0
      
      lqu=[]
      for c1,c2 in zip(self.coeff, pol2.coeff):
         if abs(c2) > 1E-12 and abs(c1) < 1E-12:
	    return 0
	 elif abs(c2) < 1E-12 and abs(c1) < 1E-12:
	    lqu.append(None)
	 else:
	    lqu.append(1.0*c2/c1)   

      lf= lqu[0]
      for ff in lqu[1:]:
         if ff is None:
	    continue
	    
	 if abs(lf-ff) > TOLERANCE:
	    return 0	    
      
      return 1.0/lf


################################################################
################################################################



class polynome_product:
   def __init__(self,lpoly, fmult=1):
      """
      lpoly is a list of tuples, each tuple consists in a 
      couple (polynome, expo), where polynome is a polynome object and expo is 
      the exponent of the factor. 
      fmult is the overall multiplicative factor
      """
      
      self.polynomes=[]
      self.fmult=fmult
      
      for tp in lpoly:
         self.polynomes.append(tp)
	 
      self.simplify()
	 
   def printout(self):
      """
      Print the current product as (polyclass_printout)^(n1)*(...)
      """
      
      if isinstance(self.fmult, complex):
         print "({0:6f})".format(self.fmult),
      elif abs(self.fmult -1) > 1E-12:
	    print "{0:6f}".format(self.fmult),
	 
      for (poly,exp) in self.polynomes:
         if exp == 0:
	    continue
	 
	       
	 print "(",
	 poly.printout()
         if exp != 1:
	    st=")^({0:2d})  ".format(int(exp))
	 else:
	    st=") "
	    
	 print st,
	 
      sys.stdout.flush()

   def eval(self,z):
      """
      Evaluate this product.
      """
      
      zprod = self.fmult*complex(1.0,0)
      zvl = complex(z)
      for (pol,expo) in self.polynomes:
         feval = pol.eval(zvl)**expo
	 zprod = zprod*feval
	 
      if abs(zprod.imag) < 1E-16 and not z is complex:
         zprod = float(zprod.real)

      return zprod
      
   def derive(self):
      """
      Make the (first) derivative of this product. 
      Return a pprod_sum object.
      """	 

      lret=[]
      for i,(fac,exp) in enumerate(self.polynomes):
	 lcur=[]
	 dp = fac.derive()
	 if dp.null():
	    continue
	    
	 lcur.append((dp,1))
	 if exp != 1:
	    lcur.append((fac,exp-1))
	    
	 lcur = lcur + self.polynomes[0:i] + self.polynomes[i+1:]
	 	 
	 lcur = polynome_product(lcur,self.fmult*exp)
	 
	 
	 lret.append(lcur)
	 
      return pprod_sum(lret)   

   def nderive(self,n):
      """
      nth order derivative. Just calls derive n times
      """
      if n < 0:
         return pprod_sum([])
      
      if n == 0:
         return self

      fder = self.derive()
      for i in range(n-1):
         fder = fder.derive()
	 
      return fder


	       
   def simplify(self):
      """
      Simplify a product by removing factors ()^0, multiplying 
      same factors (minus sign) and putting zero-th order 
      polynomes to the multiplicative factor.
      
      #note: still to do the multiplication of same factors
      """
      
      npoly=[]
      for tp in self.polynomes:
         if abs(tp[1]) < 1E-12:
	    continue
	     
	 if tp[0].order == 0:
	    self.fmult = self.fmult * tp[0].coeff[0]**tp[1]
	    
	 else:
	    npoly.append(tp)   


      nfmult = self.fmult
      
      nnpol=[]
      excl_list=[]
      for i, tp1 in enumerate(npoly):
         sexp=tp1[1]
	 sfac=1
	 if i in excl_list:
	    continue
	    
	 for j in range(i+1,len(npoly)):
      	    tp2 = npoly[j]
	    eqv = tp1[0].equal(tp2[0])
	    if eqv != 0:
               sexp = sexp + tp2[1]
	       sfac = sfac*eqv**(-tp2[1])
	       excl_list.append(j)
	    
         if abs(sexp)>1E-16:
	    nnpol.append((tp1[0],sexp))
	 
	 nfmult = nfmult*sfac   
	       
      
      
      self.polynomes = nnpol
      self.fmult = nfmult

   def equal(self, prod2):
      """
      Check whether two products are the same and return 
      the ratio of multiplicative factors r: self = r*prod2
      """
      
      if len(self.polynomes) != len(prod2.polynomes):
         return 0

      r=1.0*self.fmult/prod2.fmult
      
      for i,p1 in enumerate(self.polynomes):
         found=False
	 fac=1
	 for j,p2 in enumerate(prod2.polynomes):
	    if p1[1] != p2[1]:
	       continue
	    
	    eqv = p1[0].equal(p2[0])
	    if eqv == 0:
	       continue
	    
	    fac = eqv**p1[1]
	    found=True
	    break
	    
	 if found is False:
	    return 0
	    
	 r = r*fac
	 
      return r
            
         
############################################################
############################################################

class pprod_sum:
   def __init__(self, lsum):
      """
      Define a sum of polynomial products. 
      lsum is the only argument that is needed for this and is a 
      LIST of polynome_product objects.
      The multiplicative factors in the sum are already defined in the 
      polynome_product objects.
      """
      
      self.lsum=list(lsum)
	
      self.simplify()
      
   def __add__(self, p2):
      nlst = self.lsum + p2.lsum
      return pprod_sum(nlst)
      
   
   def simplify(self):
      """
      Sum equal terms together
      """
      nlist=[]
      
      excl=[]
      
      for i,cpol in enumerate(self.lsum):
         nfac = cpol.fmult
	 if i in excl:
	    continue
	    
	 for j in range(i+1, len(self.lsum)):
	    polj = self.lsum[j]
	    eqv = cpol.equal(polj)
	    
	    if eqv == 0:
	       continue
	       
	    excl.append(j)
	    nfac = nfac*(1 + 1.0/eqv)
	    
	 cpol.fmult=nfac
	 nlist.append(cpol)
	 
      self.lsum = nlist
         
      
   def printout(self):
      if len(self.lsum) == 0:
         print "0"
	 
      for lf in self.lsum:
         lf.printout()
	 print "  +  ",
	 
      print "\n"

   def derive(self):
      """
      Derivative of a sum of polynomials. Returns another pprod_sum 
      object
      """
      
      dersum=pprod_sum([])
      for pols in self.lsum:
         nterm = pols.derive()
	 dersum = dersum + nterm
	 
      return dersum
      
   def nderive(self,n):
      """
      nth order derivative. Just calls derive n times
      """
      if n < 0:
         return pprod_sum([])
      
      if n == 0:
         return self

      fder = self.derive()
      for i in range(n-1):
         fder = fder.derive()
	 
      return fder
            

   def eval(self,z):
      """
      Evaluate this sum
      """
      vret = complex(0,0)
      zvl=complex(z)
      
      for pols in self.lsum:
         vret = vret + pols.eval(z)
	 
      if abs(vret.imag) < 1E-16 and not z is complex:
         vret = float(vret.real)
	 
      return vret     
