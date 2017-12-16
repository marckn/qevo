"""
List of dictionaries containing keywords
"func","M","bpars","caps","factor","ID","dim"
where func is the function and M the corresponding matrix 
of the integrals. The integrator will try to see if the 
required element can be found before starting the 
Monte Carlo integration. 
"""

import numpy as np
from aux.enum import enumerator as enum
from inspect import getargspec as nargs
import sys

onebody_funcs=[]
twobody_funcs=[]


# compound functions
def compoundFunc(newidfunc, fassoc, lidfunc, lfactors):
   """
   From n functions already stored build the elements for a new one which is 
   a linear combination with coefficients lfactors.
   """ 
   
   arf = nargs(fassoc)
   body = int(len(arf.args))
   
   if body == 1:
      targ = onebody_funcs
   else:
      targ = twobody_funcs
      
   
   idxs=[]
   for i,entry in enumerate(targ):
      cid = entry["ID"]
      if cid in lidfunc:
         ilid = lidfunc.index(cid)
         idxs.append((cid,i,ilid))
	 

   if len(idxs) != len(lidfunc):
      print "Some id were not found"
      raise ValueError
      
   
   # consistency check
   caps = targ[idxs[0][1]]["caps"]
   bpars = targ[idxs[0][1]]["bpars"]
   dim = len(caps)
   for tp in idxs[1:]:
      if targ[tp[1]]["caps"] != caps:
         print "Consistency error: different matrices sizes"
	 raise ValueError
	 
      if targ[tp[1]]["bpars"] != bpars:
         print "Consistency error: different bpars"
	 raise ValueError

	 
      
   newm = np.zeros(targ[idxs[0][1]]["M"].shape)
   
   for fac,tp in zip(lfactors,idxs):
      Mtadd = targ[tp[1]]["M"]
      fcadd = targ[tp[1]]["factor"]
      newm = newm + float(fac)*float(fcadd)*Mtadd
      
   
   ttadd = {"func":fassoc,"M":newm,"bpars":bpars,"caps":caps,"factor":1.0,"ID":newidfunc,"dim":dim}
   targ.append(ttadd)      
   


def loadFunc(idfunc,fassoc,factor=1, fdir="./data", symm=True):
   """
   Scan through the data file in search for the function with ID idfunc and 
   add the relative entry to the onebody_func list.
   """
   
   print "Loading ",fassoc," elements from table entry ",idfunc
   
   arf = nargs(fassoc)
   body = int(len(arf.args))
      
   if body == 1:
      targ = onebody_funcs
      fpath = fdir+"/onebody.dat"
      dimcoord = 2
   else:
      targ = twobody_funcs
      fpath = fdir+"/twobody.dat"
      dimcoord = 4
      
   for entry in targ:
      if idfunc == entry["ID"]:
         print idfunc," already loaded"
	 return
	 
   f = open(fpath,"r")
   found=False
   datablock=[]
   
   for fline in f:
      fline=fline[:-1]
      if found is False:
         if fline.find("ID") == -1:
            continue
	 
         res = fline.split("=")
         fid = res[1]
	 if fid == idfunc:
            found = True
	    inblock=True 
	    datablock.append(fline)  

      else:
         if inblock:
            if fline.find("ID") != -1:
	       inblock=False
	       break
	 
	    datablock.append(fline)
	    
   
   if not found:
      print "Cannot find ",idfunc
      raise ValueError
      
   
   # reading header
   dim = datablock[1]
   bpar = datablock[2]
   caps = datablock[3]
  
   
   dim = dim.split("=")
   dim = int(dim[1])
   
   def getTuple(st,dtype=str):
      data= st.split("=")
      if len(data)>1:
         data=data[1]
      else:
         data = data[0]
	  
      vals = data.split(",")
      vals = map(lambda i : dtype(i), vals)
      return tuple(vals)
   
   bpar = getTuple(bpar,float)
   caps = getTuple(caps,int)
   
   maxn = reduce(lambda i,j: i*j, caps, 1)
   mdim = [maxn]*(2*body)
   
   
   matrix = np.zeros(mdim)
   
   for dataline in datablock[4:]:
      tp = getTuple(dataline,float)
      coord = tuple(map(lambda i : int(i), tp[:-1]))
      matrix[coord] = tp[-1]
      if symm is True:
         if body == 1:
	    i,j=coord[0],coord[1]
	    matrix[j,i]=tp[-1]
	 if body == 2:
	    i,j,k,l=coord[0],coord[1],coord[2],coord[3]
	    matrix[i,j,k,l]=tp[-1]
	    matrix[j,i,l,k]=tp[-1]
	    matrix[i,j,l,k]=tp[-1]
	    matrix[k,l,i,j]=tp[-1]
	    matrix[k,l,j,i]=tp[-1]
	    matrix[l,k,j,i]=tp[-1]
	    matrix[l,k,i,j]=tp[-1]
      
   
   ttadd = {"func":fassoc,"M":matrix,"bpars":bpar,"caps":caps,"factor":factor,"ID":idfunc,"dim":dim}
   targ.append(ttadd)
   
   
def writeFunc(idfunc,fassoc,bpars,caps,mcpoints,fdir,
                integrator, iszero=lambda i : False, analytical = lambda i : None, symm = True):
   """
   Compute and write the integrals to the corresponding file.
   IN: idfunc   =   the ID of the entry
       fassoc   =   the function to be used
       caps     =   a tuple defining the caps for each quantum number
       mcpoints =   the number of MC steps of integration
       iszero   =   a given function that can tell from the levels whether the integral will be 
                    null or not. This function takes a tuple as argument, structured as 
		    (t1,t2,...tx) where tx are tuples containing the levels. Note that 
		    for 1 body the arg is (t1,t2) and for 2 bodies (t1,t2,t3,t4)
     analytical =   a function that given the same arguments as iszero computes the integral analytically.
     symm       =   it's a fully symmetric tensor, that means i<=j for 1-body and i<=j<=k<=l for 2-body 
   """	 
   
   dim = len(caps)
   arf = nargs(fassoc)
   body = int(len(arf.args))
   
   if body == 1:
      fpath=fdir+"/onebody.dat"
   else:
      fpath=fdir+"/twobody.dat"
      
   # check if it already exists
   try:
      f = open(fpath, "r")
   
      for ln in f:
         ln=ln[:-1]
	 if ln.find("ID=") != -1:
	    (n,val) = ln.split("=")
	    if val == idfunc:
	       print idfunc," already present in file."
	       raise ValueError
	    
      f.close()
   except IOError:
      pass
   
   
   f = open(fpath,"a")
      
   enu = enum(caps)
   maxsize=enu.maxidx
   f.write("ID="+idfunc+"\n")
   f.write("DIM="+str(dim)+"\n")
   
   bpline="BPAR="
   for bp in bpars:
      bpline = bpline+"{0:14.10f},".format(float(bp))
   
   bpline = list(bpline)
   bpline="".join(bpline[:-1])+"\n"
   f.write(bpline)
   
   capsline="CAPS="
   for cp in caps:
      capsline=capsline+"{0:3d},".format(int(cp))
   
   capsline=list(capsline)
   capsline="".join(capsline[:-1])+"\n"
   f.write(capsline)
   
   fcaps = [maxsize]*(2*body)
   fenu = enum(fcaps)

   print "Header done. \n Writing datablock \n [ ",
   sys.stdout.flush()
   percblock = fenu.maxidx/10
   
   for i in xrange(fenu.maxidx):
      flev = fenu.tolevels(i)
      
      llist=[]
      for l in flev:
         ltadd = enu.tolevels(l)
	 llist.append(ltadd)
      
      isord=True
      if symm is True:
         if body == 1:
	    if flev[0]>flev[1]:
	       isord=False
	 
	 else:
	    if flev[0]>flev[1] or flev[0]>flev[2] or flev[2]>flev[3]:
	       isord=False
	          
      val = 0
      if not iszero(llist) and isord:
         val = analytical(llist)
	 if val is None:
            val = integrator(bpars,*llist,func=fassoc,borders=None,nint=mcpoints)
         
	 if abs(val)>1E-4:
	    ln=""
            for l in flev:
               ln=ln+"{0:3d},".format(int(l))
	 
            ln=ln+"{0:20.10f}\n".format(float(val))
            f.write(ln)
      
      if (i+1)%percblock == 0:
         print "*",
	 sys.stdout.flush()
      
   print "] \n Done."      
   
   


def find_precomputed_1b(bpars,l1,l2,V1):
   """
   Scan through the onebody_funcs list to find the 
   entry corresponding to these values.
   """
   
   
   for idx,entry in enumerate(onebody_funcs):
      if entry["dim"] != len(bpars):
         continue
	 
      if entry["func"] != V1:
         continue
      
      stop=False
      
      for b1,b2 in zip(entry["bpars"], bpars):
         if abs(float(b1) - float(b2)) > 1E-2:
	    stop=True
	    break
      
      if stop:
         continue
      
      for x1,x2,x3 in zip(entry["caps"],l1,l2):
         if x1 <= x2 or x1 <= x3:
	    stop=True
	    break
	    
      if stop:
         continue
	 
      return idx
      
   return False

def precomputed_integral_1b(l1,l2,idx):
   """
   Return the precomputed value stored at position idx in onebody_funcs with matrix 
   elements l1,l2
   """
   
   targ = onebody_funcs[idx]
   enu = enum(targ["caps"])
   n1 = enu.toidx(l1)
   n2 = enu.toidx(l2)
   return targ["M"][n1,n2]*targ["factor"]


def find_precomputed_2b(bpars,l1,l2,l3,l4,V2):
   """
   Scan through the twobody_funcs list to find the 
   entry corresponding to these values. Can be merged with the onebody version
   """
   
   for idx,entry in enumerate(twobody_funcs):
      if entry["dim"] != len(bpars):
         continue
	 
      if entry["func"] != V2:
         continue
      
      stop=False
      
      for b1,b2 in zip(entry["bpars"], bpars):
         if abs(float(b1) - float(b2)) > 1E-2:
	    stop=True
	    break
      
      if stop:
         continue
      
      for x1,x2,x3,x4,x5 in zip(entry["caps"],l1,l2,l3,l4):
         if x1 <= x2 or x1 <= x3 or x1 <= x4 or x1 <= x5:
	    stop=True
	    break
	    
      if stop:
         continue
	 
      return idx
      
   return False

def precomputed_integral_2b(l1,l2,l3,l4,idx):
   """
   Return the precomputed value stored at position idx in onebody_funcs with matrix 
   elements l1,l2
   """
   
   targ = twobody_funcs[idx]
   enu = enum(targ["caps"])
   n1 = enu.toidx(l1)
   n2 = enu.toidx(l2)
   n3 = enu.toidx(l3)
   n4 = enu.toidx(l4)
   return targ["M"][n1,n2,n3,n4]*targ["factor"]

