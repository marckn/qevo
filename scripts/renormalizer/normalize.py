from os import system as sys
from os import mkdir
import re
import math

fmax=1500
dirpath="../../applications/planewaves/mov"
outpath="outsnap/"

recols = r"\S+\b"
rex = re.compile(recols)

recomm = r"^\b#+.*"
rec = re.compile(recomm)

maxpsi=0.0
minpsi=0.0
maxrho=0.0


try:
   mkdir(outpath)
except:
   pass



for i in xrange(fmax):
   
   fsrc = dirpath+"/snap_{0}.dat".format(i) 
   print "Doing snapshot {0:0>5d}".format(i)
   
   fo = open(fsrc, "r+")
   
   cnorm=0
   cmaxrho=0.0
   cmaxpsi=0.0
   cminpsi=0.0
   dx=[]
   dy=[]
   
   
   knt=0
   prevcols=None
   for line in fo:
      if not rec.match(line) is None:
         continue
	 
      
      cols = rex.findall(line)
      cols=map(lambda i : float(i), cols)
      
      
      if len(cols) == 0:
         continue
	 
      if len(cols) == 3:
         rho = cols[2]

      else:
         rho=cols[2]**2 + cols[3]**2
	 cmaxpsi = max(cols[2],cols[3],cmaxpsi)
	 cminpsi = min(cols[2],cols[3],cminpsi)

      cnorm = cnorm + rho 
      cmaxrho=max(rho,cmaxrho)
      
      if knt > 0:
         cdx=abs(cols[0]-prevcols[0])
	 cdy=abs(cols[1]-prevcols[1])
	 dx.append(cdx)
	 dy.append(cdy)
	 dx=list(set(dx))
	 dy=list(set(dy))
      
      prevcols=cols
      knt = knt+1
      
   dx=sorted(dx)
   dy=sorted(dy)
   
   
   # when you subtract adjacent elements of an
   # evenly spaced grid there are only three
   # different outcomes.
   dx=dx[1]
   dy=dy[1]
   
   
   cnorm = cnorm*dx*dy
   cmaxrho = cmaxrho/cnorm
   cmaxpsi = cmaxpsi/math.sqrt(cnorm)
   cminpsi = cminpsi/math.sqrt(cnorm)
   
   

   
   maxrho=max(cmaxrho,maxrho)
   maxpsi=max(cmaxpsi,maxpsi)
   minpsi=min(cminpsi,minpsi)
   
   # writing normalized file
   fwr = open(outpath+"/snap_{0}.dat".format(i),"w+")
   print "Writing normalzed snapshot"
   
   fo.close()
   fo = open(fsrc, "r+")
   
   for line in fo:
      if not rec.match(line) is None:
         fwr.write(line)
	 continue
	 
      cols = rex.findall(line)
      cols=map(lambda i : float(i), cols)
      
      if len(cols) == 0:
         fwr.write(line)
	 continue
      
      lout="{0:10f} {1:10f}".format(cols[0],cols[1])
      if len(cols) == 3:
         lout = lout+" {0:15f}\n".format(cols[2]/cnorm)
      else:
         csqn=math.sqrt(cnorm)
         lout = lout+" {0:15f} {1:15f}\n".format(cols[2]/csqn,cols[3]/csqn)
	 
      fwr.write(lout)       
         



print "All done."
print "Automatically determined dx = "+str(dx)+" and dy = "+str(dy)

print "Maximum value of rho: "+str(maxrho)
print "Maximum value of psi: "+str(maxpsi)
print "Minimum value of psi: "+str(minpsi)


fres=open(outpath+"/zrange.dat", "w+")
fres.write("Maximum value of rho: "+str(maxrho)+"\n")
fres.write("Maximum value of psi: "+str(maxpsi)+"\n")
fres.write("Minimum value of psi: "+str(minpsi)+"\n")
      
