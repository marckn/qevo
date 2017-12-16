from os import system as sys
from os import mkdir
import re
import math

fmax=2000
dirpath="../../applications/tunneling_1d/mov"
outpath="outsnap/"

# NOT USED HERE
#comm = """
#gnuplot -e "filename='../snap_{0:0>5d}.dat' ; fnout='snap_{1:0>5d}.png'" script.plg
#"""

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
   
   fsrc = dirpath+"/snap_{0}.dat.dat".format(i) 
   print "Doing snapshot {0:0>5d}".format(i)
   
   fo = open(fsrc, "r+")
   
   cnorm=0
   cmaxrho=0.0
   cmaxpsi=0.0
   cminpsi=0.0
   dx=[]
   
   
   knt=0
   prevcols=None
   for line in fo:
      if not rec.match(line) is None:
         continue
	 
      
      cols = rex.findall(line)
      cols=map(lambda i : float(i), cols)
      
      
      if len(cols) == 0:
         continue
	 
      if len(cols) == 2:
         rho = cols[1]

      else:
         rho=cols[1]**2 + cols[2]**2
	 cmaxpsi = max(cols[1],cols[2],cmaxpsi)
	 cminpsi = min(cols[1],cols[2],cminpsi)

      cnorm = cnorm + rho 
      cmaxrho=max(rho,cmaxrho)
      
      if knt > 0:
         cdx=abs(cols[0]-prevcols[0])
	 dx.append(cdx)
	 dx=list(set(dx))
	 
      prevcols=cols
      knt = knt+1
      
   dx=sorted(dx)
   
   
   # when you subtract adjacent elements of an
   # evenly spaced grid there are only three
   # different outcomes.
   dx=dx[1]
   
   
   cnorm = cnorm*dx
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
      
      lout="{0:10f}".format(cols[0])
      if len(cols) == 2:
         lout = lout+" {0:15f}\n".format(cols[1]/cnorm)
      else:
         csqn=math.sqrt(cnorm)
         lout = lout+" {0:15f} {1:15f}\n".format(cols[1]/csqn,cols[2]/csqn)
	 
      fwr.write(lout)       
         



print "All done."
print "Automatically determined dx = "+str(dx)

print "Maximum value of rho: "+str(maxrho)
print "Maximum value of psi: "+str(maxpsi)
print "Minimum value of psi: "+str(minpsi)


fres=open(outpath+"/zrange.dat", "w+")
fres.write("Maximum value of rho: "+str(maxrho)+"\n")
fres.write("Maximum value of psi: "+str(maxpsi)+"\n")
fres.write("Minimum value of psi: "+str(minpsi)+"\n")
      
