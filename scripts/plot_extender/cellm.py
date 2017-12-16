fnm="psi_0.dat.csv"
lcell=1

with open(fnm,"r") as f:
   raw = f.read()
   
rbuff=[]
ln=""
for c in raw:
   if c != "\n":
      ln=ln+c
   else:
      lp=ln.split(",")
      if len(lp) > 0 and lp[0] != '':
         rbuff.append(lp)
      
      ln=""

mxm = float(rbuff[0][0])
mym = float(rbuff[0][1])
mzm = float(rbuff[0][2])

Lx= 2*abs(mxm)
Ly= 2*abs(mym)
Lz= 2*abs(mzm)


fo = open("rho_big.csv","w+")

for i in xrange(-lcell,lcell+1):
   for j in xrange(-lcell,lcell+1):
      for k in xrange(-lcell,lcell+1):
         xbase = i*Lx
	 ybase = j*Ly
	 zbase = k*Lz
	 
	 for el in rbuff:
	    if len(el) == 0:
	       continue
	       
	    cx=float(el[0])+xbase
	    cy=float(el[1])+ybase
	    cz=float(el[2])+zbase
	    
	    sto = "{0:.4f},{1:.4f},{2:.4f},{3:.8f}\n".format(cx,cy,cz,float(el[5]))
            fo.write(sto)   

fo.close()
      
