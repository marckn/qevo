lcell=1

#hcp lattice
a= 4.257
Lx= a
Ly = a*3.0**0.5
Lz = a*1.633

f=open("atompos.csv","w+")
fx=open("atompos.xyz","w+")

for i in xrange(-lcell,lcell+2):
   for j in xrange(-lcell,lcell+2):
      for k in xrange(-lcell,lcell+2):
         cx= -Lx/2 + i*Lx
	 cy= -Ly/2 + j*Ly
	 cz= -Lz/2 + k*Lz
	 
	 at1 = (cx,cy,cz)
	 at2 = (cx+Lx/2,cy+Ly/2,cz)
	 at3 = (cx+Lx/2,cy+Ly/6,cz+Lz/2)
	 at4 = (cx+Lx,cy+2*Ly/3,cz+Lz/2)
	 
	 atoms=[at1,at2,at3,at4]
	 
	 for at in atoms:
	    pl="{0:.4f},{1:.4f},{2:.4f}\n".format(*at)
	    pl2="HE {0:.4f} {1:.4f} {2:.4f}\n".format(*at)
	    
	    f.write(pl)
	    fx.write(pl2)
