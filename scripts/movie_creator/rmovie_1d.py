from os import system as sys
from os import mkdir


framerate=10
one_only=False

psirange=[-0.4,0.4]
rhorange=[0,0.1]
range_x=[-35.0,35.0]


movname="tunneling5"
input_dir="../../applications/tunneling_1d/mov/normalized_snapshots"
workdir="workdir"



fmax=2000
fjump=1



try:
   mkdir(workdir)
except:
   pass


if one_only:
   comm = """
   gnuplot -e "filename='{0}/snap_{1}.dat' ; fnout='{2}/snap_{3:0>5d}.png' ; \
   zmin='{4:.4f}' ; zmax='{5:.4f}' ;  xmin='{6:.4f}' ; xmax='{7:.4f}'" script_single_1d.plg \ 
   """
   ranges=(psirange[0], psirange[1],range_x[0], range_x[1])
else:
   comm = """
   gnuplot -e "filename='{0}/snap_{1}.dat' ; fnout='{2}/snap_{3:0>5d}.png' ; \
   rhorange='{4:.4f}' ; psimin='{5:.4f}' ; psimax='{6:.4f}' ; \
   xmin='{7:.4f}' ; xmax='{8:.4f}'" script_multi_1d.plg \
   """
   ranges=(rhorange[1],psirange[0],psirange[1],range_x[0], range_x[1])
   

cnt=0
for i in xrange(0,fmax,fjump):
   print "Doing snapshot {0:0>5d}".format(i)
   cmd=comm.format(input_dir,i,workdir,cnt,*ranges)
   print cmd
   sys(cmd)
   cnt=cnt+1
   

cmdmov= """
ffmpeg -framerate {0} -i {1}/snap_%05d.png -c:v libx264 -r 30 -pix_fmt yuv420p {2}.mp4
""".format(framerate,workdir,movname)   

sys(cmdmov)
