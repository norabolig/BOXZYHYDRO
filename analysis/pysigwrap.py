import os
DIR="./spin_destroy_slow/"
filelist=DIR+"renderlist"
bposfile="../spin_destroy_slow/b.pos"
movdir="./mov/"

fl=open(filelist,"r")

iter=0
for line in fl:
	filename=line.rstrip()
	command="python makesig.py "+DIR+filename+" "+bposfile+" "+DIR+filename+".png "+repr(iter)
	print command
	os.system(command)
	command="convert "+DIR+filename+".png "+movdir+repr(iter)+".jpg"
	os.system(command)
	iter+=1
