import os
DIR="./spin_destroy_slow/"
filelist=DIR+"rendertk"
bposfile="../spin_destroy_slow/b.pos"
movdir="./movtk/"

fl=open(filelist,"r")

iter=0
for line in fl:
	filename=line.rstrip()
	command="python maketk.py "+DIR+filename+" "+bposfile+" "+DIR+filename+".png "+repr(iter)
	print command
	os.system(command)
	command="convert "+DIR+filename+".png "+movdir+repr(iter)+".jpg"
	os.system(command)
	iter+=1
