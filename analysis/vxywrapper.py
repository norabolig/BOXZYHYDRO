import os
filelist="time_step.log"
prefix="../highopac/celldump."
outdir="./image_high_vxy/"
outpre="10mjup_high_vxy."

fl=open(filelist,"r")

iter=0
for line in fl:
	val=line.split()
	step=int(val[6])
	time=float(val[5])
	filename=prefix+repr(step).zfill(8)
	command="./vxy "+filename+ " T > "+outdir+outpre+repr(iter).zfill(8)
	print command
	os.system(command)
	iter+=1
