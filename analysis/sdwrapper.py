import os
filelist="./spin_destroy_slow/time_step.log"
prefix="../spin_destroy_slow/celldump."
outdir="./spin_destroy_slow/"
outpre="sds."

fl=open(filelist,"r")

iter=0
for line in fl:
	val=line.split()
	step=int(val[6])
	time=float(val[5])
	filename=prefix+repr(step).zfill(8)
	command="./dxy "+filename+ " T > "+outdir+outpre+repr(iter).zfill(8)
	print command
	os.system(command)
	iter+=1
