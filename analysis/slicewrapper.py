import os
filelist="../time_step.log"
prefix="../celldump."
outdir="./imageyz/"
outpre="10mjup_p025."

fl=open(filelist,"r")

iter=0
for line in fl:
	val=line.split()
	step=int(val[6])
	time=float(val[5])
	filename=prefix+repr(step).zfill(8)
	command="python get_sliceyz.py "+filename+ " > "+outdir+outpre+repr(iter).zfill(8)
	print command
	os.system(command)
	iter+=1
