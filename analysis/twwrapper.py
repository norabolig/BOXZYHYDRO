import os
filelist="../time_step.log"
prefix="../celldump."
outdir="./"
outfile="tw.log"

fl=open(filelist,"r")

iter=242
for line in fl:
	val=line.split()
	step=int(val[6])
	time=float(val[5])
	filename=prefix+repr(step).zfill(8)
	command="python get_tw.py  "+filename+ " >> "+outdir+outfile
	print command
	os.system(command)
	iter+=1
