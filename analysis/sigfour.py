import os
filelist="../time_step.log"
prefix="../celldump."
outdir="./"
outfile="fourier.log"

fl=open(filelist,"r")

iter=242
for line in fl:
	val=line.split()
	step=int(val[6])
	time=float(val[5])
	filename=prefix+repr(step).zfill(8)
	command="./sigfourier  "+filename+ " >> "+outdir+outfile
	print command
	os.system(command)
	iter+=1
