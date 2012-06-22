import os
filelist="./time_step.log"
prefix="../p025ds/celldump."
outdir="./"
outfile="pd."

fl=open(filelist,"r")

iter=0
for line in fl:
	val=line.split()
	step=int(val[6])
	time=float(val[5])
	filename=prefix+repr(step).zfill(8)
	command="./phase  "+filename+ " T > "+outdir+outfile+repr(iter).zfill(8)
	print command
	os.system(command)
	iter+=1
