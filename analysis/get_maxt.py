import math as M
import sys

if len(sys.argv)!=2:
	print "You must give the log filename."
	exit()

fh=open(sys.argv[1],'r')

iter=0
for line in fh:
	lineold=line
	if iter==0:continue
	if line.find('MaxT')<0:continue
	print lineold.rstrip()+line.rstrip()	

fh.close()
