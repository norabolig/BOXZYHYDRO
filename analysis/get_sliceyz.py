import math as M
import pylab as P
import numpy as N
import sys

filename=sys.argv[1]
deltax=0.25
deltaz=0.25
h=3*deltax
lhex=deltax*M.tan(M.pi/6.)
vol=lhex*deltax*deltaz*1.5
print "#",lhex,deltax

f=open(filename,"r")
ihead=0
iter=0
xmax=0.;xmin=0.;ymax=0.;ymin=0.;zmax=0.;zmin=0.
for line in f:
	val=line.split()
	if ihead==0:
		time=float(val[5])
		ihead=ihead+1
		continue
	if ihead==1:
		lconv=float(val[1])
		timeconv=float(val[2])
		dconv=float(val[3])
		massconv=float(val[4])
		velconv=float(val[5])
		epsconv=float(val[6])
		kelvinconv=float(val[7])
		rgascode=float(val[8])
		ihead=ihead+1
		continue

	x=(float(val[0]))
	y=(float(val[1]))
	z=(float(val[2]))
        if not x==0.:continue
        print y,z,float(val[3])

f.close()


