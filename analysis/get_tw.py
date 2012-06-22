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

x=[]
y=[]
z=[]
rho=[]
p=[]
vx=[]
vy=[]
vz=[]
phi=[]
boundary=[]
X=[]
Y=[]
image=[]

n1=301
n2=301

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
	x.append(float(val[0]))
	if x[iter]>xmax:xmax=x[iter]
	if x[iter]<xmin:xmin=x[iter]
	y.append(float(val[1]))
	if y[iter]>ymax:ymax=y[iter]
	if y[iter]<ymin:ymin=y[iter]
	z.append(float(val[2]))
	if z[iter]>zmax:zmax=z[iter]
	if z[iter]<zmin:zmin=z[iter]

	rho.append(float(val[3]))
	p.append(float(val[4]))
	vx.append(float(val[5]))
	vy.append(float(val[6]))
	vz.append(float(val[7]))
	phi.append(float(val[8]))
	boundary.append(int(val[9]))
	iter+=1
print "#read file"
f.close()

n=len(x)
T=0.;W=0.;massup=0.;masslo=0.
for i in xrange(n):
	rcyl=M.sqrt(x[i]**2+y[i]**2)
	if rcyl>0.: T+=.5*((x[i]*vy[i]-y[i]*vx[i])/rcyl)**2*rho[i]*vol
	W+=.5*rho[i]*phi[i]*vol
	if(z[i]>=0.):massup+=rho[i]*vol
	if(z[i]<=0.):masslo+=rho[i]*vol
print "T,W,T/W", time, T,W,T/W,massup,masslo


