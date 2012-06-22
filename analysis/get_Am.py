import math as M
import pylab as P
import numpy as N

filename="ascii.00009047"
MMAX=12
rmax=40.
deltax=0.5
deltaz=0.5
deltay=0.5
h=3*deltax
lhex=deltax/M.tan(M.pi/6.)
vol=deltax*deltay*deltaz

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

n1=101
n2=101

f=open(filename,"r")
ihead=0
iter=0
xmax=0.;xmin=0.;ymax=0.;ymin=0.;zmax=0.;zmin=0.
for line in f:
	val=line.split()
	if ihead==2:
		time=float(val[1])
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
	if ihead < 6:
		ihead+=1
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
	iter+=1
print "#read file"
f.close()

am=[]
bm=[]
for i in xrange(MMAX):
	am.append(0.)
	bm.append(0.)

n=len(x)

for im in xrange(MMAX):
	for i in xrange(n):
		rcyl=M.sqrt(x[i]**2+y[i]**2)
		if rcyl>rmax:continue
		if rcyl==0.:ang=0.
		else: ang=M.atan2(y[i],x[i])
		am[im]+=rho[i]*M.cos(im*ang)*vol/M.pi
		bm[im]+=rho[i]*M.sin(im*ang)*vol/M.pi
	

print am[0],bm[0]
for im in xrange(1,MMAX):
	print im, M.sqrt(am[im]**2+bm[im]**2)/am[0]


