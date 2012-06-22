import math as M
import pylab as P
import numpy as N
import sys

filename=sys.argv[1]
rmax=10.
deltax=0.25
deltay=0.25
deltaz=0.25
muc=2.33
lhex=deltax*M.tan(M.pi/6.)
vol=deltax*deltay*deltaz
print "#",lhex,deltax

#x0=-3.12023343788989937E-002; y0=-1.06965667105984995E-02
#z0=-1.31938906824096130E-002
x0=0.
y0=0.
z0=0.


nr=int(rmax/deltax)
dr=deltax

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
rshell=[]
angshell=[]
mshell=[]


f=open(filename,"r")
ihead=0
iter=0
xmax=0.;xmin=0.;ymax=0.;ymin=0.;zmax=0.;zmin=0.
tmass=0.
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
		velconv=float(val[8])
		epsconv=float(val[6])
		kelvinconv=float(val[7])
		rgascode=float(val[5])
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
	tmass+=rho[iter]*vol
	p.append(float(val[4]))
	vx.append(float(val[5]))
	vy.append(float(val[6]))
	vz.append(float(val[7]))
	phi.append(float(val[8]))
	#boundary.append(int(val[9]))
	iter+=1
print "#read file"
print "# tmass is ",tmass
f.close()

dshell=[]
pshell=[]
tshell=[]
angtot=[]
mtot=[]
nshell=[]
for i in xrange(nr):
	rshell.append(dr*(float(i)+0.5))
	angshell.append(0.)
	mshell.append(0.)
	dshell.append(0.)
	pshell.append(0.)
	tshell.append(0.)
	angtot.append(0.)
	mtot.append(0.)
	nshell.append(0.)

n=len(x)
ang=0.
for i in xrange(n):
	rsph=M.sqrt((x[i]-x0)**2+(y[i]-y0)**2+(z[i]-z0)**2)
	ishell=int(rsph/dr)
	if ishell>nr-1:continue
	angshell[ishell]+=rho[i]*(x[i]*vy[i]-y[i]*vx[i])*vol
	mshell[ishell]+=rho[i]*vol
        pshell[ishell]+=p[i]
	dshell[ishell]+=rho[i]
        tshell[ishell]+=p[i]*muc/rgascode/rho[i]
	ang+=rho[i]*(x[i]*vy[i]-y[i]*vx[i])*vol
	nshell[ishell]+=1.

for ishell in xrange(nr):
	angtot[ishell]=angshell[ishell]
	mtot[ishell]=mshell[ishell]

for ishell in xrange(nr-1):
	angtot[ishell+1]+=angtot[ishell]
	mtot[ishell+1]+=mtot[ishell]

tmass=0.
for ishell in xrange(nr):
	tmass+=mshell[ishell]
	vshell=(4.*M.pi*dr**3*(float(i+1)**3-float(i)**3)/3.)
	if nshell[ishell]>0.:
		dshell[ishell]/=nshell[ishell]
		pshell[ishell]/=nshell[ishell]
		tshell[ishell]/=nshell[ishell]
print "#time ang.mom.", time, ang,mshell[nr-1]
for ishell in xrange(nr):
	print rshell[ishell],angshell[ishell],mshell[ishell],angtot[ishell],mtot[ishell],dshell[ishell],pshell[ishell],tshell[ishell]
	

print "tmass is ",tmass,rgascode

