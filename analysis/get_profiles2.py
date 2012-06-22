import math as M
import pylab as P
import numpy as N
import sys

filename=sys.argv[1]
rmax=1.5
deltax=0.025
deltaz=0.025
muc=2.33
lhex=deltax*M.tan(M.pi/6.)
vol=lhex*deltax*deltaz*1.5
print "#",lhex,deltax

eostable="engtable.dat"
ttablemin=3.
dtk_eos=5.

au=1.496e13
msun=1.989e33
opacfac=0.1
dlimit=1e-7


def opac(T):
        if(T<80.0):
        	get_kappa=(T**2)/3200.0
        elif(T<170.0):
        	get_kappa=-2.0 + 0.050*T
        elif(T<180.0):
        	get_kappa=62.60 - 0.330*T
        elif(T<270.0):
        	get_kappa=-1.0 + 0.0233*T
        elif(T<300.0):
        	get_kappa=8.0 - 0.010*T
        elif(T<425.0):
        	get_kappa=1.88 + 0.0104*T
        elif(T<440.0):
        	get_kappa=128.13 - 0.2867*T
        elif(T<670.0):
        	get_kappa=0.57 + 0.0033*T
        elif(T<700.0):
        	get_kappa=19.50 - 0.0250*T
        elif(T<1300.0):
        	get_kappa=-0.33 + 0.0033*T
        elif(T<1350.0):
        	get_kappa=24.80 - 0.0160*T
        elif(T<1449.66):
                get_kappa=46.40 - 0.0320*T
        else:
	        get_kappa=0.01
	return get_kappa*msun/au**2*opacfac

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
gtable=[]
ttable=[]

f=open(eostable,"r")
for line in f:
	val=line.split()
	ttable.append(float(val[0]))
	gtable.append(float(val[2]))

f.close()

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
gshell=[]
beta_ad=[]
beta_rad=[]
tau=[]
dtau=[]
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
	gshell.append(0.)
	beta_ad.append(0.)
	beta_rad.append(0.)
	tau.append(0.)
	dtau.append(0.)

n=len(x)
ang=0.
for i in xrange(n):
	rsph=M.sqrt(x[i]**2+y[i]**2+z[i]**2)
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
ntable=len(ttable)
for ishell in xrange(nr):
	tmass+=mshell[ishell]
	vshell=(4.*M.pi*dr**3*(float(i+1)**3-float(i)**3)/3.)
	if nshell[ishell]>0.:
		dshell[ishell]/=nshell[ishell]
		pshell[ishell]/=nshell[ishell]
		tshell[ishell]/=nshell[ishell]
		itable=min(max(int((tshell[ishell]-ttablemin)/dtk_eos),0),ntable-1)
		itable2=min(itable+1,ntable-1)
		gshell[ishell]=gtable[itable]+(gtable[itable2]-gtable[itable])/dtk_eos*(tshell[ishell]-ttable[itable])
		beta_ad[ishell]=(gshell[ishell]-1.)/gshell[ishell]

for ishell in xrange(nr-1):
	beta_rad[ishell]=pshell[ishell]/tshell[ishell]*(tshell[ishell+1]-tshell[ishell])/(pshell[ishell+1]-pshell[ishell])
beta_rad[nr-1]=beta_rad[nr-2]

for ishell in xrange(nr-2,0,-1):
	dtau[ishell]=opac(tshell[ishell])*dshell[ishell]*dr
	if dshell[ishell]>dlimit:
		tau[ishell]=dtau[ishell]+tau[ishell+1]
	else: tau[ishell]=tau[ishell+1]

print "#time ang.mom.", time, ang,mshell[nr-1]
for ishell in xrange(nr):
	print rshell[ishell],angshell[ishell],mshell[ishell],angtot[ishell],mtot[ishell],dshell[ishell],pshell[ishell],tshell[ishell],gshell[ishell],beta_ad[ishell],beta_rad[ishell],tau[ishell],dtau[ishell]
	

print "tmass is ",tmass,rgascode

