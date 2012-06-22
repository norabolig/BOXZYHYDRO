import pylab as P
import math as M
import sys

MAX=1e-9
MIN=1e-15
NX=128
NY=128

conv=5.94e-7

FORMAT='%3.1f'
#TICKS=[-15,-14.5,-14,-13.5,-13,-12.5,-12,-11.5,-11,-10.5,-10,-9.5,-9]
TICKS=[-15,-14,-13,-12,-11,-10,-9]
x0=-15
x1=15
y0=-15
y1=15

if len(sys.argv)<5:
	print("./command file_density file_binary output binary_pos")
	sys.exit()
ifile=int(sys.argv[4])
fb=file(sys.argv[2],"r")
fd=file(sys.argv[1],"r")

xs=[];ys=[]
for line in fb:
	val=line.split()
	xs.append(float(val[3])*M.cos(float(val[4])))
	ys.append(float(val[3])*M.sin(float(val[4])))

x=[];y=[];v=[]
ix=0;iy=-1
xfirst=True
iter=0
for line in fd:
	if line[0]=="#" or line[1]=="#":continue
	val=line.split()
	if ix==0:
		iy+=1
		v.append([])
		y.append(float(val[1]))
	v[iy].append(M.log10(float(val[2])*conv))
	if xfirst:x.append(float(val[0]))
	ix+=1
	if ix==NX:
		ix=0
		xfirst=False
	iter+=1

print len(x),len(y),len(v),len(v[0])

fb.close()
fd.close()

lev=[]
logMAX=M.log10(MAX)
logMIN=M.log10(MIN)
for i in xrange(256):
	lev.append(logMIN+( (logMAX-logMIN)/255.*i ))

P.rcParams['xtick.major.size']=10
P.rcParams['ytick.major.size']=10
P.rcParams['axes.labelsize']=18
P.rcParams['xtick.labelsize']=18
P.rcParams['ytick.labelsize']=18
#P.rcParams['ytick.color']='white'
#P.rcParams['xtick.color']='white'

fig=P.figure(figsize=(10,8))
ax1=fig.add_axes([0.1,0.1,0.85*4/5.,0.85])
cp=ax1.contourf(x,y,v,levels=lev,antialiased=False,cmap=P.cm.jet)#hot,extend='both')
ax1.set_xlim(x0,x1)
ax1.set_ylim(y0,y1)
ax1.set_xlabel("AU")
ax1.set_ylabel("AU")

ax2=fig.add_axes([0.81,0.1,0.20*4/5,0.20])
ax2.scatter(0,0,c='black',marker='+')
ax2.scatter(xs,ys,c='b',edgecolor='b')
ax2.scatter(xs[ifile],ys[ifile],c='r',s=60)
ax2.set_xlim(-100,100)
ax2.set_ylim(-100,100)
ax2.set_xticks([])
ax2.set_yticks([])

ax3=fig.add_axes([0.85,0.35,0.025,0.60])
cb=P.colorbar(cp,cax=ax3,extend='both',format=FORMAT,ticks=TICKS)


P.savefig(sys.argv[3])


