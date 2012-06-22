import pylab as P
import sys

time0=0.
Select=0.

DIR='pimage/'

file_list=open("plist","r")

count=0
for file in file_list:
	print count
	handle=open(file.rstrip(),"r")
	iter=0
	x=[]
	r=[]
	for line in handle:
		if iter==0:
			iter+=1
			val=line.split()
			time=(float(val[0])-time0)/6.28
			continue
		val=line.split()
		if Select>0.:
			if not Select*.9/1.5e13 < float(val[11]) < Select*1.1/1.5e13: continue
		x.append(float(val[2]))
		r.append((float(val[2])**2+float(val[3])**2+float(val[4])**2)**0.5)
	
	handle.close()
	
	P.figure()
	P.xlabel('X (AU)')
	P.ylabel('R (AU)')
	P.scatter(x,r,marker='o',c='black',s=0.1)
	P.xlim(-1.2,1.2)
	P.ylim(0,1.0)
	#P.axis('equal')
	P.title('Time = '+'{0:0.2f}'.format(time)+' yr')
	P.savefig(DIR+'part_xrad'+repr(count)+'.png')
	count+=1
	
file_list.close()	
