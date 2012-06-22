import os
import sys

os.system('ls ../spin_destroy_slow/celldump* > tmpfile')
f=open('tmpfile','r')

iter=0
for line in f:
        filename=line.rstrip()
	command="./get_time_step "+filename+ " T "
	#print command
	os.system(command)
	iter+=1

os.system('rm -f tmpfile')
