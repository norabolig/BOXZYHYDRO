program sdxy
 implicit none

 integer,parameter::pre=8
 integer::nx=128,ny=128,nz=64,ix,iy,iz,idown,icheck=1234
 real(pre)::ds=0.025d0,dz=0.025d0

 character*72:: filename,input!="../celldump.00003750 "
 real(pre)::deltax,deltaz,h,lhex,dx,dy,xmin,xmax,ymin,ymax,zmin,zmax,pi,vol,dist,val,xx,yy,mass
 real(pre)::phase,area,val2,time
 real(pre),dimension(:),allocatable::x,y,z,rho,p,vx,vy,vz,phi,sort
 real(pre),dimension(:),allocatable::XIM,YIM,sig
 integer,dimension(:),allocatable::indx_sort,boundary
 real(pre),dimension(:,:),allocatable::image
 character junk

 type units
  real(pre)::length,time,mass,density
  real(pre)::rgas,eps,kelvin,vel
 end type

 integer::ioerr,ipix,jpix,nentry,ihead,iter,idx,ientry,isearch,ixs,jump,idx1
 integer::nentryxy,seed
 logical::binary=.false.
 
 type(units)::scl 

 call getarg(1,filename)
 idx1=iargc()
 if(idx1/=2)then
   print *,"Please give filename on command line and binary option [.true. or .false.]"
   stop
 endif
 call getarg(2,input)
 read(input,"(L)")binary

 if (binary)then

   open(unit=100,file=trim(filename),form="UNFORMATTED")
   read(100)ihead
   if (.not.ihead==icheck)then
     print *, "Corrupted file or wrong endian.  Header flag fail."
     stop
   endif
   read(100)time,ihead
   print '(A25,1pe16.8E3,I8)', "Time of snapshot at step:",time,ihead
 else
   print *, "What are you doing?"
   stop
 endif
 
 close(100)

end program


