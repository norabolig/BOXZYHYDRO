program sdxy
 implicit none

 integer,parameter::pre=8
 integer::nx=8,ny=200,nz=1,ix,iy,iz,idown,icheck=1234
 real(pre)::dy=0.01d0,dx=0.01d0,dz=0.01d0,r_vol=60d0,r,muc=2.333
 logical::use_offset=.false.

 character*72:: filename
 real(pre)::deltax,deltaz,h,lhex,xmin,xmax,ymin,ymax,zmin,zmax,pi,vol,dist,val,xx,yy,mass
 real(pre)::phase,area,val2,time,dmax,xoff,yoff,zoff
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

 pi=acos(-1d0)

 call getarg(1,filename)
 idx1=iargc()
 if(idx1/=1)then
   print *,"Please give filename on command line and binary option [.true. or .false.]"
   stop
 endif

 open(unit=100,file=trim(filename),form="UNFORMATTED")
 read(100)ihead
 if (.not.ihead==icheck)then
   print *, "Corrupted file or wrong endian.  Header flag fail."
   stop
 endif
 read(100)time,ihead
 read(100)scl
 ioerr=0
 nentry=0
 do while (ioerr>-1)
   read(100,iostat=ioerr)
   nentry=nentry+1
 enddo
 nentry=nentry-1
 print *, "#Found ",nentry," entries for step ",ihead
 print *, "#Scalings",scl
 rewind(100)
 
 allocate(x(nentry))
 allocate(sort(nentry))
 allocate(y(nentry))
 allocate(z(nentry))
 allocate(rho(nentry))
 allocate(p(nentry))
 allocate(vx(nentry))
 allocate(vy(nentry))
 allocate(vz(nentry))
 allocate(phi(nentry))
 allocate(boundary(nentry))
 allocate(indx_sort(nentry))

 if(nentry/=ny*nx*nz)then
   print *, 'ordering off. check dimensions and entry ',nx*ny*nz,nentry," diff ",nx*ny*nz-nentry
   stop
 endif
 
 dmax=0d0
 read(100)ihead
 read(100)time,ihead
 read(100)scl
 xmin=0.;xmax=0.;ymin=0.;ymax=0.;zmin=0.;zmax=0.
 iter=0
 do iz=1,nz
 do ix=1,nx
 do iy=1,ny
   iter=iter+1
   read(100)x(iter),y(iter),z(iter),rho(iter),p(iter),&
      vx(iter),vy(iter),vz(iter),phi(iter),boundary(iter)
   if(rho(iter)>dmax)then
     dmax=rho(iter)
     xoff=x(iter)
     yoff=y(iter)
     zoff=z(iter)
   endif
   indx_sort(iter)=iter
   sort(iter)=x(iter)
   if(xmax<x(iter))xmax=x(iter)
   if(xmin>x(iter))xmin=x(iter)
   if(ymax<y(iter))ymax=y(iter)
   if(ymin>y(iter))ymin=y(iter)
   if(zmax<z(iter))zmax=z(iter)
   if(zmin>z(iter))zmin=z(iter)
 enddo
 enddo
 enddo
 close(100)

 nentryxy=nentry/nz

 allocate(XIM(nx))
 allocate(YIM(ny))
 allocate(sig(nentryxy))
 allocate(image(nx,ny))

 
 print *, '#time ',time
 print *, '#nrntry, nentryxy, nx ny nz, nx ny, ', nentry, nentryxy, nx*ny*nz, nx*ny
 print *, '#maxden ',dmax,' at ',xoff,yoff,zoff
 print *, '#Using offset -> ',use_offset

 do iter=1,nentry
     write(6,"(9(1pe16.8e3,1X),I2)")x(iter),y(iter),z(iter),rho(iter),p(iter),&
        vx(iter),vy(iter),vz(iter),phi(iter),boundary(iter)
 enddo
 
 deallocate(x,y,z,rho,p,vx,vy,vz,phi,boundary)


end program

