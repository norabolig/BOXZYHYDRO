program vyz
 implicit none

 integer,parameter::pre=8
 integer::nx=128,ny=128,nz=128,ix,iy,iz,idown,icheck=1234,idy,idz
 integer::SLICE=64
 real(pre)::dx=1.00d0,dy=1.00,dz=1.00d0,vmax=1d100

 character*72:: filename
 real(pre)::deltax,deltaz,h,lhex,xmin,xmax,ymin,ymax,zmin,zmax,pi,vol,dist,val,xx,yy,mass
 real(pre)::phase,area,val2,time,vmag
 real(pre),dimension(:),allocatable::x,y,z,rho,p,vx,vy,vz,phi,sort,eps,muc_array
 real(pre),dimension(:),allocatable::XIM,YIM,sig,sig2
 integer,dimension(:),allocatable::indx_sort,boundary
 real(pre),dimension(:,:),allocatable::image
 character junk

 type units
  real(pre)::length,time,mass,density
  real(pre)::rgas,eps,kelvin,vel
 end type

 integer::ioerr,ipix,jpix,nentry,ihead,iter,idx,ientry,isearch,ixs,jump,idx1
 integer::nentryyz,seed
 logical::binary=.false.
 
 type(units)::scl 

 pi=acos(-1d0)

 call getarg(1,filename)
 idx1=iargc()
 if(idx1/=1)then
   print *,"Please give filename on command line and binary option [.true. or .false.]"
   stop
 endif

 nentry=0

 open(unit=100,file=trim(filename),form="UNFORMATTED")
 read(100)ihead
 if (.not.ihead==icheck)then
     print *, "Corrupted file or wrong endian.  Header flag fail."
     stop
 endif
 read(100)time,ihead
 read(100)scl
 ioerr=0
 do while (ioerr==0)
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
 allocate(eps(nentry))
 allocate(muc_array(nentry))
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
      vx(iter),vy(iter),vz(iter),phi(iter),eps(iter),muc_array(iter),boundary(iter)
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

 nentryyz=nentry/nx

 allocate(XIM(ny))
 allocate(YIM(nz))
 allocate(sig(nentryyz))
 allocate(sig2(nentryyz))
 allocate(image(ny,nz))

 
 print *, '#nrntry, nentryyz, nx ny nz, nx ny, ', nentry, nentryyz, nx*ny*nz, nx*ny

 iter=0
 sig=0d0
 idx=SLICE
 ixs=0
 do idz=1,nz 
   do idy=1,ny
   iter=(idz-1)*nx*ny + nx*(idy-1)+idx
   ixs=ixs+1
   sig(ixs)=vy(iter)
   sig2(ixs)=vz(iter)
  enddo
 enddo 

 print *,"#", xmin,xmax,ymin,ymax,zmin,zmax
 do iter=1,ny
  XIM(iter)=( dble(iter-1) *dy -(dble(nx)/2d0-0.5d0)*dy)
 enddo
 do iter=1,nz
  YIM(iter)=( dble(iter-1) *dz -(dble(nz)/2d0-0.5d0)*dz)
 enddo

 do jpix=1,nz,2
  do ipix=1,ny,2
    ientry=(jpix-1)*ny+ipix
    vmag=sqrt(sig(ientry)**2+sig2(ientry)**2)
!    if(jpix==ny/2.and.ipix==nx/2)cycle
!    if(jpix==ny/2.and.ipix==nx/2+1)cycle
!    if(jpix==ny/2+1.and.ipix==nx/2)cycle
!    if(jpix==ny/2+1.and.ipix==nx/2+1)cycle
    if (vmag>vmax)then 
       vmag=vmag/vmax
    else
       vmag=1d0
    endif
    print *, XIM(ipix),YIM(jpix),sig(ientry)/vmag,sig2(ientry)/vmag
 enddo
enddo
 
 deallocate(x,y,z,rho,p,vx,vy,vz,phi,boundary)

end program

