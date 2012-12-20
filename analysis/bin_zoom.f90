program sdxy
 implicit none

 integer,parameter::pre=8
 integer::nx=192,ny=192,nz=192,ix,iy,iz,idown,icheck=1234,b(8)
 integer::nx2=384,ny2=384,nz2=384,iter2,nentry2,idir,ninterp
 real(pre)::dy=0.025d0,dx=0.025d0,dz=0.025d0,r_vol=60d0,r,vol_ratio,w(8)
 real(pre)::dy2=0.0125d0,dx2=0.0125d0,dz2=0.0125d0,mtot,smallrho=1d-20,mtot0
 logical::use_offset=.false.

 character*72:: filename,input!="../celldump.00003750 "
 real(pre)::deltax,deltaz,h,lhex,xmin,xmax,ymin,ymax,zmin,zmax,pi,vol,dist,val,xx,yy,mass
 real(pre)::phase,area,val2,time,dmax,xoff,yoff,zoff,delx,dely,delz
 real(pre),dimension(:),allocatable::x,y,z,rho,p,vx,vy,vz,phi,sort
 real(pre),dimension(:),allocatable::x2,y2,z2,rho2,p2,vx2,vy2,vz2,phi2,muc2,eps2
 real(pre),dimension(:),allocatable::XIM,YIM,sig,eps,muc
 integer,dimension(:),allocatable::indx_sort,boundary,boundary2
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
 else
   open(unit=100,file=trim(filename))
   ioerr=0
   ihead=1
   nentry=0
   do while (ioerr>-1)
     if(ihead<3)then
        read(100,"(A1)")junk
        ihead=ihead+1
     endif
     read(100,"(A1)",iostat=ioerr)
     nentry=nentry+1
   enddo
   nentry=nentry-1
   print *, "#Found ",nentry," entries and ", ihead, " header entries"
   rewind(100)
 endif
 
 allocate(x(nentry))
 allocate(y(nentry))
 allocate(z(nentry))
 allocate(rho(nentry))
 allocate(p(nentry))
 allocate(vx(nentry))
 allocate(vy(nentry))
 allocate(vz(nentry))
 allocate(phi(nentry))
 allocate(muc(nentry))
 allocate(eps(nentry))
 allocate(boundary(nentry))

 nentry2=nx2*ny2*nz2
 
 allocate(x2(nentry2))
 allocate(y2(nentry2))
 allocate(z2(nentry2))
 allocate(rho2(nentry2))
 allocate(p2(nentry2))
 allocate(vx2(nentry2))
 allocate(vy2(nentry2))
 allocate(vz2(nentry2))
 allocate(phi2(nentry2))
 allocate(muc2(nentry2))
 allocate(eps2(nentry2))
 allocate(boundary2(nentry2))



 if(nentry/=ny*nx*nz)then
   print *, 'ordering off. check dimensions and entry ',nx*ny*nz,nentry," diff ",nx*ny*nz-nentry
   stop
 endif
 
 dmax=0d0
 mtot0=0d0
 if(binary)then
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
        vx(iter),vy(iter),vz(iter),phi(iter),eps(iter),muc(iter),boundary(iter)
     mtot0=mtot0+rho(iter)*dx*dy*dz
     if(rho(iter)>dmax)then
       dmax=rho(iter)
       xoff=x(iter)
       yoff=y(iter)
       zoff=z(iter)
     endif
     if(xmax<x(iter))xmax=x(iter)
     if(xmin>x(iter))xmin=x(iter)
     if(ymax<y(iter))ymax=y(iter)
     if(ymin>y(iter))ymin=y(iter)
     if(zmax<z(iter))zmax=z(iter)
     if(zmin>z(iter))zmin=z(iter)
   enddo
   enddo
   enddo
 else
   do iter=1,ihead-1
     read(100,"(A1)")junk
   enddo
   xmin=0.;xmax=0.;ymin=0.;ymax=0.;zmin=0.;zmax=0.
   iter=0
   do iz=1,nz
   do ix=1,nx
   do iy=1,ny
     iter=iter+1
     read(100,"(9(1pe16.8e3,1X),I2)")x(iter),y(iter),z(iter),rho(iter),p(iter),&
        vx(iter),vy(iter),vz(iter),phi(iter),eps(iter),muc(iter),boundary(iter)
     mtot0=mtot0+rho(iter)*dx*dy*dz
     if(rho(iter)>dmax)then
       dmax=rho(iter)
       xoff=x(iter)
       yoff=y(iter)
       zoff=z(iter)
     endif
     if(xmax<x(iter))xmax=x(iter)
     if(xmin>x(iter))xmin=x(iter)
     if(ymax<y(iter))ymax=y(iter)
     if(ymin>y(iter))ymin=y(iter)
     if(zmax<z(iter))zmax=z(iter)
     if(zmin>z(iter))zmin=z(iter)
   enddo
   enddo
   enddo
 endif
 close(100)

 x2=0d0;y2=0d0;z2=0d0;rho2=0d0;p2=0d0;vx2=0d0;vy2=0d0;vz2=0d0;phi2=0d0;muc2=0d0;eps2=0d0;boundary2=0d0

 iter=1
 do iz=1,nz2
  do iy=1,ny2
   do ix=1,nx2
    z2(iter)=dz2*(dble(iz)-(dble(nz2)+1d0)/2d0)
    y2(iter)=dy2*(dble(iy)-(dble(ny2)+1d0)/2d0)
    x2(iter)=dx2*(dble(ix)-(dble(nx2)+1d0)/2d0)
    iter=iter+1
   enddo
  enddo
 enddo

 vol_ratio=1d0!dx*dy*dz/(dx2*dy2*dz2)/dble(ninterp)**3
 do iter2=1,nentry2

    iter=get_grid_indx(x2(iter2),y2(iter2),z2(iter2),dx,dy,dz,nx,ny,nz)
    if(boundary(iter)>0)cycle
    delx=-x(iter)+x2(iter2)
    dely=-y(iter)+y2(iter2)
    delz=-z(iter)+z2(iter2)
    if(delx<0d0)iter=iter-1
    if(dely<0d0)iter=iter-nx
    if(delz<0d0)iter=iter-nx*ny
    b(1)=iter
    b(2)=iter+1
    b(3)=iter+1+nx
    b(4)=iter+nx
    b(5)=b(1)+nx*ny
    b(6)=b(2)+nx*ny
    b(7)=b(3)+nx*ny
    b(8)=b(4)+nx*ny

    call get_weights_8(b,x2(iter2),y2(iter2),z2(iter2),x,y,z,nentry,w)
  
    do iter=1,8
     rho2(iter2)=rho2(iter2)+rho(b(iter))*vol_ratio*w(iter)
     p2(iter2)=p2(iter2)+p(b(iter))*vol_ratio*w(iter)
     vx2(iter2)=vx2(iter2)+vx(b(iter))*rho(b(iter))*vol_ratio*w(iter)
     vy2(iter2)=vy2(iter2)+vy(b(iter))*rho(b(iter))*vol_ratio*w(iter)
     vz2(iter2)=vz2(iter2)+vz(b(iter))*rho(b(iter))*vol_ratio*w(iter)
     phi2(iter2)=phi2(iter2)+phi(b(iter))*vol_ratio*w(iter)
     muc2(iter2)=muc2(iter2)+muc(b(iter))*vol_ratio*w(iter)
     eps2(iter2)=eps2(iter2)+eps(b(iter))*vol_ratio*w(iter)
   enddo
 enddo
 mtot=0d0
 do iter2=1,nentry2
    mtot=mtot+rho2(iter2)*dx2*dy2*dz2
    if(rho2(iter2)>=smallrho)then
     vx2(iter2)=vx2(iter2)/rho2(iter2)
     vy2(iter2)=vy2(iter2)/rho2(iter2)
     vz2(iter2)=vz2(iter2)/rho2(iter2)
    else
     vx2(iter2)=0d0
     vy2(iter2)=0d0
     vz2(iter2)=0d0
    endif
    rho2(iter2)=max(rho2(iter2),smallrho)
 enddo

 iter=0
 open(unit=101,file='celldump.interpolated',form="UNFORMATTED")
 write(101)ihead
 write(101)time,ihead
 write(101)scl
 do iz=1,nz2
 do iy=1,ny2
 do ix=1,nx2
     iter=iter+1
     write(101)x2(iter),y2(iter),z2(iter),rho2(iter),p2(iter),&
        vx2(iter),vy2(iter),vz2(iter),phi2(iter),eps2(iter),muc2(iter),boundary2(iter)
 enddo;enddo;enddo
 close(101)

 print *, '#time ',time
 print *, '#Total mass (initial and final) ',mtot0,mtot

 deallocate(x,y,z,rho,p,vx,vy,vz,phi,boundary)
 deallocate(x2,y2,z2,rho2,p2,vx2,vy2,vz2,phi2,boundary2)

 contains

 integer function get_grid_indx(x,y,z,dx,dy,dz,nx,ny,nz)
   real(pre)::x,y,z,dx,dy,dz,half=0.5d0
   integer::ix,igrid,iz,iy,nx,ny,nz

   iz=int((z+half*dz)/dz+dble(nz+1)*half)
   if(iz>nz)iz=-1
   if(iz<1)then; igrid=-1;return;endif

   igrid=(iz-1)*nx*ny
   
   iy=int((y+half*dy)/dy+dble(ny+1)*half)
   if(iy>ny)iy=-1
   if(iy<1)then; igrid=-1;return;endif

   ix=int((x+half*dx)/dx+dble(nx+1)*half)
   if(ix>nx)ix=-1
   if(ix<1)then; igrid=-1;return;endif

   igrid=igrid+nx*(iy-1)+ix
   get_grid_indx=igrid

 end function
!
 subroutine get_weights_8(ig,x,y,z,xg,yg,zg,ng,w)
  integer,intent(in)::ng
  real(pre),intent(in)::x,y,z,xg(ng),yg(ng),zg(ng)
  real(pre),intent(out)::w(8)
  real(pre)::delx1,delx2,delx3,delx4,dely1,dely2,dely3,dely4,area1,area2,area3,area4
  real(pre)::weight
  integer,intent(out)::ig(8)
  integer::iter,flag
    delx1=abs(xg(ig(1))-x)
    dely1=abs(yg(ig(1))-y)
    delx2=abs(xg(ig(2))-x)
    dely2=abs(yg(ig(2))-y)
    delx3=abs(xg(ig(3))-x)
    dely3=abs(yg(ig(3))-y)
    delx4=abs(xg(ig(4))-x)
    dely4=abs(yg(ig(4))-y)

    area3=delx1*dely1
    area4=delx2*dely2
    area1=delx3*dely3
    area2=delx4*dely4
 
    w(1)=area1*(zg(ig(5))-z)
    w(2)=area2*(zg(ig(5))-z)
    w(3)=area3*(zg(ig(5))-z)
    w(4)=area4*(zg(ig(5))-z)
    w(5)=area1*(z-zg(ig(1)))
    w(6)=area2*(z-zg(ig(1)))
    w(7)=area3*(z-zg(ig(1)))
    w(8)=area4*(z-zg(ig(1)))

    weight=0d0
    flag=0
    do iter=1,8
     weight=weight+w(iter)
     if(w(iter)<0)flag=1
    enddo
    if (weight<=0d0)flag=1
    if(flag==1)then
      ig=-1
      w=0d0
      return
    endif
 
    do iter=1,8
     w(iter)=w(iter)/weight
    enddo

 end subroutine
 
end program
