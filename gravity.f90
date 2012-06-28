module selfgravity
 use parameters
 use derived_types
 use grid_commons
 use utils
 implicit none

 type(gravcell),dimension(:),allocatable,save::grav_grid
 real(pre),dimension(:),allocatable,save::gdy,gdx,gdz,gravrho,phi_new
 real(pre)::mbox,xcom_grid,ycom_grid,zcom_grid
 real(pre)::abin=70.,ebin=0.0d0,mbin=1d0
 integer,dimension(:),allocatable,save::gny,gnx,gnz,ngrav_grid,grav_bound_indx,ngrav_bound
 integer,save::nchild=8,nlevel=4,nmulti=1

 real*8,save::thetaopen=0.8d0

 contains

 recursive subroutine factorial(n,fac)
  logical, save::first=.true.
  real(pre)::n,fac

  if(first)then
    fac=n
    if(fac==zero)fac=one
    if(n>one)first=.false.
  endif
  if((n-one)>one)then
    fac=fac*(n-one)
    call factorial(n-one,fac)
  else
    first=.true.
  endif
  return  
 end subroutine factorial

subroutine init_grav_grid
 
 integer::nleaf,ileaf,ilevel,ngrav_tot,iskip
 integer::igrid,iz,ix,iy
 integer::ilevel2,iskip2,inode,ichild,ibound,tbound
 real(pre)::x,y,z,xg,yg,zg,dist,boxx,boxy,boxz,yoff=zero

 allocate(phi_new(ngrid))
 allocate(gdx(nlevel))
 allocate(gdy(nlevel))
 allocate(gdz(nlevel))
 allocate(gnx(nlevel))
 allocate(gny(nlevel))
 allocate(gnz(nlevel))
 allocate(ngrav_grid(nlevel))
 allocate(ngrav_bound(nlevel))

 ngrav_tot=0

 boxx= dble(nx)*dx
 boxy=dble(ny)*dy
 boxz= dble(nz)*dz

! do while (dx*two**(nlevel)>boxx)
!   nlevel=nlevel-1
! enddo
! print *, "Gravity box too big. Reseting nlevel to be ",nlevel

 do ilevel=1,nlevel
   gdx(ilevel) = dx*two**(ilevel)
   gdy(ilevel) = dy*two**(ilevel)
   gdz(ilevel) = dz*two**(ilevel)
   gnx(ilevel) = boxx/gdx(ilevel)
   gny(ilevel) = boxy/gdy(ilevel)
   gnz(ilevel) = boxz/gdz(ilevel)

   ngrav_grid(ilevel) = gnx(ilevel)*gny(ilevel)*gnz(ilevel)
   ngrav_tot=ngrav_tot+ngrav_grid(ilevel)
   print *, "#",ilevel,ngrid,ngrav_grid(ilevel),gnx(ilevel),gny(ilevel),gnz(ilevel),gdx(ilevel),gdy(ilevel),gdz(ilevel)
 enddo
 allocate(grav_grid(ngrav_tot))
 
 ilevel=1
 nleaf=nchild

 iskip=0
 do ilevel =1,nlevel
 igrid=0

 do iz = 1,gnz(ilevel)
  z = gdz(ilevel)*(dble(iz)-(dble(gnz(ilevel))+one)*half)
  do iy = 1,gny(ilevel)
   y = gdy(ilevel)*(dble(iy)-(dble(gny(ilevel))+one)*half)
   do ix = 1,gnx(ilevel)
    x = gdx(ilevel)*(dble(ix)-(dble(gnx(ilevel))+one)*half)
    igrid=igrid+1
    grav_grid(igrid+iskip)%x=x
    grav_grid(igrid+iskip)%y=y
    grav_grid(igrid+iskip)%z=z
    grav_grid(igrid+iskip)%id=igrid+iskip
    grav_grid(igrid+iskip)%rmax=zero
    !print *, igrid,iskip,igrid+iskip,x,y,z,ilevel
    !if(ilevel==3)print *,x,y,z
    grav_grid(igrid+iskip)%boundary=.false.
    do ileaf=1,nchild
     grav_grid(igrid+iskip)%ichild(ileaf)=0
    enddo
   enddo
  enddo
 enddo


 iskip=iskip+ngrav_grid(ilevel)
 
 enddo

 iskip=0
 ilevel=1
 print *,"# Building links for grav level ",ilevel
  nleaf=0
  do igrid=1,ngrid
   x=grid(igrid)%x
   y=grid(igrid)%y
   z=grid(igrid)%z

    iz = int( (z+boxz*half)/gdz(ilevel) )+1
    iy = int( (y+boxy*half)/gdy(ilevel) )+1
    ix = int( (x+boxx*half)/gdx(ilevel) )+1
 
   inode = (iz-1)*gny(ilevel)*gnx(ilevel)+ (iy-1)*gnx(ilevel) +ix 
!   print *, ix,iy,iz,inode,ngrid,ngrav_grid(ilevel)
!   print *, x,y,z,grav_grid(inode+iskip)%x,grav_grid(inode+iskip)%y,grav_grid(inode+iskip)%z
!   print *, gdz(ilevel),gdy(ilevel),gdx(ilevel)
   xg=grav_grid(inode+iskip)%x
   yg=grav_grid(inode+iskip)%y
   zg=grav_grid(inode+iskip)%z

   if(abs(xg-x)*(one)>gdx(ilevel)*half .or.  &
      abs(yg-y)*(one)>gdy(ilevel)*half .or.  &
      abs(zg-z)*(one)>gdz(ilevel)*half ) then
      print *, "Mistmatch between node and leaf for node ",inode+iskip,igrid
      print *, xg-x,yg-y,zg-z,gdx(ilevel)*half,gdy(ilevel)*half,gdz(ilevel)*half
      print *, xg,yg,zg,x,y,z,ix,iy,iz
      stop
   endif

   dist=sqrt( (xg-x)**2+(yg-y)**2+(zg-z)**2 )
   if (dist>grav_grid(inode+iskip)%rmax) grav_grid(inode+iskip)%rmax=dist
 
   do ichild=1,nchild
     ileaf= grav_grid(inode+iskip)%ichild(ichild)
     if(ileaf==0)then
       grav_grid(inode+iskip)%ichild(ichild)=igrid
       if(grid(igrid)%boundary==1)then
          grav_grid(inode+iskip)%boundary=.true.
       endif
       nleaf=nleaf+1
       !print *, inode+iskip,xg,yg,zg,x,y,z,nleaf,igrid,ichild
       exit
     else
       if(ichild>nchild)then
         print *, grav_grid(inode+iskip)%ichild(ichild),ichild,nchild,ilevel
         print *,"More children than allowed.  Linking must be broken."
         print *, xg-x,yg-y,zg-z,gdx(ilevel)*half,gdy(ilevel)*half,gdz(ilevel)*half
         print *, xg,yg,zg,x,y,z
 
         stop
       endif
     endif
   enddo 
    !print *, ileaf,ileaf+iskip,xg,yg,zg,x,y,z,nleaf,igrid
  enddo
  !print *, "For grav_grid ",ileaf," found ",nleaf, " leaves. "
 print *, "Assigned a total of ",nleaf," children"
 !stop

 iskip=0
 iskip2=ngrav_grid(1)
 do ilevel=1,nlevel-1
  ilevel2=ilevel+1
  print *,"# Building links for grav level ",ilevel2
  do igrid=1,ngrav_grid(ilevel)
    x=grav_grid(igrid+iskip)%x
    y=grav_grid(igrid+iskip)%y
    z=grav_grid(igrid+iskip)%z
     do ix=1,gnx(ilevel2)
       xg = gdx(ilevel2)*(dble(ix)-one)-boxx*half+gdx(ilevel2)
       if ( x < xg ) exit
     enddo
     do iy=1,gny(ilevel2)
       yg = gdy(ilevel2)*(dble(iy)-one)-boxy*half+gdy(ilevel2)
       if ( y < yg ) exit
     enddo
     do iz=1,gnz(ilevel2)
       zg = gdz(ilevel2)*(dble(iz)-one)-boxz*half+gdz(ilevel2)
       if ( z < zg ) exit
     enddo
    inode = (iz-1)*gny(ilevel2)*gnx(ilevel2)+ (iy-1)*gnx(ilevel2) +ix  + iskip2
 
    xg=grav_grid(inode)%x
    yg=grav_grid(inode)%y
    zg=grav_grid(inode)%z
 
   if(abs(xg-x)>gdx(ilevel2)*half .and.  &
      abs(yg-y)>gdy(ilevel2)*half .and.  &
      abs(zg-z)>gdz(ilevel2)*half ) then
      print *, "Mistmatch between node and leaf"
      stop
   endif

   dist=sqrt( (xg-x)**2+(yg-y)**2+(zg-z)**2 )
   if (dist>grav_grid(inode)%rmax) grav_grid(inode)%rmax=dist
 
   do ichild=1,nchild
     ileaf= grav_grid(inode)%ichild(ichild)
     if(ileaf==0)then
       grav_grid(inode)%ichild(ichild)=igrid+iskip
       exit
     else
       if(ichild>nchild)then
         print *, inode,grav_grid(inode)%ichild(ichild),ichild,nchild,ilevel
         print *,"More children than allowed.  Linking must be broken."
         stop
       endif
     endif
   enddo 
    !print *, ileaf,ileaf+iskip,xg,yg,zg,nleaf,igrid
  enddo
  iskip=iskip+ngrav_grid(ilevel)
  iskip2=iskip2+ngrav_grid(ilevel2)
 enddo

 tbound=0
 do ilevel=1,nlevel
  iskip=0
  do inode=1,ilevel-1
    iskip=iskip+ngrav_grid(ilevel)
  enddo
  ibound=0
  do igrid=1,ngrav_grid(ilevel)
    inode=igrid+iskip
    if (grav_grid(inode)%boundary)ibound=ibound+1
  enddo
  ngrav_bound(ilevel)=ibound
  tbound=tbound+ibound
 enddo
 allocate(grav_bound_indx(tbound))
 ibound=0
 do ilevel=1,nlevel
  iskip=0
  do inode=1,ilevel-1
    iskip=iskip+ngrav_grid(ilevel)
  enddo
  do igrid=1,ngrav_grid(ilevel)
    inode=igrid+iskip
    if (grav_grid(inode)%boundary)then
      ibound=ibound+1
      grav_bound_indx(ibound)=inode
    endif
  enddo
 enddo
 

 print *, "# Built Gravity Tree "

 end subroutine

 subroutine set_com_in_tree
  real(pre)::x,y,z,mass,xcm,ycm,zcm,mtot_loc,mnode
  real(pre)::rmax
  integer::ilevel,iskip,inode,ichild,indx

  ilevel=1
  iskip=0
  mtot_loc=zero
  xcom_grid=zero;ycom_grid=zero;zcom_grid=zero
  do inode=1,ngrav_grid(ilevel)
    mnode=zero
    xcm=zero;ycm=zero;zcm=zero
   ! xg=grav_grid(inode+iskip)%x;yg=grav_grid(inode+iskip)%y;zg=grav_grid(inode+iskip)%z
    do ichild=1,nchild
     indx=grav_grid(inode+iskip)%ichild(ichild)
     if(indx==0)exit
     x=grid(indx)%x;y=grid(indx)%y;z=grid(indx)%z
     mass=rhotot(indx)*dx*dy*dz
     xcm=xcm+mass*(x)
     ycm=ycm+mass*(y)
     zcm=zcm+mass*(z)
     mnode=mnode+mass
    enddo
    xcom_grid=xcom_grid+xcm
    ycom_grid=ycom_grid+ycm
    zcom_grid=zcom_grid+zcm
    if(mnode==zero)cycle
    xcm=xcm/mnode;ycm=ycm/mnode;zcm=zcm/mnode
    rmax=zero
    do ichild=1,nchild
     indx=grav_grid(inode+iskip)%ichild(ichild)
     if(indx==0)exit
     x=grid(indx)%x;y=grid(indx)%y;z=grid(indx)%z
     rmax=max(rmax,sqrt( (x-xcm)**2+(y-ycm)**2+(z-zcm)**2))
    enddo
    mtot_loc=mtot_loc+mnode
    grav_grid(inode+iskip)%xcm=xcm
    grav_grid(inode+iskip)%ycm=ycm
    grav_grid(inode+iskip)%zcm=zcm
    grav_grid(inode+iskip)%mass=mnode
    grav_grid(inode+iskip)%rmax=rmax
    !print *, mnode,xcm,ycm,zcm,mnode
    !print *, mnode,xcm,ycm,zcm,rmax,sqrt(gdx(ilevel)**2+gdy(ilevel)**2+gdz(ilevel)**2)/2.
  enddo
#ifdef VERBOSE
  print *, "# Mtot on level 1 ",mtot_loc
#endif
  xcom_grid=xcom_grid/mtot_loc
  ycom_grid=ycom_grid/mtot_loc
  zcom_grid=zcom_grid/mtot_loc
  mbox=mtot_loc

  iskip = ngrav_grid(1)
  do ilevel=2,nlevel
   mtot_loc=zero
   do inode=1,ngrav_grid(ilevel)
     mnode=zero
     xcm=zero;ycm=zero;zcm=zero
     !xg=grav_grid(inode+iskip)%x;yg=grav_grid(inode+iskip)%y;zg=grav_grid(inode+iskip)%z
     do ichild=1,nchild
      indx=grav_grid(inode+iskip)%ichild(ichild)
      if(indx==0)exit
      x=grav_grid(indx)%xcm;y=grav_grid(indx)%ycm;z=grav_grid(indx)%zcm
      mass=grav_grid(indx)%mass
      xcm=xcm+mass*(x)
      ycm=ycm+mass*(y)
      zcm=zcm+mass*(z)
      mnode=mnode+mass
     enddo
     if(mnode==zero)cycle
     xcm=xcm/mnode;ycm=ycm/mnode;zcm=zcm/mnode
     rmax=zero
     do ichild=1,nchild
      indx=grav_grid(inode+iskip)%ichild(ichild)
      if(indx==0)cycle
      x=grav_grid(indx)%xcm;y=grav_grid(indx)%ycm;z=grav_grid(indx)%zcm
      rmax=max(rmax,sqrt( (x-xcm)**2+(y-ycm)**2+(z-zcm)**2))
      !if(ilevel==2)print *, mnode,xcm,ycm,zcm,rmax,sqrt(gdx(ilevel)**2+gdy(ilevel)**2+gdz(ilevel)**2)/2.
     enddo
     grav_grid(inode+iskip)%xcm=xcm
     grav_grid(inode+iskip)%ycm=ycm
     grav_grid(inode+iskip)%zcm=zcm
     grav_grid(inode+iskip)%mass=mnode
     grav_grid(inode+iskip)%rmax=rmax
     mtot_loc=mtot_loc+mnode
    !if(ilevel==2)print *, mnode,xcm/mnode,ycm/mnode,zcm/mnode,xg,yg,zg
   enddo
   iskip=iskip+ngrav_grid(ilevel)
   !print *,"#total mass on level ",ilevel,mtot_loc
  enddo

 return

end subroutine

subroutine get_pot_from_tree
  integer::igrid,ilevel,indx,l,flag,inode,iskip
  real(pre)::x,y,z,xx,yy,zz,r,phi_loc,theta


!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP&PRIVATE(igrid,ilevel,indx,l,flag,inode,iskip)&
!$OMP&PRIVATE(x,y,z,xx,yy,zz,r,phi_loc,theta)
!$OMP DO SCHEDULE(dynamic)
  do igrid=1,ngrid
   phi_loc=zero
   x=grid(igrid)%x;y=grid(igrid)%y;z=grid(igrid)%z
   iskip=0
   do ilevel=1,nlevel-1
     iskip=iskip+ngrav_grid(ilevel)
   enddo
   ilevel=nlevel
   do inode=1,ngrav_grid(ilevel)
    if(.not.(grav_grid(inode+iskip)%mass>zero))cycle
    xx=grav_grid(inode+iskip)%xcm;yy=grav_grid(inode+iskip)%ycm;zz=grav_grid(inode+iskip)%zcm
    r=sqrt( (x-xx)**2+(y-yy)**2+(z-zz)**2)
    flag=0
    if(r>zero)then
     theta = sqrt(gdx(ilevel)**2+gdy(ilevel)**2+gdz(ilevel)**2)/r !grav_grid(inode+iskip)%rmax/r
     !theta = 2.*grav_grid(inode+iskip)%rmax/r
!     if(grid(igrid)%iz==nz/2.and.grid(igrid)%irow==nrow/2)&
!       print *, ilevel,theta,grav_grid(inode+iskip)%rmax,sqrt(gdx(ilevel)**2+gdy(ilevel)**2+gdz(ilevel)**2),&
!                grav_grid(inode+iskip)%mass,r
     if(theta>thetaopen)flag=1
    else
     flag=1
    endif
    if(flag==0)then
      !print *, "Flag level, ",nlevel
      phi_loc=phi_loc-grav_grid(inode+iskip)%mass/r !- &
!             (grav_grid(inode+iskip)%dpole(1)*(x-xx)    + &
!              grav_grid(inode+iskip)%dpole(2)*(y-yy)    + &
!              grav_grid(inode+iskip)%dpole(3)*(z-zz))/r**3
    else
      indx=inode+iskip !grav_grid(inode)%id
      l=ilevel
      call gather_pot(phi_loc,x,y,z,indx,l,1)
    endif
   enddo
   !print *, igrid,phi_loc
   phi(igrid)=phi_loc
 enddo
!$OMP ENDDO
!$OMP END PARALLEL

end subroutine


subroutine get_pot_bc_from_tree
  integer::igrid,ilevel,indx,l,flag,inode,iskip,ibound
  real(pre)::x,y,z,xx,yy,zz,r,phi_loc,theta

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP&PRIVATE(igrid,ilevel,indx,l,flag,inode,iskip)&
!$OMP&PRIVATE(x,y,z,xx,yy,zz,r,phi_loc,theta)
!$OMP DO SCHEDULE(dynamic)
  do ibound=1,nbound
   igrid=indx_bound(ibound)
   !print *, "Entering grid object ",igrid 
   phi_loc=zero
   x=grid(igrid)%x;y=grid(igrid)%y;z=grid(igrid)%z
   iskip=0
   do ilevel=1,nlevel-1
     iskip=iskip+ngrav_grid(ilevel)
   enddo
   ilevel=nlevel
   do inode=1,ngrav_grid(ilevel)
    xx=grav_grid(inode+iskip)%xcm;yy=grav_grid(inode+iskip)%ycm;zz=grav_grid(inode+iskip)%zcm
    r=sqrt( (x-xx)**2+(y-yy)**2+(z-zz)**2)
    flag=0
    if(r>zero)then
     !theta = 2.*grav_grid(inode+iskip)%rmax/r
     theta = sqrt(gdx(ilevel)**2+gdy(ilevel)**2+gdz(ilevel)**2)/r !grav_grid(inode+iskip)%rmax/r
     if(theta>thetaopen)flag=1
    else
     flag=1
    endif
    if(flag==0)then
      !print *, "Flag level, ",nlevel
      phi_loc=phi_loc-grav_grid(inode+iskip)%mass/r !- &
!             (grav_grid(inode+iskip)%dpole(1)*(x-xx)    + &
!              grav_grid(inode+iskip)%dpole(2)*(y-yy)    + &
!              grav_grid(inode+iskip)%dpole(3)*(z-zz))/r**3
    else
      indx=inode+iskip !grav_grid(inode)%id
      l=ilevel
      call gather_pot(phi_loc,x,y,z,indx,l,1)
    endif
   enddo
   !print *, igrid,phi_loc
   phi(igrid)=phi_loc
 enddo
!$OMP ENDDO
#if 1 == 0
!$OMP DO SCHEDULE(dynamic)
  do igrid=1,ngrid,10
   if(grid(igrid)%boundary==1)cycle
   !print *, "Entering grid object ",igrid 
   phi_loc=zero
   x=grid(igrid)%x;y=grid(igrid)%y;z=grid(igrid)%z
   iskip=0
   do ilevel=1,nlevel-1
     iskip=iskip+ngrav_grid(ilevel)
   enddo
   ilevel=nlevel
   do inode=1,ngrav_grid(ilevel)
    xx=grav_grid(inode+iskip)%xcm;yy=grav_grid(inode+iskip)%ycm;zz=grav_grid(inode+iskip)%zcm
    r=sqrt( (x-xx)**2+(y-yy)**2+(z-zz)**2)
    flag=0
    if(r>zero)then
     !theta = 2.*grav_grid(inode+iskip)%rmax/r
     theta = sqrt(gdx(ilevel)**2+gdy(ilevel)**2+gdz(ilevel)**2)/r !grav_grid(inode+iskip)%rmax/r
     if(theta>thetaopen)flag=1
    else
     flag=1
    endif
    if(flag==0)then
      !print *, "Flag level, ",nlevel
      phi_loc=phi_loc-grav_grid(inode+iskip)%mass/r !- &
!             (grav_grid(inode+iskip)%dpole(1)*(x-xx)    + &
!              grav_grid(inode+iskip)%dpole(2)*(y-yy)    + &
!              grav_grid(inode+iskip)%dpole(3)*(z-zz))/r**3
    else
      indx=inode+iskip !grav_grid(inode)%id
      l=ilevel
      call gather_pot(phi_loc,x,y,z,indx,l,1)
    endif
   enddo
   !print *, igrid,phi_loc
   phi(igrid)=phi_loc
 enddo
!$OMP ENDDO
#endif

!$OMP END PARALLEL

end subroutine

subroutine get_pot_at_higher_level(mlevel)
  integer::igrid,ilevel,indx,l,flag,inode,iskip,mlevel,id,ibottom,iskip2
  real(pre)::x,y,z,xx,yy,zz,r,phi_loc,theta

  iskip2=0
   do ilevel=1,mlevel-1
     iskip2=iskip2+ngrav_grid(ilevel)
   enddo
 
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP&PRIVATE(id,x,y,z,inode,r,flag,ilevel,theta, &
!$OMP& xx,yy,zz,indx,l,phi_loc,ibottom,iskip)
!$OMP DO SCHEDULE(dynamic) 
  do ibottom=1,ngrav_grid(mlevel)
   id=ibottom+iskip2
   !print *, "Entering grid object ",igrid 
   phi_loc=zero
   x=grav_grid(id)%x;y=grav_grid(id)%y;z=grav_grid(id)%z
   iskip=0
   do ilevel=1,nlevel-1
     iskip=iskip+ngrav_grid(ilevel)
   enddo
   ilevel=nlevel
   do inode=1,ngrav_grid(ilevel)
    xx=grav_grid(inode+iskip)%xcm;yy=grav_grid(inode+iskip)%ycm;zz=grav_grid(inode+iskip)%zcm
    r=sqrt( (x-xx)**2+(y-yy)**2+(z-zz)**2)
    flag=0
    if(r>zero)then
     !theta = 2.*grav_grid(inode+iskip)%rmax/r
     theta = sqrt(gdx(ilevel)**2+gdy(ilevel)**2+gdz(ilevel)**2)/r !grav_grid(inode+iskip)%rmax/r
     if(theta>thetaopen)flag=1
    else
     flag=1
    endif
    if(flag==0)then
      !print *, "Flag level, ",nlevel
      phi_loc=phi_loc-grav_grid(inode+iskip)%mass/r !- &
!             (grav_grid(inode+iskip)%dpole(1)*(x-xx)    + &
!              grav_grid(inode+iskip)%dpole(2)*(y-yy)    + &
!              grav_grid(inode+iskip)%dpole(3)*(z-zz))/r**3
    else
      indx=inode+iskip !grav_grid(inode)%id
      l=ilevel
      call gather_pot(phi_loc,x,y,z,indx,l,mlevel)
    endif
   enddo
   !print *, igrid,phi_loc
   grav_grid(id)%phi=phi_loc
 enddo
!$OMP ENDDO
!$OMP ENDPARALLEL

end subroutine

subroutine get_pot_bc_at_higher_level(mlevel)
  integer::igrid,ilevel,indx,l,flag,inode,iskip,mlevel,id,ibottom,iskip2
  real(pre)::x,y,z,xx,yy,zz,r,phi_loc,theta

  iskip2=0
  do ilevel=1,mlevel-1
    iskip2=iskip2+ngrav_bound(ilevel)
  enddo
 
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP& PRIVATE(id,x,y,z,inode,r,flag,ilevel,theta, &
!$OMP& xx,yy,zz,indx,l,phi_loc,ibottom,iskip) 

!$OMP DO SCHEDULE(dynamic) 
  do ibottom=1,ngrav_bound(mlevel)
   id=grav_bound_indx(ibottom+iskip2)
   !print *, "Entering grid object ",igrid 
   phi_loc=zero
   x=grav_grid(id)%x;y=grav_grid(id)%y;z=grav_grid(id)%z
   iskip=0
   do ilevel=1,nlevel-1
     iskip=iskip+ngrav_grid(ilevel)
   enddo
   ilevel=nlevel
   do inode=1,ngrav_grid(ilevel)
    xx=grav_grid(inode+iskip)%xcm;yy=grav_grid(inode+iskip)%ycm;zz=grav_grid(inode+iskip)%zcm
    r=sqrt( (x-xx)**2+(y-yy)**2+(z-zz)**2)
    flag=0
    if(r>zero)then
     !theta = 2.*grav_grid(inode+iskip)%rmax/r
     theta = sqrt(gdx(ilevel)**2+gdy(ilevel)**2+gdz(ilevel)**2)/r !grav_grid(inode+iskip)%rmax/r
     if(theta>thetaopen)flag=1
    else
     flag=1
    endif
    if(flag==0)then
      !print *, "Flag level, ",nlevel
      phi_loc=phi_loc-grav_grid(inode+iskip)%mass/r !- &
!             (grav_grid(inode+iskip)%dpole(1)*(x-xx)    + &
!              grav_grid(inode+iskip)%dpole(2)*(y-yy)    + &
!              grav_grid(inode+iskip)%dpole(3)*(z-zz))/r**3
    else
      indx=inode+iskip !grav_grid(inode)%id
      l=ilevel
      call gather_pot(phi_loc,x,y,z,indx,l,mlevel)
    endif
   enddo
   !print *, igrid,phi_loc
   grav_grid(id)%phi=phi_loc
 enddo
!$OMP ENDDO
!$OMP END PARALLEL

end subroutine

recursive subroutine gather_pot(phi_loc,x,y,z,id,ilevel,level_limit)
  real(pre)::phi_loc,x,y,z,xx,yy,zz,theta,r
  integer::id,ilevel,flag,ilm,new_indx,ichild,indx,level_limit
  ilm=ilevel-1

  !print *, "Gathering on level ",ilevel
  if(ilevel>1)then

  do ichild=1,nchild
   indx=grav_grid(id)%ichild(ichild)
   if(indx==0)exit
   !Nprint *, id,indx,ilevel
   xx=grav_grid(indx)%xcm;yy=grav_grid(indx)%ycm;zz=grav_grid(indx)%zcm
   r=sqrt( (x-xx)**2+(y-yy)**2+(z-zz)**2)
   flag=0
   if(r>zero)then
    !theta = 2.*grav_grid(indx)%rmax/r
    theta = sqrt(gdx(ilm)**2+gdy(ilm)**2+gdz(ilm)**2)/r!grav_grid(indx)%rmax/r
    if(theta>thetaopen.and.ilevel>level_limit)flag=1
   else
    flag=1
   endif
   if(flag==0)then
     phi_loc=phi_loc-grav_grid(indx)%mass/r     !- &
!            (grav_grid(indx)%dpole(1)*(x-xx)    + &
!             grav_grid(indx)%dpole(2)*(y-yy)    + &
!             grav_grid(indx)%dpole(3)*(z-zz))/r**3

     !print *, phi_loc,grav_grid(indx)%mass,r,indx
     !print *, grav_grid(indx)%mass/r,(grav_grid(indx)%dpole(1)*(x-xx)    + &
     !        grav_grid(indx)%dpole(2)*(y-yy)    + &
     !        grav_grid(indx)%dpole(3)*(z-zz))/r**3
   else
     new_indx=indx!grav_grid(indx)%id
     !print *, new_indx,indx,ilevel,theta*180./pi
     call gather_pot(phi_loc,x,y,z,new_indx,ilm,level_limit)
     !print *,"AFTER",new_indx,indx,ilevel
   endif
  enddo
   
  else

  do ichild=1,nchild
    indx=grav_grid(id)%ichild(ichild)
    !print *, id,ilevel,indx,ichild
    if(indx==0)exit
    xx=grid(indx)%x;yy=grid(indx)%y;zz=grid(indx)%z
    r=sqrt( (x-xx)**2+(y-yy)**2+(z-zz)**2 )
    if(r>zero)then
     phi_loc=phi_loc-rhotot(indx)*dx*dy*dz/r
!    else
!     phi_loc=phi_loc-rhotot(indx)*dx*dy*dz/(quarter*sqrt(dy**2+dz**2+dx**2))
    endif
  enddo

  endif
end subroutine

subroutine get_boundary_grav(ilevel,id,b)
  integer::ilevel,id,b(6)

  b(1)=id+gnx(ilevel)
  b(2)=id-gnx(ilevel)
  b(3)=id+1
  b(4)=id-1
  b(5)=id+gnx(ilevel)*gny(ilevel)
  b(6)=id-gnx(ilevel)*gny(ilevel)

end subroutine

real(pre) function project_phi_from_coarse(ilevel,id,iskip,xc,yc,zc)
  integer::id,b(8),ilevel,zskip,iadd,iskip,i
  real(pre)::x,y,z,xc,yc,zc,phi1,phi2,phi3,phi4,phi5,phi6

  x=grav_grid(id)%x;y=grav_grid(id)%y;z=grav_grid(id)%z

  zskip=gnx(ilevel)*gny(ilevel)
  if (xc < x) then
    if (yc < y) then
      if ( zc < z ) then

       b(1) = id-1-gnx(ilevel)-zskip

      else

       b(1) = id-1-gnx(ilevel)
 
      endif

    else ! x negative y positive
      if ( zc < z ) then

       b(1) = id-1-zskip

      else

       b(1) = id-1
 
      endif

    endif
    
  else ! x positive

    if (yc < y) then
      if ( zc < z ) then

        b(1) = id - gnx(ilevel)-zskip

      else

        b(1) = id - gnx(ilevel)
 
      endif
    else ! x & y positive
      if ( zc < z ) then
        b(1) = id - zskip
      else
        b(1) = id 
      endif
    endif

  endif

  b(2) = b(1)+1
  b(3) = b(1)+gnx(ilevel)
  b(4) = b(3)+1
  b(5) = b(1)+zskip
  b(6) = b(2)+zskip
  b(7) = b(3)+zskip
  b(8) = b(4)+zskip

  iadd=0
  if(xc<grav_grid(iskip+1)%x)iadd=1
  if(yc<grav_grid(iskip+1)%y)iadd=iadd+gnx(ilevel)
  if(zc<grav_grid(iskip+1)%z)iadd=iadd+zskip

  if(xc>grav_grid(iskip+ngrav_grid(ilevel))%x)iadd=iadd-1
  if(yc>grav_grid(iskip+ngrav_grid(ilevel))%y)iadd=iadd-gnx(ilevel)
  if(zc>grav_grid(iskip+ngrav_grid(ilevel))%z)iadd=iadd-zskip

  do i=1,8
    b(i)=b(i)+iadd
  enddo

  x=grav_grid(b(1))%x
  y=grav_grid(b(1))%y
  z=grav_grid(b(1))%z
  phi1 = grav_grid(b(1))%phi+(grav_grid(b(2))%phi-grav_grid(b(1))%phi)*(xc-x)/gdx(ilevel)
  phi2 = grav_grid(b(3))%phi+(grav_grid(b(4))%phi-grav_grid(b(3))%phi)*(xc-x)/gdx(ilevel)
  phi3 = grav_grid(b(5))%phi+(grav_grid(b(6))%phi-grav_grid(b(5))%phi)*(xc-x)/gdx(ilevel)
  phi4 = grav_grid(b(7))%phi+(grav_grid(b(8))%phi-grav_grid(b(7))%phi)*(xc-x)/gdx(ilevel)

  phi5=phi1+(phi2-phi1)*(yc-y)/gdy(ilevel)
  phi6=phi3+(phi4-phi3)*(yc-y)/gdy(ilevel)

  project_phi_from_coarse=phi5+(phi6-phi5)*(zc-z)/gdz(ilevel)

end function

 subroutine walk_the_V(first)
  integer::igrid,b(6),iter,ilevel,inode,id,indx,iskip,ichild,ilower
  real(pre)::denom,alpha,err_loc,gradphi,resid,mloc
  real(pre)::dxloc,dyloc,dzloc,rho_loc,phi_loc
  real(pre)::phi0,phi1,phi2,phi3,phi4,phi5,phi6,rho0
  logical::active,first

  ilower=0
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP&PRIVATE(inode,id,indx,phi0,ichild,iskip,ilevel) &
!$OMP&PRIVATE(rho_loc,phi_loc)  

  ! baselevel
  if(first)then
    ilevel=1
!$OMP DO SCHEDULE(STATIC)
    do inode=1,ngrav_grid(ilevel)
       id=inode
       rho_loc=zero
       phi_loc=zero
       do ichild=1,nchild
         indx=grav_grid(id)%ichild(ichild)
         rho_loc=rho_loc+rhotot(indx)
         phi_loc=phi_loc+phi(indx)
       enddo
       grav_grid(id)%rho=rho_loc*0.125d0
       grav_grid(id)%phi=phi_loc*0.125d0
    enddo
!$OMP ENDDO
!$OMP BARRIER

    ! other levels
    do ilevel=2,nmulti
     iskip=0
     do inode=1,ilevel-1
       iskip=iskip+ngrav_grid(inode)
     enddo
!$OMP DO SCHEDULE(STATIC)
     do inode=1,ngrav_grid(ilevel)
       id=inode+iskip
       rho_loc=zero
       phi_loc=zero
       do ichild=1,nchild
         indx=grav_grid(id)%ichild(ichild)
         rho_loc=rho_loc+grav_grid(indx)%rho
         phi_loc=phi_loc+grav_grid(indx)%phi
       enddo
       grav_grid(id)%rho=rho_loc*0.125d0
       grav_grid(id)%phi=phi_loc*0.125d0
     enddo
!$OMP ENDDO
!$OMP BARRIER
    enddo
  else
   print *,"ERROR: Trying to do multiple V passes."
   print *,"Currently set for single V pass with cascading solution."
   stop "Stopping in walk_the_V"
  endif
!$OMP END PARALLEL

  do ilevel = nmulti,1,-1
    maxerr=zero
    !if(first)then
     !call get_pot_bc_at_higher_level(ilevel)
     if(ilevel==nmulti)then
         call get_pot_at_higher_level(ilevel)
     else
         call get_pot_bc_at_higher_level(ilevel)
     endif
    !endif
      maxerr=one
      iter=0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(err_loc,iter) 
      err_loc=maxerr
      do while (err_loc>grav_err_tol_low*half**ilevel)
!$OMP BARRIER
!$OMP MASTER
        maxerr=zero
!$OMP END MASTER
        iter=iter+1
        call sor_potential_coarse(ilevel,err_loc)
!$OMP ATOMIC
        maxerr=max(maxerr,err_loc)
!$OMP BARRIER
!$OMP MASTER
!        print *, ilevel,iter,err_loc,maxerr
!$OMP END MASTER
!$OMP BARRIER
       err_loc=maxerr
      enddo
!$OMP END PARALLEL
    iskip=0
    do inode=1,ilevel-1
      iskip=iskip+ngrav_grid(inode)
    enddo

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP&PRIVATE(inode,id,indx,phi0,ichild) &
!$OMP&PRIVATE(phi_loc)  
!$OMP DO SCHEDULE(STATIC)
    do inode=1,ngrav_grid(ilevel)
      id=inode+iskip
       !print *, ilevel,id,grav_grid(id)%x,grav_grid(id)%y,grav_grid(id)%z,grav_grid(id)%phi
      if(ilevel>1)then
       !  if(ilevel==nmulti)print *, ilevel,id,grav_grid(id)%x,grav_grid(id)%y,grav_grid(id)%z,grav_grid(id)%phi
        do ichild=1,nchild
         indx=grav_grid(id)%ichild(ichild)
         if(grav_grid(indx)%boundary)cycle
         phi0=project_phi_from_coarse(ilevel,id,iskip,grav_grid(indx)%x,grav_grid(indx)%y,grav_grid(indx)%z)
         grav_grid(indx)%phi=phi0
         !print *, ilevel,indx,grav_grid(indx)%x,grav_grid(indx)%y,grav_grid(indx)%z,phi0
        enddo
      else
        do ichild=1,nchild
         indx=grav_grid(id)%ichild(ichild)
         if(grid(indx)%boundary==1)cycle
         phi0=project_phi_from_coarse(ilevel,id,iskip,grid(indx)%x,grid(indx)%y,grid(indx)%z)
         phi(indx)=phi0
        enddo
      endif
    enddo
!$OMP ENDDO
!$OMP ENDPARALLEL
  enddo
  !stop

 end subroutine

 subroutine sor_potential_coarse(ilevel,err_loc)
  integer::igrid,b(6),iter,ilevel,inode,id,indx,iskip,ichild,redblack
  real(pre)::denom,alpha,err_loc,gradphi,resid,mloc
  real(pre)::dxloc,dyloc,dzloc,rho_loc,phi_loc
  real(pre)::phi0,phi1,phi2,phi3,phi4,phi5,phi6,rho0
  
   iskip=0
   do inode=1,ilevel-1
     iskip=iskip+ngrav_grid(inode)
   enddo
   err_loc=zero

!!!$OMP PARALLEL DEFAULT(SHARED) &
!!!$OMP&PRIVATE(inode,id,rho0,phi0,b) &
!!!$OMP&PRIVATE(phi1,phi2,phi3,phi4,phi5,phi6,resid) &
!!!!$OMP&REDUCTION(max:maxerr) 
   dxloc=gdx(ilevel)
   dyloc=gdy(ilevel)
   dzloc=gdz(ilevel)
   denom=(dyloc*dzloc/dxloc+dxloc*dzloc/dyloc+dxloc*dyloc/dzloc)*two
 
!$OMP DO SCHEDULE(STATIC) PRIVATE(inode,id,rho0,phi0,b) &
!$OMP&PRIVATE(phi1,phi2,phi3,phi4,phi5,phi6,resid) 
   do inode=1,ngrav_grid(ilevel)
     id=inode+iskip

     rho0=grav_grid(id)%rho
     phi0=grav_grid(id)%phi

     if(.not.grav_grid(id)%boundary)then
      call get_boundary_grav(ilevel,id,b)

      phi1=grav_grid(b(1))%phi
      phi2=grav_grid(b(2))%phi
      phi3=grav_grid(b(3))%phi
      phi4=grav_grid(b(4))%phi
      phi5=grav_grid(b(5))%phi
      phi6=grav_grid(b(6))%phi

      resid=( (phi1+phi2)*dxloc*dzloc/dyloc+(phi3+phi4)*dyloc*dzloc/dxloc &
           +(phi5+phi6)*dxloc*dyloc/dzloc &
           - phi0*denom-four*pi*rho0*dxloc*dyloc*dzloc )/(denom)

      grav_grid(id)%phi=phi0+resid*0.5d0
      
      if(phi0/=zero)then
        err_loc=max(err_loc,abs(resid/phi0))
      else
        err_loc=max(err_loc,abs(resid))
      endif
     endif
   enddo
!$OMP ENDDO
 end subroutine
   
 subroutine sor_potential_fine(maxiter)
  integer::igrid,b(6),iter,master_iter,maxiter
  real(pre)::denom,alpha,err_loc,gradphi,maxerr_loc,resid
  real(pre)::alpha0=1.0d0,dx2,dy2,dz2
  logical::active

   maxerr=zero
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP&private(b,iter,dx2,dy2,dz2) &
!$OMP&private(denom,alpha,err_loc,active,gradphi,maxerr_loc,resid)

  alpha=alpha0 
  maxerr=one
  iter=0
  denom=(dy*dz/dx+dx*dz/dy+dx*dy/dz)*two
  maxerr_loc=maxerr
!$OMP MASTER
  maxerr=zero
!$OMP ENDMASTER
!$OMP BARRIER
  do while (iter<maxiter)
!$OMP DO SCHEDULE(STATIC) REDUCTION(max:maxerr)
   do igrid=1,ngrid
     if(grid(igrid)%boundary==1)cycle
     call get_boundary_wb(igrid,b) 

      resid=( (phi(b(1))+phi(b(2)))*dx*dz/dy+(phi(b(3))+phi(b(4)))*dy*dz/dx &
           +(phi(b(5))+phi(b(6)))*dx*dy/dz &
           - phi(igrid)*denom-four*pi*rhotot(igrid)*dx*dy*dz )/(denom)

      phi_new(igrid)=phi(igrid)+resid*0.5d0

     err_loc = abs( resid/(phi_new(igrid)))
     maxerr=max(err_loc,maxerr)
 
   enddo
!$OMP ENDDO
!$OMP BARRIER
   maxerr_loc=maxerr
!$OMP DO  SCHEDULE(STATIC) 
   do igrid=1,ngrid
     if(grid(igrid)%boundary==1)cycle ! don't change the boundary or anchor  potential!
       phi(igrid)=phi_new(igrid)!phi(igrid)+((alpha)*phi_new(igrid)) ! subtract residual
   enddo
!$OMP ENDDO
#ifdef VERBOSE
!$OMP MASTER
  print *, "#",iter,maxerr
  master_iter=iter
!$OMP END MASTER
#endif
  iter=iter+1
  enddo

!$OMP END PARALLEL

 end subroutine

 subroutine vcycle_pot()
   integer::iter
   logical first

   iter=0
   maxerr=1d6
   first=.true.
   do while (maxerr>grav_err_tol.and.iter<1000)
     call walk_the_V(first)
     iter=iter+1
     first=.false.
   print *, "Walking the V (iter,maxerr): ",iter,maxerr
   enddo
   call sor_potential2(iter)
   print *, "Refining on fine grid (iter,maxerr): ",iter,maxerr

 end subroutine

 subroutine sor_potential2(master_iter)
  integer::igrid,b(6),iter,master_iter
  real(pre)::denom,alpha,err_loc,gradphi,maxerr_loc,resid
  real(pre)::alpha0=1.0d0,dx2,dy2,dz2
  logical::active

   maxerr=zero
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP&private(b,iter,dx2,dy2,dz2) &
!$OMP&private(denom,alpha,err_loc,active,gradphi,maxerr_loc,resid)

  alpha=alpha0 
  maxerr=one
  iter=0
  denom=(dy*dz/dx+dx*dz/dy+dx*dy/dz)*two
  maxerr_loc=maxerr
  do while(maxerr_loc>grav_err_tol.and.iter<50000.or.iter<miniter)
!$OMP MASTER
  maxerr=zero
!$OMP ENDMASTER
!$OMP BARRIER
!$OMP DO SCHEDULE(STATIC) REDUCTION(max:maxerr)
   do igrid=1,ngrid
     if(grid(igrid)%boundary==1)cycle
     call get_boundary_wb(igrid,b) 

      resid=( (phi(b(1))+phi(b(2)))*dx*dz/dy+(phi(b(3))+phi(b(4)))*dy*dz/dx &
           +(phi(b(5))+phi(b(6)))*dx*dy/dz &
           - phi(igrid)*denom-four*pi*rhotot(igrid)*dx*dy*dz )/(denom)

      phi_new(igrid)=phi(igrid)+resid*0.5d0

     err_loc = abs( resid/(phi(igrid)))
     maxerr=max(err_loc,maxerr)
 
   enddo
!$OMP ENDDO
!$OMP BARRIER
   maxerr_loc=maxerr
!$OMP DO  SCHEDULE(STATIC) 
   do igrid=1,ngrid
     if(grid(igrid)%boundary==1)cycle ! don't change the boundary or anchor  potential!
       phi(igrid)=phi_new(igrid)!phi(igrid)+((alpha)*phi_new(igrid)) ! subtract residual
   enddo
!$OMP ENDDO
   iter=iter+1
   if(mod(iter,1000)==0)then
     print *, "#WARNING SOR IS GOING THROUGH MANY ITERATIONS"
     print *, "#",iter,maxerr_loc,alpha
   endif
   alpha=alpha0 !- log10(min(one,maxerr))/ten
  enddo   
!$OMP BARRIER
!#ifdef VERBOSE
!$OMP MASTER
  !print *, "#",iter,maxerr
  master_iter=iter
!$OMP END MASTER
!#endif
!$OMP END PARALLEL

 end subroutine

 subroutine external_phi()! simple Roche potential
  integer::igrid
  real(pre)::x,y,z,rbin,rcylbin,omega2,pert
  real(pre)::x0,y0,z0,mbin,mbox
  real(pre)::xcom,ycom,zcom,plan,soft,r


!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(grid,phi,ngrid,dx)
  mbin=one
  mbox=30d-3/320.
  x0=-5.2d0;y0=0d0;z0=0d0
  xcom=x0-x0*mbox/(mbin+mbox);ycom=zero;zcom=zero
  omega2=(mbin+mbox)/abs(x0)**3
  soft=2*dx
!$OMP DO SCHEDULE(STATIC)
  do igrid=1,ngrid
    x=grid(igrid)%x
    y=grid(igrid)%y
    z=grid(igrid)%z
     
    rcylbin=sqrt( (x-x0)**2+(y-y0)**2)
    rbin = sqrt(rcylbin**2+(z-z0)**2 )
    r=sqrt(x*x+y*y+z*z)
    pert=-mbin/rbin-half*omega2*( (x-xcom)**2+(y-ycom)**2 )
    !print *, xcom,ycom,zcom,omega2,rbin,rcylbin,pert,phi(igrid)
    if (rbin<soft)then
      plan=-mbox/soft*( (r/soft)*(two + (r/soft)*( r/soft-two) ) )
    else
      plan=-mbox/r
    endif
    phi(igrid)=pert+plan
  enddo
!$OMP END DO
!$OMP END PARALLEL
 
 end subroutine
 
 subroutine rochepert()
  integer::igrid
  real(pre)::x,y,z,pert
  real(pre)::xb,yb,mbin_loc,rbin_off
  real(pre)::thistime,rbin,thetabin,omegabin,mumass,rb_loc

  !mbox=2.7508338628d-3
  !thistime=time-2900d0+1105d0
  !thistime=time+1105d0-4600d0
  !thistime=time ! taken from spin moving inward
  thistime=time+1110d0-2000d0
  call getBinaryPosition(thisTime,rBin,thetaBin,omegaBin,mbox,muMass,abin,ebin,mbin)
  rbin_off=rbin!sqrt( (rbin*cos(thetabin)-xcom_grid)**2+(rbin*sin(thetabin)-ycom_grid)**2 )
  print "(A,7(1X,1pe15.8))", "#binary position ",time,rbin_off,thetabin,xcom_grid,ycom_grid,&
        zcom_grid,mbox

!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(grid,phi,ngrid,rbin,thetabin,mbin) &
!$OMP& SHARED(rbin_off,xcom_grid,ycom_grid,zcom_grid)
  mbin_loc=mbin
  xb=rbin*cos(thetabin)!-xcom_grid
  yb=rbin*sin(thetabin)!-zcom_grid
  rb_loc=rbin_off
!$OMP DO SCHEDULE(STATIC)
  do igrid=1,ngrid
    x=grid(igrid)%x
    y=grid(igrid)%y
    z=grid(igrid)%z

    pert=-mbin_loc/sqrt( (x-xb)**2+(y-yb)**2+z**2) &
         +mbin_loc/rb_loc**3*(x*xb+y*yb)
     
    !print *, xb,yb,rb_loc,thetabin,pert,phi(igrid),mbin_loc
    phi(igrid)=phi(igrid)+pert
  enddo
!$OMP END DO
!$OMP END PARALLEL

  !abin=max(abin*(one-dt*omegaBin/(two*pi)),50d0) for hold 58, 50, 81
  !abin=abin*(one-dt*omegaBin/(eight*pi))
 
 end subroutine

 subroutine test_phi()
    integer::igrid,jgrid
    real(pre)::x,y,z,xx,yy,zz,testphi
    do igrid=1,ngrid
     x=grid(igrid)%x;y=grid(igrid)%y;z=grid(igrid)%z
    if(grid(igrid)%iz/=nz/2.or.grid(igrid)%iy/=ny/2)cycle
    testphi=zero
    do jgrid=1,ngrid
       xx=grid(jgrid)%x;yy=grid(jgrid)%y;zz=grid(jgrid)%z
       if(xx==x.and.yy==y.and.zz==z)then
         testphi=testphi-cons(1,jgrid)*dx*dy*dz/sqrt((quarter*dx)**2+(quarter*dz)**2+(quarter*dy)**2)
       else
         testphi=testphi-cons(1,jgrid)*dx*dy*dz/sqrt((x-xx)**2+(y-yy)**2+(z-zz)**2)
      endif
    enddo
    print *,grid(igrid)%x,grid(igrid)%y,grid(igrid)%z,cons(1,igrid),p(igrid),&
       u(1,igrid),u(2,igrid),u(3,igrid),phi(igrid),testphi,phi(igrid)/testphi,phi_new(igrid),&
       grid(igrid)%boundary
    enddo
    stop
 end subroutine

 

end module selfgravity

