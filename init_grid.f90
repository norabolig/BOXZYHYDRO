!
! Initialize the grid.  Most of this will not need to be modified.
! The only section that needs to be touched is if an anchor-cell 
! distribution is desired.  That will eventually go into a 
! a separate file.
!
subroutine init_grid
 use parameters
 use grid_commons
 implicit none
 
 integer:: igrid,iz,ibound,ix,iy,ibound2
 integer::flag1ix,flag1iy,flag1iz,flag2ix,flag2iy,flag2iz
 integer::flag1,flag2,flag,ineigh
 real(pre)::x,y,z,r,mslope,ael,bel,length,x0,aoa,xp,yp
 real(pre)::caf,paf,taf,maf,xc,xx,yt,yc,yu,yl,xu,xl,theta,dycdx,pc

 ngrid=nx*ny*nz
 print *,"# Total grid elements = ",ngrid

 allocate(grid(ngrid))
 allocate(p  (ngrid))
 allocate(adindx  (ngrid))
 allocate(muc_array  (ngrid))
 allocate(cons(5,ngrid))
 allocate(cons_old(5,ngrid))
 allocate(fluxtmp(5,ngrid))
 allocate(u  (3,ngrid))
 allocate(qq (3,ngrid))
 allocate(gforce (3,ngrid))
 allocate(pforce (3,ngrid))
 allocate(rhotot(ngrid))
 allocate(phi(ngrid))
 allocate(phi_new(ngrid))

 allocate(slope_p(3,ngrid))
 allocate(slope_d(3,ngrid))
 allocate(slope_e(3,ngrid))
 allocate(slope_g(3,ngrid))

 allocate(slope_u(3,3,ngrid))

 allocate(state_p_p(3,ngrid))
 allocate(state_d_p(3,ngrid))
 allocate(state_e_p(3,ngrid))
 allocate(state_g_p(3,ngrid))
 allocate(state_p_m(3,ngrid))
 allocate(state_d_m(3,ngrid))
 allocate(state_e_m(3,ngrid))
 allocate(state_g_m(3,ngrid))


 allocate(state_u_p(3,3,ngrid))
 allocate(state_u_m(3,3,ngrid))

 if(fluxangmom)then
   allocate(state_rmom_m(3,ngrid))
   allocate(state_amom_m(3,ngrid))
   allocate(state_rmom_p(3,ngrid))
   allocate(state_amom_p(3,ngrid))
   allocate(slope_rmom(3,ngrid))
   allocate(slope_amom(3,ngrid))
 endif

! allocate(f_cor(3,5,ngrid))

 call first_touch() ! below

 max_den_change_old=one
 max_den_old=zero

 nbound=0
 nghost=0
 nanchor=0

 if (nz<6)then
   no_outflow_zl=.false.
   no_outflow_zr=.false.
 endif

 igrid=0
 do iz=1,nz
  flag1iz=0; if(nz>=6) then; if ((iz==1).or.(iz==nz))flag1iz=1; endif
  flag2iz=0; if(nz>=6) then; if ((iz==2).or.(iz==nz-1))flag2iz=1; endif
  z=dz*(dble(iz)-(dble(nz)+one)/two)
  do iy=1,ny
    flag1iy=0; if ((iy==1).or.(iy==ny))flag1iy=1
    flag2iy=0; if ((iy==2).or.(iy==ny-1))flag2iy=1
    y=dy*(dble(iy)-(dble(ny)+one)/two)
    do ix=1,nx
      flag1ix=0; if ((ix==1).or.(ix==nx))flag1ix=1
      flag2ix=0; if ((ix==2).or.(ix==nx-1))flag2ix=1
      x=dx*(dble(ix)-(dble(nx)+one)/two)

      igrid=igrid+1
 
      grid(igrid)%boundary=0
      flag1=flag1iz+flag1ix+flag1iy
      flag2=flag2iz+flag2ix+flag2iy
      if(flag1>0)then
        grid(igrid)%boundary=1
        nbound=nbound+1
        nghost=nghost+1
      elseif(flag2>0)then
        grid(igrid)%boundary=2
        nghost=nghost+1
      endif
 
      grid(igrid)%id=igrid
      grid(igrid)%ix=ix
      grid(igrid)%iy=iy
      grid(igrid)%iz=iz
      grid(igrid)%x=x
      grid(igrid)%y=y
      grid(igrid)%z=z
   enddo
  enddo
 enddo


 allocate(indx_bound(nbound))
 print *, "Allocated ",nbound," boundary cells"
 allocate(indx_ghost(nghost))
 print *, "Allocated ",nghost," ghost cells"

!
!***
! The following is for setting anchors and obstructions on the grid.
! The example below is for a sphere.
!***
!
!
#ifdef EXTRAANCHORS 
  paf=object_radius  !planet radius

 do igrid=1,ngrid
   x=grid(igrid)%x
   y=grid(igrid)%y
   z=grid(igrid)%z

   flag=0
   ! set object geometry here
   if ( (x+object_x_displace)**2 + (y+object_y_displace)**2 + (z+object_z_displace)**2 < paf**2 )flag=1
     !( (x+object_x_displace)**2 + y*y + z*z < paf**2 )flag=1

    if (flag==1.and.grid(igrid)%boundary==0)then
      nanchor=nanchor+1
      grid(igrid)%boundary=3
    endif    
 enddo
 print *, nanchor,dble(nanchor)/dble(ngrid)
 allocate(indx_anchor(nanchor))

 ibound=0
 do igrid=1,ngrid
    if(grid(igrid)%boundary>2)then
        ibound=ibound+1
        indx_anchor(ibound)=igrid
    endif
 enddo
!
! 
#endif /* end ifdef EXTRAANCHORS */
!
!
 ibound=0
 ibound2=0
 do igrid=1,ngrid
   if(grid(igrid)%boundary==1)then
     ibound=ibound+1
     indx_bound(ibound)=igrid
     ibound2=ibound2+1
     indx_ghost(ibound2)=igrid
 
     ix=grid(igrid)%ix
     iy=grid(igrid)%iy
     iz=grid(igrid)%iz
     grid(igrid)%ineigh(:)=0

     if(ix==1)then ! always store the boundary donor in neighbor 1
        grid(igrid)%ineigh(1)=grid(igrid)%id+2
     elseif(ix==nx)then
        grid(igrid)%ineigh(1)=grid(igrid)%id-2
     elseif(iy==1)then
       grid(igrid)%ineigh(1)=grid(igrid)%id+2*nx
     elseif(iy==ny)then
       grid(igrid)%ineigh(1)=grid(igrid)%id-2*nx
     elseif(iz==1)then
      if(nz<6)then ! if it is 2D, then just make it so.
       grid(igrid)%ineigh(1)=grid(igrid)%id
      else
       grid(igrid)%ineigh(1)=grid(igrid)%id+2*nx*ny
      endif
     elseif(iz==nz)then
      if(nz<6)then ! if 2D, make it so.
       grid(igrid)%ineigh(1)=grid(igrid)%id
      else
       grid(igrid)%ineigh(1)=grid(igrid)%id-2*nx*ny
      endif
     else
       print *," Error. Grid boundary is not really a boundary."
       stop"Forced stop in init_grid"
     endif
   elseif(grid(igrid)%boundary==2)then
     ibound2=ibound2+1
     indx_ghost(ibound2)=igrid
     ix=grid(igrid)%ix
     iy=grid(igrid)%iy
     iz=grid(igrid)%iz
     ineigh=0
     grid(igrid)%ineigh(:)=ineigh

!     if(ix==2)then 
!       ineigh=grid(igrid)%id+1
!     elseif(ix==nx-1)then
!       ineigh=grid(igrid)%id-1
!     elseif(iy==2)then
!       ineigh=grid(igrid)%id+nx
!     elseif(iy==ny-1)then
!       ineigh=grid(igrid)%id-nx
!     elseif(iz==2)then
!      if(nz<6)then ! if 2D, make it so.
!       ineigh=grid(igrid)%id
!      else
!       ineigh=grid(igrid)%id+nx*ny
!      endif
!     elseif(iz==nz-1)then
!      if(nz<6)then ! if 2D, make it so.
!       ineigh=grid(igrid)%id
!      else
!!       ineigh=grid(igrid)%id-nx*ny
 !     endif
 !    endif

  !   grid(igrid)%ineigh(1)=ineigh

    grid(igrid)%ineigh(1)=grid(igrid)%id+nx
    grid(igrid)%ineigh(2)=grid(igrid)%id-nx
    grid(igrid)%ineigh(3)=grid(igrid)%id+1
    grid(igrid)%ineigh(4)=grid(igrid)%id-1
    if(nz<6)then ! if 2D, make it so.
      grid(igrid)%ineigh(5)=grid(igrid)%id
      grid(igrid)%ineigh(6)=grid(igrid)%id
    else
      grid(igrid)%ineigh(5)=grid(igrid)%id+nx*ny
      grid(igrid)%ineigh(6)=grid(igrid)%id-nx*ny
    endif

   else

    grid(igrid)%ineigh(1)=grid(igrid)%id+nx
    grid(igrid)%ineigh(2)=grid(igrid)%id-nx
    grid(igrid)%ineigh(3)=grid(igrid)%id+1
    grid(igrid)%ineigh(4)=grid(igrid)%id-1
    if(nz<6)then ! if 2D, make it so.
      grid(igrid)%ineigh(5)=grid(igrid)%id
      grid(igrid)%ineigh(6)=grid(igrid)%id
    else
      grid(igrid)%ineigh(5)=grid(igrid)%id+nx*ny
      grid(igrid)%ineigh(6)=grid(igrid)%id-nx*ny
    endif
   endif
 
 enddo
 end subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First touch principle.  This is for accelerating OpenMP by smartly
! allocating memory. As the name implies, the trick is to touch the
! arrays as soon as possible.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 subroutine first_touch()
  use parameters
  use grid_commons
  implicit none
  integer :: igrid

!$OMP PARALLEL 
!$OMP DO SCHEDULE(STATIC)
 do igrid=1,ngrid
   grid(igrid)%x=zero
   grid(igrid)%y=zero
   grid(igrid)%z=zero
   grid(igrid)%ineigh=zero
   grid(igrid)%id=zero
   grid(igrid)%ix=zero
   grid(igrid)%iy=zero
   grid(igrid)%iz=zero
   grid(igrid)%boundary=0
   phi(igrid)=zero
   p(igrid)=zero
   adindx(igrid)=zero
   muc_array(igrid)=zero
   cons(:,igrid)=zero
   cons_old(:,igrid)=zero
   u(:,igrid)=zero
   qq(:,igrid)=zero
   gforce(:,igrid)=zero
   pforce(:,igrid)=zero
   rhotot(igrid)=zero
   fluxtmp(:,igrid)=zero
 enddo
!$OMP ENDDO 
!$OMP END PARALLEL 

 end subroutine
 
