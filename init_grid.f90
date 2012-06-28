subroutine init_grid
 use parameters
 use grid_commons
 implicit none
 
 integer:: igrid,iz,ibound,ix,iy,ibound2
 integer::flag1ix,flag1iy,flag1iz,flag2ix,flag2iy,flag2iz
 integer::flag1,flag2,flag
 real(pre)::x,y,z,r

 ngrid=nx*ny*nz
 print *,"# Total grid elements = ",ngrid

 allocate(grid(ngrid))
 allocate(phi(ngrid))
 allocate(p  (ngrid))
 allocate(adindx  (ngrid))
 allocate(muc_array  (ngrid))
 allocate(cons(5,ngrid))
 allocate(cons_old(5,ngrid))
 allocate(cons_new(5,ngrid))
 allocate(u  (3,ngrid))
 allocate(qq (3,ngrid))
 allocate(gforce (3,ngrid))
 allocate(pforce (3,ngrid))
 allocate(rhotot(ngrid))

 call first_touch()

 max_den_change_old=one
 max_den_old=zero

 nbound=0
 nghost=0
 nanchor=0

 igrid=0
 do iz=1,nz
  flag1iz=0; if ((iz==1).or.(iz==nz))flag1iz=1
  flag2iz=0; if ((iz==2).or.(iz==nz-1))flag2iz=1
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

#ifdef EXTRAANCHORS ! create a bullet as an example
 anchorradius=5*dx

 nanchor=0
 do igrid=1,ngrid
    x=grid(igrid)%x-10*dx
   y=grid(igrid)%y
    z=grid(igrid)%z
    r=sqrt(x*x+y*y+z*z)
    flag=0
    if (r<anchorradius)then
      if (y>-anchorradius)flag=1
    else
      if (y>-anchorradius)then
        if (x>-20.*dx .and. x<0.)then
          if (y< anchorradius)then
              flag=1
          endif
        endif
      endif
    endif
    if (flag==1)then
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
 
#endif


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
       grid(igrid)%ineigh(1)=grid(igrid)%id+2*nx*ny
     elseif(iz==nz)then
       grid(igrid)%ineigh(1)=grid(igrid)%id-2*nx*ny
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
     grid(igrid)%ineigh(:)=0

     if(ix==2)then ! always store the boundary donor in neighbor 1
        grid(igrid)%ineigh(1)=grid(igrid)%id+1
     elseif(ix==nx-1)then
        grid(igrid)%ineigh(1)=grid(igrid)%id-1
     elseif(iy==2)then
       grid(igrid)%ineigh(1)=grid(igrid)%id+nx
     elseif(iy==ny-1)then
       grid(igrid)%ineigh(1)=grid(igrid)%id-nx
     elseif(iz==2)then
       grid(igrid)%ineigh(1)=grid(igrid)%id+nx*ny
     elseif(iz==nz-1)then
       grid(igrid)%ineigh(1)=grid(igrid)%id-nx*ny
     else
       print *," Error. Grid boundary is not really a boundary."
       stop"Forced stop in init_grid"
     endif

 
   else

    grid(igrid)%ineigh(1)=grid(igrid)%id+nx
    grid(igrid)%ineigh(2)=grid(igrid)%id-nx
    grid(igrid)%ineigh(3)=grid(igrid)%id+1
    grid(igrid)%ineigh(4)=grid(igrid)%id-1
    grid(igrid)%ineigh(5)=grid(igrid)%id+nx*ny
    grid(igrid)%ineigh(6)=grid(igrid)%id-nx*ny
   endif
   !if(grid(igrid)%boundary>0)print *, "BOUNDARY NEIGHBOR",igrid,ix,grid(igrid)%ineigh(1)
 
 enddo
 

 end subroutine

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
   cons_new(:,igrid)=zero
   cons_old(:,igrid)=zero
   u(:,igrid)=zero
   qq(:,igrid)=zero
   gforce(:,igrid)=zero
   pforce(:,igrid)=zero
   rhotot(igrid)=zero
 enddo
!$OMP ENDDO 
!$OMP END PARALLEL 

 end subroutine
 
