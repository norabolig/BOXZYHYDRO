subroutine init_grid
 use parameters
 use grid_commons
 implicit none
 
 integer:: igrid,iz,ibound,ix,iy
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

 igrid=1
 id_grid_center=0
 do iz=1,nz
  z=dz*(dble(iz)-(dble(nz)+one)/two)
  do iy=1,ny
    y=dy*(dble(iy)-(dble(ny)+one)/two)
    do ix=1,nx
      x=dx*(dble(ix)-(dble(nx)+one)/two)
 
      if(iz==1.or.iz==nz)then
        grid(igrid)%boundary=.true.
        grid(igrid)%anchor=.true.
        nbound=nbound+1
      elseif(iy==1.or.iy==ny)then
        grid(igrid)%boundary=.true.
        grid(igrid)%anchor=.true.
        nbound=nbound+1
      elseif(ix==1.or.ix==nx)then
        grid(igrid)%boundary=.true.
        grid(igrid)%anchor=.true.
        nbound=nbound+1
      else
        grid(igrid)%boundary=.false.
        grid(igrid)%anchor=.false.
      endif
      grid(igrid)%id=igrid
      grid(igrid)%ix=ix
      grid(igrid)%iy=iy
      grid(igrid)%iz=iz
      grid(igrid)%x=x
      grid(igrid)%y=y
      grid(igrid)%z=z
!      if(mod(igrid,10)==0)then
 !       grid(igrid)%anchor=.true.
 !     else
!        grid(igrid)%anchor=.false.
 !     endif
      igrid=igrid+1
   enddo
  enddo
 enddo

 allocate(indx_bound(nbound))
 print *, "Allocated ",nbound," boundary cells"



#ifdef EXTRAANCHORS
 anchorradius=max(nx*dx,ny*dy)
 anchorradius=max(anchorradius,nz*dz)*half

 nanchor=0
 do igrid=1,ngrid
    x=grid(igrid)%x
    y=grid(igrid)%y
    z=grid(igrid)%z
    r=sqrt(x*x+y*y+z*z)
    if (r>anchorradius)then
      nanchor=nanchor+1
      grid(igrid)%anchor=.true.
    endif    
 enddo
 print *, nanchor,dble(nanchor)/dble(ngrid)
 allocate(indx_anchor(nanchor))

 ibound=0
 do igrid=1,ngrid
    if(grid(igrid)%anchor)then
        ibound=ibound+1
        indx_anchor(ibound)=igrid
    endif
 enddo
 
#endif

 ibound=0
 do igrid=1,ngrid
   if(grid(igrid)%boundary)then
     ibound=ibound+1
     indx_bound(ibound)=igrid
     ix=grid(igrid)%ix
     iy=grid(igrid)%iy
     iz=grid(igrid)%iz
     grid(igrid)%ineigh(:)=0

     if(ix==1)then
       grid(igrid)%ineigh(1)=grid(igrid)%id+nx
       grid(igrid)%ineigh(2)=grid(igrid)%id-nx
       grid(igrid)%ineigh(4)=grid(igrid)%id-1
       if(iy==1)grid(igrid)%ineigh(2)=0
       if(iy==ny)grid(igrid)%ineigh(1)=0
     elseif(ix==nx)then
       grid(igrid)%ineigh(1)=grid(igrid)%id+nx
       grid(igrid)%ineigh(2)=grid(igrid)%id-nx
       grid(igrid)%ineigh(3)=grid(igrid)%id+1
       if(iy==1)grid(igrid)%ineigh(2)=0
       if(iy==ny)grid(igrid)%ineigh(1)=0
     elseif(iy==1)then
       grid(igrid)%ineigh(1)=grid(igrid)%id+nx
       grid(igrid)%ineigh(3)=grid(igrid)%id+1
       grid(igrid)%ineigh(4)=grid(igrid)%id-1
       if(ix==1)grid(igrid)%ineigh(3)=0
       if(ix==nx)grid(igrid)%ineigh(4)=0
     elseif(iy==ny)then
       grid(igrid)%ineigh(2)=grid(igrid)%id-nx
       grid(igrid)%ineigh(3)=grid(igrid)%id+1
       grid(igrid)%ineigh(4)=grid(igrid)%id-1
       if(ix==1)grid(igrid)%ineigh(3)=0
       if(ix==nx)grid(igrid)%ineigh(4)=0
     endif
     if(iz/=1)grid(igrid)%ineigh(5)=grid(igrid)%id+nx*ny
     if(iz/=nz)grid(igrid)%ineigh(6)=grid(igrid)%id-nx*ny
     cycle
   endif

   grid(igrid)%ineigh(1)=grid(igrid)%id+nx
   grid(igrid)%ineigh(2)=grid(igrid)%id-nx
   grid(igrid)%ineigh(3)=grid(igrid)%id+1
   grid(igrid)%ineigh(4)=grid(igrid)%id-1
   grid(igrid)%ineigh(5)=grid(igrid)%id+nx*ny
   grid(igrid)%ineigh(6)=grid(igrid)%id-nx*ny
 
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
   grid(igrid)%boundary=.false. 
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
 
