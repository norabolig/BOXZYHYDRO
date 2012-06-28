subroutine source
 use parameters
 use derived_types
 use grid_commons
 use utils, only : get_boundary, get_boundary_wb,get_polar,calc_gforce,set_ghost_cells
#ifdef PARTICLE
 use particle
#endif
 implicit none

 integer :: igrid
 integer :: b(6)
 logical::active
 real(pre)::x,y,r,ux,uy,uz,consold
 real(pre)::ekin,divphi,eps


 type(units)::scl


#ifdef ROTATE
  real(pre)::mbin,mbox,x0,omega2,omega
      mbin=one
      mbox=30d-3/320.
      x0=-5.2d0
      omega2=(mbin+mbox)/abs(x0)**3
      omega=sqrt(omega2)
#endif

 call get_units(scl)

#ifndef NOHYDRO
 
!  call avisc()

!$OMP DO SCHEDULE(STATIC) &
!$OMP private(ekin,b,active) &
!$OMP private(divphi) &
!$OMP private(x,y,r) 

 do igrid=1,ngrid

  call get_boundary(igrid,b)
  if(grid(igrid)%boundary>0)cycle

    x=grid(igrid)%x;y=grid(igrid)%y+yoffset;r=sqrt(x*x+y*y)

    pforce(1,igrid)=-(p(b(3))-p(b(4)))/(cons(1,igrid)*two*dx)
    pforce(2,igrid)=-(p(b(1))-p(b(2)))/(cons(1,igrid)*two*dy)
    pforce(3,igrid)=-(p(b(5))-p(b(6)))/(cons(1,igrid)*two*dz)

    !Y direction
  
    cons(3,igrid)=cons(3,igrid)+pforce(2,igrid)*dt*cons(1,igrid)
    if(fluxangmom)then
      cons(3,igrid)=cons(3,igrid)+cons(1,igrid)*(x*u(2,igrid)-y*u(1,igrid))**2*y/r**4*dt
    endif
#ifdef ROTATE
      cons(3,igrid)=cons(3,igrid)+cons(1,igrid)*(-u(1,igrid)*omega)*dt*two
#endif

    !AV
!    cons(3,igrid)=cons(3,igrid)+cons(1,igrid)*(qq(2,b(1))-qq(2,b(2)))*dt/(two*dy)
 
    ! X

    cons(2,igrid)=cons(2,igrid)+pforce(1,igrid)*dt*cons(1,igrid)
    if(fluxangmom)then
      cons(2,igrid)=cons(2,igrid)+cons(1,igrid)*(x*u(2,igrid)-y*u(1,igrid))**2*x/r**4*dt
    endif
#ifdef ROTATE
      cons(2,igrid)=cons(2,igrid)+cons(1,igrid)*(u(2,igrid)*omega)*dt*two
#endif


    !AV
!    cons(2,igrid)=cons(2,igrid)+cons(1,igrid)*(qq(1,b(3))-qq(1,b(4)))*dt/(two*dx)

    ! Z DIRECTION

    cons(4,igrid)=cons(4,igrid)+pforce(3,igrid)*dt*cons(1,igrid)
   !AV
!    cons(4,igrid)=cons(4,igrid)-cons(1,igrid)*(qq(3,b(5))-qq(3,b(6)))*dt/(two*dz)

 enddo
!$OMP ENDDO

! call velocity()
#endif
#ifdef NOHYDRO
  if(npart>0.and.use_pic) call calc_gforce()
#else
  call calc_gforce()
#endif
!$OMP BARRIER
#ifdef PARTICLE
 call kick_particles(dt)
#ifndef NOHYDRO
 call add_direct_togrid()  ! kicked particles, so now add potential of particles for gas
#endif
#endif
!$OMP MASTER
  rate_expand=zero
!$OMP END MASTER
!$OMP BARRIER

#ifndef NOHYDRO
#ifndef NOGRAVITY
!$OMP DO SCHEDULE(STATIC) REDUCTION(max:rate_expand) &
!$OMP private(ekin,b,active) &
!$OMP private(divphi,eps,ux,uy,uz)

 do igrid=1,ngrid

    if(grid(igrid)%boundary>0)cycle

    ux=u(1,igrid)
    uy=u(2,igrid)
    uz=u(3,igrid)

!    ekin = half*(cons(2,igrid)**2+cons(3,igrid)**2+cons(4,igrid)**2)/cons(1,igrid)

    !Y direction
  
    cons(3,igrid)=cons(3,igrid)+(cons(1,igrid)*(gforce(2,igrid)))*dt

    !AV
 
    ! X

    cons(2,igrid)=cons(2,igrid)+(cons(1,igrid)*(gforce(1,igrid)))*dt

 
    !AV
 
    ! Z DIRECTION

    cons(4,igrid)=cons(4,igrid)+((gforce(3,igrid))*cons(1,igrid))*dt
   !AV

 
    ekin = half*(cons(2,igrid)**2+cons(3,igrid)**2+cons(4,igrid)**2)/cons(1,igrid)

    consold=cons(5,igrid)
    cons(5,igrid)=max(cons(5,igrid)-ekin+dt*cons(1,igrid)*((&
                  (ux)*(gforce(1,igrid))+(uy)* &
                  (gforce(2,igrid))+&
                  (uz)*(gforce(3,igrid)))),small_eps)+ekin

!    rate_expand=max(rate_expand,abs(consold-cons(5,igrid))/(cons(5,igrid)*dt))
 
 enddo
!$OMP ENDDO
#endif

 call set_ghost_cells()

#endif

end subroutine

