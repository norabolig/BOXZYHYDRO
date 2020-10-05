!
! Calculate the force on the fluid. Update the momenta and energy based on 
! these changes and divergences of flows. 
!
subroutine source(tphase,siter)
 use parameters
 use derived_types
 use grid_commons
 use utils, only : get_boundary, get_boundary_wb,get_polar,calc_gforce,&
     set_ghost_cells,get_boundary_noanchor
 use selfgravity, only : external_grav
!
!
#ifdef PARTICLE
 use particle, only : kick_particles,add_direct_togrid,zero_forces_on_direct, &
     zero_forces_on_indirect ! particle.f90
#endif
!
!
 implicit none

 integer :: igrid,siter
 integer :: b(6)
 logical::active
 real(pre)::x,y,r,ux,uy,uz,consold
 real(pre)::ekin,divphi,eps,tphase

 type(units)::scl
!
!
#ifdef ROTATE
!
!
  real(pre)::mbin,mbox,x0,omega2,omega
      !print *, "YOU MUST DEFINE THE CONDITIONS FOR THE ROTATING FRAME"
      !print *, "THIS CAN BE DONE IN SOURCE"
      !stop
      mbin=object_mass
      x0=object_x_displace
      omega2=(mbin)/abs(x0)**3
      omega=sqrt(omega2)
!
!
#endif /* end ifdef ROTATE */
!
!
 call get_units(scl) ! units.f90
!
!
#ifndef NOHYDRO
!
! 
!$OMP DO SCHEDULE(STATIC) &
!$OMP private(ekin,b) &
!$OMP private(divphi) &
!$OMP private(x,y,r) 

 do igrid=1,ngrid

  if(grid(igrid)%boundary>0)cycle
    call get_boundary_noanchor(igrid,b)

    x=grid(igrid)%x;y=grid(igrid)%y;r=sqrt(x*x+y*y)

    pforce(1,igrid)=-(p(b(3))-p(b(4)))/(cons(1,igrid)*two*dx)
    pforce(2,igrid)=-(p(b(1))-p(b(2)))/(cons(1,igrid)*two*dy)
    pforce(3,igrid)=-(p(b(5))-p(b(6)))/(cons(1,igrid)*two*dz)
!
!***
!Y direction
!***
!
    if(fluxangmom)then
      cons(3,igrid)=cons(3,igrid)+cons(1,igrid)*(x*u(2,igrid)-y*u(1,igrid))**2*y/r**4*dt
    endif
!
!
#ifdef ROTATE
      cons(3,igrid)=cons(3,igrid)+cons(1,igrid)*(-u(1,igrid)*omega)*dt*two &
                   +cons(1,igrid)*(1.5d0*omega*x/x0*x0dot)*dt
#endif
!
!***
! X DIRECTION
!***
!
    if(fluxangmom)then
      cons(2,igrid)=cons(2,igrid)+cons(1,igrid)*(x*u(2,igrid)-y*u(1,igrid))**2*x/r**4*dt
    endif
!
!
#ifdef ROTATE
      cons(2,igrid)=cons(2,igrid)+cons(1,igrid)*(u(2,igrid)*omega)*dt*two&
                   +cons(1,igrid)*(-1.5d0*omega*y/x0*x0dot)*dt
#endif
!
!
!***
! Z DIRECTION
!***
!
! DO NOTHING

 enddo
!$OMP ENDDO
!
!
#endif /* end ifdef NOHYDRO */
!
!
!
!
#ifdef PARTICLE
!
!
 call zero_forces_on_direct()
 call zero_forces_on_indirect()
!
!
#ifndef NOHYDRO
 call add_direct_togrid() ! particle.f90: kicked particles, so now add potential of particles for gas
#endif
!
! Kick called after add to grid in case gravity backreaction of gas is included.
!
 call kick_particles(dt) ! particle.f90
!
!
#endif /* end ifdef PARTICLE */
!
!
!$OMP MASTER
  rate_expand=zero
!$OMP END MASTER
!$OMP BARRIER
!
!
#ifndef NOHYDRO
!
!
#ifndef NOGRAVITY
!
!
!$OMP DO SCHEDULE(STATIC) REDUCTION(max:rate_expand) &
!$OMP private(ekin,b,active) &
!$OMP private(divphi,eps,ux,uy,uz)

 do igrid=1,ngrid

    if(grid(igrid)%boundary>0)cycle

    ux=u(1,igrid)
    uy=u(2,igrid)
    uz=u(3,igrid)

!
!***
!Y direction
!***
!
    cons(3,igrid)=cons(3,igrid)+(cons(1,igrid)*(gforce(2,igrid)))*dt

!
!*** 
! X
!***
!
    cons(2,igrid)=cons(2,igrid)+(cons(1,igrid)*(gforce(1,igrid)))*dt
!
!***
! Z DIRECTION
!***
!
    cons(4,igrid)=cons(4,igrid)+((gforce(3,igrid))*cons(1,igrid))*dt
 
    ekin = half*(cons(2,igrid)**2+cons(3,igrid)**2+cons(4,igrid)**2)/cons(1,igrid)

    consold=cons(5,igrid)
    cons(5,igrid)=max(cons(5,igrid)-ekin+dt*cons(1,igrid)*((&
                  (ux)*(gforce(1,igrid))+(uy)* &
                  (gforce(2,igrid))+&
                  (uz)*(gforce(3,igrid)))),small_eps)+ekin

!
!
#ifdef EXPANSION_LIMITED
    rate_expand=max(rate_expand,abs(consold-cons(5,igrid))/(cons(5,igrid)*dt))
#endif
!
!
 enddo
!$OMP ENDDO
!
!
#endif /* end ifndef NOGRAVITY */
!
!
!$OMP BARRIER
 call set_ghost_cells() ! utils.f90
!
!
#endif /* end ifndef NOHYDRO */
!
!
end subroutine

