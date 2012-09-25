!
! Set up initial conditions for the simulations.
! These ICs create a 4 AU poltyrope witn n=1.5
! and total mass 10 Mjup.  Modify as needed.
!
subroutine init_conditions()
 use parameters
 use derived_types
 use grid_commons
 use eos
 implicit none

 integer:: igrid,idim
 
 type(units)::scale

 real(pre)::x,y,z,r,zi,ekin,angle
 real(pre)::eps,rho,tk,gam,eng,r2
 real(pre)::rpoly=4.0,kpoly=8.055e13
 real(pre)::structure_gamma=1.4d0

 real(pre)::hfitd= 5.25d-10, &
            gfitd= -3.56d-12, &
            ffitd=-5.167d-10, &
            efitd= 3.3d-10, &
            dfitd= -7.89d-11,& 
            cfitd= 6.67d-12, &
            bfitd=zero, &
            afitd=zero

 call get_units(scale)
!
!***
! set your initial conditions here
!***
!
!$OMP DO SCHEDULE(STATIC) PRIVATE(x,y,z,tk,eps) 
 do igrid=1,ngrid
  x=grid(igrid)%x;y=grid(igrid)%y;z=abs(grid(igrid)%z)
  phi(igrid)=zero
  adindx(igrid)=gammafix
  r=sqrt(x*x+y*y+z*z)
  r2=r
  angle=atan2(y,x)
  if (r < rpoly) then
     cons(1,igrid)= max( (((((((afitd*r2+bfitd)*r2+cfitd)*r2+dfitd)*r2+efitd)&
                      *r2+ffitd)*r2+gfitd)*r2+hfitd)/scale%density,small_rho)
  else
   cons(1,igrid)=small_rho
  endif
  u(1:3,igrid)=zero
  ekin=half*(u(1,igrid)**2+u(2,igrid)**2+u(3,igrid)**2)*cons(1,igrid)
  p(igrid)=kpoly*(cons(1,igrid)*scale%density)**structure_gamma/scale%eps
  tk=p(igrid)/(cons(1,igrid)*scale%rgas)*muc
  call get_gamma_from_tk(eps,cons(1,igrid),tk,muc_array(igrid),adindx(igrid))
  cons(5,igrid)=eps+ekin
 
  cons(2,igrid)=cons(1,igrid)*u(1,igrid)
  cons(3,igrid)=cons(1,igrid)*u(2,igrid)
  cons(4,igrid)=cons(1,igrid)*u(3,igrid)

  do idim=1,5
    cons_old(idim,igrid)=cons(idim,igrid)
  enddo

 enddo
!$OMP ENDDO
!$OMP DO SCHEDULE(STATIC)
 do igrid=1,nbound
  cons(1,indx_bound(igrid))=small_rho
  cons(2,indx_bound(igrid))=zero 
  cons(3,indx_bound(igrid))=zero 
  cons(4,indx_bound(igrid))=zero 
  cons(5,indx_bound(igrid))=small_eps
 enddo
!$OMP ENDDO NOWAIT

 call state() ! state.f90

 print *, "#Done with ICs"

end subroutine 
   
