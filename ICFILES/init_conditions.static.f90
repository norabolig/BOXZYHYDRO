subroutine init_conditions()
 use parameters
 use derived_types
 use grid_commons
 use eos
 implicit none

 integer:: igrid,idim
 
 type(units)::scale

 real(pre)::x,y,z,r,zi,kpoly,rtrope=12d0,ekin,angle
 real(pre)::eps,rho,tk,gam,eng

 
 call get_units(scale)
! set your initial conditions here

!$OMP DO SCHEDULE(STATIC) PRIVATE(x,y,z) 
 do igrid=1,ngrid
  x=grid(igrid)%x;y=grid(igrid)%y;z=grid(igrid)%z
  phi(igrid)=zero
  adindx(igrid)=two
  cons(1,igrid)=small_rho
  cons(2,igrid)=zero
  cons(3,igrid)=zero
  cons(4,igrid)=zero
  u(1,igrid)=zero
  u(2,igrid)=zero
  u(3,igrid)=zero
  if(abs(z)<ten)cons(1,igrid)=1d-6
   cons(2,igrid)=u(1,igrid)*cons(1,igrid)
   cons(3,igrid)=u(2,igrid)*cons(1,igrid)
   ekin=half*(u(1,igrid)**2+u(2,igrid)**2+u(3,igrid)**2)*cons(1,igrid)

   kpoly=four*pi*(rtrope/pi)**2*half
   p(igrid)= cons(1,igrid)*scl%rgas/muc*200d0 !kpoly*cons(1,igrid)**2
   cons(5,igrid)=max(p(igrid)/(adindx(igrid)-one),small_eps)+ekin
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

 print *, "#Done with ICs"

end subroutine 
   
