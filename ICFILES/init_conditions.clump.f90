subroutine init_conditions()
 use parameters
 use derived_types
 use grid_commons
 use eos
 implicit none

 integer:: igrid,idim
 
 type(units)::scale

 real(pre)::x,y,z,r,zi,kpoly,rtrope=1d0,ekin,angle
 real(pre)::eps,rho,tk,gam,eng

 
 call get_units(scale)
! set your initial conditions here

!$OMP DO SCHEDULE(STATIC) PRIVATE(x,y,z) 
 do igrid=1,ngrid
  x=grid(igrid)%x;y=grid(igrid)%y;z=grid(igrid)%z
  phi(igrid)=zero
  p(igrid)=1e-20
  adindx(igrid)=two
  cons(1,igrid)=1e-10
  cons(2,igrid)=zero
  cons(3,igrid)=zero
  cons(4,igrid)=zero
  cons(5,igrid)=max(p(igrid)/(adindx(igrid)-one),small_eps)
  u(1,igrid)=zero
  u(2,igrid)=zero
  u(3,igrid)=zero
!  if(x==zero.and.y==zero.and.z==zero)then
   r=sqrt(x*x+z*z+y*y)
  !if( (x/8.)**2+(y/8.)**2+(z/2.)**2<1.)then
  if(r<rtrope)then
   zi=pi*r/rtrope
   if(zi>zero)then
    angle=atan2(y,x)
    cons(1,igrid)=0.200e-5*sin(zi)/zi
   else
    cons(1,igrid)=0.200e-5!+cons(1,igrid)
   endif
   u(1,igrid)=-2.0e-3*y! was 2.5
   u(2,igrid)=2.0e-3*x ! was 2.5
   cons(2,igrid)=u(1,igrid)*cons(1,igrid)
   cons(3,igrid)=u(2,igrid)*cons(1,igrid)
  endif
   kpoly=four*pi*(rtrope/pi)**2*half
   ekin=half*(u(1,igrid)**2+u(2,igrid)**2+u(3,igrid)**2)*cons(1,igrid)

   kpoly=four*pi*(rtrope/pi)**2*half
   p(igrid)=kpoly*cons(1,igrid)**2
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
   
