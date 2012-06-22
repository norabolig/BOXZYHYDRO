subroutine init_conditions()
 use parameters
 use derived_types
 use grid_commons
 use eos
 implicit none

 integer:: igrid,idim
 
 type(units)::scale

 integer,parameter::ntble=62
 integer::itble,jtble,iter
 real(pre)::x,y,z,r,zi,ekin,angle,dx
 real(pre)::eps,rho,tk,gam,eng,val0,val1,jn
 real(pre)::rhotble(ntble,ntble),jntble(ntble,ntble),dum1,dum2,rtble(ntble)
 real(pre)::rpoly=3.6,kpoly=29d15,rhopoly=4.4e-10,jnpoly=3.10e18,r2
 real(pre)::afit=1.57365e-12 , &
            bfit=-8.31884e-12, &
            cfit=-1.76319e-12, &
            dfit=4.10389e-11

 call get_units(scale)
! set your initial conditions here
 dx=ds*cospi6
 print *, dx

!$OMP DO SCHEDULE(STATIC) PRIVATE(x,y,z) 
 do igrid=1,ngrid
  x=grid(igrid)%x;y=grid(igrid)%y;z=abs(grid(igrid)%z)
  phi(igrid)=zero
  adindx(igrid)=gammafix
  r=sqrt(x*x+y*y+z*z)
  r2=r
  angle=atan2(y,x)

  if (r < rpoly) then
     cons(1,igrid)= max((((afit*r2+bfit)*r2+cfit)*r2+dfit)/scale%density,small_rho)
     u(1:3,igrid)=zero
   !print *, r,rtble(itble+1),rtble(itble)
  else
   cons(1,igrid)=small_rho
   u(1:3,igrid)=zero
  endif
  ekin=half*(u(1,igrid)**2+u(2,igrid)**2+u(3,igrid)**2)*cons(1,igrid)
  p(igrid)=kpoly*(cons(1,igrid)*scale%density)**adindx(igrid)/scale%eps
  cons(5,igrid)=max(p(igrid)/(adindx(igrid)-one),small_eps)+ekin
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

 print *, "#Done with ICs"

end subroutine 
   
