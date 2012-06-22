subroutine init_conditions()
 use parameters
 use derived_types
 use grid_commons
 use eos
 implicit none

 integer:: igrid,idim
 
 type(units)::scale

 real(pre)::x,y,z,r,zi,rtrope=4d0,ekin,angle
 real(pre)::eps,rho,tk,gam,eng
 real(pre)::rpoly=4.50,kpoly=1.934d14,rhopoly=4.4e-10,jnpoly=3.10e18,r2
 real(pre)::structure_gamma=1.435d0
 real(pre)::hfitd= 0.000331409, &
            gfitd= 2.64063e-05, &
            ffitd=-0.000204844, &
            efitd= 9.52241e-05, &
            dfitd= -1.71391e-05,& 
            cfitd= 1.10999e-06, &
            bfitd=zero, &
            afitd=zero


! real(pre)::afitd=0.0d0 , &
!            bfitd=0.0d0, &
!            cfitd=0.0d0,&
!            dfitd=-0.146658d00, &
!            efitd=0.332d0,&
!            ffitd=-0.21174d0, &
!            gfitd=-0.00206944d0, &
!            hfitd=0.0277743d0
!
 real(pre)::afitp=-0.0006826d0 , &
            bfitp=0.00395d0, &
            cfitp=-0.008202d0, &
            dfitp=0.0080d0, &
            efitp=-0.00376d0, &
            ffitp=0.00065d0


! fit is for r=1 and mass=10

 call get_units(scale)
! set your initial conditions here

!$OMP DO SCHEDULE(STATIC) PRIVATE(x,y,z,tk,eps) 
 do igrid=1,ngrid
  x=grid(igrid)%x;y=grid(igrid)%y;z=abs(grid(igrid)%z)
  phi(igrid)=zero
  adindx(igrid)=gammafix
  r=sqrt(x*x+y*y+z*z)
  r2=r
  angle=atan2(y,x)
  if (r < rpoly) then
     cons(1,igrid)= max(((((((afitd*r2+bfitd)*r2+cfitd)*r2+dfitd)*r2+efitd)&
                      *r2+ffitd)*r2+gfitd)*r2+hfitd,small_rho)
     u(1:3,igrid)=zero
  else
   cons(1,igrid)=small_rho
   cons(5,igrid)=small_eps
   u(1:3,igrid)=zero
  endif
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

 call state()

 print *, "#Done with ICs"

end subroutine 
   
