module coolinglaws
 use parameters
 use derived_types
 use grid_commons
 use utils
 implicit none

 real(pre)::scale_kappa

 contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 real(pre) function get_kappa(T)
   real(pre),intent(in)::T
   get_kappa=(T/64d0)**.5d0*scale_kappa
   !get_kappa=.33d0*scale_kappa
 end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine get_local_divflux()
   real(pre):: sigmaSBcode,dl,vol,opacfac,fluxfac,tirr,ebgrnd,dflux
   real(pre)::T,kappa,dtau,coolTime,lum
   type(units)::scl

   call get_units(scl)

   scale_kappa=scl%mass/scl%length**2
 
   lum=zero
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(dx,dy,dz,dt) 
!$OMP&SHARED(cons,p,muc,scl,ngrid) REDUCTION(+:lum)
!$OMP call get_units(scl)
   dl=sqrt(dx*dx+dy*dy+dz*dz)
   vol=dx,dy,dz
   opacfac=sqrt(3.)/2.
   fluxfac=8./sqrt(3.)
   sigmaSBcode=5.67d-5/(scl%mass/scl%time**3)
   Tirr=tk_bgrnd
   ebgrnd=tirr*scl%rgas/(muc)
!$OMP DO SCHEDULE(STATIC)
   do igrid=1,ngrid
     T=p(igrid)/(scl%rgas*cons(1,igrid))*muc
     kappa=get_kappa(T)
     dtau=kappa*cons(1,igrid)*dl
        
    ekin=half*(cons(2,igrid)**2+cons(3,igrid)**2+cons(4,igrid)**2)/cons(1,igrid)
    eps=max(cons(5,igrid)-ekin,small_eps)
 
     tau_fac=dtau+one/dtau
     dflux=fluxfac*sigmaSBcode*((T)**4-(tirr)**4)/(dl*tau_fac)
     coolTime = eps/abs(dflux)
     dflux=( eps-ebgrnd*cons(1,igrid) &
                        /(adindx(1)-one) )/dt &
                        *(one-exp(-(dt/coolTime)**2))  &
                        + dflux*exp(-(dt/coolTime)**2)
          
     lum=lum+dflux*vol
     cons(5,igrid)=max(cons(5,igrid)-dflux*dt,ekin)
   enddo
!$OMP ENDDO
!$OMP END PARALLEL
   print *, "Luminosity at time ",lum,time
   
  end subroutine
end module
  
 
