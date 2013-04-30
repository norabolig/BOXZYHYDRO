!
! Routine for relating energy temperature and pressure
! It also allows simple EOSs to be used. POLY is set
! as an example for a polytropic EOS.
!
subroutine state
 use parameters
 use derived_types
 use grid_commons
 use eos
 implicit none

 real(pre)::ekin,rtrope=1d0,coolfac,tk,eps,x,y,z,kpolycgs=3.61352628338d+17,kpoly
 integer :: igrid

 type(units)::scale

 call get_units(scale)

#ifdef POLYEOS
 kpoly=four*pi*(rtrope/pi)**2*half
!
!***
! The following should be kept for an example of a slowing cooling Polytropic EOS
!***
!
! kpoly=29d15*scale%density**gammafix/scale%eps - &
!   (gammafix-one)*1.67d-3**(one-gammafix)*1.9e-1/scale%vel**2*scale%time*time ! scales to 1d-9 g/cc. ! 1e-4 Lsun
!!$OMP MASTER
!  print *,"POLY",time,kpoly
!!$OMP END MASTER
!
!
#endif /* end ifdef POLYEOS */
!
!
!$OMP DO SCHEDULE(STATIC) private(ekin,eps,tk)
 do igrid=1,ngrid
   x=grid(igrid)%x;y=grid(igrid)%y;z=grid(igrid)%z
   ekin=cons(1,igrid)*half*(u(1,igrid)**2+u(2,igrid)**2+u(3,igrid)**2)
!
!
#ifdef POLYEOS
  p(igrid)=(kpoly*cons(1,igrid)**gammafix)
  adindx(igrid)=gammafix
  muc_array(igrid)=muc
  cons(5,igrid)=p(igrid)/(adindx(igrid)-one)+ekin 
!
!
#else
!
!
  ekin=cons(1,igrid)*half*(u(1,igrid)**2+u(2,igrid)**2+u(3,igrid)**2)
  eps=max(cons(5,igrid)-ekin,small_eps)
  if(H2STAT==-1)then
    adindx(igrid)=gammafix
    muc_array(igrid)=muc
    p(igrid)=eps*(gammafix-one)
  else
    call get_gamma2(eps,cons(1,igrid),tk,muc_array(igrid),adindx(igrid))
    p(igrid)=max((cons(1,igrid)*tk*scale%rgas/muc_array(igrid)),scale%rgas*tk_bgrnd/muc*cons(1,igrid))
    if (tk<tk_bgrnd)then
      call get_gamma_from_p(eps,cons(1,igrid),p(igrid),muc_array(igrid),adindx(igrid))
      cons(5,igrid)=eps+ekin
    endif
  endif
!
!
#endif /* end endif POLYEOS */
!
!

 enddo
!$OMP ENDDO

end subroutine
