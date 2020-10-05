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

 real(pre)::ekin,rtrope=1d0,tk,rho,eps,x,y,z,kpoly,r
 real(pre)::hscale
 integer :: igrid

 type(units)::scale

 call get_units(scale)

#ifdef POLYEOS
!***
!
  kpoly=0.011d0**(one-gammafix)*300d0*scale%rgas/(muc)/gammafix**2/1.15
#endif /* end ifdef POLYEOS */
!
!
!
!$OMP DO SCHEDULE(STATIC) private(ekin,eps,tk,x,y,z,r)
 do igrid=1,ngrid
   x=grid(igrid)%x;y=grid(igrid)%y;z=grid(igrid)%z
   ekin=cons(1,igrid)*half*(u(1,igrid)**2+u(2,igrid)**2+u(3,igrid)**2)
!
!
#ifdef POLYEOS
  r=sqrt(x*x+y*y)
  p(igrid)=cons(1,igrid)*min(300d0/sqrt(r),1d4)*scale%rgas/muc
  !p(igrid)=(kpoly*cons(1,igrid)**gammafix)
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
    if (cons(1,igrid)<1d-4) then
        tk = min(5d3,p(igrid)/(cons(1,igrid)*scale%rgas)*muc)
        p(igrid)=(cons(1,igrid)*tk*scale%rgas/muc_array(igrid))
        eps = p(igrid)/(gammafix-one)
        cons(5,igrid)=eps+ekin
    else 
        p(igrid)=eps*(gammafix-one)
    endif
  else
    if(nz==1)then
      x=grid(igrid)%x;y=grid(igrid)%y
      hscale = sqrt( ( adindx(igrid)-1)*eps/(two*ekin) * (x*x+y*y) )
      eps=eps/(two*hscale)
      rho=cons(1,igrid)/(two*hscale)
!      if (cons(1,igrid)>1.)then
!        print *, rho*scale%density,hscale, hscale/sqrt(x*x+y*y)
!      endif
      call get_gamma2(eps,rho,tk,muc_array(igrid),adindx(igrid))
    else
       call get_gamma2(eps,cons(1,igrid),tk,muc_array(igrid),adindx(igrid))
    endif
    p(igrid)=(cons(1,igrid)*tk*scale%rgas/muc_array(igrid))
    if (cons(1,igrid)<1d-4)then
      p(igrid)=scale%rgas*tk_bgrnd/muc*cons(1,igrid)
      eps = p(igrid)/(gammafix-one)
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
