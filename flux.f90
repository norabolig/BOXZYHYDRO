!
! Fluxing routine that uses centeral differencing or upwind schemes
! At the moment, variables f2,f1,f0,fm1 are calculated but not used
! This will likely be removed, but it may stick around for a bit
!
subroutine flux(avg)
 use parameters
 use derived_types
 use grid_commons
 use utils
 implicit none


 integer :: igrid,b(6),idim,isten(5),avg,flag
 real(pre)::fx(2,5),fy(2,5),fz(2,5),areaxy,areayz,areaxz
 real(pre)::vol
 real(pre)::x,y,z,r,angmom,rmom

 type(units)::scale

 call get_units(scale)

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP&PRIVATE(areaxy,areayz,areaxz,vol) &
!$OMP&private(fz,fx,fy,flag) &
!$OMP&private(b,angmom,rmom,x,y,z,r,isten)

 areaxy=dy*dx
 areayz=dy*dz
 areaxz=dx*dz
 vol=dx*dy*dz

!
!***
! Will be solving for the total energy, excluding gravity. 
! Add pressure to this quantity.
!***
!
!
!***
! If angular momentum fluxing is used, make a variable switch.
!***
!
 if(fluxangmom)then
!$OMP DO SCHEDULE(STATIC)
  do igrid=1,ngrid
     x=grid(igrid)%x;y=grid(igrid)%y+yoffset
     r=sqrt(x*x+y*y)
     angmom=x*cons(3,igrid)-y*cons(2,igrid)
     rmom=cons(2,igrid)*x/r+cons(3,igrid)*y/r
     cons(2,igrid)=rmom
     cons(3,igrid)=angmom
  enddo
!$OMP ENDDO
 endif
!
!***
! Enter main loop
!***
!
! First get slopes
!
!$OMP DO SCHEDULE(static)
 do igrid=1,ngrid

  call clear_slope(igrid)
  fluxtmp(:,igrid)=zero

  if(grid(igrid)%boundary>0)then
      b(:)=igrid
  else
      call get_boundary_wb(igrid,b)
  endif

  isten(1)  = b(4)
  isten(2)  = igrid
  isten(3)  = b(3)
  isten(5)  = 1

  call calculate_slopes(isten) ! only indices 1:3 matter here

  isten(1)  = b(2)
  isten(2)  = igrid
  isten(3)  = b(1)
  isten(5)  = 2

  call calculate_slopes(isten) ! only indices 1:3 matter here

  isten(1)  = b(6)
  isten(2)  = igrid
  isten(3)  = b(5)
  isten(5)  = 3

  call calculate_slopes(isten) ! only indices 1:3 matter here

 enddo
!$OMP ENDDO 
!$OMP BARRIER
!
!
! next construct left and right states
!
!$OMP DO SCHEDULE(static)
 do igrid=1,ngrid

  if(grid(igrid)%boundary>0)then
      b(:)=igrid
  else
      call get_boundary_wb(igrid,b)
  endif


  isten(1)  = b(4)
  isten(2)  = igrid
  isten(3)  = b(3)

  call calculate_states(isten) ! only indices 1:3 matter here

  isten(1)  = b(2)
  isten(2)  = igrid
  isten(3)  = b(1)

  call calculate_states(isten) ! only indices 1:3 matter here

  isten(1)  = b(6)
  isten(2)  = igrid
  isten(3)  = b(5)

  call calculate_states(isten) ! only indices 1:3 matter here

 enddo
!$OMP ENDDO 
!$OMP BARRIER
!
! next calculate fluxes
!
!$OMP DO SCHEDULE(static)
 do igrid=1,ngrid

  if(grid(igrid)%boundary>0)cycle

  call get_boundary_wb(igrid,b)

  isten(1)  = b(4)
  isten(2)  = igrid
  isten(3)  = b(3)
  isten(5)  = 1

  call get_fluxes(isten,fx) 
  flag=0
  call adjust_for_boundary(isten,3,flag)
  if(flag==1)then
    fx(2,:)=zero
    fx(2,2)=p(igrid)
  endif
  flag=0
  call adjust_for_boundary(isten,1,flag)
  if(flag==1)then
    fx(1,:)=zero
    fx(1,2)=p(igrid)
  endif

  isten(1)  = b(2)
  isten(2)  = igrid
  isten(3)  = b(1)
  isten(5)  = 2

  call get_fluxes(isten,fy) 
  flag=0
  call adjust_for_boundary(isten,3,flag)
  if(flag==1)then
    fy(2,:)=zero
    fy(2,3)=p(igrid)
  endif
  flag=0
  call adjust_for_boundary(isten,1,flag)
  if(flag==1)then
    fy(1,:)=zero
    fy(1,3)=p(igrid)
  endif


  isten(1)  = b(6)
  isten(2)  = igrid
  isten(3)  = b(5)
  isten(5)  = 3

  call get_fluxes(isten,fz) 
  flag=0
  call adjust_for_boundary(isten,3,flag)
  if(flag==1)then
    fz(2,:)=zero
    fz(2,4)=p(igrid)
  endif
  flag=0
  call adjust_for_boundary(isten,1,flag)
  if(flag==1)then
    fz(1,:)=zero
    fz(1,4)=p(igrid)
  endif


  do idim=1,5
     fluxtmp(idim,igrid)=-(areayz*(fx(2,idim)-fx(1,idim))+&
                           areaxz*(fy(2,idim)-fy(1,idim))+&
                           areaxy*(fz(2,idim)-fz(1,idim)))*dt/vol
  enddo
 enddo
!$OMP ENDDO 
!$OMP BARRIER
!
!***
! Done with the big loop
! Now update the conserved quantities
!***
!
if(avg==0)then
!$OMP DO SCHEDULE(STATIC)
 do igrid=1,ngrid
  do idim=1,5
   cons(idim,igrid)=cons(idim,igrid)+fluxtmp(idim,igrid)
  enddo
 enddo
!$OMP ENDDO
else
!$OMP DO SCHEDULE(STATIC)
 do igrid=1,ngrid
  do idim=1,5
   cons(idim,igrid)=half*(cons_old(idim,igrid)+cons(idim,igrid)+fluxtmp(idim,igrid))
  enddo
 enddo
!$OMP ENDDO
endif
!
!***
! Put things back in order
!***
!
!$OMP DO SCHEDULE(STATIC)
 do igrid=1,ngrid
  cons(1,igrid)=max(cons(1,igrid),small_rho)
  if(cons(5,igrid)<zero)then
    cons(5,igrid)=half*(cons(2,igrid)**2+cons(3,igrid)**2+cons(4,igrid)**2)/cons(1,igrid)+small_eps
  endif
 enddo
!$OMP ENDDO
!
!***
! If fluxing angular momentum, put it back into linear momentum.
!***
!
 if(fluxangmom)then
!$OMP DO SCHEDULE(STATIC) private(x,y,r,angmom,rmom) !!!private(ekin)
 do igrid=1,ngrid
    x=grid(igrid)%x;y=grid(igrid)%y+yoffset
    r=sqrt(x*x+y*y)
    angmom=cons(3,igrid)
    rmom=cons(2,igrid)
    cons(3,igrid)=(rmom*y+x*angmom/r)/(x*x/r+y*y/r) 
    if(y==zero)then
      cons(2,igrid)=(rmom-y/r*cons(3,igrid))*r/x
    else
      cons(2,igrid)=(x*cons(3,igrid)-angmom)/y
    endif
 enddo
!$OMP ENDDO
 endif ! fluxangmom
 call set_ghost_cells() ! utils.f90

!$OMP END PARALLEL

end subroutine

subroutine clear_slope(id)
 use parameters
 use grid_commons
 implicit none
!
 integer, intent(in)::id
!
 slope_u(:,:,id)=zero
 slope_p(:,id)=zero
 slope_d(:,id)=zero
 slope_e(:,id)=zero
 slope_g(:,id)=zero
end subroutine
!
!
!
subroutine calculate_slopes(idx)
 use parameters
 use grid_commons
 use utils,only:slope
 implicit none

 integer,intent(in):: idx(5)
 integer:: idim,idir

 real(pre)::eps(3)

 idir=idx(5)
!
! slopes for velocity
!
 do idim=1,3
   slope_u(idir,idim,idx(2))=slope(u(idim,idx(1)),u(idim,idx(2)),u(idim,idx(3)))
 enddo
!
! slopes for pressure,density,energy,and gamma
!
  do idim=1,3
    eps(idim)=max(cons(5,idx(idim))-half*(cons(2,idx(idim))**2 &
                                   +cons(3,idx(idim))**2       &
                                   +cons(4,idx(idim))**2       &
                                   )/cons(1,idx(idim)),small_eps)
  enddo
!
  slope_p(idir,idx(2))=slope(p(idx(1)),p(idx(2)),p(idx(3)))
  slope_d(idir,idx(2))=slope(cons(1,idx(1)),cons(1,idx(2)),cons(1,idx(3)))
  slope_e(idir,idx(2))=slope(eps(1),eps(2),eps(3))
  slope_g(idir,idx(2))=slope(adindx(idx(1)),adindx(idx(2)),adindx(idx(3)))
 
end subroutine
!
!
!
subroutine interface_muscl(q,lm,lp,s,minus,plus)
 use parameters
 real(pre),intent(in)::q,lp,lm,s
 real(pre)::minus,plus
!
 minus=q+half*(-one-lm)*s
 plus =q+half*( one-lp)*s
!
!
end subroutine
!
!
!
subroutine calculate_states(isten)
 use parameters
 use grid_commons
 use eos, only :  get_gamma_from_p
 implicit none
!
 integer,intent(in)::isten(5)
 integer::idim,jdim,idir1,idir2,idx
!
 real(pre)::sound,lambda_p,lambda_m,ds(3),transverse
 real(pre)::lambda_p1, lambda_p2, lambda_m1,lambda_m2,m,ekin,eps
 real(pre)::cc,rho,uvel,vvel,wvel,drx,dry,drz,dux,dvx,dvz,duy,dvy,dwy,duz,dwz,dwx
 real(pre)::sr0,su0,sv0,sw0,sp0,dpx,dpy,dpz,dex,dey,dez,dgx,dgy,dgz,se0,sg0
!
 idx=isten(2)
 sound=sqrt(p(idx)*adindx(idx)/cons(1,idx))
 ds=[dx,dy,dz]
 do idim=1,3
  lambda_p=max(u(idim,idx)+sound,zero)*dt/ds(idim)
  lambda_m=min(u(idim,idx)-sound,zero)*dt/ds(idim)
  !lambda_p=u(idim,idx)*dt/ds(idim)
  !lambda_m=u(idim,idx)*dt/ds(idim)
!
  do jdim=1,3
     call interface_muscl(u(jdim,idx),lambda_m,lambda_p,slope_u(idim,jdim,idx),state_u_m(idim,jdim,idx),&
          state_u_p(idim,jdim,idx))
!  
  enddo

  eps = max(cons(5,idx)-half*(cons(2,idx)**2+cons(3,idx)**2+cons(4,idx)**2)/cons(1,idx),small_eps)

  call interface_muscl(p(idx),lambda_m,lambda_p,slope_p(idim,idx),state_p_m(idim,idx),&
       state_p_p(idim,idx))
  call interface_muscl(cons(1,idx),lambda_m,lambda_p,slope_d(idim,idx),state_d_m(idim,idx),&
      state_d_p(idim,idx))
  call interface_muscl(eps,lambda_m,lambda_p,slope_e(idim,idx),state_e_m(idim,idx),&
      state_e_p(idim,idx))
  call interface_muscl(adindx(idx),lambda_m,lambda_p,slope_g(idim,idx),state_g_m(idim,idx),&
      state_g_p(idim,idx))
!
 enddo
!
! now add corrections from transverse flux
!

 rho=cons(1,idx)
 uvel=u(1,idx)
 vvel=u(2,idx)
 wvel=u(3,idx)

 drx=slope_d(1,idx)
 dux=slope_u(1,1,idx)
 dvx=slope_u(1,2,idx)
 dwx=slope_u(1,3,idx)
 dpx=slope_p(1,idx)
 dex=slope_e(1,idx)
 dgx=slope_g(1,idx)

 dry=slope_d(2,idx)
 duy=slope_u(2,1,idx)
 dvy=slope_u(2,2,idx)
 dwy=slope_u(2,3,idx)
 dpy=slope_p(2,idx)
 dey=slope_e(2,idx)
 dgy=slope_g(2,idx)

 drz=slope_d(3,idx)
 duz=slope_u(3,1,idx)
 dvz=slope_u(3,2,idx)
 dwz=slope_u(3,3,idx)
 dpz=slope_p(3,idx)
 dez=slope_e(3,idx)
 dgz=slope_g(3,idx)

 su0=  zero    -vvel*duy -wvel*duz
 sv0=-uvel*dvx + zero    -wvel*dvz
 sw0=-uvel*dwx -vvel*dwy + zero 

 sr0=zero
 sp0=zero
 se0=zero
 sg0=zero

 cc=dt*half/dx
!
! Tried several variations of source terms.  
! Standard approach does not seem to be appropriate
! because not assuming a perfect gas (p=(g-1)eps).
! Blast waves look good at the moment, so will keep
! only the velocity cross terms for now until I
! figure out something better to do.  
!
! This can also be done with loops, etc., but I wrote it
! out so the cross terms become more apparent. 
!
! sr0=  zero   -vvel*dry -wvel*drz
! sp0=  zero   -vvel*dpy -wvel*dpz
! se0=  zero   -vvel*dey -wvel*dez
! sg0=  zero   -vvel*dgy -wvel*dgz
! x dir
 state_d_p(1,idx)  =state_d_p(1,idx)+sr0*cc
 state_u_p(1,1,idx)=state_u_p(1,1,idx)+su0*cc
 state_u_p(1,2,idx)=state_u_p(1,2,idx)+sv0*cc
 state_u_p(1,3,idx)=state_u_p(1,3,idx)+sw0*cc
 state_p_p(1,idx)  =state_p_p(1,idx)+sp0*cc
 state_g_p(1,idx)  =state_g_p(1,idx)+sg0*cc
 state_e_p(1,idx)  =state_e_p(1,idx)+se0*cc
!
 state_d_m(1,idx)  =state_d_m(1,idx)+sr0*cc
 state_u_m(1,1,idx)=state_u_m(1,1,idx)+su0*cc
 state_u_m(1,2,idx)=state_u_m(1,2,idx)+sv0*cc
 state_u_m(1,3,idx)=state_u_m(1,3,idx)+sw0*cc
 state_p_m(1,idx)  =state_p_m(1,idx)+sp0*cc
 state_e_m(1,idx)  =state_e_m(1,idx)+se0*cc
 state_g_m(1,idx)  =state_g_m(1,idx)+sg0*cc

 cc=dt*half/dy
! sr0=  -uvel*drx   +zero -wvel*drz 
! sp0=  -uvel*dpx   +zero -wvel*dpz 
! se0=  -uvel*dex   +zero -wvel*dez 
! sg0=  -uvel*dgx   +zero -wvel*dgz
!
! y dir
 state_d_p(2,idx)  =state_d_p(2,idx)+sr0*cc
 state_u_p(2,1,idx)=state_u_p(2,1,idx)+su0*cc
 state_u_p(2,2,idx)=state_u_p(2,2,idx)+sv0*cc
 state_u_p(2,3,idx)=state_u_p(2,3,idx)+sw0*cc
 state_p_p(2,idx)  =state_p_p(2,idx)+sp0*cc
 state_e_p(2,idx)  =state_e_p(2,idx)+se0*cc
 state_g_p(2,idx)  =state_g_p(2,idx)+sg0*cc
!
 state_d_m(2,idx)  =state_d_m(2,idx)+sr0*cc
 state_u_m(2,1,idx)=state_u_m(2,1,idx)+su0*cc
 state_u_m(2,2,idx)=state_u_m(2,2,idx)+sv0*cc
 state_u_m(2,3,idx)=state_u_m(2,3,idx)+sw0*cc
 state_p_m(2,idx)  =state_p_m(2,idx)+sp0*cc
 state_g_m(2,idx)  =state_g_m(2,idx)+sg0*cc
 state_e_m(2,idx)  =state_e_m(2,idx)+se0*cc

 cc=dt*half/dz

! sr0=  -uvel*drx   -vvel*dry  +zero  
! sp0=  -uvel*dpx   -vvel*dpy  +zero   
! se0=  -uvel*dex   -vvel*dey  +zero  
! sg0=  -uvel*dgx   -vvel*dgy  +zero
!
! z dir
 state_d_p(3,idx)  =state_d_p(3,idx)+sr0*cc
 state_u_p(3,1,idx)=state_u_p(3,1,idx)+su0*cc
 state_u_p(3,2,idx)=state_u_p(3,2,idx)+sv0*cc
 state_u_p(3,3,idx)=state_u_p(3,3,idx)+sw0*cc
 state_p_p(3,idx)  =state_p_p(3,idx)+sp0*cc
 state_g_p(3,idx)  =state_g_p(3,idx)+sg0*cc
 state_e_p(3,idx)  =state_e_p(3,idx)+se0*cc
!
 state_d_m(3,idx)  =state_d_m(3,idx)+sr0*cc
 state_u_m(3,1,idx)=state_u_m(3,1,idx)+su0*cc
 state_u_m(3,2,idx)=state_u_m(3,2,idx)+sv0*cc
 state_u_m(3,3,idx)=state_u_m(3,3,idx)+sw0*cc
 state_p_m(3,idx)  =state_p_m(3,idx)+sp0*cc
 state_g_m(3,idx)  =state_g_m(3,idx)+sg0*cc
 state_e_m(3,idx)  =state_e_m(3,idx)+se0*cc


end subroutine
!
!
!

subroutine get_fluxes(isten,fhll)
 use parameters
 use grid_commons
 
 implicit none

  integer::isten(5),idir,m_or_p,idx,imp
  real(pre)::fluxes_l(5),fluxes_r(5),states_l(5),states_r(5),fhll(2,5)
  real(pre)::vxl,vxr,vyl,vyr,vzl,vzr,gl,gr,pl,pr,dl,dr,el,er,vr,vl,fac

  idir=isten(5)
  idx =isten(2)

! complete minus state

 do m_or_p=0,1

  dl=state_d_p(idir,isten(1+m_or_p)) ! plus state of the minus cell
  pl=state_p_p(idir,isten(1+m_or_p))
  el=state_e_p(idir,isten(1+m_or_p))
  gl=state_g_p(idir,isten(1+m_or_p))
  vxl=state_u_p(idir,1,isten(1+m_or_p))
  vyl=state_u_p(idir,2,isten(1+m_or_p))
  vzl=state_u_p(idir,3,isten(1+m_or_p))
   
  states_l(1)=dl
  states_l(2)=dl*vxl
  states_l(3)=dl*vyl
  states_l(4)=dl*vzl
  states_l(5)=el + half*dl*(vxl**2+vyl**2+vzl**2)

  dr=state_d_m(idir,isten(2+m_or_p)) ! minus state of active cell
  pr=state_p_m(idir,isten(2+m_or_p))
  er=state_e_m(idir,isten(2+m_or_p))
  gr=state_g_m(idir,isten(2+m_or_p))
  vxr=state_u_m(idir,1,isten(2+m_or_p))
  vyr=state_u_m(idir,2,isten(2+m_or_p))
  vzr=state_u_m(idir,3,isten(2+m_or_p))
   
  states_r(1)=dr
  states_r(2)=dr*vxr
  states_r(3)=dr*vyr
  states_r(4)=dr*vzr
  states_r(5)=er + half*dr*(vxr**2+vyr**2+vzr**2)

  fluxes_l=zero
  fluxes_r=zero
  if(isten(5)==1)then
    vr=vxr
    vl=vxl
    fluxes_l(2)=pl
    fluxes_r(2)=pr
  elseif(isten(5)==2)then
    vr=vyr
    vl=vyl
    fluxes_l(3)=pl
    fluxes_r(3)=pr
  else
    vr=vzr
    vl=vzl
    fluxes_l(4)=pl
    fluxes_r(4)=pr
  endif

!
! left states
!
  fluxes_l(1)=states_l(1)*vl
  fluxes_l(2)=fluxes_l(2)+states_l(2)*vl
  fluxes_l(3)=fluxes_l(3)+states_l(3)*vl
  fluxes_l(4)=fluxes_l(4)+states_l(4)*vl
  fluxes_l(5)=(states_l(5)+pl)*vl
!
! right states  
!
  fluxes_r(1)=states_r(1)*vr
  fluxes_r(2)=fluxes_r(2)+states_r(2)*vr
  fluxes_r(3)=fluxes_r(3)+states_r(3)*vr
  fluxes_r(4)=fluxes_r(4)+states_r(4)*vr
  fluxes_r(5)=(states_r(5)+pr)*vr

  imp=m_or_p+1
  call flux_hll(fluxes_l,fluxes_r,states_l,states_r,gl,gr,dl,dr,pl,pr,vl,vr,fhll(imp,:))

 enddo ! m_or_p

end subroutine

subroutine flux_hll(fl,fr,ul,ur,gl,gr,dl,dr,pl,pr,vl,vr,fhll)
    use parameters
    implicit none
    integer::i
    real(pre)::fl(5),fr(5),fhll(5)
    real(pre)::ul(5),ur(5)
    real(pre)::gl,pl,dl,gr,pr,dr,cr,cl,ap,am,vr,vl
    !print *,"In flux_hll",gr,gl,pl,pr,dr,dl
    cr=sqrt(gr*pr/dr)
    cl=sqrt(gl*pl/dl)
    ap=max(zero,vr+cr,vl+cl)
    am=max(zero,-(vr-cr),-(vl-cl))
    do i=1,5
       fhll(i)=(ap*fl(i)+am*fr(i)-ap*am*(ur(i)-ul(i)))/(ap+am)
    enddo
end subroutine

subroutine adjust_for_boundary(isten,ib,flag)
    use parameters
    use grid_commons
    implicit none
    integer::flag,ib,isten(5),i
    flag=0
    if(grid(isten(ib))%boundary>2)then
      flag=1
    endif
end subroutine 

