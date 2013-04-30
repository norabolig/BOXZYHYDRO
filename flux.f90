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
 real(pre)::fxl(5),fxr(5),fyt(5),fyb(5),fzt(5),fzb(5),areaxy,areayz,areaxz
 real(pre)::vol,left,right
 real(pre)::x,y,z,r,angmom,rmom

 type(units)::scale

 call get_units(scale)

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP&PRIVATE(areaxy,areayz,areaxz) &
!$OMP&private(fzt,fzb,fxr,fxl,fyt,fyb,left,right) &
!$OMP&private(b,angmom,rmom,x,y,z,r,isten,flag)

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
!$OMP DO SCHEDULE(static)
 do igrid=1,ngrid

  fluxtmp(:,igrid)=zero
  if(grid(igrid)%boundary>0)cycle

  call get_boundary_wb(igrid,b)

  isten(1)  = b(4)
  isten(2)  = igrid
  isten(3)  = b(3)
  isten(4)  = grid(b(3))%ineigh(3)
  isten(5)  = 1
  call adjust_for_boundary(isten,3,2,flag)
  call get_states(isten,fxr)
  if(flag==1)then
     fxr=zero
     fxr(2)=p(igrid)
  endif

  isten(1)  = grid(b(4))%ineigh(4)
  isten(2)  = b(4)
  isten(3)  = igrid
  isten(4)  = b(3)
  isten(5)  = 1
  call adjust_for_boundary(isten,2,3,flag)
  call get_states(isten,fxl)
  if(flag==1)then
     fxl=zero
     fxl(2)=p(igrid)
  endif


  isten(1)  = b(2)
  isten(2)  = igrid
  isten(3)  = b(1)
  isten(4)  = grid(b(1))%ineigh(1)
  isten(5)  = 2
  call adjust_for_boundary(isten,3,2,flag)
  call get_states(isten,fyt)
  if(flag==1)then
     fyt=zero
     fyt(3)=p(igrid)
  endif

  isten(1)  = grid(b(2))%ineigh(2)
  isten(2)  = b(2)
  isten(3)  = igrid
  isten(4)  = b(1)
  isten(5)  = 2
  call adjust_for_boundary(isten,2,3,flag)
  call get_states(isten,fyb)
  if(flag==1)then
     fyb=zero
     fyb(3)=p(igrid)
  endif

  isten(1)  = b(6)
  isten(2)  = igrid
  isten(3)  = b(5)
  isten(4)  = grid(b(5))%ineigh(5)
  isten(5)  = 3
  call adjust_for_boundary(isten,3,2,flag)
  call get_states(isten,fzt)
  if(flag==1)then
     fzt=zero
     fzt(4)=p(igrid)
  endif

  isten(1)  = grid(b(6))%ineigh(6)
  isten(2)  = b(6)
  isten(3)  = igrid
  isten(4)  = b(5)
  isten(5)  = 3
  call adjust_for_boundary(isten,2,3,flag)
  call get_states(isten,fzb)
  if(flag==1)then
     fzb=zero
     fzb(4)=p(igrid)
  endif

  do idim=1,5
     fluxtmp(idim,igrid)=-( areaxy*(fzt(idim)-fzb(idim))+areayz*(fxr(idim)-fxl(idim))&
                        +areaxz*(fyt(idim)-fyb(idim)))/(vol)*dt
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

subroutine get_states(idx,fhll)
 use parameters
 use grid_commons
 use utils,only : left_right_states
 implicit none

  integer::idx(5)
  real(pre)::fluxes_l(5),fluxes_r(5),states_l(5),states_r(5),fhll(5)
  real(pre)::vxl,vxr,vyl,vyr,vzl,vzr,gl,gr,pl,pr,dl,dr,el,er,vr,vl,fac

  call left_right_states(u(1,idx(1)),u(1,idx(2)),u(1,idx(3)),u(1,idx(4)),vxl,vxr)
  call left_right_states(u(2,idx(1)),u(2,idx(2)),u(2,idx(3)),u(2,idx(4)),vyl,vyr)
  call left_right_states(u(3,idx(1)),u(3,idx(2)),u(3,idx(3)),u(3,idx(4)),vzl,vzr)
  call left_right_states(p(idx(1))  ,p(idx(2))  ,p(idx(3))  ,p(idx(4))  ,pl,pr)
  call left_right_states(adindx(idx(1)),adindx(idx(2)),adindx(idx(3)),adindx(idx(4)),gl,gr)
  call left_right_states(cons(1,idx(1)),cons(1,idx(2)),cons(1,idx(3)),cons(1,idx(4)),dl,dr)
  call left_right_states(cons(5,idx(1)),cons(5,idx(2)),cons(5,idx(3)),cons(5,idx(4)),el,er)

  if(idx(5)==1)then
    vr=vxr
    vl=vxl
! left states
    fluxes_l(1)=dl*vxl
    fluxes_l(2)=dl*vxl*vxl+pl
    fluxes_l(3)=dl*vyl*vxl
    fluxes_l(4)=dl*vzl*vxl
    fluxes_l(5)=vxl*(el+pl)
! right states  
    fluxes_r(1)=dr*vxr
    fluxes_r(2)=dr*vxr*vxr+pr
    fluxes_r(3)=dr*vyr*vxr
    fluxes_r(4)=dr*vzr*vxr
    fluxes_r(5)=vxr*(er+pr)
  elseif(idx(5)==2)then
    vr=vyr
    vl=vyl
! left states
    fluxes_l(1)=dl*vyl
    fluxes_l(2)=dl*vxl*vyl
    fluxes_l(3)=dl*vyl*vyl+pl
    fluxes_l(4)=dl*vzl*vyl
    fluxes_l(5)=vyl*(el+pl)
! right states  
    fluxes_r(1)=dr*vyr
    fluxes_r(2)=dr*vxr*vyr
    fluxes_r(3)=dr*vyr*vyr+pr
    fluxes_r(4)=dr*vzr*vyr
    fluxes_r(5)=vyr*(er+pr)
  else
    vr=vzr
    vl=vzl
 ! left states
    fluxes_l(1)=dl*vzl
    fluxes_l(2)=dl*vxl*vzl
    fluxes_l(3)=dl*vyl*vzl
    fluxes_l(4)=dl*vzl*vzl+pl
    fluxes_l(5)=vzl*(el+pl)
! right states  
    fluxes_r(1)=dr*vzr
    fluxes_r(2)=dr*vxr*vzr
    fluxes_r(3)=dr*vyr*vzr
    fluxes_r(4)=dr*vzr*vzr+pr
    fluxes_r(5)=vzr*(er+pr)
  endif
!
  states_l(1)=dl
  states_l(2)=dl*vxl
  states_l(3)=dl*vyl
  states_l(4)=dl*vzl
  states_l(5)=el
!
  states_r(1)=dr
  states_r(2)=dr*vxr
  states_r(3)=dr*vyr
  states_r(4)=dr*vzr
  states_r(5)=er

  call flux_hll(fluxes_l,fluxes_r,states_l,states_r,gl,gr,dl,dr,pl,pr,vl,vr,fhll)

end subroutine

subroutine flux_hll(fl,fr,ul,ur,gl,gr,dl,dr,pl,pr,vl,vr,fhll)
    use parameters
    implicit none
    integer::i
    real(pre)::fl(5),fr(5),fhll(5)
    real(pre)::ul(5),ur(5)
    real(pre)::gl,pl,dl,gr,pr,dr,cr,cl,ap,am,vr,vl
    cr=sqrt(gr*pr/dr)
    cl=sqrt(gl*pl/dl)
    ap=max(zero,vr+cr,vl+cl)
    am=max(zero,-(vr-cr),-(vl-cl))
    do i=1,5
       fhll(i)=(ap*fl(i)+am*fr(i)-ap*am*(ur(i)-ul(i)))/(ap+am)
    enddo
end subroutine

subroutine adjust_for_boundary(isten,ib,ig,flag)
    use parameters
    use grid_commons
    implicit none
    integer::flag,ib,ig,isten(5),i
    flag=0
    if(grid(isten(ib))%boundary>2)then
      flag=1
    endif
end subroutine 

