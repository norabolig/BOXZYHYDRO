!
! Fluxing routine that uses centeral differencing or upwind schemes
! At the moment, variables f2,f1,f0,fm1 are calculated but not used
! This will likely be removed, but it may stick around for a bit
!
subroutine flux()
 use parameters
 use derived_types
 use grid_commons
 use utils
 implicit none

 real(pre),allocatable,dimension(:,:)::fluxtmp
 real(pre),allocatable,dimension(:)::rhophi,rhophi_pt

 integer :: igrid,b(6),idim,flag
 real(pre)::fxl,fxr,fyt,fyb,fzt,fzb,areaxy,areayz,areaxz
 real(pre)::vol,v,f2,f1,f0,fm1,slopef,slopeb,left,right
 real(pre)::x,y,z,r,angmom,rmom,cfast,alocal(6),clocal
 real(pre)::u1,u2,u0,um1,fleft,fright,theta=1.0d0

 logical :: assoc_con
 type(units)::scale

 call get_units(scale)

 allocate(fluxtmp(5,ngrid))
! allocate(rhophi(ngrid))
! allocate(rhophi_pt(ngrid))

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP&PRIVATE(areaxy,areayz,areaxz,u1,u2,u0,um1,fleft,fright) &
!$OMP&private(v,f2,f1,f0,fm1,slopef,slopeb,flag,cfast,clocal) &
!$OMP&private(fzt,fzb,fxr,fxl,fyt,fyb,left,right,alocal) &
!$OMP&private(b,angmom,rmom,x,y,z,r,assoc_con)

 areaxy=dy*dx
 areayz=dy*dz
 areaxz=dx*dz
 vol=dx*dy*dz
 assoc_con=associated(cons_pt,target=cons) ! if target is cons, don't double operate

!
!***
! Will be solving for the total energy, excluding gravity. 
! Add pressure to this quantity.
!***
!
!!$OMP DO SCHEDULE(STATIC)
!  do igrid=1,ngrid
!   rhophi(igrid)=cons(1,igrid)*phi(igrid)
!   rhophi_pt(igrid)=cons_pt(1,igrid)*phi(igrid)
!  enddo  
!!$OMP ENDDO
!$OMP DO SCHEDULE(STATIC)
 do igrid=1,ngrid
  cons(5,igrid)=cons(5,igrid)+p(igrid)!+rhophi(igrid)
 enddo
!$OMP ENDDO
 if(.not.assoc_con)then
!$OMP DO SCHEDULE(STATIC)
 do igrid=1,ngrid
  cons_pt(5,igrid)=cons_pt(5,igrid)+p(igrid)!+rhophi_pt(igrid)
 enddo
!$OMP ENDDO
 endif
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
  if(.not.assoc_con)then
!$OMP DO SCHEDULE(STATIC) private(x,y,r,angmom,rmom)
   do igrid=1,ngrid
     x=grid(igrid)%x;y=grid(igrid)%y+yoffset
     r=sqrt(x*x+y*y)
     angmom=x*cons_pt(3,igrid)-y*cons_pt(2,igrid)
     rmom=cons_pt(2,igrid)*x/r+cons_pt(3,igrid)*y/r
     cons_pt(2,igrid)=rmom
     cons_pt(3,igrid)=angmom
   enddo
!$OMP ENDDO
  endif
 endif
!
!***
! Enter main loop
!***
!
!$OMP DO SCHEDULE(static)
 do igrid=1,ngrid

  if(grid(igrid)%boundary>0)cycle

  call get_boundary_wb(igrid,b)
  x=grid(igrid)%x;y=grid(igrid)%y+yoffset;z=grid(igrid)%z
  r=sqrt(x*x+y*y)

  clocal=sqrt(max(adindx(igrid)*p(igrid)/cons(1,igrid),zero))
  do idim=1,6
     alocal(idim)=sqrt(max(adindx(b(idim))*p(b(idim))/cons(1,b(idim)),zero))
     alocal(idim)=half*(alocal(idim)+clocal)
  enddo

  do idim=1,5

      v=half*(u(2,igrid)+u(2,b(1)))

      u1=cons(idim,b(1))
      u0=cons(idim,igrid)
      um1=cons(idim,b(2))
      f1=u1*u(2,b(1))
      f0=u0*u(2,igrid)
      fm1=um1*u(2,b(2))

     flag=0
     if(grid(b(1))%boundary>0)then
       if(no_out_flow_y.or.(grid(b(1))%boundary>2))flag=1
     else
         u2=cons(idim,grid(b(1))%ineigh(1))
         f2=u2*u(2,grid(b(1))%ineigh(1))
     endif
 
     if(flag==0)then
!
!
#ifdef UPWIND
     fyt=v*upwind(u2,u1,u0,um1,v,dy,dt)
#else
!
!
#ifdef TCDIFFERENCE
!
!
     cfast=abs(v)*theta
     call left_right_states(u2,u1,u0,um1,left,right)
     fyt=v*half*(left+right)-cfast*(right-left)
!
!
#else
     stop"Unknown flux type"
#endif /* end ifdef TCDIFFERENCE */
!
!
#endif /* end ifdef UPWIND */
!
!
     else
       fyt=zero
     endif

     v=-half*(u(2,igrid)+u(2,b(2)))

     u1=cons(idim,b(2))
     u0=cons(idim,igrid)
     um1=cons(idim,b(1))
     f1=-u1*u(2,b(2))
     f0=-u0*u(2,igrid)
     fm1=-um1*u(2,b(1))
     flag=0

     if(grid(b(2))%boundary>0)then
       if(no_out_flow_y.or.(grid(b(2))%boundary>2))flag=1
     else
       u2=cons(idim,grid(b(2))%ineigh(2))
       f2=-u2*u(2,grid(b(2))%ineigh(2))
     endif
      
     if(flag==0)then 
!
!
#ifdef UPWIND
     fyb=v*upwind(u2,u1,u0,um1,v,dy,dt)
#else
!
!
#ifdef TCDIFFERENCE
!
!
     cfast=abs(v)*theta
     call left_right_states(u2,u1,u0,um1,left,right)
     fyb=v*half*(left+right)-cfast*(right-left)
!
!
#endif /* end ifdef TCDIFFERENCE */
!
!
#endif /* end ifdef UPWIND */
!
!
     else
      fyb=zero
     endif
!
!***
!X 3 and 4
!***
!
      v=half*(u(1,igrid)+u(1,b(3)))

      u1=cons(idim,b(3))
      u0=cons(idim,igrid)
      um1=cons(idim,b(4))
      f1=u1*u(1,b(3))
      f0=u0*u(1,igrid)
      fm1=um1*u(1,b(4))
 
     flag=0

     if(grid(b(3))%boundary>0)then
       if(no_out_flow_x.or.(grid(b(3))%boundary>2))flag=1
     else
         u2=cons(idim,grid(b(3))%ineigh(3))
         f2=u2*u(1,grid(b(3))%ineigh(3))
     endif
 
     if(flag==0)then
!
!
#ifdef UPWIND
     fxr=v*upwind(u2,u1,u0,um1,v,dx,dt)
#else
!
!
#ifdef TCDIFFERENCE
!
!
     cfast=abs(v)*theta
     call left_right_states(u2,u1,u0,um1,left,right)
     fxr=v*half*(left+right)-cfast*(right-left)
!
!
#endif /* end ifdef TCDIFFERENCE */
!
!
#endif /* end ifdef UPWIND */
!
!
     else
      fxr=zero
     endif

     v=-half*(u(1,igrid)+u(1,b(4)))

     u1=cons(idim,b(4))
     u0=cons(idim,igrid)
     um1=cons(idim,b(3))
     f1=-u1*u(1,b(4))
     f0=-u0*u(1,igrid)
     fm1=-um1*u(1,b(3))

     flag=0

     if(grid(b(4))%boundary>0)then
       if(no_out_flow_x.or.(grid(b(4))%boundary>2))flag=1
     else
       u2=cons(idim,grid(b(4))%ineigh(4))
       f2=-u2*u(1,grid(b(4))%ineigh(4))
     endif
       

     if(flag==0)then
!
!
#ifdef UPWIND
     fxl=v*upwind(u2,u1,u0,um1,v,dx,dt)
#else
!
!
#ifdef TCDIFFERENCE
!
!
     cfast=abs(v)*theta
     call left_right_states(u2,u1,u0,um1,left,right)
     fxl=v*half*(left+right)-cfast*(right-left)
!
!
#endif /* end ifdef TCDIFFERENCE */
!
!
#endif /* end ifdef UPWIND */
!
!
     else
      fxl=zero
     endif
!
!***
!Z 5 and 6
!***
!
      v=half*(u(3,igrid)+u(3,b(5)))

      u1=cons(idim,b(5))
      u0=cons(idim,igrid)
      um1=cons(idim,b(6))
      f1=u1*u(3,b(5))
      f0=u0*u(3,igrid)
      fm1=um1*u(3,b(6))
 
     flag=0

     if(grid(b(5))%boundary>0)then
       if(no_out_flow_z.or.(grid(b(5))%boundary>2))flag=1
     else
         u2=cons(idim,grid(b(5))%ineigh(5))
         f2=u2*u(3,grid(b(5))%ineigh(5))
     endif

     if(flag==0)then
!
!
#ifdef UPWIND
     fzt=v*upwind(u2,u1,u0,um1,v,dz,dt)
#else
!
!
#ifdef TCDIFFERENCE
!
!
     cfast=abs(v)*theta
     call left_right_states(u2,u1,u0,um1,left,right)
     fzt=v*half*(left+right)-cfast*(right-left)
!
!
#endif /* end ifdef TCDIFFERENCE */
!
!
#endif /* end ifdef UPWIND */
!
!
     else
       fzt=zero
     endif

     v=-half*(u(3,igrid)+u(3,b(6)))

     u1=cons(idim,b(6))
     u0=cons(idim,igrid)
     um1=cons(idim,b(5))
     f1=-u1*u(3,b(6))
     f0=-u0*u(3,igrid)
     fm1=-um1*u(3,b(5))

     flag=0
     if(grid(b(6))%boundary>0)then
       if(no_out_flow_z.or.(grid(b(6))%boundary>2))flag=1
     else
         u2=cons(idim,grid(b(6))%ineigh(6))
         f2=-u2*u(3,grid(b(6))%ineigh(6))
     endif
       
     if(flag==0)then
!
!
#ifdef UPWIND
     fzb=v*upwind(u2,u1,u0,um1,v,dz,dt)
#else
!
!
#ifdef TCDIFFERENCE
!
!
     cfast=abs(v)*theta
     call left_right_states(u2,u1,u0,um1,left,right)
     fzb=v*half*(left+right)-cfast*(right-left)
!
!
#endif /* end ifdef TCDIFFERENCE */
!
!
#endif /* end ifdef UPWIND */
!
!
     else
       fzb=zero
     endif

     fluxtmp(idim,igrid)=-( areaxy*(fzt+fzb)+areayz*(fxr+fxl)+areaxz*(fyt+fyb))&
                     /(vol)*dt
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
!$OMP DO SCHEDULE(STATIC)
 do igrid=1,ngrid
  do idim=1,5
   cons_pt(idim,igrid)=cons_pt(idim,igrid)+fluxtmp(idim,igrid)
  enddo
 enddo
!$OMP ENDDO
!
!***
! Put things back in order
!***
!
!$OMP DO SCHEDULE(STATIC)
 do igrid=1,ngrid
  cons(5,igrid)=cons(5,igrid)-p(igrid)!-rhophi(igrid)
  if(cons(5,igrid)<zero)then
    cons(5,igrid)=half*(cons(2,igrid)**2+cons(3,igrid)**2+cons(4,igrid)**2)/cons(1,igrid)
  endif
 enddo
!$OMP ENDDO
 if(.not.assoc_con)then
!$OMP DO SCHEDULE(STATIC)
 do igrid=1,ngrid
  cons_pt(1,igrid)=max(cons_pt(1,igrid),small_rho)
  cons_pt(5,igrid)=cons_pt(5,igrid)-p(igrid)!-rhophi_pt(igrid)
  if(cons_pt(5,igrid)<zero)then
    cons_pt(5,igrid)=half*(cons_pt(2,igrid)**2+cons_pt(3,igrid)**2+cons_pt(4,igrid)**2) &
                    /cons_pt(1,igrid)
  endif
 enddo
!$OMP ENDDO
 endif
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
 if(.not.assoc_con)then
!$OMP DO SCHEDULE(STATIC) private(x,y,r,angmom,rmom) !!!private(ekin)
 do igrid=1,ngrid
    x=grid(igrid)%x;y=grid(igrid)%y+yoffset
    r=sqrt(x*x+y*y)
    cons_pt(1,igrid)=max(cons_pt(1,igrid),small_rho)
     angmom=cons_pt(3,igrid)
     rmom=cons_pt(2,igrid)
     cons_pt(3,igrid)=(rmom*y+x*angmom/r)/(x*x/r+y*y/r) 
     if(y==zero)then
       cons_pt(2,igrid)=(rmom-y/r*cons_pt(3,igrid))*r/x
     else
       cons_pt(2,igrid)=(x*cons_pt(3,igrid)-angmom)/y
     endif
 enddo
!$OMP ENDDO
 endif ! assoc_con
 endif ! fluxangmom
 call set_ghost_cells_pt() ! utils.f90

!$OMP END PARALLEL

 deallocate(fluxtmp)!,rhophi,rhophi_pt)

end subroutine
