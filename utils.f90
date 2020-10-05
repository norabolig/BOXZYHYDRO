
!
! A collection of a variety of functions and routines for the code
!
module utils
 use parameters
 use derived_types
 implicit  none
!
!***
! These variables need to be shared for reduction methods in parallelized loops
!***
!
 real(pre),save::mass,eint,ekin,etot,egra,amom,momx,momy,momz,xcom,ycom,zcom
 real(pre),save::xcom0,ycom0,zcom0
 logical,save::comfirst=.true.

 contains
!
!
real(pre) function set_inverse(xx)
  real(pre)::xx
  if(xx/=zero)then
    set_inverse=one/xx
  else
    set_inverse=zero
  endif
end function
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Find the grid index for a given x,y,z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 integer function get_grid_indx(x,y,z)
   use grid_commons
   real(pre)::x,y,z
   integer::ix,igrid,iz,iy

   iz=int((z+half*dz)/dz+dble(nz+1)*half)
   if(iz<1)iz=1
   if(iz>nz)iz=nz

   igrid=(iz-1)*nx*ny
   
   iy=int((y+half*dy)/dy+dble(ny+1)*half)
   if(iy<1)iy=1
   if(iy>ny)iy=ny

   ix=int((x+half*dx)/dx+dble(nx+1)*half)
   if(ix<1)ix=1
   if(ix>nx)ix=nx 

   igrid=igrid+nx*(iy-1)+ix
   get_grid_indx=igrid

 end function
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get cell neighbors, excluding boundary cells.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 subroutine get_boundary(igrid,b)
  use derived_types
  use grid_commons
  integer, intent(in)::igrid
  integer, intent(out)::b(6)
  integer::idx

   do idx=1,6
    b(idx)=grid(igrid)%ineigh(idx)
    if(grid(b(idx))%boundary>0)b(idx)=igrid
   enddo

 end subroutine get_boundary
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get cell neighbors, excluding anchor cells.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 subroutine get_boundary_noanchor(igrid,b)
  use derived_types
  use grid_commons
  integer, intent(in)::igrid
  integer, intent(out)::b(6)
  integer::idx

   do idx=1,6
    b(idx)=grid(igrid)%ineigh(idx)
    if(grid(b(idx))%boundary>2)b(idx)=igrid
   enddo

 end subroutine get_boundary_noanchor
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get cell neighbors, including boundary cells.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 subroutine get_boundary_wb(igrid,b)
  use derived_types
  use grid_commons
  integer, intent(in)::igrid
  integer, intent(out)::b(6)
  integer::idx
  logical::active_grid

  do idx=1,6
   b(idx)=grid(igrid)%ineigh(idx)
  enddo

 end subroutine get_boundary_wb
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get nearest 8 cells for a random x,y,z position
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 subroutine get_nearest_8(x,y,z,neigh)
  use grid_commons
  integer::igrid,igridt
  integer,intent(out)::neigh(8)
  integer::bb(6)
  integer::bt(6)
  real(pre),intent(in)::x,y,z
  real(pre)::delz

  igrid=get_grid_indx(x,y,z)
  if(grid(igrid)%boundary>0)then
      neigh=-1
      return
  endif 
  call get_boundary_wb(igrid,bb)

  delz=z-grid(igrid)%z
  if(delz<zero)then
     igrid=bb(6)
     call get_boundary_wb(igrid,bb)
  endif
  igridt=bb(5)
  call get_boundary_wb(igridt,bt)
  
  neigh(5)=bb(5)
  neigh(6)=igrid

  if(x-grid(igrid)%x<zero )then
    if(y-grid(igrid)%y<zero) then
      neigh(1)=grid(bb(2))%ineigh(4)
      neigh(2)=bb(2)
      neigh(3)=igrid
      neigh(4)=grid(igrid)%ineigh(4)
      neigh(5)=grid(bt(2))%ineigh(4)
      neigh(6)=bt(2)
      neigh(7)=igridt
      neigh(8)=grid(igridt)%ineigh(4)
    else
      neigh(1)=grid(igrid)%ineigh(4)
      neigh(2)=igrid
      neigh(3)=bb(1)
      neigh(4)=grid(bb(1))%ineigh(4)
      neigh(5)=grid(igridt)%ineigh(4)
      neigh(6)=igridt
      neigh(7)=bt(1)
      neigh(8)=grid(bt(1))%ineigh(4)
    endif
  else
    if(y-grid(igrid)%y<zero) then
      neigh(1)=bb(2)
      neigh(2)=grid(bb(2))%ineigh(3)
      neigh(3)=grid(igrid)%ineigh(3)
      neigh(4)=igrid
      neigh(5)=bt(2)
      neigh(6)=grid(bt(2))%ineigh(3)
      neigh(7)=grid(igridt)%ineigh(3)
      neigh(8)=igridt
    else
      neigh(1)=igrid
      neigh(2)=bb(3)
      neigh(3)=grid(bb(3))%ineigh(1)
      neigh(4)=bb(1)
      neigh(5)=igridt
      neigh(6)=bt(3)
      neigh(7)=grid(bt(3))%ineigh(1)
      neigh(8)=bt(1)
    endif
   endif
  
 end subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get grid ID for a cell with indexes ix,iy,iz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 integer function get_id(ix,iy,iz)
  integer,intent(in)::ix,iy,iz
  get_id=ix+(iy-1)*nx+(iz-1)*nx*ny
 end function
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Check whether a boundary cell is near a corner of the domain.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 integer function check_corner(ix,iy,iz,ineigh)
  integer,intent(in)::ix,iy,iz
  integer::ineigh

  if (ix<3) then
    if (iy<3) then
      if(nz>=6)then
        if (iz<3)then; ineigh=get_id(3,3,3)
        elseif(iz>nz-2)then; ineigh=get_id(3,3,nz-2)
        else; ineigh=get_id(3,3,iz)
        endif
      else
        ineigh=get_id(3,3,iz)
      endif
    else if (iy>ny-2) then
      if(iz>=6)then
        if (iz<3)then; ineigh=get_id(3,ny-2,3)
        elseif(iz>nz-2)then; ineigh=get_id(3,ny-2,nz-2)
        else; ineigh=get_id(3,ny-2,iz)
        endif
      else
        ineigh=get_id(3,ny-2,iz)
      endif
    else 
      if(nz>=6)then
        if (iz<3)then; ineigh=get_id(3,iy,3)
        elseif(iz>nz-2)then; ineigh=get_id(3,iy,nz-2)
        endif
      else
        ineigh=get_id(3,iy,iz)
      endif
    endif
  else if (ix>nx-2) then
    if (iy<3) then
      if(nz>=6)then
        if (iz<3)then; ineigh=get_id(nx-2,3,3)
        elseif(iz>nz-2)then; ineigh=get_id(nx-2,3,nz-2)
        else; ineigh=get_id(nx-2,3,iz)
        endif
      else
        ineigh=get_id(nx-2,3,iz)
      endif
    else if (iy>ny-2) then
      if(nz>=6)then
        if (iz<3)then; ineigh=get_id(nx-2,ny-2,3)
        elseif(iz>nz-2)then; ineigh=get_id(nx-2,ny-2,nz-2)
        else; ineigh=get_id(nx-2,ny-2,iz)
        endif
      else
        ineigh=get_id(nx-2,ny-2,iz)
      endif
    else 
      if(nz>=6)then
        if (iz<3)then; ineigh=get_id(nx-2,iy,3)
        elseif(iz>nz-2)then; ineigh=get_id(nx-2,iy,nz-2)
        endif
      else
        ineigh=get_id(nx-2,iy,iz)
      endif
    endif
  else 
    if (iy<3) then
      if(nz>=6)then
        if (iz<3)then; ineigh=get_id(ix,3,3)
        elseif(iz>nz-2)then; ineigh=get_id(ix,3,nz-2)
        endif
      else
        ineigh=get_id(ix,3,iz)
      endif
    else if (iy>ny-2) then
      if(nz>=6)then
        if (iz<3)then; ineigh=get_id(ix,ny-2,3)
        elseif(iz>nz-2)then; ineigh=get_id(ix,ny-2,nz-2)
        endif
      else
        ineigh=get_id(ix,ny-2,iz)
      endif
    endif
  endif

  check_corner=ineigh
 end function
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set ghost cells 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 subroutine set_ghost_helper(idx,ip,fac)
  use grid_commons
  integer::idx,ip
  real(pre)::fac(3)
    cons(1,idx)= cons(1,ip)
    cons(2,idx)= cons(2,ip)*fac(1)
    cons(3,idx)= cons(3,ip)*fac(2)
    cons(4,idx)= cons(4,ip)*fac(3)
    cons(5,idx)= cons(5,ip)
 end subroutine
!
 subroutine zero_flow_dir(idx,idir)
  use grid_commons
  integer::idx,idir
    cons(5,idx)=max(cons(5,idx)-half*cons(idir,idx)**2/cons(1,idx),small_eps)
    cons(idir+1,idx)=zero
 end subroutine
!
 subroutine set_ghost_cells()
  use grid_commons
  use eos, only : get_gamma_from_tk
  integer::igrid,idx,idim,ineigh,ix,iy,iz,ip,ib
  real(pre)::v,eps,ploc,mucloc,adloc,fac(3)
  type(units)::scale

  call get_units(scale) ! units.f90

!$OMP DO SCHEDULE(STATIC) PRIVATE(idx,ineigh,eps,mucloc,adloc,ix,iy,iz,fac,ip,ib)
  do igrid=1,nghost
    idx=indx_ghost(igrid)
    ix=grid(idx)%ix
    iy=grid(idx)%iy
    iz=grid(idx)%iz
    ib=grid(idx)%boundary
    if (ib==1)then
       ineigh=grid(idx)%ineigh(1)
    else
       if    (ix==2   )then; ineigh=grid(idx)%ineigh(3)
       elseif(ix==nx-1)then; ineigh=grid(idx)%ineigh(4)
       elseif(iy==2   )then; ineigh=grid(idx)%ineigh(1)
       elseif(iy==ny-1)then; ineigh=grid(idx)%ineigh(2)
       elseif(iz==2   )then; ineigh=grid(idx)%ineigh(5)
       elseif(iz==nz-1)then; ineigh=grid(idx)%ineigh(6)    
       endif
    endif
    ineigh=check_corner(ix,iy,iz,ineigh)
 
       do idim=1,5
         cons(idim,idx)=cons(idim,ineigh)
       enddo

       if(no_outflow_xl)then
         fac(1:3)=one
         fac(1)=-one 
         if(ix==1)then
           ip=nx*ny*(iz-1)+nx*(iy-1)+ix+3
           call set_ghost_helper(idx,ip,fac)
         elseif(ix==2)then
           ip=nx*ny*(iz-1)+nx*(iy-1)+ix+1
           call set_ghost_helper(idx,ip,fac)
         endif
       endif
       if(no_outflow_xr)then
         fac(1:3)=one
         fac(1)=-one 
         if(ix==nx)then
           ip=nx*ny*(iz-1)+nx*(iy-1)+ix-3
           call set_ghost_helper(idx,ip,fac)
         elseif(ix==nx-1)then
           ip=nx*ny*(iz-1)+nx*(iy-1)+ix-1
           call set_ghost_helper(idx,ip,fac)
         endif
       endif
       if(no_outflow_yl)then
         fac(1:3)=one
         fac(2)=-one 
         if(iy==1)then
           ip=nx*ny*(iz-1)+nx*((iy+3)-1)+ix
           call set_ghost_helper(idx,ip,fac)
         elseif(iy==2)then
           ip=nx*ny*(iz-1)+nx*((iy+1)-1)+ix
           call set_ghost_helper(idx,ip,fac)
         endif
       endif
       if(no_outflow_yr)then
         fac(1:3)=one
         fac(2)=-one 
         if(iy==ny)then
           ip=nx*ny*(iz-1)+nx*((iy-3)-1)+ix
           call set_ghost_helper(idx,ip,fac)
         elseif(iy==ny-1)then
           ip=nx*ny*(iz-1)+nx*((iy-1)-1)+ix
           call set_ghost_helper(idx,ip,fac)
         endif
       endif
       if(no_outflow_zl)then
         fac(1:3)=one
         fac(3)=-one 
         if(iz==1)then
           ip=nx*ny*(iz+3-1)+nx*(iy-1)+ix
           call set_ghost_helper(idx,ip,fac)
         elseif(iz==2)then
           ip=nx*ny*(iz+1-1)+nx*(iy-1)+ix
           call set_ghost_helper(idx,ip,fac)
         endif
       endif
       if(no_outflow_zr)then
         fac(1:3)=one
         fac(3)=-one 
         if(iz==nz)then
           ip=nx*ny*(iz-3-1)+nx*(iy-1)+ix
           call set_ghost_helper(idx,ip,fac)
         elseif(iz==nz-1)then
           ip=nx*ny*(iz-1-1)+nx*(iy-1)+ix
           call set_ghost_helper(idx,ip,fac)
         endif
       endif
!! START for control valve
#ifdef CONTROLVALVES
       if(ix==1.or.ix==2)then        !negative x
         idim=1
         if(u(idim,idx)>zero)call zero_flow_dir(idx,idim)
       endif
       if(ix==nx-1.or.ix==nx)then   !positve x
         idim=1
         if(u(idim,idx)<zero)call zero_flow_dir(idx,idim)
       endif
       if(iy==1.or.iy==2)then       !negative y
         idim=2
         if(u(idim,idx)>zero)call zero_flow_dir(idx,idim)
       endif
       if(iy==ny.or.iy==ny-1)then    !positive y
         idim=2
         if(u(idim,idx)<zero)call zero_flow_dir(idx,idim)
       endif
!       if(iz==1.or.iz==2)then       !negative z
!         idim=3
!         if(u(idim,idx)>zero)call zero_flow_dir(idx,idim)
!       endif
!       if(iz==nz.or.iz==nz-1)then    !positive z
!         idim=3
!         if(u(idim,idx)<zero)call zero_flow_dir(idx,idim)
!       endif
#endif /* end ifdef CONTROLVALVES */
!! FINISH for control valve

#ifdef WINDTUNNEL
       if(ix==nx.or.ix==nx-1)then
         u(1,idx)=-vflow
         u(2,idx)=zero
         u(3,idx)=zero
         cons(1,idx)=rhoflow
         cons(2,idx)=u(1,idx)*cons(1,idx)
         cons(3,idx)=zero
         cons(4,idx)=zero
         call get_gamma_from_tk(eps,cons(1,idx),tkflow,mucloc,adloc)
         adindx(idx)=adloc
         p(idx)=scale%rgas*cons(1,idx)*tkflow/mucloc
         cons(5,idx)=eps+half*cons(1,idx)*u(1,idx)**2
       endif
#endif /* end ifdef WINDTUNNEL */

#ifdef SHOCKDIFF
       if(iy>ny/2.and.(ix==1.or.ix==2))then
         cons(1,idx)=5.0293806477913563d0
         p(idx)=21.471035714285712d0
         u(1,idx)=4.0779469548133598d0
         cons(2,idx)=u(1,idx)*cons(1,idx) 
         cons(5,idx)=p(idx)/(adindx(idx)-1.0)+0.5*cons(2,idx)**2/cons(1,idx)
       endif
#endif /* end ifdef SHOCKDIFF */


  enddo
!$OMP ENDDO
!$OMP BARRIER
!
!
#ifdef EXTRAANCHORS
!$OMP DO SCHEDULE(STATIC) PRIVATE(idx)
  do igrid=1,nanchor
    idx=indx_anchor(igrid)
    cons(1,idx)=small_rho
    cons(2,idx)=zero
    cons(3,idx)=zero
    cons(4,idx)=zero
    cons(5,idx)=small_eps
  enddo
!$OMP ENDDO
#endif /* end ifdef EXTRAANCHORS */
!
!

 end subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate cubic spline for softening potential.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 real(pre) function wspline3(x)
  real(pre), intent(in)::x
   if (x<one)then
     wspline3= x*(1.5d0*x-0.75d0*x*x)
   elseif (x<two) then
     wspline3=one-0.25d0*(two-x)**3
   else
     wspline3=one
   endif
 end function
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate cubic spline for softening potential.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 real(pre) function derivwspline3(x)
  real(pre), intent(in)::x
   if (x<one)then
     derivwspline3= x*(3d0-2.25d0*x)
   elseif (x<two) then
     derivwspline3= 0.75*(two-x)**2
   else
     derivwspline3=zero
   endif
 end function
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate exponetial Plummer for softening potential.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 real(pre) function wexpplum(x)
  real(pre), intent(in)::x
  
  wexpplum = x/(x+exp(-x))

 end function
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate derivative of exponetial Plummer for softening potential.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 real(pre) function derivwexpplum(x)
  real(pre), intent(in)::x
  
  derivwexpplum = one/(x+exp(-x)) -x/(x+exp(-x))**2 * (one-exp(-x))

 end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate gravitational force from grid.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 subroutine calc_gforce_helper(isten,b,igrid,im,ip)
  use grid_commons, only : grid
  integer::isten(5),b(6),igrid,im,ip
  isten(3)=igrid
  isten(4)=b(ip)
  isten(5)=grid(b(ip))%ineigh(ip)
  isten(2)=b(im)
  isten(1)=grid(b(im))%ineigh(im)
 end subroutine

 subroutine calc_gforce()
  use grid_commons
  integer::igrid,b(6),sten(5)

!$OMP DO SCHEDULE(STATIC) &
!$OMP&PRIVATE(b,sten)
  do igrid=1,ngrid
    if(grid(igrid)%boundary>0)cycle
    call get_boundary_wb(igrid,b)
    
!    call calc_gforce_helper(sten,b,igrid,4,3)
!    gforce(1,igrid)=-(-phi(sten(5))+eight*phi(sten(4))-eight*phi(sten(2))+phi(sten(1)))/(12d0*dx)
!    call calc_gforce_helper(sten,b,igrid,2,1)
!    gforce(2,igrid)=-(-phi(sten(5))+eight*phi(sten(4))-eight*phi(sten(2))+phi(sten(1)))/(12d0*dy)
!    call calc_gforce_helper(sten,b,igrid,6,5)
!    gforce(3,igrid)=-(-phi(sten(5))+eight*phi(sten(4))-eight*phi(sten(2))+phi(sten(1)))/(12d0*dz)

    gforce(1,igrid)=-(phi(b(3))-phi(b(4)))/(two*dx)
    gforce(2,igrid)=-(phi(b(1))-phi(b(2)))/(two*dy)
    gforce(3,igrid)=-(phi(b(5))-phi(b(6)))/(two*dz)
  enddo
!$OMP ENDDO
 end subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get polar angle. Like atan2.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 subroutine get_polar(r,an,x,y)
  real(pre),intent(in)::x,y
  real(pre),intent(out)::r,an
  
  r=sqrt(x*x+y*y)
  if(x==zero)then
     an=pi*half 
     if(y<zero)an=1.5d0*pi
  elseif(y==zero)then
     an=zero
     if(x<zero)an=pi
  else
     an=atan(y/x)
     if(x<zero)then
       an=an+pi
     elseif(y<zero)then
       an=an+two*pi
     endif
  endif

  return
 end subroutine get_polar
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Flux limiter for hydro fluxing.  It is set to MINMOD always,
! but can be changed.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 real(pre) function flux_limiter(a,b,c,d)
   real(pre),intent(in)::a,b,c,d
   real(pre)::r,den
!
!
#define MINMOD
!
!
   flux_limiter=zero
!
!
#ifdef VANLEER
     den=(c-d)
     if(den==zero)then 
       flux_limiter=two
       return
     endif
     r=(a-b)/den
     flux_limiter = (r+abs(r))/(one+abs(r)) 
#endif /* end ifdef VANLEER */
!
!
#ifdef VANALDABADA1
     den=(c-d)
     if(den==zero)then 
       flux_limiter=one
       return
     endif
     r=(a-b)/den
     if(r<zero)r=zero
     flux_limiter = (r*r + r)/(r*r + one)
#endif /* end ifdef VANALDABADA1 */
!
!
#ifdef SUPERBEE
     den=(c-d)
     if(den==zero)then 
       flux_limiter=two
       return
     endif
     r=(a-b)/den
     if(r<zero)then
       flux_limiter=zero
     elseif(r<half)then
       flux_limiter=two*r
     elseif(r<one)then
       flux_limiter=one
     elseif(r<two)then
       flux_limiter=r
     else
       flux_limiter=two
     endif
#endif /* end ifdef SUPERBEE */
!
!
#ifdef MINMOD 
     den=(c-d)
     if(den==zero)then 
       flux_limiter=one
       return
     endif
     r=(a-b)/den
     flux_limiter = max(zero,min(r,one))
#endif /* end ifdef MINMOD (DEFAULT */
!
!
 end function
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Upwind fluxing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 real(pre) function upwind(f2,f1,f0,fm1,v,ds,t)
   real(pre),intent(in)::f2,f1,f0,fm1,v,ds,t
   real(pre)::limit,nu

   nu=v*t/ds
   if (v<zero)then
     limit=flux_limiter(f2,f1,f1,f0)
     upwind = f1-limit*half*(one+nu)*(f1-f0)
   else
     limit=flux_limiter(f0,fm1,f1,f0)
     upwind = f0+limit*half*(one-nu)*(f1-f0)
   endif

 end function
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Left and Right states for central differencing.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 real(pre) function genminmod(x1,x2,x3)
    real(pre)::x1,x2,x3,tmp
    if((x1>zero) .and. ((x2>zero) .and. (x3>zero)))then
        tmp=min(x2,x3)
        genminmod=min(x1,tmp)
    else if((x1<zero) .and. ((x2<zero) .and. (x3<zero)))then
        tmp=max(x2,x3)
        genminmod=max(x1,tmp)
    else
        genminmod=zero
    endif
 end function

 real(pre) function vanleer(a,b)
    real(pre)::a,b,absa,absb,signa,signb
    absa=abs(a)
    absb=abs(b)
    if(absa==0)then
      signa=0d0
    else
      signa=a/absa
    endif
    if(absb==0)then
      signb=0d0
    else
      signb=b/absb
    endif
    vanleer=absa*absb/(absa+absb+1d-9)*(signa+signb)
 end function

 subroutine left_right_states(cm1,c0,c1,c2,c_l,c_r)
    real(pre)::cm1,c0,c1,c2,c_l,c_r,tt=SLOPE_THETA,ufaceLL,ufaceRR
    c_l=c0+half*genminmod(tt*(c0-cm1),half*(c1-cm1),tt*(c1-c0))
    c_r=c1-half*genminmod(tt*(c1-c0),half*(c2-c0),tt*(c2-c1))
 end subroutine

 subroutine left_right_states_v(cm1,c0,c1,c2,ccL,ccR,c_l,c_r)
    real(pre)::cm1,c0,c1,c2,c_l,c_r,tt=SLOPE_THETA,ccL,ccR
    c_l=c0+half*(one-ccL)*genminmod(tt*(c0-cm1),half*(c1-cm1),tt*(c1-c0))
    c_r=c1-half*(one+ccR)*genminmod(tt*(c1-c0),half*(c2-c0),tt*(c2-c1))
 end subroutine

 real(pre) function slope(c0,c1,c2)
    real(pre)::c0,c1,c2,tt=SLOPE_THETA
    slope=genminmod(tt*(c1-c0),half*(c2-c0),tt*(c2-c1))
 end function slope


!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Another function like atan2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      real(pre) function getAngle(x,y) 
       real(pre),intent(in):: x,y
       if(x==0d0)then
        if(y==0d0)then
          getAngle=0d0
        elseif(y>0d0)then
          getAngle=pi*0.5d0
        else
          getAngle=pi*1.5d0
        endif 
       else
        getAngle=atan(y/x)
        if (y>0d0)then
         if(x<0d0)getAngle=getAngle+pi 
        else
         if(x<0d0)then
          getAngle=getAngle+pi
         else
          getAngle=getAngle+2.d0*pi
         endif
        endif
       endif  
       return 
      end function getAngle
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! radial part of Kepler
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      real(pre) function radialCoorKep(a,ee,tAnom)
       real(pre),intent(in)::a,ee,tAnom
       radialCoorKep=a*(1d0-ee**2)/(1d0+ee*cos(tAnom))
       return
      end function radialCoorKep
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
! Solve Kepler's problem
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine trueAnomaly(tt,mu,a,ee,theta,eAnom)
       real(pre),intent(in)::tt,mu,a,ee
       real(pre),intent(out)::theta,eAnom
       real(pre)::calc,errMag,err,orbits,val
       real(pre)::tol_kep=1d-3

       val=tt*sqrt(mu/a**3)
       orbits=dble(int(val/(2.d0*pi)))
       !print *,"Anomaly r,val,orbits,mu", r,val,orbits,mu
       val=val-orbits*2d0*pi
       if(val<1d-32)then ! don't do work unless we must
        theta=0d0;eAnom=0d0;return
       endif

       eAnom=1d0
       errMag=1d0
       do while(errMag>tol_kep)
        calc=eAnom-ee*sin(eAnom)
        err=(val-calc)/val
        errMag=abs(err)
        eAnom=eAnom+0.6d0*min(errMag,1d0)*eAnom*err/errMag
        !print *, eAnom,val,err,time
       end do
       theta=two*getAngle(sqrt(one-ee)*cos(eAnom*half),sqrt(one+ee)*sin(eAnom*half))
       return
      end subroutine trueAnomaly
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This routine calculates the position of a wide-orbit binary.
! For the position of the binary, only consider r and phi because...
! Work in progress and is not used.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
subroutine getBinaryPosition(tt,rBin,thetaBin,omegaBin, &
     &  massDisk,muBin,abin,ebin,massbin)
      implicit none
      real(pre),intent(in)::tt,massDisk,abin,ebin
      real(pre),intent(out)::rBin,thetaBin,omegaBin,muBin
      real(pre)::eAnom,massbin,tBin

      muBin=massDisk+massBin
      omegaBin=sqrt((muBin)/(aBin**3))
      tBin=tt

      call trueAnomaly(tBin,muBin,aBin,eBin,thetaBin,eAnom)
      rBin=radialCoorKep(aBin,eBin,thetaBin)
      return
      end subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Maximum magnitude function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 real(pre) function magmax(a,b)
   real(pre),intent(in)::a,b

     if(abs(a)>abs(b))then
        magmax=a
     else
        magmax=b
     endif
 end function
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Maximum magnitude of three variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 real(pre) function magmax3(a,b,c)
   real(pre),intent(in)::a,b,c
     magmax3=magmax(magmax(a,b),c)
 end function
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Spectral radius function for Jacobian solver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 real(pre) function get_spectral_radius_from_flux_jacobian(var,v,igrid,idim)
   real(pre),intent(in)::var(:,:),v(:,:)
   real(pre)::jac(3,3),l(3)
   integer::igrid,b(6),idim
   logical::flag,active

   call get_boundary(igrid,b)
   call flux_jacobian(var,v,idim,b,jac)
   call matrix3_eigen(l,jac,flag)
   if(flag)then
     get_spectral_radius_from_flux_jacobian=abs(l(1))
   else
     get_spectral_radius_from_flux_jacobian=abs(magmax3(l(1),l(2),l(3)))
   endif
 end function
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Eigen values of a 3X3 matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 subroutine matrix3_eigen(l,a,flag)
   real(pre),intent(in)::a(3,3)
   real(pre),intent(out)::l(3)
   real(pre)::a1,a2,a3,q,r,d,s,t
   logical::flag

   a1=-(a(3,3)+a(2,2)+a(1,1))
   a2=(a(2,2)*a(3,3)+a(1,1)*a(3,3)+a(1,1)*a(2,2)-a(2,1)*a(1,2)-a(3,2)*a(2,3)-a(3,1)*a(1,3))
   a3=-(a(1,1)*a(2,2)*a(3,3)+a(2,1)*a(3,2)*a(1,3)+a(3,1)*a(1,2)*a(2,3)-a(1,1)*a(3,2)*a(2,3) &
       -a(2,1)*a(1,2)*a(3,3)-a(3,1)*a(2,2)*a(1,3))


   q=(three*a2-a1**2)/nine
   r=(nine*a1*a2-27d0*a3-two*a1**3)/54d0
   d=q**3+r**3
   flag=.false.
   
   if(d>zero)then
     s=(r+sqrt(d))**(one/three)
     t=(r+sqrt(d))**(one/three)
     flag=.true.
     l(1)=s+t-a1/three
     l(2)=zero
     l(3)=zero
   elseif(d==zero)then
     s=r**(one/three);t=s
     l(1)=s+t-a1/three
     l(2)=-half*(s+t)-a1/three
     l(3)=l(2)
   else
      t=acos(max(min(r/sqrt(-q**3),one),-one))
      l(1)=two*sqrt(-q)*cos(t/three)-a1/three
      l(2)=two*sqrt(-q)*cos(t/three+two*pi/three)-a1/three
      l(3)=two*sqrt(-q)*cos(t/three+four*pi/three)-a1/three
   endif

 end subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Jacobian
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 subroutine flux_jacobian(var,v,idim,b,jac)
   real(pre)::var(:,:),v(:,:)
   integer,intent(in)::idim
   real(pre),intent(out)::jac(3,3)
   real(pre)::df,du
   integer::b(6)
   
   df=var(idim,b(3))*v(1,b(3))-var(idim,b(4))*v(1,b(4))
   du=var(idim,b(3))-var(idim,b(4))
   jac(1,1)=zero
   if(.not.du==zero)jac(1,1)=df/du

   du=var(idim,b(1))-var(idim,b(2))
   jac(1,2)=zero
   if(.not.du==zero)jac(1,2)=df/du
 
   du=var(idim,b(5))-var(idim,b(6))
   jac(1,3)=zero
   if(.not.du==zero)jac(1,3)=df/du
   
   df=var(idim,b(1))*v(2,b(1))-var(idim,b(2))*v(2,b(2))
   du=var(idim,b(3))-var(idim,b(4))
   jac(2,1)=zero
   if(.not.du==zero)jac(2,1)=df/du

   du=var(idim,b(1))-var(idim,b(2))
   jac(2,2)=zero
   if(.not.du==zero)jac(2,2)=df/du
 
   du=var(idim,b(5))-var(idim,b(6))
   jac(2,3)=zero
   if(.not.du==zero)jac(2,3)=df/du

   df=var(idim,b(5))*v(3,b(5))-var(idim,b(6))*v(3,b(6))
   du=var(idim,b(3))-var(idim,b(4))
   jac(3,1)=zero
   if(.not.du==zero)jac(3,1)=df/du

   du=var(idim,b(1))-var(idim,b(2))
   jac(3,2)=zero
   if(.not.du==zero)jac(3,2)=df/du
 
   du=var(idim,b(5))-var(idim,b(6))
   jac(3,3)=zero
   if(.not.du==zero)jac(3,3)=df/du

 end subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Suppress grid drift by forcing the COM of the entire domain to be
! at the center of the grid.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 subroutine suppress_drift(xcom,ycom,zcom)
  use grid_commons
  use parameters
  real(pre),intent(in)::xcom,ycom,zcom
  real(pre)::vxcom,vycom,vzcom
  integer::igrid
!$OMP MASTER
  if(comfirst)then
    xcom0=xcom
    ycom0=ycom
    zcom0=zcom
    comfirst=.false.
  endif
!$OMP END MASTER
!$OMP BARRIER
  if(dt>zero)then
   vxcom=zero;vycom=zero;vzcom=zero
   if(abs(xcom)>abs(xcom0)) vxcom=xcom/dt*0.1d0 ! v*com are private
   if(abs(ycom)>abs(ycom0)) vycom=ycom/dt*0.1d0 ! v*com are private
   if(abs(zcom)>abs(zcom0)) vzcom=zcom/dt*0.1d0 ! v*com are private
!$OMP DO SCHEDULE(STATIC)
   do igrid=1,ngrid
      cons(2,igrid)=cons(2,igrid)-vxcom*cons(1,igrid)
      cons(3,igrid)=cons(3,igrid)-vycom*cons(1,igrid)
      cons(4,igrid)=cons(4,igrid)-vzcom*cons(1,igrid)
   enddo
!$OMP ENDDO
  endif
!$OMP MASTER
 xcom0=xcom
 ycom0=ycom
 zcom0=zcom
!$OMP END MASTER

 end subroutine
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Check how we are doing with conserving quantities when applicable. 
! Check drift of COM.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 subroutine conservation_diagnostic()
   use grid_commons
   real(pre)::x,y,z,angle
   integer::igrid
!$OMP MASTER
   mass=zero
   eint=zero
   ekin=zero
   etot=zero
   egra=zero
   amom=zero
   momx=zero
   momy=zero
   xcom=zero
   ycom=zero
   zcom=zero
!$OMP END MASTER 
!$OMP BARRIER
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:mass,ekin,eint,amom,momx)        &
!$OMP&REDUCTION(+:momy,momz,etot,egra,xcom,ycom,zcom) PRIVATE(angle,x,y,z)
   do igrid=1,ngrid
    if(grid(igrid)%boundary>0)cycle
    x=grid(igrid)%x
    y=grid(igrid)%y
    z=grid(igrid)%z
    angle=atan2(y,x)
    amom=amom+(cons(3,igrid)*cos(angle)-cons(2,igrid)*sin(angle))&
        *sqrt(x**2+y**2)*dx*dy*dz
    momx=momx+cons(2,igrid)*dz*dx*dy
    momy=momy+cons(3,igrid)*dz*dx*dy
    momz=momz+cons(4,igrid)*dz*dx*dy
    mass=mass+cons(1,igrid)*dz*dx*dy
    xcom=xcom+x*cons(1,igrid)*dz*dx*dy
    ycom=ycom+y*cons(1,igrid)*dz*dx*dy
    zcom=zcom+z*cons(1,igrid)*dz*dx*dy
    etot=etot+(cons(5,igrid))*dz*dx*dy
    ekin=ekin+(cons(2,igrid)**2+cons(3,igrid)**2+cons(4,igrid)**2) &
        /(two*cons(1,igrid))*dx*dy*dz
    egra=egra+half*(cons(1,igrid)*phi(igrid))*dx*dy*dz
   enddo
!$OMP ENDDO
!$OMP BARRIER
!
!
#ifdef VERBOSE
!$OMP MASTER
  print *," Total mass   is ",mass,time
  print *," Total momx,y,z is ",momx,momy,momz,time
  print *," Total amom   is ",amom,time
  print *," Total energy is ",etot+egra,etot,etot-ekin,ekin,egra,time
  print *,"COM ",time,xcom/mass,ycom/mass,zcom/mass
!$OMP END MASTER
#endif /* end ifdef VERBOSE */
!
!
 end subroutine

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate total energy. Used for correction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 subroutine get_etot(e,m,f)
   use grid_commons
   real(pre)::e,m,f
   integer::igrid
!
!$OMP DO SCHEDULE(STATIC)REDUCTION(+:e,m,f)
   do igrid=1,ngrid
    if(grid(igrid)%boundary>0)cycle
    m=m+cons(1,igrid)*dx*dy*dz
    f=f+cons(1,igrid)*dx*dy*dz*abs(u(1,igrid)*gforce(1,igrid)+u(2,igrid)*gforce(2,igrid)+&
      u(3,igrid)*gforce(3,igrid))
    e=e+(cons(5,igrid)+half*(cons(1,igrid)*phi(igrid)))*dx*dy*dz
   enddo
!$OMP ENDDO
!
 end subroutine

end module

