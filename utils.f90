module utils
 use parameters
 use derived_types
 implicit  none

 real(pre),save::xcom0,ycom0,zcom0
 logical,save::comfirst=.true.

 contains

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

 subroutine get_boundary_anchor(igrid,b,active_grid)
  use derived_types
  use grid_commons
  integer, intent(in)::igrid
  integer, intent(out)::b(6)
  integer::idx
  logical::active_grid

  active_grid=.true.
  if(grid(igrid)%anchor)then
    active_grid=.false.
    return
  endif
  
  do idx=1,6
   b(idx)=grid(igrid)%ineigh(idx)
   if(grid(b(idx))%anchor)b(idx)=igrid
  enddo

 end subroutine get_boundary_anchor


 subroutine get_boundary(igrid,b,active_grid)
  use derived_types
  use grid_commons
  integer, intent(in)::igrid
  integer, intent(out)::b(6)
  integer::idx
  logical::active_grid

  active_grid=.true.
  if(grid(igrid)%boundary)then
    active_grid=.false.
    return
  endif
  
  do idx=1,6
   b(idx)=grid(igrid)%ineigh(idx)
   if(grid(b(idx))%boundary)b(idx)=igrid
  enddo

 end subroutine get_boundary

 subroutine get_boundary_wb_loc(g,b,active_grid)
  use derived_types
  integer, intent(out)::b(6)
  integer::idx
  type(gridcell), intent(in)::g
  logical::active_grid

  active_grid=.true.
  if(g%boundary)then
    active_grid=.false.
    return
  endif
  
  do idx=1,6
   b(idx)=g%ineigh(idx)
  enddo

 end subroutine get_boundary_wb_loc


 subroutine get_boundary_wb(igrid,b,active_grid)
  use derived_types
  use grid_commons
  integer, intent(in)::igrid
  integer, intent(out)::b(6)
  integer::idx
  logical::active_grid

  active_grid=.true.
  if(grid(igrid)%boundary)then
    active_grid=.false.
    return
  endif
  
  do idx=1,6
   b(idx)=grid(igrid)%ineigh(idx)
  enddo

 end subroutine get_boundary_wb

 subroutine get_nearest_8(x,y,z,neigh)
  use grid_commons
  integer::igrid,igridt
  integer,intent(out)::neigh(8)
  integer::bb(6)
  integer::bt(6)
  real(pre),intent(in)::x,y,z
  real(pre)::delz
  logical::active 

  igrid=get_grid_indx(x,y,z)
  call get_boundary(igrid,bb,active)

  delz=z-grid(igrid)%z
  if(delz<zero)then
     igrid=bb(6)
     call get_boundary(igrid,bb,active)
  endif
  igridt=bb(5)
  call get_boundary(igridt,bt,active)
  
  neigh(5)=bb(5)
  neigh(6)=igrid

  if(x-grid(igrid)%x<zero )then
    if(y-grid(igrid)%y<zero) then
      neigh(1)=igrid
      neigh(2)=bb(4)
      neigh(3)=grid(bb(4))%ineigh(2)
      neigh(4)=bb(2)
      neigh(5)=igridt
      neigh(6)=bt(4)
      neigh(7)=grid(bb(4))%ineigh(2)
      neigh(8)=bt(2)
    else
      neigh(1)=igrid
      neigh(2)=bb(1)
      neigh(3)=grid(bb(1))%ineigh(4)
      neigh(4)=bb(4)
      neigh(5)=igridt
      neigh(6)=bt(1)
      neigh(7)=grid(bb(1))%ineigh(4)
      neigh(8)=bt(4)
    endif
  else
    if(y-grid(igrid)%y<zero) then
      neigh(1)=igrid
      neigh(2)=bb(2)
      neigh(3)=grid(bb(2))%ineigh(3)
      neigh(4)=bb(3)
      neigh(5)=igridt
      neigh(6)=bt(2)
      neigh(7)=grid(bb(2))%ineigh(3)
      neigh(8)=bt(3)
    else
      neigh(1)=igrid
      neigh(2)=bb(3)
      neigh(3)=grid(bb(3))%ineigh(1)
      neigh(4)=bb(1)
      neigh(5)=igridt
      neigh(6)=bt(3)
      neigh(7)=grid(bb(3))%ineigh(1)
      neigh(8)=bt(1)
    endif
   endif
  
 end subroutine

 subroutine calc_gforce()
  use grid_commons
  integer::igrid,b(6)
  logical::active

!$OMP DO SCHEDULE(STATIC) &
!$OMP&PRIVATE(b,active)
  do igrid=1,ngrid
    call get_boundary_wb(igrid,b,active)
    if(.not.active)cycle

    gforce(1,igrid)=-(phi(b(3))-phi(b(4)))/(two*dx)
    gforce(2,igrid)=-(phi(b(1))-phi(b(2)))/(two*dy)
    gforce(3,igrid)=-(phi(b(5))-phi(b(6)))/(two*dz)
  enddo
!$OMP ENDDO

 end subroutine

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

 real(pre) function flux_limiter(a,b,c,d)
   real(pre),intent(in)::a,b,c,d
   real(pre)::r,den
#define MINMOD

   flux_limiter=zero
#ifdef VANLEER
     den=(c-d)
     if(den==zero)then 
       flux_limiter=two
       return
     endif
     r=(a-b)/den
     flux_limiter = (r+abs(r))/(one+abs(r)) 
#endif
#ifdef VANALDABADA1
     den=(c-d)
     if(den==zero)then 
       flux_limiter=one
       return
     endif
     r=(a-b)/den
     if(r<zero)r=zero
     flux_limiter = (r*r + r)/(r*r + one)
#endif
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
#endif
#ifdef MINMOD 
     den=(c-d)
     if(den==zero)then 
       flux_limiter=one
       return
     endif
     r=(a-b)/den
     flux_limiter = max(zero,min(r,one))
#endif
 end function

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

 subroutine left_right_states(f2,f1,f0,fm1,left,right)
   real(pre),intent(in)::f2,f1,f0,fm1
   real(pre),intent(out)::left,right
   real(pre)::limit

   limit=flux_limiter(f1,f0,f2,f1)
   right = f1-limit*half*(f2-f1)

   limit=flux_limiter(f0,fm1,f1,f0)
   left =  f0+limit*half*(f1-f0)

 end subroutine

!.....Translate cvmgp and cvmgz into equivalent if statements
 real(pre) FUNCTION SLOPE(a,b,c)
  real(pre):: a,b,c,tmp,tmm
  slope=zero
  if ((b-c)*(a-b) < zero) return
   tmm=a-c
  if (tmm.eq.0.0) return
   slope=(b-c)*(a-b)/tmm*two
  return
 end  function 

!*******************************************************************************
!.........and van leer interpolation
 real(pre) FUNCTION VLI(qf,qb,slopf,slopb,vel,dt,ds)
  real(pre):: qf,qb,slopf,slopb,vel,dt,ds
  if(vel.lt.zero) then 
    vli=qf-(ds+vel*dt)*slopf/(two*ds)
  else 
    vli=qb+(ds-vel*dt)*slopb/(two*ds)
  endif
  return
 end  function

!!!! From CHYMERA
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

      !!!!!!!!!!!!!!!!!!!!!!!
      ! radial part of Kepler
      !!!!!!!!!!!!!!!!!!!!!!!
 
      real(pre) function radialCoorKep(a,ee,tAnom)
       real(pre),intent(in)::a,ee,tAnom
       radialCoorKep=a*(1d0-ee**2)/(1d0+ee*cos(tAnom))
       return
      end function radialCoorKep
    
      !!!!!!!!!!!!!!!!!!!!!!!!
      ! Solve Kepler's problem
      !!!!!!!!!!!!!!!!!!!!!!!! 
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This routine calculates the position of a wide-orbit binary.
!!!
!!! For the position of the binary, only consider r and phi because
!!! we are still in mirror symmetry land for the disk. 
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
      !!print *, "GetBINARY",tBin,muBin,omegaBin,rBin,thetaBin,massDisk
      return
      end subroutine

!
 real(pre) function magmax(a,b)
   real(pre),intent(in)::a,b

     if(abs(a)>abs(b))then
        magmax=a
     else
        magmax=b
     endif
 end function


 real(pre) function magmax3(a,b,c)
   real(pre),intent(in)::a,b,c
     magmax3=magmax(magmax(a,b),c)
 end function

 real(pre) function get_spectral_radius_from_flux_jacobian(var,v,igrid,idim)
   real(pre),intent(in)::var(:,:),v(:,:)
   real(pre)::jac(3,3),l(3)
   integer::igrid,b(6),idim
   logical::flag,active

   call get_boundary(igrid,b,active)
   call flux_jacobian(var,v,idim,b,jac)
   call matrix3_eigen(l,jac,flag)
   if(flag)then
     get_spectral_radius_from_flux_jacobian=abs(l(1))
   else
     get_spectral_radius_from_flux_jacobian=abs(magmax3(l(1),l(2),l(3)))
   endif
 end function

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

end module

