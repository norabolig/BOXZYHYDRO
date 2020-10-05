!
! Guts for the particle-gas coupling.
!
module pdrag
 use parameters
 use derived_types
 use grid_commons
 use utils
 use eos, only: get_gamma_from_tk
 implicit none
!
 contains
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Simple function for the sign
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 integer function get_sign(a)
    real(pre)::a
    if (a==0)then
       get_sign=0
    else if(a<0) then
       get_sign=-1
    else
       get_sign=1
    endif
 end function
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate the amount of drift in the analytic limit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 real(pre) function calc_drift(f,ts)
    real(pre)::f,ts
    calc_drift=f*ts
    return
 end function
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Main work horse.  Find the drag and return that forces.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 subroutine get_drag(dfx,dfy,dfz,tg,pg,dg,x,y,z,vx,vy,vz,d, &
                           rhoa,asize,ig,w,t,fx,fy,fz,alpha,de)
  integer,intent(in)::ig(8)
  real(pre),intent(in)::t,rhoa,asize,vx,vy,vz,x,y,z,d,w(8),fx,fy,fz
  real(pre),intent(out)::tg,pg,dg,de
  real(pre)::azimuth,dumeps,dumm,momx,momy,momz,vxn,vyn,vzn,weight
  real(pre)::vgx,vgy,vgz,delvx,delvy,delvz,cs,kn,dfx,dfy,dfz,vdrift,e0,e1
  real(pre)::delvx_new,delvy_new,delvz_new,ggam,alpha,exp_talpha,ts,tdrag
  real(pre)::vcirc,r,dg0
  type(units)::scl

  integer::iter

  call get_units(scl)

  vgx=zero;vgy=zero;vgz=zero
  dg=zero;tg=zero;pg=zero
  momx=zero
  momy=zero
  momz=zero
  ggam=zero
  do iter=1,8
   dg =dg +cons(1,ig(iter))*w(iter)
   vgx=vgx+cons(2,ig(iter))*w(iter)
   vgy=vgy+cons(3,ig(iter))*w(iter)
   vgz=vgz+cons(4,ig(iter))*w(iter)
   tg =tg +muc_array(ig(iter))*w(iter)
   pg =pg +p(ig(iter))*w(iter)
   ggam=ggam+adindx(ig(iter))*w(iter)
  enddo
  momx=vgx+vx*d
  momy=vgy+vy*d
  momz=vgz+vz*d
  vgx=vgx/dg
  vgy=vgy/dg
  vgz=vgz/dg
  tg=pg/(scl%rgas*dg)*tg
  delvx=vgx-vx
  delvy=vgy-vy
  delvz=vgz-vz
  cs=sqrt(ggam*pg/dg)
  e0=half*( (vgx**2+vgy**2+vgz**2)*dg + (vx**2+vy**2+vz**2)*d )

  if (nz<6)then
     vcirc = abs(x*vgy-y*vgx)/sqrt(x*x+y*y)
     dg0=dg*vcirc/(two*cs*sqrt(x*x+y*y))
     !print * , d,dg
  else
    dg0=dg
  endif

  call get_gamma_from_tk(dumeps,dg,tg,dumm,ggam)
  kn=half*dumm*mp/scl%mass/(dg0*pi*4.43d-43*asize) !4.44d-43->1d-16cm^2

  if(tg<tg_immediate_couple.or.asize>=a_sublimate_limit)then

    call update_dv(delvx_new,delvx,kn, &
                   cs,dg0,rhoa,asize,t,alpha)
    call update_dv(delvy_new,delvy,kn, &
                   cs,dg0,rhoa,asize,t,alpha)
    call update_dv(delvz_new,delvz,kn, &
                   cs,dg0,rhoa,asize,t,alpha)
!
!
#ifdef TURN_ON_ONLY_FOR_DEBUGGING_PURPOSES
    if(IsNan(delvx_new))then 
     print *,"dg,vgx,vgy,vgz,tg,pg,vx,vy,vz,cs,delvx_new,delvy_new,delvz_new"
     print *, dg,vgx,vgy,vgz,tg,pg,vx,vy,vz,cs,delvx_new,delvy_new,delvz_new
     do iter=1,8
      print *,iter,cons(1,ig(iter))*w(iter),w(iter),ig(iter),grid(ig(iter))%boundary
     enddo
     do iter=1,8
      print *,iter,cons(1:5,ig(iter))
     print *,"Check location ",grid(ig(iter))%ix,grid(ig(iter))%iy,grid(ig(iter))%iz
     enddo
      stop
    endif
#endif
!
!
    ts=one/alpha
    exp_talpha=exp(-two)
    if(t>ts*exp_talpha)then
       weight=(ts*exp_talpha/t)**2
       delvx_new=delvx_new*weight-calc_drift(fx,ts)*(one-weight)
       delvy_new=delvy_new*weight-calc_drift(fy,ts)*(one-weight)
       delvz_new=delvz_new*weight-calc_drift(fz,ts)*(one-weight)
    endif
!
  else
    delvx_new=zero
    delvy_new=zero
    delvz_new=zero
  endif
!
  vxn=(momx+delvx_new*d)/(d+dg)-delvx_new
  vyn=(momy+delvy_new*d)/(d+dg)-delvy_new
  vzn=(momz+delvz_new*d)/(d+dg)-delvz_new


  e1=half*( ((vxn+delvx_new)**2+(vyn+delvy_new)**2+(vzn+delvz_new)**2)*dg &
    + ((vxn)**2+(vyn)**2+(vzn)**2)*d )

  dfx=(vxn-vx)/t
  dfy=(vyn-vy)/t
  dfz=(vzn-vz)/t

  de=(e0-e1)/(d+dg) ! assuming the dust takes all of the missing energy
!  de=vxn+delvx_new-vgx

 end subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update the velocity difference between particles and gas due to drag.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 subroutine update_dv(dvnew,dv,kn,c,dg,rho0,asize,dt,alpha)
  real(pre),intent(in)::dv,kn,dg,c,rho0,asize,dt
  real(pre),intent(out)::dvnew
  real(pre)::magdv,mach,Re,beta2,beta,alpha,bb,dv_eps
  real(pre)::dv_sto,pq,pqinv

  magdv=abs(dv)

  if(dv==zero)then
    dvnew=zero
    return
  endif

  mach=magdv/c
  Re=three*sqrt(pi/eight)*mach/kn
  beta2=128d0*c**2/(nine*pi)
  beta =sqrt(beta2)
  alpha=dg*c/(rho0*asize)*sqrt(eight/pi)
  bb=beta+sqrt(dv**2+beta2)
!
!***
! First do the Epstein limit
!***
!
  dv_eps=two*bb*dv*beta*exp(-alpha*dt) &
        /(bb*bb-dv*dv*exp(-two*alpha*dt))
!
!***
! Now work on the Stokes regime
!***
!
  if(Re<=500d0)then
     pq=0.687d0
     pqinv=one/pq
     beta=0.15d0*(three*sqrt(pi/eight)/(c*kn))**pq
     dv_sto=exp(-three*alpha*kn*dt)
     if(dv_sto>zero)then
       dv_sto=dv_sto/( (beta*magdv**pq+one)*magdv**(-pq) &
             - beta*exp(-three*pq*alpha*kn*dt) )**pqinv
     endif
     dv_sto=dv_sto*magdv/dv
  elseif(Re<=1500d0)then
     pq=2.4d0
     pqinv=one/pq
     beta=3.96d-6*(three*sqrt(pi/eight)/(c*kn) )**pq
     dv_sto=(magdv**(-pq)+7.2d0*kn*alpha*beta*dt)**(-pqinv)*dv/magdv
  else
     dv_sto=dv/(one+magdv*0.99d0*sqrt(pi/eight)*alpha/c*dt)
  endif

   dvnew=dv_eps*(three*kn)**2/((three*kn)**2+one)+dv_sto/((three*kn)**2+one)

 end subroutine

end module

