module pdrag
 use parameters
 use derived_types
 use grid_commons
 use utils
 use eos, only: get_gamma_from_tk
 implicit none


 contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

 real(pre) function calc_drift(f,ts)
    real(pre)::f,ts
    calc_drift=f*ts
    return
 end function

 subroutine get_drag(dfx,dfy,dfz,tg,pg,dg,x,y,z,vx,vy,vz,d, &
                           rhoa,asize,ig,w,t,fx,fy,fz,alpha)
  integer,intent(in)::ig(8)
  real(pre),intent(in)::t,rhoa,asize,vx,vy,vz,x,y,z,d,w(8),fx,fy,fz
  real(pre),intent(out)::tg,pg,dg
  real(pre)::azimuth,dumeps,dumm,momx,momy,momz,vxn,vyn,vzn,weight
  real(pre)::vgx,vgy,vgz,delvx,delvy,delvz,cs,kn,dfx,dfy,dfz,vdrift
  real(pre)::delvx_new,delvy_new,delvz_new,ggam,alpha,exp_talpha,ts,tdrag
  type(units)::scl

  integer::iter

  call get_units(scl)

  vgx=zero;vgy=zero;vgz=zero
  dg=zero;tg=zero;pg=zero
  momx=zero
  momy=zero
  momz=zero
  do iter=1,8
   !print *,w(iter)
   dg =dg +cons(1,ig(iter))*w(iter)
   vgx=vgx+cons(2,ig(iter))*w(iter)
   vgy=vgy+cons(3,ig(iter))*w(iter)
   vgz=vgz+cons(4,ig(iter))*w(iter)
   tg =tg +muc_array(ig(iter))*w(iter)
   pg =pg +p(ig(iter))*w(iter)
  enddo
  momx=vgx+vx*d
  momy=vgy+vy*d
  momz=vgz+vz*d
  vgx=vgx/dg
  vgy=vgy/dg
  vgz=vgz/dg
  tg =pg/(scl%rgas*dg)*tg
  delvx=vgx-vx
  delvy=vgy-vy
  delvz=vgz-vz
  call get_gamma_from_tk(dumeps,dg,tg,dumm,ggam)
  cs=sqrt(ggam*pg/dg)
  !print *, "TG =", tg,pg,dg,dumeps,dumm,ggam
  kn=half*dumm*mp/scl%mass/(dg*pi*4.43d-43*asize) !4.44d-43->1d-16cm^2

  

  if(tg<tg_immediate_couple.or.asize>=a_sublimate_limit)then

!    ts=one/(dg*cs/(rhoa*asize)*sqrt(eight/pi))
!    exp_talpha=exp(-two)
!    tdrag=min(t,ts*exp_talpha)
 
    call update_dv(delvx_new,delvx,kn, &
                   cs,dg,rhoa,asize,t,alpha)
    call update_dv(delvy_new,delvy,kn, &
                   cs,dg,rhoa,asize,t,alpha)
    call update_dv(delvz_new,delvz,kn, &
                   cs,dg,rhoa,asize,t,alpha)

    ts=one/alpha
    exp_talpha=exp(-two)
    if(t>ts*exp_talpha)then
!       exp_talpha=exp(-(t-ts*exp_talpha)*alpha)
!       delvx_new=delvx_new*exp_talpha-(one-exp_talpha)*calc_drift(fx,ts)
!       delvy_new=delvy_new*exp_talpha-(one-exp_talpha)*calc_drift(fy,ts)
!       delvz_new=delvz_new*exp_talpha-(one-exp_talpha)*calc_drift(fz,ts)

       weight=(ts*exp_talpha/t)**2

       delvx_new=delvx_new*weight-calc_drift(fx,ts)*(one-weight)
       delvy_new=delvy_new*weight-calc_drift(fy,ts)*(one-weight)
       delvz_new=delvz_new*weight-calc_drift(fz,ts)*(one-weight)
       !print *, "DVEL",time,delvz_new,-calc_drift(fz,ts)
       !if (delvz_new>-calc_drift(fz,ts))delvz_new=-calc_drift(fz,ts)
       
 
    endif

!      else
!       delvx_new=-calc_drift(fx,ts)
!       delvy_new=-calc_drift(fy,ts)
!       delvz_new=-calc_drift(fz,ts)
!      endif
      
!    endif

  else
    delvx_new=zero
    delvy_new=zero
    delvz_new=zero
  endif

  !print *, vx,vy,vz
  !print *, vgx,vgy,vgz
  !print *, "GAS ",delvx,delvy,delvz,d,dg
  !print *, "PAR ",delvx_new,delvy_new,delvz_new

  vxn=(momx+delvx_new*d)/(d+dg)-delvx_new
  vyn=(momy+delvy_new*d)/(d+dg)-delvy_new
  vzn=(momz+delvz_new*d)/(d+dg)-delvz_new

  !print *, "OLD P Velocity ",vx,vy,vz
  !print *, "New Velocity ",vxn,vyn,vzn,vxn+delvx_new,vyn+delvy_new,vzn+delvz_new

  dfx=(vxn-vx)/t
  dfy=(vyn-vy)/t
  dfz=(vzn-vz)/t

 end subroutine

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
  dv_eps=two*bb*dv*beta*exp(-alpha*dt) &
        /(bb*bb-dv*dv*exp(-two*alpha*dt))
  !print *, "TSTOP ",one/alpha,kn,Re,mach,c

  !print *, mach,Re,beta2,beta,alpha,bb,dv_eps,dt
!!!! epstein done.  now stokes

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

  !print *, beta,dv_sto,kn

   dvnew=dv_eps*(three*kn)**2/((three*kn)**2+one)+dv_sto/((three*kn)**2+one)
   !print *,"TIMING ts,t,ts/t",one/alpha,dt,dt*alpha

 end subroutine

end module


