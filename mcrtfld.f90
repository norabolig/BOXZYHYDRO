module mcrtfld
 use parameters
 use derived_types
 use grid_commons
 use utils
 use eos, only: eng_table
 implicit none

 integer::nphottot=6
 integer::RANSEED(1)=314159
 integer::RANLENGTH=1000000

 integer::currentidx
 integer::iter_cool,maxiter=10

 real(pre)::taulimit=half
 real(pre)::taufac=one
 real(pre)::cfrac=0.25d0
 real(pre)::opac_scale=1.00d0
 real(pre),dimension(:),allocatable::divflux,ran_colat,ran_azimuth
 real(pre)::scale_kappa,sigmaSBcode,coolrate,cooltime,totraden,tcoolold
 real(pre)::rho_divflux_limit=1d-10,rho_timestep=1d-10,tau_stream_limit=1d-12
 real(pre)::lumlimit=1d-20
 real(pre)::r_follow_limit=10.0d0,r_follow_limit2
 real(pre)::maxT,ds=zero
 logical:: central_print=.true.
 
 type(units)::scl

 contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 integer function nextidx(i)
  integer,intent(in)::i
  nextidx=i+1
  if(nextidx>RANLENGTH)nextidx=1
 end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 real(pre) function get_kappa(T)
   real(pre),intent(in)::T
!   get_kappa=(T/16d0)**.5d0*scale_kappa
!   get_kappa=1d-1*scale_kappa
!   return
!     Pollack et al. (1994) rosseland opacities

      if(T.lt.80.0) then
         get_kappa=(T**2)/3200.0
      else if(T.lt.170.0) then
         get_kappa=-2.0 + 0.050*T
      else if(T.lt.180.0) then
         get_kappa=62.60 - 0.330*T
      else if(T.lt.270.0) then
         get_kappa=-1.0 + 0.0233*T
      else if(T.lt.300.0) then
         get_kappa=8.0 - 0.010*T
      else if(T.lt.425.0) then
         get_kappa=1.88 + 0.0104*T
      else if(T.lt.440.0) then
         get_kappa=128.13 - 0.2867*T
      else if(T.lt.670.0) then
         get_kappa=0.57 + 0.0033*T
      else if(T.lt.700.0) then
         get_kappa=19.50 - 0.0250*T
      else if(T.lt.1300.0) then
         get_kappa=-0.33 + 0.0033*T
      else if(T.lt.1350.0) then
         get_kappa=24.80 - 0.0160*T
      else if(T.lt.1449.66) then
         get_kappa=46.40 - 0.0320*T
      else
         get_kappa=0.01
      endif
      get_kappa=get_kappa*opac_scale*scale_kappa

 end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 real(pre) function get_luminosity(rho,kappa,T)
  real(pre),intent(in)::rho,kappa,T
  get_luminosity=four*sigmaSBcode*rho*kappa*T**4/dble(nphottot)
 end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine initialize_mcrtfld()
  integer::i
  
  allocate(divflux(ngrid))
  allocate(ran_azimuth(RANLENGTH))
  allocate(ran_colat(RANLENGTH))
  call get_units(scl)
  scale_kappa=scl%mass/scl%length**2
  sigmaSBcode=5.67d-5/(scl%mass/scl%time**3)
  ds=sqrt(dx*dx+dy*dy+dz*dz)
  !call random_seed(put=RANSEED)
  !call ramdom_seed()
  call random_number(ran_colat)
  call random_number(ran_azimuth)
  do i=1,RANLENGTH
    ran_azimuth(i)=ran_azimuth(i)*two*pi
    ran_colat(i)=acos(one-two*ran_colat(i))
  enddo
  currentidx=1
  r_follow_limit2=r_follow_limit**2

 end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 real(pre) function fluxlimiter(rk,T,dtds_mag)
   real(pre)::rk,T,dtds_mag,y
   y=four/(rk*T)*dtds_mag
   fluxlimiter=(two+y)/(six+y*(three+y))
 end function

 subroutine fld(igrid,boundary_e)
  integer,intent(in)::igrid
  real(pre)::Ta,kappa0,rho0,dtds,dtds_mag,rgas_loc
  real(pre)::T0,T1,kappa1,rho1,rk,beta,vol
  real(pre)::fzt,fzb,fxl,fxr,fyt,fyb
  real(pre),dimension(6)::boundary_e
  integer::b(6)
  logical::active

  rgas_loc=scl%rgas 

  vol=dx*dy*dz

   call get_boundary(igrid,b)

   rho0=cons(1,igrid)

   T0=p(igrid)/(rgas_loc*rho0)*muc_array(igrid)
   kappa0=get_kappa(T0)

   rho1=cons(1,b(1))
   T1=p(b(1))/(rgas_loc*rho1)*muc_array(b(1))
   kappa1=get_kappa(T1)

   rk=half*(rho0*kappa0+rho1*kappa1)
   Ta  =half*(T0+t1)

   dtds=(T1-T0)/dy
   dtds_mag=abs(dtds)
   
   beta=fluxlimiter(rk,ta,dtds_mag)

   fyt=16d0*sigmaSBcode*beta*Ta**3*dtds/rk

!in2

   rho1=cons(1,b(2))
   T1=p(b(2))/(rgas_loc*rho1)*muc_array(b(2))
   kappa1=get_kappa(T1)

   rk=half*(rho0*kappa0+rho1*kappa1)
   Ta  =half*(T0+t1)

   dtds=(T1-T0)/dy
   dtds_mag=abs(dtds)
   
   beta=fluxlimiter(rk,ta,dtds_mag)

   fyb=16d0*sigmaSBcode*beta*Ta**3*dtds/rk


!b(3)

   rho1=cons(1,b(3))
   T1=p(b(3))/(rgas_loc*rho1)*muc_array(b(3))
   kappa1=get_kappa(T1)

   rk=half*(rho0*kappa0+rho1*kappa1)
   Ta  =half*(T0+t1)

   dtds=(T1-T0)/dx
   dtds_mag=abs(dtds)
   
   beta=fluxlimiter(rk,ta,dtds_mag)

   fxr=16d0*sigmaSBcode*beta*Ta**3*dtds/rk

!b(4)

   rho1=cons(1,b(4))
   T1=p(b(4))/(rgas_loc*rho1)*muc_array(b(4))
   kappa1=get_kappa(T1)

   rk=half*(rho0*kappa0+rho1*kappa1)
   Ta  =half*(T0+t1)

   dtds=(T1-T0)/dx
   dtds_mag=abs(dtds)
   
   beta=fluxlimiter(rk,ta,dtds_mag)

   fxl=16d0*sigmaSBcode*beta*Ta**3*dtds/rk

!b(5)

   rho1=cons(1,b(5))
   T1=p(b(5))/(rgas_loc*rho1)*muc_array(b(5))
   kappa1=get_kappa(T1)

   rk=half*(rho0*kappa0+rho1*kappa1)
   Ta  =half*(T0+t1)

   dtds=(T1-T0)/dz
   dtds_mag=abs(dtds)
   
   beta=fluxlimiter(rk,ta,dtds_mag)

   fzt=16d0*sigmaSBcode*beta*Ta**3*dtds/rk

!b(6)

   rho1=cons(1,b(6))
   T1=p(b(6))/(rgas_loc*rho1)*muc_array(b(6))
   kappa1=get_kappa(T1)

   rk=half*(rho0*kappa0+rho1*kappa1)
   Ta  =half*(T0+t1)

   dtds=(T1-T0)/dz
   dtds_mag=abs(dtds)
   
   beta=fluxlimiter(rk,ta,dtds_mag)

   fzb=16d0*sigmaSBcode*beta*Ta**3*dtds/rk

!!$OMP ATOMIC
!     divflux(igrid)=divflux(igrid)+( area_side*((ftt+fbb+ftr+fbl+ftl+fbr))/(vol) &
!                     + (area_top*(ftz+fbz)/vol))

   boundary_e(1)=dz*dx*max(-fyt,zero)/vol
   boundary_e(2)=dz*dx*max(-fyb,zero)/vol
   boundary_e(3)=dz*dy*max(-fxr,zero)/vol
   boundary_e(4)=dz*dy*max(-fxl,zero)/vol
   boundary_e(5)=dx*dy*max(-fzt,zero)/vol
   boundary_e(6)=dx*dy*max(-fzb,zero)/vol
!   print *, boundary_e
!   print *, rho0,t0,rho1,t1,p(igrid),p(b(1)),p(b(2)),p(b(3)),p(b(4)),p(b(5)),p(b(6)),p(in7),p(in8)
!   print *, cons(1,igrid),cons(1,b(1)),cons(1,b(2)),cons(1,b(3)),cons(1,b(4)),cons(1,b(5)),cons(1,b(6)), &
! cons(1,in7),cons(1,in8)
!   print *, muc,rgas_loc
!   stop


 end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine transfer_radiation()
  integer::idx,igrid,iphot,nphot
  real(pre)::lum,lum0,x,y,z,x0,y0,z0,azimuth,colat,azimuthp,colatp,dcol,daz,T,kappa,vol
  real(pre)::area_top,area_side,area_tot,dtau,dl,dtaucell,r
  real(pre),dimension(6)::boundary_e

  vol=dx*dy*dz
  dl=ds*taufac

  idx=currentidx
  maxT=zero

!$OMP DO SCHEDULE(STATIC)
  do igrid=1,ngrid
   divflux(igrid)=zero
  enddo
!$OMP ENDDO
!$OMP MASTER
  central_print=.true.
!$OMP END MASTER
!$OMP BARRIER
!$OMP DO SCHEDULE(DYNAMIC) PRIVATE(T,azimuthp,colatp,lum,lum0,x,y,z,x0,y0,z0) &
!$OMP&PRIVATE(azimuth,colat,dcol,daz,kappa,area_top,area_side,area_tot,boundary_e) &
!$OMP&PRIVATE(dtau,dtaucell,r) &
!$OMP&REDUCTION(max:maxT)
  do igrid=1,ngrid

   if(grid(igrid)%boundary>0)cycle
   x0=grid(igrid)%x;y0=grid(igrid)%y;z0=grid(igrid)%z

   T=p(igrid)/(scl%rgas*cons(1,igrid))*muc_array(igrid)
   if(T<zero)then
    print *, "WTF???",T,p(igrid)
   endif
   maxT=max(maxT,T)
   kappa=get_kappa(T)
   dtau=kappa*cons(1,igrid)*dl
   dtaucell=dtau*ds/dl

   if(x0**2+y0**2+z0**2<dx**2)then
!!!!!$OMP ATOMIC
!!!!!      divflux(igrid)=divflux(igrid)+1d-4*3.8e33/(scl%mass*scl%length**2/scl%time**3)/(eight*dx*dy*dz)
      if(central_print)print *, "Temperature Central ",time,T,dtau,cons(1,igrid),time,muc_array(igrid),p(igrid)
      central_print=.false.
   endif


   if(dtaucell>taulimit)then
     call fld(igrid,boundary_e)
     lum0=zero
     do iphot=1,6
       lum0=lum0+boundary_e(iphot)
     enddo
     nphot=6
     lum0=lum0/dble(nphot)
   else
     lum0=get_luminosity(cons(1,igrid),kappa,T)
     nphot=nphottot
   endif

!$OMP ATOMIC
   divflux(igrid)=divflux(igrid)-lum0*dble(nphot)



   if(dtaucell>tau_stream_limit.and.T>tk_bgrnd*1.1.and.cons(1,igrid)>rho_divflux_limit)then
   do iphot=1,nphot

!     azimuth=ran_azimuth(idx)
!     colat=ran_colat(idx)
!     idx=nextidx(idx) 

     x=x0;y=y0;z=z0
     lum=lum0

     if(dtaucell>taulimit)then
       colatp=zero
       azimuthp=zero
       daz=zero
       dcol=zero
       select case(iphot)

       case(1)
         colatp=half*pi
         azimuthp=half*pi
         lum=boundary_e(1)
       case(2)
         colatp=half*pi
         azimuthp=1.5*pi
         lum=boundary_e(2)
       case(3)
         colatp=pi*half
         azimuthp=zero
         lum=boundary_e(3)
       case(4)
         colatp=pi*half
         azimuthp=pi
         lum=boundary_e(4)
       case(5)
         colatp=zero
         azimuthp=zero
         lum=boundary_e(5)
       case(6)
         colatp=pi
         azimuthp=zero
         lum=boundary_e(6)
       case(7)
         print *, "Only 6 photons allowed for now."
         stop
       end select
  
       x=x+ds*1.01d0*sin(colatp)*cos(azimuthp)*half
       y=y+ds*1.01d0*sin(colatp)*sin(azimuthp)*half
       z=z+dz*1.01d0*cos(colatp)*half

       colat=colatp
       azimuth=azimuthp
     else
      azimuth=ran_azimuth(idx)
      colat=ran_colat(idx)
      idx=nextidx(idx) 
     endif

     call propogate_photons(colat,azimuth,lum,x,y,z)
   enddo
   endif

  enddo
!$OMP ENDDO 
!$OMP MASTER
  currentidx=idx
!$OMP END MASTER

 end subroutine
     
 subroutine propogate_photons(colat,azimuth,lum0,x0,y0,z0)
     integer::igrid,igrid0
     real(pre),intent(in)::colat,azimuth,x0,y0,z0,lum0
     real(pre)::lum,x,y,z,dtau,kappa,T,absorb,dl,tau
     dl=ds*taufac

     x=x0;y=y0;z=z0
     igrid=get_grid_indx(x,y,z)
     T=p(igrid)/(scl%rgas*cons(1,igrid))*muc_array(igrid)
     igrid0=igrid
     lum=lum0
     tau=zero
     do while(lum>lumlimit)
       igrid=get_grid_indx(x,y,z)
       if(grid(igrid)%boundary>0)exit
       if(igrid>ngrid)then
         print *, "Major problem.  igrid>ngrid",igrid,ngrid
         stop
       endif
       T=p(igrid)/(scl%rgas*cons(1,igrid))*muc_array(igrid)
       kappa=get_kappa(T)
       dtau=cons(1,igrid)*dl*kappa
       !tau=tau+dtau
       absorb=(one-exp(-dtau))*lum

!$OMP ATOMIC
       divflux(igrid)=divflux(igrid)+absorb

       lum=lum-absorb
       x=x+dl*sin(colat)*cos(azimuth)
       y=y+dl*sin(colat)*sin(azimuth)
       z=z+dl*cos(colat)
       if(x*x+y*y+z*z>r_follow_limit2)lum=zero
     enddo
     !if(t0>0..and.t0<450.)then
!     if(tau<10.)then
!       print *,x0,y0,z0,x,y,z,lum,lum0,dtau,tau,igrid,igrid0,grid(igrid0)%x,grid(igrid0)%y,grid(igrid0)%z,T,T0, &
!       cons(1,igrid0),p(igrid0),muc_array(igrid0)
!     endif
 end subroutine

 subroutine heat_and_cool(iterate)
   logical::iterate
   integer::igrid
   real(pre)::tcool,eps,ekin,vol
  

  vol=dz*dx*dy

   call transfer_radiation()
!$OMP MASTER
  coolrate=zero
  totraden=zero
  iter_cool=iter_cool+1
!$OMP END MASTER
!$OMP BARRIER
!$OMP DO SCHEDULE(STATIC) REDUCTION(max:coolrate) PRIVATE(ekin,eps)
  do igrid=1,ngrid
   if(grid(igrid)%boundary>0)cycle
   if(cons(1,igrid)<rho_timestep)cycle
   
   ekin=half*(cons(2,igrid)**2+cons(3,igrid)**2+cons(4,igrid)**2)/cons(1,igrid)
   eps=max(cons(5,igrid)-ekin,small_eps)
   if(eps<=eng_table(1,1)*cons(1,igrid))cycle
   !if(eps<=small_eps)cycle
   coolrate=max(coolrate,abs(divflux(igrid))/(cfrac*eps))
  enddo
!$OMP ENDDO
!$OMP BARRIER
  tcool=one/coolrate
  !if(tcoolold>zero)tcool=max(tcool,1.1d0*tcoolold)
  if(tcool+cooltime>dt.or.iter_cool==maxiter)tcool=dt-cooltime
!$OMP MASTER
  tcoolold=tcool
  cooltime=cooltime+tcool
  print *, "maxT ",time,maxT
  print *, "tcool,cooltime,dt",tcool,cooltime,dt
!$OMP END MASTER
!$OMP BARRIER
!$OMP DO SCHEDULE(STATIC) PRIVATE(ekin,eps) REDUCTION(+:totraden)
  do igrid=1,ngrid
   if(grid(igrid)%boundary>0)cycle
   ekin=half*(cons(2,igrid)**2+cons(3,igrid)**2+cons(4,igrid)**2)/cons(1,igrid)
   eps=max(cons(5,igrid)-ekin+divflux(igrid)*tcool,small_eps)
   totraden=totraden+divflux(igrid)*vol
   cons(5,igrid)=eps+ekin
  enddo
!$OMP ENDDO
!$OMP BARRIER
  if(cooltime<dt)then 
     iterate= .true.
  else
     iterate= .false.
  endif
  call state()
!$OMP MASTER
   print *,"#Total Luminosity", time,totraden*scl%mass*scl%length**2/scl%time**3/3.8d33
!$OMP END MASTER
!stop

 end subroutine

 subroutine mcrtfld_transfer
   logical::iterate

!$OMP MASTER
   cooltime=zero
   tcoolold=zero
   iter_cool=0
!$OMP END MASTER
!$OMP BARRIER
   iterate=.true.
   do while(iterate)
      call  heat_and_cool(iterate)
   enddo
  
 end subroutine

end module



