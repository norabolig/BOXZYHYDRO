!
! Module for radiative transfer. This uses a combination of flux limited diffusion
! and a particle/ray-tracing approach.
!
module mcrtfld
 use parameters
 use derived_types
 use grid_commons
 use utils
 use eos, only: eng_table
 implicit none

 integer::nphottot=6
 integer::RANSEED(1)=314159
 integer::RANLENGTH=50000

 integer::currentidx,ran60_id
 integer::iter_cool,maxiter=1000
!
!
#ifdef FLDONLY
 real(pre)::taulimit=0.0d0
#else
 real(pre)::taulimit=1.0d20
#endif /* end ifdef FLDONLY */
!
!
 real(pre)::taufac=1d0
 real(pre)::cfrac=0.1d0
 real(pre)::opac_scale=1d0,taum=1d-1
 real(pre),dimension(:),allocatable::divflux,ran_colat,ran_azimuth,ran_pm45
 real(pre)::scale_kappa,sigmaSBcode,coolrate,cooltime,totraden,tcoolold,eold,total_cool
 real(pre)::rho_divflux_limit=1d-7,rho_timestep=1d-7,tau_stream_limit=1d-100
 real(pre)::lumlimit=-1d0
 real(pre)::r_follow_limit=1d1,r_follow_limit2
 real(pre)::maxT,ds=zero
 logical:: central_print,print_iter=.true.
 
 type(units)::scl

 contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Move to the next id. Used for the Monte Carlo-like portion.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 integer function nextidx(i)
  integer,intent(in)::i
  nextidx=i+1
  if(nextidx>RANLENGTH)nextidx=1
 end function
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interpolation function for the FLD and the free-streaming limits
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 real(pre) function expinterp(r,k,dy)
   real(pre),intent(in)::r,k,dy
      expinterp=exp(-three*(r*k*dy)**2)!*0.1d0)
 end function
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Given an alignment along an x,y,z axis (idir), find a perturbed
! angle for propagating a photon packet.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

 subroutine get_angle_pert(xp,yp,zp,idx,idir)
   real(pre)::az,co,x,y,z,xp,yp,zp,cosaz,sinaz,cosco,sinco
   integer::idx,idir

   az=ran_azimuth(idx);idx=nextidx(idx)
   co=ran_pm45(idx);idx=nextidx(idx)

   cosco=cos(co)
   sinco=sin(co)

   cosaz=cos(az)
   sinaz=sin(az)

   x=sinco*cosaz
   y=sinco*sinaz
   z=cosco

   xp=z;yp=y;zp=z
   
   select case(idir)

   case(1)

    xp=x
    yp=z
    zp=-y

   case(2)

    xp=x
    yp=-z
    zp=y

   case(3)

    xp=z
    yp=y
    zp=-x

   case(4)

    xp=-z
    yp=y
    zp=x
 
   case(5)

    xp=x
    yp=y
    zp=z

   case(6)

    xp=x
    yp=-y
    zp=-z
      
   end select 

   return

 end subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get the opacity.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 real(pre) function get_kappa(T)
   real(pre),intent(in)::T
!
!***
!     Pollack et al. (1994) rosseland opacities
!***
!
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
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get emission from a cell in free-streaming limit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 real(pre) function get_luminosity(rho,kappa,T)
  real(pre),intent(in)::rho,kappa,T
  get_luminosity=four*sigmaSBcode*rho*kappa*T**4/dble(nphottot)
 end function
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Similar to above, but integrated through the cell. This is still in the
! test phase and is not used.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 real(pre) function get_luminosity_2(rho,kappa,T,ds)
  real(pre),intent(in)::rho,T,ds,kappa
  real(pre)::exptau_loc
  exptau_loc=exp(-rho*kappa*ds)
  get_luminosity_2=four*sigmaSBcode*(T**4-tk_bgrnd**4)*(one-exptau_loc) &
      / (dble(nphottot)*ds)
 end function
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize variables and arrays for mcrtfld routine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 subroutine initialize_mcrtfld()
  integer::i
  
  allocate(divflux(ngrid))
  allocate(ran_azimuth(RANLENGTH))
  allocate(ran_colat(RANLENGTH))
  allocate(ran_pm45(RANLENGTH))
  call get_units(scl)
  scale_kappa=scl%mass/scl%length**2
  sigmaSBcode=5.67d-5/(scl%mass/scl%time**3)
  ds=sqrt(dx*dx+dy*dy+dz*dz)
  call random_number(ran_colat)
  call random_number(ran_azimuth)
  do i=1,RANLENGTH
    ran_azimuth(i)=ran_azimuth(i)*two*pi
    ran_pm45(i)=acos(ran_colat(i))
    ran_colat(i)=acos(one-two*ran_colat(i))
  enddo
  currentidx=1
  r_follow_limit2=r_follow_limit**2
  if (lumlimit<zero)lumlimit=four*sigmaSBcode*Tk_bgrnd**4/(dble(nphottot)*(dx+dy+dz)/three)
  print *,"Luminosity limiter is set to ",lumlimit

 end subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Bodenheimer et al. 1990 flux limiter for FLD contribution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 real(pre) function fluxlimiter(rk,T,dtds_mag)
   real(pre)::rk,T,dtds_mag,y
   y=four/(rk*T)*dtds_mag
   fluxlimiter=(two+y)/(six+y*(three+y))
 end function
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Flux-limited diffusion solution.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
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
!
!***
!in2
!***
!
   rho1=cons(1,b(2))
   T1=p(b(2))/(rgas_loc*rho1)*muc_array(b(2))
   kappa1=get_kappa(T1)

   rk=half*(rho0*kappa0+rho1*kappa1)
   Ta  =half*(T0+t1)

   dtds=(T1-T0)/dy
   dtds_mag=abs(dtds)
   
   beta=fluxlimiter(rk,ta,dtds_mag)

   fyb=16d0*sigmaSBcode*beta*Ta**3*dtds/rk
!
!***
!b(3)
!***
!
   rho1=cons(1,b(3))
   T1=p(b(3))/(rgas_loc*rho1)*muc_array(b(3))
   kappa1=get_kappa(T1)

   rk=half*(rho0*kappa0+rho1*kappa1)
   Ta  =half*(T0+t1)

   dtds=(T1-T0)/dx
   dtds_mag=abs(dtds)
   
   beta=fluxlimiter(rk,ta,dtds_mag)

   fxr=16d0*sigmaSBcode*beta*Ta**3*dtds/rk
!
!***
!b(4)
!***
!
   rho1=cons(1,b(4))
   T1=p(b(4))/(rgas_loc*rho1)*muc_array(b(4))
   kappa1=get_kappa(T1)

   rk=half*(rho0*kappa0+rho1*kappa1)
   Ta  =half*(T0+t1)

   dtds=(T1-T0)/dx
   dtds_mag=abs(dtds)
   
   beta=fluxlimiter(rk,ta,dtds_mag)

   fxl=16d0*sigmaSBcode*beta*Ta**3*dtds/rk
!
!***
!b(5)
!***
!
   rho1=cons(1,b(5))
   T1=p(b(5))/(rgas_loc*rho1)*muc_array(b(5))
   kappa1=get_kappa(T1)

   rk=half*(rho0*kappa0+rho1*kappa1)
   Ta  =half*(T0+t1)

   dtds=(T1-T0)/dz
   dtds_mag=abs(dtds)
   
   beta=fluxlimiter(rk,ta,dtds_mag)

   fzt=16d0*sigmaSBcode*beta*Ta**3*dtds/rk
!
!***
!b(6)
!***
!
   rho1=cons(1,b(6))
   T1=p(b(6))/(rgas_loc*rho1)*muc_array(b(6))
   kappa1=get_kappa(T1)

   rk=half*(rho0*kappa0+rho1*kappa1)
   Ta  =half*(T0+t1)

   dtds=(T1-T0)/dz
   dtds_mag=abs(dtds)
   
   beta=fluxlimiter(rk,ta,dtds_mag)

   fzb=16d0*sigmaSBcode*beta*Ta**3*dtds/rk
!
!
#ifdef FLDONLY
!$OMP ATOMIC
     divflux(igrid)=divflux(igrid)+( dx*dz*(fyt+fyb)+dy*dz*(fxl+fxr)+dx*dy*(fzt+fzb))/(vol) 
     boundary_e=zero
#else
   boundary_e(1)=dz*dx*max(-fyt,zero)/vol
   boundary_e(2)=dz*dx*max(-fyb,zero)/vol
   boundary_e(3)=dz*dy*max(-fxr,zero)/vol
   boundary_e(4)=dz*dy*max(-fxl,zero)/vol
   boundary_e(5)=dx*dy*max(-fzt,zero)/vol
   boundary_e(6)=dx*dy*max(-fzb,zero)/vol
#endif /* end ifdef FLDONLY */
!
!
 end subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Main work function for photon propagation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 subroutine transfer_radiation()
  integer::idx,igrid,iphot,nphot,id60
  real(pre)::lum,lum0,x,y,z,x0,y0,z0,azimuth,colat,azimuthp,colatp,dcol,daz,T,kappa,vol
  real(pre)::area_top,area_side,area_tot,dtau,dl,dtaucell,r,atten,exptau
  real(pre)::xp,yp,zp,delx,dely,delz
  real(pre),dimension(6)::boundary_e

  vol=dx*dy*dz
  dl=ds*taufac

  idx=currentidx
  id60=currentidx
  maxT=zero

!$OMP DO SCHEDULE(STATIC) PRIVATE(kappa,dtau,T,x0,y0,z0,r)
  do igrid=1,ngrid
!
!***
! Below is a bit of a mess.  The commented code is used for
! dissipation, which is used for testing. Leaving it in for now.
!***
!
!   r=sqrt(x0*x0+y0*y0+z0*z0)
!   if(r>1.4*dz)then
      divflux(igrid)=zero
!   else
      !T=p(igrid)/(scl%rgas*cons(1,igrid))*muc_array(igrid)
      !if(T<zero)then
      !  print *, "WTF???",T,p(igrid)
      !endif
      !kappa=get_kappa(T)
      !dtau=kappa*cons(1,igrid)*ds
!      divflux(igrid)=1d-4*3.8e33/(scl%mass*scl%length**2/scl%time**3)/(eight*dx*dy*dz)
      !divflux(igrid)=0.14*dtau/(scl%mass/scl%time**3)/(taum*two*dz)
      !divflux(igrid)=1d-4*3.8e33/(scl%mass*scl%length**2/scl%time**3)/(4d0*pi*16d0)*(4d0/r)
!   endif
  enddo
!$OMP ENDDO
!!!$OMP MASTER
!!  central_print=.true.
!!!$OMP END MASTER
!$OMP BARRIER
!
!
!$OMP DO SCHEDULE(DYNAMIC) PRIVATE(T,azimuthp,colatp,lum,lum0,x,y,z,x0,y0,z0) &
!$OMP&PRIVATE(azimuth,colat,dcol,daz,kappa,area_top,area_side,area_tot,boundary_e) &
!$OMP&PRIVATE(dtau,dtaucell,r,atten,exptau,delx,dely,delz,xp,yp,zp) &
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
   dtaucell=max(dtau*ds/dl,1d-20)

!
!
#ifdef VERBOSE
   if(x0**2+y0**2+z0**2<two*dx**2)then
      if(central_print)print *, "Temperature Central ",time,T,dtau,cons(1,igrid),time,muc_array(igrid),p(igrid)
      central_print=.false.
   endif
#endif /* VERBOSE */
!
!
     call fld(igrid,boundary_e)
     nphot=nphottot

   do iphot=1,nphot

     x=x0;y=y0;z=z0
     lum=zero
     exptau=zero
     delx=zero
     dely=zero
     delz=zero

       select case(iphot)
       case(1)
         delz=zero
         delx=zero
         dely= half*dy*1.01d0
         lum=get_luminosity_2(cons(1,igrid),kappa,T,dy)
       exptau=expinterp(cons(1,igrid),kappa,dy)
       case(2)
         delz=zero
         delx=zero
         dely=-half*dy*1.01d0
         lum=get_luminosity_2(cons(1,igrid),kappa,T,dy)
       exptau=expinterp(cons(1,igrid),kappa,dy)
       case(3)
         delz=zero
         delx= half*dx*1.01d0
         dely=zero
         lum=get_luminosity_2(cons(1,igrid),kappa,T,dx)
       exptau=expinterp(cons(1,igrid),kappa,dx)
       case(4)
         delz=zero
         delx=-half*dx*1.01d0
         dely=zero
         lum=get_luminosity_2(cons(1,igrid),kappa,T,dx)
       exptau=expinterp(cons(1,igrid),kappa,dx)
       case(5)
         delz= half*dz*1.01d0
         delx=zero
         dely=zero
         lum=get_luminosity_2(cons(1,igrid),kappa,T,dz)
       exptau=expinterp(cons(1,igrid),kappa,dz)
       case(6)
         delz=-half*dz*1.01d0
         delx=zero
         dely=zero
         lum=get_luminosity_2(cons(1,igrid),kappa,T,dz)
       exptau=expinterp(cons(1,igrid),kappa,dz)
 
       case(7)
         print *, "Only 6 photons allowed for now."
         stop
       end select
!
       lum=lum*exptau+boundary_e(iphot)*(one-exptau)
!
!
#ifdef FLDONLY
       lum=zero
#endif
!
!
!$OMP ATOMIC
      divflux(igrid)=divflux(igrid)-lum

   if(cons(1,igrid)>rho_divflux_limit)then
       x=x0+delx
       y=y0+dely
       z=z0+delz
       call get_angle_pert(xp,yp,zp,idx,iphot)
       call propogate_photons(xp,yp,zp,lum,x,y,z)
   endif
   enddo
  enddo
!$OMP ENDDO 
!
!
!$OMP MASTER
  currentidx=idx
!$OMP END MASTER
!
!
 end subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
! Follow the photons along a given direction.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 subroutine propogate_photons(xp,yp,zp,lum0,x0,y0,z0)
     integer::igrid,igrid0
     real(pre),intent(in)::xp,yp,zp,x0,y0,z0,lum0
     real(pre)::lum,x,y,z,dtau,kappa,T,absorb,dl,tau,delx,dely,delz
     dl=min(min(dx,dy),dz)*half
     delx=dl*xp
     dely=dl*yp
     delz=dl*zp

     x=x0;y=y0;z=z0
     igrid=get_grid_indx(x,y,z)
     T=p(igrid)/(scl%rgas*cons(1,igrid))*muc_array(igrid)
     igrid0=igrid
     lum=lum0
     tau=zero

     do while(lum>lumlimit)
       igrid=get_grid_indx(x,y,z)
       if(grid(igrid)%boundary>0)exit
       if(igrid>ngrid.or.igrid<1)then
         print *, "Major problem.  igrid>ngrid",igrid,ngrid
         stop
       endif
       T=p(igrid)/(scl%rgas*cons(1,igrid))*muc_array(igrid)
       kappa=get_kappa(T)
       dtau=cons(1,igrid)*dl*kappa
!
!
#ifdef FLDONLY
       if(cons(1,igrid)<rho_divflux_limit)lum=zero
       absorb=lum
#else
       absorb=(one-exp(-dtau))*lum
#endif /* end ifdef FLDONLY */
!
!
!$OMP ATOMIC
       divflux(igrid)=divflux(igrid)+absorb

       lum=lum-absorb
       x=x+delx
       y=y+dely
       z=z+delz
       if(x*x+y*y+z*z>r_follow_limit2)lum=zero
     enddo
 end subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Apply divergence of the flux and determine whether a subcycle is 
! necessary.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 subroutine heat_and_cool(iterate)
   logical::iterate
   integer::igrid
   real(pre)::tcool,eps,ekin,vol,rate_temp,T

   vol=dz*dx*dy

   call transfer_radiation()
!
!$OMP MASTER
  coolrate=zero
  totraden=zero
  iter_cool=iter_cool+1
!$OMP END MASTER
!$OMP BARRIER
!
!
!$OMP DO SCHEDULE(STATIC) REDUCTION(max:coolrate) PRIVATE(ekin,eps)
  do igrid=1,ngrid
   if(grid(igrid)%boundary>0)cycle
   if(cons(1,igrid)<rho_timestep)cycle
   
   ekin=half*(cons(2,igrid)**2+cons(3,igrid)**2+cons(4,igrid)**2)/cons(1,igrid)
   eps=max(cons(5,igrid)-ekin,small_eps)
   if(eps<=eng_table(1,1)*cons(1,igrid))cycle
   coolrate=max(coolrate,abs(divflux(igrid))/(cfrac*eps))
  enddo
!$OMP ENDDO
!
!
!$OMP BARRIER
  tcool=one/coolrate
  if(tcool+cooltime>dt.or.iter_cool==maxiter)tcool=dt-cooltime
!
!$OMP MASTER
  tcoolold=tcool
  cooltime=cooltime+tcool
  if(cooltime>=dt)print_iter=.true.
  if(print_iter)then
   print *, iter_cool,"maxT ",time,maxT, "tcool,cooltime,dt",tcool,cooltime,dt
  endif
!$OMP END MASTER
!
!
!$OMP BARRIER
!
!
!$OMP DO SCHEDULE(STATIC) PRIVATE(ekin,eps,eold) REDUCTION(+:totraden)
  do igrid=1,ngrid
   if(grid(igrid)%boundary>0)cycle
   ekin=half*(cons(2,igrid)**2+cons(3,igrid)**2+cons(4,igrid)**2)/cons(1,igrid)
   eold=max(cons(5,igrid)-ekin,small_eps)
   eps=max(eold+divflux(igrid)*tcool,small_eps)
!
!
!#ifdef FLDONLY
!   if((eps-eold)/eold>cfrac)then
!       divflux(igrid)=eold*cfrac/tcool
!       eps=eold*(one+cfrac)
!   elseif( -(eps-eold)/eold>cfrac) then
!       divflux(igrid)=-eold*cfrac/tcool
!       eps=eold*(one-cfrac)
!   endif
!#endif /* end ifdef FLDONLY */
!
!
   totraden=totraden+divflux(igrid)*vol
   cons(5,igrid)=eps+ekin
  enddo
!$OMP ENDDO
!$OMP BARRIER
!
!
  if(cooltime<dt)then 
     iterate= .true.
  else
     iterate= .false.
  endif
  call state() ! state.f90
!$OMP MASTER
!
!
  total_cool=total_cool+totraden*tcool
  if(print_iter)then
   print *,"#Total Luminosity", time,totraden*scl%mass*scl%length**2/scl%time**3/3.8d33,totraden*dt
   print_iter=.false.
  endif
!$OMP END MASTER
 end subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Main driving routine.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 subroutine mcrtfld_transfer()
   logical::iterate

!$OMP MASTER
   cooltime=zero
   tcoolold=zero
   total_cool=zero
   iter_cool=0
   print_iter=.true.
   central_print=.true.
!$OMP END MASTER
!$OMP BARRIER
   iterate=.true.
   do while(iterate)
      call  heat_and_cool(iterate)
   enddo
 end subroutine

end module
