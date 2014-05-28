!
! Module for adding particles to the hydrodynamics calculations
! Particles can be integrated with or without drag.
! All drag calculations include feedback with the gas.
! Particles can be used as particle-in-cell, and there can also
! be particle-particle interactions using direct sum.
!
module particle
!$ use omp_lib
 use parameters
 use derived_types
 use grid_commons
 use utils
!
!
#ifdef WITHDRAG
 use pdrag
#endif
!
!
 implicit none

 type(particle_type),dimension(:),allocatable::part,part_direct
 real(pre)::prate,ds
 integer::nact1,nact2
 type(units)::scl

 contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Weights for cloud in particle-in-cell routine.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 subroutine get_weights_8(ig,x,y,z,w)
  real(pre),intent(in)::x,y,z
  real(pre),intent(out)::w(8)
  real(pre)::delx1,delx2,delx3,delx4,dely1,dely2,dely3,dely4,area1,area2,area3,area4
  real(pre)::weight
  integer,intent(out)::ig(8)
  integer::iter,flag,id
    flag=0
    call get_nearest_8(x,y,z,ig) ! 1-3 are low
    do iter=1,8
      if(ig(iter)<1.or.ig(iter)>ngrid)flag=1
      if(grid(ig(iter))%boundary==1.or.grid(ig(iter))%boundary==3)flag=1
!      print *, "BOUNDARY VALUE = ",grid(ig(iter))%boundary," for neighbor ",iter
!      print *, x/dx,y/dy,z/dz,flag,ig(iter),ngrid
    enddo
    if(flag==1)then
      ig=-1
      w=zero
      return
    endif
    delx1=abs(grid(ig(1))%x-x)
    dely1=abs(grid(ig(1))%y-y)
    delx2=abs(grid(ig(2))%x-x)
    dely2=abs(grid(ig(2))%y-y)
    delx3=abs(grid(ig(3))%x-x)
    dely3=abs(grid(ig(3))%y-y)
    delx4=abs(grid(ig(4))%x-x)
    dely4=abs(grid(ig(4))%y-y)

    area3=delx1*dely1
    area4=delx2*dely2
    area1=delx3*dely3
    area2=delx4*dely4

    if (area1<0.0.or.area2<0.0.or.area3<0.0.or.area4<0.0)then
       print *, "AREA PROBLEM ",area1,area2,area3,area4
    endif
 
    if(.not.z==zero)then
      w(1)=area1*(grid(ig(5))%z-z)
      w(2)=area2*(grid(ig(5))%z-z)
      w(3)=area3*(grid(ig(5))%z-z)
      w(4)=area4*(grid(ig(5))%z-z)
      w(5)=area1*(z-grid(ig(1))%z)
      w(6)=area2*(z-grid(ig(1))%z)
      w(7)=area3*(z-grid(ig(1))%z)
      w(8)=area4*(z-grid(ig(1))%z)
    else
      w(1)=area1
      w(2)=area2
      w(3)=area3
      w(4)=area4
      w(5)=area1
      w(6)=area2
      w(7)=area3
      w(8)=area4
    endif

    weight=zero
    flag=0
    do iter=1,8
     weight=weight+w(iter)
     if(w(iter)<0)then
        print *, "WEIGHT I ",iter,w(iter)
        print *, grid(ig(5))%z, grid(ig(1))%z,z,z/dz+nz/2
     endif
     if(w(iter)<0)flag=1
    enddo
    if (weight<=zero)flag=1
    if(flag==1)then
      ig=-1
      w=zero
      return
    endif
 
    do iter=1,8
     w(iter)=w(iter)/weight
    enddo
 
 end subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Intialize variables and arrays for particle integration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 subroutine initialize_particles()
  integer::ipart,id,iseed=314,istat,step,isize
  real(pre)::dphi,r,z,rpoly=12.,momx,momy,mass,afit,bfit,cfit,r_ran,theta_ran,phi_ran
  real(pre)::vcyl,vr,fexp,vz,vx,vy,x,y,m,soft,dtheta,theta,vt,a,ecc
!
!
#ifdef WITHDRAG
  real(pre)::tg,pg,dg,tm,pm,dm,asize,rhoa
#endif
!
!
  logical::active
  character*8::cindx
  character*17::filename

  call get_units(scl)
  
  ds=sqrt(dx*dx+dy*dy+dz*dz)

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First check whether we need to read in a file.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  if(irestart>0.and. (.not.initialize_particles_now))then
   filename=""
   write(cindx,'(I8.8)')irestart
   filename="pdump."//cindx 
   open(unit=101,file=filename)
   read(101,'(1pe16.9,1X,3(I9))')starttime,step,npart,npart_direct
   print*, npart,npart_direct

   time=starttime
 
   allocate(part(npart)) !  particle in cell
   allocate(part_direct(npart_direct)) ! direct
   call srand(iseed)
   id=0

   print *, "Restarting from particle file ",filename
   do ipart=1,npart+npart_direct
!
!
#ifdef WITHDRAG
     read(101,'(I9,1X,I9,16(1X,1pe16.9),1X,L)')istat,id,x,y,z,vx,vy,vz,m,soft,rhoa,asize,dg,tg,pg,dm,tm,pm,active
#else
     read(101,'(I9,1X,I9,8(1X,1pe16.9),1X,L)')istat,id,x,y,z,vx,vy,vz,m,soft,active
#endif /* end ifdef WITHDRAG */
!
!
     if(istat==1)then
       part_direct(id-npart+1)%x=x
       part_direct(id-npart+1)%y=y
       part_direct(id-npart+1)%z=z
       part_direct(id-npart+1)%vx=vx
       part_direct(id-npart+1)%vy=vy
       part_direct(id-npart+1)%vz=vz
       part_direct(id-npart+1)%m=m
       part_direct(id-npart+1)%soft=soft
       part_direct(id-npart+1)%active=active
       part_direct(id-npart+1)%fx=zero
       part_direct(id-npart+1)%fy=zero
       part_direct(id-npart+1)%fz=zero
!
!
#ifdef WITHDRAG
       part_direct(id-npart+1)%rho0=rhoa
       part_direct(id-npart+1)%r=asize
       part_direct(id-npart+1)%d=dg
       part_direct(id-npart+1)%t=tg
       part_direct(id-npart+1)%p=pg
       part_direct(id-npart+1)%dm=dm
       part_direct(id-npart+1)%tm=tm
       part_direct(id-npart+1)%pm=pm
#endif /* end ifdef WITHDRAG */
!
!
       part_direct(id-npart+1)%id=id
     else 
       part(id+1)%x=x
       part(id+1)%y=y
       part(id+1)%z=z
       part(id+1)%vx=vx
       part(id+1)%vy=vy
       part(id+1)%vz=vz
       part(id+1)%m=m
       part(id+1)%soft=soft
       part(id+1)%active=active
       part(id+1)%fx=zero
       part(id+1)%fy=zero
       part(id+1)%fz=zero
!
!
#ifdef WITHDRAG
       part(id+1)%rho0=rhoa
       part(id+1)%r=asize
       part(id+1)%d=dg
       part(id+1)%t=tg
       part(id+1)%p=pg
       part(id+1)%dm=dm*zero
       part(id+1)%tm=tm*zero
       part(id+1)%pm=pm*zero
#endif /* end ifdef WITHDRAG */
!
!
       part(id+1)%id=id
     endif

   enddo
   close(101)


  else

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! There is no particle file to read.  We need to set up the particles
! according to user specifications.  This is done in particle_setup.inc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    if(npart>0)allocate(part(npart))
    if(npart_direct>0)allocate(part_direct(npart_direct))
    call srand(iseed)
    id=0
    isize=1

#include "particle_user_init.inc"
#include "particle_direct_user_init.inc"

  endif
 end subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Modify timestep based on particles. Take the smallest of the courant
! timestep and particle cell-crossing timestep.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 subroutine timestep_particle()
  integer::ipart
  real(pre)::fx,fy,fz,vx,vy,vz,tmp,ff,vv
!
!$OMP MASTER
  prate=zero
!$OMP END MASTER
!
!$OMP BARRIER
  if(npart>0)then
!
!
!$OMP DO SCHEDULE(STATIC) REDUCTION(max:prate)
  do ipart=1,npart
   if(.not.part(ipart)%active)cycle
   fx=part(ipart)%fx;fy=part(ipart)%fy;fz=part(ipart)%fz
   vx=part(ipart)%vx;vy=part(ipart)%vy;vz=part(ipart)%vz
   ff=sqrt(fx*fx+fy*fy+fz*fz)
   vv=sqrt(vx*vx+vy*vy+vz*vz)
   if(.not.(vv*ff)==zero)then
     prate=max(prate,ff/vv*1d2)
   else
     prate=max(prate,max(sqrt(vx**2+vy**2)/(0.25d0*dx), &
         abs(vz)/(0.25d0*dz)))
   endif
  enddo
!$OMP ENDDO
!
!$OMP MASTER
   dt=min(dt,one/prate)
   prate=zero
!$OMP END MASTER
!
!$OMP BARRIER
!
!
  endif
!
!
  if(npart_direct>0)then
!
!
!$OMP DO SCHEDULE(STATIC) REDUCTION(max:prate) private(fx,fy,fz,vx,vy,vz,ff,vv)
  do ipart=1,npart_direct
   if(.not.part_direct(ipart)%active)cycle
   fx=part_direct(ipart)%fx;fy=part_direct(ipart)%fy;fz=part_direct(ipart)%fz
   vx=part_direct(ipart)%vx;vy=part_direct(ipart)%vy;vz=part_direct(ipart)%vz
   ff=sqrt(fx*fx+fy*fy+fz*fz)
   vv=sqrt(vx*vx+vy*vy+vz*vz)
   if(.not.(vv*ff)==zero)then
     prate=max(prate,ff/vv*1d2)
   else
     prate=max(prate,max(sqrt(vx**2+ &
       vy**2)/(0.25d0*dx),abs(vz)/(0.25d0*dz)))
   endif
  enddo
!$OMP ENDDO
!
!$OMP MASTER
   dt=min(dt,one/prate)
!$OMP END MASTER
!
!$OMP BARRIER
!
!
  endif
 end subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Add particle-in-cell particles to the total cell density using a
! cloud interpolation.  This is used for self-gravity.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 subroutine add_particle_density()
  integer::ipart,ig(8),iter,np_reset_limit,nthreads
  real(pre)::x,y,z,d,vol,r
  real(pre)::w(8),step_release

  vol=dx*dy*dz
  nthreads=1
  np_reset_limit=0
  step_release=zero
!$ nthreads=omp_get_num_threads()

!$OMP DO SCHEDULE(STATIC) 
  do ipart=1,npart
    if(.not.part(ipart)%active)then
#include "particle_false_user.inc"
    endif
    x=part(ipart)%x
    y=part(ipart)%y
    z=part(ipart)%z
    d=part(ipart)%m/vol ! *100 ! 100 is for dust to gas ratio

    call get_weights_8(ig,x,y,z,w)
    if(ig(1)==-1)then
#include "particle_oob_user.inc"
    endif 
    do iter=1,8
!$OMP ATOMIC
      rhotot(ig(iter))=rhotot(ig(iter))+d*w(iter)
    enddo

  enddo
!$OMP ENDDO  

 end subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Drift step for all particles.  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 subroutine drift_particles(t)
  integer::ipart
  real(pre),intent(in)::t
  real(pre)::x,y,z
!
!$OMP MASTER
  nact1=0
  nact2=0
!$OMP ENDMASTER
!
!$OMP BARRIER
!
!
  if(npart>0)then
!
!
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:nact1)
  do ipart=1,npart
    if(.not.part(ipart)%active)cycle
    x=part(ipart)%x+part(ipart)%vx*t
    y=part(ipart)%y+part(ipart)%vy*t
    z=part(ipart)%z+part(ipart)%vz*t
    part(ipart)%x=x
    part(ipart)%y=y
    part(ipart)%z=z
    nact1=nact1+1
  enddo
!$OMP ENDDO
!
!
  endif
  if(npart_direct>0)then
!
!
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:nact2)
  do ipart=1,npart_direct
    if(.not.part_direct(ipart)%active)cycle
    x=part_direct(ipart)%x+part_direct(ipart)%vx*t
    y=part_direct(ipart)%y+part_direct(ipart)%vy*t
    z=part_direct(ipart)%z+part_direct(ipart)%vz*t
    part_direct(ipart)%x=x
    part_direct(ipart)%y=y
    part_direct(ipart)%z=z
    nact2=nact2+1
  enddo
!$OMP ENDDO
!
!
  endif
!
!
#ifdef VERBOSE
!
!
!$OMP MASTER
  print *, "Direct active particles = ",nact2
  print *, "Cloud  active particles = ",nact1
!$OMP END MASTER
!
!
#endif /* end ifdef VERBOSE */
!
! 
 end subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Kick step for all particles.  Drag terms are defined regardless
! whether drag is in use.  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 subroutine kick_particles(t)
  integer::ipart,ig(8),iter,nthreads,np_reset_limit
  real(pre),intent(in)::t
  real(pre)::x,y,z,r,vx,vy,vz,fx,fy,fz,de,ekin,t0,tR
  real(pre)::w(8),azimuth
  real(pre)::dfx,dfy,dfz,d,pg,tg,dg,vol,rhoa,asize,d_loc,tdrag
  real(pre)::pfx,pfy,pfz,alpha_loc,omega,omega2,step_release
  integer::niter,idrag
  nthreads=1
  np_reset_limit=0
  step_release=zero
!$ nthreads=omp_get_num_threads()

#ifdef ROTATE
  real(pre)::x0
      x0=object_x_displace
      omega2=(object_mass)/abs(object_x_displace)**3
      omega=sqrt(omega2)
#endif
  vol=dx*dy*dz
!
!
  if(npart>0)then
!
!
  call particle_direct_interp()
!
!$OMP BARRIER
!
!
!$OMP DO SCHEDULE(STATIC)                                               &
!$OMP&PRIVATE(x,y,z,vx,vy,vz,fx,fy,fz,ig,w,pfx,pfy,pfz,rhoa,d)          &
!$OMP&PRIVATE(tg,pg,dg,alpha_loc,asize,de,ekin,t0,tR,iter,r)
  do ipart=1,npart
    if(.not.part(ipart)%active)then
#include "particle_false_user.inc"
    endif

    x=part(ipart)%x
    y=part(ipart)%y
    z=part(ipart)%z
    vx=part(ipart)%vx
    vy=part(ipart)%vy
    vz=part(ipart)%vz
    fx=part(ipart)%fx
    fy=part(ipart)%fy
    fz=part(ipart)%fz

    call get_weights_8(ig,x,y,z,w)
    if(ig(1)==-1)then
#include "particle_oob_user.inc"
    endif 
 
    do iter=1,8
     fx=fx+gforce(1,ig(iter))*w(iter)
     fy=fy+gforce(2,ig(iter))*w(iter)
     fz=fz+gforce(3,ig(iter))*w(iter)
    enddo

#ifdef ROTATE
     fx=fx+two*vy*omega
     fx=fx-1.5d0*omega*y/x0*x0dot
     fy=fy-two*vx*omega
     fy=fy+1.5d0*omega*x/x0*x0dot
#endif

    vx=vx+fx*t
    vy=vy+fy*t
    vz=vz+fz*t
!
!
#ifdef WITHDRAG
!
!
    pfx=zero;pfy=zero;pfz=zero
    do iter=1,8
     pfx=pfx-pforce(1,ig(iter))*w(iter)
     pfy=pfy-pforce(2,ig(iter))*w(iter)
     pfz=pfz-pforce(3,ig(iter))*w(iter)
    enddo


    d=part(ipart)%m/vol
    rhoa=part(ipart)%rho0
    asize=part(ipart)%r
    call get_drag(dfx,dfy,dfz,tg,pg,dg,x,y,z,vx,vy,vz, &
                  d,rhoa,asize,ig,w,t,pfx,pfy,pfz,alpha_loc,de) ! pdrag.f90
 
    de=max(de,zero)
    do iter=1,8
     ekin=half*(cons(2,ig(iter))**2+cons(3,ig(iter))**2+cons(4,ig(iter))**2)/cons(1,ig(iter))
!$OMP ATOMIC
     cons(5,ig(iter))=cons(5,ig(iter))-ekin
!$OMP ATOMIC
     cons(2,ig(iter))=cons(2,ig(iter))-d*dfx*w(iter)*t! set to zero for step in calc_force
!$OMP ATOMIC
     cons(3,ig(iter))=cons(3,ig(iter))-d*dfy*w(iter)*t ! set to zero for step in calc_force
!$OMP ATOMIC
     cons(4,ig(iter))=cons(4,ig(iter))-d*dfz*w(iter)*t ! set to zero for step in calc_force
!
     ekin=half*(cons(2,ig(iter))**2+cons(3,ig(iter))**2+cons(4,ig(iter))**2)/cons(1,ig(iter))
!$OMP ATOMIC
     cons(5,ig(iter))=cons(5,ig(iter))+ekin+(dg)*de*w(iter)
    enddo

    vx=vx+(dfx)*t
    vy=vy+(dfy)*t
    vz=vz+(dfz)*t

    part(ipart)%t=tg

#ifdef THERMALHIST
    part(ipart)%t=tg
    part(ipart)%p=pg
    part(ipart)%d=dg
    if(part(ipart)%tm<tg)part(ipart)%tm=tg
    if(part(ipart)%pm<pg)part(ipart)%pm=pg
    if(part(ipart)%dm<dg)part(ipart)%dm=dg
!
!
#endif /* end ifdef THERMALHIST */
!
!
#endif /* end ifdef WITHDRAG */
!
!
#ifdef ZERO_OUT_Z
    vz=zero
#endif
!
    part(ipart)%vx=vx
    part(ipart)%vy=vy
    part(ipart)%vz=vz

  enddo
!$OMP ENDDO
!
!
  endif
!
!
  if(npart_direct>0)then
!
!
  call particle_direct_direct()
!
!
!$OMP BARRIER
!$OMP DO SCHEDULE(STATIC)                                               &
!$OMP&PRIVATE(x,y,z,vx,vy,vz,fx,fy,fz,ig,w,pfx,pfy,pfz,rhoa,d)          &
!$OMP&PRIVATE(tg,pg,dg,alpha_loc,asize,de)
  do ipart=1,npart_direct
    if(.not.part_direct(ipart)%active)cycle
    vy=part_direct(ipart)%vy
    vx=part_direct(ipart)%vx
    vz=part_direct(ipart)%vz
    fy=part_direct(ipart)%fy
    fx=part_direct(ipart)%fx
    fz=part_direct(ipart)%fz

  if(npart>0)then
    x=part_direct(ipart)%x
    y=part_direct(ipart)%y
    z=part_direct(ipart)%z
 
    call get_weights_8(ig,x,y,z,w)
    if(ig(1)==-1)then
#include "particle_direct_oob_user.inc"
    endif 

    do iter=1,8
     fx=fx+gforce(1,ig(iter))*w(iter)
     fy=fy+gforce(2,ig(iter))*w(iter)
     fz=fz+gforce(3,ig(iter))*w(iter)
    enddo

#ifdef ROTATE
     fx=fx+two*vy*omega
     fx=fx-1.5d0*omega*y/x0*x0dot
     fy=fy-two*vx*omega
     fy=fy+1.5d0*omega*x/x0*x0dot
#endif



    vx=vx+fx*t
    vy=vy+fy*t
    vz=vz+fz*t
!
!
#ifdef WITHDRAG
!
!
    pfx=zero;pfy=zero;pfz=zero
    do iter=1,8
     pfx=pfx-pforce(1,ig(iter))*w(iter)
     pfy=pfy-pforce(2,ig(iter))*w(iter)
     pfz=pfz-pforce(3,ig(iter))*w(iter)
    enddo

    d=part_direct(ipart)%m/vol
    rhoa=part_direct(ipart)%rho0
    asize=part_direct(ipart)%r
    call get_drag(dfx,dfy,dfz,tg,pg,dg,x,y,z,vx,vy,vz, &
                  d,rhoa,asize,ig,w,t,pfx,pfy,pfz,alpha_loc,de)

    do iter=1,8
!$OMP ATOMIC
     cons(2,ig(iter))=cons(2,ig(iter))-d*dfx*w(iter)*t! set to zero for step in calc_force
!$OMP ATOMIC
     cons(3,ig(iter))=cons(3,ig(iter))-d*dfy*w(iter)*t ! set to zero for step in calc_force
!$OMP ATOMIC
     cons(4,ig(iter))=cons(4,ig(iter))-d*dfz*w(iter)*t ! set to zero for step in calc_force
    enddo

    vx=vx+(dfx)*t
    vy=vy+(dfy)*t
    vz=vz+(dfz)*t
!
!
#ifdef THERMALHIST
    part_direct(ipart)%t=tg
    part_direct(ipart)%p=pg
    part_direct(ipart)%d=dg
    if(part_direct(ipart)%tm<tg)part_direct(ipart)%tm=tg
    if(part_direct(ipart)%pm<pg)part_direct(ipart)%pm=pg
    if(part_direct(ipart)%dm<dg)part_direct(ipart)%dm=dg
#endif /* end ifdef THERMALHIST */
!
!
#endif /* end ifdef WITHDRAG */
!
!
else
    vx=vx+fx*t
    vy=vy+fy*t
    vz=vz+fz*t
endif
    part_direct(ipart)%vy=vy
    part_direct(ipart)%vx=vx
    part_direct(ipart)%vz=vz
  enddo
!$OMP ENDDO
endif
 end subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Find gravitational force using direct sum among selected particles.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 subroutine particle_direct_direct()
  integer::ipart,jpart
  real(pre)::xi,yi,zi,mi,si,xj,yj,zj,mj,sj,r,fxi,fyi,fzi,fxj,fyj,fzj
 
  if(npart_direct>0)then
!
!
!$OMP DO SCHEDULE(STATIC)
  do ipart=1,npart_direct
   part_direct(ipart)%fx=zero
   part_direct(ipart)%fy=zero
   part_direct(ipart)%fz=zero
  enddo
!$OMP ENDDO
!$OMP BARRIER
!
!$OMP DO SCHEDULE(STATIC)
  do ipart=1,npart_direct
   if(.not.part_direct(ipart)%active)cycle
   xi=part_direct(ipart)%x
   yi=part_direct(ipart)%y
   zi=part_direct(ipart)%z
   mi=part_direct(ipart)%m
   si=part_direct(ipart)%soft
   do jpart=ipart+1,npart_direct
     if(.not.part_direct(jpart)%active)cycle
     xj=part_direct(jpart)%x
     yj=part_direct(jpart)%y
     zj=part_direct(jpart)%z
     mj=part_direct(jpart)%m
     sj=part_direct(jpart)%soft
     r=sqrt( (xi-xj)**2+(yi-yj)**2+(zi-zj)**2 )
     if(r<sj)then
      fxi=mj/sj**3*( three*r/sj-four)*(xi-xj)
      fyi=mj/sj**3*( three*r/sj-four)*(yi-yj)
      fzi=mj/sj**3*( three*r/sj-four)*(zi-zj)

      if(r<si)then
       fxj=-mi/si**3*( three*r/si-four)*(xi-xj)
       fyj=-mi/si**3*( three*r/si-four)*(yi-yj)
       fzj=-mi/si**3*( three*r/si-four)*(zi-zj)
      else
       fxj=mi/r**3*(xi-xj)
       fyj=mi/r**3*(yi-yj)
       fzj=mi/r**3*(zi-zj)
      endif
     else
      fxi=-mj/r**3*(xi-xj)
      fyi=-mj/r**3*(yi-yj)
      fzi=-mj/r**3*(zi-zj)

      if(r<si)then
       fxj=-mi/si**3*( three*r/si-four)*(xi-xj)
       fyj=-mi/si**3*( three*r/si-four)*(yi-yj)
       fzj=-mi/si**3*( three*r/si-four)*(zi-zj)
      else
       fxj=mi/r**3*(xi-xj)
       fyj=mi/r**3*(yi-yj)
       fzj=mi/r**3*(zi-zj)
      endif
     endif
!$OMP ATOMIC
     part_direct(ipart)%fx=part_direct(ipart)%fx+fxi
!$OMP ATOMIC
     part_direct(ipart)%fy=part_direct(ipart)%fy+fyi
!$OMP ATOMIC
     part_direct(ipart)%fz=part_direct(ipart)%fz+fzi

!$OMP ATOMIC
     part_direct(jpart)%fx=part_direct(jpart)%fx+fxj
!$OMP ATOMIC
     part_direct(jpart)%fy=part_direct(jpart)%fy+fyj
!$OMP ATOMIC
     part_direct(jpart)%fz=part_direct(jpart)%fz+fzj
   enddo
  enddo
!$OMP ENDDO
!
!
  endif
 end subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Add force of direct sum particles on PIC particles.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 subroutine particle_direct_interp()
  integer::ipart,jpart
  real(pre)::xi,yi,zi,xj,yj,zj,mi,si,sj,mj,r,fxi,fyi,fzi
 
!$OMP DO SCHEDULE(STATIC)
  do ipart=1,npart
   part(ipart)%fx=zero
   part(ipart)%fy=zero
   part(ipart)%fz=zero
  enddo
!$OMP ENDDO
!$OMP BARRIER
  if(npart_direct>0)then
!$OMP DO SCHEDULE(STATIC)
  do ipart=1,npart_direct
   if(.not.part_direct(ipart)%active)cycle
   xi=part_direct(ipart)%x
   yi=part_direct(ipart)%y
   zi=part_direct(ipart)%z
   mi=part_direct(ipart)%m
   si=part_direct(ipart)%soft
   do jpart=1,npart
     if(.not.part(jpart)%active)cycle
     xj=part(jpart)%x
     yj=part(jpart)%y
     zj=part(jpart)%z
     r=sqrt( (xi-xj)**2+(yi-yj)**2+(zi-zj)**2 )
     if(r<si)then
      fxi=-mi/si**3*( three*r/si-four)*(xi-xj)
      fyi=-mi/si**3*( three*r/si-four)*(yi-yj)
      fzi=-mi/si**3*( three*r/si-four)*(zi-zj)

     else
      fxi=mi/r**3*(xi-xj)
      fyi=mi/r**3*(yi-yj)
      fzi=mi/r**3*(zi-zj)

     endif
!$OMP ATOMIC
     part(jpart)%fx=part(jpart)%fx+fxi
!$OMP ATOMIC
     part(jpart)%fy=part(jpart)%fy+fyi
!$OMP ATOMIC
     part(jpart)%fz=part(jpart)%fz+fzi

   enddo
  enddo
!$OMP ENDDO
  endif
 end subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Add potential from direct particles to the grid.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 subroutine add_direct_togrid()
  integer::jpart,igrid
  real(pre)::xi,yi,zi,xj,yj,zj,sj,mj,r,fxi,fyi,fzi

  if(npart_direct>0)then
!
!
#ifdef NOHYDRO
 if(.not.npart>0.or..not.use_pic)return
#endif 
!
!
!$OMP DO SCHEDULE(STATIC)
  do igrid=1,ngrid
   xi=grid(igrid)%x
   yi=grid(igrid)%y
   zi=grid(igrid)%z
   do jpart=1,npart_direct
     if(.not.part_direct(jpart)%active)cycle
     xj=part_direct(jpart)%x
     yj=part_direct(jpart)%y
     zj=part_direct(jpart)%z
     mj=part_direct(jpart)%m
     sj=part_direct(jpart)%soft
     r=sqrt( (xi-xj)**2+(yi-yj)**2+(zi-zj)**2 )
     if(r<sj)then
      fxi=mj/sj**3*( three*r/sj-four)*(xi-xj)
      fyi=mj/sj**3*( three*r/sj-four)*(yi-yj)
      fzi=mj/sj**3*( three*r/sj-four)*(zi-zj)
     else
      fxi=-mj/r**3*(xi-xj)
      fyi=-mj/r**3*(yi-yj)
      fzi=-mj/r**3*(zi-zj)
     endif
!$OMP ATOMIC
     gforce(1,igrid)=gforce(1,igrid)+fxi
!$OMP ATOMIC
     gforce(2,igrid)=gforce(2,igrid)+fyi
!$OMP ATOMIC
     gforce(3,igrid)=gforce(3,igrid)+fzi
   enddo
  enddo
!$OMP ENDDO
  endif
 end subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! logfile diagnostic for particle tracking at every time step
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 subroutine print_select_particles()
 integer::ipart
!$OMP MASTER
 if(npart>0)then
  do ipart=1,1!5!npart
!   if (ipart==51) then
!   if (ipart<=npart) then
!
!
#ifdef THERMALHIST
      print "(A10,I6,1X,10(1pe16.8))", "PARTICLE:",ipart,time,part(ipart)%x,part(ipart)%y, &
        part(ipart)%z,part(ipart)%vx,part(ipart)%vy,part(ipart)%vz,&
        part(ipart)%d, part(ipart)%t,part(ipart)%p
#else
      print "(A10,I6,1X,7(1pe16.8))", "PARTICLE:",ipart,time,part(ipart)%x,part(ipart)%y, &
        part(ipart)%z,part(ipart)%vx,part(ipart)%vy,part(ipart)%vz
#endif /* end ifdef THERMALHIST */
!
!
!   endif
  enddo
 endif
!$OMP END MASTER
 end subroutine

end module

