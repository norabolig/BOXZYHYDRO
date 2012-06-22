program boxzyhydro
 use parameters
 use derived_types
 use grid_commons
 use eos
 use input
#ifdef RADTRAN
 use mcrtfld
#endif
#ifdef SELFGRAVITY
 use selfgravity
#endif
#ifdef EXTERNALPHI
 use selfgravity
#endif
#ifdef PARTICLE
 use particle
#endif
 implicit none

 integer::igrid,idim,step,idx,ibd,ipart
 real(pre)::mass,eint,ekin,etot,amom,angle,x,y,z,momx,momy,egra,phian,r
 real(pre)::nextstop,cons_swap,f1,max_den_change_loc,allow=0.01,timer_start,timer_stop
 real(pre)::max_den_loc,max_den_old_loc,dtold,xcom,ycom,zcom,momz
 type(units)::scale
 character*8::cindx
 character*17::filename
 character*80::namelist_file
 integer::readtest=1234

 logical:: checkden=.false.

! debugging declarations
 real(pre)::testphi,xx,yy,zz
 integer::jgrid

 logical:: first=.true.,f1_cycle

 call getarg(1,namelist_file)
 call read_params(namelist_file)
 call init_grid()
 call initialize_eos
 call calc_eos_table
#ifdef RADTRAN
 call initialize_mcrtfld()
#endif
!$OMP PARALLEL DEFAULT(SHARED)
#ifndef NOHYDRO
 if(irestart>0)then
  call read_hydro()
 else
  call init_conditions()
 endif
#else
  call init_conditions()
#endif
  call avisc()
!$OMP END PARALLEL
#ifdef PARTICLE
 call initialize_particles()
#endif
 time=starttime
 nextstop=starttime
#ifdef SELFGRAVITY
 call init_grav_grid
!$OMP PARALLEL
!$OMP DO SCHEDULE(STATIC)
 do igrid=1,ngrid
   rhotot(igrid)=cons(1,igrid)
 enddo
!$OMP ENDDO
#ifdef PARTICLE
 call add_particle_density()
#endif
!$OMP END PARALLEL
 call set_com_in_tree
#ifdef FASTGRAVITY
 call get_pot_bc_from_tree
 call vcycle_pot() !sor_potential2(sor_iter)
#else
 call get_pot_from_tree
print *, "Got pot using expansion"
#endif
#if VERBOSE
#endif
#ifdef USEPERT
 rhophi=phi
#if VERBOSE
 print *,"Should call rochepert"
#endif
 call rochepert()
#endif
#endif

#ifdef EXTERNALPHI
  call external_phi()
#endif

#ifdef NOHYDRO
  if(.not.npart>0.or..not.use_pic)phi=zero
#endif

#if 1 == 0
  call test_phi()
#endif

 call get_units(scale)

 f1=one
 dtold=zero
 step=0
 if(irestart>0)step=irestart
#ifdef VERBOSE
print *, "#Starting time integration"
#endif
call cpu_time(timer_start)
do 
 if(time>=endtime.or.step>=1000000)exit
 step=step+1

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
 max_den_loc=zero
 max_den_old_loc=zero
!$OMP  PARALLEL DEFAULT(SHARED) private(x,y,z,angle)
#ifndef NOHYDRO
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:mass,ekin,eint,amom,momx,momy,momz,etot) &
!$OMP REDUCTION(+:egra,xcom,ycom,zcom) 
 do igrid=1,ngrid
  x=grid(igrid)%x
  y=grid(igrid)%y
  z=grid(igrid)%z
  angle=atan2(y,x)
  amom=amom+(cons(3,igrid)*cos(angle)-cons(2,igrid)*sin(angle))*sqrt(x**2+y**2)*dx*dy*dz
  momx=momx+cons(2,igrid)*dz*dx*dy
  momy=momy+cons(3,igrid)*dz*dx*dy
  momz=momz+cons(4,igrid)*dz*dx*dy
  mass=mass+cons(1,igrid)*dz*dx*dy
  xcom=xcom+x*cons(1,igrid)*dz*dx*dy
  ycom=ycom+y*cons(1,igrid)*dz*dx*dy
  zcom=zcom+z*cons(1,igrid)*dz*dx*dy
  etot=etot+(cons(5,igrid))*dz*dx*dy
  ekin=ekin+(cons(2,igrid)**2+cons(3,igrid)**2+cons(4,igrid)**2)/(two*cons(1,igrid))*dx*dy*dz
  egra=egra+half*(cons(1,igrid)*phi(igrid))*dx*dy*dz
 enddo
!$OMP ENDDO

#ifdef SUPPRESSDRIFT
  call suppress_drift(xcom,ycom,zcom)
#endif

!$OMP DO SCHEDULE(STATIC) REDUCTION(max:max_den_old_loc)
 do igrid=1,ngrid
  max_den_old_loc=max(max_den_old_loc,cons(1,igrid))
 enddo
!$OMP ENDDO
!$OMP DO SCHEDULE(STATIC)
 do igrid=1,ngrid
 do idim=1,5
   cons_old(idim,igrid)=cons(idim,igrid)
 enddo
 enddo
!$OMP ENDDO

#ifdef VERBOSE
!$OMP MASTER
 print *,"# Total mass   is ",mass,time
 print *,"# Total momx,y,z is ",momx,momy,momz,time
 print *,"# Total amom   is ",amom,time
 print *,"# Total energy is ",etot+egra,etot,etot-ekin,ekin,egra,time
 print *,"#COM ",time,xcom/mass,ycom/mass,zcom/mass
!$OMP END MASTER
#endif

 call state()

!$omp barrier
 
 call courant()
!$OMP BARRIER
#endif 
! NOHYDRO
#ifdef PARTICLE
#ifdef VERBOSE
!$OMP MASTER
 if(npart>0)then
 print *, "PARTICLE1: ",time,part(1)%x,part(1)%y,part(1)%z, &
        part(1)%vx,part(1)%vy,part(1)%vz!,part(1)%d,part(1)%t
 if(npart_direct>0)then
 print *, "PARTICLED1: ",time,part_direct(1)%x,part_direct(1)%y,part_direct(1)%z, &
        part_direct(1)%vx,part_direct(1)%vy,part_direct(1)%vz
 endif
 endif
!$OMP END MASTER
#endif
#ifdef NOHYDRO
!$OMP MASTER
 dt=endtime
!$OMP END MASTER
!$OMP BARRIER
#endif
 call timestep_particle()
#endif 
! end for PARTICLE
!$OMP MASTER
 dt=dt*f1
 if(step==1)dt=1d-1*dt
 if(step>1.and.dtold>zero)dt=min(dt,dtold*1.1d0)
 dtold=dt
#ifdef VERBOSE
 print *,"#dt and Time",dt,time," at step ",step,f1
#endif
 if(dt+time>endtime)dt=endtime-time
 dt=dt*half
!$OMP END MASTER
!$OMP BARRIER
 call source()
!$OMP BARRIER
 
#ifndef NOHYDRO
#ifdef RADTRAN
 call mcrtfld_transfer
#endif
!$OMP BARRIER
 call velocity()
!$OMP BARRIER
 
!$OMP DO SCHEDULE(STATIC)
 do igrid=1,ngrid
 do idim=1,5
   cons_new(idim,igrid)=cons(idim,igrid)
 enddo
 enddo
!$OMP ENDDO

#endif 
! NOHYDRO
!$OMP END PARALLEL

#ifndef NOHYDRO

 cons_pt=>cons

 call flux()

!$OMP PARALLEL 
!$OMP DO SCHEDULE(STATIC) REDUCTION(max:max_den_loc)
 do igrid=1,ngrid
  max_den_loc=max(max_den_loc,cons(1,igrid))
 enddo
!$OMP ENDDO
!$OMP MASTER
 f1_cycle=.false.
 max_den_change_loc=abs(one-max_den_loc/max_den_old_loc)
 if(checkden.and.max_den_change_loc>allow)then
   f1=allow/max_den_change_loc*f1*half
   print *, "# Density change of ",  max_den_change_loc, " detected. Allowing ", allow,".",f1
   f1_cycle=.true.
 else
   f1=min(f1*1.1d0,one)
 endif
!$OMP END MASTER
!$OMP BARRIER
 if(f1_cycle)then
!$OMP DO SCHEDULE(STATIC)
 do igrid=1,ngrid
 do idim=1,5
   cons(idim,igrid)=cons_old(idim,igrid)
 enddo
 enddo
!$OMP ENDDO
 call velocity()
 endif
!$OMP END PARALLEL
#endif 
!NOHYDRO

 if(f1_cycle)cycle

!$OMP PARALLEL DEFAULT(SHARED) private(cons_swap)
#ifndef NOHYDRO
 call velocity()

#endif
!NOHYDRO
!$OMP END PARALLEL

 dt=two*dt

#ifndef NOHYDRO
 cons_pt=>cons_new
 call flux()

!$OMP PARALLEL DEFAULT(SHARED) 
!$OMP DO SCHEDULE(STATIC)
  do igrid=1,ngrid
   do idim=1,5
    cons(idim,igrid)=cons_pt(idim,igrid)
   enddo
  enddo
!$OMP ENDDO
!$OMP END PARALLEL
!NOHYDRO
#endif


#ifdef SELFGRAVITY
!$OMP PARALLEL
#ifdef PARTICLE
 call drift_particles(dt)
#ifdef NOHYDRO
 if(npart>0.and.use_pic)then
#endif
#endif
!$OMP DO SCHEDULE(STATIC)
 do igrid=1,ngrid
   rhotot(igrid)=cons(1,igrid)
 enddo
!$OMP ENDDO
#ifdef PARTICLE
!$OMP BARRIER
 call add_particle_density()
#endif
#ifdef NOHYDRO
 endif ! do not bother with gravity if npart !>0 if NOHYDRO defined
#endif
!$OMP END PARALLEL
#ifdef NOHYDRO
 if(npart>0.and.use_pic)then
#endif
 call set_com_in_tree
#ifdef USEPERT
!$OMP PARALLEL 
!$OMP DO SCHEDULE(STATIC)
 do igrid=1,ngrid
   phi(igrid)=rhophi(igrid)
 enddo
!$OMP END DO
!$OMP END PARALLEL
#endif
#ifdef FASTGRAVITY
 call get_pot_bc_from_tree
 call vcycle_pot() !sor_potential2(sor_iter)
#else
 call get_pot_from_tree
#endif
#ifdef USEPERT
!$OMP PARALLEL DO SCHEDULE(STATIC)
 do igrid=1,ngrid
   rhophi(igrid)=phi(igrid)
 enddo
!$OMP END PARALLEL DO
 call rochepert()
#endif
#ifdef NOHYDRO
 endif
#endif 
#endif
! selfgravity

#ifdef EXTERNALPHI
  call external_phi()
#endif

 dt=half*dt

!$OMP PARALLEL DEFAULT(SHARED) 
#ifndef NOHYDRO
 call velocity() ! update from last flux.  Done after self-gravity.
 call state()
#endif
 call source()
#ifndef NOHYDRO
#ifdef RADTRAN
 call mcrtfld_transfer
#endif
 call velocity()
#endif
!$OMP END PARALLEL

 dt=dt*two
 time=time+dt


 if((time>=nextstop.or.time==endtime))then
#ifndef NOHYDRO
!$OMP PARALLEL DEFAULT(SHARED)
  call state
!$OMP END PARALLEL
  write(cindx,'(I8.8)')step  
  filename="celldump."//cindx
  print *, filename
  if(write_ascii)then
     open(unit=100,file=filename)
     write(100,'(A25,1pe16.8E3,I8)')"#Time of snapshot at step:",time,step
     write(100,'(A,1X,8(1pe15.8,1X))')'#',scale
     do igrid=1,ngrid
       ibd=0
       if(grid(igrid)%boundary)ibd=1
       !if(grid(igrid)%boundary)cycle
        write(100,'(11(1pe16.8E3,1X),I2)'),grid(igrid)%x,grid(igrid)%y,grid(igrid)%z,&
           cons(1,igrid),p(igrid),&
           u(1,igrid),u(2,igrid),u(3,igrid),phi(igrid),cons(5,igrid),muc_array(igrid),ibd
     enddo
   else
     open(unit=100,file=filename,form='UNFORMATTED')
     write(100)readtest
     write(100)time,step
     write(100)scale
     do igrid=1,ngrid
       ibd=0
       if(grid(igrid)%boundary)ibd=1
       !if(grid(igrid)%boundary)cycle
        write(100)grid(igrid)%x,grid(igrid)%y,grid(igrid)%z,&
           cons(1,igrid),p(igrid), &
           u(1,igrid),u(2,igrid),u(3,igrid),phi(igrid),cons(5,igrid),muc_array(igrid),ibd
     enddo
   endif
   close(100)
 
#endif
#ifdef PARTICLE
   filename=""
   write(cindx,'(I8.8)')step  
   filename="pdump."//cindx
   print *, filename
   open(unit=100,file=filename)
   write(100,'(1pe16.9,1X,3(I9))')time,step,npart,npart_direct
   if(npart_direct>0)then
     do ipart=1,npart_direct
#ifdef WITHDRAG
       write(100,'(I9,1X,I9,16(1X,1pe16.9),1X,L)')1,part_direct(ipart)%id, &
         part_direct(ipart)%x,part_direct(ipart)%y, &
         part_direct(ipart)%z,part_direct(ipart)%vx,part_direct(ipart)%vy, &
         part_direct(ipart)%vz,part_direct(ipart)%m,part_direct(ipart)%soft, &
         part_direct(ipart)%rho0, &
         part_direct(ipart)%r, &
         part_direct(ipart)%d, &
         part_direct(ipart)%t, &
         part_direct(ipart)%p, &
         part_direct(ipart)%dm, &
         part_direct(ipart)%tm, &
         part_direct(ipart)%pm, &
         part_direct(ipart)%active
#else
       write(100,'(I9,1X,I9,8(1X,1pe16.9),1X,L)')1,part_direct(ipart)%id, &
         part_direct(ipart)%x,part_direct(ipart)%y, &
         part_direct(ipart)%z,part_direct(ipart)%vx,part_direct(ipart)%vy, &
         part_direct(ipart)%vz,part_direct(ipart)%m,part_direct(ipart)%soft, &
         part_direct(ipart)%active
#endif
     enddo
   endif
   if(npart>0)then
     do ipart=1,npart
#ifdef WITHDRAG
       write(100,'(I9,1X,I9,16(1X,1pe16.9),1X,L)')0,part(ipart)%id,part(ipart)%x,&
         part(ipart)%y,part(ipart)%z,part(ipart)%vx,part(ipart)%vy,part(ipart)%vz, &
         part(ipart)%m,part(ipart)%soft,&
         part(ipart)%rho0, &
         part(ipart)%r, &
         part(ipart)%d, &
         part(ipart)%t, &
         part(ipart)%p, &
         part(ipart)%dm, &
         part(ipart)%tm, &
         part(ipart)%pm, &
         part(ipart)%active
#else
       write(100,'(I9,1X,I9,8(1X,1pe16.9),1X,L)')0,part(ipart)%id,part(ipart)%x,&
         part(ipart)%y,part(ipart)%z,part(ipart)%vx,part(ipart)%vy,part(ipart)%vz, &
         part(ipart)%m,part(ipart)%soft,&
         part(ipart)%active
#endif 
     enddo
   endif
   close(100)
#endif
   nextstop=nextstop+dtout
  endif

enddo
call cpu_time(timer_stop)
print *, "Elapsed Time is : ",timer_stop-timer_start

#ifndef NOHYDRO
 call cleanup
#endif

end program
   
