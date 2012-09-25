!
! boxzy hydro.  A Cartesian FVM hydrodynamics code with:
!  o radiation transport
!  o self-gravity
!  o a detailed equation of state, plus simplified versions
!  o particles
!    * gas-drag coupling, including feedback reactions
!    * particle-in-cell + particle-particle N-body
!  o flux linear or angular momentum 
!  o Tadmore & Kalgnor fluxing with min-mod flux limiter
!  o arbitrary boundaries for creating wind tunnels, etc.
!
program boxzyhydro
!
! basic modules
!
 use parameters
 use derived_types
 use grid_commons
 use eos
 use input
 use utils
! 
! extra modules
!
#ifdef RADTRAN
 use mcrtfld
#endif
!
!
#if defined SELFGRAVITY || defined EXTERNALPHI
 use selfgravity
#endif
!
!
#ifdef PARTICLE
 use particle
#endif
!
!
 implicit none

 integer::igrid,idim,step,idx,ibd,ipart
 real(pre)::angle,x,y,z,phian,r
 real(pre)::nextstop,cons_swap,f1,max_den_change_loc,allow=0.01,timer_start,timer_stop
 real(pre)::max_den_loc,min_den_old_loc,min_den_loc,max_den_old_loc,dtold
 real(pre)::saback,safron,satop,sabot,pback,pfron,ptop,pbot,rhoav,vxav,avcount,chord=150d0
 type(units)::scale
 character*8::cindx
 character*17::filename
 character*80::namelist_file
 integer::iter

 logical:: checkden=.false.

! debugging declarations
 real(pre)::testphi,xx,yy,zz
 integer::jgrid

 logical:: first=.true.,f1_cycle

 call getarg(1,namelist_file) ! read command line for namelist file
 call read_params(namelist_file) ! input.f90
 call init_grid() ! init_grid.f90
 call initialize_eos ! eos.f90
 call calc_eos_table  ! eos.f90
!
!
#ifdef RADTRAN
 call initialize_mcrtfld() ! mcrtfld.f90
#endif
!
!
!$OMP PARALLEL DEFAULT(SHARED)
!
!
#ifndef NOHYDRO
 if(irestart>0)then
  call read_hydro() ! read_hydro.f90
 else
  call init_conditions() ! init_conditions.f90
 endif
#else
  call init_conditions()
#endif /* NOHYDRO */ 
! 
!
!$OMP END PARALLEL
!
!
#ifdef PARTICLE
 call initialize_particles() ! particles.f90
#endif
!
!
 time=starttime
 nextstop=starttime
 if (irestart>0)nextstop=nextstop+dtout ! write next output in time=time0+dtout
!
!
#ifdef SELFGRAVITY
!
!
 call init_grav_grid() ! gravity.f90
!
!$OMP PARALLEL
!$OMP DO SCHEDULE(STATIC)
 do igrid=1,ngrid
   rhotot(igrid)=cons(1,igrid)
!  rhotot is necessary for the self-gravity solver with particles treated as clouds. 
 enddo
!$OMP ENDDO
!
!
#ifdef PARTICLE
 call add_particle_density() ! particle.f90
#endif
!
!
!$OMP END PARALLEL
 call set_com_in_tree ! gravity.f90
!
!
#ifdef FASTGRAVITY
!
 call get_pot_bc_from_tree ! gravity.f90
 call vcycle_pot() ! gravity.f90
!
#else
!
 call get_pot_from_tree ! gravity.f90
!
!
#ifdef VERBOSE
print *, "Got pot using expansion"
#endif
!
#endif /* endif FASTGRAVITY */
!
!
#ifdef USEPERT
!
!
#ifdef VERBOSE
 print *,"Should call rochepert"
#endif
!
!
 call rochepert() ! gravity.f90
!
!
#endif /* end USEPERT */
!
!
#endif /* endif SELFGRAVITY */
!
!
#ifdef EXTERNALPHI
  call external_phi() ! gravity.f90, user must alter routine appropriately.
#endif
!
!
#ifdef NOHYDRO
  if(.not.npart>0.or..not.use_pic)phi=zero
#endif
!
!
#ifdef RUN_PHI_TEST
  call test_phi() ! test for gravity solver. Routine in gravity.f90
#endif
!
!
 call get_units(scale)

 f1=one
 dtold=zero
 step=0
 if(irestart>0)step=irestart
#ifdef VERBOSE
 print *, "Starting time integration"
#endif
call cpu_time(timer_start)
!
!***
! 
! This is the main loop.
!
!***
!
do 
 if(time>=endtime.or.step>=1000000)exit ! master loop control
 step=step+1

 max_den_old_loc=zero
!$OMP PARALLEL DEFAULT(SHARED) private(x,y,z,angle)                    
#ifndef NOHYDRO
  call conservation_diagnostic() ! in utils.f90 xcom,ycom,zcom defined in utils
!
!
#ifdef SUPPRESSDRIFT
  call suppress_drift(xcom,ycom,zcom) ! in utils.f90
#endif 
!
!
!$OMP DO SCHEDULE(STATIC) REDUCTION(max:max_den_old_loc)
 do igrid=1,ngrid
  max_den_old_loc=max(max_den_old_loc,cons(1,igrid))
 enddo
!$OMP ENDDO
!$OMP MASTER
  min_den_old_loc=max_den_old_loc
!$OMP END MASTER
!$OMP DO SCHEDULE(STATIC) REDUCTION(min:min_den_old_loc)
 do igrid=1,ngrid
  min_den_old_loc=min(min_den_old_loc,cons(1,igrid))
 enddo
!$OMP ENDDO
!
!$OMP DO SCHEDULE(STATIC)
 do igrid=1,ngrid
  if(mod(step,25)==zero.and.time<250)then;cons(2:4,igrid)=zero
  elseif(mod(step,50)==zero.and.time<500)then;cons(2:4,igrid)=zero
  elseif(mod(step,100)==zero.and.time<1000)then;cons(2:4,igrid)=zero
  elseif(mod(step,200)==zero.and.time<1500)then;cons(2:4,igrid)=zero
  endif
 do idim=1,5
   cons_old(idim,igrid)=cons(idim,igrid)
 enddo
 enddo
!$OMP ENDDO

 call state() ! state.f90

!$OMP BARRIER
!
!***
! Equation of state and velocities are now current. 
! It is time to find the courant condition.
!***
!
 call courant() ! courant.f90
!$OMP BARRIER
!
!
#endif /* endif NOHYDRO */
!
!
#ifdef PARTICLE
!
!
 call print_select_particles() ! particle.f90
!
!
#ifdef NOHYDRO /* If no hydro, limit the starting point for the timestep solver. */
!$OMP MASTER
 dt=endtime
!$OMP END MASTER
!$OMP BARRIER
#endif /* endif NOHYDRO */
!
!
 call timestep_particle() ! particle.f90
#endif /* endif PARTICLE */
!
!
!$OMP MASTER  
 dt=dt*f1
 if(step==1)dt=1d-3*dt ! use this with caution.  It is best for models that are not quite in equilibrium
                       ! and need a kick-start for stability.
 if(step>1.and.dtold>zero)dt=min(dt,dtold*1.1d0)
 dtold=dt
!
!
#ifdef VERBOSE
 print *,"#dt and Time",dt,time," at step ",step,f1
#endif
!
!
 if(dt+time>endtime)dt=endtime-time
 dt=dt*half
!$OMP END MASTER
!
!***
! Time step is now updated for all cases, and cut in half to do the 
! hydro halfstep for the second-order scheme.
!***
!
!$OMP BARRIER
 call source() ! source.f90, apply forces
!$OMP BARRIER
!
!
#ifndef NOHYDRO
!
!
#ifdef RADTRAN
 call mcrtfld_transfer ! mcrtfld.f90
!$OMP BARRIER
#endif
!
!
 call velocity() ! velocity.f90
!$OMP BARRIER
!
!***
! Conservative variables and primitives are now updated over the halfstep. 
! State is called at the end of sourcing and radiative transfer, so it is not
! included here. First, store the conservative variables at the half step.
!***
!
!$OMP DO SCHEDULE(STATIC)
 do igrid=1,ngrid
 do idim=1,5
   cons_new(idim,igrid)=cons(idim,igrid)
 enddo
 enddo
!$OMP ENDDO
!
#endif  /* end ifndef NOHYDRO */
!$OMP END PARALLEL
!
!
#ifndef NOHYDRO
!
!***
! Point to the flux variable that will be used as the base for the updates
! during fluxing, i.e., cons = flux(cons_pt). 
!***
!
 cons_pt=>cons

 call flux() ! flux.f90

!$OMP PARALLEL 
!
 call set_ghost_cells() ! utils.f90, user can set special boundary conditions.
!
!$OMP BARRIER
!
!***
! half fluxing complete.  Now check for changes in density field, and limit
! evolution accordingly. This helps with code stability when self-gravity
! is used. Perhaps this will be moved into a separate function later.
!***
!
!$OMP DO SCHEDULE(STATIC) REDUCTION(max:max_den_loc)
 do igrid=1,ngrid
  max_den_loc=max(max_den_loc,cons(1,igrid))
 enddo
!$OMP ENDDO
!
!$OMP MASTER
  min_den_loc=max_den_loc
!$OMP END MASTER
!
!$OMP BARRIER
!
!$OMP DO SCHEDULE(STATIC) REDUCTION(min:min_den_loc)
 do igrid=1,ngrid
  min_den_loc=min(min_den_loc,cons(1,igrid))
 enddo
!$OMP ENDDO
!
!$OMP MASTER
 f1_cycle=.false.
 max_den_change_loc=abs(one-max_den_loc/max_den_old_loc)
 max_den_change_loc=max(max_den_change_loc,abs(one-min_den_loc/min_den_old_loc))
!
 if(checkden.and.max_den_change_loc>allow)then
   f1=f1*half !allow/max_den_change_loc*f1*half
   print *, "# Density change of ",  max_den_change_loc, " detected. Allowing ", allow,".",f1
   f1_cycle=.true.
 else
   f1=min(f1*1.1d0,one)
 endif
!
!$OMP END MASTER
!
!$OMP BARRIER
!
 if(f1_cycle)then
!
!$OMP DO SCHEDULE(STATIC)
 do igrid=1,ngrid
  do idim=1,5
    cons(idim,igrid)=cons_old(idim,igrid)
  enddo
 enddo
!$OMP ENDDO
!
 call velocity() ! velocity.f90
!
 endif
!
!$OMP END PARALLEL
!
if(f1_cycle)cycle

!$OMP PARALLEL DEFAULT(SHARED) 
!
 call velocity() ! velocity.f90
!
!$OMP END PARALLEL
!
#endif  /* end ifndef NOHYDRO */
!
!
 dt=two*dt 
!
!***
! Now take a full step based on quantities at the half step.
!***
!
!
#ifndef NOHYDRO
!
!
 cons_pt=>cons_new
 call flux() ! flux.f90
!
!$OMP PARALLEL DEFAULT(SHARED) 
!$OMP DO SCHEDULE(STATIC)
  do igrid=1,ngrid
   do idim=1,5
    cons(idim,igrid)=cons_pt(idim,igrid)
   enddo
  enddo
!$OMP ENDDO
 call set_ghost_cells() ! utils.f90
!$OMP END PARALLEL
!
!
#endif /* end ifndef NOHYDRO */
!
!
!$OMP PARALLEL
!
!
#ifdef PARTICLE
!
!
 call drift_particles(dt) ! particle.f90
!
!
#ifdef NOHYDRO
 if(npart>0.and.use_pic)then
#endif /* A matching conditional will be used below for a FORTRAN endif */
!
!
#endif /* end ifdef PARTICLE */
!
!
#ifdef SELFGRAVITY
!
!
!$OMP DO SCHEDULE(STATIC)
 do igrid=1,ngrid
   rhotot(igrid)=cons(1,igrid)
 enddo
!$OMP ENDDO
!
!
#endif  /* end ifdef SELFGRAVITY */
!
!
#ifdef PARTICLE
!$OMP BARRIER
 call add_particle_density() ! particle.f90
!
!
#ifdef NOHYDRO
 endif ! do not bother with gravity if npart !>0 &&  if NOHYDRO defined
#endif /* Matching conditional for above FORTRAN if (npart>0.and.use_pic) */
!
!
#endif /* end ifdef PARTICLE */
!
!
!$OMP END PARALLEL
!
!
#ifdef SELFGRAVITY
!
!
#ifdef NOHYDRO
 if(npart>0.and.use_pic)then
#endif /* Conditional for skipping self-gravity if no hydro and no particle in cell. */
!
!
 call set_com_in_tree() ! gravity.f90
!
!
#ifdef FASTGRAVITY
 call get_pot_bc_from_tree() ! gravity.f90
 call vcycle_pot() ! gravity.f90
#else
 call get_pot_from_tree() ! gravity.f90
#endif /* endif FASTGRAVITY */
!
!
#ifdef USEPERT
 call rochepert() ! gravity.f90
#endif
!
!
#ifdef NOHYDRO
 endif
#endif  /* Matching conditional for skipping self-gravity */
!
!
#endif /* endif SELFGRAVITY */
!
!
#ifdef EXTERNALPHI
  call external_phi() ! gravity.f90
#endif
!
!***
! Full source and gravity updates complete. Final half source.
!***
!
 dt=half*dt
!
!
!$OMP PARALLEL DEFAULT(SHARED) 
#ifndef NOHYDRO
 call velocity() ! update from last flux.  Done after self-gravity.
 call state() ! state.f90
#endif /* end ifndef NOHYDRO */
!
!
 call source() ! source.f90
!
!
#ifndef NOHYDRO
!
!
#ifdef RADTRAN
 call mcrtfld_transfer() ! mcrtfld.f90
#endif
!
!
 call velocity() ! velocity.f90
!
!
#endif /* end ifndef NOHYDRO */
!
!
!$OMP END PARALLEL
!
!
 dt=dt*two
 time=time+dt
!
!
 if((time>=nextstop.or.time==endtime))then
   call write_files(step) ! write_files.f90
   nextstop=nextstop+dtout
 endif

enddo ! end of the major loop. Continue loop until final time or max steps is reached.
call cpu_time(timer_stop)
print *, "Elapsed Time is : ",timer_stop-timer_start
!
!#ifndef NOHYDRO
! call cleanup ! do not know why this is failing.
!#endif
!
end program
