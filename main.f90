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
 real(pre)::nextstop,f1,timer_start,timer_stop
 real(pre)::dtold
 real(pre)::saback,safron,satop,sabot,pback,pfron,ptop,pbot,rhoav,vxav,avcount,chord=150d0
 real(pre)::etot_loc,mass_loc
 real(pre)::etot_loc0,mass_loc0,eloss_tot,ediff,factor,f_loc,f_loc0
 type(units)::scale
 character*8::cindx
 character*17::filename
 character*80::namelist_file
 integer::iter

 logical:: checkden=.false.,perturb=.false.

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

 etot_loc=zero;mass_loc=zero
 etot_loc0=zero;mass_loc0=zero
 f_loc=zero;f_loc0=zero
 eloss_tot=zero
!$OMP PARALLEL DEFAULT(SHARED) private(x,y,z,angle)                    
#ifndef NOHYDRO
  call conservation_diagnostic() ! in utils.f90 xcom,ycom,zcom defined in utils
  call get_etot(etot_loc0,mass_loc0,f_loc0)
!
!
#ifdef SUPPRESSDRIFT
  call suppress_drift(xcom,ycom,zcom) ! in utils.f90
#endif 
!
!
!
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
! if(step==1)dt=1d-3*dt ! use this with caution.  It is best for models that are not quite in equilibrium
!                       ! and need a kick-start for stability.
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
#ifndef NOHYDRO
!
!
#ifdef RADTRAN
 call mcrtfld_transfer ! mcrtfld.f90
!$OMP BARRIER
#endif
!$OMP MASTER
#ifdef MCFLDRT
 eloss_tot=eloss_tot+total_cool
#else
 eloss_tot=zero
#endif
!$OMP END MASTER
!
!
 call set_ghost_cells() ! utils.f90, user can set special boundary conditions.
!$OMP BARRIER
 call state()
!$OMP BARRIER
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
   cons_old(idim,igrid)=cons(idim,igrid)
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

 dt=dt*two
 call flux(0) ! flux.f90

!$OMP PARALLEL 
!
 call set_ghost_cells() ! utils.f90, user can set special boundary conditions.
!$OMP BARRIER
 call state()
!$OMP BARRIER
 call velocity()
!
!***
! first fluxing complete.  Now check for changes in density field, and limit
! evolution accordingly. This helps with code stability when self-gravity
! is used. Perhaps this will be moved into a separate function later.
!***
!
!$OMP END PARALLEL
!
 call flux(1)

!$OMP PARALLEL DEFAULT(SHARED) 
!
 call set_ghost_cells() ! utils.f90, user can set special boundary conditions.
!$OMP BARRIER
 call state()
!$OMP BARRIER
 call velocity()
!
!$OMP END PARALLEL
!
#endif  /* end ifndef NOHYDRO */
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
!$OMP MASTER
#ifdef MCFLDRT
 eloss_tot=eloss_tot+total_cool
#else
 eloss_tot=zero
#endif
!$OMP END MASTER
!
!
 call set_ghost_cells() ! utils.f90, user can set special boundary conditions.
!$OMP BARRIER
 call state()
!$OMP BARRIER
 call velocity() ! velocity.f90
!
#ifdef SELFGRAVITY
#ifdef ENERGYCOR
 call get_etot(etot_loc,mass_loc,f_loc)
!$OMP BARRIER
!$OMP MASTER
 print *, etot_loc0,etot_loc,mass_loc0,mass_loc,eloss_tot
 ediff=(etot_loc-etot_loc0-eloss_tot)/f_loc
!$OMP END MASTER
!$OMP BARRIER
!$OMP DO SCHEDULE(STATIC) PRIVATE(factor)
 do igrid=1,ngrid
    factor=cons(1,igrid)*abs(u(1,igrid)*gforce(1,igrid)+u(2,igrid)*gforce(2,igrid)+&
      u(3,igrid)*gforce(3,igrid))*mass_loc/mass_loc0
    cons(5,igrid)=max(cons(5,igrid)-ediff*factor,small_eps)
 enddo
!$OMP ENDDO
#endif /* end ifdef ENERGYCOR */
#endif /* end ifdef SELFGRAVITY */
!
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
#ifndef NOWRITEOUT
 if((time>=nextstop.or.time==endtime))then
   call write_files(step) ! write_files.f90
   nextstop=nextstop+dtout
 endif
#endif /* end ifndef NOWRITEOUT */

enddo ! end of the major loop. Continue loop until final time or max steps is reached.
call cpu_time(timer_stop)
print *, "Elapsed Time is : ",timer_stop-timer_start
!
!#ifndef NOHYDRO
! call cleanup ! do not know why this is failing.
!#endif
!
end program
