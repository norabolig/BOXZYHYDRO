!
! Get input using the namelist file
!
module input
 use parameters
 implicit none
 namelist /io_control/write_ascii
 namelist /grid_params/nx,ny,nz,dx,dy,dz,yoffset,fluxangmom, &
   no_outflow_xl,no_outflow_xr,no_outflow_yl,no_outflow_yr, &
   no_outflow_zl,no_outflow_zr
 namelist /hydro/ small_rho,small_eps,vlimit,endtime,den_change_tol, &
                  avmagx,avmagy,avmagz,dtout,starttime,irestart,cfl
 namelist /gravity/ yml_max,nrad_yml,anchor_space,grav_err_tol,miniter,grav_err_tol_low
 namelist /eos_input/mu_z,brot,vib,dTk_eos,tk_bgrnd,ac,bc,gammafix,&
               xabun,yabun,zabun,NEOS_RHO,H2STAT
 namelist /particle_input/npart,npart_direct,tg_immediate_couple,a_sublimate_limit,use_pic,initialize_particles_now


 contains

 subroutine read_params(filename)
  character filename*80
  open(102,file=trim(filename))
  read(102,nml=io_control)
  rewind(102)
  read(102,nml=grid_params)
  rewind(102)
  read(102,nml=hydro)
  rewind(102)
  read(102,nml=gravity)
  rewind(102)
  read(102,nml=eos_input)
  rewind(102)
  read(102,nml=particle_input)
  close(102)
 end subroutine read_params

end module

