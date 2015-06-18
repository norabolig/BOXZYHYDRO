!
! Basic run parameters and simulation control defaults
!
module parameters
 integer,parameter::pre=8

 real(pre)::gcgs=6.67d-8
 real(pre)::msun=1.989d33
 real(pre)::rgasCGS=8.254d7
 real(pre)::mp=1.67d-24
 real(pre)::pi=3.1415926535897931d0

 integer::irestart=0

!
!***
!io
!***
!
 logical::write_ascii=.true.
!
!***
!grid parameters
!***
!
 integer::nx=50,ny=50,nz=50
 real(pre)::dx=0.25d0
 real(pre)::dy=0.25d0
 real(pre)::dz=0.25d0
 real(pre)::yoffset=0.d0
 logical::fluxangmom=.false.
 logical::no_outflow_xl=.true.
 logical::no_outflow_xr=.true.
 logical::no_outflow_yl=.true.
 logical::no_outflow_yr=.true.
 logical::no_outflow_zl=.true.
 logical::no_outflow_zr=.true.
!
!***
!eos parameters
!***
!
 real(pre)::mu_z=16.78d0
 real(pre)::brot=85.4d0
 real(pre)::vib=5987d0
 real(pre)::diss=52000d0
 real(pre)::drho_eos=0d0
 real(pre)::dTk_eos=5d0
 real(pre)::deng_eos=0d0
 real(pre)::tk_bgrnd=1d-6
 real(pre)::ac=1d0 ! parahydrogen mixture
 real(pre)::bc=3d0 ! orthohydrogen
 real(pre)::gammafix=1.66666666667d0 ! used only for H2STAT=-1
 real(pre)::xabun=0.73d0,yabun=0.25d0,zabun=0.02d0
 real(pre)::rho_eos_high=1d-4
 real(pre)::rho_eos_low=1d-15
 real(pre)::tk_eos_cutoff=1d4
!
 integer::NEOS_T=600
 integer::NEOS_RHO=500
 integer::H2STAT=0 ! 0 mixture, 1 equilibrium, -1 fixed gamma
!
 logical::extend_table=.true.
!
!***
!hydro
!***
!
 real(pre)::small_rho=1d-10
 real(pre)::small_eps=1d-40
 real(pre)::vlimit=1d1
 real(pre)::starttime=0d0
 real(pre)::endtime=0d0
 real(pre)::dtout=10d0
 real(pre)::cfl=0.5d0
!
 real(pre)::den_change_tol=0.1d0
 real(pre)::avmagx=0.0d0,avmagy=0.d0,avmagz=0.00d0
!
 real(pre)::anchorradius=100d0
 real(pre)::tkflow=300d0
 real(pre)::rhoflow=1d-9
 real(pre)::vflow=9e5
 real(pre)::object_radius=3.1
 real(pre)::object_x_displace=-15d0
 real(pre)::object_mass=1d0
 real(pre)::x0dot=0d0
 real(pre)::pmass_factor=3.375d-5
 real(pre)::psize1=0.1d0
 real(pre)::psize2=0.1d0
 real(pre)::turbulence=0.01d0


 character(32)::flux_limit_type='minmod'
#define MINMOD
!
!***
! Above once was used as an on-the-fly switch, but it is a HUGE performance hit.
! better to just define it using the preprocessor
!***
!
!
!***
!gravity
!***
 integer::yml_max=10
 integer::nrad_yml=100
 integer::anchor_space=10
 integer::anchor_number
 integer::miniter=10
 real(pre)::grav_err_tol=1e-5
 real(pre)::grav_err_tol_low=1e-5
!
!***
!particles
!***
!
 integer::npart=100000
 integer::npart_direct=0
 integer::ntrace=1
 integer::nstep_print_part=10
 real(pre)::tg_immediate_couple=1450d0
 real(pre)::a_sublimate_limit=6.67e-20
 logical::use_pic=.true.
 logical::initialize_particles_now=.true.
!
!***
! constants
!***
!
 real(pre),parameter:: &
  zero=0d0,&
  one=1.d0,&
  two=2d0,&
  three=3d0,&
  four=4d0,&
  five=5d0,&
  six=6d0,&
  seven=7d0,&
  eight=8d0,&
  nine=9d0,&
  ten=10d0,&
  half=one/two,&
  quarter=one/four
  
end module
 
