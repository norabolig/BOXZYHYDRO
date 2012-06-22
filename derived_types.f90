module derived_types
 use parameters 

 type gridcell
  integer:: id,ix,iy,iz
  integer,dimension(6):: ineigh
  logical::boundary,anchor
  real(pre)::x,y,z
 end type

 type units
  real(pre)::length,time,mass,density
  real(pre)::rgas,eps,kelvin,vel
 end type

 type bidx
   integer::i(6)
 end type

 type gravcell
  integer::id
  !integer,dimension(27)::ichild
  integer,dimension(8)::ichild
  real(pre)::x,y,z,xcm,ycm,zcm,mass,rmax
  real(pre)::rho,phi
  logical::boundary
!  real(pre)::dpole(3)
!  real(pre)::qpole(3,3)
 end type

 type particle_type
  integer::id
  real(pre)::x,y,z,vx,vy,vz,fx,fy,fz,soft,r,m,rho0
#ifdef THERMALHIST
  real(pre)::t,p,d,tm,pm,dm
#endif
  logical::active
 end type

 type flux_data
   real(pre)::f2,f1,f0,fm1,v,s2,s1,s0,sm1,ds,dt
 end type

end module

