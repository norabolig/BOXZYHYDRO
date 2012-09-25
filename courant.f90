!
! Calculate courant time for hydrodynamics time stepping.
!
subroutine courant
 use parameters
 use derived_types
 use grid_commons
 implicit none

 integer::igrid
 real(pre)::sound2,s1,s2,s3,s4
 real(pre)::r1,r2,r0

!$OMP MASTER
 dtrate=zero
!$OMP END MASTER
!$OMP BARRIER
!$OMP DO SCHEDULE(STATIC) REDUCTION(max:dtrate)  &
!$OMP PRIVATE(s1,s2,s3,s4,r1,r2,r0,sound2)
 do igrid=1,ngrid
  sound2=(adindx(igrid)*p(igrid)/cons(1,igrid))/(dx*dx+dy*dy+dz*dz)
  s1=(u(1,igrid)/dx)**2
  s2=(u(2,igrid)/dy)**2
  s3=(u(3,igrid)/dz)**2
  s4=16d0*((qq(1,igrid)/dx)**2+(qq(2,igrid)/dy)**2+(qq(3,igrid)/dz)**2)
  r0=sqrt(s1+s2+s3+s4+sound2)
  dtrate=max(dtrate,r0)
 enddo
!$OMP ENDDO 
!$OMP BARRIER

!$OMP MASTER
 dtrate=max(dtrate,rate_expand)
 dt=cfl/dtrate!*min(f1*f2,one)
!
!
#ifdef VERBOSE
 print *, "max speeds ",dtrate
#endif
!
!
!$OMP END MASTER

end subroutine
  
