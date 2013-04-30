!
! Calculate courant time for hydrodynamics time stepping.
!
subroutine courant
 use parameters
 use derived_types
 use grid_commons
 implicit none

 integer::igrid
 real(pre)::r0,sound

!$OMP MASTER
 dtrate=zero
!$OMP END MASTER
!$OMP BARRIER
!$OMP DO SCHEDULE(STATIC) REDUCTION(max:dtrate)  &
!$OMP PRIVATE(r0,sound)
 do igrid=1,ngrid
   sound=sqrt(adindx(igrid)*p(igrid)/cons(1,igrid))
   r0=(sound+abs(u(1,igrid)))/dx+(sound+abs(u(2,igrid)))/dy+(sound+abs(u(3,igrid)))/dz
   dtrate=max(dtrate,r0)
 enddo
!$OMP ENDDO 
!$OMP BARRIER

!$OMP MASTER
! dtrate=max(dtrate,rate_expand)
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
  
