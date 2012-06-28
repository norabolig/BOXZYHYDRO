subroutine courant
 use parameters
 use derived_types
 use grid_commons
 implicit none

 integer::igrid
 real(pre)::sound2,s1,s2,s3,s4
 real(pre)::allow,max_den_change,f1,f2,r1,r2,r0

 allow=0.01d0
!$OMP MASTER
 dtrate=zero
 max_den=zero
!$OMP END MASTER
!$OMP BARRIER
!$OMP DO SCHEDULE(STATIC) REDUCTION(max:max_den,dtrate)  &
!$OMP PRIVATE(s1,s2,s3,s4,r1,r2,r0,sound2)
 do igrid=1,ngrid
  sound2=(adindx(igrid)*p(igrid)/cons(1,igrid))/(dx*dx+dy*dy+dz*dz)

  s1=(u(1,igrid)/dx)**2
  s2=(u(2,igrid)/dy)**2
  s3=(u(3,igrid)/dz)**2
  s4=16d0*((qq(1,igrid)/dx)**2+(qq(2,igrid)/dy)**2+(qq(3,igrid)/dz)**2)

  if(cons(1,igrid)>max_den) max_den=cons(1,igrid)

  r0=sqrt(s1+s2+s3+s4+sound2)
  dtrate=max(dtrate,r0)
 enddo
!$OMP ENDDO 
!$OMP BARRIER

!$OMP MASTER
 max_den_change=abs(one-max_den_old/max_den)
 if(max_den_old>zero)then
   f1=allow/(max(allow,max_den_change))
   f2=allow/(max(allow,max_den_change_old))
   max_den_change_old=max_den_change
   max_den_old=max_den
 else
   f1=one;f2=one
   max_den_old=max_den
   max_den_change_old=one
 endif
 dtrate=max(dtrate,rate_expand)
 dt=cfl/dtrate!*min(f1*f2,one)
#ifdef VERBOSE
 print *, "#max speeds ",dtrate,f1,f2
#endif
!$OMP END MASTER

end subroutine
  
