subroutine velocity
 use parameters
 use derived_types
 use grid_commons
 implicit none

 real(pre)::ekin,speed
 integer :: igrid,idim,flag

 maxv=zero
!$OMP DO SCHEDULE(STATIC) PRIVATE(ekin,flag)
 do igrid=1,ngrid
    do idim=1,3
       u(idim,igrid)=cons(idim+1,igrid)/cons(1,igrid)
    enddo
 enddo
!$OMP ENDDO
!$OMP DO SCHEDULE(STATIC) REDUCTION(max:maxv) PRIVATE(ekin,flag,speed)
 do igrid=1,ngrid
  flag=0
  do idim=1,3
   maxv=max(maxv,abs(u(idim,igrid)))
   if(abs(u(idim,igrid))>vlimit)flag=1
  enddo
  if(flag==1)then
     ekin=half*(cons(2,igrid)**2+cons(3,igrid)**2+cons(4,igrid)**2)/cons(1,igrid)
     speed=zero
     do idim=1,3
       speed=speed+u(idim,igrid)**2
     enddo
     speed=sqrt(speed)
     do idim=1,3
       u(idim,igrid)=u(idim,igrid)*vlimit/speed
       cons(idim+1,igrid)=u(idim,igrid)*cons(1,igrid)
     enddo
     cons(5,igrid)=cons(5,igrid)-ekin
     ekin=half*(cons(2,igrid)**2+cons(3,igrid)**2+cons(4,igrid)**2)/cons(1,igrid)
     cons(5,igrid)=max(cons(5,igrid),small_eps)+ekin
  endif
 enddo
!$OMP ENDDO
#ifdef VERBOSE
!$OMP MASTER
 print *, "# maxv: ",maxv
!$OMP END MASTER
#endif

end subroutine
