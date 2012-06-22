subroutine avisc
 use parameters
 use derived_types
 use grid_commons
 use utils, only : get_boundary
 implicit none

 integer :: igrid,b(6)
 logical :: active

!$OMP DO SCHEDULE(STATIC)
 do igrid=1,ngrid
   call get_boundary(igrid,b,active)
   if(active)then
       qq(1,igrid)=avmagx*min(zero,u(1,b(3))-u(1,b(4)))**2
       qq(2,igrid)=avmagy*min(zero,u(2,b(1))-u(2,b(2)))**2
       qq(3,igrid)=avmagz*min(zero,u(3,b(5))-u(3,b(6)))**2
   else
     qq(1:3,igrid)=zero
   endif
 enddo
!$OMP ENDDO NOWAIT

end subroutine
