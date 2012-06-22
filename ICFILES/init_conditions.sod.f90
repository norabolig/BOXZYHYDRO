subroutine init_conditions()
 use parameters
 use derived_types
 use grid_commons
 implicit none

 integer:: igrid
 
 type(units)::scale

 real(pre)::x,y,z
 
 call get_units(scale)
! set your initial conditions here
! these initial conditions are for a sod shock tube test
! directions x, y, and z can be chosen below where indicated 

!$OMP DO SCHEDULE(STATIC) PRIVATE(x,y,z)
 do igrid=1,ngrid
  x=grid(igrid)%x;y=grid(igrid)%y;z=grid(igrid)%z
  phi(igrid)=zero
  
  if(y>zero)then ! change flow direction here
   cons(1,igrid)=one
   p(igrid)=one
   adindx(igrid)=five/three
   cons(5,igrid)=p(igrid)/(adindx(igrid)-one)
   cons(2,igrid)=zero
   cons(3,igrid)=zero
   cons(4,igrid)=zero
   u(1,igrid)=zero
   u(2,igrid)=zero
   u(3,igrid)=zero
  else
   cons(1,igrid)=0.125
   p(igrid)=one/ten
   adindx(igrid)=five/three
   cons(5,igrid)=p(igrid)/(adindx(igrid)-one)
   cons(2,igrid)=zero
   cons(3,igrid)=zero
   cons(4,igrid)=zero
   u(1,igrid)=zero
   u(2,igrid)=zero
   u(3,igrid)=zero
  endif

 enddo
!$OMP ENDDO
!$OMP DO SCHEDULE(STATIC)
 do igrid=1,nbound
  cons(1,indx_bound(igrid))=small_rho
  cons(2,indx_bound(igrid))=zero 
  cons(3,indx_bound(igrid))=zero 
  cons(4,indx_bound(igrid))=zero 
  cons(5,indx_bound(igrid))=small_eps
 enddo
!$OMP ENDDO
!$OMP DO SCHEDULE(STATIC)
do igrid=1,ngrid
 cons_old(1,igrid)=cons(1,igrid)
 cons_old(2,igrid)=cons(2,igrid)
 cons_old(3,igrid)=cons(3,igrid)
 cons_old(4,igrid)=cons(4,igrid)
 cons_old(5,igrid)=cons(5,igrid)
enddo
!$OMP ENDDO
 print *, "#Done with ICs"

end subroutine 
   

