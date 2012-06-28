subroutine init_conditions()
 use parameters
 use derived_types
 use grid_commons
 use utils, only : set_ghost_cells
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
  
   cons(1,igrid)=rhoflow/scale%density
   p(igrid)=cons(1,igrid)*scale%rgas/muc*300
   adindx(igrid)=five/three
   u(1,igrid)=-vflow/scale%vel
   cons(2,igrid)=u(1,igrid)*cons(1,igrid)
   cons(3,igrid)=zero
   cons(4,igrid)=zero
   cons(5,igrid)=p(igrid)/(adindx(igrid)-one)+cons(2,igrid)**2/cons(1,igrid)*half
   u(2,igrid)=zero
   u(3,igrid)=zero

 enddo
!$OMP ENDDO

 call set_ghost_cells()

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
   

