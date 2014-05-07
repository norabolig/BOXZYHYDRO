!
! This is used to treat out of bounds.  
! Default is to remove the particle completely
! Recall that variables are defined in the actual
! .f90 source code (particle.f90).  
! 
! Be sure to clean before making, as included files
! are not generally re-made with a regular make.
!
! DEFAULT BEGIN
       print *, "Removing Particle :", ipart,x/dx,y/dy,z/dz
       part(ipart)%active=.false.
       cycle
! DEFAULT END
!
! USER BEGIN
!  EXAMPLE CODE
! This will reset particles to the +x boundary with a given -vx
!
!
!       part(ipart)%active=.true.
!       part(ipart)%x=dble(nx/2-2.001)*dx
!       call random_number(y)
!       part(ipart)%y=(one-two*y)*dble(ny/2-2)*dy
!       part(ipart)%z=zero
!
!       part(ipart)%vx=-vflow
!       part(ipart)%vy=zero
!       part(ipart)%vz=zero
