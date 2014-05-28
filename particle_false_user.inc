!
! This is used to treat non active particles.
! Recall that variables are defined in the actual
! .f90 source code (particle.f90).  
! 
! Be sure to clean before making, as included files
! are not generally re-made with a regular make.
!
! DEFAULT BEGIN
!  this simply keeps removed particles removed from the simulation

       cycle

! DEFAULT END
!
! USER BEGIN
!  EXAMPLE CODE
! This will reset particles to the +x boundary with a given -vx
! at a given rate, which is compared with a random variable between 0 and 1
!
! 
!         call random_number(x)
!         if(vflow/(dx*nx)*dt>x)then
!           part(ipart)%active=.true.
!           part(ipart)%x=(dble(nx/2)-2.001)*dx
!           call random_number(y)
!           part(ipart)%y=(one-two*y)*(dble(ny/2)-2.001)*dy
!           call random_number(z)
!           part(ipart)%z=(one-two*z)*(dble(nz/2)-2.001)*dz
!  
!           part(ipart)%vx=-vflow
!           part(ipart)%vy=zero
!           part(ipart)%vz=zero
!         endif
!         cycle

