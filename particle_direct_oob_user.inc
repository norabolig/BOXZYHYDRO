!
! This is used to treat out of bounds.  Default is to remove the particle completely
! The file is supposed to be for direct particles
! See also particle_oob_user.inc
!
! BEGIN DEFAULT
       print *, "Removing Particle :", ipart,x/dx,y/dy,z/dz
       part_direct(ipart)%active=.false.
       cycle
! END DEFAULT
!
! BEGIN USER
