  if(npart_direct>0)then
    momy=zero
    momx=zero
    mass=zero
    do ipart=1,npart_direct
!
!***
! The following is special use. Please redefine as necessary.
! WIll eventually go into its own file.
!***
!
       part_direct(ipart)%active=.true.
       part_direct(ipart)%x=rand()
       part_direct(ipart)%y=rand()
       part_direct(ipart)%z=rand()
       part_direct(ipart)%m=rand()
       part_direct(ipart)%vy=zero
       part_direct(ipart)%vx=zero
       part_direct(ipart)%vz=zero
       part_direct(ipart)%fx=zero
       part_direct(ipart)%fy=zero
       part_direct(ipart)%fz=zero
       part_direct(ipart)%soft=0.1d0

       momx=momx+part_direct(ipart)%vx*part_direct(ipart)%m
       momy=momy+part_direct(ipart)%vy*part_direct(ipart)%m
       mass=mass+part_direct(ipart)%m
       part_direct(ipart)%id=id
       id=id+1
     enddo
     do ipart=1,npart_direct
       part_direct(ipart)%vx=part_direct(ipart)%vx-momx/mass
       part_direct(ipart)%vy=part_direct(ipart)%vy-momy/mass
     enddo
  endif
