    if(npart>0)then
      do ipart=1,npart
         part(ipart)%active=.true.
!
!***
! The following is special use, and should be setup according to the user.
! Will make this a separate file.
!***

    r=zero
    do while (r < object_radius)
        call random_number(x)
        part(ipart)%x=dble(nx/2-2)*dx*(one-two*x)
        call random_number(y)
        part(ipart)%y=dble(ny/2-2)*dy*(one-two*y)
        part(ipart)%z=0
        r=sqrt(part(ipart)%x**2+part(ipart)%y**2)
    enddo
       part(ipart)%m=rhoflow*(dble(nx)*dble(ny)*dble(nz))*dx*dy*dz*pmass_factor/dble(NPART)
       part(ipart)%vy=zero
       part(ipart)%vx=-vflow
       part(ipart)%vz=zero
       part(ipart)%fx=zero
       part(ipart)%fy=zero
       part(ipart)%fz=zero
       part(ipart)%soft=zero
!
!
#ifdef WITHDRAG
       part(ipart)%rho0=3./scl%density
       asize=0.03/scl%length
       part(ipart)%r=asize
       part(ipart)%d=zero
       part(ipart)%t=zero
       part(ipart)%p=zero
       part(ipart)%dm=zero
       part(ipart)%tm=zero
       part(ipart)%pm=zero
#endif /* end ifdef WITHDRAG */
!
!
      part(ipart)%id=id
      id=id+1
    enddo
  endif 

