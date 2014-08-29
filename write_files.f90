!
! Write hydro and particle files, celldump and pdump, respectively.
!
!
subroutine write_files(step)
 use parameters
 use derived_types
 use grid_commons
!
!
#ifdef PARTICLE
 use particle
#endif
!
!
 implicit none
!
 integer, intent(in)::step
 integer:: igrid,ibd,readtest=1234,ipart
 real(pre)::x,y,z,eps,tk
 character::cindx*8,filename*80
 type(units)::scale
!
!
 call get_units(scale)
!
!
#ifndef NOHYDRO
!
!
!$OMP PARALLEL DEFAULT(SHARED)
  call state() ! state.f90, make sure pressure is fully updated.
!$OMP END PARALLEL
!
!
  write(cindx,'(I8.8)')step  
  filename="celldump."//cindx
!
!
#ifdef VERBOSE
  print *, filename
#endif
!
!
  if(write_ascii)then
     open(unit=100,file=filename)
     write(100,'(A25,1pe16.8E3,I8)')"#Time of snapshot at step:",time,step
     write(100,'(A,1X,8(1pe15.8,1X))')'#',scale
     do igrid=1,ngrid
        write(100,'(11(1pe16.8E3,1X),I2)'),grid(igrid)%x,grid(igrid)%y,grid(igrid)%z,&
           cons(1,igrid),p(igrid),&
           u(1,igrid),u(2,igrid),u(3,igrid),phi(igrid),cons(5,igrid),muc_array(igrid),grid(igrid)%boundary
     enddo
   else
     open(unit=100,file=filename,form='UNFORMATTED')
     write(100)readtest
     write(100)time,step
     write(100)scale
     do igrid=1,ngrid
        write(100)grid(igrid)%x,grid(igrid)%y,grid(igrid)%z,&
           cons(1,igrid),p(igrid), &
           u(1,igrid),u(2,igrid),u(3,igrid),phi(igrid),cons(5,igrid),muc_array(igrid),grid(igrid)%boundary
     enddo
   endif
   close(100)
!
!
#endif /* end ifndef NOHYDRO */
!
!
#ifdef PARTICLE
!
!
   filename=""
   write(cindx,'(I8.8)')step  
   filename="pdump."//cindx
   print *, filename
   open(unit=100,file=filename)
   write(100,'(1pe16.9,1X,3(I9))')time,step,npart,npart_direct
   if(npart_direct>0)then
     do ipart=1,npart_direct
!
!
#ifdef WITHDRAG
!
!
       write(100,'(I9,1X,I9,16(1X,1pe16.9),1X,L1)')1,part_direct(ipart)%id, &
         part_direct(ipart)%x,part_direct(ipart)%y, &
         part_direct(ipart)%z,part_direct(ipart)%vx,part_direct(ipart)%vy, &
         part_direct(ipart)%vz,part_direct(ipart)%m,part_direct(ipart)%soft, &
         part_direct(ipart)%rho0, &
         part_direct(ipart)%r, &
         part_direct(ipart)%d, &
         part_direct(ipart)%t, &
         part_direct(ipart)%p, &
         part_direct(ipart)%dm, &
         part_direct(ipart)%tm, &
         part_direct(ipart)%pm, &
         part_direct(ipart)%active
!
!
#else
!
!
       write(100,'(I9,1X,I9,8(1X,1pe16.9),1X,L1)')1,part_direct(ipart)%id, &
         part_direct(ipart)%x,part_direct(ipart)%y, &
         part_direct(ipart)%z,part_direct(ipart)%vx,part_direct(ipart)%vy, &
         part_direct(ipart)%vz,part_direct(ipart)%m,part_direct(ipart)%soft, &
         part_direct(ipart)%active
!
!
#endif /* end ifdef WITHDRAG */
!
!
     enddo
   endif
   if(npart>0)then
     do ipart=1,npart
!
!
#ifdef WITHDRAG
!
!
       write(100,'(I9,1X,I9,16(1X,1pe16.9),1X,L1)')0,part(ipart)%id,part(ipart)%x,&
         part(ipart)%y,part(ipart)%z,part(ipart)%vx,part(ipart)%vy,part(ipart)%vz, &
         part(ipart)%m,part(ipart)%soft,&
         part(ipart)%rho0, &
         part(ipart)%r, &
         part(ipart)%d, &
         part(ipart)%t, &
         part(ipart)%p, &
         part(ipart)%dm, &
         part(ipart)%tm, &
         part(ipart)%pm, &
         part(ipart)%active
!
!
#else
!
!
       write(100,'(I9,1X,I9,8(1X,1pe16.9),1X,L1)')0,part(ipart)%id,part(ipart)%x,&
         part(ipart)%y,part(ipart)%z,part(ipart)%vx,part(ipart)%vy,part(ipart)%vz, &
         part(ipart)%m,part(ipart)%soft,&
         part(ipart)%active
!
!
#endif /* end ifdef WITHDRAG */
!
!
     enddo
   endif
   close(100)
!
!
#endif /* end ifdef PARTICLE */
!
!
end subroutine 
