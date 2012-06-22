subroutine read_hydro()
 use parameters
 use derived_types
 use grid_commons
 use eos
 implicit none

 integer:: igrid,ibd,step
 real(pre)::x,y,z,eps,tk
 character::junk*25,cindx*8,filename*80
 type(units)::scl_read
 
!$OMP MASTER
 write(cindx,'(I8.8)')irestart
 filename="celldump."//cindx
 print *, "#RESTART FILE ",filename
if(write_ascii)then
 open(unit=100,file=filename)
 read(100,'(A25,1pe16.8E3,I8)')junk,time,step
 starttime=time
 read(100,'(A1)')Junk
 igrid=1
#if OLDREAD
 do while(igrid<=ngrid)
      read(100,'(9(1pe16.8E3,1X),I2)'),x,y,z,&
         cons(1,igrid),p(igrid),&
         u(1,igrid),u(2,igrid),u(3,igrid),phi(igrid),ibd
      if (abs(z)>dz*dble(nz-1)*half)cycle
      igrid=igrid+1
 enddo
else
 open(unit=100,file=filename,form="UNFORMATTED")
 read(100)step
 read(100)time,step
 starttime=time
 read(100)scl_read
 igrid=1
 do while(igrid<=ngrid)
      read(100),x,y,z,&
         cons(1,igrid),p(igrid),&
         u(1,igrid),u(2,igrid),u(3,igrid),phi(igrid),ibd
      if (abs(z)>dz*dble(nz-1)*half)cycle
      igrid=igrid+1
 enddo
#else
 do while(igrid<=ngrid)
      read(100,'(9(1pe16.8E3,1X),I2)'),x,y,z,&
         cons(1,igrid),p(igrid),&
         u(1,igrid),u(2,igrid),u(3,igrid),phi(igrid),cons(5,igrid),muc_array(igrid),ibd
      if (abs(z)>dz*dble(nz-1)*half)cycle
      igrid=igrid+1
 enddo
else
 open(unit=100,file=filename,form="UNFORMATTED")
 read(100)step
 read(100)time,step
 starttime=time
 read(100)scl_read
 igrid=1
 do while(igrid<=ngrid)
      read(100),x,y,z,&
         cons(1,igrid),p(igrid),&
        u(1,igrid),u(2,igrid),u(3,igrid),phi(igrid),cons(5,igrid),muc_array(igrid),ibd
      if (abs(z)>dz*dble(nz-1)*half)cycle
      igrid=igrid+1
 enddo
#endif
endif
 close(100)
!$OMP END MASTER
!$OMP BARRIER
 call get_units(scl_read)

! set your initial conditions here

!$OMP DO SCHEDULE(STATIC) PRIVATE(x,y,z,tk,eps) 
 do igrid=1,ngrid
  !cons(5,igrid)=p(igrid)
#if OLDREAD
  tk=p(igrid)/scl%rgas/cons(1,igrid)*muc
  !if(tk>4d0)print *, igrid,tk,cons(1,igrid),eps
  call get_gamma_from_tk(eps,cons(1,igrid),tk,muc_array(igrid),adindx(igrid))
  !call get_gamma_from_p(eps,cons(1,igrid),p(igrid),muc_array(igrid),adindx(igrid))
   x=grid(igrid)%x;y=grid(igrid)%y;z=grid(igrid)%z

   if(x**2+y**2+z**2<dx**2+dy**2+dz**2)then
      print *, "IC Central ",cons(1,igrid),p(igrid), &
           muc_array(igrid),adindx(igrid),eps,&
           p(igrid)*muc_array(igrid)/scl_read%rgas/cons(1,igrid),scl_read%rgas
   endif
  !call get_gamma_from_p(eps,cons(1,igrid),p(igrid),muc_array(igrid),adindx(igrid))

  cons(2,igrid)=u(1,igrid)*cons(1,igrid)
  cons(3,igrid)=u(2,igrid)*cons(1,igrid)
  cons(4,igrid)=u(3,igrid)*cons(1,igrid)
  cons(5,igrid)=eps+half*cons(1,igrid)*(u(1,igrid)**2+u(2,igrid)**2+u(3,igrid)**2)
#else
  cons(2,igrid)=u(1,igrid)*cons(1,igrid)
  cons(3,igrid)=u(2,igrid)*cons(1,igrid)
  cons(4,igrid)=u(3,igrid)*cons(1,igrid)
#endif
 
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
!$OMP ENDDO NOWAIT

 call state()

 print *, "#Done with ICs"

end subroutine 
   
