!
! Basic utility for reading celldump files.
!
subroutine read_hydro()
 use parameters
 use derived_types
 use grid_commons
 use eos
 implicit none

 integer:: igrid,ibd,step
 real(pre)::x,y,z,eps,tk,rho,hscale,ekin,ploc
 character::cindx*8,filename*80
 type(units)::scl_read
 
#ifdef DOUBLEGRID
 integer:: ngridin,nxin,nyin,nzin,jgrid,in,ix,iy
 real(pre),dimension(:),allocatable::pin,phiin,muc_arrayin
 real(pre),dimension(:,:),allocatable::consin,uin
#endif
 
!$OMP MASTER
 write(cindx,'(I8.8)')irestart
 filename="celldump."//cindx
 
#ifdef DOUBLEGRID
 nxin=1600
 nyin=1600
 nzin=1
 ngridin=nxin*nyin*nzin
 allocate(pin(ngridin))
 allocate(phiin(ngridin))
 allocate(muc_arrayin(ngridin))
 allocate(consin(2,ngridin))
 allocate(uin(3,ngridin))
#endif
 
!
 print *, "RESTART FILE ",filename
!
!***
! Output can be written in ascii form.
! This is good for very small runs and tests.
!***
!
#ifdef DOUBLEGRID
!
 open(unit=100,file=filename,form="UNFORMATTED")
 read(100)step
 read(100)time,step
 starttime=time
 read(100)scl_read
 igrid=1
 do while(igrid<=ngridin)
      read(100)x,y,z,&
         consin(1,igrid),pin(igrid),&
        uin(1,igrid),uin(2,igrid),uin(3,igrid),phiin(igrid),consin(2,igrid),muc_arrayin(igrid),ibd
      igrid=igrid+1
 enddo

 do iy = 1,nyin
   do ix = 1,nxin

     in = nxin*(iy-1) + ix

     jgrid = nx*( (iy-1)*2 ) + (ix-1)*2 + 1

     cons(1,jgrid)=consin(1,in)
     cons(5,jgrid)=consin(2,in)
     u(1,jgrid)=uin(1,in)
     u(2,jgrid)=uin(2,in)
     u(3,jgrid)=uin(3,in)
     phi(jgrid)=phiin(in)
     muc_array(jgrid)=muc_arrayin(in)
     p(jgrid)=pin(in)

     jgrid = nx*( (iy-1)*2 ) + 1 + (ix-1)*2 + 1

     cons(1,jgrid)=consin(1,in)
     cons(5,jgrid)=consin(2,in)
     u(1,jgrid)=uin(1,in)
     u(2,jgrid)=uin(2,in)
     u(3,jgrid)=uin(3,in)
     phi(jgrid)=phiin(in)
     muc_array(jgrid)=muc_arrayin(in)
     p(jgrid)=pin(in)


     jgrid = nx*( (iy)*2-1 ) + (ix-1)*2 +1

     cons(1,jgrid)=consin(1,in)
     cons(5,jgrid)=consin(2,in)
     u(1,jgrid)=uin(1,in)
     u(2,jgrid)=uin(2,in)
     u(3,jgrid)=uin(3,in)
     phi(jgrid)=phiin(in)
     muc_array(jgrid)=muc_arrayin(in)
     p(jgrid)=pin(in)


     jgrid = nx*( (iy)*2-1 ) + (ix-1)*2+1 +1

     cons(1,jgrid)=consin(1,in)
     cons(5,jgrid)=consin(2,in)
     u(1,jgrid)=uin(1,in)
     u(2,jgrid)=uin(2,in)
     u(3,jgrid)=uin(3,in)
     phi(jgrid)=phiin(in)
     muc_array(jgrid)=muc_arrayin(in)
     p(jgrid)=pin(in)

!     print* ,jgrid,in,ngrid,ngridin 
 enddo
 enddo

 deallocate(consin,phiin,pin,uin,muc_arrayin)

!
#else
!

!
 open(unit=100,file=filename,form="UNFORMATTED")
 read(100)step
 read(100)time,step
 starttime=time
 read(100)scl_read
 igrid=1
 do while(igrid<=ngrid)
      read(100)x,y,z,&
         cons(1,igrid),p(igrid),&
        u(1,igrid),u(2,igrid),u(3,igrid),phi(igrid),cons(5,igrid),muc_array(igrid),ibd
      igrid=igrid+1
 enddo
#endif
!
!
!
 close(100)
!
!$OMP END MASTER
!$OMP BARRIER
!
 call get_units(scl_read)
!
!***
! set your initial conditions here
!***
!

!$OMP DO SCHEDULE(STATIC) PRIVATE(x,y,z,rho,tk,hscale,eps,ekin,ploc) 
 do igrid=1,ngrid
!
!
  cons(2,igrid)=u(1,igrid)*cons(1,igrid)
  cons(3,igrid)=u(2,igrid)*cons(1,igrid)
  cons(4,igrid)=u(3,igrid)*cons(1,igrid)
  adindx(igrid)=gammafix ! temporary re-initialization.  

#ifdef RESETENG
  if(igrid==1)print *,"Resetting ENG in read"
  eps = p(igrid)/(gammafix-one)
  cons(5,igrid) = eps + half*cons(1,igrid)*( u(1,igrid)**2 + u(2,igrid)**2 + u(3,igrid)**2 )
#endif

!
!
 enddo
!$OMP ENDDO


!
!$OMP DO SCHEDULE(STATIC)
 do igrid=1,nbound
  cons(1,indx_bound(igrid))=small_rho
  cons(2,indx_bound(igrid))=zero 
  cons(3,indx_bound(igrid))=zero 
  cons(4,indx_bound(igrid))=zero 
  cons(5,indx_bound(igrid))=small_eps
 enddo
!$OMP ENDDO NOWAIT

 call state() ! state.f90

 print *, "Done with ICs"

end subroutine 
   
