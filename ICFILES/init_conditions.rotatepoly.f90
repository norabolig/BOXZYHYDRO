subroutine init_conditions()
 use parameters
 use derived_types
 use grid_commons
 use eos
 implicit none

 integer:: igrid,idim
 
 type(units)::scale

 integer,parameter::ntble=62
 integer::itble,jtble,iter
 real(pre)::x,y,z,r,zi,ekin,angle
 real(pre)::eps,rho,tk,gam,eng,rpoly,val0,val1,jn
 real(pre)::rhotble(ntble,ntble),jntble(ntble,ntble),dum1,dum2,rtble(ntble)
 real(pre)::polyr=16.,kpoly=29d15,rhopoly=4.4e-10,jnpoly=3.10e18

!$OMP CRITICAL
 open(unit=3,file="init.hachisu",form="FORMATTED")
 read(3,"((1pe15.8))")dum1
 do itble=1,ntble
  rtble(itble)=dble(itble-1)*polyr/dble(ntble)
  do jtble=1,ntble
   iter=iter+1
   read(3,"(4(1pe15.8,1X))")dum1,dum2,rhotble(itble,jtble),jntble(itble,jtble)
  enddo
 enddo
 close(3)
!$OMP END CRITICAL

 
 call get_units(scale)
! set your initial conditions here

!$OMP DO SCHEDULE(STATIC) PRIVATE(x,y,z) 
 do igrid=1,ngrid
  x=grid(igrid)%x;y=grid(igrid)%y;z=abs(grid(igrid)%z)
  phi(igrid)=zero
  adindx(igrid)=gammafix
  muc_array(igrid)=muc
  r=sqrt(x*x+y*y)
  angle=atan2(y,x)
  rpoly=zero
  itble=1
  do while(rpoly<r)
    itble=itble+1
    if(itble>ntble-1)then
      itble=ntble+1
      exit
    endif
    rpoly=rtble(itble)
  enddo
  if(itble>1)itble=itble-1
  rpoly=zero
  jtble=1
  do while(rpoly<(z))
    jtble=jtble+1
    if(jtble>ntble-1)then
      jtble=ntble+1
      exit
     endif
    rpoly=rtble(jtble)
    !print *, z,rpoly
  enddo
  if(jtble>1)jtble=jtble-1
  if((itble<ntble.and.jtble<ntble))then
  !print *,"START", r,z,r-rtble(itble),rtble(itble+1),rhotble(itble,jtble),itble,jtble
  if ((rtble(itble+1)-rtble(itble))==zero)stop
   val0 = rhotble(itble,jtble)+(rhotble(itble+1,jtble)-rhotble(itble,jtble))*(r-rtble(itble))&
        / (rtble(itble+1)-rtble(itble))
   val1 = rhotble(itble,jtble+1)+(rhotble(itble+1,jtble+1)-rhotble(itble,jtble+1))*(r-rtble(itble))&
        / (rtble(itble+1)-rtble(itble))
   cons(1,igrid)=(val0+(val1-val0)*(z-rtble(jtble))/(rtble(jtble+1)-rtble(jtble)))*rhopoly/scale%density
   cons(1,igrid)=max(cons(1,igrid),small_rho)*(1.+0.01*cos(2.*angle))

   val0 = jntble(itble,jtble)+(jntble(itble+1,jtble)-jntble(itble,jtble))*(r-rtble(itble))&
        / (rtble(itble+1)-rtble(itble))
   val1 = jntble(itble,jtble+1)+(jntble(itble+1,jtble+1)-jntble(itble,jtble+1))*(r-rtble(itble))&
        / (rtble(itble+1)-rtble(itble))
   jn=(val0+(val1-val0)*(z-rtble(jtble))/(rtble(jtble+1)-rtble(jtble)))*jnpoly/(scale%length**2/scale%time)
   jn=max(jn,zero)

   if(r>zero)then
    u(1,igrid)=-jn*y/r**2
    u(2,igrid)=jn*x/r**2
   else
    u(1:2,igrid)=zero
   endif
   u(3,igrid)=zero
   !print *, r,rtble(itble+1),rtble(itble)
  else
   cons(1,igrid)=small_rho
   u(1:3,igrid)=zero
  endif
  ekin=half*(u(1,igrid)**2+u(2,igrid)**2+u(3,igrid)**2)*cons(1,igrid)
  p(igrid)=kpoly*(cons(1,igrid)*scale%density)**adindx(igrid)/scale%eps
  cons(5,igrid)=max(p(igrid)/(adindx(igrid)-one),small_eps)+ekin
  cons(2,igrid)=cons(1,igrid)*u(1,igrid)
  cons(3,igrid)=cons(1,igrid)*u(2,igrid)
  cons(4,igrid)=cons(1,igrid)*u(3,igrid)
  do idim=1,5
    cons_old(idim,igrid)=cons(idim,igrid)
  enddo

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

 print *, "#Done with ICs"

end subroutine 
   
