module drag
 use parameters
 use derived_types
 use grid_commons
 use utils
 use particles
 implicit none

 type(units)::scl
 real(pre)::prate
 integer::nact1,nact2

 contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine get_weights_8(ipart,x,y,z,w)
  real(pre)::intent(in)::x,y,z
  real(pre)::intent(out)::w(8)
  real(pre)::delx1,delx2,delx3,delx4,dely1,dely2,dely3,dely4
  integer,intent(in)::ipart
  integer::ig(8)
    call get_nearest_8(x,y,z,ig) ! 1-3 are low
    if(ig(1)<1.or.ig(1)>ngrid.or. &
       grid(ig(1))%boundary>0.or. &
       grid(ig(2))%boundary>0.or. &
       grid(ig(3))%boundary>0.or. &
       grid(ig(4))%boundary>0.or. &
       grid(ig(5))%boundary>0.or. & 
       grid(ig(6))%boundary>0.or. &
       grid(ig(7))%boundary>0.or. &
       grid(ig(8))%boundary>0)then
       print *, "Removing Particle :", ipart,ig(1),x,y,z
       part(ipart)%active=.false.
       cycle
    endif

    delx1=abs(grid(ig(1))%x-x)
    dely1=abs(grid(ig(1))%y-y)
    delx2=abs(grid(ig(2))%x-x)
    dely2=abs(grid(ig(2))%y-y)
    delx3=abs(grid(ig(3))%x-x)
    dely3=abs(grid(ig(3))%y-y)
    delx4=abs(grid(ig(4))%x-x)
    dely4=abs(grid(ig(4))%y-y)

    area3=delx1*dely1
    area4=delx2*dely2
    area1=delx3*dely3
    area2=delx4*dely4
    
    w(1)=area1*(grid(ig(5))%z-z)
    w(2)=area2*(grid(ig(5))%z-z)
    w(3)=area3*(grid(ig(5))%z-z)
    w(4)=area4*(grid(ig(5))%z-z)
    w(5)=area1*(z-grid(ig(1))%z)
    w(6)=area2*(z-grid(ig(1))%z)
    w(7)=area3*(z-grid(ig(1))%z)
    w(8)=area4*(z-grid(ig(1))%z)

    weight=zero
    do iter=1,8
     weight=weight+w(iter)
    enddo
    do iter=1,8
     w(iter)=w(iter)/weight
    enddo
 
 end subroutine

 subroutine add_particle_density()
  integer::ig(8)
  real(pre)::x,y,z,dist1,d,vol
  real(pre)::delx1,delx2,delx3,delx4,dely1,dely2,dely3,dely4
  real(pre)::w1,w2,w3,w4,w5,w6,w7,w8,weight,azimuth,a,b,c,sinb

  vol=dx*dy*dz

!$OMP DO SCHEDULE(STATIC) 
  do ipart=1,npart
    if(.not.part(ipart)%active)cycle
    x=part(ipart)%x
    y=part(ipart)%y
    z=part(ipart)%z
    d=part(ipart)%m/vol

    call get_nearest_8(x,y,z,ig) ! 1-3 are low
    if(ig(1)<1.or.ig(1)>ngrid.or. &
       grid(ig(1))%boundary>0.or. &
       grid(ig(2))%boundary>0.or. &
       grid(ig(3))%boundary>0.or. &
       grid(ig(4))%boundary>0.or. &
       grid(ig(5))%boundary>0.or. & 
       grid(ig(6))%boundary>0.or. &
       grid(ig(7))%boundary>0.or. &
       grid(ig(8))%boundary>0)then
       print *, "Removing Particle :", ipart,ig(1),x,y,z
       part(ipart)%active=.false.
       cycle
    endif

    delx1=abs(grid(ig(1))%x-x)
    dely1=abs(grid(ig(1))%y-y)
    delx2=abs(grid(ig(2))%x-x)
    dely2=abs(grid(ig(2))%y-y)
    delx3=abs(grid(ig(3))%x-x)
    dely3=abs(grid(ig(3))%y-y)
    delx4=abs(grid(ig(4))%x-x)
    dely4=abs(grid(ig(4))%y-y)

    area3=delx1*dely1
    area4=delx2*dely2
    area1=delx3*dely3
    area2=delx4*dely4
 
    if (abs((area1+area2+area3+area4)-(dx*dy))/(dx*dy)>1d-6) then
       print *, x,y,z
       print *, grid(ig(1))%x,grid(ig(1))%y,grid(ig(1))%z
       print *, grid(ig(2))%x,grid(ig(2))%y,grid(ig(2))%z
       print *, grid(ig(3))%x,grid(ig(3))%y,grid(ig(3))%z
       print *, grid(ig(4))%x,grid(ig(4))%y,grid(ig(4))%z
       print *, "Areas do not add up PD",area1+area2+area3+area4,dx*dy
       print *, area1,area2,area3,area4
       stop
    endif


    w1=area1*(grid(ig(5))%z-z)
    w2=area2*(grid(ig(5))%z-z)
    w3=area3*(grid(ig(5))%z-z)
    w4=area4*(grid(ig(5))%z-z)
    w5=area1*(z-grid(ig(1))%z)
    w6=area2*(z-grid(ig(1))%z)
    w7=area3*(z-grid(ig(1))%z)
    w8=area4*(z-grid(ig(1))%z)

    weight=(w1+w2+w3+w4+w5+w6+w7+w8)
    w1=w1/weight
    w2=w2/weight
    w3=w3/weight
    w4=w4/weight
    w5=w5/weight
    w6=w6/weight
    w7=w7/weight
    w8=w8/weight
 
!$OMP ATOMIC
      rhotot(ig(1))=rhotot(ig(1))+d*w1

!$OMP ATOMIC
      rhotot(ig(2))=rhotot(ig(2))+d*w2

!$OMP ATOMIC
      rhotot(ig(3))=rhotot(ig(3))+d*w3

!$OMP ATOMIC
      rhotot(ig(4))=rhotot(ig(4))+d*w4

!$OMP ATOMIC
      rhotot(ig(5))=rhotot(ig(5))+d*w5

!$OMP ATOMIC
      rhotot(ig(6))=rhotot(ig(6))+d*w6

!$OMP ATOMIC
      rhotot(ig(7))=rhotot(ig(7))+d*w7
 
!$OMP ATOMIC
      rhotot(ig(8))=rhotot(ig(8))+d*w8
 
  enddo
!$OMP ENDDO  

 end subroutine

 subroutine drift_particles(t)
  integer::ipart
  real(pre),intent(in)::t
  real(pre)::x,y,z
!$OMP MASTER
  nact1=0
  nact2=0
!$OMP ENDMASTER

  if(npart>0)then
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:nact1)
  do ipart=1,npart
    if(.not.part(ipart)%active)cycle
    x=part(ipart)%x+part(ipart)%vx*t
    y=part(ipart)%y+part(ipart)%vy*t
    z=part(ipart)%z+part(ipart)%vz*t
    part(ipart)%x=x
    part(ipart)%y=y
    part(ipart)%z=z
    nact1=nact1+1
  enddo
!$OMP ENDDO
  endif
  if(npart_direct>0)then
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:nact2)
  do ipart=1,npart_direct
    if(.not.part_direct(ipart)%active)cycle
    x=part_direct(ipart)%x+part_direct(ipart)%vx*t
    y=part_direct(ipart)%y+part_direct(ipart)%vy*t
    z=part_direct(ipart)%z+part_direct(ipart)%vz*t
    part_direct(ipart)%x=x
    part_direct(ipart)%y=y
    part_direct(ipart)%z=z
    nact2=nact2+1
  enddo
!$OMP ENDDO
  endif

!$OMP MASTER
  print *, "Direct active particles = ",nact2
  print *, "Cloud  active particles = ",nact1
!$OMP END MASTER
 
 end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine kick_particles(t)
  integer::ipart,b(6)
  real(pre),intent(in)::t
  real(pre)::x,y,z,vx,vy,vz,fx,fy,fz
  real(pre)::dist1,sinb,weight
  real(pre)::delx1,delx2,delx3,delx4,dely1,dely2,dely3,dely4
  real(pre)::side1,side2,area1,area2,area3,area4
  real(pre)::w1,w2,w3,w4,w5,w6,w7,w8,a,b,c,azimuth

  if(npart>0)then
  call particle_direct_interp()
!$OMP BARRIER
!$OMP DO SCHEDULE(STATIC) 
  do ipart=1,npart
    if(.not.part(ipart)%active)cycle
    x=part(ipart)%x
    y=part(ipart)%y
    z=part(ipart)%z
    vx=part(ipart)%vx
    vy=part(ipart)%vy
    vz=part(ipart)%vz
    fx=part(ipart)%fx
    fy=part(ipart)%fy
    fz=part(ipart)%fz

    call get_nearest_8(x,y,z,ig) ! 1-3 are low
    if(ig(1)<1.or.ig(1)>ngrid.or. &
       grid(ig(1))%boundary>0.or. &
       grid(ig(2))%boundary>0.or. &
       grid(ig(3))%boundary>0.or. &
       grid(ig(4))%boundary>0.or. &
       grid(ig(5))%boundary>0.or. & 
       grid(ig(6))%boundary>0.or. &
       grid(ig(7))%boundary>0.or. &
       grid(ig(8))%boundary>0)then
       print *, "Removing Particle :", ipart,ig(1),x,y,z
       part(ipart)%active=.false.
       cycle
    endif
    delx1=abs(grid(ig(1))%x-x)
    dely1=abs(grid(ig(1))%y-y)
    delx2=abs(grid(ig(2))%x-x)
    dely2=abs(grid(ig(2))%y-y)
    delx3=abs(grid(ig(3))%x-x)
    dely3=abs(grid(ig(3))%y-y)
    delx4=abs(grid(ig(4))%x-x)
    dely4=abs(grid(ig(4))%y-y)

    area3=delx1*dely1
    area4=delx2*dely2
    area1=delx3*dely3
    area2=delx4*dely4

    if (abs((area1+area2+area3+area4)-(dx*dy))/(dx*dy)>1d-6) then
       print *, x,y,z
       print *, grid(ig(1))%x,grid(ig(1))%y,grid(ig(1))%z
       print *, grid(ig(2))%x,grid(ig(2))%y,grid(ig(2))%z
       print *, grid(ig(3))%x,grid(ig(3))%y,grid(ig(3))%z
       print *, grid(ig(4))%x,grid(ig(4))%y,grid(ig(4))%z
       print *, "Areas do not add up PD",area1+area2+area3+area4,dx*dy
       print *, area1,area2,area3,area4
       stop
    endif


    w1=area1*(grid(ig(5))%z-z)
    w2=area2*(grid(ig(5))%z-z)
    w3=area3*(grid(ig(5))%z-z)
    w4=area4*(grid(ig(5))%z-z)
    w5=area1*(z-grid(ig(1))%z)
    w6=area2*(z-grid(ig(1))%z)
    w7=area3*(z-grid(ig(1))%z)
    w8=area4*(z-grid(ig(1))%z)

    weight=(w1+w2+w3+w4+w5+w6+w7+w8)
    w1=w1/weight
    w2=w2/weight
    w3=w3/weight
    w4=w4/weight
    w5=w5/weight
    w6=w6/weight
    w7=w7/weight
    w8=w8/weight
    
    fx=fx+gforce(1,ig(1))*w1
    fx=fx+gforce(1,ig(2))*w2
    fx=fx+gforce(1,ig(3))*w3
    fx=fx+gforce(1,ig(4))*w4
    fx=fx+gforce(1,ig(5))*w5
    fx=fx+gforce(1,ig(6))*w6
    fx=fx+gforce(1,ig(7))*w7
    fx=fx+gforce(1,ig(8))*w8

    fy=fy+gforce(2,ig(1))*w1
    fy=fy+gforce(2,ig(2))*w2
    fy=fy+gforce(2,ig(3))*w3
    fy=fy+gforce(2,ig(4))*w4
    fy=fy+gforce(2,ig(5))*w5
    fy=fy+gforce(2,ig(6))*w6
    fy=fy+gforce(2,ig(7))*w7
    fy=fy+gforce(2,ig(8))*w8
 
    fz=fz+gforce(3,ig(1))*w1
    fz=fz+gforce(3,ig(2))*w2
    fz=fz+gforce(3,ig(3))*w3
    fz=fz+gforce(3,ig(4))*w4
    fz=fz+gforce(3,ig(5))*w5
    fz=fz+gforce(3,ig(6))*w6
    fz=fz+gforce(3,ig(7))*w7
    fz=fz+gforce(3,ig(8))*w8
#if 1==0
    PRINT *,"CHECK FORCES",ipart,time,x,y,z,fx,fy,fz
    PRINT *,"CHECK WEIGHTS",ipart,time,w1,w2,w3,w4,w5,w6,w7,w8,weight
    PRINT *,"CHECK AREA ",area1,area2,area3,area4,dist1,side1,side2,ds
    PRINT *,"CHECK angles ",a,b,c,a+b+c,sinb
    print *, "GRID CHECK"
    print *, x,y,z
    print *, grid(ig(1))%x,grid(ig(1))%y,grid(ig(1))%z
    print *, grid(ig(2))%x,grid(ig(2))%y,grid(ig(2))%z
    print *, grid(ig(3))%x,grid(ig(3))%y,grid(ig(3))%z
    print *, grid(ig(4))%x,grid(ig(4))%y,grid(ig(4))%z
    print *, grid(ig(5))%x,grid(ig(5))%y,grid(ig(5))%z
    print *, grid(ig(6))%x,grid(ig(6))%y,grid(ig(6))%z
    print *, grid(ig(7))%x,grid(ig(7))%y,grid(ig(7))%z
    print *, grid(ig(8))%x,grid(ig(8))%y,grid(ig(8))%z
#endif
  
    vx=vx+fx*t
    vy=vy+fy*t
    vz=vz+fz*t

    part(ipart)%vy=vy
    part(ipart)%vx=vx
    part(ipart)%vz=vz

  enddo
!$OMP ENDDO
  endif
  if(npart_direct>0)then
  call particle_direct_direct()
!$OMP BARRIER
#ifdef NOHYDRO
  if(.not.npart>0)then
!$OMP DO SCHEDULE(STATIC)
  do ipart=1,npart_direct
    if(.not.part_direct(ipart)%active)cycle
    fy=part_direct(ipart)%fy
    fx=part_direct(ipart)%fx
    fz=part_direct(ipart)%fz
    vy=part_direct(ipart)%vy
    vx=part_direct(ipart)%vx
    vz=part_direct(ipart)%vz
 
    vx=vx+fx*t
    vy=vy+fy*t
    vz=vz+fz*t

    part_direct(ipart)%vy=vy
    part_direct(ipart)%vx=vx
    part_direct(ipart)%vz=vz
  enddo
!$OMP ENDDO
  else
#endif
!$OMP DO SCHEDULE(STATIC)
  do ipart=1,npart_direct
    if(.not.part_direct(ipart)%active)cycle
    x=part_direct(ipart)%x
    y=part_direct(ipart)%y
    z=part_direct(ipart)%z
    vy=part_direct(ipart)%vy
    vx=part_direct(ipart)%vx
    vz=part_direct(ipart)%vz
    fy=part_direct(ipart)%fy
    fx=part_direct(ipart)%fx
    fz=part_direct(ipart)%fz

    call get_nearest_8(x,y,z,ig) ! 1-3 are low
    if(ig(1)<1.or.ig(1)>ngrid.or. &
       grid(ig(1))%boundary>0.or. &
       grid(ig(2))%boundary>0.or. &
       grid(ig(3))%boundary>0.or. &
       grid(ig(4))%boundary>0.or. &
       grid(ig(5))%boundary>0.or. & 
       grid(ig(6))%boundary>0.or. &
       grid(ig(7))%boundary>0.or. &
       grid(ig(8))%boundary>0)then
       print *, "Removing Particle :", ipart,ig(1),x,y,z
       part_direct(ipart)%active=.false.
       cycle
    endif
    delx1=abs(grid(ig(1))%x-x)
    dely1=abs(grid(ig(1))%y-y)
    delx2=abs(grid(ig(2))%x-x)
    dely2=abs(grid(ig(2))%y-y)
    delx3=abs(grid(ig(3))%x-x)
    dely3=abs(grid(ig(3))%y-y)
    delx4=abs(grid(ig(4))%x-x)
    dely4=abs(grid(ig(4))%y-y)

    area3=delx1*dely1
    area4=delx2*dely2
    area1=delx3*dely3
    area2=delx4*dely4
 

    if (abs((area1+area2+area3+area4)-(dx*dy))/(dx*dy)>1d-6) then
       print *, x,y,z
       print *, grid(ig(1))%x,grid(ig(1))%y,grid(ig(1))%z
       print *, grid(ig(2))%x,grid(ig(2))%y,grid(ig(2))%z
       print *, grid(ig(3))%x,grid(ig(3))%y,grid(ig(3))%z
       print *, grid(ig(4))%x,grid(ig(4))%y,grid(ig(4))%z
       print *, "Areas do not add up PD",area1+area2+area3+area4,dx*dy
       print *, area1,area2,area3,area4
       stop
    endif

    w1=area1*(grid(ig(5))%z-z)
    w2=area2*(grid(ig(5))%z-z)
    w3=area3*(grid(ig(5))%z-z)
    w4=area4*(grid(ig(5))%z-z)
    w5=area1*(z-grid(ig(1))%z)
    w6=area2*(z-grid(ig(1))%z)
    w7=area3*(z-grid(ig(1))%z)
    w8=area4*(z-grid(ig(1))%z)

    weight=(w1+w2+w3+w4+w5+w6+w7+w8)
    w1=w1/weight
    w2=w2/weight
    w3=w3/weight
    w4=w4/weight
    w5=w5/weight
    w6=w6/weight
    w7=w7/weight
    w8=w8/weight
   
    fx=fx+gforce(1,ig(1))*w1
    fx=fx+gforce(1,ig(2))*w2
    fx=fx+gforce(1,ig(3))*w3
    fx=fx+gforce(1,ig(4))*w4
    fx=fx+gforce(1,ig(5))*w5
    fx=fx+gforce(1,ig(6))*w6
    fx=fx+gforce(1,ig(7))*w7
    fx=fx+gforce(1,ig(8))*w8

    fy=fy+gforce(2,ig(1))*w1
    fy=fy+gforce(2,ig(2))*w2
    fy=fy+gforce(2,ig(3))*w3
    fy=fy+gforce(2,ig(4))*w4
    fy=fy+gforce(2,ig(5))*w5
    fy=fy+gforce(2,ig(6))*w6
    fy=fy+gforce(2,ig(7))*w7
    fy=fy+gforce(2,ig(8))*w8
 
    fz=fz+gforce(3,ig(1))*w1
    fz=fz+gforce(3,ig(2))*w2
    fz=fz+gforce(3,ig(3))*w3
    fz=fz+gforce(3,ig(4))*w4
    fz=fz+gforce(3,ig(5))*w5
    fz=fz+gforce(3,ig(6))*w6
    fz=fz+gforce(3,ig(7))*w7
    fz=fz+gforce(3,ig(8))*w8
    !PRINT *,"CHECK FORCES",ipart,time,x,y,z,fx,fy,fz
    !PRINT *,"CHECK WEIGHTS",ipart,time,w1,w2,w3,w4,w5,w6

    vx=vx+fx*t
    vy=vy+fy*t
    vz=vz+fz*t

    part_direct(ipart)%vy=vy
    part_direct(ipart)%vx=vx
    part_direct(ipart)%vz=vz
  enddo
!$OMP ENDDO
#ifdef NOHYDRO
  endif
#endif
  endif
 end subroutine

 subroutine particle_direct_direct()
  integer::ipart,jpart
  real(pre)::xi,yi,zi,mi,si,xj,yj,zj,mj,sj,r,fxi,fyi,fzi,fxj,fyj,fzj
 
  if(npart_direct>0)then
!$OMP DO SCHEDULE(STATIC)
  do ipart=1,npart_direct
   part_direct(ipart)%fx=zero
   part_direct(ipart)%fy=zero
   part_direct(ipart)%fz=zero
  enddo
!$OMP ENDDO
!$OMP BARRIER
!$OMP DO SCHEDULE(STATIC)
  do ipart=1,npart_direct
   if(.not.part_direct(ipart)%active)cycle
   xi=part_direct(ipart)%x
   yi=part_direct(ipart)%y
   zi=part_direct(ipart)%z
   mi=part_direct(ipart)%m
   si=part_direct(ipart)%soft
   do jpart=ipart+1,npart_direct
     if(.not.part_direct(jpart)%active)cycle
     xj=part_direct(jpart)%x
     yj=part_direct(jpart)%y
     zj=part_direct(jpart)%z
     mj=part_direct(jpart)%m
     sj=part_direct(jpart)%soft
     r=sqrt( (xi-xj)**2+(yi-yj)**2+(zi-zj)**2 )
     if(r<sj)then
      fxi=mj/sj**3*( three*r/sj-four)*(xi-xj)
      fyi=mj/sj**3*( three*r/sj-four)*(yi-yj)
      fzi=mj/sj**3*( three*r/sj-four)*(zi-zj)

      if(r<si)then
       fxj=-mi/si**3*( three*r/si-four)*(xi-xj)
       fyj=-mi/si**3*( three*r/si-four)*(yi-yj)
       fzj=-mi/si**3*( three*r/si-four)*(zi-zj)
      else
       fxj=mi/r**3*(xi-xj)
       fyj=mi/r**3*(yi-yj)
       fzj=mi/r**3*(zi-zj)
      endif
     else
      fxi=-mj/r**3*(xi-xj)
      fyi=-mj/r**3*(yi-yj)
      fzi=-mj/r**3*(zi-zj)

      if(r<si)then
       fxj=-mi/si**3*( three*r/si-four)*(xi-xj)
       fyj=-mi/si**3*( three*r/si-four)*(yi-yj)
       fzj=-mi/si**3*( three*r/si-four)*(zi-zj)
      else
       fxj=mi/r**3*(xi-xj)
       fyj=mi/r**3*(yi-yj)
       fzj=mi/r**3*(zi-zj)
      endif
     endif
!$OMP ATOMIC
     part_direct(ipart)%fx=part_direct(ipart)%fx+fxi
!$OMP ATOMIC
     part_direct(ipart)%fy=part_direct(ipart)%fy+fyi
!$OMP ATOMIC
     part_direct(ipart)%fz=part_direct(ipart)%fz+fzi

!$OMP ATOMIC
     part_direct(jpart)%fx=part_direct(jpart)%fx+fxj
!$OMP ATOMIC
     part_direct(jpart)%fy=part_direct(jpart)%fy+fyj
!$OMP ATOMIC
     part_direct(jpart)%fz=part_direct(jpart)%fz+fzj
   enddo
  enddo
!$OMP ENDDO
  endif
 end subroutine

 subroutine particle_direct_interp()
  integer::ipart,jpart
  real(pre)::xi,yi,zi,xj,yj,zj,mi,si,sj,mj,r,fxi,fyi,fzi
 
!$OMP DO SCHEDULE(STATIC)
  do ipart=1,npart
   part(ipart)%fx=zero
   part(ipart)%fy=zero
   part(ipart)%fz=zero
  enddo
!$OMP ENDDO
!$OMP BARRIER
  if(npart_direct>0)then
!$OMP DO SCHEDULE(STATIC)
  do ipart=1,npart_direct
   if(.not.part_direct(ipart)%active)cycle
   xi=part_direct(ipart)%x
   yi=part_direct(ipart)%y
   zi=part_direct(ipart)%z
   mi=part_direct(ipart)%m
   si=part_direct(ipart)%soft
   do jpart=1,npart
     if(.not.part(jpart)%active)cycle
     xj=part(jpart)%x
     yj=part(jpart)%y
     zj=part(jpart)%z
     r=sqrt( (xi-xj)**2+(yi-yj)**2+(zi-zj)**2 )
     if(r<si)then
      fxi=-mi/si**3*( three*r/si-four)*(xi-xj)
      fyi=-mi/si**3*( three*r/si-four)*(yi-yj)
      fzi=-mi/si**3*( three*r/si-four)*(zi-zj)

     else
      fxi=mi/r**3*(xi-xj)
      fyi=mi/r**3*(yi-yj)
      fzi=mi/r**3*(zi-zj)

     endif
!$OMP ATOMIC
     part(jpart)%fx=part(jpart)%fx+fxi
!$OMP ATOMIC
     part(jpart)%fy=part(jpart)%fy+fyi
!$OMP ATOMIC
     part(jpart)%fz=part(jpart)%fz+fzi

   enddo
  enddo
!$OMP ENDDO
  endif
 end subroutine

 subroutine add_direct_togrid()
  integer::jpart,igrid
  real(pre)::xi,yi,zi,xj,yj,zj,sj,mj,r,fxi,fyi,fzi

  print *, "IN add direct to grid"
 
  if(npart_direct>0)then
#ifdef NOHYDRO
 if(.not.npart>0)return
#endif
!$OMP DO SCHEDULE(STATIC)
  do igrid=1,ngrid
   xi=grid(igrid)%x
   yi=grid(igrid)%y
   zi=grid(igrid)%z
   do jpart=1,npart_direct
     if(.not.part_direct(jpart)%active)cycle
     xj=part_direct(jpart)%x
     yj=part_direct(jpart)%y
     zj=part_direct(jpart)%z
     mj=part_direct(jpart)%m
     sj=part_direct(jpart)%soft
     r=sqrt( (xi-xj)**2+(yi-yj)**2+(zi-zj)**2 )
     if(r<sj)then
      fxi=mj/sj**3*( three*r/sj-four)*(xi-xj)
      fyi=mj/sj**3*( three*r/sj-four)*(yi-yj)
      fzi=mj/sj**3*( three*r/sj-four)*(zi-zj)
     else
      fxi=-mj/r**3*(xi-xj)
      fyi=-mj/r**3*(yi-yj)
      fzi=-mj/r**3*(zi-zj)
     endif
!$OMP ATOMIC
     pforce(1,igrid)=pforce(1,igrid)+fxi
!$OMP ATOMIC
     pforce(2,igrid)=pforce(2,igrid)+fyi
!$OMP ATOMIC
     pforce(3,igrid)=pforce(3,igrid)+fzi
   enddo
  enddo
!$OMP ENDDO
  endif
 end subroutine

end module


