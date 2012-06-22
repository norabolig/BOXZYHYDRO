module eos
! Written by A. C. Boley (updated 19 Nov 2010)
! See Pathria for questions.  This uses E=NkT^2 d ln Z /d T to calculate internal energies.
! However, the zero point energies are subtracted out (relevant for orthohydrogen
! and for the vibrational states).  The 16.78 in the metals
! term is from cameron 1968. This takes into account rotational states for 
! hydrogen and 1 vibrational state. zp is the parahydrogen partition function
! and dzpdt and ddzpdtt are its derivatives, *e for equilibrium, and *o for ortho.
 use parameters
 use derived_types
 use grid_commons
 implicit none

 type(units)::scl
 real(pre)::eul,log_rho_eos_low
 real(pre),parameter::kb=1.38065d-16
 real(pre),parameter::hplanck=6.626d-27
 integer::NEOS

 contains

   subroutine initialize_eos()

     call get_units(scl)
     NEOS_T=int((tk_eos_cutoff-tk_bgrnd)/dTk_eos)+1
     drho_eos=(log10(rho_eos_high)-log10(rho_eos_low))/(NEOS_RHO-1)
     print *,"NEOS_T for large table  = ",NEOS_T
     print *,"drho_eos for large table = ",drho_eos
     NEOS=NEOS_T

     allocate(gamma_table(NEOS_T,NEOS_RHO))
     allocate(gamma_table2(NEOS_T,NEOS_RHO))
     allocate(eng_table(NEOS_T,NEOS_RHO))
     allocate(eng_table2(NEOS_T,NEOS_RHO))
     allocate(tk_table(NEOS_T))
     allocate(p_table(NEOS_T,NEOS_RHO))
     allocate(tk_table2(NEOS_T,NEOS_RHO))
     allocate(rho_table(NEOS_RHO))
     allocate(muc_table(NEOS_T,NEOS_RHO))
     allocate(muc_table2(NEOS_T,NEOS_RHO))
     allocate(deng_eos_array(NEOS_RHO))

   end subroutine

   subroutine clean_eos()
     deallocate(gamma_table,tk_table,eng_table,gamma_table2,eng_table2,tk_table2,muc_table,muc_table2,rho_table)
   end subroutine

   subroutine h2_partition(t,zp,zo,ze,zoprime,dzpdt,dzodt,dzedt,dzoprimedt)
     real(pre),intent(in)::t
     real(pre),intent(out)::zp,zo,ze,zoprime,dzpdt,dzodt,dzedt,dzoprimedt
     real(pre)::f1,f2,f3,f4
     integer::j
     zp=zero
     zo=zero
     ze=zero
     zoprime=zero 
     dzpdt=zero
     dzodt=zero
     dzedt=zero
     dzoprimedt=zero
     do j=0,24,2
       f1=dble(2*j+1)
       f2=dble(j*(j+1))
       f3=dble(2*(j+1)+1)
       f4=dble((j+1)*(j+2))
       zp=zp+f1*exp(-f2*brot/t)
       zo=zo+three*f3*exp(-f4*brot/t)
       zoprime=zoprime+three*f3*exp(-(f4-two)*brot/t)
       dzpdt=dzpdt+f1*f2*brot*exp(-f2*brot/t)/t**2
       dzodt=dzodt+three*f3*f4*brot*exp(-f4*brot/t)/t**2
       dzoprimedt=dzoprimedt+three*f3*(f4-two)*brot*exp(-(f4-two)*brot/t)/t**2
     enddo
     ze=zp+zo
     dzedt=dzpdt+dzodt
     return
   end subroutine

   real(pre) function translate(m,t,n)
     real(pre)::m,t,n
     translate=eul*(sqrt(two*pi*m*mp*kB*t/hplanck**2))**3/n
     return
   end function

   real(pre) function debroglie(m,t)
     real(pre)::m,t
     debroglie=hplanck/sqrt(two*pi*m*mp*kB*t)
     return
   end function

   subroutine calc_eos_table()
     real(pre)::t0,tm,tp,denm,denp,den0,logdenp,logdenm,logden0
     real(pre)::logt0,logtm,gamu=zero,gaml=zero,xdiss,apot
     real(pre)::zp,zo,ze,dzpdt,dzodt,dzedt,zoprime,dzoprimedt
     real(pre)::nh,nh2,nhe,nz,trans_h,trans_he,trans_h2,trans_z,zh2_int,rhs
     real(pre)::den,gam,s0,s1,s2,s3,tk,log_eng_low

     real(pre),dimension(:),allocatable::xdiss_a
     real(pre),dimension(:,:),allocatable::sent

     integer i,irho,itk,irhop,irhom

     call get_units(scl)
     eul=exp(one)
   
     muc = one/(xabun*half + yabun*quarter + zabun/mu_z)
     print *, "muc is:", muc

     allocate(sent(NEOS_T,NEOS_RHO))
     allocate(xdiss_a(NEOS_T))

     do i = 1,NEOS_T
       tk_table(i)=tk_bgrnd + dble(i-1)*dTk_eos
     enddo
     log_rho_eos_low=log10(rho_eos_low)
     do i = 1,NEOS_RHO
       rho_table(i)=ten**(log_rho_eos_low+dble(i-1)*drho_eos)
     enddo
     !log_rho_eos_low=log_rho_eos_low-log10(scl%density)
     !drho_eos=drho_eos-log10(scl%density)

! we've done the easy part.  Now we need to make the large table. Same idea, but must be done for each rho and for each

     do irho=1,NEOS_RHO
       den=rho_table(irho)
       nh=xabun*den/(mp)
       nh2=nh*half
       nhe=yabun*den/(four*mp)
       nz=zabun*den/(mu_z*mp)
    
       select case(H2STAT)
       case(0)
         do itk=1,NEOS_T
           tk=tk_table(itk) 
           trans_he=one
           trans_z=one
      
           trans_h=translate(one,tk,nh)
           trans_h2=translate(two,tk,nh2)
           if(yabun>zero)trans_he=translate(four,tk,nhe)
           if(zabun>zero)trans_z=translate(mu_z,tk,nz)
            
           call h2_partition(tk,zp,zo,ze,zoprime,dzpdt,dzodt,dzedt,dzoprimedt)

           zh2_int=zp**(ac/(ac+bc))*zoprime**(bc/(ac+bc))/(one-exp(-vib/tk))
           rhs=(debroglie(two,tk)/debroglie(one,tk)**2)**3/(two*nh*zh2_int)*exp(-diss/tk)
           xdiss_a(itk)=rhs*half+one-half*sqrt((rhs+two)**2-four)
           xdiss=xdiss_a(itk)
           
           eng_table(itk,irho)=rgasCGS*tk*(xabun* ( (one-xdiss)*(1.5d0+half*diss/tk ) &
               + half*xdiss*(1.5d0+tk*(ac/(ac+bc)*dzpdt/zp+bc/(ac+bc)*(dzoprimedt/zoprime) &
               + vib/tk**2*exp(-vib/tk)/(one-exp(-vib/tk)) ) ) ) &
               + yabun*0.25d0*1.5d0 + zabun*1.5d0/mu_z)
           apot=-rgasCGS*tk*( xabun * ( (one-xdiss)*(log(trans_h)-half*diss/tk) + &
              half*xdiss*( log(trans_h2) + ac/(ac+bc)*log(zp) + bc/(ac+bc)*log(zoprime) &
              - log(one-exp(-vib/tk))) ) + yabun*0.25d0*log(trans_he) + zabun*log(trans_z)/mu_z)
           sent(itk,irho)=(eng_table(itk,irho)-apot)/tk
         enddo

       case(1)

         do itk=1,NEOS_T
           tk=tk_table(itk) 
           trans_he=one
           trans_z=one
      
           trans_h=translate(one,tk,nh)
           trans_h2=translate(two,tk,nh2)
           if(yabun>zero)trans_he=translate(four,tk,nhe)
           if(zabun>zero)trans_z=translate(mu_z,tk,nz)
            
           call h2_partition(tk,zp,zo,ze,zoprime,dzpdt,dzodt,dzedt,dzoprimedt)

           zh2_int=ze/(one-exp(-vib/tk))
           rhs=(debroglie(two,tk)/debroglie(one,tk)**2)**3/(two*nh*zh2_int)*exp(-diss/tk)
           xdiss_a(itk)=rhs*half+one-half*sqrt((rhs+two)**2-four)
           xdiss=xdiss_a(itk)
           
           eng_table(itk,irho)=rgasCGS*tk*(xabun* ( (one-xdiss)*(1.5d0+half*diss/tk ) &
               + half*xdiss*(1.5d0+tk*(dzedt/ze &
               + vib/tk**2*exp(-vib/tk)/(one-exp(-vib/tk)) ) ) ) &
               + yabun*0.25d0*1.5d0 + zabun*1.5d0/mu_z)
           apot=-rgasCGS*tk*( xabun * ( (one-xdiss)*(log(trans_h)-half*diss/tk) + &
              half*xdiss*( log(trans_h2) + log(ze)     &
              - log(one-exp(-vib/tk))) ) + yabun*0.25d0*log(trans_he) + zabun*log(trans_z)/mu_z)
           sent(itk,irho)=(eng_table(itk,irho)-apot)/tk
         enddo
       case(-1)
         !print *, "Single gamma.  Generating table for completeness, but it is not used."
         do itk =1,NEOS_T
           tk=tk_table(itk) 
           eng_table(itk,irho)=rgasCGS*tk/(gammafix-one)/muc
           apot=-rgasCGS*tk/muc
           xdiss_a(itk)=one
           sent(itk,irho)=(eng_table(itk,irho)-apot)/tk
           gamma_table(itk,irho)=gammafix
         enddo
       end select
       do itk=1,NEOS_T
         muc_table(itk,irho)=(one/(xabun*(one-half*xdiss_a(itk)) + yabun*0.25d0 + zabun/mu_z))
       enddo
     enddo ! loop over density 

     ! This part is a pain, but it is a straight-forward way to derive the adiabatic index
     if(H2STAT>=0)then
     do irho=1,NEOS_RHO
       irhop=irho+1
       irhom=irho-1
       den0=rho_table(irho)
       logden0=log(den0)
       denm=rho_table(irhom)
       logdenm=log(denm)
       denp=rho_table(irhop)
       logdenp=log(denp)
       
       if(irho==1)then
         do itk=1,NEOS_T
           if(itk<NEOS_T)then
             s0=sent(itk,irho)
             s1=sent(itk,irhop)
             s2=sent(itk+1,irho)
             s3=sent(itk+1,irhop)
             t0=tk_table(itk)
             logt0=log(t0)
             tp=(s0-s1)*( log(tk_table(itk+1))-logt0)/(s3-s1)+logt0
             if (s1 < s0 .and. s0 < s3)then
              gam=one+(tp-logt0)/(logdenp-logden0)
             else
              tp=tk_table(itk+1)
              den=(s0-s2)*(logdenp-logden0)/(s3-s2)+logden0
              gam=one+(log(tp)-logt0)/(den-logden0)
             endif
           else
             s0=sent(itk-1,irho)
             s1=sent(itk-1,irhop)
             s2=sent(itk,irho)
             s3=sent(itk,irhop)
             t0=tk_table(itk-1)
             logt0=log(t0)
             tp=(s0-s1)*( log(tk_table(itk))-logt0)/(s3-s1)+logt0
             if (s1 < s0 .and. s0 < s3)then
              gam=one+(tp-logt0)/(logdenp-logden0)
             else
              tp=tk_table(itk)
              den=(s0-s2)*(logdenp-logden0)/(s3-s2)+logden0
              gam=one+(log(tp)-logt0)/(den-logden0)
             endif
           endif
           gamma_table(itk,irho)=gam
         enddo     
         cycle ! cycle rho
       endif
       if(irho==NEOS_RHO)then
         do itk=1,NEOS_T
          if(itk<NEOS_T)then
            s0=sent(itk,irhom)
            s1=sent(itk,irho)
            s2=sent(itk+1,irhom)
            s3=sent(itk+1,irho)
            t0=tk_table(itk+1)
            tm=tk_table(itk)
            logt0=log(t0)
            logtm=log(tm)
            tp=logtm-(s1-s0)*(logt0-logtm)/(s2-s0)
            !if(s0<s3.and.s3<s2)then
            if(s1<s0.and.s0<s3)then
              gam=one+(tp-logtm)/(logden0-logdenm)
            else
              den=logdenm-(s3-s1)*(logdenm-logden0)/(s0-s1)
              gam=one+(logt0-logtm)/(den-logdenm)
            endif
          else
            s3=sent(itk,irho)
            s1=sent(itk-1,irho)
            s0=sent(itk-1,irhom)
            s2=sent(itk,irhom)
            t0=tk_table(itk)
            tm=tk_table(itk-1)
            logt0=log(t0)
            logtm=log(tm)
            tp=logtm-(s3-s2)*(logtm-logt0)/(s0-s2)
            if(s0<s3.and.s3<s2)then
              gam=one+(tp-logtm)/(logden0-logdenm)
            else
              den=logdenm-(s3-s1)*(logdenm-logden0)/(s0-s1)
              gam=one+(logt0-logtm)/(den-logdenm)
            endif
          endif
          gamma_table(itk,irho)=gam
         enddo
         cycle ! cycle rho
       endif
       do itk=1,NEOS_T
          if(itk<NEOS_T)then
            s0=sent(itk,irho)
            s1=sent(itk,irhop)
            s2=sent(itk+1,irho)
            s3=sent(itk+1,irhop)
            t0=tk_table(itk)
            logt0=log(t0)
            tp=(s0-s1)*(log(tk_table(itk+1))-logt0)/(s3-s1)+logt0
            if(s1<s0.and.s0<s3)then
              gamu=one+(tp-logt0)/(logdenp-logden0) 
            else
              tp=tk_table(itk+1)
              den=(s0-s2)*(logdenp-logden0)/(s3-s2)+logden0
              gamu=one+(log(tp)-logt0)/(den-logden0)
            endif
          endif
          if(itk>1)then
            s3=sent(itk,irho)
            s1=sent(itk-1,irho)
            s0=sent(itk-1,irhom)
            s2=sent(itk,irhom)
            t0=tk_table(itk)
            tm=tk_table(itk-1)
            logt0=log(t0)
            logtm=log(tm)
            tp=logtm-(s3-s2)*(logtm-logt0)/(s0-s2)
            if(s1<s0.and.s0<s3)then
              gaml=one+(tp-logtm)/(logden0-logdenm) 
            else
              den=logdenm-(s3-s1)*(logdenm-logden0)/(s0-s1)
              gaml=one+(logt0-logtm)/(den-logdenm)
            endif
          endif
          if(itk<NEOS_T.and.itk>1)then 
            gam=half*(gamu+gaml)
          elseif(itk<NEOS_T)then
            gam=gamu
          else
            gam=gaml
          endif
          gamma_table(itk,irho)=gam 
       enddo
     enddo
     endif

#ifdef VERBOSE
  print *, "#EOS table initialized.  XABUN, YABUN, ZABUN, MU_Z, MUC ",xabun,yabun,zabun,mu_z,muc
#endif

     gamma_table(:,1)=gamma_table(:,2)
     gamma_table(:,NEOS_RHO)=gamma_table(:,NEOS_RHO-1)
     do irho=1,NEOS_RHO
      deng_eos=(log10(eng_table(NEOS_T,irho))-log10(eng_table(1,irho)))/dble(NEOS_T-1)
      deng_eos_array(irho)=deng_eos
      log_eng_low=log10(eng_table(1,irho))
      do itk=1,NEOS_T
        eng_table2(itk,irho)=ten**(log_eng_low+dble(itk-1)*deng_eos)
        call get_gamma_norho(eng_table2(itk,irho),tk_table2(itk,irho), &
                 muc_table2(itk,irho),gamma_table2(itk,irho),irho)
      enddo
     enddo

     open(unit=102,file="eostable.dat")
     do irho=1,NEOS_RHO
      write(102,"(A,1X,1pe15.8)")"#den: ",rho_table(irho)
      do itk=1,NEOS_T
        p_table(itk,irho)=rho_table(irho)/scl%density*scl%rgas*tk_table(itk)/muc_table(itk,irho)
        write(102,"(10(1X,1pe15.8))")rho_table(irho),tk_table(itk), &
              eng_table(itk,irho),muc_table(itk,irho),gamma_table(itk,irho), &
              tk_table2(itk,irho),eng_table2(itk,irho),muc_table2(itk,irho),gamma_table2(itk,irho),p_table(itk,irho)
        eng_table(itk,irho)=eng_table(itk,irho)*(scl%time/scl%length)**2
        eng_table2(itk,irho)=eng_table2(itk,irho)*(scl%time/scl%length)**2
      enddo
      rho_table(irho)=rho_table(irho)/scl%density
     enddo
     close(102)

     deallocate(sent)

   end subroutine calc_eos_table

   subroutine get_gamma_norho(eng,tk,m,gam,irho)
    real(pre)::eng,tk,gam,m
    integer::ientry,jump,flag,inext,irho
    
    ientry=1
    jump=NEOS/4 ! we are usually at low T, so take a small jump.
    flag=0
    do
      inext=min(ientry+1,NEOS)
      if(eng_table(ientry,irho)<=eng.and.eng<eng_table(inext,irho))exit
      if(eng_table(ientry,irho)>eng)then
        ientry=ientry-jump
        jump=max(int(jump*.75),1)
        if(ientry<1)then
          ientry=1
          if(eng<=eng_table(ientry,irho))then
            flag=1
            exit
          endif
        endif
      else
        ientry=ientry+jump
        jump=max(int(jump*.75),1)
        if(ientry>NEOS-1)then
         ientry=NEOS-1
         if(eng>=eng_table(ientry,irho))then
            flag=2
            exit
         endif
        endif
      endif
    enddo
   if(flag>1)then
     eng=eng_table(NEOS,irho)
   elseif(flag>0)then
     eng=eng_table(1,irho)
   endif
   tk=tk_table(ientry)+(tk_table(ientry+1)-tk_table(ientry)) &
     /(eng_table(ientry+1,irho)-eng_table(ientry,irho))*(eng-eng_table(ientry,irho))
   gam=gamma_table(ientry,irho)+(gamma_table(ientry+1,irho)-gamma_table(ientry,irho)) &
     /(eng_table(ientry+1,irho)-eng_table(ientry,irho))*(eng-eng_table(ientry,irho))
   m=muc_table(ientry,irho)+(muc_table(ientry+1,irho)-muc_table(ientry,irho)) &
     /(eng_table(ientry+1,irho)-eng_table(ientry,irho))*(eng-eng_table(ientry,irho))


   end subroutine

   subroutine get_gamma(eps,rho,tk,m,gam)
    real(pre)::eng,tk,gam,eps,rho,m
    integer::ientry,jump,flag,inext,irho
    
    irho=min(int( (log10(rho*scl%density)-log_rho_eos_low)/drho_eos) + 1 ,NEOS_RHO)
    if(irho<1)irho=1
 
    eng=eps/rho
    ientry=1
    jump=NEOS/4 ! we are usually at low T, so take a small jump.
    flag=0
    do
      inext=min(ientry+1,NEOS)
      if(eng_table(ientry,irho)<=eng.and.eng<eng_table(inext,irho))exit
      if(eng_table(ientry,irho)>eng)then
        ientry=ientry-jump
        jump=max(int(jump*.75),1)
        if(ientry<1)then
          ientry=1
          if(eng<=eng_table(ientry,irho))then
            flag=1
            exit
          endif
        endif
      else
        ientry=ientry+jump
        jump=max(int(jump*.75),1)
        if(ientry>NEOS-1)then
         ientry=NEOS-1
         if(eng>=eng_table(ientry,irho))then
            flag=2
            exit
         endif
        endif
      endif
    enddo
   if(flag>1)then
     eng=eng_table(NEOS,irho)
     eps=eng*rho ! note that this does change the energy
   elseif(flag>0)then
     eng=eng_table(1,irho)
     eps=eng*rho ! note that this does change the energy
   endif
   tk=tk_table(ientry)+(tk_table(ientry+1)-tk_table(ientry)) &
     /(eng_table(ientry+1,irho)-eng_table(ientry,irho))*(eng-eng_table(ientry,irho))
   gam=gamma_table(ientry,irho)+(gamma_table(ientry+1,irho)-gamma_table(ientry,irho)) &
     /(eng_table(ientry+1,irho)-eng_table(ientry,irho))*(eng-eng_table(ientry,irho))
   m=muc_table(ientry,irho)+(muc_table(ientry+1,irho)-muc_table(ientry,irho)) &
     /(eng_table(ientry+1,irho)-eng_table(ientry,irho))*(eng-eng_table(ientry,irho))


   end subroutine

   subroutine get_gamma2(eps,rho,tk,m,gam)
    real(pre)::eng,tk,gam,eps,rho,m
    integer::ientry,irho
 
    irho=min(int( (log10(rho*scl%density)-log_rho_eos_low)/drho_eos) + 1 ,NEOS_RHO)
    if(irho<1)irho=1

    eng=eps/rho
    !print *, rho,eng,eps,irho
    ientry=int( (log10(eng)-log10(eng_table2(1,irho)))/deng_eos_array(irho)) + 1
    if (ientry>NEOS-1)then
         ientry=NEOS-1
         eng=eng_table2(NEOS,irho)
         eps=eng*rho
    endif
    if(ientry<1.or.eng<eng_table2(1,irho))then
      !print *, rho,eng,eps,irho
        ientry=1
        eng=eng_table2(1,irho)
        eps=eng*rho
    endif
   
    gam=gamma_table2(ientry,irho)+(gamma_table2(ientry+1,irho)-gamma_table2(ientry,irho))&
            /(eng_table2(ientry+1,irho)-eng_table2(ientry,irho))*(eng-eng_table2(ientry,irho))
    m=muc_table2(ientry,irho)+(muc_table2(ientry+1,irho)-muc_table2(ientry,irho))&
            /(eng_table2(ientry+1,irho)-eng_table2(ientry,irho))*(eng-eng_table2(ientry,irho))
    tk=tk_table2(ientry,irho)+(tk_table2(ientry+1,irho)-tk_table2(ientry,irho))&
            /(eng_table2(ientry+1,irho)-eng_table2(ientry,irho))*(eng-eng_table2(ientry,irho))

   end subroutine

   subroutine get_gamma_from_tk(eps,rho,tk,m,gam)
    real(pre)::tk,gam,eps,rho,m
    integer::ientry,irho
  
    irho=min(int( (log10(rho*scl%density)-log_rho_eos_low)/drho_eos) + 1 ,NEOS_RHO)
    if(irho<1)irho=1
   
    ientry=int( (tk-tk_bgrnd)/dTk_eos)+1
    if (ientry>NEOS-1)ientry=NEOS-1
    if(ientry<1)ientry=1

    eps=rho*(eng_table(ientry,irho)+(eng_table(ientry+1,irho)-eng_table(ientry,irho))/dTK_eos*(tk-tk_table(ientry)))
    gam=(gamma_table(ientry,irho)+(gamma_table(ientry+1,irho)-gamma_table(ientry,irho))/dTK_eos*(tk-tk_table(ientry)))
    m=(muc_table(ientry,irho)+(muc_table(ientry+1,irho)-muc_table(ientry,irho))/dTK_eos*(tk-tk_table(ientry)))

   end subroutine

   subroutine get_gamma_from_p(eps,rho,p_loc,m,gam)
    real(pre)::p_loc,gam,eps,rho,m
    integer::ientry,jump,flag,inext,irho,i,irhop,irho0
    real(pre)::eps0,eps1,m0,m1,gam0,gam1
    
    irho0=min(int( (log10(rho*scl%density)-log_rho_eos_low)/drho_eos) + 1 ,NEOS_RHO)
    if(irho0<1)irho0=1
    irhop=min(irho0+1,NEOS_RHO)
    m0=zero
    gam0=zero
    eps0=zero
    m1=zero
    gam1=zero
    eps1=zero
 
    do i=1,2
    irho=irho0
    if(i==2)irho=irhop
    ientry=1
    jump=NEOS/4 ! we are usually at low T, so take a small jump.
    flag=0
    do
      inext=min(ientry+1,NEOS)
      if(p_table(ientry,irho)<=p_loc.and.p_loc<p_table(inext,irho))exit
      if(p_table(ientry,irho)>p_loc)then
        ientry=ientry-jump
        jump=max(int(jump*.75),1)
        if(ientry<1)then
          ientry=1
          if(p_loc<=p_table(ientry,irho))then
            flag=1
            exit
          endif
        endif
      else
        ientry=ientry+jump
        jump=max(int(jump*.75),1)
        if(ientry>NEOS-1)then
         ientry=NEOS-1
         if(p_loc>=p_table(ientry,irho))then
            flag=2
            exit
         endif
        endif
      endif
    enddo
   if(flag>1)then
     !print *, p_loc,p_table(NEOS,irho),rho*scl%density,irho,log_rho_eos_low,(log10(rho*scl%density)-log_rho_eos_low)/drho_eos
     p_loc=p_table(NEOS,irho)
   elseif(flag>0)then
     p_loc=p_table(1,irho)
   endif
   select case(i)
   case(1)
   eps0=eng_table(ientry,irho)+(eng_table(ientry+1,irho)-eng_table(ientry,irho)) &
     /(p_table(ientry+1,irho)-p_table(ientry,irho))*(p_loc-p_table(ientry,irho))
   gam0=gamma_table(ientry,irho)+(gamma_table(ientry+1,irho)-gamma_table(ientry,irho)) &
     /(p_table(ientry+1,irho)-p_table(ientry,irho))*(p_loc-p_table(ientry,irho))
   m0=muc_table(ientry,irho)+(muc_table(ientry+1,irho)-muc_table(ientry,irho)) &
     /(p_table(ientry+1,irho)-p_table(ientry,irho))*(p_loc-p_table(ientry,irho))

   case(2)
   eps1=eng_table(ientry,irho)+(eng_table(ientry+1,irho)-eng_table(ientry,irho)) &
     /(p_table(ientry+1,irho)-p_table(ientry,irho))*(p_loc-p_table(ientry,irho))
   gam1=gamma_table(ientry,irho)+(gamma_table(ientry+1,irho)-gamma_table(ientry,irho)) &
     /(p_table(ientry+1,irho)-p_table(ientry,irho))*(p_loc-p_table(ientry,irho))
   m1=muc_table(ientry,irho)+(muc_table(ientry+1,irho)-muc_table(ientry,irho)) &
     /(p_table(ientry+1,irho)-p_table(ientry,irho))*(p_loc-p_table(ientry,irho))
   end select
   enddo

   if (irhop>irho0)then
   eps=(eps0+(eps1-eps0)*(rho-rho_table(irho0))/(rho_table(irhop)-rho_table(irho0)))*rho
   gam=gam0+(gam1-gam0)*(rho-rho_table(irho0))/(rho_table(irhop)-rho_table(irho0))
   m=m0+(m1-m0)*(rho-rho_table(irho0))/(rho_table(irhop)-rho_table(irho0))
   else
     eps=eps0*rho
     m=m0
     gam=gam0
   endif 

   end subroutine



end module
