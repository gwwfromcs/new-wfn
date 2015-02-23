!*******************************************
!*** This program using (v_in+dv, rho) -> (v_out+dv) to fixed the
!*** Hamiltonian: vion (so V_hxc(rho)+vion=v_out+dv)
!***              dvp: (so v_in+dv & rho, Ef0 satisfy the solution eq.)
!***  Then, it resolve rho', so that v_in'=v_out'
!***  Finally, the v_in' is used for small k 
!***  and kerker mixing is used for large k
!****************************************************************
!*** The final equation for psi=dsqrt(rho)
!*** -0.5 \nabla^2 psi +( rho^(2/3)+ v_in)*psi + dvp = Ef*psi
!***  and v_in= vion + v_Hxc(rho)
!****************************************************************

      subroutine Thomas_Fermi(vin,vout,vion,Ef0,max_iter, &
                       ffts,crys,gs, den,den_core,vnl_ave)
  !
  !     2001 by L.W.Wang from PEtot and adapted into Paratec by Dave Raczkowski
  !
!*****************************************************************
!*** in return, v_in will be replaced by new v_in
!*** rho, v_out, dv are not changed
!*****************************************************************
  use all_to_all_module
  include 'use.h'  
  IMPLICIT NONE
  include 'interface.h' 
  include 'all_to_all.h'  
  include 'flibcalls.ph'  
  !
  !     INPUT:
  !     -----
  ! 
  type(crystal), intent(in) :: crys            ! for vcell, nspin, ztot
  type(parallel_gspace), intent(in) ::  gs  
  type(fft_struc), intent(in) :: ffts 
  real(dp), intent(in) :: den(2*gs%length,crys%nspin),den_core(2*gs%length)
  real(dp) Ef0      ! fermi level
  real(dp), intent(in) :: vnl_ave(2*gs%length)
  integer, intent(in) :: max_iter
  !
  !     INPUT/OUTPUT:
  !     -------------
  !
  real(dp), intent(inout) :: &
       vin(2*gs%length),&    ! screening pot. which the eigvectors were calc.
       vout(2*gs%length),&   ! screening pot. which the eigvectors were calc.  
       vion(2*gs%length)
  !
  !
  !     --------------- local variables -------------------------
  !
  real(dp) work1_c(max(gs%length*2,2*gs%r_size) ),&
            work2_c(max(gs%length*2,2*gs%r_size))
  real(dp) den0(max(gs%length*2,gs%r_size)),&
           ugr(max(gs%length*2,gs%r_size)),&
            ughr(max(gs%length*2,2*gs%r_size)),&
           pgr(max(gs%length*2,gs%r_size)),&
           pgr_old(max(gs%length*2,gs%r_size)),&
          dvp(gs%r_size),&
          vion_r(gs%r_size)

  integer i,nr,nint
  real(dp) cos_th,sin_th,etot1,theta0,y,beta ,e_upg,&
            ecut0,theta ,x ,fa, UxcCA,s1,s2,s3,alpha,err,etot,rr0,&
            rr1,Ek,E_thth,predE,Ef,s,uxc2,rr2,Etot_old
  logical ireset

!********* vion, dvp are derived from input v_in,v_out,rho

  nr=gs%fftsize(1)*gs%fftsize(2)*gs%fftsize(3)
  alpha=1d0
  pgr_old = 0.0
  ireset = .false.
  Etot_old =1d8

  call mdaxpy(2*gs%length,done,vion(1),1,vin(1),1)  
  call mdaxpy(2*gs%length,done,vion(1),1,vout(1),1)  

! need to remove overall volume scaling in charge density
  call mdcopy(2*gs%length,den(1,1),1,den0(1),1)
  if (crys%nspin > 1)   call mdaxpy(2*gs%length,done,den(1,2),1,den0(1),1)

  call mdscal(2*gs%length,done/crys%vcell,den0(1),1)

  call fourier_transform(-1, ffts, gs, den0(1), work1_c(1), 1)
! den0=den
  call mdcopy(gs%r_size,work1_c(1),2,den0(1),1)

!********************************************
!*** make rho(i) = workr(i)**2 and workr(i) is in the Ecut2
!********************************************
  s1=dzero
  work1_c=dzero
  do i=1,gs%r_size
    s1=s1+den0(i)
    work1_c(2*i-1)=dsqrt(abs(den0(i)))
  enddo

  call all_sum_all(s1) 

  call fourier_transform(1, ffts, gs, work2_c(1), work1_c(1), 1)
  call fourier_transform(-1, ffts, gs, work2_c(1), work1_c(1), 1)

  s2=dzero

  do i=1,gs%r_size
    den0(i)=work1_c(2*i-1)**2
    s2=s2+den0(i)
  enddo

  call all_sum_all(s2) 

  s=s1/s2

  call mdscal(gs%r_size,s,den0(1),1)

!*********************************************
!   PREPARE VION
!*********************************************
  work1_c=dzero
  call mdcopy(gs%r_size,den0(1),1,work1_c(1),2)

  call fourier_transform(1, ffts, gs, work2_c(1), work1_c(1), 1)
! pgr=fft of charge density
  call mdcopy(2*gs%length,work2_c(1),1,pgr(1),1)

! put coulomb interaction in pgr
  do i=1,gs%length
    if (gs%ekin(i) .gt. dzero) then 
      pgr(2*i)=4*pgr(2*i)*2*pi/gs%ekin(i)   
      pgr(2*i-1)=4*pgr(2*i-1)*2*pi/gs%ekin(i)  
    else           
      pgr(2*i)=dzero
      pgr(2*i-1)=dzero
    end if
  end do

! work1_c=pgr
  call mdcopy(2*gs%length,pgr(1),1,work2_c(1),1)
  call fourier_transform(-1, ffts, gs, work2_c(1), work1_c(1), 1)
! pgr=real( work1_c )
  call mdcopy(gs%r_size,work1_c(1),2,pgr(1),1)  ! pgr has coulomb in rsp now

! inverse fft vout and core charge density
  call fourier_transform(-1, ffts, gs, den_core(1), work2_c(1), 1)
  call mdscal(2*gs%r_size,dzero, work2_c(1),1)
  call fourier_transform(-1, ffts, gs, vout(1), work1_c(1), 1)

  call fourier_transform(-1, ffts, gs, vnl_ave(1), ughr(1), 1)

  do i=1,gs%r_size
    pgr(i)=pgr(i)+2*UxcCA(den0(i)+work2_c(2*i-1),uxc2)  ! add in Uxc
! add in local Vps and new vout. seems screwy
    vion_r(i)=(work1_c(2*i-1) + ughr(2*i-1))-pgr(i) 
  enddo

!*********************************************
! PREPARE dvp 
!*********************************************

! inverse fft vin
  call fourier_transform(-1, ffts, gs, vin(1), work2_c(1), 1)
 
! pgr now equals exchange? + vin +vloc
  fa=(3*pi**2)**(2.d0/3.d0)
  do i=1,gs%r_size
    pgr(i)=fa*den0(i)**(2.d0/3.d0)+work2_c(2*i-1) + ughr(2*i-1) 
  enddo

! PUT (KINETIC ENERGY)*(charge density) INTO ugr
  work1_c=dzero
  do i=1,gs%r_size
    work1_c(2*i-1)=dsqrt(den0(i))
  enddo

  call fourier_transform(1, ffts, gs, ugr(1), work1_c(1), 1)

  do i=1,gs%length
    ugr(i*2)=alpha*gs%ekin(i)*ugr(i*2)
    ugr(i*2-1)=alpha*gs%ekin(i)*ugr(i*2-1)
   enddo

  call fourier_transform(-1, ffts, gs, ugr(1), work1_c(1), 1)

! get difference (Hdft - Htf) * psi^1/2
  do i=1,gs%r_size
    dvp(i)=(Ef0-pgr(i))*dsqrt(den0(i)) -work1_c(2*i-1)   
    dvp(i)=dvp(i)/dsqrt(done*crys%ztot)         
  enddo

!ccccccccccccccccccccccccccccccccccccccccccccccc
!ccc finished obtaining vion(r) and dvp(r)
!ccc Now, the Hamiltonian has be determined, so just need
!ccc to solve the Hamiltonian selfconsistently
!ccccccccccccccccccccccccccccccccccccccccccccccc
  s1=dzero
  s2=dzero
  s3=dzero

  do i=1,gs%r_size
    s1=s1+dvp(i)**2
    s2=s2+den0(i)
!    s3=s3+dabs(work1_c(2*i-1))  !vionT_n
   enddo

  call all_sum_all(s1) 
  call all_sum_all(s2) 
!  call all_sum_all(s3) 

  s1=dsqrt(s1/s2)
!      s3=s3/nr

  write(9,*) "average dvp (Hartree)", s1
!  write(6,*) "average v_in+dv", s3 
  
  
!******************************************************
!*****, Now, the potential vion(i) and dvp(i) has been calculated
!******************************************************


  work1_c=dzero
  do i=1,gs%r_size
    work1_c(2*i-1)=dsqrt(den0(i)/crys%ztot)
  enddo
  call fourier_transform(1, ffts, gs, ugr(1), work1_c(1), 1)

! ugr is averaged CD^1/2. It is like our wavefunction for TF theory

  do 1000 nint=1,max_iter

    if(nint.eq.1) rr0=1.D+40

! NORMALIZE ugr
    s=mddot(2*gs%length, ugr(1), 1, ugr(1), 1)    
    call all_sum_all(s)
    s=s*crys%vcell  
    s=dsqrt(done/s)
    call mdscal(2*gs%length,s,ugr(1),1)  !ugr=s*ugr

! CALCULATE H*u=ughr and Etot
    call Hpsi3(ugr,ughr,vion_r,dvp,Etot,alpha,work1_c,&
                 ffts,crys,gs,den_core) 

! pgr=ughr-H*CD
    call mdcopy(2*gs%length,ughr(1),1,pgr(1),1)
!Ef=(pgr.ugr)  it is eigenvalue
    Ef=mddot(2*gs%length, pgr(1), 1, ugr(1), 1)
    call all_sum_all(Ef)
    Ef=Ef*crys%vcell  

!pgr=pgr-Ef*ugr   H*CD - (CD'*H*CD)CD    THIS IS THE GRADIENT
    call mdaxpy(2*gs%length,-Ef,ugr(1),1,pgr(1),1)
  
!err=(pgr.pgr)  NORM OF GRADIENT
    err=mddot(2*gs%length, pgr(1), 1,pgr(1), 1)
    call all_sum_all(err)
    err=dsqrt(err*crys%vcell) 

    if((nint.eq.1.or.nint.eq.max_iter) ) then
      write(9,*) "Ef,err",Ef,err,Etot,nint,dcos(theta)
    endif

     if (err .lt. 1d-12) go to 2000

!   Ek=0.5d0 
     Ek=done   ! now in Rydbergs was 0.5d0 in hartree
    rr1=dzero
    rr2=dzero

    if (Etot .gt. Etot_old+1d-6 .and. .not. ireset) then
!     NORMALIZE AGAIN
!s=(pgr.pgr)  s=(Z.Z)
   call mdcopy(2*gs%length,pgr_old(1),1,pgr(1),1)
       s=mddot(2*gs%length, pgr(1), 1,pgr(1), 1)
       call all_sum_all(s)
       s=1.d0/dsqrt(s*crys%vcell) 
!pgr=s*pgr  Z=(Z.Z)Z
       call mdscal(2*gs%length,s,pgr(1),1)

       call mdaxpy(2*gs%length,-sin_th,pgr(1),1,ugr(1),1)
       call mdscal(2*gs%length,1/cos_th,ugr(1),1)
       write(9,*) 'reset'
      write(9,*) "Ef,err",Ef,err,Etot,nint,dcos(theta)
       ireset = .true.
       go to 1000
!       go to 2000
    end if 
    Etot_old=Etot

    if (nint .gt. 1 .and. .not. ireset) then

! calculate preconditioned u*u
    do i=1,gs%length
      x=gs%ekin(i)/Ek
      y=1/(1+x**4)**0.25d0
      rr1=rr1+y*(pgr(i*2-1)**2+pgr(i*2)**2)
    enddo

    call all_sum_all(rr1) 
    rr1=rr1*crys%vcell  

    beta=rr1/rr0
    rr0=rr1


    else

    do i=1,gs%length
      x=gs%ekin(i)/Ek
      y=1/(1+x**4)**0.25d0
      rr1=rr1+y*(pgr(i*2-1)**2+pgr(i*2)**2)
    enddo
   call all_sum_all(rr1) 
    rr1=rr1*crys%vcell  
     beta=dzero
     rr0=rr1

    end if

     ireset = .false.

! GET INITIAL PRECONDITIONED CONJUGATE SEARCH DIRECTION
    do i=1,gs%length
      x=gs%ekin(i)/Ek
      y=1/(1+x**4)**0.25d0
      pgr(i*2)=-pgr(i*2)*y+beta*pgr_old(i*2)
      pgr(i*2-1)=-pgr(i*2-1)*y+beta*pgr_old(i*2-1)
    enddo

!       PRJOECT ONTO SUBSPACE AGAIN
!s=(pgr.ugr)  (Z.U)   
    s=mddot(2*gs%length, pgr(1), 1,ugr(1), 1)
    call all_sum_all(s)     
    s=s*crys%vcell  
!pgr=pgr-s*ugr   Z=Z-(Z.U)U
    call mdaxpy(2*gs%length,-s,ugr(1),1,pgr(1),1)

!pgr_old=pgr  Z_old=Z
    call mdcopy(2*gs%length,pgr(1),1,pgr_old(1),1)

!     NORMALIZE AGAIN
!s=(pgr.pgr)  s=(Z.Z)
    s=mddot(2*gs%length, pgr(1), 1,pgr(1), 1)
    call all_sum_all(s)
    s=1.d0/dsqrt(s*crys%vcell) 
!pgr=s*pgr  Z=(Z.Z)Z
    call mdscal(2*gs%length,s,pgr(1),1)

!*********************************************
!    UPDATE UGR 
!*********************************************

!E_upg=(pgr.pgr)  (Z.Z)
    E_upg=mddot(2*gs%length, pgr(1), 1,ughr(1), 1)
    call all_sum_all(E_upg)
    E_upg=E_upg*crys%vcell  
    E_upg=2*E_upg*crys%ztot

!**************************************************
!ccccc  theta0=0.2D-3, might cause slightly different results for
!ccccc  parallel and serial code due to round off error.
!cccccc for debug, use theta0=0.2D-2,but for actual run, use 0.2D-03.
!**************************************************
    theta0=0.2D-03
!      theta0=0.2D-02

    cos_th=dcos(theta0)
    sin_th=dsin(theta0)

!ugr=cos_th*ugr   STAYING NORMALIZED Unew=cos*U+sin*Z
    call mdscal(2*gs%length,cos_th,ugr(1),1)
!ugr=ugr+sin_th*pgr
    call mdaxpy(2*gs%length,sin_th,pgr(1),1,ugr(1),1)

!ughr=H*Utrial  to get Etot1
    call Hpsi3(ugr,ughr,vion_r,dvp,Etot1,alpha,work1_c,&
                  ffts,crys,gs,den_core)

!   GET OLD U BACK
!ugr=ugr-sin_th*pgr    U=(U-Z*sin)/cos
    call mdaxpy(2*gs%length,-sin_th,pgr(1),1,ugr(1),1)
!ugr=cos_th*ugr
    call mdscal(2*gs%length,1/cos_th,ugr(1),1)

    E_thth=2*(Etot1-Etot-E_upg*theta0)/theta0**2

    theta=0.5d0*dabs(datan(E_upg/E_thth))
    cos_th=dcos(theta)
    sin_th=dsin(theta)

    predE=0.5d0*E_upg*dsin(2*theta)- 0.25d0*E_thth*dcos(2*theta)+ &
        (Etot+0.25d0*E_thth)

500 format(2(f20.15,1x), 2(E10.5,1x))
!**************************************************

!ugr=cos_th*ugr   STAYING NORMALIZED U=cos*U+sin*Z
    call mdscal(2*gs%length,cos_th,ugr(1),1)
!ugr=ugr+sin_th*pgr
    call mdaxpy(2*gs%length,sin_th,pgr(1),1,ugr(1),1)

1000 continue

2000 continue

!*********************************************
!    GET NEW VIN FROM UGR
!*********************************************

    call fourier_transform(-1, ffts, gs, ugr(1), work1_c(1), 1)

    do i=1,gs%r_size
       work1_c(2*i-1)=crys%ztot*work1_c(2*i-1)**2
    enddo

!*******************************************

  call fourier_transform(-1, ffts, gs, den_core(1), work2_c(1), 1)
  call mdscal(2*gs%r_size,dzero, work2_c(1),1)

  do i=1,gs%r_size
    pgr(i)=vion_r(i) +2*UxcCA(work1_c(2*i-1)+work2_c(2*i-1),uxc2)
  enddo

!*******************************************

  call fourier_transform(1, ffts, gs, ugr(1), work1_c(1), 1)

! ADD IN COULOMB
  do i=1,gs%length
    if (gs%ekin(i) .gt.0.d0) then  
      ugr(2*i)=4*ugr(2*i)*2*pi/gs%ekin(i)  
      ugr(2*i-1)=4*ugr(2*i-1)*2*pi/gs%ekin(i) 
!    else           
!      ugr(2*i)=0.d0
!      ugr(2*i-1)=0.d0
    end if
  end do

  call fourier_transform(-1, ffts, gs, ugr(1), work1_c(1), 1)
  call fourier_transform(-1, ffts, gs, vion(1), work2_c(1), 1)

!*******************************************

  call fourier_transform(-1, ffts, gs, vnl_ave(1), ughr(1), 1)

  s=dzero
  do i=1,gs%r_size
    pgr(i)=pgr(i)+work1_c(2*i-1) -ughr(2*i-1)  
    s=s+pgr(i)
  enddo

  call all_sum_all(s)     
  s=s/nr                    ! renormalize 

  do i=1,gs%r_size
    pgr(i)=pgr(i)-s
  enddo

!*******************************************
!*** Now, the pgr is the new v_in potential
!*******************************************

!ugr=vout
  call mdcopy(2*gs%length,vout(1),1,ugr(1),1)
!ughr=vin
  call mdcopy(2*gs%length,vin(1),1,ughr(1),1)

  work1_c=dzero
  call mdcopy(gs%r_size,pgr(1),1,work1_c(1),2)
  call fourier_transform(1, ffts, gs, pgr(1), work1_c(1), 1)

  Ecut0=1.0 
!   write(9,*) y,' yhuh'
!   y=done 
!   if (err .gt. .1) then
!    Ecut0=y*0.05   ! now in Ryd was 0.5 in hartree
!  else  if (err .gt. 1d-2) then
!      Ecut0=y*0.075
!   else  if (err .gt. 1d-3) then
!      Ecut0=y*0.1
!   else  if (err .gt. 1d-4) then
!      Ecut0=y*0.125
!!   else  if (err .gt. 1d-4) then
!!      Ecut0=1.0
!   else
!      Ecut0=y*0.15
!   end if          

  do i=1,gs%length

!    if (err .gt. 1d-2) then
!    x=0.4*gs%ekin(i)/(gs%ekin(i)+0.5)
!    x=0.01*gs%ekin(i)/(gs%ekin(i)+0.5)
!    else
    x=0.8*gs%ekin(i)/(gs%ekin(i)+0.5)
!    end if

    y=dexp(-gs%ekin(i)/Ecut0)
    pgr(i*2)=pgr(i*2)*y+(1-y)*(ughr(i*2)*(1-x)+ugr(i*2)*x)
    pgr(i*2-1)=pgr(i*2-1)*y+(1-y)*(ughr(i*2-1)*(1-x)+ugr(i*2-1)*x)
  enddo

  call mdcopy(gs%length*2,pgr(1),1,vin(1),1) 

 call mdaxpy(2*gs%length,dmone,vion(1),1,vin(1),1)  

  call mdcopy(gs%length*2,pgr(1),1,vion(1),1) 



  call fourier_transform(-1, ffts, gs, pgr(1), work1_c(1), 1)

   
      return
      end

!**************************************************
!         HPSI3
!**************************************************

      subroutine Hpsi3(ugr,ur,vion_r,dvp,Etot,alpha,workr_n,&
                        ffts,crys,gs,den_core)
  use all_to_all_module
  include 'use.h'  
  IMPLICIT NONE
  include 'interface.h' 
  include 'all_to_all.h'  
  include 'flibcalls.ph'  
  !
  !     INPUT:
  !     -----
  ! 
  type(crystal), intent(in) :: crys            ! for vcell, nspin, ztot
  type(parallel_gspace), intent(in) ::  gs  
  type(fft_struc), intent(in) :: ffts 
  real(dp), intent(in) :: den_core(2*gs%length),&
      dvp(gs%r_size) ,&
   vion_r(max(gs%length*2,gs%r_size) ),&
  ugr(max(gs%length*2,gs%r_size) )
  !
  !     WORK:
  !     -----
  ! 
  real(dp)    workr_n(max(gs%length*2,2*gs%r_size) ),&
              ur(max(gs%length*2,gs%r_size) )
  !
  !
  !     --------------- local variables -------------------------
  !
  integer ng2_n,i,nr_n,nr   

  real(dp) fa,uxc2,d,Ekk,UxcCA,vtot,Etot,alpha,s
  real(dp) rho_d(max(gs%length*2,gs%r_size) )
  !
  !
  !     --------------- CODE -------------------------
  !
  ng2_n=gs%length
  nr_n=gs%r_size
  nr=gs%fftsize(1)*gs%fftsize(2)*gs%fftsize(3)
  fa=(3*pi**2)**(2.d0/3.d0)

  call fourier_transform(-1, ffts, gs, ugr(1), workr_n(1), 1)

!  ur = workr_n
  call mdcopy(gs%r_size,workr_n(1),2,ur(1),1)

!***************************************

  workr_n=dzero
  do i=1,gs%r_size
    workr_n(2*i-1)=crys%ztot*ur(i)**2 
  enddo

  call fourier_transform(1, ffts, gs, rho_d(1), workr_n(1), 1)

! put coulomb interaction in pgr
  do i=1,gs%length
    if (gs%ekin(i) .gt. dzero) then 
      rho_d(2*i)=4*rho_d(2*i)*2*pi/gs%ekin(i)  
      rho_d(2*i-1)=4*rho_d(2*i-1)*2*pi/gs%ekin(i)  
    else           
      rho_d(2*i)=dzero
      rho_d(2*i-1)=dzero
    end if
  end do

  call fourier_transform(-1, ffts, gs, rho_d(1), workr_n(1), 1)

!  rho_d = workr_n
  call mdcopy(gs%r_size,workr_n(1),2,rho_d(1),1)

!***************************************

  call fourier_transform(-1, ffts, gs, den_core, workr_n(1), 1)
  call mdscal(2*gs%r_size,dzero, workr_n(1),1)

  Etot=dzero

  workr_n=dzero
  do i=1,gs%r_size

    d=crys%ztot*ur(i)**2
    vtot=vion_r(i)+2*UxcCA(d+workr_n(2*i-1),uxc2)+fa*d**(2.d0/3.d0)+ rho_d(i)

    Etot=Etot+vion_r(i)*d +2*uxc2*(d+workr_n(2*i-1))+fa*(3.d0/5.d0)* &
     d**(5.d0/3.d0) +  0.5d0*rho_d(i)*d &
     +2*dvp(i)*crys%ztot*ur(i)
      
    workr_n(2*i-1)=vtot*ur(i)+dvp(i)
  enddo

  call all_sum_all(Etot)

  Etot=Etot*crys%vcell/nr

  call fourier_transform(1, ffts, gs, ur(1), workr_n(1), 1)

  Ekk=dzero

  do i=1,gs%length

    ur(i*2)=ur(i*2)+alpha*gs%ekin(i)*ugr(i*2)
    ur(i*2-1)=ur(i*2-1)+alpha*gs%ekin(i)*ugr(i*2-1)

    Ekk=Ekk+alpha*gs%ekin(i)*(ugr(i*2-1)**2+ugr(i*2)**2)

  enddo

  call all_sum_all(Ekk)

  Etot=Etot+Ekk*crys%vcell*crys%ztot  

  return
  end

      Function UxcCA(rho,Uxc2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)


      DATA CEX/-1.969490099D0/
      DATA ONE/1.0D0/,HALF/0.5D0/,THALF/1.5D0/
      DATA ZERO/0.0D0/,THRD/0.333333333333333D0/
      DATA TFT/0.75D0/,TWO/2.0D0/
      DATA RC/0.02258D0/,EPS/1.0D-20/
!C-----------------------------------------------------------------------
!C
      RH3=max(rho,1.d-16)

!C.s
!C. KOHN-SHAM EXCHANGE
!C. CEPERLEY-ALDER / RPA CORRELATION
!C.
      DATA GAMMA,BETA1,BETA2 /-0.2846d0, 1.0529d0, 0.3334d0/
      DATA A,B,C,D /0.0622d0,-0.096d0,0.0040d0,-0.0232d0 /
      DATA PI /3.141592654D0/
      DATA BETA11,BETA22 /1.228383333D0,0.44453333333D0/
      DATA B1,C1,D1 /-0.11673333D0,0.0026666667D0,-0.0168D0/
!C
700   RHO=RH3
      RH3=RH3**THRD
      RS=(PI*RHO/TFT)**(-THRD)
      IF(RS .GE. ONE) THEN
       ROOTRS=SQRT(RS)
       VC=GAMMA*(ONE+BETA11*ROOTRS+BETA22*RS) &
         /(ONE+BETA1*ROOTRS+BETA2*RS)**2
       EC=GAMMA/(ONE+BETA1*ROOTRS+BETA2*RS)
      ELSE
       XLNRS=LOG(RS)
       VC=A*XLNRS+B1+C1*RS*XLNRS+D1*RS
       EC=A*XLNRS+B+C*RS*XLNRS+D*RS
      END IF
      VEXCOR=CEX*RH3+VC
      Uxc2=(TFT*CEX*RH3+EC)/TWO
      UxcCA=VEXCOR/TWO
      RETURN
      END


      
