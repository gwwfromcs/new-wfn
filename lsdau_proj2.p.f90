! m4undef.m4
!
! resets various m4 commands to have a prefix m4_
! this can be achieved using option -P but only on some architectures
! any file to be preprocessed with m4 should use this set of macro definitions
!
! David Prendergast, June 6, 2006



! fft_macros.m4
!
!     fft_aux_space(aux,naux,aux2,naux2,dec1,dec2)
!     fft_local_init(N,A,LDA,LOT,DIREC)
!     fft_perm_init(N,AUX)
!
!     fft_multiple_backward(N,A,LOTSTRIDE,LOT,LOOPDUMMY,VECOUT,OFFSET)
!     fft_multiple_forward (N,A,LOTSTRIDE,LOT,LOOPDUMMY,VECOUT,OFFSET)
!
!     fft_backward(N,A,LOTSTRIDE)
!     fft_forward (N,A,LOTSTRIDE)
!
!     fft_convol(N,A,B,POS)
!
!     fft_multiple_scale(N,A,LOTSTRIDE,LOT,LOOPDUMMY,SCALE)
!
!     timeget (t0)
!     timediff(t1,t0)
!
!
!     fft_fcblock     defines blocking factor for fast convolute
!     fft_local_free  release memory
!

!-------------------------------------------------------------------------
        subroutine angular_wfnpz_d2(wavefn,nkpt,kxyz,tau,rc,cart,  &
                          kkx,kky,kkz,vcell,myproc,ngrid,step,Rl,cnk,Rlm,proj_flag)


        use all_to_all_module
        implicit none

      include 'mpif.h'


        integer ngrid,myproc,proj_flag,nkpt,kxyz(3,nkpt)
        double precision tau(3),rc,kkx,kky,kkz,cart(3,3),step,vcell,Rl(ngrid)
        double complex:: wavefn(nkpt),Rlm(ngrid,5),cnk(5)
          
! local   
          
        integer ii,jj,kk,ll,ierr
        double precision bess(4),r1,r2,r3,r4,r5,r6,r7,r8,kx,ky,kz
        double precision k1,k2,k3,kgr,kgr0,k10,k20,k30,ekin0
        double precision SBESSJ,kgdotr0,scale,norm,step2

        double precision ekin,const1,const2,const3
        double complex::  phase,ylm(5),ctem1,ctem2,ione,itwo

        double complex, allocatable::Rlm_sum(:,:)
        double complex, allocatable::Rlm2(:,:)

!============================================================================
!
         
       ione=cmplx(0.d0,1.d0,kind=8)
       itwo=cmplx(0.d0,2.d0,kind=8)
      
       const1=dsqrt(15.d0/(8.d0*pi4))
       const2=dsqrt(15.d0/(2.d0*pi4))
       const3=dsqrt(5.d0/(4.d0*pi4))

       step2=step*step

       allocate(Rlm2(ngrid,5))
      
       Rlm(:,:)=(0.d0,0.d0)

! sum over G points

! Rlm(r)=4\pi i^l e(ikr_0)\sum_G C^k(G) e(iGr_0) j_l(|k+G|r)Y_lm(\Omega_{k+G})
!
! where k is the wavevector in BZ, r_0 is the position of the atom of interest.
!  


           k10 = cart(1,1)*kkx+cart(1,2)*kky+cart(1,3)*kkz
           k20 = cart(2,1)*kkx+cart(2,2)*kky+cart(2,3)*kkz
           k30 = cart(3,1)*kkx+cart(3,2)*kky+cart(3,3)*kkz


      do jj=1,nkpt

           k1=real(kxyz(1,jj),kind=8)
           k2=real(kxyz(2,jj),kind=8)
           k3=real(kxyz(3,jj),kind=8)

           kx = cart(1,1)*k1+cart(1,2)*k2+  &
                cart(1,3)*k3+k10
           ky = cart(2,1)*k1+cart(2,2)*k2+  &
                cart(2,3)*k3+k20
           kz = cart(3,1)*k1+cart(3,2)*k2+  &
                cart(3,3)*k3+k30
           
           kgdotr0=(kx*tau(1)+ky*tau(2)+kz*tau(3))

           ekin = kx*kx+ky*ky+kz*kz
 
           if(ekin.lt.(1.d-10)) then
              ylm(1)=(0.d0,0.d0)
              ylm(2)=(0.d0,0.d0)
              ylm(3)=(0.d0,0.d0)
              ylm(4)=(0.d0,0.d0)
              ylm(5)=(0.d0,0.d0)
           else
!              kx=kx/ekin
!              ky=ky/ekin
!              kz=kz/ekin

!   <Ylm|Psi>=\int Ylm^*\Psi

!              ylm(1)= const21*kx*ky/ekin
!              ylm(2)= const21*ky*kz/ekin
!              ylm(3)= const21*kz*kx/ekin
!              ylm(4)= const21*(kx*kx-ky*ky)/2.d0/ekin
!              ylm(5)= const22*(1.5d0*kz*kz-0.5d0)/ekin


              ylm(1)= const1*(kx*kx-ky*ky+itwo*kx*ky)/ekin       !Y2(-2)^*
              ylm(2)= const2*(kx*kz+ione*ky*kz)/ekin             !Y2(-1)^*
              ylm(3)= const3*(3.d0*kz*kz/ekin-1.d0)              !Y2(0)^*
              ylm(4)= -DCONJG(ylm(2))                            !Y2(1)^*
              ylm(5)= DCONJG(ylm(1))                             !Y2(2)^*
          end if
         

! loop over radial points

          phase=exp(cmplx(0.0d0,kgdotr0,kind=8))
          ctem1=phase*wavefn(jj)
          ekin0=dsqrt(ekin)*step

          do ii=1,ngrid
              kgr=ekin0*real(ii,kind=8)

!              bess(3)=SBESSJ(2,kgr)
!              ctem2=ctem1*bess(3)

              ctem2=ctem1*SBESSJ(2,kgr)

              Rlm(ii,1)=Rlm(ii,1)+ctem2*ylm(1)
              Rlm(ii,2)=Rlm(ii,2)+ctem2*ylm(2)
              Rlm(ii,3)=Rlm(ii,3)+ctem2*ylm(3)
              Rlm(ii,4)=Rlm(ii,4)+ctem2*ylm(4)
              Rlm(ii,5)=Rlm(ii,5)+ctem2*ylm(5)
          end do
       end do

! 4*pi is the Bessel-Fourier transform factor, vcell is the wavefunction
! normalization fortor.
!

       allocate(Rlm_sum(ngrid,5))
        Rlm_sum(:,:)=(0.d0,0.d0)
        do jj=1,5
        call MPI_ALLREDUCE(Rlm(1,jj),Rlm_sum(1,jj),ngrid,MPI_DOUBLE_COMPLEX,  &
         MPI_SUM,MPI_COMM_WORLD,ierr)
        end do
        Rlm(:,:)=Rlm_sum(:,:)
        deallocate(Rlm_sum)

       
       Rlm(:,:)=Rlm(:,:)*pi4/dsqrt(vcell)

       do kk=1,5
       if(proj_flag.eq.1) then
       do ii=1,ngrid
          Rlm2(ii,kk)=Rl(ii)*Rlm(ii,kk)
       end do
       else
       do ii=1,ngrid
          Rlm2(ii,kk)=Rlm(ii,kk)*DCONJG(Rlm(ii,kk))
       end do
       end if

       end do

! do integration in real space from the atomic position to cutoff r_C

       cnk=(0.d0,0.d0)

       do kk=1,5
       do ii=1,ngrid-7,7
          r1=real(ii*ii,kind=8)
          r2=real((ii+1)*(ii+1) ,kind=8)
          r3=real((ii+2)*(ii+2) ,kind=8)
          r4=real((ii+3)*(ii+3) ,kind=8)
          r5=real((ii+4)*(ii+4) ,kind=8)
          r6=real((ii+5)*(ii+5) ,kind=8)
          r7=real((ii+6)*(ii+6) ,kind=8)
          r8=real((ii+7)*(ii+7) ,kind=8)

          cnk(kk)=cnk(kk)+    &
          751.d0*(r1*Rlm2(ii,kk)+ r8*Rlm2(ii+7,kk))+    &
          3577.d0*(r2*Rlm2(ii+1,kk)+r7*Rlm2(ii+6,kk))+    &
          1323.d0*(r3*Rlm2(ii+2,kk)+r6*Rlm2(ii+5,kk))+    &
          2989.d0*(r4*Rlm2(ii+3,kk)+r5*Rlm2(ii+4,kk))
       end do
       end do

      scale=step*step2*(7.d0/17280.d0)

      cnk=cnk*scale
      if(proj_flag.ne.1) cnk=sqrt(cnk)

 101  format(4f15.8)
        
      deallocate(Rlm2)

      return
      end

!-------------------------------------------------------------------------
        subroutine angular_wfnpz_p2(wavefn,nkpt,kxyz,tau,rc,cart,  &
                          kkx,kky,kkz,vcell,myproc,ngrid,step,Rl,cnk,Rlm,proj_flag)



        use all_to_all_module
        implicit none

      include 'mpif.h'


      
        integer ngrid,myproc,proj_flag,nkpt,kxyz(3,nkpt)
        double precision tau(3),rc,kkx,kky,kkz,cart(3,3),step,vcell,Rl(ngrid)
        double complex:: wavefn(nkpt),Rlm(ngrid,3),cnk(3)
          
! local   
          
        integer ii,jj,kk,ll,ierr
        double precision bess(4),r1,r2,r3,r4,r5,r6,r7,r8
        double precision k1,k2,k3,kgr,kgr0,k10,k20,k30,ekin0
        double precision SBESSJ,kgdotr0,scale,norm,step2

        double precision ekin,const1,const0
        double complex::  phase,ylm(3),ione,ctem1,ctem2

        double complex, allocatable::Rlm_sum(:,:)
        double complex, allocatable::Rlm2(:,:)

!============================================================================
!
       ione=(0.d0,1.d0)
       const1=dsqrt(3.d0/(2.d0*pi4))
       const0=dsqrt(3.d0/(pi4))
         
       step2=step*step
       allocate(Rlm2(ngrid,3))
      
       Rlm(:,:)=(0.d0,0.d0)

! sum over G points

! Rlm(r)=4\pi i^l e(ikr_0)\sum_G C^k(G) e(iGr_0) j_l(|k+G|r)Y_lm(\Omega_{k+G})
!
! where k is the wavevector in BZ, r_0 is the position of the atom of interest.
!  

           k10 = cart(1,1)*kkx+cart(1,2)*kky+cart(1,3)*kkz
           k20 = cart(2,1)*kkx+cart(2,2)*kky+cart(2,3)*kkz
           k30 = cart(3,1)*kkx+cart(3,2)*kky+cart(3,3)*kkz


      do jj=1,nkpt

           k1 = cart(1,1)*kxyz(1,jj)+cart(1,2)*kxyz(2,jj)+  &
                cart(1,3)*kxyz(3,jj)+k10
           k2 = cart(2,1)*kxyz(1,jj)+cart(2,2)*kxyz(2,jj)+  &
                cart(2,3)*kxyz(3,jj)+k20
           k3 = cart(3,1)*kxyz(1,jj)+cart(3,2)*kxyz(2,jj)+  &
                cart(3,3)*kxyz(3,jj)+k30

           kgdotr0=(k1*tau(1)+k2*tau(2)+k3*tau(3))

           ekin = dsqrt(k1*k1+k2*k2+k3*k3)
 
           if(ekin.lt.(1.d-10)) then
              ylm(1)=(0.d0,0.d0)
              ylm(2)=(0.d0,0.d0)
              ylm(3)=(0.d0,0.d0)
           else
              k1=k1/ekin
              k2=k2/ekin
              k3=k3/ekin

!   <Ylm|Psi>=\int Ylm^*\Psi

              ylm(1)= const1*(k1+ione*k2)        !Y1(-1)^*
              ylm(2)= const0*k3                  !Y1(0)^*
              ylm(3)= -const1*(k1-ione*k2)       !Y1(1)^*

          end if
         

! loop over radial points

          phase=exp(cmplx(0.0d0,kgdotr0,kind=8))
          ctem1=phase*wavefn(jj)

          ekin0=ekin*step
          do ii=1,ngrid
              kgr=ekin0*ii

!              bess(2)=SBESSJ(1,kgr)

              ctem2=ctem1*SBESSJ(1,kgr)

              Rlm(ii,1)=Rlm(ii,1)+ctem2*ylm(1)
              Rlm(ii,2)=Rlm(ii,2)+ctem2*ylm(2)
              Rlm(ii,3)=Rlm(ii,3)+ctem2*ylm(3)

          end do

       end do

!
! 4*pi is the Bessel-Fourier transform factor, vcell is the wavefunction
! normalization fortor.
!


       allocate(Rlm_sum(ngrid,3))
       Rlm_sum(:,:)=(0.d0,0.d0)

       do jj=1,3
       call MPI_ALLREDUCE(Rlm(1,jj),Rlm_sum(1,jj),ngrid,MPI_DOUBLE_COMPLEX,  &
        MPI_SUM,MPI_COMM_WORLD,ierr)
       end do

       Rlm(:,:)=Rlm_sum(:,:)
       deallocate(Rlm_sum)

       
       Rlm(:,:)=Rlm(:,:)*pi4/dsqrt(vcell)

       do kk=1,3

       if(proj_flag.eq.1) then
       do ii=1,ngrid
          Rlm2(ii,kk)=Rl(ii)*Rlm(ii,kk)
       end do
       else
       do ii=1,ngrid
          Rlm2(ii,kk)=Rlm(ii,kk)*DCONJG(Rlm(ii,kk))
       end do
       end if

       end do


! do integration in real space from the atomic position to cutoff r_C

       cnk=(0.d0,0.d0)

       do kk=1,3
       do ii=1,ngrid-7,7
          r1=ii*ii
          r2=(ii+1)*(ii+1) 
          r3=(ii+2)*(ii+2) 
          r4=(ii+3)*(ii+3) 
          r5=(ii+4)*(ii+4) 
          r6=(ii+5)*(ii+5) 
          r7=(ii+6)*(ii+6) 
          r8=(ii+7)*(ii+7) 

          cnk(kk)=cnk(kk)+    &
          751*(r1*Rlm2(ii,kk)+ r8*Rlm2(ii+7,kk))+    &
          3577*(r2*Rlm2(ii+1,kk)+r7*Rlm2(ii+6,kk))+    &
          1323*(r3*Rlm2(ii+2,kk)+r6*Rlm2(ii+5,kk))+    &
          2989*(r4*Rlm2(ii+3,kk)+r5*Rlm2(ii+4,kk))
       end do
       end do

      scale=step*step2*(7.d0/17280.d0)

      cnk=cnk*scale

 101  format(4f15.8)


! uncomment the following line if use decomposition scheme
! \int |<Ylm|Psi>|^2
! if use
! \int <Ylm*Rl|Psi>
! then do not do sqrt.
!

      if(proj_flag.ne.1) cnk(:)=sqrt(cnk(:))
        
      deallocate(Rlm2)
 

      return
      end


!-------------------------------------------------------------------------
        subroutine angular_wfnpz_s2(wavefn,nkpt,kxyz,tau,rc,cart,  &
                          kkx,kky,kkz,vcell,myproc,ngrid,step,Rl,cnk,Rlm,proj_flag)


        use all_to_all_module
        implicit none

      include 'mpif.h'


	integer ngrid,myproc,proj_flag,nkpt,kxyz(3,nkpt)
        double precision tau(3),rc,kkx,kky,kkz,cart(3,3),step,vcell,Rl(ngrid)
        double complex:: wavefn(nkpt),Rlm(ngrid,1),cnk(1)

! local

	integer ii,jj,kk,ll,ierr
	double precision bess(4),r1,r2,r3,r4,r5,r6,r7,r8
	double precision k1,k2,k3,kgr,kgr0,k10,k20,k30,ekin0
        double precision SBESSJ,kgdotr0,scale,norm,step2

        double precision ekin,const1,const0
	double complex::  phase,ylm(1),ione,ctem1,ctem2

        double complex, allocatable::Rlm_sum(:,:)
        double complex, allocatable::Rlm2(:,:)


!============================================================================
!
       ione=(0.d0,1.d0)
      
       const0=1.d0/dsqrt(pi4)



         
       step2=step*step
       allocate(Rlm2(ngrid,1))
      

       Rlm(:,:)=(0.d0,0.d0)


! sum over G points

! Rlm(r)=4\pi i^l e(ikr_0)\sum_G C^k(G) e(iGr_0) j_l(|k+G|r)Y_lm(\Omega_{k+G})
!
! where k is the wavevector in BZ, r_0 is the position of the atom of interest.
!  

           k10 = cart(1,1)*kkx+cart(1,2)*kky+cart(1,3)*kkz
           k20 = cart(2,1)*kkx+cart(2,2)*kky+cart(2,3)*kkz
           k30 = cart(3,1)*kkx+cart(3,2)*kky+cart(3,3)*kkz
           ylm(1)=const0

      do jj=1,nkpt

           k1 = cart(1,1)*kxyz(1,jj)+cart(1,2)*kxyz(2,jj)+  &
                cart(1,3)*kxyz(3,jj)+k10
           k2 = cart(2,1)*kxyz(1,jj)+cart(2,2)*kxyz(2,jj)+  &
                cart(2,3)*kxyz(3,jj)+k20
           k3 = cart(3,1)*kxyz(1,jj)+cart(3,2)*kxyz(2,jj)+  &
                cart(3,3)*kxyz(3,jj)+k30

           kgdotr0=(k1*tau(1)+k2*tau(2)+k3*tau(3))

           ekin0 = dsqrt(k1*k1+k2*k2+k3*k3)*step

! loop over radial points

          phase=exp(cmplx(0.0d0,kgdotr0,kind=8))
          ctem1=phase*wavefn(jj)

          do ii=1,ngrid
              kgr=ekin0*ii

!              bess(1)=SBESSJ(0,kgr)

              ctem2=ctem1*SBESSJ(0,kgr)

              Rlm(ii,1)=Rlm(ii,1)+ctem2*ylm(1)

          end do

       end do

!
! 4*pi is the Bessel-Fourier transform factor, vcell is the wavefunction
! normalization fortor.
!

       allocate(Rlm_sum(ngrid,1))
       Rlm_sum(:,:)=(0.d0,0.d0)
       jj=1
       call MPI_ALLREDUCE(Rlm(1,jj),Rlm_sum(1,jj),ngrid,MPI_DOUBLE_COMPLEX,  &
        MPI_SUM,MPI_COMM_WORLD,ierr)
       Rlm(:,:)=Rlm_sum(:,:)
       deallocate(Rlm_sum)

       
       Rlm(:,:)=Rlm(:,:)*pi4/dsqrt(vcell)

       kk=1
       if(proj_flag.eq.1) then
       do ii=1,ngrid
          Rlm2(ii,kk)=Rl(ii)*Rlm(ii,kk)
       end do
       else
       do ii=1,ngrid
          Rlm2(ii,kk)=Rlm(ii,kk)*DCONJG(Rlm(ii,kk))
       end do
       end if


! do integration in real space from the atomic position to cutoff r_C

       cnk=(0.d0,0.d0)

       kk=1
       do ii=1,ngrid-7,7
          r1=ii*ii
          r2=(ii+1)*(ii+1) 
          r3=(ii+2)*(ii+2) 
          r4=(ii+3)*(ii+3) 
          r5=(ii+4)*(ii+4) 
          r6=(ii+5)*(ii+5) 
          r7=(ii+6)*(ii+6) 
          r8=(ii+7)*(ii+7) 

          cnk(kk)=cnk(kk)+    &
          751*(r1*Rlm2(ii,kk)+ r8*Rlm2(ii+7,kk))+    &
          3577*(r2*Rlm2(ii+1,kk)+r7*Rlm2(ii+6,kk))+    &
          1323*(r3*Rlm2(ii+2,kk)+r6*Rlm2(ii+5,kk))+    &
          2989*(r4*Rlm2(ii+3,kk)+r5*Rlm2(ii+4,kk))
       end do

      scale=step*step2*(7.d0/17280.d0)

      cnk=cnk*scale


! uncomment the following line if use decomposition scheme
! \int |<Ylm|Psi>|^2
! if use
! \int <Ylm*Rl|Psi>
! then do not do sqrt.
!

      if(proj_flag.ne.1) cnk(:)=sqrt(cnk(:))
        
      deallocate(Rlm2)
 

      return
      end

        subroutine angular_wfnpz(wavefn,nkpt,kxyz,  &
                          tau,rc,cart,kkx,kky,kkz,vcell,ngrid,myproc)


        use all_to_all_module
        implicit none

      include 'mpif.h'


        integer nkpt,kxyz(3,nkpt),ngrid
        double precision tau(3),rc,kkx,kky,kkz,cart(3,3),vcell
        complex (kind=8) wavefn(nkpt)
                
        integer ii,jj,kk,ll,ierr,myproc
        double precision k1,k2,k3,kgr,kgr0
        double precision ylm(9),bess(4),r1,r2,r3,r4,r5,r6,r7,r8
        double precision ekin,const0,const1,const21,const22
        double precision SBESSJ,kgdotr0,scale,step,norm
        double precision, allocatable::Rlm2(:,:)
        complex (kind=8)  phase,alpha(9)
        complex (kind=8), allocatable::Rlm(:,:),Rlm_sum(:,:)



        const0=1.d0/dsqrt(pi4)
        const1=dsqrt(3.d0)*const0
        const21=dsqrt(15.d0)*const0
        const22=dsqrt(5.d0)*const0

!============================================================================
!       ngrid=150

       step=rc/ngrid
       allocate(Rlm(ngrid,9))
       allocate(Rlm2(ngrid,9))

       Rlm(:,:)=(0.d0,0.d0)

! sum over G points

! Rlm(r)=4\pi i^l e(ikr_0)\sum_G C^k(G) e(iGr_0) j_l(|k+G|r)Y_lm(\Omega_{k+G})
!
! where k is the wavevector in BZ, r_0 is the position of the atom of interest.
!

      do jj=1,nkpt

           k1 = cart(1,1)*(kxyz(1,jj)+kkx)+cart(1,2)*(kxyz(2,jj)+kky)+  &
                cart(1,3)*(kxyz(3,jj)+kkz)
           k2 = cart(2,1)*(kxyz(1,jj)+kkx)+cart(2,2)*(kxyz(2,jj)+kky)+  &
                cart(2,3)*(kxyz(3,jj)+kkz)
           k3 = cart(3,1)*(kxyz(1,jj)+kkx)+cart(3,2)*(kxyz(2,jj)+kky)+  &
                cart(3,3)*(kxyz(3,jj)+kkz)


           kgdotr0=(k1*tau(1)+k2*tau(2)+k3*tau(3))

           ekin = dsqrt(k1*k1+k2*k2+k3*k3)

           if(ekin.lt.(1.d-10)) then
              ylm(1)=const0
              ylm(2)=0.d0
              ylm(3)=0.d0
              ylm(4)=0.d0
              ylm(5)=0.d0
              ylm(6)=0.d0
              ylm(7)=0.d0
              ylm(8)=0.d0
              ylm(9)=0.d0
           else
              k1=k1/ekin
              k2=k2/ekin
              k3=k3/ekin

              ylm(1)=const0
              ylm(2)=const1*k1
              ylm(3)=const1*k2
              ylm(4)=const1*k3
              ylm(5)= const21*k1*k2
              ylm(6)= const21*k2*k3
              ylm(7)= const21*K3*k1
              ylm(8)= const21*(k1*k1-k2*k2)/2.d0
              ylm(9)= const22*(1.5d0*k3*k3-0.5d0)
          end if

! loop over radial points

          do ii=1,ngrid
              kgr=ekin*step*ii
              phase=exp(cmplx(0.0d0,kgdotr0))

              bess(1)=SBESSJ(0,kgr)
              bess(2)=SBESSJ(1,kgr)
              bess(3)=SBESSJ(2,kgr)
!        write(239,*) sin(kgr)/kgr,(sin(kgr)/(kgr*kgr)-cos(kgr)/kgr),
!     $ (3.d0/(kgr*kgr*kgr)-1.d0/kgr)*sin(kgr)-3.d0*cos(kgr)/(kgr*kgr)
              Rlm(ii,1)=Rlm(ii,1)+wavefn(jj)*bess(1)*phase*ylm(1)
              Rlm(ii,2)=Rlm(ii,2)+wavefn(jj)*bess(2)*phase*ylm(2)
              Rlm(ii,3)=Rlm(ii,3)+wavefn(jj)*bess(2)*phase*ylm(3)
              Rlm(ii,4)=Rlm(ii,4)+wavefn(jj)*bess(2)*phase*ylm(4)
              Rlm(ii,5)=Rlm(ii,5)+wavefn(jj)*bess(3)*phase*ylm(5)
              Rlm(ii,6)=Rlm(ii,6)+wavefn(jj)*bess(3)*phase*ylm(6)
              Rlm(ii,7)=Rlm(ii,7)+wavefn(jj)*bess(3)*phase*ylm(7)
              Rlm(ii,8)=Rlm(ii,8)+wavefn(jj)*bess(3)*phase*ylm(8)
              Rlm(ii,9)=Rlm(ii,9)+wavefn(jj)*bess(3)*phase*ylm(9)

          end do

       end do

!
! 4*pi is the Bessel-Fourier transform factor, vcell is the wavefunction
! normalization fortor.
!

       allocate(Rlm_sum(ngrid,9))
       Rlm_sum(:,:)=(0.d0,0.d0)
       do jj=1,9
       call MPI_REDUCE(Rlm(1,jj),Rlm_sum(1,jj),ngrid,MPI_DOUBLE_COMPLEX,  &
        MPI_SUM,0,MPI_COMM_WORLD,ierr)
       end do
       if(myproc.eq.0) then
       Rlm(:,:)=Rlm_sum(:,:)
       end if
       deallocate(Rlm_sum)



       if(myproc.eq.0) then


       Rlm(:,:)=Rlm(:,:)*pi4/dsqrt(vcell)

       do jj=1,9
       do ii=1,ngrid
          Rlm2(ii,kk)=Rlm(ii,kk)*DCONJG(Rlm(ii,kk))
       end do

       end do

! do integration in real space from the atomic position to cutoff r_C

       do jj=1,9
          alpha(jj)=0.d0
       do ii=1,ngrid-7,7
          r1=step*ii*step*ii
          r2=step*(ii+1)*step*(ii+1)
          r3=step*(ii+2)*step*(ii+2)
          r4=step*(ii+3)*step*(ii+3)
          r5=step*(ii+4)*step*(ii+4)
          r6=step*(ii+5)*step*(ii+5)
          r7=step*(ii+6)*step*(ii+6)
          r8=step*(ii+7)*step*(ii+7)

          alpha(jj)=alpha(jj)+    &
          751*(r1*Rlm2(ii,jj)+ r8*Rlm2(ii+7,jj))+    &
          3577*(r2*Rlm2(ii+1,jj)+r7*Rlm2(ii+6,jj))+    &
          1323*(r3*Rlm2(ii+2,jj)+r6*Rlm2(ii+5,jj))+    &
          2989*(r4*Rlm2(ii+3,jj)+r5*Rlm2(ii+4,jj))
       end do
       end do

      scale=step*(7.d0/17280.d0)

      alpha(:)=alpha(:)*scale
      norm=0.d0


      do jj=1,9
         norm=norm+alpha(jj)
      end do

      do jj=1,9
         write(239,*) jj,alpha(jj),alpha(jj)/norm
      end do
      end if

      deallocate(Rlm)
      deallocate(Rlm2)


      return
      end


!-------------------------------------------------------------------------
! real spherical harmonic
!
        subroutine angular_wfnpz_f2r(wavefn,nkpt,kxyz,tau,rc,cart,  &
                          kkx,kky,kkz,vcell,myproc,ngrid,step,Rl,cnk,Rlm,proj_flag)


        use all_to_all_module
        implicit none

      include 'mpif.h'


        integer ngrid,myproc,proj_flag,nkpt,kxyz(3,nkpt)
        double precision tau(3),rc,kkx,kky,kkz,cart(3,3),step,vcell,Rl(ngrid)
        double complex:: wavefn(nkpt),Rlm(ngrid,7),cnk(7)
          
! local   
          
        integer ii,jj,kk,ll,ierr
        double precision bess(4),r1,r2,r3,r4,r5,r6,r7,r8,kx,ky,kz,kx2,ky2,kz2
        double precision kkxyz,kxy2,kxz2,kyz2,kyx2,kzx2,kzy2,kx3,ky3,kz3
        double precision k1,k2,k3,kgr,kgr0,k10,k20,k30,ekin0
        double precision SBESSJ,kgdotr0,scale,norm,step2,const31,const32,const33,const34

        double precision ekin,const1,const2,const3
        double complex::  phase,ylm(7),ctem1,ctem2,ione,itwo

        double complex, allocatable::Rlm_sum(:,:)
        double complex, allocatable::Rlm2(:,:)

!============================================================================
!
         
       ione=cmplx(0.d0,1.d0,kind=8)
       itwo=cmplx(0.d0,2.d0,kind=8)
      
       const31=dsqrt(7.d0/(16.d0*pi4))
       const32=dsqrt(35.d0/(32.d0*pi4))
       const33=dsqrt(105.d0/(16.d0*pi4))
       const34=dsqrt(21.d0/(32.d0*pi4))

       step2=step*step

       allocate(Rlm2(ngrid,7))
      
       Rlm(:,:)=(0.d0,0.d0)

! sum over G points

! Rlm(r)=4\pi i^l e(ikr_0)\sum_G C^k(G) e(iGr_0) j_l(|k+G|r)Y_lm(\Omega_{k+G})
!
! where k is the wavevector in BZ, r_0 is the position of the atom of interest.
!  


           k10 = cart(1,1)*kkx+cart(1,2)*kky+cart(1,3)*kkz
           k20 = cart(2,1)*kkx+cart(2,2)*kky+cart(2,3)*kkz
           k30 = cart(3,1)*kkx+cart(3,2)*kky+cart(3,3)*kkz


      do jj=1,nkpt

           k1=real(kxyz(1,jj),kind=8)
           k2=real(kxyz(2,jj),kind=8)
           k3=real(kxyz(3,jj),kind=8)

           kx = cart(1,1)*k1+cart(1,2)*k2+  &
                cart(1,3)*k3+k10
           ky = cart(2,1)*k1+cart(2,2)*k2+  &
                cart(2,3)*k3+k20
           kz = cart(3,1)*k1+cart(3,2)*k2+  &
                cart(3,3)*k3+k30
           
           kgdotr0=(kx*tau(1)+ky*tau(2)+kz*tau(3))

           ekin = dsqrt(kx*kx+ky*ky+kz*kz)
           kx=kx/ekin
           ky=ky/ekin
           kz=kz/ekin

           kx2=kx*kx
           ky2=ky*ky
           kz2=kz*kz
         
           kx3=kx2*kx
           ky3=ky2*ky
           kz3=kz2*kz
 
           kkxyz=kx*ky*kz
           kxy2=kx*ky2
           kyx2=ky*kx2
           kzx2=kz*kx2
           kxz2=kx*kz2
           kzy2=kz*ky2
           kyz2=ky*kz2
           
           
           if(ekin.lt.(1.d-10)) then
              ylm(1)=(0.d0,0.d0)
              ylm(2)=(0.d0,0.d0)
              ylm(3)=(0.d0,0.d0)
              ylm(4)=(0.d0,0.d0)
              ylm(5)=(0.d0,0.d0)
              ylm(6)=(0.d0,0.d0)
              ylm(7)=(0.d0,0.d0)
           else
!              kx=kx/ekin
!              ky=ky/ekin
!              kz=kz/ekin

!   <Ylm|Psi>=\int Ylm^*\Psi

              ylm(1)= const31*(2*kz3-3*kzx2-3*kzy2)
              ylm(2)= const32*(3*kyx2-ky3)
              ylm(3)= const32*(kx3-3*kxy2)
              ylm(4)= const33*(kzx2-kzy2)
              ylm(5)= const33*kkxyz
              ylm(6)= const34*(4*kyz2-kyx2-ky3)
              ylm(7)= const34*(4*kxz2-kx3-kxy2)

          end if
         

! loop over radial points

          phase=exp(cmplx(0.0d0,kgdotr0,kind=8))
          ctem1=phase*wavefn(jj)
          ekin0=dsqrt(ekin)*step

          do ii=1,ngrid
              kgr=ekin0*real(ii,kind=8)

!              bess(3)=SBESSJ(2,kgr)
!              ctem2=ctem1*bess(3)

              ctem2=ctem1*SBESSJ(3,kgr)

              Rlm(ii,1)=Rlm(ii,1)+ctem2*ylm(1)
              Rlm(ii,2)=Rlm(ii,2)+ctem2*ylm(2)
              Rlm(ii,3)=Rlm(ii,3)+ctem2*ylm(3)
              Rlm(ii,4)=Rlm(ii,4)+ctem2*ylm(4)
              Rlm(ii,5)=Rlm(ii,5)+ctem2*ylm(5)
              Rlm(ii,6)=Rlm(ii,6)+ctem2*ylm(6)
              Rlm(ii,7)=Rlm(ii,7)+ctem2*ylm(7)
          end do
       end do

! 4*pi is the Bessel-Fourier transform factor, vcell is the wavefunction
! normalization fortor.
!

       allocate(Rlm_sum(ngrid,7))
        Rlm_sum(:,:)=(0.d0,0.d0)
        do jj=1,7
        call MPI_ALLREDUCE(Rlm(1,jj),Rlm_sum(1,jj),ngrid,MPI_DOUBLE_COMPLEX,  &
         MPI_SUM,MPI_COMM_WORLD,ierr)
        end do
        Rlm(:,:)=Rlm_sum(:,:)
        deallocate(Rlm_sum)

       
       Rlm(:,:)=Rlm(:,:)*pi4/dsqrt(vcell)

       do kk=1,7
       if(proj_flag.eq.1) then
       do ii=1,ngrid
          Rlm2(ii,kk)=Rl(ii)*Rlm(ii,kk)
       end do
       else
       do ii=1,ngrid
          Rlm2(ii,kk)=Rlm(ii,kk)*DCONJG(Rlm(ii,kk))
       end do
       end if

       end do

! do integration in real space from the atomic position to cutoff r_C

       cnk=(0.d0,0.d0)

       do kk=1,7
       do ii=1,ngrid-7,7
          r1=real(ii*ii,kind=8)
          r2=real((ii+1)*(ii+1) ,kind=8)
          r3=real((ii+2)*(ii+2) ,kind=8)
          r4=real((ii+3)*(ii+3) ,kind=8)
          r5=real((ii+4)*(ii+4) ,kind=8)
          r6=real((ii+5)*(ii+5) ,kind=8)
          r7=real((ii+6)*(ii+6) ,kind=8)
          r8=real((ii+7)*(ii+7) ,kind=8)

          cnk(kk)=cnk(kk)+    &
          751.d0*(r1*Rlm2(ii,kk)+ r8*Rlm2(ii+7,kk))+    &
          3577.d0*(r2*Rlm2(ii+1,kk)+r7*Rlm2(ii+6,kk))+    &
          1323.d0*(r3*Rlm2(ii+2,kk)+r6*Rlm2(ii+5,kk))+    &
          2989.d0*(r4*Rlm2(ii+3,kk)+r5*Rlm2(ii+4,kk))
       end do
       end do

      scale=step*step2*(7.d0/17280.d0)

      cnk=cnk*scale
      if(proj_flag.ne.1) cnk=sqrt(cnk)

 101  format(4f15.8)
        
      deallocate(Rlm2)

      return
      end

!-------------------------------------------------------------------------
! real spherical harmonic
!
        subroutine angular_wfnpz_d2r(wavefn,nkpt,kxyz,tau,rc,cart,  &
                          kkx,kky,kkz,vcell,myproc,ngrid,step,Rl,cnk,Rlm,proj_flag)


        use all_to_all_module
        implicit none

      include 'mpif.h'


        integer ngrid,myproc,proj_flag,nkpt,kxyz(3,nkpt)
        double precision tau(3),rc,kkx,kky,kkz,cart(3,3),step,vcell,Rl(ngrid)
        double complex:: wavefn(nkpt),Rlm(ngrid,5),cnk(5)
          
! local   
          
        integer ii,jj,kk,ll,ierr
        double precision bess(4),r1,r2,r3,r4,r5,r6,r7,r8,kx,ky,kz
        double precision k1,k2,k3,kgr,kgr0,k10,k20,k30,ekin0
        double precision SBESSJ,kgdotr0,scale,norm,step2

        double precision ekin,const1,const2,const3
        double complex::  phase,ylm(5),ctem1,ctem2,ione,itwo

        double complex, allocatable::Rlm_sum(:,:)
        double complex, allocatable::Rlm2(:,:)

!============================================================================
!
         
       ione=cmplx(0.d0,1.d0,kind=8)
       itwo=cmplx(0.d0,2.d0,kind=8)
      
       const1=dsqrt(5.d0/(4.d0*pi4)) 
       const2=dsqrt(15.d0/pi4)
       const3=dsqrt(15.d0/(4.d0*pi4))

       step2=step*step

       allocate(Rlm2(ngrid,5))
      
       Rlm(:,:)=(0.d0,0.d0)

! sum over G points

! Rlm(r)=4\pi i^l e(ikr_0)\sum_G C^k(G) e(iGr_0) j_l(|k+G|r)Y_lm(\Omega_{k+G})
!
! where k is the wavevector in BZ, r_0 is the position of the atom of interest.
!  


           k10 = cart(1,1)*kkx+cart(1,2)*kky+cart(1,3)*kkz
           k20 = cart(2,1)*kkx+cart(2,2)*kky+cart(2,3)*kkz
           k30 = cart(3,1)*kkx+cart(3,2)*kky+cart(3,3)*kkz


      do jj=1,nkpt

           k1=real(kxyz(1,jj),kind=8)
           k2=real(kxyz(2,jj),kind=8)
           k3=real(kxyz(3,jj),kind=8)

           kx = cart(1,1)*k1+cart(1,2)*k2+  &
                cart(1,3)*k3+k10
           ky = cart(2,1)*k1+cart(2,2)*k2+  &
                cart(2,3)*k3+k20
           kz = cart(3,1)*k1+cart(3,2)*k2+  &
                cart(3,3)*k3+k30
           
           kgdotr0=(kx*tau(1)+ky*tau(2)+kz*tau(3))

           ekin = kx*kx+ky*ky+kz*kz
 
           if(ekin.lt.(1.d-10)) then
              ylm(1)=(0.d0,0.d0)
              ylm(2)=(0.d0,0.d0)
              ylm(3)=(0.d0,0.d0)
              ylm(4)=(0.d0,0.d0)
              ylm(5)=(0.d0,0.d0)
           else
!              kx=kx/ekin
!              ky=ky/ekin
!              kz=kz/ekin

!   <Ylm|Psi>=\int Ylm^*\Psi

              ylm(1)= const2*kx*ky/ekin         ! i/sqrt(2)Y(-2)-Y(2)
              ylm(2)= const2*ky*kz/ekin         ! i/sqrt(2)Y(-1)+Y(1)
              ylm(3)= const2*kz*kx/ekin         ! 1/sqrt(2)Y(-1)-Y(1)
              ylm(4)= const3*(kx*kx-ky*ky)/ekin ! 1/sqrt(2)Y(-2)+Y(2)
              ylm(5)= const1*(3.d0*kz*kz/ekin-1.0d0)  ! Y(0)


!              ylm(1)= const1*(kx*kx-ky*ky+itwo*kx*ky)/ekin       !Y2(-2)^*
!              ylm(2)= const2*(kx*kz+ione*ky*kz)/ekin             !Y2(-1)^*
!              ylm(3)= const3*(3.d0*kz*kz/ekin-1.d0)              !Y2(0)^*
!              ylm(4)= -DCONJG(ylm(2))                            !Y2(1)^*
!              ylm(5)= DCONJG(ylm(1))                             !Y2(2)^*
          end if
         

! loop over radial points

          phase=exp(cmplx(0.0d0,kgdotr0,kind=8))
          ctem1=phase*wavefn(jj)
          ekin0=dsqrt(ekin)*step

          do ii=1,ngrid
              kgr=ekin0*real(ii,kind=8)

!              bess(3)=SBESSJ(2,kgr)
!              ctem2=ctem1*bess(3)

              ctem2=ctem1*SBESSJ(2,kgr)

              Rlm(ii,1)=Rlm(ii,1)+ctem2*ylm(1)
              Rlm(ii,2)=Rlm(ii,2)+ctem2*ylm(2)
              Rlm(ii,3)=Rlm(ii,3)+ctem2*ylm(3)
              Rlm(ii,4)=Rlm(ii,4)+ctem2*ylm(4)
              Rlm(ii,5)=Rlm(ii,5)+ctem2*ylm(5)
          end do
       end do

! 4*pi is the Bessel-Fourier transform factor, vcell is the wavefunction
! normalization fortor.
!

       allocate(Rlm_sum(ngrid,5))
        Rlm_sum(:,:)=(0.d0,0.d0)
        do jj=1,5
        call MPI_ALLREDUCE(Rlm(1,jj),Rlm_sum(1,jj),ngrid,MPI_DOUBLE_COMPLEX,  &
         MPI_SUM,MPI_COMM_WORLD,ierr)
        end do
        Rlm(:,:)=Rlm_sum(:,:)
        deallocate(Rlm_sum)

       
       Rlm(:,:)=Rlm(:,:)*pi4/dsqrt(vcell)

       do kk=1,5
       if(proj_flag.eq.1) then
       do ii=1,ngrid
          Rlm2(ii,kk)=Rl(ii)*Rlm(ii,kk)
       end do
       else
       do ii=1,ngrid
          Rlm2(ii,kk)=Rlm(ii,kk)*DCONJG(Rlm(ii,kk))
       end do
       end if

       end do

! do integration in real space from the atomic position to cutoff r_C

       cnk=(0.d0,0.d0)

       do kk=1,5
       do ii=1,ngrid-7,7
          r1=real(ii*ii,kind=8)
          r2=real((ii+1)*(ii+1) ,kind=8)
          r3=real((ii+2)*(ii+2) ,kind=8)
          r4=real((ii+3)*(ii+3) ,kind=8)
          r5=real((ii+4)*(ii+4) ,kind=8)
          r6=real((ii+5)*(ii+5) ,kind=8)
          r7=real((ii+6)*(ii+6) ,kind=8)
          r8=real((ii+7)*(ii+7) ,kind=8)

          cnk(kk)=cnk(kk)+    &
          751.d0*(r1*Rlm2(ii,kk)+ r8*Rlm2(ii+7,kk))+    &
          3577.d0*(r2*Rlm2(ii+1,kk)+r7*Rlm2(ii+6,kk))+    &
          1323.d0*(r3*Rlm2(ii+2,kk)+r6*Rlm2(ii+5,kk))+    &
          2989.d0*(r4*Rlm2(ii+3,kk)+r5*Rlm2(ii+4,kk))
       end do
       end do

      scale=step*step2*(7.d0/17280.d0)

      cnk=cnk*scale
      if(proj_flag.ne.1) cnk=sqrt(cnk)

 101  format(4f15.8)
        
      deallocate(Rlm2)

      return
      end

!-------------------------------------------------------------------------
! real spherical harmonic
        subroutine angular_wfnpz_p2r(wavefn,nkpt,kxyz,tau,rc,cart,  &
                          kkx,kky,kkz,vcell,myproc,ngrid,step,Rl,cnk,Rlm,proj_flag)



        use all_to_all_module
        implicit none

      include 'mpif.h'


      
        integer ngrid,myproc,proj_flag,nkpt,kxyz(3,nkpt)
        double precision tau(3),rc,kkx,kky,kkz,cart(3,3),step,vcell,Rl(ngrid)
        double complex:: wavefn(nkpt),Rlm(ngrid,3),cnk(3)
          
! local   
          
        integer ii,jj,kk,ll,ierr
        double precision bess(4),r1,r2,r3,r4,r5,r6,r7,r8
        double precision k1,k2,k3,kgr,kgr0,k10,k20,k30,ekin0
        double precision SBESSJ,kgdotr0,scale,norm,step2

        double precision ekin,const1,const0
        double complex::  phase,ylm(3),ione,ctem1,ctem2

        double complex, allocatable::Rlm_sum(:,:)
        double complex, allocatable::Rlm2(:,:)

!============================================================================
!
       ione=(0.d0,1.d0)
       const0=dsqrt(3.d0/(pi4))
         
       step2=step*step
       allocate(Rlm2(ngrid,3))
      
       Rlm(:,:)=(0.d0,0.d0)

! sum over G points

! Rlm(r)=4\pi i^l e(ikr_0)\sum_G C^k(G) e(iGr_0) j_l(|k+G|r)Y_lm(\Omega_{k+G})
!
! where k is the wavevector in BZ, r_0 is the position of the atom of interest.
!  

           k10 = cart(1,1)*kkx+cart(1,2)*kky+cart(1,3)*kkz
           k20 = cart(2,1)*kkx+cart(2,2)*kky+cart(2,3)*kkz
           k30 = cart(3,1)*kkx+cart(3,2)*kky+cart(3,3)*kkz


      do jj=1,nkpt

           k1 = cart(1,1)*kxyz(1,jj)+cart(1,2)*kxyz(2,jj)+  &
                cart(1,3)*kxyz(3,jj)+k10
           k2 = cart(2,1)*kxyz(1,jj)+cart(2,2)*kxyz(2,jj)+  &
                cart(2,3)*kxyz(3,jj)+k20
           k3 = cart(3,1)*kxyz(1,jj)+cart(3,2)*kxyz(2,jj)+  &
                cart(3,3)*kxyz(3,jj)+k30

           kgdotr0=(k1*tau(1)+k2*tau(2)+k3*tau(3))

           ekin = dsqrt(k1*k1+k2*k2+k3*k3)
 
           if(ekin.lt.(1.d-10)) then
              ylm(1)=(0.d0,0.d0)
              ylm(2)=(0.d0,0.d0)
              ylm(3)=(0.d0,0.d0)
           else
              k1=k1/ekin
              k2=k2/ekin
              k3=k3/ekin

! real spherical harmonic

              ylm(1)= const0*k1
              ylm(2)= const0*k2                  
              ylm(3)= const0*k3

!   <Ylm|Psi>=\int Ylm^*\Psi
!              ylm(1)= const1*(k1+ione*k2)        !Y1(-1)^*
!              ylm(2)= const0*k3                  !Y1(0)^*
!              ylm(3)= -const1*(k1-ione*k2)       !Y1(1)^*

          end if
         

! loop over radial points

          phase=exp(cmplx(0.0d0,kgdotr0,kind=8))
          ctem1=phase*wavefn(jj)

          ekin0=ekin*step
          do ii=1,ngrid
              kgr=ekin0*ii


              bess(2)=SBESSJ(1,kgr)

              ctem2=ctem1*bess(2)

              Rlm(ii,1)=Rlm(ii,1)+ctem2*ylm(1)
              Rlm(ii,2)=Rlm(ii,2)+ctem2*ylm(2)
              Rlm(ii,3)=Rlm(ii,3)+ctem2*ylm(3)

          end do

       end do

!
! 4*pi is the Bessel-Fourier transform factor, vcell is the wavefunction
! normalization fortor.
!


       allocate(Rlm_sum(ngrid,3))
       Rlm_sum(:,:)=(0.d0,0.d0)

       do jj=1,3
       call MPI_ALLREDUCE(Rlm(1,jj),Rlm_sum(1,jj),ngrid,MPI_DOUBLE_COMPLEX,  &
        MPI_SUM,MPI_COMM_WORLD,ierr)
       end do

       Rlm(:,:)=Rlm_sum(:,:)
       deallocate(Rlm_sum)

       
       Rlm(:,:)=Rlm(:,:)*pi4/dsqrt(vcell)

       do kk=1,3

       if(proj_flag.eq.1) then
       do ii=1,ngrid
          Rlm2(ii,kk)=Rl(ii)*Rlm(ii,kk)
       end do
       else
       do ii=1,ngrid
          Rlm2(ii,kk)=Rlm(ii,kk)*DCONJG(Rlm(ii,kk))
       end do
       end if

       end do


! do integration in real space from the atomic position to cutoff r_C

       cnk=(0.d0,0.d0)

       do kk=1,3
       do ii=1,ngrid-7,7
          r1=ii*ii
          r2=(ii+1)*(ii+1) 
          r3=(ii+2)*(ii+2) 
          r4=(ii+3)*(ii+3) 
          r5=(ii+4)*(ii+4) 
          r6=(ii+5)*(ii+5) 
          r7=(ii+6)*(ii+6) 
          r8=(ii+7)*(ii+7) 

          cnk(kk)=cnk(kk)+    &
          751*(r1*Rlm2(ii,kk)+ r8*Rlm2(ii+7,kk))+    &
          3577*(r2*Rlm2(ii+1,kk)+r7*Rlm2(ii+6,kk))+    &
          1323*(r3*Rlm2(ii+2,kk)+r6*Rlm2(ii+5,kk))+    &
          2989*(r4*Rlm2(ii+3,kk)+r5*Rlm2(ii+4,kk))
       end do
       end do

      scale=step*step2*(7.d0/17280.d0)

      cnk=cnk*scale

 101  format(4f15.8)


! uncomment the following line if use decomposition scheme
! \int |<Ylm|Psi>|^2
! if use
! \int <Ylm*Rl|Psi>
! then do not do sqrt.
!

      if(proj_flag.ne.1) cnk(:)=sqrt(cnk(:))
        
      deallocate(Rlm2)
 

      return
      end

