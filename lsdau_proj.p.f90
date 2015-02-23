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

        subroutine lsdau_proj_init(k_gspace,kpoints,crys)

        use lsdau_shared

        use constants
        use crystal_module
        use parallel_gspace_module
        use kpoint_module
        implicit none

        type(kpoint), intent(in) :: kpoints  
        type(parallel_gspace), intent(in) :: k_gspace(*)
        type(crystal) :: crys       ! crystal structure

!-------------------------------------------------------------------------

        integer::it,i2,i1,i0,G1,G2,G3,ia,im
        integer::len_gk_max,ngrid_max
	integer ii,jj,kk,irk

        double precision::k_1,k_2,k_3,k_x,k_y,k_z,kgr,kg_x,kg_y,kg_z,ekin,ekin0
        double precision const21,const22,const23,const11,const12,const0,kgr0
        double precision SBESSJ

        double precision,allocatable::tau(:,:)

	double complex:: ione,itwo

!-------------------------------------------------------------
       i2=0
       i1=0
       i0=0
       len_gk_max=MAXVAL(k_gspace(1:kpoints%nrk)%length)
       ngrid_max=MAXVAL(ngrid1)

       allocate(tau(3,Nsite))

       kk=0
       do it=1,crys%ntype
          if(lproj_tem(it).eq.2) then
             i2=1
             do ia=1, crys%natom(it)
                 kk=kk+1
                 tau(:,kk)=matmul(crys%avec,crys%rat(:,ia,it))/pi2
             end do
           end if

          if(lproj_tem(it).eq.1) then
             i1=1
             do ia=1, crys%natom(it)
                 kk=kk+1
                 tau(:,kk)=matmul(crys%avec,crys%rat(:,ia,it))/pi2
             end do
          end if

          if(lproj_tem(it).eq.0) then
             i0=1
             do ia=1, crys%natom(it)
                 kk=kk+1
                 tau(:,kk)=matmul(crys%avec,crys%rat(:,ia,it))/pi2
             end do
          end if
       end do


       if(i2.eq.1) then
         allocate(ylm2(5,len_gk_max,kpoints%nrk))

! this matrix could be huge
         allocate(sb2(ngrid_max,len_gk_max,kpoints%nrk))
       end if
 
       if(i1.eq.1) then
         allocate(ylm1(3,len_gk_max,kpoints%nrk))
         allocate(sb1(ngrid_max,len_gk_max,kpoints%nrk))
       end if

       if(i0.eq.1) then
         allocate(sb0(ngrid_max,len_gk_max,kpoints%nrk))
       end if

       allocate(phase_ikgr0(len_gk_max,kpoints%nrk,Nsite))
       allocate(ikg(3,len_gk_max,kpoints%nrk))


!---------------------------------------------------------------
        ione=cmplx(0.d0,1.d0,kind=8)
        itwo=cmplx(0.d0,2.d0,kind=8)
      
        const21=dsqrt(15.d0/(8.d0*pi4))
        const22=dsqrt(15.d0/(2.d0*pi4))
        const23=dsqrt(5.d0/(4.d0*pi4))

        const12=dsqrt(3.d0/(2.d0*pi4))
        const11=dsqrt(3.d0/(pi4))

        const0=1.d0/dsqrt(pi4)

        ylm0=const0


!------------------------------------------------------------------------------

! Rlm(r)=4*pi i^l e(ikr_0)\sum_G C^k(G) e(iGr_0) j_l(|k+G|r)Y_lm(\Omega_{k+G})
!
! where k is the wavevector in BZ, r_0 is the position of the atom of interest.
! 
! here we calculate j_l and Y_lm and store them
!  
!------------------------------------------------------------------------------

        do irk = 1, kpoints%nrk
           k_1=kpoints%rk(1,irk)
           k_2=kpoints%rk(2,irk)
           k_3=kpoints%rk(3,irk)

           k_x = crys%bvec(1,1)*k_1+crys%bvec(1,2)*k_2+crys%bvec(1,3)*k_3
           k_y = crys%bvec(2,1)*k_1+crys%bvec(2,2)*k_2+crys%bvec(2,3)*k_3
           k_z = crys%bvec(3,1)*k_1+crys%bvec(3,2)*k_2+crys%bvec(3,3)*k_3

        do jj=1,k_gspace(irk)%length

           G1= k_gspace(irk)%gvec(1,jj)
           G2= k_gspace(irk)%gvec(2,jj)
           G3= k_gspace(irk)%gvec(3,jj)

           kg_x = crys%bvec(1,1)*G1+crys%bvec(1,2)*G2+crys%bvec(1,3)*G3+k_x
           kg_y = crys%bvec(2,1)*G1+crys%bvec(2,2)*G2+crys%bvec(2,3)*G3+k_y
           kg_z = crys%bvec(3,1)*G1+crys%bvec(3,2)*G2+crys%bvec(3,3)*G3+k_z

           ekin = kg_x*kg_x+kg_y*kg_y+kg_z*kg_z

           ikg(1,jj,irk)=ione*kg_x
           ikg(2,jj,irk)=ione*kg_y
           ikg(3,jj,irk)=ione*kg_z

           do ia=1,Nsite
              kgr0= kg_x*tau(1,ia)+kg_y*tau(2,ia)+kg_z*tau(3,ia)
              phase_ikgr0(jj,irk,ia)=exp(cmplx(0.0d0,kgr0,kind=8))
           end do

!--------------------------------------------------------------------------------------- 
!        ylm1(3,len_gk_max,kpoints%nrk)

           if(ekin.lt.(1.d-10)) then

             if(i1.eq.1) then
              do im=1,3
                 ylm1(im,jj,irk)=(0.d0,0.d0)
              end do
             end if

             if(i2.eq.1) then
              do im=1,5
                 ylm2(im,jj,irk)=(0.d0,0.d0)
              end do
             end if

           else
             
              if(i2.eq.1) then
              ylm2(1,jj,irk)= const21*(kg_x*kg_x-kg_y*kg_y+itwo*kg_x*kg_y)/ekin       !Y2(-2)^*
              ylm2(2,jj,irk)= const22*(kg_x*kg_z+ione*kg_y*kg_z)/ekin             !Y2(-1)^*
              ylm2(3,jj,irk)= const23*(3.d0*kg_z*kg_z/ekin-1.d0)              !Y2(0)^*
              ylm2(4,jj,irk)= -DCONJG(ylm2(2,jj,irk))                            !Y2(1)^*
              ylm2(5,jj,irk)= DCONJG(ylm2(1,jj,irk))                             !Y2(2)^*
              end if

              if(i1.eq.1) then
              ylm1(1,jj,irk)= const12*(kg_x+ione*kg_y)/dsqrt(ekin)        !Y1(-1)^*
              ylm1(2,jj,irk)= const11*kg_z/dsqrt(ekin)                  !Y1(0)^*
              ylm1(3,jj,irk)= -const12*(kg_x-ione*kg_y)/dsqrt(ekin)       !Y1(1)^*
              end if

          end if
!--------------------------------------------------------------------------------------- 
         
!         sb1(ngrid_max,len_gk_max,kpoints%nrk)

! loop over radial points


          ekin0=dsqrt(ekin)*step

          do ii=1,ngrid_max
             kgr=ekin0*ii

             if(i2.eq.1) sb2(ii,jj,irk)=SBESSJ(2,kgr)
             if(i1.eq.1) sb1(ii,jj,irk)=SBESSJ(1,kgr)
             if(i0.eq.1) sb0(ii,jj,irk)=SBESSJ(0,kgr)
          end do
        end do
        end do

        deallocate(tau)

        return
        end
!-------------------------------------------------------------------------


        subroutine angular_wfnpz_dr(wavefn,k_gspace,crys,kpoints,idr)

        use lsdau_shared
        include 'use.h'
        implicit none


      include 'mpif.h'



	integer::idr

        type(crystal),intent(in) :: crys       ! crystal structure
        type(kpoint), intent(in) :: kpoints  
        type(parallel_gspace), intent(in) :: k_gspace(*)
        type(complex_gspace_array), intent(in) ::  wavefn

!-------------------------------------------------------------------------


        integer:: ia,irk,is,ib,ig,im,ii,len_g,twolp1,ierr,ix

        double precision::scale

	double precision r1,r2,r3,r4,r5,r6,r7,r8

        double complex::ione,itwo,ctem1,ctem2,ctem3(3)
        double complex, allocatable::Rlm_sum(:),Rlm(:,:),Rlm2(:,:),cnk(:)
        double complex, allocatable::drRlm(:,:,:),drRlm2(:,:,:),drcnk(:,:)


!-------------------------------------------------------------



        scale=step*step*step*(7.d0/17280.d0)*pi4/dsqrt(crys%vcell)

        ione=cmplx(0.d0,1.d0,kind=8)
        itwo=cmplx(0.d0,2.d0,kind=8)
      
!------------------------------------------------------------------------------
! sum over G points

! Rlm(r)=4\pi i^l \sum_G C^k(G) ei(k+G)r_0 j_l(|k+G|r)Y_lm(\Omega_{k+G})
!
! where k is the wavevector in BZ, r_0 is the position of the atom of interest.
!
!------------------------------------------------------------------------------

!       cnk_LSDA(1:5,ib,irk,is,ia)=cnk(1:5)

        do ia=1,Nsite

        twolp1=2*lproj(ia)+1

        allocate(Rlm(ngrid(ia),twolp1))
        allocate(Rlm2(ngrid(ia),twolp1))
        allocate(cnk(twolp1))
        if(idr.eq.1) then
          allocate(drcnk(twolp1,3))
          allocate(drRlm(ngrid(ia),twolp1,3))
          allocate(drRlm2(ngrid(ia),twolp1,3))
        end if

        do is=1,crys%nspin
        do irk = 1, kpoints%nrk

        len_g=k_gspace(irk)%length

        do ib=nband0+1,nband1

        Rlm=(0.d0,0.d0)
        if(idr.eq.1) drRlm=(0.d0,0.d0)

!*******************************
        if(lproj(ia).eq.2) then
        do ig=1,len_g

           ctem1=wavefn%data((ib-1)*len_g+ig,irk,is)* phase_ikgr0(ig,irk,ia) 

           do im=1,twolp1
           ctem2=ctem1*ylm2(im,ig,irk)

           do ii=1,ngrid(ia)
                Rlm(ii,im)=Rlm(ii,im)+ctem2*sb2(ii,ig,irk)
           end do

          
           if(idr.eq.1) then
           do ix=1,3
           ctem3(ix)=ctem2*ikg(ix,ig,irk)

           do ii=1,ngrid(ia)
                drRlm(ii,im,ix)=drRlm(ii,im,ix)+ctem3(ix)*sb2(ii,ig,irk)
           end do
           end do
           end if
           end do

        end do
        end if

        if(lproj(ia).eq.1) then
        do ig=1,len_g

           ctem1=wavefn%data((ib-1)*len_g+ig,irk,is)* phase_ikgr0(ig,irk,ia) 

           do im=1,twolp1
           ctem2=ctem1*ylm1(im,ig,irk)

           do ii=1,ngrid(ia)
                Rlm(ii,im)=Rlm(ii,im)+ctem2*sb1(ii,ig,irk)
           end do

          
           if(idr.eq.1) then
           do ix=1,3
           ctem3(ix)=ctem2*ikg(ix,ig,irk)

           do ii=1,ngrid(ia)
                drRlm(ii,im,ix)=drRlm(ii,im,ix)+ctem3(ix)*sb1(ii,ig,irk)
           end do
           end do
           end if
           end do

        end do
        end if

        if(lproj(ia).eq.0) then
        do ig=1,len_g

           ctem1=wavefn%data((ib-1)*len_g+ig,irk,is)* phase_ikgr0(ig,irk,ia) 

           do im=1,twolp1
           ctem2=ctem1*ylm0

           do ii=1,ngrid(ia)
                Rlm(ii,im)=Rlm(ii,im)+ctem2*sb0(ii,ig,irk)
           end do

          
           if(idr.eq.1) then
           do ix=1,3
           ctem3(ix)=ctem2*ikg(ix,ig,irk)

           do ii=1,ngrid(ia)
                drRlm(ii,im,ix)=drRlm(ii,im,ix)+ctem3(ix)*sb0(ii,ig,irk)
           end do
           end do
           end if
           end do

        end do
        end if



! collect all results



        allocate(Rlm_sum(ngrid(ia)))

        call mpi_barrier(MPI_COMM_WORLD, ierr)
        do im=1,twolp1

        call MPI_REDUCE(Rlm(1,im),Rlm_sum,ngrid(ia),MPI_DOUBLE_COMPLEX,  &
        MPI_SUM,0,MPI_COMM_WORLD,ierr)

        if(myproc.eq.0) then
        Rlm(1:ngrid(ia),im)=Rlm_sum(1:ngrid(ia))
        end if

        end do


        if(idr.eq.1) then
        do ix=1,3
        do im=1,twolp1
        call MPI_REDUCE(drRlm(1,im,ix),Rlm_sum,ngrid(ia),MPI_DOUBLE_COMPLEX,  &
        MPI_SUM,0,MPI_COMM_WORLD,ierr)


        if(myproc.eq.0) then
        drRlm(1:ngrid(ia),im,ix)=Rlm_sum(1:ngrid(ia))
        end if
        end do
        end do
        end if

        deallocate(Rlm_sum)


        if(myproc.eq.0) then

        do im=1,twolp1
        do ii=1,ngrid(ia)
           Rlm2(ii,im)=Rnl_intp(ii,ntp(ia))*Rlm(ii,im)
        end do
        end do

        if(idr.eq.1) then
        do ix=1,3
        do im=1,twolp1
        do ii=1,ngrid(ia)
           drRlm2(ii,im,ix)=Rnl_intp(ii,ntp(ia))*drRlm(ii,im,ix)
        end do
        end do
        end do
        end if

! do integration in real-space

        cnk(:)=(0.d0,0.d0)
        if(idr.eq.1) drcnk=(0.d0,0.d0)

        do im=1,twolp1
        do ii=1,ngrid(ia)-7,7
           r1=ii*ii
           r2=(ii+1)*(ii+1)
           r3=(ii+2)*(ii+2)
           r4=(ii+3)*(ii+3)
           r5=(ii+4)*(ii+4)
           r6=(ii+5)*(ii+5)
           r7=(ii+6)*(ii+6)
           r8=(ii+7)*(ii+7)

           cnk(im)=cnk(im)+    &
            751.d0*(r1*Rlm2(ii, im)+ r8*Rlm2(ii+7,im))+    &
           3577.d0*(r2*Rlm2(ii+1,im)+r7*Rlm2(ii+6,im))+    &
           1323.d0*(r3*Rlm2(ii+2,im)+r6*Rlm2(ii+5,im))+    &
           2989.d0*(r4*Rlm2(ii+3,im)+r5*Rlm2(ii+4,im))

           if(idr.eq.1) then
           do ix=1,3
           drcnk(im,ix)=drcnk(im,ix)+    &
           751.d0*(r1* drRlm2(ii,  im,ix)+r8*drRlm2(ii+7,im,ix))+    &
           3577.d0*(r2*drRlm2(ii+1,im,ix)+r7*drRlm2(ii+6,im,ix))+    &
           1323.d0*(r3*drRlm2(ii+2,im,ix)+r6*drRlm2(ii+5,im,ix))+    &
           2989.d0*(r4*drRlm2(ii+3,im,ix)+r5*drRlm2(ii+4,im,ix))
           end do
           end if

 
        end do
        end do
        
        cnk=cnk*scale
        if(idr.eq.1) drcnk=drcnk*scale

        end if ! myproc.eq.1



      call mpi_barrier(MPI_COMM_WORLD, ierr)
      call mpi_bcast(cnk,twolp1,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
      if(idr.eq.1) then
        do ix=1,3
        call mpi_bcast(drcnk(1,ix),twolp1,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
        end do
      end if


        cnk_LSDA(1:twolp1,ib-nband0,irk,is,ia)=cnk(1:twolp1)
        if(idr.eq.1) then
           do ix=1,3
           drcnk_LSDA(1:twolp1,ib-nband0,irk,is,ia,ix)=drcnk(1:twolp1,ix)
           end do
        end if


!*******************************

        end do ! lopp over band
        end do ! k-point
        end do ! spin
        deallocate(Rlm)
        deallocate(Rlm2)
        deallocate(cnk)
        if(idr.eq.1) then
        deallocate(drRlm)
        deallocate(drRlm2)
        deallocate(drcnk)
        end if

        end do ! atomic sites


      return
      end



       REAL*8 FUNCTION SBESSJ(N,X)
       IMPLICIT REAL*8 (A-H,O-Z)
!
!
!
       PARAMETER(ONE=1.0d0,TWO=2.0d0,THREE=3.0d0,ZERO=0.0d0)
       PARAMETER(FIVE = 5.0d0,TEN = 10.0d0,FOURTN = 14.0d0)
!
!
!      SPHERICAL BESSEL FUNCTION OF THE FIRST KIND
!
       IF(ABS(X) .GT. 0.001) THEN
         SB0 = SIN(X)/X
       ELSE
         X2 = X*X/TWO
         SB0 = ONE - (X2/THREE)*(ONE - X2/TEN)
       ENDIF
       IF(N .EQ. 0) THEN
         SBESSJ = SB0
       ELSE
         IF(ABS(X) .GT. 0.001) THEN
           SB1 = (SIN(X)/X - COS(X)) / X
         ELSE
           X2 = X*X/TWO
           SB1 = (X/THREE)*(ONE - (X2/FIVE)*(1.0d0 - X2/FOURTN))
         ENDIF
         IF(N .EQ. 1) THEN
           SBESSJ = SB1
         ELSEIF(X .EQ. ZERO) THEN
           SBESSJ = ZERO
         ELSE
           BY = SB1
           BYM = SB0
           UX = ONE / X
           DO 10 J=1,N-1
             BYP = REAL(2*J+1)*UX*BY - BYM
             BYM = BY
             BY = BYP
 10        CONTINUE
           SBESSJ = BY
         ENDIF
       ENDIF
       RETURN
       END


!-------------------------------------------------------------------------
        subroutine angular_wfnpz_p(wavefn,nkpt,kxyz,tau,rc,cart,  &
                          kkx,kky,kkz,vcell,myproc,ngrid,step,Rl,cnk, &
                          cnk_f, Rlm,proj_flag,syms,reduced)

        use all_to_all_module
        use symmetry_module

        implicit none

      include 'mpif.h'


	integer ii,jj,kk,ll,ngrid,ierr,myproc,proj_flag,ifrac,im
        integer nkpt,kxyz(3,nkpt)
        logical::reduced
        double precision tau(3),rc,kkx,kky,kkz,cart(3,3),vcell,Rl(ngrid)
	complex (kind=8) wavefn(nkpt)
        double complex::Rlm(ngrid,3)

	double precision k1,k2,k3,kx,ky,kz,kgr

        double complex::cnk(5),cnk_f(5,48)

	double precision r1,r2,r3,r4,r5,r6,r7,r8,k10,k20,k30,ekin0
        double precision ekin,const0,const1
        double precision SBESSJ,kgdotr0,scale,step,norm,step2,gdotts
	double complex::  phase,ylm(3),ione,itwo,ctem1

        double complex, allocatable::Rlm_sum(:,:)
        double complex, allocatable::Rlm2(:,:)
        double complex, allocatable::Rlmf(:,:,:),Rlm2f(:,:,:)
        double complex, allocatable::ctem2(:),ctem3(:),ctemf(:)

        type(symmetry),intent(in):: syms


        ifrac=0
        if(syms%nfrac.ne.0) then
          ifrac=1
          allocate(Rlmf(ngrid,3,syms%nfrac))
          allocate(Rlm2f(ngrid,3,syms%nfrac))
          allocate(ctemf(ngrid))
          allocate(ctem3(ngrid))
          Rlmf=(0.d0,0.d0)
        end if

        ione=cmplx(0.d0,1.d0,kind=8)
        itwo=cmplx(0.d0,2.d0,kind=8)
      

        const1=dsqrt(3.d0/(2.d0*pi4))
        const0=dsqrt(3.d0/(pi4))

!============================================================================
!
       allocate(ctem2(ngrid))
         
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

           k1=kxyz(1,jj)
           k2=kxyz(2,jj)
           k3=kxyz(3,jj)

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
           else
!   <Ylm|Psi>=\int Ylm^*\Psi

              ylm(1)= const1*(kx+ione*ky)/dsqrt(ekin)        !Y1(-1)^*
              ylm(2)= const0*kz/dsqrt(ekin)                  !Y1(0)^*
              ylm(3)= -const1*(kx-ione*ky)/dsqrt(ekin)       !Y1(1)^*

          end if
         

! loop over radial points

          phase=exp(cmplx(0.0d0,kgdotr0,kind=8))
          ctem1=phase*wavefn(jj)
          ekin0=dsqrt(ekin)*step

          do ii=1,ngrid
             kgr=ekin0*ii
             ctem2(ii)=ctem1*SBESSJ(1,kgr)
          end do

          do ii=1,ngrid
             Rlm(ii,1)=Rlm(ii,1)+ctem2(ii)*ylm(1)
             Rlm(ii,2)=Rlm(ii,2)+ctem2(ii)*ylm(2)
             Rlm(ii,3)=Rlm(ii,3)+ctem2(ii)*ylm(3)
          end do

!----------------------------------------
! for systems with fractional translations

          if(ifrac.eq.1) then

          do kk=1,syms%nfrac
             gdotts=k1*syms%tnp(1,syms%indfrac(kk))+k2*syms%tnp(2,syms%indfrac(kk))+ &
                    k3*syms%tnp(3,syms%indfrac(kk))
             ctemf(kk)=exp(cmplx(0.0d0,gdotts,kind=8))
          end do

          do kk=1,syms%nfrac

             ctem3(:)=ctem2(:)*ctemf(kk)

          do ii=1,ngrid
             do im=1,3
             Rlmf(ii,im,kk)=Rlmf(ii,im,kk)+ctem3(ii)*ylm(im)
             end do
          end do
          end do

          end if
!----------------------------------------

       end do

! 4*pi is the Bessel-Fourier transform factor, vcell is the wavefunction
! normalization fortor.
!

       allocate(Rlm_sum(ngrid,3))
        Rlm_sum(:,:)=(0.d0,0.d0)
        do jj=1,3
        call MPI_ALLREDUCE(Rlm(1,jj),Rlm_sum(1,jj),ngrid,MPI_DOUBLE_COMPLEX,  &
         MPI_SUM,MPI_COMM_WORLD,ierr)
        end do
        Rlm(:,1:3)=Rlm_sum(:,1:3)

! systems with fractional translations
        if(ifrac.eq.1) then

        do ii=1,syms%nfrac
           Rlm_sum(:,:)=(0.d0,0.d0)
           do jj=1,3
              call MPI_ALLREDUCE(Rlmf(1,jj,ii),Rlm_sum(1,jj),ngrid,MPI_DOUBLE_COMPLEX,  &
                                 MPI_SUM,MPI_COMM_WORLD,ierr)
           end do
           Rlmf(:,:,ii)=Rlm_sum(:,:)
        end do
        end if

        deallocate(Rlm_sum)

       
       Rlm(:,1:3)=Rlm(:,1:3)*pi4/dsqrt(vcell)
       if(ifrac.eq.1) Rlmf(:,:,:)=Rlmf(:,:,:)*pi4/dsqrt(vcell)

!---------------------------------------------------
       if(proj_flag.eq.1) then

       do kk=1,3
       do ii=1,ngrid
          Rlm2(ii,kk)=Rl(ii)*Rlm(ii,kk)
       end do
       end do


       if(ifrac.eq.1) then

         do jj=1,syms%nfrac
         do kk=1,3
            do ii=1,ngrid
               Rlm2f(ii,kk,jj)=Rl(ii)*Rlmf(ii,kk,jj)
            end do
         end do
         end do

       end if

       else

       do kk=1,3
       do ii=1,ngrid
          Rlm2(ii,kk)=Rlm(ii,kk)*DCONJG(Rlm(ii,kk))
       end do
       end do

       if(ifrac.eq.1) then
         do jj=1,syms%nfrac
         do kk=1,3
            do ii=1,ngrid
               Rlm2f(ii,kk,jj)=Rlmf(ii,kk,jj)*DCONJG(Rlmf(ii,kk,jj))
            end do
         end do
         end do
       end if

       end if

!---------------------------------------------------

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
          751.d0*(r1*Rlm2(ii,kk)+ r8*Rlm2(ii+7,kk))+    &
          3577.d0*(r2*Rlm2(ii+1,kk)+r7*Rlm2(ii+6,kk))+    &
          1323.d0*(r3*Rlm2(ii+2,kk)+r6*Rlm2(ii+5,kk))+    &
          2989.d0*(r4*Rlm2(ii+3,kk)+r5*Rlm2(ii+4,kk))
       end do
       end do

       if(ifrac.eq.1) then

       cnk_f=(0.d0,0.d0)

       do jj=1,syms%nfrac

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

          cnk_f(kk,jj)=cnk_f(kk,jj)+    &
          751.d0*(r1*Rlm2f(ii,kk,jj)+ r8*Rlm2f(ii+7,kk,jj))+    &
          3577.d0*(r2*Rlm2f(ii+1,kk,jj)+r7*Rlm2f(ii+6,kk,jj))+    &
          1323.d0*(r3*Rlm2f(ii+2,kk,jj)+r6*Rlm2f(ii+5,kk,jj))+    &
          2989.d0*(r4*Rlm2f(ii+3,kk,jj)+r5*Rlm2f(ii+4,kk,jj))
       end do
       end do
       end do
       end if

      scale=step*step2*(7.d0/17280.d0)

      cnk=cnk*scale

      if(ifrac.eq.1) cnk_f=cnk_f*scale

      if(proj_flag.ne.1) then
         cnk=sqrt(cnk)
         if(ifrac.eq.1) cnk_f=sqrt(cnk_f)
      end if
        
      deallocate(Rlm2)
      deallocate(ctem2)

      if(ifrac.eq.1) then
        deallocate(Rlmf)
        deallocate(Rlm2f)
        deallocate(ctem3)
        deallocate(ctemf)
      end if


      return
      end




!-------------------------------------------------------------------------
        subroutine angular_wfnpz_s(wavefn,nkpt,kxyz,tau,rc,cart,  &
                          kkx,kky,kkz,vcell,myproc,ngrid,step,Rl,cnk,Rlm,proj_flag)


        use all_to_all_module
        implicit none

      include 'mpif.h'


	integer ii,jj,kk,ll,ngrid,ierr,myproc,proj_flag
        integer nkpt,kxyz(3,nkpt)
        double precision tau(3),rc,kkx,kky,kkz,cart(3,3),vcell,Rl(ngrid)
	complex (kind=8) wavefn(nkpt)
        double complex::Rlm(ngrid,1)

	double precision k1,k2,k3,kgr,kgr0,k10,k20,k30,ekin0

        double complex::cnk(1)

	double precision r1,r2,r3,r4,r5,r6,r7,r8
        double precision ekin,const1,const0
        double precision SBESSJ,kgdotr0,scale,step,norm,step2
	double complex::  phase,ylm(1),ione,ctem1,ctem2

!        double complex, allocatable::Rlm(:,:),Rlm_sum(:,:)
        double complex, allocatable::Rlm_sum(:,:)
        double complex, allocatable::Rlm2(:,:)

        ione=(0.d0,1.d0)
      
       const0=1.d0/dsqrt(pi4)



!============================================================================
!
         
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

!-------------------------------------------------------------------------
        subroutine angular_wfnpz_d(wavefn,nkpt,kxyz,tau,rc,cart,  &
                          kkx,kky,kkz,vcell,myproc,ngrid,step,Rl,cnk, &
                          cnk_f, Rlm,proj_flag,syms,reduced)

        use all_to_all_module
        use symmetry_module

        implicit none

      include 'mpif.h'


	integer ii,jj,kk,ll,ngrid,ierr,myproc,proj_flag,ifrac,im
        integer nkpt,kxyz(3,nkpt)
        logical::reduced
        double precision tau(3),rc,kkx,kky,kkz,cart(3,3),vcell,Rl(ngrid)
	complex (kind=8) wavefn(nkpt)
        double complex::Rlm(ngrid,5)

	double precision k1,k2,k3,kx,ky,kz,kgr

        double complex::cnk(5),cnk_f(5,48)

	double precision r1,r2,r3,r4,r5,r6,r7,r8,k10,k20,k30,ekin0
        double precision ekin,const1,const2,const3,const0,const21,const22
        double precision SBESSJ,kgdotr0,scale,step,norm,step2,gdotts
	double complex::  phase,ylm(5),ione,itwo,ctem1,ctem4,ctem11(3)

        double complex, allocatable::Rlm_sum(:,:)
        double complex, allocatable::Rlm2(:,:)
        double complex, allocatable::Rlmf(:,:,:),Rlm2f(:,:,:)
        double complex, allocatable::ctem2(:),ctem3(:),ctemf(:),ctem22(:,:)

        type(symmetry),intent(in):: syms


        ifrac=0
        if(syms%nfrac.ne.0) then
          ifrac=1
          allocate(Rlmf(ngrid,5,syms%nfrac))
          allocate(Rlm2f(ngrid,5,syms%nfrac))
          allocate(ctemf(ngrid))
          allocate(ctem3(ngrid))
          Rlmf=(0.d0,0.d0)
        end if

        ione=cmplx(0.d0,1.d0,kind=8)
        itwo=cmplx(0.d0,2.d0,kind=8)
      
        const1=dsqrt(15.d0/(8.d0*pi4))
        const2=dsqrt(15.d0/(2.d0*pi4))
        const3=dsqrt(5.d0/(4.d0*pi4))

!============================================================================
!
       allocate(ctem2(ngrid))
         
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

           k1=kxyz(1,jj)
           k2=kxyz(2,jj)
           k3=kxyz(3,jj)

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
             kgr=ekin0*ii
             ctem2(ii)=ctem1*SBESSJ(2,kgr)
          end do

          do ii=1,ngrid
             do im=1,5
             Rlm(ii,im)=Rlm(ii,im)+ctem2(ii)*ylm(im)
             end do
          end do

!----------------------------------------
! for systems with fractional translations

          if(ifrac.eq.1) then

          do kk=1,syms%nfrac
             gdotts=k1*syms%tnp(1,syms%indfrac(kk))+k2*syms%tnp(2,syms%indfrac(kk))+ &
                    k3*syms%tnp(3,syms%indfrac(kk))
             ctemf(kk)=exp(cmplx(0.0d0,gdotts,kind=8))
          end do

          do kk=1,syms%nfrac

             ctem3(:)=ctem2(:)*ctemf(kk)

          do ii=1,ngrid
             Rlmf(ii,1,kk)=Rlmf(ii,1,kk)+ctem3(ii)*ylm(1)
             Rlmf(ii,2,kk)=Rlmf(ii,2,kk)+ctem3(ii)*ylm(2)
             Rlmf(ii,3,kk)=Rlmf(ii,3,kk)+ctem3(ii)*ylm(3)
             Rlmf(ii,4,kk)=Rlmf(ii,4,kk)+ctem3(ii)*ylm(4)
             Rlmf(ii,5,kk)=Rlmf(ii,5,kk)+ctem3(ii)*ylm(5)
          end do
          end do

          end if
!----------------------------------------

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


! systems with fractional translations
        if(ifrac.eq.1) then

        do ii=1,syms%nfrac
           Rlm_sum(:,:)=(0.d0,0.d0)
           do jj=1,5
              call MPI_ALLREDUCE(Rlmf(1,jj,ii),Rlm_sum(1,jj),ngrid,MPI_DOUBLE_COMPLEX,  &
                                 MPI_SUM,MPI_COMM_WORLD,ierr)
           end do
           Rlmf(:,:,ii)=Rlm_sum(:,:)
        end do
        end if

        deallocate(Rlm_sum)

       
       Rlm(:,:)=Rlm(:,:)*pi4/dsqrt(vcell)
       if(ifrac.eq.1) Rlmf(:,:,:)=Rlmf(:,:,:)*pi4/dsqrt(vcell)

!---------------------------------------------------
       if(proj_flag.eq.1) then

       do kk=1,5
       do ii=1,ngrid
          Rlm2(ii,kk)=Rl(ii)*Rlm(ii,kk)
       end do
       end do


       if(ifrac.eq.1) then

         do jj=1,syms%nfrac
         do kk=1,5
            do ii=1,ngrid
               Rlm2f(ii,kk,jj)=Rl(ii)*Rlmf(ii,kk,jj)
            end do
         end do
         end do

       end if

       else

       do kk=1,5
       do ii=1,ngrid
          Rlm2(ii,kk)=Rlm(ii,kk)*DCONJG(Rlm(ii,kk))
       end do
       end do

       if(ifrac.eq.1) then
         do jj=1,syms%nfrac
         do kk=1,5
            do ii=1,ngrid
               Rlm2f(ii,kk,jj)=Rlmf(ii,kk,jj)*DCONJG(Rlmf(ii,kk,jj))
            end do
         end do
         end do
       end if

       end if

!---------------------------------------------------

! do integration in real space from the atomic position to cutoff r_C

       cnk=(0.d0,0.d0)

       do kk=1,5
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
          751.d0*(r1*Rlm2(ii,kk)+ r8*Rlm2(ii+7,kk))+    &
          3577.d0*(r2*Rlm2(ii+1,kk)+r7*Rlm2(ii+6,kk))+    &
          1323.d0*(r3*Rlm2(ii+2,kk)+r6*Rlm2(ii+5,kk))+    &
          2989.d0*(r4*Rlm2(ii+3,kk)+r5*Rlm2(ii+4,kk))


       end do
       end do

       if(ifrac.eq.1) then

       cnk_f=(0.d0,0.d0)

       do jj=1,syms%nfrac

       do kk=1,5
       do ii=1,ngrid-7,7
          r1=ii*ii
          r2=(ii+1)*(ii+1)
          r3=(ii+2)*(ii+2)
          r4=(ii+3)*(ii+3)
          r5=(ii+4)*(ii+4)
          r6=(ii+5)*(ii+5)
          r7=(ii+6)*(ii+6)
          r8=(ii+7)*(ii+7)

          cnk_f(kk,jj)=cnk_f(kk,jj)+    &
          751.d0*(r1*Rlm2f(ii,kk,jj)+ r8*Rlm2f(ii+7,kk,jj))+    &
          3577.d0*(r2*Rlm2f(ii+1,kk,jj)+r7*Rlm2f(ii+6,kk,jj))+    &
          1323.d0*(r3*Rlm2f(ii+2,kk,jj)+r6*Rlm2f(ii+5,kk,jj))+    &
          2989.d0*(r4*Rlm2f(ii+3,kk,jj)+r5*Rlm2f(ii+4,kk,jj))
       end do
       end do
       end do
       end if

      scale=step*step2*(7.d0/17280.d0)

      cnk=cnk*scale

      if(ifrac.eq.1) cnk_f=cnk_f*scale

      if(proj_flag.ne.1) then
         cnk=sqrt(cnk)
         if(ifrac.eq.1) cnk_f=sqrt(cnk_f)
      end if
        
      deallocate(Rlm2)
      deallocate(ctem2)

      if(ifrac.eq.1) then
        deallocate(Rlmf)
        deallocate(Rlm2f)
        deallocate(ctem3)
        deallocate(ctemf)
      end if


      return
      end


        subroutine angular_wfnpz_dw(wavefn,nkpt,kxyz,tau,rc,cart,  &
                          kkx,kky,kkz,vcell,myproc,ngrid,step,Rl,cnk,Rlm,proj_flag)

        use all_to_all_module
        implicit none

      include 'mpif.h'


	integer ii,jj,kk,ll,ngrid,ierr,myproc,proj_flag
        integer nkpt,kxyz(3,nkpt)
        double precision tau(3),rc,kkx,kky,kkz,cart(3,3),vcell,Rl(ngrid)
	complex (kind=8) wavefn(nkpt)
        double complex::Rlm(ngrid,5)

	double precision k1,k2,k3,kgr,kgr0,kx,ky,kz

        double complex::cnk(5)

	double precision bess(4),r1,r2,r3,r4,r5,r6,r7,r8,k10,k20,k30,ekin0
        double precision ekin,const1,const2,const3,const0,const21,const22
        double precision SBESSJ,kgdotr0,scale,step,norm,step2
	double complex::  phase,ylm(5),ione,itwo,ctem1,ctem2

        double complex, allocatable::Rlm_sum(:,:)
        double complex, allocatable::Rlm2(:,:)

        ione=cmplx(0.d0,1.d0,kind=8)
        itwo=cmplx(0.d0,2.d0,kind=8)
      
        const1=dsqrt(15.d0/(8.d0*pi4))
        const2=dsqrt(15.d0/(2.d0*pi4))
        const3=dsqrt(5.d0/(4.d0*pi4))

!        const0=1.d0/dsqrt(pi4)
!        const1=dsqrt(3.d0)*const0
!        const21=dsqrt(15.d0)*const0
!        const22=dsqrt(5.d0)*const0

!============================================================================
!
         
       call myflush(240)
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
        subroutine angular_wfnpz_pw(wavefn,nkpt,kxyz,tau,rc,cart,  &
                          kkx,kky,kkz,vcell,myproc,ngrid,step,Rl,cnk,Rlm,proj_flag)


        use all_to_all_module
        implicit none

      include 'mpif.h'


	integer ii,jj,kk,ll,ngrid,ierr,myproc,proj_flag
        integer nkpt,kxyz(3,nkpt)
        double precision tau(3),rc,kkx,kky,kkz,cart(3,3),vcell,Rl(ngrid)
	complex (kind=8) wavefn(nkpt)
        double complex::Rlm(ngrid,3)

	double precision k1,k2,k3,kgr,kgr0,k10,k20,k30,ekin0

        double complex::cnk(3)

	double precision bess(4),r1,r2,r3,r4,r5,r6,r7,r8
        double precision ekin,const1,const0
        double precision SBESSJ,kgdotr0,scale,step,norm,step2
	double complex::  phase,ylm(3),ione,ctem1,ctem2

!        double complex, allocatable::Rlm(:,:),Rlm_sum(:,:)
        double complex, allocatable::Rlm_sum(:,:)
        double complex, allocatable::Rlm2(:,:)

        ione=(0.d0,1.d0)
      
        const1=dsqrt(3.d0/(2.d0*pi4))
        const0=dsqrt(3.d0/(pi4))


!============================================================================
!
         
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
!     if(myproc.eq.0) then
!     do jj=1,3
!         write(233,*) jj,REAL(cnk(jj)),REAL(sqrt(cnk(jj)))
!      end do
!      end if

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
        subroutine angular_wfnpz_sw(wavefn,nkpt,kxyz,tau,rc,cart,  &
                          kkx,kky,kkz,vcell,myproc,ngrid,step,Rl,cnk,Rlm,proj_flag)


        use all_to_all_module
        implicit none

      include 'mpif.h'


	integer ii,jj,kk,ll,ngrid,ierr,myproc,proj_flag
        integer nkpt,kxyz(3,nkpt)
        double precision tau(3),rc,kkx,kky,kkz,cart(3,3),vcell,Rl(ngrid)
	complex (kind=8) wavefn(nkpt)
        double complex::Rlm(ngrid,1)

	double precision k1,k2,k3,kgr,kgr0,k10,k20,k30,ekin0

        double complex::cnk(1)

	double precision bess(4),r1,r2,r3,r4,r5,r6,r7,r8
        double precision ekin,const1,const0
        double precision SBESSJ,kgdotr0,scale,step,norm,step2
	double complex::  phase,ylm(1),ione,ctem1,ctem2

!        double complex, allocatable::Rlm(:,:),Rlm_sum(:,:)
        double complex, allocatable::Rlm_sum(:,:)
        double complex, allocatable::Rlm2(:,:)

        ione=(0.d0,1.d0)
      
       const0=1.d0/dsqrt(pi4)



!============================================================================
!
         
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

              bess(1)=SBESSJ(0,kgr)

              ctem2=ctem1*bess(1)

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
