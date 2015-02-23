!-------------------------------------------------------------------------

        subroutine lsdau_proj_init(k_gspace,kpoints,crys)

        use lsdau_shared

        use constants
        use crystal_module
        use parallel_gspace_module
        use kpoint_module

        type(kpoint), intent(in) :: kpoints  
        type(parallel_gspace), intent(in) :: k_gspace(*)
        type(crystal) :: crys       ! crystal structure

!-------------------------------------------------------------------------

        integer::it,i2,i1,i0,G1,G2,G3
        integer::len_gk_max,ngrid_max
	integer ii,jj,irk

        double precision::k_1,k_2,k_3,k_x,k_y,k_z,kgr,kg_x,kg_y,kg_z,ekin,ekin0

        double precision const21,const22,const23,const11,const12,const0
        double precision SBESSJ
	double complex:: ione,itwo

!-------------------------------------------------------------
       i2=0
       i1=0
       i0=0
       len_gk_max=MAXVAL(k_gspace(1:kpoints%nrk)%length)
       ngrid_max=MAXVAL(ngrid1)

       do it=1,crys%ntype
          if(lproj_tem(it).eq.2) i2=1+1
          if(lproj_tem(it).eq.1) i1=1+1
          if(lproj_tem(it).eq.0) i0=1+1
       end do

       if(i2.ge.1) then
         allocate(ylm2(5,len_gk_max,kpoints%nrk))

! this matrix could be huge
         allocate(sb2(ngrid_max,len_gk_max,kpoints%nrk))
       end if
 
       if(i1.ge.1) then
         allocate(ylm1(3,len_gk_max,kpoints%nrk))
         allocate(sb1(ngrid_max,len_gk_max,kpoints%nrk))
       end if

       if(i0.ge.1) then
         allocate(sb0(ngrid_max,len_gk_max,kpoints%nrk))
       end if

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

! Rlm(r)=4\pi i^l e(ikr_0)\sum_G C^k(G) e(iGr_0) j_l(|k+G|r)Y_lm(\Omega_{k+G})
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

           kg_x = crys%bvec(1,1)*G1+crys%bvec(1,2)*G2+  &
                crys%bvec(1,3)*G3+k_x
           kg_y = crys%bvec(2,1)*G1+crys%bvec(2,2)*G2+  &
                crys%bvec(2,3)*G3+k_y
           kg_z = crys%bvec(3,1)*G1+crys%bvec(3,2)*G2+  &
                crys%bvec(3,3)*G3+k_z

           ekin = kg_x*kg_x+kg_y*kg_y+kg_z*kg_z

!--------------------------------------------------------------------------------------- 
!        ylm1(3,len_gk_max,kpoints%nrk)

           if(ekin.lt.(1.d-10)) then

             if(i1.ge.1) then
              do im=1,twolp1
                 ylm1(im,jj,irk)=(0.d0,0.d0)
              end do
             end if
             if(i1.ge.2) then
              do im=1,twolp1
                 ylm2(im,jj,irk)=(0.d0,0.d0)
              end do
             end if

           else
             
              if(i2.ge.1) then
              ylm2(1,jj,irk)= const21*(kg_x*kg_x-kg_y*kg_y+itwo*kg_x*kg_y)/ekin       !Y2(-2)^*
              ylm2(2,jj,irk)= const22*(kg_x*kg_z+ione*kg_y*kg_z)/ekin             !Y2(-1)^*
              ylm2(3,jj,irk)= const23*(3.d0*kg_z*kg_z/ekin-1.d0)              !Y2(0)^*
              ylm2(4,jj,irk)= -DCONJG(ylm2(2,jj,irk))                            !Y2(1)^*
              ylm2(5,jj,irk)= DCONJG(ylm2(1,jj,irk))                             !Y2(2)^*
              end if

              if(i1.ge.1) then
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

             if(i2.ge.1) sb2(ii,jj,irk)=SBESSJ(2,kgr)
             if(i1.ge.1) sb1(ii,jj,irk)=SBESSJ(1,kgr)
             if(i0.ge.1) sb0(ii,jj,irk)=SBESSJ(0,kgr)
          end do
        end do
        end do

      return
      end
