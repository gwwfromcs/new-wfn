!-*-Fortran-*-
!
subroutine setup_nonlocal_potential(ham, pspot, crys,iter)

  include 'use.h'   
  implicit none                    ! implicit? Just say no!
  include 'interface.h'
  include 'flibcalls.ph'
  integer iter  
  type(hamiltonian), intent(inout) :: ham  
  type(pseudo_potential), intent(inout), target :: pspot    
  type(crystal), intent(in) :: crys  
  !
  !     Set up the nonlocal part of the Hamiltonian for a given k-point.
  !     Kleinman-Bylander Pseudopotential.
  !
  !     all it does is setting up the arrays ham%vnloc and ham%xnorm
  !     at the kpoint specified by ham%gspace%rk
  !
  !     2001 added real space Non-Loc stuff and if statement for rspace option
  !
  !     1996 Bernd Pfrommer (parallel version)
  !          I changed the JLM representation "anlga" to be the conjugate
  !          complex "vnloc", because it allows a more efficient use of
  !          BLAS calls later in the diagonalization routine.
  !
  !     1990 JLM
  !     adapted from Sverre Froyen plane wave program
  !     --------------------- local variables ----------------------------
  !
  !     ----- variables for the gspace loop ---------
  !
  integer :: igv(4), fftn(4), ffth(4), igv3, irod, iord, igs, nt,ico1,ia,nr  
  real(dp) :: gv(3)  
  integer :: ind, i, k, kk, l, ni, ni1, lmmin, lmmax, lm, j, jcar,neg,&
             jabs  
  real(dp) :: qi, qk(3), qcar(3), flm(16), qinv,qinv1,qinv2,qinv3, fi, xni, &
       vq, xi, xdum, xa, xb,rkmod(3,4),ek,ek1,ek2
  real(dp), parameter :: eps = 1.0d-8
  complex(dp) :: sum  
  complex(dp), allocatable :: st(:)     ! work array for structure factors
  !
  !     ------------------------------------------------------------------
  !
  if(iter .eq. 1 .or. .not. pspot%NL_rspace(1)) then

  ind = 0  
  do k = 1, crys%ntype  
     do l = 1, 5  
        if (pspot%nkb(l, k) /= 0) then  
           lmmin = pspot%lo(l, k) * pspot%lo(l, k) + 1  
           lmmax = (pspot%lo(l, k) + 1) * (pspot%lo(l, k) + 1)  
           do lm = lmmin, lmmax  
              do kk = 1, crys%natom(k)  
                 ind = ind + 1  
                 ham%xnorm(ind) = real(pspot%nkb(l, k), dp)  
               end do
           end do
        end if
     end do
  end do 

  allocate(st(crys%mxdatm))  
  fftn(1:3) = ham%gspace%fftsize(1:3)  
  fftn(4) = ham%gspace%fftsize(3)  
  ffth(:) = fftn(:) / 2  
  !
  !      starts loop over g-vectors in small gspace
  !
  igs = 0  
  do iord = 1, ham%gspace%lorder                  ! loop through x/y gspace
     irod = ham%gspace%order(1, iord)  
     igv(1) = irod / ham%gspace%fftsize(2)  
     igv(2) = mod(irod, ham%gspace%fftsize(2))  
     igv(3) = ham%gspace%order(2, iord)  
     igv(4) = ham%gspace%order(3, iord)  
     igv(:) = mod(igv(:) + ffth(:), fftn(:)) - ffth(:)  
     gv(1:2) = real(igv(1:2), dp)  
     do igv3 = igv(3), igv(4)                     ! loop over z axis
        gv(3) = real(igv3, dp)  
        igs = igs + 1  
        qi = sqrt(ham%gspace%ekin(igs))           ! get kinetic energy from
        qk(:) = ham%gspace%rk(:) + gv(:)  
        !
        !           cartesian of k+G
        !
        qcar(1) = crys%bvec(1, 1) * qk(1) + crys%bvec(1, 2) * qk(2) + &
             crys%bvec(1, 3) * qk(3)
        qcar(2) = crys%bvec(2, 1) * qk(1) + crys%bvec(2, 2) * qk(2) + &
             crys%bvec(2, 3) * qk(3)
        qcar(3) = crys%bvec(3, 1) * qk(1) + crys%bvec(3, 2) * qk(2) + &
             crys%bvec(3, 3) * qk(3)
        !
        !           ----- compute angular functions
        !
        flm (1) = done  
        if (qi > eps) then  
           qinv1 = done / qi  
           flm(2) = qcar(1) * qinv1  
           flm(3) = qcar(2) * qinv1  
           flm(4) = qcar(3) * qinv1  
           qinv2 = qinv1 * qinv1 
           flm(5) = drt3 * qcar(1) * qcar(2) * qinv2 
           flm(6) = drt3 * qcar(2) * qcar(3) * qinv2 
           flm(7) = drt3 * qcar(3) * qcar(1) * qinv2 
           flm(8) = dtrhalf * qcar(3) * qcar(3) * qinv2- dhalf
           flm(9) = drt3 * (qcar(1) * qcar(1) - qcar(2) * qcar(2)) * qinv2 * &
                dhalf

           qinv3 = qinv1 * qinv2
           flm(10) = qinv3*0.5d0*drt5*oort2* (3.d0*qcar(1)*qcar(1)*qcar(2) - qcar(2)* qcar(2)*qcar(2) )
           flm(11) = qinv3*drt15*qcar(1)*qcar(2)*qcar(3)
           flm(12) = 0.5d0*drt3*oort2*(5.d0*qcar(3)*qcar(3)*qinv2-1.d0)*qcar(2)*qinv1
           flm(13) = 0.5d0*(5.d0*qcar(3)*qcar(3)*qcar(3)*qinv3-3.d0*qcar(3)*qinv1)
           flm(14) = 0.5d0*drt3*oort2*(5.d0*qcar(3)*qcar(3)*qinv2-1.d0)*qcar(1)*qinv1
           flm(15) = qinv3*0.5d0*drt15*(qcar(1)*qcar(1)-qcar(2)*qcar(2))*qcar(3)
           flm(16) =-qinv3*0.5d0*drt5*oort2* (3.d0*qcar(2)*qcar(2)*qcar(1) - qcar(1)* qcar(1)*qcar(1) )

        else  
           do lm = 2, 16  
              flm(lm) = dzero  
           end do
        end if
        !
        !        starts loop over second index
        !
        ind = 0  
        do k = 1, crys%ntype  
           do kk = 1, crys%natom(k)            ! compute complex phase factor
              fi = gv(1) * crys%rat(1, kk, k) + gv(2) * crys%rat(2, kk, k) + &
                   gv(3) * crys%rat(3, kk, k)
              st(kk) = exp(cmplx(dzero, fi, dp))  
           end do
           !
           ! loop over lo is replaced by loop over l
           ! and l quantum number is given  explicitely by referencing psp
           !
           do l = 1, 5                         ! loop over angular momenta l
              if (pspot%nkb(l, k) /= 0) then   ! found pot of that ang. mom.
                 ! interpolate potential
                 xni = qi / pspot%delqnl(k) + dtwo  
                 !
                 ! this was the old extrapolation. changed 1996 pfrommer/ma
                 !
                 !      ni  = xni + half
                 !      if(ni .le. 3) ni = 4
                 !      vq = zero
                 !      if(ni .lt. pspot%nqnl(k)) then
                 !      xi  = xni - dble(ni)
                 !      vq = pspot%vkb(ni,l,k) * (one+xi) * (one-xi)
                 !     $  + half * (pspot%vkb(ni+1,l,k)*(one+xi)
                 !     $  -pspot%vkb(ni-1,l,k)*(one-xi)) * xi
                 !      endif
                 !
                 ! cubic spline interpolation
                 !
                 vq = dzero  
                 ni = xni  
                 if (ni <= 2) ni = 3  
                 ni1 = ni + 1  
                 if (ni < pspot%nqnl(k)) then  
                    xa = real(ni1, dp) - xni  
                    xb = xni - real(ni, dp)  
                    vq = xa * pspot%vkb(ni, l, k) + xb * &
                         pspot%vkb(ni1, l, k) + ((xa**3 - xa) * &
                         pspot%d2vkbdq2(ni - 2, l, k) + (xb**3 - xb) * &
                         pspot%d2vkbdq2(ni1 - 2, l, k)) * dsixth
                 end if
                 lmmin = pspot%lo(l, k) * pspot%lo(l, k) + 1  
                 lmmax = (pspot%lo(l, k) + 1) * (pspot%lo(l, k) + 1)
                 !
                 do lm = lmmin, lmmax     ! loop over m quantum number
                    do kk = 1, crys%natom(k)  
                       ind = ind + 1  
                       !
                       !  conjugate here is for structure factor.
                       !  we use a -G*r convention that inverts the axes
                       !
                       ham%vnloc%data(igs, ind, 1) =conjg(st(kk) * &
                            flm(lm) * vq)
                    end do
                 end do
              end if
           end do
        end do           ! end of loop over atomic types
     end do              ! end of loop over 3rd dimension
  end do                 ! end of loop over small gspace
  !      do k=1,crys%ntype
  !         do l=1,5
  !            if(pspot%nkb(l,k) .ne. 0) then ! found pot of that ang. mom
  !               sum =0
  !               do i=1, size(pspot%vkb(:,l,k))
  !                  sum = sum + pspot%vkb(i,l,k)*pspot%vkb(i,l,k)
  !     $            *sqrt(crys%vcell)
  !     $            *pspot%delqnl(k)*(2*4*atan(one))**(-3)*4*4*atan(one)
  !               enddo
  !            endif
  !         end do
  !      end do
  !      check for norm
  !
  !     do ind=1, pspot%nanl
  !        sum =0
  !        do i=1, ham%gspace%length
  !           sum = sum +
  !    $           ham%vnloc%data(i,ind,1)*
  !    $           conjg(ham%vnloc%data(i,ind,1))
  !        end do
  !        write(9,*) 'second a2norm/anorm for projector',ind,'is:', sum
  !     end do
  !
  !      normalization matrix
  !
  deallocate(st) 
  end if
  !
  !       write(9,200)
  !       do ind=1,pspot%nanl
  !          write(9,300) (ham%vnloc%data(k,ind,1), k=1,ham%gspace%length)
  !       end do
  if (pspot%NL_rspace(1)) then

    ind = 0  
    do k = 1, crys%ntype  
      do kk = 1, crys%natom(k)
        do l = 1, 5  
          if (pspot%nkb(l, k) /= 0) then  
            lmmin = pspot%lo(l, k) * pspot%lo(l, k) + 1  
            lmmax = (pspot%lo(l, k) + 1) * (pspot%lo(l, k) + 1)  
            do lm = lmmin, lmmax  
              ind = ind + 1  
              ham%rsp_norm(ind) = real(pspot%nkb(l, k), dp)  
            end do
          end if
        end do
      end do
    end do 

    nr=ham%fftstruc%fftsize(1)*ham%fftstruc%fftsize(2)*&
         ham%fftstruc%fftsize(3)
    call mdscal(pspot%nanl,crys%vcell/nr,ham%rsp_norm(1),1)

    qk(:) = ham%gspace%rk(:)

    qcar(1) = crys%bvec(1, 1) * qk(1) + crys%bvec(1, 2) * qk(2) + &
             crys%bvec(1, 3) * qk(3)
    qcar(2) = crys%bvec(2, 1) * qk(1) + crys%bvec(2, 2) * qk(2) + &
             crys%bvec(2, 3) * qk(3)
    qcar(3) = crys%bvec(3, 1) * qk(1) + crys%bvec(3, 2) * qk(2) + &
             crys%bvec(3, 3) * qk(3)    

    ico1=0
    do ia=1,pspot%natom_tot
      do i=1,pspot%nmap(ia)
        ico1=ico1+1
        ham%cphase(ico1)=exp(-zi*(pspot%xyzmap(ico1*3-2)*qcar(1)+ &
          pspot%xyzmap(ico1*3-1)*qcar(2)+pspot%xyzmap(ico1*3)*qcar(3))) 
      enddo
    enddo


  end if 

  return  

200 format (/' Vectors: ')  

300 format (2000('(',f10.4,',',f10.4,')'))  

end subroutine setup_nonlocal_potential
