!     @process extchk
!
subroutine forstressnloc(ipr, crys, bands, kgs, ekinmod, pspot, &
     rk, aadot, eigenvec, fsum, ssum)
  !
  !     computes the kinetic energy and nonlocal potential contributions
  !     to the Hellman-Feynman forces and the stress.
  !     Adapted from a code by J.L. Martins 1988
  !
  !     parallel version 1995/1996 Bernd Pfrommer
  !
  use all_to_all_module  
  include 'use.h'  
  implicit none                   ! never do implicit !
  include 'interface.h'  
  include 'all_to_all.h'  
  include 'flibcalls.ph'  
  !
  !     INPUT:
  !     -----
  !
  type(parallel_gspace), intent(in) :: &
       kgs(*)                            ! gspaces for the different k-points
  type(pseudo_potential), intent(in) :: pspot  
  type(band), intent(in) :: bands  
  type(crystal), intent(in) :: crys  
  type(complex_gspace_array), intent(in) :: &
       eigenvec                               ! the groundstate wavefunctions
  integer, intent(in) :: ipr                                     ! print flag
  real(dp), intent(in) :: &
       rk(3, bands%nrk), &          ! k-points in units of reciprocal lattice
       ekinmod(3), &             ! parameters for the modified kinetic energy
       aadot(6)                           ! metric in realspace:
                                          ! adot(11,22,33,12,23,31)/(2*pi)**2
  !
  !     INPUT/OUTPUT:  (the nonlocal+kinetic part is added)
  !     -------------
  !
  real(dp), intent(inout) :: &
       fsum(3, crys%mxdatm, crys%ntype), & ! unsymm forces [latt. coordinates]
                                                ! from local part on each atom
       ssum(6)                     ! unsymmetrized, metric stress tensor [Ryd]
  !
  !     ---------------- local variables ---------------------------------
  !     dynamic arrays:
  !     --------------
  complex(dp), allocatable :: anlga(:,:,:), dhd(:,:,:)  
  real(dp), allocatable :: xnlkb(:)
  !  
  !     other simple variables:
  !     -----------------------
  integer :: i, j, k, l, n, kk, kl, lm, ind, lo, ist, is2, is3, irk, &
       neig, ndim, jmin, jmax, mxdpdim, nanl, nanl10, nanlga, nanlga10, &
       ni, lmmin, lmmax, is
  real(dp) :: &
       delta, rkpt(3), qki(3), qk(3), qi, qinv, dqi(6), qcar(3), &
       dqcar(6, 3), flm(16), dflm(6, 16), forefac, zsumr, fi, xni, &
       vq, dvq, xi, enl, temp(6)
  logical :: isp  
  complex(dp) :: st
  real(dp), parameter :: eps = 1.0d-12 

  isp = .false.  
  if (crys%nspin == 2) isp = .true.  

  forefac = dtwo / (sqrt(pi) * ekinmod(3)) * ekinmod(1)  
  !     ------- first allocate dynamic arrays ------------
  mxdpdim = 0  
  do i = 1, bands%nrk  
     if (kgs(i)%length > mxdpdim) mxdpdim = kgs(i)%length  
  end do
  !
  !     count the max number of nonlocal components nanlga for array dimen
  !
  nanlga = 0  
  do k = 1, crys%ntype  
     ind = 0  
     do l = 1, 5  
        if (pspot%nkb(l, k) /= 0) then  
           lmmin = pspot%lo(l, k) * pspot%lo(l, k) + 1  
           lmmax = (pspot%lo(l, k) + 1) *(pspot%lo(l, k) + 1)  
           do lm = lmmin, lmmax  
              ind = ind + 1  
           end do
        end if
     end do
     nanlga = max(nanlga, ind)
  end do

  nanlga10 = 10 * nanlga
  allocate(anlga(10, nanlga, mxdpdim)) ; anlga = zzero
  allocate(dhd(10, nanlga, bands%max)) ; dhd = zzero

  if (pspot%nanl > 0) allocate(xnlkb(pspot%nanl))  
  do is = 1, crys%nspin  
     do irk = 1, bands%nrk                           ! loop over all k-points
        neig = bands%nband(irk, is)  
        ndim = kgs(irk)%length  
        rkpt(1) = rk(1, irk)  
        rkpt(2) = rk(2, irk)  
        rkpt(3) = rk(3, irk)  
        !        ----------- kinetic energy contribution to the stress ---------
        !        find min and max occupied band
        jmin = 0  
        jmax = 0  
        do i = 1, neig  
           if (bands%occup(i, irk, is) /= dzero .and. jmin == 0) jmin = i  
           if (bands%occup(i, irk, is) /= dzero .and. jmin /= 0) jmax = i  
        end do
        !
        !      start loop over g(i)
        !
        do i = 1, ndim  
           !
           !           compute |g+k|
           !
           qki(1) = rkpt(1) + real(kgs(irk)%gvec(1, i), dp)  
           qki(2) = rkpt(2) + real(kgs(irk)%gvec(2, i), dp)  
           qki(3) = rkpt(3) + real(kgs(irk)%gvec(3, i), dp)  
           !
           !           do eigenvector multiplication and band sum
           !
           zsumr = dzero  
           do k = jmin, jmax  
              zsumr = zsumr + bands%occup(k, irk, is) / &
                   real(crys%nspin, dp) * &
                   (real(eigenvec%data((k - 1) * ndim + i, irk, is)) * &
                   real(eigenvec%data((k - 1) * ndim + i, irk, is)) + &
                   aimag(eigenvec%data((k - 1) * ndim + i, irk, is)) * &
                   aimag(eigenvec%data((k - 1) * ndim + i, irk, is)))
           end do
           !
           !           add the kinetic energy functional modifications
           !
           delta = (kgs(irk)%ekin(i) - ekinmod(2)) / ekinmod(3)  

           zsumr = zsumr * (done + forefac * exp(-delta * delta))  

           zsumr = dtwo * zsumr  
           ssum(1) = ssum(1) + zsumr * qki(1) * qki(1)  
           ssum(2) = ssum(2) + zsumr * qki(2) * qki(2)  
           ssum(3) = ssum(3) + zsumr * qki(3) * qki(3)  
           ssum(4) = ssum(4) + zsumr * qki(1) * qki(2)  
           ssum(5) = ssum(5) + zsumr * qki(2) * qki(3)  
           ssum(6) = ssum(6) + zsumr * qki(3) * qki(1)  
        end do                                     ! end of loop over k-points
     end do                                           ! end of loop over spins
  end do
  !     ----------- nonlocal potential contribution to forces and stress -
  !
  !     Consolidate across processors by summing. THIS CAN ONLY BE DONE ON
  !
  call all_sum_all(fsum, 3 * crys%mxdatm * crys%ntype)  
  call all_sum_all(ssum, 6)  

  temp=dzero

  do irk = 1, bands%nrk                               ! loop over all k-points
     ndim = kgs(irk)%length  
     rkpt(1) = rk(1, irk)  
     rkpt(2) = rk(2, irk)  
     rkpt(3) = rk(3, irk)  
     do k = 1, crys%ntype  
        !
        nanl = 0  
        do l = 1, 5  
           if (pspot%nkb(l, k) /= 0) then  
              lmmin = pspot%lo(l, k) * pspot%lo(l, k) + 1  
              lmmax = (pspot%lo(l, k) + 1) * (pspot%lo(l, k) + 1)  
              do lm = lmmin, lmmax  
                 nanl = nanl + 1  
              end do
           end if
        end do
        nanl10 = 10 * nanl  
        !
        do kk = 1, crys%natom(k)  
           !
           !                calculates the matrix elements
           !
           do i = 1, ndim  
              qk(1) = rkpt(1) + real(kgs(irk)%gvec(1, i), dp)  
              qk(2) = rkpt(2) + real(kgs(irk)%gvec(2, i), dp)  
              qk(3) = rkpt(3) + real(kgs(irk)%gvec(3, i), dp)  
              qi = (qk(1) * crys%bdot(1, 1) + qk(2) * crys%bdot(2, 1) + &
                   qk(3) * crys%bdot(3, 1)) * qk(1) + &
                   (qk(1) * crys%bdot(1, 2) + qk(2) * crys%bdot(2, 2) + &
                   qk(3) * crys%bdot(3, 2)) * qk(2) + &
                   (qk(1) * crys%bdot(1, 3) + qk(2) * crys%bdot(2, 3) + &
                   qk(3) * crys%bdot(3, 3)) * qk(3)
              if (qi > eps) then  
                 qi = sqrt(qi)  
                 qinv = done / qi  
                 !
                 !         derivatives of qi with respect to strain
                 !         in lattice coordinates
                 !
                 dqi(1) = -qk(1) * qk(1) * qinv  
                 dqi(2) = -qk(2) * qk(2) * qinv  
                 dqi(3) = -qk(3) * qk(3) * qinv  
                 dqi(4) = -qk(1) * qk(2) * qinv  
                 dqi(5) = -qk(2) * qk(3) * qinv  
                 dqi(6) = -qk(3) * qk(1) * qinv  
                 !
                 !         qvector in cartesian coordinates and normalized
                 !
                 qcar(1) = crys%bvec(1, 1) * qk(1) + &
                      crys%bvec(1, 2) * qk(2) + crys%bvec(1, 3) * qk(3)
                 qcar(2) = crys%bvec(2, 1) * qk(1) + &
                      crys%bvec(2, 2) * qk(2) + crys%bvec(2, 3) * qk(3)
                 qcar(3) = crys%bvec(3, 1) * qk(1) + &
                      crys%bvec(3, 2) * qk(2) + crys%bvec(3, 3) * qk(3)
                 qcar(1) = qcar(1) * qinv  
                 qcar(2) = qcar(2) * qinv  
                 qcar(3) = qcar(3) * qinv  
                 !
                 !         derivatives of qcar with respect to strain
                 !         in lattice coordinates
                 !
                 dqcar(1, 1) = -crys%avec(1, 1) * qk(1) / pi2  
                 dqcar(2, 1) = -crys%avec(1, 2) * qk(2) / pi2  
                 dqcar(3, 1) = -crys%avec(1, 3) * qk(3) / pi2  
                 dqcar(4, 1) = -crys%avec(1, 1) * qk(2) / pi2  
                 dqcar(5, 1) = -crys%avec(1, 2) * qk(3) / pi2  
                 dqcar(6, 1) = -crys%avec(1, 3) * qk(1) / pi2  
                 dqcar(1, 2) = -crys%avec(2, 1) * qk(1) / pi2  
                 dqcar(2, 2) = -crys%avec(2, 2) * qk(2) / pi2  
                 dqcar(3, 2) = -crys%avec(2, 3) * qk(3) / pi2  
                 dqcar(4, 2) = -crys%avec(2, 1) * qk(2) / pi2  
                 dqcar(5, 2) = -crys%avec(2, 2) * qk(3) / pi2  
                 dqcar(6, 2) = -crys%avec(2, 3) * qk(1) / pi2  
                 dqcar(1, 3) = -crys%avec(3, 1) * qk(1) / pi2  
                 dqcar(2, 3) = -crys%avec(3, 2) * qk(2) / pi2  
                 dqcar(3, 3) = -crys%avec(3, 3) * qk(3) / pi2  
                 dqcar(4, 3) = -crys%avec(3, 1) * qk(2) / pi2  
                 dqcar(5, 3) = -crys%avec(3, 2) * qk(3) / pi2  
                 dqcar(6, 3) = -crys%avec(3, 3) * qk(1) / pi2  
                 !
                 do kl = 1, 6
                    dqcar(kl, 1) = (dqcar(kl, 1) - dqi(kl) * qcar(1)) * qinv
                    dqcar(kl, 2) = (dqcar(kl, 2) - dqi(kl) * qcar(2)) * qinv
                    dqcar(kl, 3) = (dqcar(kl, 3) - dqi(kl) * qcar(3)) * qinv
                 end do
                 !
                 !         angular functions
                 !
                 flm(1) = done  
                 flm(2) = qcar(1)  
                 flm(3) = qcar(2)  
                 flm(4) = qcar(3)  
                 flm(5) = drt3 * qcar(1) * qcar(2)  
                 flm(6) = drt3 * qcar(2) * qcar(3)  
                 flm(7) = drt3 * qcar(3) * qcar(1)  
                 flm(8) = dthree * qcar(3) * qcar(3) * dhalf - dhalf
                 flm(9) = drt3 * (qcar(1) * qcar(1) - qcar(2) * qcar(2)) * dhalf
                 flm(10) = 0.5d0*drt5*oort2* (3.d0*qcar(1)*qcar(1)*qcar(2) - qcar(2)* qcar(2)*qcar(2) )
                 flm(11) = drt15*qcar(1)*qcar(2)*qcar(3)
                 flm(12) = 0.5d0*drt3*oort2*(5.d0*qcar(3)*qcar(3)-1.d0)*qcar(2)
                 flm(13) = 0.5d0*(5.d0*qcar(3)*qcar(3)*qcar(3)-3.d0*qcar(3))
                 flm(14) = 0.5d0*drt3*oort2*(5.d0*qcar(3)*qcar(3)-1.d0)*qcar(1)
                 flm(15) = 0.5d0*drt15*(qcar(1)*qcar(1)-qcar(2)*qcar(2))*qcar(3)
                 flm(16) =-0.5d0*drt5*oort2* (3.d0*qcar(2)*qcar(2)*qcar(1) - qcar(1)* qcar(1)*qcar(1) )

                 !
                 !                    strain derivatives of angular functions
                 !                    in lattice coordinates
                 !
                 dflm=dzero
                 do kl = 1, 6  
                    dflm(kl, 1) = dzero  
                    dflm(kl, 2) = dqcar(kl, 1)  
                    dflm(kl, 3) = dqcar(kl, 2)  
                    dflm(kl, 4) = dqcar(kl, 3)  
                    dflm(kl, 5) = drt3 * (dqcar(kl, 1) * qcar(2) + &
                         qcar(1) * dqcar(kl, 2))
                    dflm(kl, 6) = drt3 * (dqcar(kl, 2) * qcar(3) + &
                         qcar(2) * dqcar(kl, 3))
                    dflm(kl, 7) = drt3 * (dqcar(kl, 3) * qcar(1) + &
                         qcar(3) * dqcar(kl, 1))
                    dflm(kl, 8) = dthree * dqcar(kl, 3) * qcar(3)  
                    dflm(kl, 9) = drt3 * (dqcar(kl, 1) * qcar(1) - &
                         dqcar(kl, 2) * qcar(2))


                 dflm(kl,10) = 0.5d0*drt5*oort2*  &
                              (6.d0*dqcar(kl,1)*qcar(1)*qcar(2)+3.d0*qcar(1)*qcar(1)*dqcar(kl,2) &
                               - 3.d0*qcar(2)* qcar(2)*dqcar(kl,2) )

                 dflm(kl,11) = drt15*(dqcar(kl,1)*qcar(2)*qcar(3)+ &
                                      dqcar(kl,2)*qcar(1)*qcar(3)+ &
                                      dqcar(kl,3)*qcar(1)*qcar(2))

                 dflm(kl,12) = 0.5d0*drt3*oort2* &
                               (10.d0*dqcar(kl,3)*qcar(3)*qcar(2)+ &
                                (5.d0*qcar(3)*qcar(3)-1.d0)*dqcar(kl,2))

                 dflm(kl,13) = 0.5d0*(15.d0*qcar(3)*qcar(3)-3.d0)*dqcar(kl,3)

                 dflm(kl,14) = 0.5d0*drt3*oort2* &
                              ( (5.d0*qcar(3)*qcar(3)-1.d0)*dqcar(kl,1)+ &
                                10.d0*qcar(3)*dqcar(kl,3)*qcar(1) )

                 dflm(kl,15) = 0.5d0*drt15* &
                               ( (2.d0*dqcar(kl,1)*qcar(1)-2.d0*dqcar(kl,2)*qcar(2))*qcar(3) + &
                                 (qcar(1)*qcar(1)-qcar(2)*qcar(2))*dqcar(kl,3) )

                 dflm(kl,16) =-0.5d0*drt5*oort2* &
                              ( 6.d0*dqcar(kl,2)*qcar(2)*qcar(1) + 3.d0*qcar(2)*qcar(2)*dqcar(kl,1) - &
                                3.d0*dqcar(kl,1)*qcar(1)*qcar(1) )

                 end do
              else                                  ! treat the G=0 component
                 qi = dzero  
                 dqi(1) = dzero  
                 dqi(2) = dzero  
                 dqi(3) = dzero  
                 dqi(4) = dzero  
                 dqi(5) = dzero  
                 dqi(6) = dzero  
                 do lm = 1, 16 
                    flm(lm) = dzero  
                    dflm(1, lm) = dzero  
                    dflm(2, lm) = dzero  
                    dflm(3, lm) = dzero  
                    dflm(4, lm) = dzero  
                    dflm(5, lm) = dzero  
                    dflm(6, lm) = dzero  
                 end do
                 flm(1) = done  
              end if
              !
              !                 structure factor
              !
              fi = real(kgs(irk)%gvec(1, i), dp) * crys%rat(1, kk, k) + &
                   real(kgs(irk)%gvec(2, i), dp) * crys%rat(2, kk, k) + &
                   real(kgs(irk)%gvec(3, i), dp) * crys%rat(3, kk, k)
              st = exp(cmplx(dzero, fi, dp))  
              !
              ind = 0  
              do l = 1, 5  
                 if (pspot%nkb(l, k) /= 0) then  
                    xni = qi / pspot%delqnl(k) + dtwo  
                    ni = int(xni + dhalf)  
                    if (ni <= 3) ni = 4  
                    vq = dzero  
                    dvq = dzero  
                    if (ni < pspot%nqnl(k)) then  
                       xi = xni - real(ni, dp)  
                       vq = pspot%vkb(ni, l, k) * (done + xi) * (done - xi) + &
                            dhalf * (pspot%vkb(ni + 1, l, k) * (done + xi) - &
                            pspot%vkb(ni - 1, l, k) * (done - xi)) * xi
                       dvq = dmtwo * xi * pspot%vkb(ni, l, k) + &
                            (xi + dhalf) * pspot%vkb(ni + 1, l, k) + &
                            (xi - dhalf) * pspot%vkb(ni - 1, l, k)
                       dvq = dvq / pspot%delqnl(k)  
                    end if
                    lmmin = pspot%lo(l, k) * pspot%lo(l, k) + 1  

                    lmmax = (pspot%lo(l, k) + 1) * (pspot%lo(l, k) + 1)
                    !
                    do lm = lmmin, lmmax  
                       ind = ind + 1  
                       anlga(1, ind, i) = st * flm(lm) * vq  
                       xnlkb(ind) = real(pspot%nkb(l, k), dp)  
                       !
                       !                          forces
                       !
                       anlga(2, ind, i) = cmplx(dzero, dmone, dp) * &
                            st * flm(lm) * vq * real(kgs(irk)%gvec(1, i), dp)
                       anlga(3, ind, i) = cmplx(dzero, dmone, dp) * &
                            st * flm(lm) * vq * real(kgs(irk)%gvec(2, i), dp)
                       anlga(4, ind, i) = cmplx(dzero, dmone, dp) * &
                            st * flm(lm) * vq * real(kgs(irk)%gvec(3, i), dp)
                       !
                       !                          stresses
                       !
                       do ist = 1, 6
                          anlga(4 + ist, ind, i) = st * ( dflm(ist, lm)  * &
                               vq   + &
                     flm(lm) *  dvq  * dqi(ist) )
                       end do
                    end do
                 end if
                 !
              end do                     ! loop over lo
              !
           end do                        ! loop over gspace
           !
           !     calculates the matrix vector product
           !
           do is = 1, crys%nspin  
              neig = bands%nband(irk, is)  
              if (nanl10 > 0) call mzgemm('N', 'N', nanl10, neig, ndim, &
                   zone, anlga(1, 1, 1) ,nanlga10, eigenvec%data(1, irk, is), &
                   ndim, zzero, dhd(1, 1, 1), nanlga10)
              !
              !              cross-processor consolidation for dhd
              !
              call all_sum_all(dhd, nanlga10 * bands%max)  


!            write(9,*) dhd(5,1,1),dhd(5,2,1),dhd(5,3,1)

              do n = 1, neig  

!       do j=1,kgs(irk)%length
!          write(9,*) eigenvec%data((n-1)*ndim+j,irk,is),j
!        end do
                 do i = 1, nanl  
                    enl =bands%occup(n, irk, is) / real(crys%nspin, dp) * &
                         xnlkb(i) * (real(dhd(1, i, n), dp) * &
                         real(dhd(1, i, n), dp) + aimag(dhd(1, i, n)) * &
                         aimag(dhd(1, i, n)))
!              write(9,*) enl,dhd(1,i,n),n,ssum(1)
                    do ist = 1, 3  
                       fsum(ist, kk, k) = fsum(ist, kk, k) + dtwo * &
                            bands%occup(n, irk, is) / real(crys%nspin, dp) * &
                            xnlkb(i) * (real(dhd(1 + ist, i, n), dp) * &
                            real(dhd(1, i, n), dp) + &
                            aimag(dhd(1 + ist, i, n)) * aimag(dhd(1, i, n)))
                    end do
                    do ist = 1, 6  
                       ssum(ist) = ssum(ist) + enl * aadot(ist) & 
                          - dtwo * &
                            bands%occup(n, irk, is) / real(crys%nspin, dp) * &
                            xnlkb(i) * (real(dhd(4 + ist, i, n), dp) * &
                            real(dhd(1, i, n), dp) + &
                            aimag(dhd(4 + ist, i, n)) * aimag(dhd(1, i, n)))
                       temp(ist)=temp(ist)  - dtwo * &
                            bands%occup(n, irk, is) / real(crys%nspin, dp) * &
                            xnlkb(i) * (real(dhd(4 + ist, i, n), dp) * &
                            real(dhd(1, i, n), dp) + &
                            aimag(dhd(4 + ist, i, n)) * aimag(dhd(1, i, n)))
                    end do
                 end do
              end do
!      call mystop
           end do               ! end of loop over spins
        end do                  ! end of loop over atoms(type)
     end do                     ! end of loop over atomic types
  end do                        ! end of loop over k-points

  !     free the dynamic arrays
  deallocate(anlga)  
  deallocate(dhd)  

  if (pspot%nanl > 0) deallocate (xnlkb)  

  return  

end subroutine forstressnloc
