!-*-F90-*-
!
!     @process extchk optimize(2)
!
subroutine setup_nonlocal_derivative(qmag, rk, pspot, &
     k_gspace, nloc_deriv_proj, crys)

  include 'use.h'  
  implicit none              ! implicit? Just say no!
  include 'interface.h'  
  !
  !     INPUT:
  !     -----
  !
  type(pseudo_potential), intent(in) :: pspot  
  type(crystal), intent(in) :: crys  
  type(parallel_gspace), intent(in) :: k_gspace  
  real(dp), intent(in) :: rk(3), &                      ! the origin kpoint
       qmag                        ! magnitude of delta q to compute derivative
  !
  !     OUTPUT:
  !     ------
  !
  !     derivative operator
  !
  complex(dp), intent(out) :: &
       nloc_deriv_proj(k_gspace%length + pspot%nanl, pspot%nanl, 2, 3)
  !
  !
  !     DESCRIPTION:
  !     -----------
  !
  !
  !     Sets up the nonlocal potential part of the hamiltonian at
  !
  !         k+qmag*e^  for all three cartesian directions directions e^
  !                    and for the negative -e^.
  !
  !     This is required to compute the derivative dV_nloc/dk.
  !
  !     dV_nloc/dk_x |phi> =
  !       (nloc_deriv_proj(:,2,x)-nloc_deriv_proj(:,1,x)) |phi>
  !
  !     notice that there is no xnorm application necessary: it is already
  !     incorporated into the projectors
  !
  !
  !     --------------------- local variables ----------------------------
  !
  !     ----- variables for the gspace loop ---------
  !
  integer :: igv(4), fftn(4), ffth(4), igv3, irod, iord, igs  
  real(dp) :: gv(3)  
  integer :: ind, i, k, kk, l, ni, ni1, lmmin, lmmax, lm, iq, iqsign
  real(dp) :: dq(3, 3), &              ! delta q in lattice coordinates for
                                       ! cartesian displacement in x,y,z

       qi, qk(3), qcar(3), flm(16), &
       qinv,qinv1,qinv2,qinv3, fi, xni, oosq2q, t0, vq, &
       xi, xdum, xa, xb
  real(dp), parameter :: eps = 1.0d-8 
  complex(dp), allocatable :: st(:)    ! work array for structure factors
  !
  !     ------------------------------------------------------------------
  !
  allocate(st(crys%mxdatm))  

  fftn(1:3) = k_gspace%fftsize(1:3)  
  fftn(4) = k_gspace%fftsize(3)  
  ffth(:) = fftn(:) / 2  
  !
  !
  !
  dq = qmag * transpose(crys%avec) / pi2  
  oosq2q = done / sqrt(dfour * qmag)    ! 2.0 for Hart<-Ryd, 2.0 for +-q deriv
  !
  !      starts loop over g-vectors in small gspace
  !
  igs = 0  
  do iord = 1, k_gspace%lorder                 ! loop through x/y gspace
     irod = k_gspace%order(1, iord)  
     igv(1) = irod / k_gspace%fftsize(2)  
     igv(2) = mod(irod, k_gspace%fftsize(2))  
     igv(3) = k_gspace%order(2, iord)  
     igv(4) = k_gspace%order(3, iord)  
     igv(:) = mod(igv(:) + ffth(:), fftn(:)) - ffth(:)  
     gv(1:2) = real(igv(1:2), dp)  
     do igv3 = igv(3), igv(4)       ! loop over z axis
        gv(3) = real(igv3, dp)  
        igs = igs + 1  
        do iq = 1, 3                ! loop over all q-points
           do iqsign = 1, -1, -2    ! sign of q
              qk = gv(:) + rk(:) + real(iqsign, dp) * dq(:, iq)  
              qcar(1) = crys%bvec(1, 1) * qk(1) + crys%bvec(1, 2) * qk(2) + &
                   crys%bvec(1, 3) * qk(3)
              qcar(2) = crys%bvec(2, 1) * qk(1) + crys%bvec(2, 2) * qk(2) + &
                   crys%bvec(2, 3) * qk(3)
              qcar(3) = crys%bvec(3, 1) * qk(1) + crys%bvec(3, 2) * qk(2) + &
                   crys%bvec(3, 3) * qk(3)
              qi = sqrt(qcar(1) * qcar(1) + qcar(2) * qcar(2) + &
                   qcar(3) * qcar(3))
              !
              !                 ----- compute angular functions
              !
              flm(1) = done  
              if (qi > eps) then  
                 qinv1 = done / qi  
                 flm(2) = qcar(1) * qinv1 
                 flm(3) = qcar(2) * qinv1 
                 flm(4) = qcar(3) * qinv1 
                 qinv2 = qinv1* qinv1 
                 flm(5) = drt3 * qcar(1) * qcar(2) * qinv2 
                 flm(6) = drt3 * qcar(2) * qcar(3) * qinv2 
                 flm(7) = drt3 * qcar(3) * qcar(1) * qinv2 
                 flm(8) = dtrhalf * qcar(3) * qcar(3) * qinv2 - dhalf
                 flm(9) = drt3 * (qcar(1) * qcar(1) - qcar(2) * qcar(2)) * &
                      qinv2 * dhalf
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
              !                 starts loop over second index
              !
              ind = 0  
              do k = 1, crys%ntype  
                 do kk = 1, crys%natom(k)      !  compute complex phase factor
                    fi = gv(1) * crys%rat(1, kk, k) + &
                         gv(2) * crys%rat(2, kk, k) + gv(3) * crys%rat(3, kk, k)
                    st(kk) = exp(cmplx(dzero, fi, dp))  
                 end do
                 !
                 ! loop over angular momenta l
                 do l = 1, 5  
                    if (pspot%nkb(l, k) /= 0) then  
                       xni = qi / pspot%delqnl(k) + dtwo  
                       !
                       !                          cubic spline interpolation
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
                               pspot%d2vkbdq2(ni - 2, l, k) + &
                               (xb**3 - xb) * pspot%d2vkbdq2(ni1 - 2, l, k)) * &
                               dsixth
                       end if
                       lmmin = pspot%lo(l, k) * pspot%lo(l, k) + 1  
                       lmmax = (pspot%lo(l, k) + 1) * (pspot%lo(l, k) + 1)
                       !
                       ! loop over m quantum number
                       do lm = lmmin, lmmax  
                          do kk = 1, crys%natom(k)  
                             ind = ind + 1  
                             nloc_deriv_proj(igs, ind, (iqsign + 3) / 2, iq) = &
                                  conjg(st(kk) * flm(lm) * vq) * oosq2q
                          end do
                       end do
                    end if
                 end do
              end do                         ! end of loop over atomic types
           end do                            ! loop over iqsign
        enddo                                ! loop over iq
     end do                                  ! end of loop over 3rd dimension
  end do                                     ! end of loop over small gspace

  deallocate(st)  
  !
  !     attach the normalization factors to the end of the nonloc_deriv
  !     array
  ind = 0  
  do k = 1, crys%ntype  
     do l = 1, 5  
        if (pspot%nkb(l, k) /= 0) then  
           lmmin = pspot%lo(l, k) * pspot%lo(l, k) + 1  
           lmmax = (pspot%lo(l, k) + 1) * (pspot%lo(l, k) + 1)  
           do lm = lmmin, lmmax  
              do kk = 1, crys%natom(k)  
                 ind = ind + 1  
                 nloc_deriv_proj(k_gspace%length + ind, 1, 1, 1) = &
                      cmplx(real(pspot%nkb(l, k), dp), dzero, dp)
              end do
           end do
        end if
     end do
  end do

end subroutine setup_nonlocal_derivative
