!
subroutine excorr(ipr, icorr, ispn, rsize, norm, chd, vcell, exc, afmagm)

  use constants
  use all_to_all_module  
  implicit none                 ! implicit? Just say no!
  include 'all_to_all.h'  
  !
  !     ------------------------------------------------------------------
  !
  !     INPUT:
  !     ------
  !
  integer, intent(in) :: ipr, &
       ispn, &                                    ! number of spin components
       norm, &                       ! normalization: nfft(1)*nfft(2)*nfft(3)
       rsize                            ! declared dimension of the array chd
  character(len=2), intent(in) :: &
       icorr                         ! determines the type of correlation pot
  !                                       that is used.
  !                                 xa  x-alpha method.
  !                                     x-alpha must be set in the subroutine
  !                                 wi  wigners interpolation formula is used
  !                                 hl  hedin-lundqvist.
  !                                 ca  ceperly-alder.
  real(dp), intent(in) :: vcell                                 ! cell volume
  !
  !     INPUT/OUTPUT:
  !
  complex(dp), intent(inout) :: &
       chd(rsize, ispn)            ! @  charge density on input,
  !                                      in electrons/unit cell,
  !                                  @  v_xc on output, in rydberg
  !     OUTPUT:
  !
  real(dp), intent(out) :: afmagm, &               ! antiferromagnetic moment
       exc  
  !
  !     DESCRIPTION:
  !     -----------
  !
  !
  !     computes the exchange correlation potential
  !     (in rydberg) given the charge density den (in electrons
  !     per unit cell).
  !
  !     the potential replaces the density in the array den.
  !     adapted from sverre froyen plane wave program
  !     written june 11 1987. jlm
  !
  !     exc        the total exchange correlation energy given by the
  !                integral of the density times epsilon xc (rydberg).
  !
  !     parallel version, explicit typing, cmplx data type,  1995
  !
  !
  !     --------------------- local variables -------------------------
  !
  real(dp), parameter :: eps = 1.0d-9          ! maximum cmplx cd to tolerate
  real(dp) :: dmax, &
       dmin, cmax, abschd, a0, alp, ro, ro1, ro2, &
       rs, a, b, c, x, aln, ec, alrs, xdum, sqrs, &
       cmin, fa, zeta, fzeta, dfdz, &
       ecu, ecp, ect, vc, vcp, vcu, vcd, &
       g, b1, b2, aca, bca, cca, dca, &
       gp, b1p, b2p, acap, bcap, ccap, dcap, &
       ax, rs_kf, tf_kin(2), lda_exchange(2), gga_exchange(2), &
       lda_correlation, gga_correlation
  logical :: isp                                          ! spin polarization
  integer :: i, j, k, ierr, idum, is  
  !
  !     ------------------------------------------------------------------

  isp = .false.  
  if (ispn == 2) isp = .true.  
  !
  ! find max and min value for real space charge density, check for complex
  !
  do is = 1, ispn  
     call check_array(chd(1, is), rsize, dmax, dmin, cmax, cmin, eps, &
          ierr)
     if (ierr /= 0) then  
        write(9, 200) ierr, cmax  
     end if
     if (dmin < dzero) then  
        if (ipr >= 1) call warn(272, dmin, idum, 'excorr')  
     end if
     if (ipr /= 0) then  
        if (isp) then  
           write(9, 304) is, dmax, dmin  
           write(9, 106) cmax, cmin  
        else  
           write(9, 104) dmax, dmin  
           write(9, 106) cmax, cmin  
        end if
     end if
  end do

  exc = dzero  

  if (rsize <= 0) goto 999  
  !
  !     compute the antiferromagnetic moment. Core charge is assumed to be
  !     not polarized, so it need not be subtracted.
  !

  afmagm = dzero                         ! calculate antiferromagnetic moment
  if (isp) then  
     do i = 1, rsize  
        afmagm = afmagm + abs(chd(i, 1) - chd(i, 2))  
     end do
  end if
  call all_sum_all(afmagm)
  
  if (norm > 0) then  
     afmagm = afmagm / real(norm, dp)                             ! normalize
  else  
     write(9, *) 'excorr: norm=0. Something is wrong.'  
     call mystop  
  end if

  a0 = (dfour / (dthree * dthree * pi))**dthird  
  fa = dtwo * (dsix / pi)**dthird

  ect = dzero
  if ((icorr == 'wi') .or. (icorr == 'xa') .or. (icorr == 'hl')) then
     if (isp) then  
        write(9, *) 'excorr: No spin polarization for icorr=', icorr
        write(9, *) '        Implement it if you need it!'  
        call mystop  
     else  
        write(9, *) 'excorr: Correlation type icorr=', icorr  
        write(9, *) '        is implemented, but not debugged!'  
        call mystop  
     end if
  end if

  if (icorr == 'xa') then                        !        x-alpha correlation
     alp = done  
     write(9, 100) alp  
     do i = 1, rsize  
        ro = dhalf * real(chd(i, 1) + conjg(chd(i, 1)), dp) / vcell  
        chd(i, 1) = zzero
        if (ro > dzero) then  
           rs = (dthree / (dfour * pi * ro))**dthird  
           chd(i, 1) = cmplx(-dthree * alp / (pi * a0 * rs), dzero, dp)
           exc = exc + dthree * ro * real(chd(i, 1), dp) * dqtr  
        end if
     end do
  else if (icorr == 'wi') then                    !        wigner correlation
     a = -0.875529d0  
     b = 7.8d0
     do i = 1, rsize  
        ro = dhalf * real(chd(i, 1) + conjg(chd(i, 1)), dp) / vcell  
        chd(i, 1) = zzero
        if (ro > dzero) then  
           rs = (dthree / (dfour * pi * ro))**dthird  
           chd(i, 1) = cmplx(-dtwo / (pi * a0 * rs), dzero, dp)  
           exc = exc + dthree * ro * real(chd(i, 1), dp) * dqtr  
           chd(i, 1) = chd(i, 1) + cmplx(a * ((dfour * rs) * dthird + b) / &
                ((rs + b) * (rs + b)), dzero, dp)
           exc = exc + ro * a / (rs + b)  
        end if
     end do
  else if (icorr == 'hl') then                       !        hedin-lundqvist
     a = 2.1d1
     c = 0.045d0  
     do i = 1, rsize  
        ro = dhalf * real(chd(i, 1) + conjg(chd(i, 1)), dp) / vcell  
        chd(i, 1) = zzero
        if (ro > dzero) then  
           rs = (dthree / (dfour * pi * ro))**dthird  
           chd(i, 1) = cmplx(-dtwo / (pi * a0 * rs), dzero, dp)  
           exc = exc + dthree * ro * real(chd(i, 1), dp) * dqtr 
           x = rs / a  
           aln = log(done + done / x)  
           chd(i, 1) = chd(i, 1) - cmplx(c * aln, dzero, dp)  
           exc = exc - ro * c * ((done + x * x * x) * aln + x * dhalf - &
                x * x - dthird)
        end if
     end do
  else if (icorr == 'ca') then                        !        ceperly alder.
     ! unpolarized constants
     g = -0.2846d0  
     b1 = 1.0529d0  
     b2 = 0.3334d0  
     aca = 0.0311d0 * dtwo  
     bca = -0.048d0 * dtwo  
     cca = 0.0020d0 * dtwo  
     dca = -0.0116d0 * dtwo  
     ! polarized constants
     gp = -0.0843d0 * dtwo 
     b1p = 1.3981d0  
     b2p = 0.2611d0  
     acap = 0.01555d0 * dtwo 
     bcap = -0.0269d0 * dtwo  
     ccap = 0.0007d0 * dtwo  
     dcap = -0.0048d0 * dtwo  
     if (.not.isp) then  
        do i = 1, rsize  
           ro = real(chd(i, 1)) / vcell  
           if (ro > dzero) then  
              rs = (dthree / (dfour * pi * ro))**dthird  
              !     exchange potential
              chd(i, 1) = cmplx(-dtwo / (pi * a0 * rs), dzero, dp)
              exc = exc + dthree * ro * real(chd(i, 1), dp) * dqtr
              if (rs >= done) then  
                 sqrs = sqrt(rs)  
                 ec = g / (done + b1 * sqrs + b2 * rs)  
                 chd(i, 1) = chd(i, 1) + cmplx(ec * ec * (done + &
                      7.0d0 * b1 * sqrs / dsix + dfour * b2 * rs * dthird) / &
                      g, dzero, dp)
                 ect = ect + ro * ec  
                 exc = exc + ro * ec  
              else  
                 alrs = log(rs)  
                 ec = 0.0622d0 * alrs - 0.096d0 + &
                      (0.0040d0 * alrs - 0.0232d0) * rs
                 chd(i, 1) = chd(i, 1) + cmplx(ec - (0.0622d0 + &
                      (0.0040d0 * alrs - 0.0192d0) * rs) * dthird, dzero, dp)
                 ect = ect + ro * ec
                 exc = exc + ro * ec
              end if
           else  
              chd(i, 1) = zzero
           end if
        end do
     else                                          ! do spin polarized version
        do i = 1, rsize  
           ro1 = real(chd(i, 1), dp) / vcell  
           ro2 = real(chd(i, 2), dp) / vcell  
           ro = ro1 + ro2  
           if ((ro > dzero) .and. (ro1 >= dzero) .and. (ro2 >= dzero)) then  
              zeta = (ro1 - ro2) / ro                     ! spin polarization
              fzeta = ((done + zeta)**dftrd + (done - zeta)**dftrd - dtwo) / &
                   (dtwo**dftrd - dtwo)
              dfdz = dftrd * ((done + zeta)**dthird - &
                   (done - zeta)**dthird) / (dtwo**dftrd - dtwo)
              !                 exchange potentials
              chd(i, 1) = cmplx(-fa * ro1**dthird, dzero, dp)  
              chd(i, 2) = cmplx(-fa * ro2**dthird, dzero, dp)  
              rs = (dthree / (dfour * pi * ro))**dthird  
              !                 exchange energy
              exc = exc + 0.75d0 * (chd(i, 1) * ro1 + chd(i, 2) * ro2)
              !                  correlation potential and energy.
              if (rs >= done) then  
                 sqrs = sqrt(rs)  
                 ecu = g / (done + b1 * sqrs + b2 * rs)         ! unpolarized
                 vc = ecu * ecu * (done + 7.0d0 * b1 * sqrs / dsix + &
                      dftrd * b2 * rs) / g
                 ecp = gp / (done + b1p * sqrs + b2p * rs)        ! polarized
                 vcp = ecp * ecp * (done + 7.0d0 * b1p * sqrs / dsix + &
                      dftrd * b2p * rs) / gp
              else  
                 alrs = log(rs)  
                 ecu = aca * alrs + bca + (cca * alrs + dca) * rs  
                 vc = ecu - (aca + (cca * alrs + cca + dca) * rs) * dthird
                 ecp = acap * alrs + bcap + (ccap * alrs + dcap) * rs
                 vcp = ecp - (acap + (ccap * alrs + ccap + dcap) * rs) * dthird
              end if
              !                 now take spin into acoount for correlation
              vcu = vc + fzeta * (vcp - vc) + (ecp - ecu) * &
                   (done - zeta) * dfdz
              vcd = vc + fzeta * (vcp - vc) + (ecp - ecu) * &
                   (dmone - zeta) * dfdz
              ec = ecu + fzeta * (ecp - ecu)  
              chd(i, 1) = chd(i, 1) + vcu  
              chd(i, 2) = chd(i, 2) + vcd  
              exc = exc + ro * ec  
              ect = ect + ro * ec                   ! total correlation energy
           else  
              chd(i, 1) = zzero
              chd(i, 2) = zzero
           end if
        end do
     end if
  else  
     call warn(280, xdum, idum, 'excorr')  
  endif

999 continue  
  exc = exc * vcell / real(norm, dp)  
  ect = ect * vcell / real(norm, dp)  
  call all_sum_all(exc)  

  return

100 FORMAT(/,' X-ALPHA CORRELATION, X-ALPHA =',F5.3,/)  
104 FORMAT(/,'  MAX AND MIN OF CHARGE DENSITY ',2F10.4)  
106 FORMAT('  MAX AND MIN COMPLEX COMPONENT   ',2g10.2)  
200 FORMAT(' *** WARNING: EXCORR FINDS ',i6, &
       &     ' POINTS WITH COMPLEX CHARGE, MAX=',g14.4)

304 FORMAT(/,'  MAX AND MIN OF CHARGE DENSITY SPIN (',i1,') ' &
       &     ,2F10.4)

end subroutine excorr
