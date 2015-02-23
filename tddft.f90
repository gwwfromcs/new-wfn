!
subroutine tddft(ipr, pot_gspace, k_gspace, pw_params, crys, &
     syms, ffts, kpoints, bands, wavefn, denc, den)

  use all_to_all_module
  include 'use.h'
  implicit none             ! implicit ? ....
  include 'interface.h'
  include 'all_to_all.h'
  include 'flibcalls.ph'

  !     ---------------- arguments ----------------------------------------
  !
  !     INPUT:
  !     -----
  !
  integer, intent(in) :: &
       ipr                  ! print flag
  type(symmetry), intent(in) :: &
       syms                 ! symmetry operations
  type(crystal), intent(in) :: &
       crys                 ! crystal structure
  type(fft_struc), intent(in) :: &
       ffts                 ! information for the fast fourier transform
  type(pw_parameter), intent(in) :: &
       pw_params            ! plane wave parameters
  type(kpoint), intent(in) :: &
       kpoints              ! BZ integration kpoints
  type(complex_gspace_array), intent(in) :: &
       denc                 ! core charge density
  type(parallel_gspace), intent(in) :: &
       pot_gspace, k_gspace(*)
  !
  !     INPUT/OUTPUT:
  !     -------------
  !
  type(complex_gspace_array), intent(inout) :: &
       den                  ! valence charge density on input, becomes total
  type(band), intent(inout) :: &
       bands                ! eigenvalues, occupation numbers etc
  type(complex_gspace_array), intent(inout) :: &
       wavefn               ! wave functions for all kpoints, bands
  !  
  !     DESCRIPTION:
  !     -----------
  !
  !     This subroutine performs a time-dependent density-functional
  !     theory calculation based upon the results of an SCF calculation.
  !     
  !     The current approximation for exchange-correlation kernel is
  !     adiabatic local density approximation.
  !
  !     --------------- local variables -------------------------
  !
  real(dp) :: &
       t0, tstart, &        ! start time
       rho, rs, sqrs, a0, fa, ecu, ecp, alrs, norm, rnorm, exc, rhotot, &
       fac, afmagm, rho1, rho2, dmax, dmin, cmax, cmin, g2, &
       vc, vcp, vcu, vcd, ec, zeta, fzeta, dfdz, d2fdz2, &
       fc, fcp, fcuu, fcud, fcdu, fcdd, f, &
       rvec(3), gdotr, rk(3), rkinv(3), rkdiff(3)
  integer :: &
       i, j, k, &
       ierr, &              ! error flag
       is, js, &            ! spin counters
       irk, &               ! kpoint counter
       istate, &            ! pair state counter
       istred, &            ! reduced pair state counter
       nstates, &           ! total number of pair states
       nstred, &            ! reduced number of pair states (time-inversion)
       is1, is2, ib1, ib2, &
       iskip, jskip, jstate, &
       nval, ncon, &
       ndim, &
       nb, &
       ioccup, &
       idum, &
       lwork, &
       info
  logical :: isp            ! spin polarization
  integer, allocatable :: &
       ireal(:), &          ! indicates whether a pair state derives from
                            !  real k-point states (0) or not (1)
       ikinv(:), &          ! indicates whether the inverse of the k-point
                            !  should be included (1) or not (0)
       indx(:), &           ! index for sorting eigenvalue differences
       ispin(:)             ! spin of pair state
  real(dp), allocatable :: &
       eig(:), &            ! eigenvalue difference for pair state
       occ(:), &            ! occupancy difference for pair state
       fxc(:,:,:), &        ! xc kernel
       vcoul(:), &          ! Coulomb potential in reciprocal space
       exe(:), &            ! excitation energies
       rwork(:)
  complex(dp) :: ctmp(4, 4), csti, cstj
  complex(dp), allocatable :: &
       wavefnr(:,:,:,:), &  ! realspace wavefunctions
       chd(:,:), &          ! charge density
       stater(:,:), &       ! realspace pair states
       stateg(:,:), &       ! reciprocal space pair states
       tdmatrix(:,:),     & ! matrix elements of (v_c + f_xc)
       tdsinglet(:,:), &    ! singlet matrix elements
       tdtriplet(:,:), &    ! triplet matrix elements
       tdcoulomb(:,:), &    ! coulomb matrix elements
       work(:)
  real(dp), parameter :: eps = 1.0d-9
  real(dp), parameter :: &
       g = -0.2846d0, &             ! unpolarized constants
       b1 = 1.0529d0, &
       b2 = 0.3334d0, &
       aca = 0.0311d0 * dtwo, &  
       bca = -0.048d0 * dtwo, &  
       cca = 0.0020d0 * dtwo, &  
       dca = -0.0116d0 * dtwo, &  
       gp = -0.0843d0 * dtwo, &     ! polarized constants
       b1p = 1.3981d0, &
       b2p = 0.2611d0, & 
       acap = 0.01555d0 * dtwo, & 
       bcap = -0.0269d0 * dtwo, &  
       ccap = 0.0007d0 * dtwo, & 
       dcap = -0.0048d0 * dtwo  
  real(dp), external :: gimmetime
  !
  !     -------------- here begins the code itself --------------------------
  !
  write(9, 100)
  write(9, 900)
  !
  isp = .false.
  if (bands%nspin == 2) isp = .true.
  !
  a0 = (dfour / (dthree * dthree * pi))**dthird
  fa = dtwo * (dsix / pi)**dthird
  d2fdz2 = 8.0d0 / (9.0d0 * (dtwo**dftrd - dtwo))
  rnorm = done / real(pot_gspace%fftsize(1) * pot_gspace%fftsize(2) * &
       pot_gspace%fftsize(3), dp)
  !
  !    check list of k-points for points whose inverses differ by
  !      a G-vector
  !
  allocate(ikinv(kpoints%nrk))
  ikinv = 1
  do irk = 1, kpoints%nrk
     rk = kpoints%rk(:, irk)
     rkinv = -rk
     rkdiff = rk - rkinv
     if (abs(rkdiff(1) - real(nint(rkdiff(1)), dp)) < eps .and. &
          abs(rkdiff(2) - real(nint(rkdiff(2)), dp)) < eps .and. &
          abs(rkdiff(3) - real(nint(rkdiff(3)), dp)) < eps) &
          ikinv(irk) = 0
  end do
  !
  !    write out list of k-points
  !
  write(9, 500)
  do irk = 1, kpoints%nrk
     if (ikinv(irk) == 0) then
        write(9, 510) irk, kpoints%rk(1:3, irk), ' NO'
     else
        write(9, 510) irk, kpoints%rk(1:3, irk), 'YES'
     end if
  end do
  !
  !    spin-unpolarized case
  !
  if (.not. isp) then
     !
     ! count number of pair states for 'vertical' transitions
     !
     nstred = 0 ; nstates = 0
     do irk = 1, bands%nrk
        nval = bands%ifmax(irk, 1)
        ncon = bands%min(1) - bands%ifmax(irk, 1)
        nstred = nstred + nval * ncon
        if (ikinv(irk) == 1) then                  ! time-inversion symmetry
           nstates = nstates + 2 * nval * ncon
        else
           nstates = nstates + nval * ncon
        end if
     end do
     write(9, 300) nstates * 2                     ! include spin degeneracy
     !
     !     allocate tddft matrix
     !
     allocate(ireal(nstred),indx(nstates))
     allocate(eig(nstates), exe(nstates), occ(nstates))
     allocate(tdsinglet(nstates, nstates), tdtriplet(nstates, nstates), &
          tdcoulomb(nstates, nstates))
     !
     ! allocate real-space wavefunctions
     !
     allocate(wavefnr(ffts%r_size, bands%min(1), bands%nrk, 1), stat = info)
     call alccheck('wavefnr', ffts%r_size * bands%min(1) * bands%nrk, info)
     !
     !     workspace for diagonalisation
     !
     lwork = 64 * nstates
     allocate(rwork(3 * nstates), work(lwork))
     !
     !     Fourier transform wave functions to real space
     !
     !     units are bohr^-1.5
     !
     fac = done / sqrt(crys%vcell)
     do irk = 1, bands%nrk
        ndim = k_gspace(irk)%length
        do i = 0, bands%min(1) - 1, pw_params%nbandsfft
           nb = min(pw_params%nbandsfft, bands%min(1) - i)
           call fourier_transform(-1, ffts, k_gspace(irk), &
                wavefn%data(i * ndim + 1, irk, 1), &
                wavefnr(1, i + 1, irk, 1), nb)
           call mzdscal(ffts%r_size * nb, fac, wavefnr(1, i + 1, irk, 1), 1)
        end do
     end do
     !
     !      finished with g-space wavefunctions
     !
     call destroy(wavefn)
     !
     !     add core density to valence
     !     inverse Fourier transform the charge density to real space
     !
     !     units are electrons per cell
     !
     !     allocate density, Coulomb potential and exchange-correlation kernel
     !
     allocate(fxc(ffts%r_size, 1, 2), stat = info)
     call alccheck('fxc', ffts%r_size * 2, info)
     !
     allocate(chd(ffts%r_size, 1), stat = info)
     call alccheck('chd', ffts%r_size, info)
     !
     call mzaxpy(pot_gspace%length, zone, denc%data(1, 1, 1), 1, &
          den%data(1, 1, 1), 1)
     !
     call fourier_transform(-1, ffts, pot_gspace, den%data(1, 1, 1), &
          chd(1, 1), 1)
     !
     ! find max and min value for real space charge density, check for complex
     !
     call check_array(chd(1, 1), ffts%r_size, dmax, dmin, cmax, cmin, &
          eps, ierr)
     if (ierr /= 0) then  
        write(9, 200) ierr, cmax  
     end if
     if (dmin < dzero) then
        if (ipr >= 1) call warn(272, dmin, idum, 'tddft')  
     end if
     if (ipr /= 0) then  
        write(9, 220) dmax, dmin  
        write(9, 230) cmax, cmin  
     end if
     !
     !     calculate exchange-correlation kernel within ALDA
     !
     !                            xc
     !                          dv     (r)
     !                            sigma
     !     fxc(r, sigma, tau) = ----------
     !                          drho   (r)
     !                              tau
     !
     !     units are Ry bohr^3
     !
     exc = dzero
     rhotot = dzero
     afmagm = dzero
     !
     do i = 1, ffts%r_size
        rho = real(chd(i, 1), dp) / crys%vcell    ! in bohr^-3
        rhotot = rhotot + rho
        if (rho > dzero) then
           rho1 = dhalf * rho
           rs = (dthree / (dfour * pi * rho))**dthird
           !     exchange kernel
           fxc(i, 1, 1) = -fa * dthird * rho1**dmttrd
           fxc(i, 1, 2) = dzero
           exc = exc - dtrhalf * rho / (pi * a0 * rs)
           !     correlation kernel
           if (rs >= done) then ! ---- fit to Monte Carlo data
              sqrs = sqrt(rs)
              ecu = g / (done + b1 * sqrs + b2 * rs)         ! unpolarized
              vc = ecu * ecu * (done + 7.0d0 * b1 * sqrs / dsix + &
                   dftrd * b2 * rs) / g
              fc = rs * ecu * ecu * ecu * dthird * (dfive * b1 / (1.2d1 * &
                   sqrs) + dtwo * dthird * b2 - 7.0d0 * b1 * b1 / 1.2d1 - &
                   2.3d1 * b1 * b2 * sqrs / 1.2d1 - dftrd * b2 * b2 * rs) / &
                   (rho * g * g)
              ecp = gp / (done + b1p * sqrs + b2p * rs)        ! polarized
              vcp = ecp * ecp * (done + 7.0d0 * b1p * sqrs / dsix + &
                   dftrd * b2p * rs) / gp
              fcp = rs * ecp * ecp * ecp * dthird * (dfive * b1p / (1.2d1 * &
                   sqrs) + dttrd * b2p - 7.0d0 * b1p * b1p / 1.2d1 - &
                   2.3d1 * b1p * b2p * sqrs / 1.2d1 - &
                   dftrd * b2p * b2p * rs) / (rho * gp * gp)
           else             ! ---- Gell-Mann & Brueckner
              alrs = log(rs)
              ecu = aca * alrs + bca + (cca * alrs + dca) * rs  
              vc = ecu - (aca + (cca * alrs + cca + dca) * rs) * dthird
              fc = (vc - ecu + (cca * alrs + dtwo * cca + dca) * rs * &
                   dthird * dthird) / rho
              ecp = acap * alrs + bcap + (ccap * alrs + dcap) * rs
              vcp = ecp - (acap + (ccap * alrs + ccap + dcap) * rs) * dthird
              fcp = (vcp - ecp + (ccap * alrs + dtwo * ccap + dcap) * rs * &
                   dthird * dthird) / rho
           end if
           fcuu = fc + (ecp - ecu) * d2fdz2 / rho
           fcud = fc - (ecp - ecu) * d2fdz2 / rho
           exc = exc + rho * ecu
           fxc(i, 1, 1) = fxc(i, 1, 1) + fcuu
           fxc(i, 1, 2) = fxc(i, 1, 2) + fcud
        else
           fxc(i, 1, 1) = dzero
           fxc(i, 1, 2) = dzero
        end if
     end do
     !
     call all_sum_all(rhotot)
     call all_sum_all(afmagm)
     call all_sum_all(exc)
!       write(9, 130) rhotot * crys%vcell * rnorm
!       write(9, 120) exc * crys%vcell * rnorm
     !
     deallocate(chd)
     !
     !     allocate memory for pair states
     !
     allocate(stater(ffts%r_size, nstred), stat = info)
     call alccheck('stater', ffts%r_size * nstred, info)
     !
     allocate(stateg(pot_gspace%length, nstred), stat = info)
     call alccheck('stateg', pot_gspace%length * nstred, info)
     !
     !     set up pair states in real space
     !
     !     units are bohr^-3
     !
     istate = 1 ; istred = 1
     do irk = 1, bands%nrk                              !--- loop over k-points
        iskip = ikinv(irk)
        do ib1 = 1, bands%ifmax(irk, 1)                 !--- occupied bands
           do ib2 = bands%ifmax(irk, 1) + 1, bands%min(1)  !--- unoccupied bands
              ireal(istred) = iskip
              eig(istate) = bands%energy(ib2, irk, 1) - &
                   bands%energy(ib1, irk, 1)
              occ(istate) = dhalf * (bands%occup(ib1, irk, 1) - &
                   bands%occup(ib2, irk, 1))
              eig(istate + iskip) = eig(istate)
              occ(istate + iskip) = occ(istate)
              do i = 1, ffts%r_size
                 stater(i, istred) = wavefnr(i, ib1, irk, 1) * &
                      conjg(wavefnr(i, ib2, irk, 1))
              end do
              istred = istred + 1
              istate = istate + 1 + iskip
           end do
        end do
     end do
     !
     !     sort eigenvalue differences
     !
     call sort(nstates, eig(1), indx(1))
     !
     !     matrix elements of exchange-correlation kernel
     !
     istate = 1
     do i = 1, nstred              ! loop over kl
        iskip = ireal(i)
        jstate = 1
        do j = 1, nstred           ! loop over ij
           jskip = ireal(j)
           ctmp = zzero
           do k = 1, ffts%r_size   ! loop over r
              csti = stater(k, i)
              cstj = stater(k, j)
              ctmp(1:2, 1) = ctmp(1:2, 1) + conjg(cstj) * fxc(k, 1, :) * csti
              ctmp(1:2, 2) = ctmp(1:2, 2) + cstj * fxc(k, 1, :) * csti
              ctmp(1:2, 3) = ctmp(1:2, 3) + conjg(cstj) * fxc(k, 1, :) * &
                   conjg(csti)
              ctmp(1:2, 4) = ctmp(1:2, 4) + cstj * fxc(k, 1, :) * conjg(csti)
           end do
           tdsinglet(jstate + jskip, istate) = ctmp(1, 2) - ctmp(2, 2)
           tdtriplet(jstate + jskip, istate) = ctmp(1, 2) + ctmp(2, 2)
           tdsinglet(jstate, istate + iskip) = ctmp(1, 3) - ctmp(2, 3)
           tdtriplet(jstate, istate + iskip) = ctmp(1, 3) + ctmp(2, 3)
           tdsinglet(jstate + jskip, istate + iskip) = ctmp(1, 4) - ctmp(2, 4)
           tdtriplet(jstate + jskip, istate + iskip) = ctmp(1, 4) + ctmp(2, 4)
           tdsinglet(jstate, istate) = ctmp(1, 1) - ctmp(2, 1)
           tdtriplet(jstate, istate) = ctmp(1, 1) + ctmp(2, 1)
           jstate = jstate + 1 + jskip
        end do
        istate = istate + 1 + iskip
     end do
     fac = crys%vcell * rnorm
     tdsinglet = tdsinglet * fac
     tdtriplet = tdtriplet * fac
     !
     !     Fourier transform pair states to reciprocal space
     !
     !     dimensionless
     !
     do istred = 1, nstred, pw_params%nbandsfft
        nb = min(pw_params%nbandsfft, nstred - istred + 1)
        call fourier_transform(1, ffts, pot_gspace, stateg(1, istred), &
             stater(1, istred), nb)
        call mzdscal(pot_gspace%length * nb, crys%vcell, stateg(1, istred), 1)
     end do
     !
     !     Coulomb potential
     !
     !     units are Ry
     !
     allocate(vcoul(pot_gspace%length), stat = info)
     call alccheck('vcoul', pot_gspace%length, info)
     fac = 8.0d0 * pi / crys%vcell
     do k = 1, pot_gspace%length
        g2 = pot_gspace%ekin(k)
        if (g2 > dzero) then
           vcoul(k) = fac / g2
        else
           vcoul(k) = dzero
        end if
     end do
     !
     !     matrix elements of coulomb potential
     !
     istate = 1
     do i = 1, nstred
        iskip = ireal(i)
        jstate = 1
        do j = 1, nstred
           jskip = ireal(j)
           ctmp = zzero
           do k = 1, pot_gspace%length
              csti = stateg(k, i)
              cstj = stateg(k, j)
              ctmp(1, 1) = ctmp(1, 1) + conjg(cstj) * vcoul(k) * csti
              ctmp(2, 1) = ctmp(2, 1) + cstj * vcoul(k) * csti
              ctmp(3, 1) = ctmp(3, 1) + conjg(cstj) * vcoul(k) * &
                   conjg(csti)
              ctmp(4, 1) = ctmp(4, 1) + cstj * vcoul(k) * conjg(csti)
           end do
           tdcoulomb(jstate + jskip, istate) = ctmp(2, 1)
           tdcoulomb(jstate, istate + iskip) = ctmp(3, 1)
           tdcoulomb(jstate + jskip, istate + iskip) = ctmp(4, 1)
           tdcoulomb(jstate, istate) = ctmp(1, 1)
           jstate = jstate + 1 + jskip
        end do
        istate = istate + 1 + iskip
     end do
     tdtriplet = tdtriplet + dtwo * tdcoulomb
     !
     !     sum across all processors
     !
     call all_sum_all_complex2(tdtriplet, nstates * nstates)
     call all_sum_all_complex2(tdsinglet, nstates * nstates)
     !
     ! Multiply by square root factors
     !
     do i = 1, nstates
        do j = 1, nstates
           fac = dtwo * sqrt(abs(occ(i) * occ(j) * eig(i) * eig(j)))
           tdtriplet(j, i) = tdtriplet(j, i) * fac
           tdsinglet(j, i) = tdsinglet(j, i) * fac
        end do
     end do
     !
     !     add diagonal elements
     !
     do i = 1, nstates
        tdtriplet(i, i) = tdtriplet(i, i) + cmplx(eig(i) * eig(i), dzero, dp)
        tdsinglet(i, i) = tdsinglet(i, i) + cmplx(eig(i) * eig(i), dzero, dp)
     end do
     !
     !     solve triplet
     !
     info = 0
     call mzheev('N', 'U', nstates, tdtriplet(1, 1), nstates, &
          exe(1), work(1), lwork, rwork(1), info)
     if (info /= 0) then
        write(9, *) 'ERROR in zheev: info=', info
        call mystop
     end if
     !
     !     output triplet
     !
     write(9, 400)
     write(9, 430)
     do i = 1, nstates
        if (exe(i) < dzero) then
           write(9, 450) i, 'i', sqrt(-exe(i)), eig(indx(i))
        else
           write(9, 440) i, sqrt(exe(i)), eig(indx(i))
        end if
     end do
     !
     !     solve singlet
     !
     info = 0
     call mzheev('N', 'U', nstates, tdsinglet(1, 1), nstates, &
          exe(1), work(1), lwork, rwork(1), info)
     if (info /= 0) then
        write(9, *) 'ERROR in zheev: info=', info
        call mystop
     end if
     !
     !     output singlet
     !
     write(9, 410)
     write(9, 430)
     do i = 1, nstates
        if (exe(i) < dzero) then
           write(9, 450) i, 'i', sqrt(-exe(i)), eig(indx(i))
        else
           write(9, 440) i, sqrt(exe(i)), eig(indx(i))
        end if
     end do
     !
     deallocate(tdsinglet, tdtriplet, tdcoulomb)
     !
  else
     !
     ! spin-polarized case
     !
     !
     ! count number of pair states for 'vertical' transitions
     !
     nstred = 0 ; nstates = 0
     do is = 1, bands%nspin
        do irk = 1, bands%nrk
           nval = bands%ifmax(irk, is)
           ncon = bands%min(is) - bands%ifmax(irk, is)
           nstred = nstred + nval * ncon
           if (ikinv(irk) == 1) then                    ! time-inversion symmetry
              nstates = nstates + 2 * nval * ncon
           else                                       ! except at gamma point
              nstates = nstates + nval * ncon
           end if
        end do
     end do
     write(9, 300) nstates
     !
     !     allocate tddft matrix
     !
     allocate(ireal(nstred), ispin(nstred), indx(nstates))
     allocate(eig(nstates), exe(nstates), occ(nstates))
     allocate(tdmatrix(nstates, nstates))
     !
     ! allocate real-space wavefunctions
     !
     allocate(wavefnr(ffts%r_size, bands%min(is), bands%nrk, bands%nspin), &
          stat = info)
     call alccheck('wavefnr', ffts%r_size * bands%min(is) * bands%nrk * &
          bands%nspin, info)
     !
     !     workspace for diagonalisation
     !
     lwork = 64 * nstates
     allocate(rwork(3 * nstates), work(lwork))
     !
     !     Fourier transform wave functions to real space
     !
     !     units are bohr^-1.5
     !
     fac = done / sqrt(crys%vcell)
     do is = 1, bands%nspin
        do irk = 1, bands%nrk
           ndim = k_gspace(irk)%length
           do i = 0, bands%min(is) - 1, pw_params%nbandsfft
              nb = min(pw_params%nbandsfft, bands%min(is) - i)
              call fourier_transform(-1, ffts, k_gspace(irk), &
                   wavefn%data(i * ndim + 1, irk, is), &
                   wavefnr(1, i + 1, irk, is), nb)
              call mzdscal(ffts%r_size * nb, fac, wavefnr(1, i + 1, irk, is), 1)
           end do
        end do
     end do
     !
     !      finished with g-space wavefunctions
     !
     call destroy(wavefn)
     !
     !     add core density to valence
     !     inverse Fourier transform the charge density to real space
     !
     !     units are electrons per cell
     !
     !     allocate density, Coulomb potential and exchange-correlation kernel
     !
     allocate(fxc(ffts%r_size, bands%nspin, 2), stat = info)
     call alccheck('fxc', ffts%r_size * bands%nspin * 2, info)
     !
     allocate(chd(ffts%r_size, bands%nspin), stat = info)
     call alccheck('chd', ffts%r_size * bands%nspin, info)
     !
     do is = 1, bands%nspin
        call mzaxpy(pot_gspace%length, zhalf, denc%data(1, 1, 1), 1, &
          den%data(1, 1, is), 1)
        !
        call fourier_transform(-1, ffts, pot_gspace, den%data(1, 1, is), &
             chd(1, is), 1)
        !
        ! find max and min value for real space charge density, 
        !    check for complex
        !
        call check_array(chd(1, is), ffts%r_size, dmax, dmin, cmax, cmin, &
             eps, ierr)
        if (ierr /= 0) then  
           write(9, 200) ierr, cmax  
        end if
        if (dmin < dzero) then
           if (ipr >= 1) call warn(272, dmin, idum, 'tddft')  
        end if
        if (ipr /= 0) then  
           write(9, 210) is, dmax, dmin  
           write(9, 230) cmax, cmin  
        end if
     end do
     !
     !     calculate exchange-correlation kernel within ALDA
     !
     !                            xc
     !                          dv     (r)
     !                            sigma
     !     fxc(r, sigma, tau) = ----------
     !                          drho   (r)
     !                              tau
     !
     !     units are Ry bohr^3
     !
     exc = dzero
     rhotot = dzero
     afmagm = dzero
     !
     do i = 1, ffts%r_size
        rho1 = real(chd(i, 1), dp) / crys%vcell
        rho2 = real(chd(i, 2), dp) / crys%vcell
        rho = rho1 + rho2
        rhotot = rhotot + rho
        if (rho > dzero .and. rho1 >= dzero .and. rho2 >= dzero) then
           rs = (dthree / (dfour * pi * rho))**dthird
           zeta = (rho1 - rho2) / rho
           fzeta = ((done + zeta)**dftrd + (done - zeta)**dftrd - dtwo) / &
                (dtwo**dftrd - dtwo)
           dfdz = dftrd * ((done + zeta)**dthird - (done - zeta)**dthird) / &
                (dtwo**dftrd - dtwo)
           d2fdz2 = dfour * ((done + zeta)**dmttrd + &
                (done - zeta)**dmttrd) / (9.0d0 * (dtwo**dftrd - dtwo))
           !     exchange kernel
           fxc(i, 1, 1) = -fa * dthird * rho1**dmttrd
           fxc(i, 1, 2) = dzero
           fxc(i, 2, 1) = dzero
           fxc(i, 2, 2) = -fa * dthird * rho2**dmttrd
           exc = exc - dtrhalf * rho / (pi * a0 * rs)
           !     correlation kernel
           if (rs >= done) then ! ---- fit to Monte Carlo data
              sqrs = sqrt(rs)
              ecu = g / (done + b1 * sqrs + b2 * rs)         ! unpolarized
              vc = ecu * ecu * (done + 7.0d0 * b1 * sqrs / dsix + &
                   dftrd * b2 * rs) / g
              fc = rs * ecu * ecu * ecu * dthird * (dfive * b1 / (1.2d1 * &
                   sqrs) + dtwo * dthird * b2 - 7.0d0 * b1 * b1 / 1.2d1 - &
                   2.3d1 * b1 * b2 * sqrs / 1.2d1 - dftrd * b2 * b2 * rs) / &
                   (rho * g * g)
              ecp = gp / (done + b1p * sqrs + b2p * rs)        ! polarized
              vcp = ecp * ecp * (done + 7.0d0 * b1p * sqrs / dsix + &
                   dftrd * b2p * rs) / gp
              fcp = rs * ecp * ecp * ecp * dthird * (dfive * b1p / (1.2d1 * &
                   sqrs) + dttrd * b2p - 7.0d0 * b1p * b1p / 1.2d1 - &
                   2.3d1 * b1p * b2p * sqrs / 1.2d1 - &
                   dftrd * b2p * b2p * rs) / (rho * gp * gp)
           else             ! ---- Gell-Mann & Brueckner
              alrs = log(rs)
              ecu = aca * alrs + bca + (cca * alrs + dca) * rs  
              vc = ecu - (aca + (cca * alrs + cca + dca) * rs) * dthird
              fc = (vc - ecu + (cca * alrs + dtwo * cca + dca) * rs * &
                   dthird * dthird) / rho
              ecp = acap * alrs + bcap + (ccap * alrs + dcap) * rs
              vcp = ecp - (acap + (ccap * alrs + ccap + dcap) * rs) * dthird
              fcp = (vcp - ecp + (ccap * alrs + dtwo * ccap + dcap) * rs * &
                   dthird * dthird) / rho
           end if
           fcuu = fc + fzeta * (fcp - fc) + &
                (vcp - vc) * dfdz * (dtwo - zeta) / rho + &
                (ecp - ecu) * ((done - zeta) * d2fdz2 - &
                (dtwo - zeta) * dfdz) / rho
           fcud = fc + fzeta * (fcp - fc) - &
                (vcp - vc) * dfdz * zeta / rho + &
                (ecp - ecu) * (zeta * dfdz - (done - zeta) * d2fdz2) / rho
           fcdu = fc + fzeta * (fcp - fc) - &
                (vcp - vc) * dfdz * zeta / rho + &
                (ecp - ecu) * (zeta * dfdz - (done + zeta) * d2fdz2) / rho
           fcdd = fc + fzeta * (fcp - fc) - &
                (vcp - vc) * dfdz * (dtwo + zeta) / rho + &
                (ecp - ecu) * ((done + zeta) * d2fdz2 + &
                (dtwo + zeta) * dfdz) / rho
           ec = ecu + fzeta * (ecp - ecu)
           exc = exc + rho * ec
           fxc(i, 1, 1) = fxc(i, 1, 1) + fcuu
           fxc(i, 1, 2) = fxc(i, 1, 2) + fcud
           fxc(i, 2, 1) = fxc(i, 2, 1) + fcdu
           fxc(i, 2, 2) = fxc(i, 2, 2) + fcdd
        else
           fxc(i, 1, 1) = dzero
           fxc(i, 1, 2) = dzero
           fxc(i, 2, 1) = dzero
           fxc(i, 2, 2) = dzero
        end if
     end do
     !
     call all_sum_all(rhotot)
     call all_sum_all(afmagm)
     call all_sum_all(exc)
!       write(9, 130) rhotot * crys%vcell * rnorm
!       write(9, 120) exc * crys%vcell * rnorm
     !
     deallocate(chd)
     !
     !     allocate memory for pair states
     !
     allocate(stater(ffts%r_size, nstred), stat = info)
     call alccheck('stater', ffts%r_size * nstred, info)
     !
     allocate(stateg(pot_gspace%length, nstred), stat = info)
     call alccheck('stateg', pot_gspace%length * nstred, info)
     !
     !     set up pair states in real space
     !
     !     units are bohr^-3
     !
     istate = 1 ; istred = 1
     do is = 1, bands%nspin
        do irk = 1, bands%nrk                         !--- loop over k-points
           iskip = ikinv(irk)
           do ib1 = 1, bands%ifmax(irk, is)                !--- occupied bands
              do ib2 = bands%ifmax(irk, is) + 1, bands%min(is) ! unoccupied bands
                 ireal(istred) = iskip
                 ispin(istred) = is
                 eig(istate) = bands%energy(ib2, irk, is) - &
                      bands%energy(ib1, irk, is)
                 occ(istate) = bands%occup(ib1, irk, is) - &
                      bands%occup(ib2, irk, is)
                 eig(istate + iskip) = eig(istate)
                 occ(istate + iskip) = occ(istate)
                 do i = 1, ffts%r_size
                    stater(i, istred) = wavefnr(i, ib1, irk, is) * &
                         conjg(wavefnr(i, ib2, irk, is))
                 end do
                 istred = istred + 1
                 istate = istate + 1 + iskip
              end do
           end do
        end do
     end do
     !
     !     sort eigenvalue differences
     !
     call sort(nstates, eig(1), indx(1))
     !
     !     matrix elements of exchange-correlation kernel
     !
     istate = 1
     do i = 1, nstred              ! loop over kl
        iskip = ireal(i)
        is = ispin(i)
        jstate = 1
        do j = 1, nstred           ! loop over ij
           jskip = ireal(j)
           js = ispin(j)
           ctmp = zzero
           do k = 1, ffts%r_size   ! loop over r
              csti = stater(k, i)
              cstj = stater(k, j)
              f = fxc(k, js, is)
              ctmp(1, 1) = ctmp(1, 1) + conjg(cstj) * f * csti
              ctmp(2, 1) = ctmp(2, 1) + cstj * f * csti
              ctmp(3, 1) = ctmp(3, 1) + conjg(cstj) * f * conjg(csti)
              ctmp(4, 1) = ctmp(4, 1) + cstj * f * conjg(csti)
           end do
           tdmatrix(jstate + jskip, istate) = ctmp(2, 1)
           tdmatrix(jstate, istate + iskip) = ctmp(3, 1)
           tdmatrix(jstate + jskip, istate + iskip) = ctmp(4, 1)
           tdmatrix(jstate + jskip, istate + iskip) = ctmp(1, 1)
           jstate = jstate + 1 + jskip
        end do
        istate = istate + 1 + iskip
     end do
     fac = crys%vcell * rnorm
     tdmatrix = tdmatrix * fac
     !
     !     Fourier transform pair states to reciprocal space
     !
     !     dimensionless
     !
     do istred = 1, nstred, pw_params%nbandsfft
        nb = min(pw_params%nbandsfft, nstred - istred + 1)
        call fourier_transform(1, ffts, pot_gspace, stateg(1, istred), &
             stater(1, istred), nb)
        call mzdscal(pot_gspace%length * nb, crys%vcell, stateg(1, istred), 1)
     end do
     !
     !     Coulomb potential
     !
     !     units are Ry
     !
     allocate(vcoul(pot_gspace%length), stat = info)
     call alccheck('vcoul', pot_gspace%length, info)
     fac = 8.0d0 * pi / crys%vcell
     do k = 1, pot_gspace%length
        g2 = pot_gspace%ekin(k)
        if (g2 > dzero) then
           vcoul(k) = fac / g2
        else
           vcoul(k) = dzero
        end if
     end do
     !
     !     matrix elements of coulomb potential
     !
     istate = 1
     do i = 1, nstred
        iskip = ireal(i)
        jstate = 1
        do j = 1, nstred
           jskip = ireal(j)
           ctmp = zzero
           do k = 1, pot_gspace%length
              csti = stateg(k, i)
              cstj = stateg(k, j)
              ctmp(1, 1) = ctmp(1, 1) + conjg(cstj) * vcoul(k) * csti
              ctmp(2, 1) = ctmp(2, 1) + cstj * vcoul(k) * csti
              ctmp(3, 1) = ctmp(3, 1) + conjg(cstj) * vcoul(k) * &
                   conjg(csti)
              ctmp(4, 1) = ctmp(4, 1) + cstj * vcoul(k) * conjg(csti)
           end do
           tdmatrix(jstate + jskip, istate) = &
                tdmatrix(jstate + jskip, istate) + ctmp(2, 1)
           tdmatrix(jstate, istate + iskip) = &
                tdmatrix(jstate, istate + iskip) + ctmp(3, 1)
           tdmatrix(jstate + jskip, istate + iskip) = &
                tdmatrix(jstate + jskip, istate + iskip) + ctmp(4, 1)
           tdmatrix(jstate, istate) = tdmatrix(jstate, istate) + ctmp(1, 1)
           jstate = jstate + 1 + jskip
        end do
        istate = istate + 1 + iskip
     end do
     !
     !     sum across all processors
     !
     call all_sum_all_complex2(tdmatrix, nstates * nstates)
     !
     ! Multiply by square root factors
     !
     do i = 1, nstates
        do j = 1, nstates
           fac = dtwo * sqrt(abs(occ(i) * occ(j) * eig(i) * eig(j)))
           tdmatrix(j, i) = tdmatrix(j, i) * fac
        end do
     end do
     !
     !     add diagonal elements
     !
     do i = 1, nstates
        tdmatrix(i, i) = tdmatrix(i, i) + cmplx(eig(i) * eig(i), dzero, dp)
     end do
     !
     !     solve matrix
     !
     info = 0
     call mzheev('N', 'U', nstates, tdmatrix(1, 1), nstates, &
          exe(1), work(1), lwork, rwork(1), info)
     if (info /= 0) then
        write(9, *) 'ERROR in zheev: info=', info
        call mystop
     end if
     !
     !     output energies
     !
     write(9, 420)
     write(9, 430)
     do i = 1, nstates
        if (exe(i) < dzero) then
           write(9, 450) i, 'i', sqrt(-exe(i)), eig(indx(i))
        else
           write(9, 440) i, sqrt(exe(i)), eig(indx(i))
        end if
     end do
     !
     deallocate(tdmatrix)
     deallocate(ispin)
     !
  end if
  !
  !    Deallocate workspace
  !
  deallocate(ikinv, ireal, indx)
  deallocate(eig, exe, occ)
  deallocate(stateg)
  deallocate(stater)
  deallocate(vcoul)
  deallocate(fxc)
  deallocate(work, rwork)
  deallocate(wavefnr)

  return

  !-------------------------------------------------------------------------

100 format(/'=== ENTERED TDDFT SUBROUTINE ==='/)
120 format('Exchange-correlation energy check :',f17.10)
130 format('Total number of electrons per cell:',f17.10)

200 FORMAT(' *** WARNING: TDDFT FINDS ',i6, &
       &     ' POINTS WITH COMPLEX CHARGE, MAX=',g14.4)
210 FORMAT(/,'  MAX AND MIN OF CHARGE DENSITY SPIN (',i1,') ' &
       &     ,2F10.4)
220 FORMAT(/,'  MAX AND MIN OF CHARGE DENSITY ',2F10.4)
230 FORMAT('  MAX AND MIN COMPLEX COMPONENT   ',2g10.2)  

300 format(/,' Number of pair states in calculation: ',I6)

400 format(/,' Triplet(?) states:',/)
410 format(/,' Singlet(?) states:',/)
420 format(/,' Excited states:',/)
430 format('   I    EXCITATION ENERGY  KS EIGENVALUE DIFF. (Ry)')
440 format(i4,4x,f15.12,4x,f15.12)
450 format(i4,3x,a1,f15.12,4x,f15.12)

500 format(/,' k-point summary:',/, &
         '    I                RK              DISTINCT INVERSE?')
510 format(I5,3(2X,F8.5),6X,A)

900 format('WARNING: this routine is still under development!'/)

end subroutine tddft
