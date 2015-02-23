! ----------------------------------------------------------------------
!
!   computes the magnetic susceptibility in the LDA by perturbation theo
!
!   (1996)      Bernd Pfrommer and Francesco Mauri
!               UC Berkeley
!
!   --------------------------------------------------------------------
!     @process extchk
!
subroutine magnetic_suscept(qmag, g0mask, pw_params, ffts, bands, &
     pot_gspace, vion, crys, syms, pspot, kpoints)

  use all_to_all_module  
  include 'use.h'  
  implicit none  
  include 'interface.h'  
  include 'all_to_all.h'
 
  type(kpoint), intent(in) :: kpoints            ! the kpoints
  type(pw_parameter), intent(in) :: pw_params    ! energy cutoffs etc
  type(parallel_gspace), intent(in) :: pot_gspace  
  type(fft_struc), intent(in) :: ffts  
  type(crystal), intent(in) :: crys  
  type(band), intent(inout) :: bands  
  type(symmetry), intent(in) :: syms  
  type(pseudo_potential), intent(in) :: pspot  
  type(complex_gspace_array), intent(in) :: &
       vion                          ! the screened ionic potential in
  real(dp), intent(in) :: &
       g0mask(3, 3), &               ! tensor to mask G=0 component
       qmag                          ! magnitude of q
  !
  !     OUTPUT:
  !     ------
  !
  !
  !
  !     DESCRIPTION:
  !     -----------
  !
  !     computes the diamagnetic susceptibility of insulators.
  !     the loop structure is as follows
  !
  !
  !    rho = 0
  !    loop over k (kpoints)
  !	generate Gspace for k
  !       setup nonlocal derivative dV/dk @ k
  !*	compute u_k
  !	pack    u_k
  !	loop over q^  (wavevector of the applied magnetic field)
  !	    loop over +-q
  !		setup H_k+q
  !               setup nonlocal derivative dV/dk @ k+q
  !		compute extended starting guess for u_k+q
  !
  !*		compute u_k+q
  !		loop over the two directions b0 perpendicular to q
  !                    p0= q^ x b0
  !		     solve for u~_k+q,p0
  !		     loop over i ( occupied bands)
  !		        Ftrans u_k,i
  !			loop over p1
  !			   apply p1 to u~_k+q,p0,i
  !			   fourier transform u~_k+q,p0,p1,i
  !			   rho(r,b0,p1)= rho+sign(q)*u_k,i(r)*u~_k+q(r),p0,p1,i
  !			end loop over p1
  !		     end loop over i
  !		loop over b0
  !	    end loop over +-q
  !	end loop over q
  !    end loop over k
  !     -------------- local variables --------------------------------
  !
  integer :: p0, p1, alf, &  ! counters for directions of momentum operator
       iq, &                 ! qpoint counter
       nsub, &               ! number of plane waves in subspace
       i, j, &               ! counter for occupied bands
       irkstart, &           ! which kpoint to start with
       irk                   ! kpoint counter
  integer, parameter :: dq(3, 3) = (/ 1, 0, 0, 0, 1, 0, 0, 0, 1 /)  
  type(parallel_gspace) :: &
       k_gspace, &           ! the gspace for a given k
       ksub_gspace           ! the subspace for a given k
  type(hamiltonian) :: ham   ! the subspace hamiltonians
  type(dense_hamiltonian) :: &
       densham, &            ! for k
       denshamkq             ! for k+q
  real(dp) :: t0, tk, chi_avg, &
       rk(3), &              ! unshifted kpoint
       kq(3), &              ! shifted kpoint k+q
       chi_macro(3, 3), &    ! g(G=0, q.ne.0), in ppm/mol [cgs]
       binverse(3, 3)        ! inverse of the lattice vectors
  complex(dp), pointer :: &
       u_kq(:,:), &          ! the wave functions for a given k+q
       u_k(:,:), &           ! the wave functions for a given k
       nloc_deriv_k(:,:,:,:), &  !  derivative of nonlocal potential at k
       nloc_deriv_kq(:,:,:,:), & ! same, but at k+q
       phi(:), phi1(:), &    ! temporary storage
       vloc(:), &            ! local potential in realspace
       rhom_r(:,:,:), &      ! rho_magnetic in realspace rhom_r(r,a0,a1)
       u_ktil(:,:,:), &      ! u_ktilde in gspace
       u_kqtil(:,:), &       ! u_kqtilde in gspace
       graduk(:,:), &        ! the gradient of uk
       u_kqtil_r(:), &       ! u_kqtilde in rspace
       u_k_rspace(:)         ! u_k in real space
  complex(dp) :: delchi, &
       chivv(3, 3), &        ! g(G=0, q.ne.0)
       chihh(3, 3), &        ! g(G=0, q.ne.0)
       kk_termvv(3, 3), &    ! g(G=0, q=0)
       kk_termhh(3, 3)       ! g(G=0, q=0)
  real(dp), allocatable :: &
       e_kq(:)               ! the eigenvalues at u_k+q
  logical :: idiaconv        ! if diagonalization has converged
  integer :: b0, b1, iperm0, iperm1, iqsign, iguess, &
       n_ukq                 ! number of computed n_ukq. Used as dummy
  real(dp), external :: gimmetime  

 
  ! CHI in gspace
  !
  Complex(dp), Allocatable ::     chi_g(:,:,:)
                                                      


  !     ----------------------------------------------------------------
  !*$*OPTIMIZE(0)

  ham%ekinmod = pw_params%ekinmod  
  iguess = 1  

  if (iand(pw_params%optimize, 4) == 4) iguess = 0  

  if (pw_params%epsdiag > 1.0d-12) write(9, 150) pw_params%epsdiag  
  binverse = transpose(crys%avec) / pi2  
  !
  !     set up the local potential, which is the same for all kpoints
  !
  allocate(vloc(ffts%r_size))  
  call setup_local_potential(1, ffts, pot_gspace, vion%data(1, 1, 1), vloc(1))
  !
  !     allocate some other arrays
  !
  allocate(rhom_r(ffts%r_size, 3, 3))  
  rhom_r = zzero
  chivv = zzero
  chihh = zzero

  allocate(e_kq(bands%max))  
  irk = 0  
  !
  !     read checkpoint file if requested
  !
  if (iand(pw_params%input, 4) == 4) call chkpt_chi(irk, rhom_r(1, 1, 1), &
       9 * pot_gspace%r_size, pot_gspace%myproc, chivv, chihh, &
       pw_params%output(1), 2)

  irkstart = irk + 1  
  do irk = irkstart, kpoints%nrk  
     tk = gimmetime()  
     rk = kpoints%rk(:, irk)  
     write(9, 800) irk, rk  
     call myflush(9)  
     !
     !     set up the sub-gspace and the gspace for a given kpoint irk
     !
     k_gspace%gmax = sqrt(pw_params%emax)  
     k_gspace%rk = kpoints%rk(:, irk)  
     k_gspace%nproc = pot_gspace%nproc                ! truly parallel
     k_gspace%name = 'kmag'  
     k_gspace%myproc = pot_gspace%myproc
     k_gspace%fftsize = pot_gspace%fftsize  
     k_gspace%istar = .false.
     k_gspace%igvec = .false.  
     k_gspace%iekin = .true.
     k_gspace%iinverse = .false.  
     k_gspace%istruc = .false.  
     k_gspace%ipackinfo = .true.  
!     call generate_gspace(pw_params%output(1), k_gspace, crys, syms, 1)  
     call generate_gspace(pw_params%output(1), k_gspace, crys, syms)  
     !
     !     allocate space for wave function u_k
     !
     allocate(u_k(k_gspace%length, bands%max))  
     allocate(u_kq(k_gspace%length, bands%max))  
     allocate(phi(k_gspace%length))  
     allocate(phi1(k_gspace%length))  
     ksub_gspace%gmax = sqrt(pw_params%emaxsub)  
     ksub_gspace%rk = k_gspace%rk     ! inherit kpoint from big brother

     ksub_gspace%nproc = 1            ! scatter across all processors
     ksub_gspace%name = 'kmagsub'  
     ksub_gspace%myproc = 0
     ksub_gspace%fftsize = (/ 0, 0, 0, 0, 0 /)  
     ksub_gspace%istar = .false.
     ksub_gspace%igvec = .true.  
     ksub_gspace%iekin = .true.
     ksub_gspace%iinverse = .false.  
     ksub_gspace%istruc = .false.  
     ksub_gspace%ipackinfo = .false.  
!     call generate_gspace(pw_params%output(1), ksub_gspace, crys, syms, 1)
     call generate_gspace(pw_params%output(1), ksub_gspace, crys, syms)
     !
     !        set up the nonlocal part of H and the derivative for given k
     !
     call create_hamiltonian(ham, k_gspace, pspot, vloc, ffts)  
     t0 = gimmetime()  
     call setup_nonlocal_potential(ham, pspot, crys,1)  
     if (iand(pw_params%output(1), 8) == 8) write(9, 922) gimmetime() - t0
     !
     !        generate starting guess for u_k
     !
     call create_dense_hamiltonian(densham, ksub_gspace)  

     call start_vectors(0, pw_params, ham, densham, bands%energy(1, irk, 1), &
          bands%max, u_k(1, 1), pot_gspace, ksub_gspace, vion%data(1, 1, 1), &
          crys)
     !
     !        compute the u_k
     !
     bands%nband(irk, 1) = bands%nband(irk, 1) - 1  
     if (pw_params%bandgap > dzero) then  

!P. Zhang
!        bands%nband(irk, 1) = bands%min
!        bands%energy(bands%min + 1, irk, 1) = &
!             bands%energy(bands%min, irk, 1) + pw_params%bandgap

        bands%nband(irk, 1) = bands%min(1)  

        bands%energy(bands%min(1) + 1, irk, 1) = &
             bands%energy(bands%min(1), irk, 1) + pw_params%bandgap
     end if
     call diagonalize_hamiltonian(pw_params, pw_params%maxitdiag, &
          idiaconv, bands%min(1), bands%nband(irk, 1), pw_params%epsdiag, &
          bands%energy(1, irk, 1), ham, u_k(1, 1))
     !
     !        setup nonlocal derivative at k
     !
     if (iand(pw_params%optimize, 1) /= 1) then  
        allocate(nloc_deriv_k(k_gspace%length + pspot%nanl, pspot%nanl, 2, 3))
        call setup_nonlocal_derivative(pw_params%output(1), qmag, rk(1), &
             pspot, k_gspace, nloc_deriv_k(1, 1, 1, 1), crys)
     else  
        allocate(nloc_deriv_k(1, 1, 1, 1))  
     end if
     !
     !        compute for the case of q=0
     !
     if (iand(pw_params%optimize, 4) == 4) then  
        allocate(u_ktil(ham%gspace%length, bands%ifmax(irk, 1), 3))
     else  
        allocate(u_ktil(ham%gspace%length, bands%ifmax(irk, 1), 1))
     end if
     call magnetic_kkterm(ham, densham, nloc_deriv_k(1, 1, 1, 1), &
          qmag, rk(1), pspot, u_k(1, 1), ham%gspace%length, ksub_gspace, &
          bands%energy(1, irk, 1), bands%ifmax(irk, 1), pw_params, crys, &
          done, u_ktil(1, 1, 1), kk_termvv(1, 1), kk_termhh(1, 1))
     !
     if (iand(pw_params%optimize, 4) /= 4) deallocate(u_ktil)  
     !
     !        free to optimize for memory. then recompute below
     !
     if (iand(pw_params%optimize, 1) == 1) deallocate(nloc_deriv_k)  
     !
     !        now compute u_k+q
     !
     call create_dense_hamiltonian(denshamkq, ksub_gspace)  
     do iq = 1, 3       ! loop over all q-points
        do iqsign = 1, -1, -2          ! sign of q
           kq = kpoints%rk(:, irk) + real(iqsign, dp) * qmag * &
                matmul(binverse, dq(:, iq))
           !
           ! Generate starting guess for u_k+q. Use same g-space as for u_k
           ! but modify its kinetic energy appropriately. Since the hamiltonian
           ! refers to that g-space, it is altered implicitly, too.

           k_gspace%rk = kq
           ksub_gspace%rk = kq
           call compute_ekin_gspace(k_gspace, crys)  
           call compute_ekin_gspace(ksub_gspace, crys)  
           t0 = gimmetime()  
           call setup_nonlocal_potential(ham, pspot, crys,1)  

           if (iand(pw_params%output(1), 8) == 8) write(9, 922) gimmetime() - t0

           if (iand(pw_params%optimize, 2) == 2) then   ! don't compute star
              u_kq = u_k  
              nsub = densham%gspace%length  
              denshamkq%neig = densham%neig  
              if (.not. associated(denshamkq%matrix)) then  
                 allocate(denshamkq%matrix(size(densham%matrix)))  
              end if
              denshamkq%matrix(:) = densham%matrix(:)  
              denshamkq%energy(1:nsub) = densham%energy(1:nsub)  
              e_kq(1:bands%max) = densham%energy(1:bands%max)  
              n_ukq = bands%min(1)
           else  
              call start_vectors(0, pw_params, ham, denshamkq, e_kq(1), &
                   bands%max, u_kq(1, 1), pot_gspace, ksub_gspace, &
                   vion%data(1, 1, 1), crys)
              n_ukq = bands%max - 1  
           end if
           !
           !           compute the u_k+q
           !
           if (pw_params%bandgap > dzero) then  
              n_ukq = bands%min(1)  
              e_kq(bands%min(1) + 1) = e_kq(bands%min(1)) + pw_params%bandgap  
           end if

           call diagonalize_hamiltonian(pw_params, pw_params%maxitdiag, &
                idiaconv, bands%min(1), n_ukq, pw_params%epsdiag, e_kq(1), ham, &
                u_kq(1, 1))
           !
           !           setup nonlocal potentials and derivatives
           !
           allocate(nloc_deriv_kq(k_gspace%length + pspot%nanl, pspot%nanl, &
                2, 3))

           call setup_nonlocal_derivative(pw_params%output(1), qmag, kq(1), &
                pspot, k_gspace, nloc_deriv_kq(1, 1, 1, 1), crys)
           !
           ! if memory optimization is switched on, we have to recompute
           ! the nonlocal derivative at k.
           if (iand(pw_params%optimize, 1) == 1) then  
              allocate(nloc_deriv_k(k_gspace%length + pspot%nanl, pspot%nanl, &
                   2, 3))
              call setup_nonlocal_derivative(pw_params%output(1), qmag, rk(1), &
                   pspot, k_gspace, nloc_deriv_k(1, 1, 1, 1), crys)
           end if
           !
           !           compute b0, p0 from q
           !
           allocate(u_kqtil(k_gspace%length, bands%ifmax(irk, 1)))  
           do iperm0 = 1, -1, -2  
              b0 = mod(iq + iperm0 + 2, 3) + 1  ! gives different b0 for iperm
              p0 = mod(iq - iperm0 + 2, 3) + 1  ! gives p0 = abs(q x b0)
              !
              !              compute the gradient of uk
              !
              allocate(graduk(k_gspace%length, bands%ifmax(irk, 1)))  
              do i = 1, bands%ifmax(irk, 1)  
                 call take_nonloc_deriv(p0, k_gspace, rk(1), u_k(1, i), &
                      ham%vnloc%nvecs, nloc_deriv_k(1, 1, 1, 1), &
                      graduk(1, i), crys)

                 call take_nonloc_deriv(p0, ham%gspace, kq(1), u_k(1, i), &
                      ham%vnloc%nvecs, nloc_deriv_kq(1, 1, 1, 1), phi(1), crys)
                 graduk(:, i) = dhalf * (graduk(:, i) + phi(:))  
              end do
              !
              !              compute u_kq~ for a given p0
              !
              if (iand(pw_params%optimize, 1) == 1) then  
                 deallocate(nloc_deriv_kq)  
                 deallocate(nloc_deriv_k)  
              end if

              if (iand(pw_params%optimize, 4) == 4) u_kqtil(:,:) = &
                   u_ktil(:,:, p0)
              ham%shift = dzero                ! no shift here
              call cg_blocksolve(pw_params%output(1), crys, p0, ham, &
                   pw_params%maxitcgsol, pw_params%epsmag, ksub_gspace, &
                   rk(1), bands%ifmax(irk, 1), ksub_gspace%length, &
                   bands%energy(1, irk, 1), densham%energy(1), e_kq(1), &
                   denshamkq%energy(1), denshamkq%matrix(1), &
                   u_kq(1, 1), u_k(1, 1), u_kqtil(1, 1), iguess, graduk(1, 1), &
                   pw_params%nbandsfft)
              deallocate(graduk)  
              !
              !           re-setup nonlocal potentials and derivatives
              !
              if (iand(pw_params%optimize, 1) == 1) then  
                 allocate(nloc_deriv_kq(k_gspace%length + pspot%nanl, &
                      pspot%nanl, 2, 3))
                 call setup_nonlocal_derivative(pw_params%output(1), qmag, kq(1), &
                      pspot, k_gspace, nloc_deriv_kq(1, 1, 1, 1), crys)
                 allocate(nloc_deriv_k(k_gspace%length + pspot%nanl, &
                      pspot%nanl, 2, 3))
                 call setup_nonlocal_derivative(pw_params%output(1), qmag, rk(1), &
                      pspot, k_gspace, nloc_deriv_k(1, 1, 1, 1), crys)
              end if

              allocate(u_k_rspace(ffts%r_size))  
              do i = 1, bands%ifmax(irk, 1)         ! loop over occupied bands
                 !  tsum1(iq)= tsum1(iq) + e_kq(i)*dble(iqsign)/(4.d0*qmag
                 !
                 !  Fourier transform u_(k,i)
                 !
                 call fourier_transform(-1, ffts, k_gspace, u_k(1, i),  &
                      u_k_rspace(1), 1)
                 do alf = 1, 3  
                    !
                    !                    apply operator p1 to u_k+q~
                    !
                    call take_nonloc_deriv(alf, k_gspace, kq(1), &
                         u_kqtil(1, i), 0, nloc_deriv_kq(1, 1, 1, 1), &
                         phi(1), crys)
                    !
                    !                    add to the "charge" density
                    !
                    call add_to_charge(ffts, k_gspace, u_k_rspace(1), phi(1), &
                         real(iperm0, dp) * real(iqsign, dp) * kpoints%w(irk), &
                         rhom_r(1, b0, alf))
                 end do

                 allocate(u_kqtil_r(ffts%r_size))  

                 call fourier_transform(-1, ffts, k_gspace, u_kqtil(1, i), &
                      u_kqtil_r(1), 1)
                 do alf = 1, 3  
                    !
                    !                    apply operator p1 to u_k
                    !
                    call take_nonloc_deriv(alf, k_gspace, rk(1), u_k(1, i), &
                         0, nloc_deriv_k(1, 1, 1, 1), phi(1), crys)
                    !
                    !                    add to the "charge" density
                    !
                    call add_to_charge(ffts, k_gspace, u_kqtil_r(1), phi(1), &
                         dmone * real(iperm0, dp) * real(iqsign, dp) * &
                         kpoints%w(irk), rhom_r(1, b0, alf))
                 end do

                 deallocate(u_kqtil_r)  
                 !
                 !                 treat the G=0 component
                 !
                 do iperm1 = 1, -1, -2  
                    b1 = mod(iq + iperm1 + 2, 3) + 1   ! gives different b1 for
                    p1 = mod(iq - iperm1 + 2, 3) + 1   ! gives p1 = abs(q x b1)
                    !
                    !                    --- half-half for G=0 ----
                    !
                    call take_nonloc_deriv(p1, k_gspace, rk(1), u_kqtil(1, i), &
                         0, nloc_deriv_k(1, 1, 1, 1), phi(1), crys)
                    call take_nonloc_deriv(p1, k_gspace, kq(1), u_kqtil(1, i), &
                         0, nloc_deriv_kq(1, 1, 1, 1), phi1(1), crys)

                    phi(:) = dhalf * (phi(:) + phi1(:))  
                    delchi = real(iperm1, dp) * real(iperm0, dp) * &
                         kpoints%w(irk) * (parallel_zdotc2(k_gspace%length, &
                         u_k(1, i), 1, phi(1), 1) - kk_termhh(p0, p1))
                    if (b0 == b1) then  
                       chihh(b0, b1) = chihh(b0, b1) + delchi  
                    else  
                       chihh(b0, b1) = chihh(b0, b1) + dtwo * delchi  
                    end if
                    !
                    !                    --- dH/dk  dH/dk for G=0 ----
                    !
                    call take_nonloc_deriv(p1, k_gspace, rk(1), u_kqtil(1, i), &
                         pspot%nanl, nloc_deriv_k(1, 1, 1, 1), phi(1), crys)
                    call take_nonloc_deriv(p1, k_gspace, kq(1), u_kqtil(1, i), &
                         pspot%nanl, nloc_deriv_kq(1, 1, 1, 1), phi1(1), crys)

                    phi(:) = dhalf * (phi(:) + phi1(:))  

                    delchi = real(iperm1, dp) * real(iperm0, dp) * &
                         kpoints%w(irk) * (parallel_zdotc2(k_gspace%length, &
                         u_k(1, i), 1, phi(1), 1) - kk_termvv(p0, p1))
                    if (b0 == b1) then  
                       chivv(b0, b1) = chivv(b0, b1) + delchi  
                    else  
                       chivv(b0, b1) = chivv(b0, b1) + dtwo * delchi  
                    end if
                 end do                    ! loop over permutation iperm1
              end do                 ! end of loop over bands i
              deallocate(u_k_rspace)  
              call myflush(9)  
           end do              ! loop over permutation iperm0
           deallocate(u_kqtil)  
           deallocate(nloc_deriv_kq)  
           if (iand(pw_params%optimize, 1) == 1) then  
              deallocate(nloc_deriv_k)  
           end if
        end do           ! end of iqsign
     end do        ! end of loop over q-points

     deallocate(phi)  
     deallocate(phi1)  
     if (iand(pw_params%optimize, 1) /= 1) then  
        deallocate(nloc_deriv_k)  
     end if

     if (iand(pw_params%optimize, 4) == 4) deallocate(u_ktil)  
     call destroy(densham)  
     call destroy(denshamkq)  
     call destroy(ksub_gspace)  
     call destroy(k_gspace)  
     deallocate(u_k)  
     deallocate(u_kq)  

     call destroy(ham)  
     !         write(9,*) tsum
     !         write(9,*) 'tsum1:'
     !         write(9,*) tsum1
     !         write(9,*) 'tsum1-tsum:'
     !         write(9,*) tsum1-tsum
     write(9, 926) gimmetime() - tk  
     !
     !        do checkpointing if required
     !
     if (iand(pw_params%miscflag, 2) == 2) call chkpt_chi(irk, &
          rhom_r(1, 1, 1), 9 * pot_gspace%r_size, pot_gspace%myproc, chivv, &
          chihh, pw_params%output(1), 1)

  end do     ! end of loop over kpoints
  deallocate(vloc)  
  deallocate(e_kq)  
  Allocate(chi_g(pot_gspace%length,3,3))  
  !
  !        spin average q  Hartree
  !
  chivv = dmtwo * dhalf * dtwo * real(chivv(:,:), dp) / (qmag * qmag)
  !
  !     first chi with v-v
  !
  chi_macro = real(chivv, dp) * 0.529177d-8**3 * 0.6022d24 / &
       137.036**2 * 1.0d6

  call symmetrize_tensor(chi_macro(1, 1), crys, syms)  

  write(9, 190)
  write(9, 200)
  call printdmatrix(chi_macro(1, 1), 3, 3)

  chi_avg = (chi_macro(1, 1) + chi_macro(2, 2) + chi_macro(3, 3)) * dthird
  write(9, 400) chi_avg  

  chi_avg = (chi_macro(1, 1) * g0mask(1, 1) + chi_macro(2, 2) * &
       g0mask(2, 2) + chi_macro(3, 3) * g0mask(3, 3)) * dthird
  write(9, 410) chi_avg / (crys%vcell * 0.529177d-8**3 * 0.6022d24) * dtwo * pi2
  !
  !     now chi computed with half-half
  !
  chihh = dmtwo * dhalf * dtwo * real(chihh(:,:), dp) / (qmag * qmag)

  chi_macro = real(chihh, dp) * 0.529177d-8**3 * 0.6022d24 / &
       137.036**2 * 1.0d6

  call symmetrize_tensor(chi_macro(1, 1), crys, syms)  

  write(9, 210)
  call printdmatrix(chi_macro(1, 1), 3, 3)  

  chi_avg = (chi_macro(1, 1) + chi_macro(2, 2) + chi_macro(3, 3)) * dthird
  write(9, 400) chi_avg  

  chi_avg = (chi_macro(1, 1) * g0mask(1, 1) + chi_macro(2, 2) * &
       g0mask(2, 2) + chi_macro(3, 3) * g0mask(3, 3)) * dthird
  write(9, 415) chi_avg / (crys%vcell * 0.529177d-8**3 * 0.6022d24) * dtwo * pi2
  !
  !     compute the  magnetic susceptibilty
  !
  write(9, 350) (g0mask(i, 1:3), i = 1, 3)  

  call wrapup_chi(crys, ffts, pot_gspace, rhom_r(1, 1, 1), &
       qmag, g0mask(1, 1), chi_g(1, 1, 1), chihh(1, 1))

  call symmetrize_tensor_global(pw_params%output(1), chi_g(1, 1, 1), &
       syms, pot_gspace)

  call writegsdat(1, pot_gspace, chi_g(1, 1, 1), 9, 9, 1, 'CHI', 3)

  call nmr_shift(pw_params, crys, pot_gspace, chi_g(1, 1, 1))  
  write(9, 420) qmag  

  write(9, 300)  

  deallocate(rhom_r)  
  Deallocate(chi_g)   
  return  

150 format(/' *** WARNING: ACCURACY OF DIAGONALIZATION ', &
       &     'MIGHT NOT BE HIGH ENOUGH:',g12.6)

190 format(/' *** WARNING: CHI IS ONLY CORRECT IF SYMMETRY-', &
       &       /' *** OPERATIONS MAP CARTESIAN AXIS INTO EACH OTHER')
200 format(/' MACROSCOPIC CHI (nonlocal velocity operator)', &
       &     /' [10^-6 cm^3/mole]:'/)
210 format(/' MACROSCOPIC CHI (half-half)', &
       &     /' [10^-6 cm^3/mole]:'/)
220 format(/' G=0 CORRECTION TO THE NMR SHIFT', &
       &     '(nonlocal velocity operator):')
222 format(/' G=0 CORRECTION TO THE NMR SHIFT', &
       &     '(half-half):')
300 format(/' THE NMR SHIFT IS COMPUTED WITH HALF/HALF AND', &
       &     ' THE G=0 ALSO WITH HH')
350 format(/' THE G=0 COMPONENT FOR THE SHIFT IS MASKED WITH:'/, &
       &     3(3f12.6/))

400 format(/' AVERAGE CHI:', f12.4,' *10^-6 cm^3/mole')  
410 format(/' VV SHIFT CONTRIBUTION:', f12.4)  
415 format(/' HH SHIFT CONTRIBUTION:', f12.4)  

420 format(/' NMR Q:', f10.6)  

800 format(//' COMPUTING CONTRIBUTION FROM KPOINT',i4,':',3f9.4)  
910 format(' TIME FOR DIAGONALIZATION:', f12.3)  
922 format(' TIME FOR NONLOCAL POTENTIAL:', f12.3)  
926 format(' TIME FOR KPOINT:', f12.3)  
930 format(' TIME [SECONDS] FOR ADD_CHARGE:',f12.3)  

end subroutine magnetic_suscept
