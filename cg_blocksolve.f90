!     @process extchk
!
subroutine cg_blocksolve(imonitor, crys, p0, ham, maxit, delta, &
      rk, n_ukq, e_k, e_kq, ukq, uk, ukqtil, iguess, phi, nbandsfft)
  !
  !     1996 by Bernd Pfrommer and Francesco Mauri
  !
  use all_to_all_module  
  include 'use.h'  
  implicit none  
  include 'interface.h'  
  include 'all_to_all.h'  
  include 'flibcalls.ph'  
  !
  !     INPUT:
  !     -----
  !
  type(hamiltonian), intent(in) :: ham                ! the hamiltonian at k+q
  type(crystal), intent(in) :: crys   ! need the bvec information for gradient
  integer, intent(in) :: iguess, &       ! if iguess=1, compute starting guess
       p0, &                         ! direction in which the p0 operator acts
       imonitor, &               ! monitor flag indicating additional printout
       n_ukq, &                       ! number of u_k+q and u_k (must be same)
       maxit, &                                     ! max number of iterations
       nbandsfft                   ! max number of bands to FFT simultaneously
  real(dp), intent(in) :: rk(3), &                        ! the kpoint for u_k
       e_k(n_ukq)                  , &             ! energy eigenvalues of u_k
       e_kq(n_ukq)                   , &         ! energy eigenvalues of u_k+q
       delta                                           ! convergence criterion
  complex(dp), intent(in) ::  phi(ham%dim, n_ukq), &   ! contains grad|u_k>.
                                                       ! Will be destroyed!
       ukq(ham%dim, n_ukq), &         ! the full occupied wave functions u_k+q
       uk(ham%dim, n_ukq)                      ! the full occupied states at k
  !
  !     INPUT/OUTPUT:
  !     -------------
  !
  complex(dp), intent(inout) :: ukqtil(ham%dim, n_ukq)      ! the desired wave
                                                           ! function ukqtilde
  !
  !     DESCRIPTION:
  !     -----------
  !
  !     solves the equation
  !
  !     (H_k+q - e_k*I + shift* Q) |u~_k+q> = -(1-Q)grad|u_k>
  !
  !     by means of a conjugate gradient algorithm
  !
  !     Note: The ham%gspace is expected to have rk=k+q
  !     Note: If iguess /= 1, then the ukqtil should contain starting
  !           guess on input
  !
  !
  !     --------------------- local variables ----------------------
  !
  integer :: i, j, k, l, iwork, &
       nunconv, &                               ! number of unconverged states
       nconvp1, &                        ! number of converged states plus one
       nnconvp1, &                 ! new number of unconverged states plus one
       iter, &                                          ! the iteration number
       len, &                                         ! length of full vectors
       noccup, &                                   ! number of occupied states
       ilumo, &                                                ! index of lumo
       ihomo                                                   ! index of homo
  real(dp) :: dgamma, &                                                ! gamma
       rkq(3), &                                                         ! k+q
       t0, &
       a, c, &
       kq(3),&                                     ! temporary storage for k+q
       time(7)
  real(dp), allocatable :: rhoold(:), &              ! previous iterations rho
       rho(:), &                                    ! absolute of the gradient
       eprecond(:), &
       shift(:)                        ! energy shift of the Hamiltonian [Ryd]
  complex(dp), allocatable ::  g(:,:), &                         ! the gradient
       gp(:,:), &                                               ! the gradient
       h(:,:), &                                              ! step direction
       hold(:,:), &                       ! previous iterations step direction
       t(:, :), &                                           ! some work vector
       work(:), &                                                ! work arrays
       mat(:), &                                         ! intermediate matrix
       me_k(:)                                ! minus e_k, in complex(dp) form
  complex(dp) :: ssum, &                                       ! debugging sum
       tkinav, &                           ! work variable for preconditioning
       dcgamma, &                                          ! conjugation gamma
       dclambda                           ! step length lambda, in complex(dp)
  real(dp), external :: gimmetime, average_ekin  
  !
  !     ------------------------------------------------------
  !
  !      write(9,*) 'ekq:',e_kq
  !      write(9,*) 'uk on input',e_k
  !      do i=1, ham%gspace%length
  !         write(9,*) i,uk(i)
  !      end do

  t0 = gimmetime()  
  len = ham%gspace%length  
  ilumo = n_ukq + 1             ! index of lumo
  ihomo = n_ukq                 ! index of homo
  noccup = n_ukq  
  rkq = ham%gspace%rk           ! k+q
  !
  !     allocate a work array that is big enough for all purposes
  !
  !This line delated, don't think it is necc jry
  !iwork = max(gapp%fftsize(1) * gapp%fftsize(2) * gapp%fftsize(3), &
  !
  !
  iwork = max(ham%vnloc%nvecs * noccup, &        ! for apply_ham_wfn
       noccup * noccup)                          ! for apply_proj_wfn
  allocate(work(iwork))  
  work = zzero

  !
  !     -------- compute the rhs of the equation to solve. this is -------
  !
  !               psi = -(1-Q) * phi
  !
  call apply_proj_wavefn(len, dmone, 1, ukq(1, 1), noccup, phi(1, 1), &
       noccup, phi(1, 1), work(1))        ! compute (1-Q) grad(u_k)
  call mzdscal(len * noccup, dmone, phi(1, 1), 1)  
  !
  !     compute starting guess. Actually the cg_startguess routine
  !     didn't improve results so we just initialize the array
  !
  if (iguess == 1) then  
       ukqtil=zzero
  end if
  !
  !     now ukqtil contains the full trial vectors
  !
  !
  !     ------------ here starts the conjugate gradient part --------
  !
  !
  !     figure out how much the Hamiltonian must be shifted to make it
  !     positive  definite.
  !
  allocate (shift (noccup) )  
  shift(1:noccup) = (e_k(1:noccup) - e_kq(1)) + done

  allocate(g(len, noccup))           ! used as scratch array here
  allocate(gp(len, noccup))  
  allocate(h(len, noccup))  
  allocate(hold(len, noccup))  
  allocate(t(len, noccup))  
  allocate(eprecond(noccup))  
  allocate(rho(noccup))
  allocate(rhoold(noccup))  

  allocate(me_k(noccup))  
  me_k(1:noccup) = cmplx( -e_k(1:noccup), dzero, dp)  
  do i = 1, noccup  
     eprecond(i) = 1.35d0 * average_ekin(ukq(1, i), tkinav, &
          ham%gspace%ekin(1), 1, len)
  end do
  nconvp1 = 1  
  nunconv = noccup - nconvp1 + 1  
  do iter = 1, maxit  

     !
     !        compute the gradient. can reuse information from previous step
     !
     !
     if (iter == 1) then  
        !           g=H*ukqtil
        call apply_ham_wavefn(time,imonitor, g(1, nconvp1), ham, &
             ukqtil(1, nconvp1), nunconv, work(1), nbandsfft)
        !           g=g + shift *Q   * uktil
        call apply_proj_wavefn(len, shift(nconvp1), nunconv, &
             ukq(1, 1), noccup, ukqtil(1, nconvp1), nunconv, &
             g(1, nconvp1), work(1))
        !           g=g -e_k * uktil
        do i = nconvp1, noccup  
           call mzaxpy(len, me_k(i), ukqtil(1, i), 1, g(1, i), 1)
        end do
        call mzaxpy(len * nunconv, zmone, phi(1, nconvp1), 1, &
             g(1, nconvp1), 1)           !    g = g - b
     end if
     !
     !        compute residual
     !
     gp(1:len, nconvp1:noccup) = g(1:len, nconvp1:noccup)  
     nnconvp1 = nconvp1  
     do i = nconvp1, noccup  
        call precondition(gp(1, i), eprecond(i), ham%ekinmod(1), &
             ham%gspace%ekin(1), 1, len)
        rho(i) = real(parallel_zdotc(len, gp(1, i), 1, g(1, i), 1),dp)
        if ((rho(i) < delta) .and. (i == nnconvp1)) nnconvp1 = nnconvp1 + 1
     end do
     !         if((iand(imonitor,8).eq.8).and.(nconvp1.lt.nnconvp1))
     !     $        write(9,920) nconvp1,nnconvp1-1,iter
     nconvp1 = nnconvp1  

     nunconv = noccup - nconvp1 + 1  

     if (nconvp1 == noccup + 1) goto 100  
     !
     !        compute the step direction h. Conjugate it to previous step
     !
     if (iter == 1) then  
        h(1:len, nconvp1:noccup) = -gp(1:len, nconvp1:noccup)  
     else  
        do i = nconvp1, noccup  
           dgamma = rho(i) / rhoold(i)  
           h(1:len, i) = -gp(1:len, i) + dgamma * hold(1:len, i)  
        end do
     end if

     rhoold(nconvp1:noccup) = rho(nconvp1:noccup)  
     hold(1:len, nconvp1:noccup) = h(1:len, nconvp1:noccup)  
     !
     !        compute t = A*h
     !
     call apply_ham_wavefn(time,imonitor, t(1, nconvp1), ham, h(1, nconvp1), &
          nunconv, work (1), nbandsfft)     ! t = H*h
     call apply_proj_wavefn(len, shift(nconvp1), nunconv, ukq(1, 1), &
          noccup, h(1, nconvp1), nunconv, t(1, nconvp1), work(1))  ! t=t+sh
     !
     !        compute the coefficients a and c for the line minimization
     !        compute step length lambda
     do i = nconvp1, noccup  

        call mzaxpy(len, me_k(i), h(1, i), 1, t(1, i), 1)     !t=t-e_k*h
        a = real(parallel_zdotc(len, h(1, i), 1, g(1, i), 1), dp)  
        c = real(parallel_zdotc(len, h(1, i), 1, t(1, i), 1), dp)  
        dclambda = cmplx( -a / c, dzero, dp)  
        !
        !           move to new position
        !
        !           ukqtil = ukqtil+lambda*h
        !
        call mzaxpy(len, dclambda, h(1, i), 1, ukqtil(1, i), 1)  
        !           update to get the gradient
        ! g = g + lambd
        call mzaxpy(len, dclambda, t(1, i), 1, g(1, i), 1)  
     end do
     !         write(9,900) iter,rho(nconvp1:noccup)
     !         call myflush(9)
  end do
  write(9, *) '*** WARNING: cg_blocksolve not converged in', &
       maxit, ' iterations'

  write(9, *) 'residual errors:', rho  

100 continue  
  deallocate(me_k)  
  deallocate(rho)
  deallocate(rhoold)  
  deallocate(shift)
  deallocate(eprecond)  
  deallocate(g)
  deallocate(gp)  
  deallocate(h)  
  deallocate(hold)  
  deallocate(t)  
  deallocate(work)  

  if (iand (imonitor, 8) == 8) then  
     write(9, 930) iter, gimmetime() - t0  
     call myflush (9)  
  end if

900 format('iter:',i3,' rho=', 1000g12.6)  
910 format('cg_solve converged at iteration:',i3, &
       &     ' with residual error rho=', g12.6)
920 format(' STATES ',i3,'-',i3,' TOOK ',i3,' ITERATIONS ')  

930 format(' CGSOLVE TOOK ',i3,' ITERATIONS AND ', &
       &     f12.3, ' SECONDS')

end subroutine cg_blocksolve
!
!     ==================================================================
!
!     @process extchk

subroutine apply_proj_wavefn(len, mu, nff, proj, nproj, wfn, &
     nwfn, pwfn, work)
  !
  !     1996 Bernd Pfrommer
  !
  use constants
  use all_to_all_module
  implicit none  
  include 'all_to_all.h'  
  include 'flibcalls.ph'  
  !
  !     INPUT:
  !     -----
  integer, intent(in) :: len, &      ! length of wave functions and projector
       nproj, &                                        ! number of projectors
       nff, &                                        ! number of fore factors
       nwfn                                        ! number of wave functions
  real(dp), intent(in) :: mu(nff)                               ! forefactors
  complex(dp), intent(in) :: proj(len, nproj), wfn(len, nwfn)  
  !
  !     OUTPUT:
  !     ------
  !
  ! projector applied to the wavefunctions
  !
  complex(dp), intent(out) :: pwfn(len, nwfn)  
  !
  !     WORK:
  !     ----
  !
  complex(dp) :: work(nproj * nwfn)  
  !
  !     DESCRIPTION:
  !     -----------
  !
  !     applies mu*projector to a wave function wfn and ADDS ! the result
  !     to pwfn:
  !
  !     pwfn = pwfn + mu*proj*wfn
  !
  !     The work array should be set to reasonable numbers,
  !     because NaN*zero = NaN still!
  !
  !
  !     ------------------------------------------------------------------
  !
  integer :: i, j  

  call mzgemm('C', 'N', nproj, nwfn, len, zone, proj(1, 1), len, &
       wfn(1, 1), len, zzero, work(1), nproj)
  call all_sum_all (work, nproj * nwfn)  
  if (nff == nwfn) then  
     do i = 1, nwfn  
        do j = 1, nproj  
           work((i - 1) * nproj + j) = work((i - 1) * nproj + j) * mu (i)
        end do
     end do
  else if (nff == 1) then  
     work = work * mu(1)  
  else  
     write(9, *) 'apply_proj_wavefn: nff wrong:', nff  
     call mystop  
  end if

  call mzgemm('N', 'N', len, nwfn, nproj, zone, proj(1, 1) , len, &
       work(1) , nproj, zone, pwfn(1, 1), len)

  return  

end subroutine apply_proj_wavefn
