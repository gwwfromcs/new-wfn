!
!     find the lowest eigenvectors and eigenvalues for a given
!     Hamiltonian with a conjugate gradient minimization scheme
!
subroutine diagonalize_hamiltonian(time,iter,num_HPsi,pw_params,neigmax,maxit,&
      iconv, neigmin, neig, delta, eval, ham, x)

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
  type(hamiltonian), intent(inout) :: ham    ! the hamiltonian to diagonalize
  type(pw_parameter), intent(in) :: pw_params            ! various parameters

  integer, intent(in) :: maxit, &                  ! max number of iterations
       neigmin, &               ! minimum number of eigenvalues to be computed
       neigmax, &               ! maximum number of eigenvalues to be computed
       iter                     ! for orthogonalizing on first step
                                ! needed for wavefunction extrapolation
  real(dp), intent(in) :: delta                 ! accuracy of diagonalization
  !
  !     INPUT/OUTPUT:
  !     ------------
  !
  integer, intent(inout) :: neig,&    ! max/real number of eigenvalues computed
            num_HPsi                  ! op count of H*Psi  

  complex(dp), intent(inout) :: &
       x(*)                      ! start/final wave functions on input/output
  real(dp), intent(inout) :: &
       eval(neig + 1)                 ! start/final eigenvalues, eval(neig+1)
  real(dp), intent(inout):: time(num_time)
  !
  !     OUTPUT:
  !     -------
  !
  logical, intent(out) :: iconv               ! true if converged, false else
  !
  !
  !     --------------------- local variables ----------------------
  !
  real(dp) :: t0, t1, t2,t3, &                              ! timing variables
       oldrho, a, b, c, d, demax, de, lambda, eprecond, delrho,rho_init, &
       rho, &                                                 ! current error
       nops, &                                         ! number of operations
       egap                ! gap between highest computed and next eigenvalue
  real(dp)   num,energy_old,energy1,denom !DBR
  integer :: myproc, ierrmall, mnew, mmax,nover_cnt,&
       diaginfo, &                  ! return value of diagonalization routine
       m, &                                       ! current number of vectors
       m_squared, &                                             ! square of m
       i, j, &                                                      ! dummies
       imax, &                                            ! largest gap index
       n, &                                                ! iteration number
       nreset, &                                 ! number of resets performed
       len                                                    ! gspace length
  complex(dp) :: dclambda, dclambdasq, dcgamma
  !
  !     work arrays
  !
  complex(dp), allocatable :: &                               ! of size m x m
       eps(:), s(:), xsi(:), eta(:), nu(:), mu(:)
  complex(dp), allocatable :: &                            ! of size nanl x m
       ews(:,:)  
  complex(dp), allocatable :: &                             ! of size len x m
       g(:,:), gold(:,:), gp(:,:), h(:,:), y(:,:)
  integer, allocatable :: &                                       ! size of m
       iord(:)
  real(dp), allocatable :: rwork(:)  

  logical :: ishift, ireset, iconjug 
  !
  !          local routines
  real(dp), external :: trace, contract, average_ekin, gimmetime  

  t0 = gimmetime()  
  if (iand(pw_params%output(1), 8) == 8) t3 = t0
  m = neig                                       ! start with maximum allowed
  m_squared = m * m  
  len = ham%gspace%length  
  myproc = ham%gspace%myproc  
  !
  !     allocate all the work spaces
  !
  allocate(eps(m * m)) ; eps = zzero 
  allocate(s(m * m)) ; s = zzero
  allocate(xsi(m * m)) ; xsi = zzero  
  allocate(eta(m * m)) ; eta = zzero  
  allocate(nu(m * m)) ; nu = zzero  
  allocate(mu(m * m)) ; mu = zzero  
  allocate(rwork(max(1, 3*m-2))) ; rwork = dzero  
  if (ham%vnloc%nvecs > 0 .and. .not. ham%pspot%NL_rspace(1) ) then  
     allocate(ews(ham%vnloc%nvecs, m))  
     ews = zzero
  end if
  allocate(g(len, m), stat = ierrmall)
  if (ierrmall /= 0) then  
     write(9, 930) 'g', len, m, ierrmall  
     call mystop  
  end if
  g = zzero  
  allocate(gold(len, m), stat = ierrmall)
  if (ierrmall /= 0) then  
     write(9, 930) 'gold', len, m, ierrmall  
     call mystop  
  end if
  gold = zzero
  allocate(h(len, m), stat = ierrmall)
  if (ierrmall /= 0) then  
     write(9, 930) 'h', len, m, ierrmall  
     call mystop  
  end if
  h = zzero
  allocate(y(len, m), stat = ierrmall)
  if (ierrmall /= 0) then  
     write(9, 930) 'y', len, m, ierrmall  
     call mystop  
  end if
  y = zzero  
  allocate(gp(len, m), stat = ierrmall)
  if (ierrmall /= 0) then  
     write(9, 930) 'gp', len, m, ierrmall  
     call mystop  
  end if
  gp = zzero  

  allocate(iord(m))
  iord = 0  

  ireset = .true.                                ! first step is like a reset
  iconjug = .false.                      ! don't make conjugate in first step
  if (neig .gt. neigmin) then 
    ishift = .true.               ! want to shift if error has decreased enough
  else
    ishift = .false.
  end if
  rho = 1.0d4                                       ! give some initial error
  nops = 0  

  t1 = dzero
  !      write(9,*) 'eigenvalues from submatrix:', eval(1:neig+1)
  !      write(9,*) 'starting vectors:'
  !      call printmatrix(x,len,m)
  !      call overlap_matrix(s,x,x,m,len) !     s = x^T*x
  !      write(9,*) 's matrix:'
  !      call printmatrix(s,m,m)
  !      call apply_ham_wavefn(pw_params%output,y(1,1),ham,x(1),
  !     $     m,ews(1,1),pw_params%nbandsfft)          ! y = H * x
  !      call overlap_matrix(eps,x,y,m,len) ! eps = x^T*y
  !      write(9,*) 'eps matrix:'
  !      call printmatrix(eps,m,m)
  !      m = 4
  !      write(9,*) 'number of states computed:', m
  !      egap=0.283134674988261281-(-0.588457289928866700)
  !      ham%shift = -0.588457289928866700 +
  !     $     egap/4 + 0.2

  if (iter .eq. 1) then
        call cheap_overlap_matrix(s(1), x(1), m, len) 
        call occsp_sqrtinv_scal('U','n',myproc,s,m,ham%gspace%nproc)

        call mztrmm('R','U','N','N',len,m,zone,s(1),m,x(1),len)   
  end if

  egap = eval(neig + 1) - eval(neig)  

  if (egap < dzero) then                        ! no good gap as input, guess
     write(9, *) ' *** DIAGONALIZATION: INPUT EIGENVALUES BAD!'  
     ham%shift = dzero 
     call apply_ham_wavefn(time,pw_params%output(1), y(1, 1), ham, x(1), &
          m, ews(1, 1), pw_params%nbandsfft)                      ! y = H * x 
     num_HPsi=num_HPsi+1
     call overlap_matrix(eps(1), x(1), y(1, 1), m, len)         ! eps = x^T*y
     egap = dmfour * eps(neigmin + (neigmin - 1) * m)  
     write(9, *) ' *** GUESSED GAP:', egap  
  end if
  !
  !     shift hamiltonian just enough to make it negative definite
  !
  ham%shift = eval(neig) + egap * dqtr + pw_params%shiftsafety  

  nreset = 0
  nover_cnt=0
  do n = 1, maxit                   ! ----------- start of big iteration loop
100  continue                                                       ! restart
     if (ireset) then  
       if (iand(pw_params%output(1), 8) == 8) call get_timing(time(4),t3) 
        call apply_ham_wavefn(time,pw_params%output(1), y(1, 1), ham, x(1), &
             m, ews(1, 1), pw_params%nbandsfft)                     ! y=H * x  
        num_HPsi=num_Hpsi+1
        if (iand(pw_params%output(1), 1) == 1) t2 = gimmetime()
        if (iand(pw_params%output(1), 8) == 8) t3 = gimmetime()      
        call cheap_overlap_matrix(s(1), x(1), m, len)             ! s = x^T*x
        call overlap_matrix(eps(1), x(1), y(1, 1), m, len)      ! eps = x^T*y
        if (iand(pw_params%output(1), 1) == 1 ) then  
           t1 = t1 + gimmetime() - t2  
           nops = nops + 1.2d1 * real(m * m * len, dp)
        end if
       if (iand(pw_params%output(1), 8) == 8) call get_timing(time(6),t3)
        ireset = .false.                              ! no resetting any more

     else  

        dclambdasq = dclambda * dclambda  
        call matrix_complete(eta(1), m)  
        call matrix_complete(nu(1), m)  
        call mzaxpy(m * m, dclambda, eta(1), 1, s(1), 1)  
        call mzaxpy(m * m, dclambdasq, nu(1), 1, s(1), 1)  
        call matrix_complete(xsi(1), m)  
        call matrix_complete(mu(1), m)  
        call mzaxpy(m * m, dclambda, mu(1), 1, eps(1), 1)  
        call mzaxpy(m * m, dclambdasq, xsi(1), 1, eps(1), 1)  

        num=dzero
        do i=1,m
          num=num+real( eps((i-1)*m+i),dp )
        end do     

     end if
     !
     !        adjust shift, find gap. In principle, this is a hack,
     !        because the energy surface changes, and we should
     !        restart the algorithm. In practice, we find that it
     !        hardly has any impact on the performance.
     !
     if (ishift .and. (rho < 0.1d0) .and. &                 ! adjust shift of
          iand(pw_params%optimize, 256) == 256) then            ! Hamiltonian
        xsi = eps  
        call matrix_polish(s, xsi, mu, nu, m)  
        mu = xsi  
        call mzheev('V', 'L', m, mu(1), m, eval(1), nu(1), m_squared, &
             rwork(1), diaginfo)
        if (diaginfo /= 0) then  
           write(9, *) 'diagonalization failed:', diaginfo  
           call printmatrix(mu, m, m)  
           call mystop  
        end if
        do i = 1, m
           iord(i) = i
        end do
        call sort(m, eval(1), iord(1))          ! sort with numerical recipes
        if (neigmin > m) then  
           write(9, *) 'diagonalize_hamiltonian: neigmin > m!'  
           write(9, *) neigmin, m  
           call mystop  
        else if (neigmin < m) then  
           !
           !               find largest gap, and throw other vectors away
           !
           demax = dzero 
           imax = neigmin + 1  

           if(m .eq. neigmax-1) then  ! necessary since eval(neigmax) 
              mmax=m-1                ! is not accurate
           else 
              mmax=m
           end if

           do i = neigmin + 1, m  
              de = (eval(i) - eval(i - 1)) / (real(i, dp)**2.5d0)  

              if (demax < de) then  
                 demax = de
                 imax = i - 1  
              end if
           end do
           mnew = imax  
           !
           !              superimpose x to form eigenvectors, and rearrange
           if (iand(pw_params%output(1), 1) == 1) t2 = gimmetime()  
           call mzgemm('N', 'N', len, m, m, zone, x(1), len, mu(1), &
                m, zzero, g(1, 1), len)
           if (iand (pw_params%output(1), 1) == 1) then  
              t1 = t1 + gimmetime() - t2  
              nops = nops + 8.0d0 * real(m * m * len, dp)
           end if
           call wavefn_rearrange(len, m, g(1, 1), x(1), iord(1), mnew)
           egap = eval(mnew + 1) - eval(mnew)  
           ham%shift = ham%shift + eval(mnew) + egap * dqtr + &
                pw_params%shiftsafety
           iconjug = .false.             ! don't reuse old conjugate vector
           ireset = .true.                   ! recompute H psi from scratch
           m = mnew  
           ishift = .false.  
           goto 100  
        else                               ! neigmin.eq.m, no rearrangement
           ham%shift = ham%shift + eval(m) + egap * dqtr + &
                pw_params%shiftsafety
           ireset = .true.                   ! recompute H psi from scratch
           ishift = .false.  
        end if
     end if
     !
     !        calculate the gradient
     !
     if (iand(pw_params%output(1), 8) == 8) call get_timing(time(4),t3)
     if (iand(pw_params%output(1), 1) == 1) t2 = gimmetime()  
     call mzgemm('N', 'N', len, m, m, zone, x(1), len, eps(1), m, &
          zzero, g(1, 1), len)
     call mzgemm('N', 'N', len, m, m, zone, y(1, 1), len, s(1), m, &
          zone, g(1, 1), len)
     if (iand(pw_params%output(1), 1) == 1) then  
        t1 = t1 + gimmetime() - t2  
        nops = nops + 1.6d1 * real(m * m * len)
     end if
      if (iand(pw_params%output(1), 8) == 8) call get_timing(time(7),t3)
     call mzaxpy(len * m, zmtwo, y(1, 1), 1, g(1, 1), 1)
     !
     !        make conjugate
     !
     !prepare for precondi
     call mzcopy(len * m, g(1, 1), 1, gp(1, 1), 1)  
     !
     !        change the preconditioning energy during the iteration.
     !        this is not good style, since it effectively changes the
     !        matrix under consideration, and we lose the property of
     !        conjugation. we tried it out, and it doesn't hurt performance.
     !        the factor of 1.35 was found by numerical experiment to give
     !        best performance. Smaller factors are dangerous, as they kill
     !        good components, and that hurts a lot.
     !
     eprecond = average_ekin(x(1), xsi(1), ham%gspace%ekin(1), m, len) * 1.35d0
     call precondition(gp(1, 1), eprecond, ham%ekinmod(1), &
          ham%gspace%ekin(1), m, len)
 
     rho = real(parallel_zdotc(len * m, g(1, 1), 1, gp(1, 1), 1), dp)
     !         if(.not.ireset) then
     if ((n > 1) .and. (.not.ireset)) then  
        !     step for Polak-Ribiere. Fletcher-Reeves would not require gold
        !        but is performing terribly in many cases.
        delrho = real(parallel_zdotc(len * m, gold(1, 1), 1, gp(1, 1), 1), dp)
        dcgamma = cmplx((rho - delrho) / oldrho, dzero, dp)  
        !   dcgamma = cmplx(rho/oldrho) ! -------- do this for fletcher-Reeves
        call mzaxpy(len * m, dcgamma, h(1, 1), 1, gp(1, 1), 1)  
     end if

     if (iand(pw_params%output(1), 536870912 ) == 536870912 ) &
              write(9,*) rho,n+nover_cnt,'   rho'
     if (n .eq. 1) rho_init=rho      

     call mzcopy(len * m, gp(1, 1), 1, h(1, 1), 1)  
     call mzcopy(len * m, g(1, 1), 1, gold(1, 1), 1)  ! save gradient into gold

     oldrho = rho  
     if (rho < delta ) then  
        if (n.eq.1) return
        iconv = .true.  
        goto 200  
     end if
     !
     !        compute some of the matrices
     !
     if (iand(pw_params%output(1), 1) == 1) t2 = gimmetime()  
     if (iand(pw_params%output(1), 8) == 8) call get_timing(time(4),t3)
     call hermit_overlap_matrix(eta(1), gp(1, 1), x(1), m, len)  
     call cheap_overlap_matrix(nu(1), gp(1, 1), m, len)  
     call hermit_overlap_matrix(mu(1), gp(1, 1), y(1, 1), m, len)  
     if (iand(pw_params%output(1), 1) == 1) then  
        t1 = t1 + gimmetime() - t2  
        nops = nops + 1.2d1 * real(m * m * len, dp)  
     end if
     if (iand(pw_params%output(1), 8) == 8) call get_timing(time(6),t3)
     call apply_ham_wavefn(time,pw_params%output(1), g(1, 1), ham, h(1, 1), &
          m, ews(1, 1), pw_params%nbandsfft)                        ! g=H * h
     num_HPsi=num_Hpsi+1
     if (iand(pw_params%output(1), 1) == 1) t2 = gimmetime() 
     if (iand(pw_params%output(1), 8) == 8  ) t3 = gimmetime()  
     call overlap_matrix(xsi(1), gp(1, 1), g(1, 1), m, len)  
     if (iand(pw_params%output(1), 1) == 1) then  
        t1 = t1 + gimmetime() - t2  
        nops = nops + 8.0d0 * real(m * m * len, dp)
     end if
     !
     !        calculate the coefficients for the minimization routine
     !
     if (iand(pw_params%output(1), 8) == 8) call get_timing(time(6),t3)
     a = dtwo * trace(mu(1), m) - contract(eta(1), eps(1), m) - &
          contract(mu(1), s(1), m)
     b = dtwo * trace(xsi(1), m) - contract(mu(1), eta(1), m) - &
          contract(eps(1), nu(1), m) - contract(xsi(1), s(1), m)
     c = -contract(xsi(1), eta(1), m) - contract(nu(1), mu(1), m)
     d = -contract(xsi(1), nu(1), m)  
     !
     !        do the line minimization
     !
     call findroot(lambda, dfour * d, dthree * c, dtwo * b, a)  

     if (iand(pw_params%output(1), 8) == 8) call get_timing(time(3),t3)
     if (lambda < dzero) then  
        ham%shift = ham%shift + dhalf  
        nreset = nreset + 1  
        write(9, *) '*** diagonalize_hamiltonian: resetting!'  
        ireset = .true.  
        if (nreset + n > maxit) goto 199        ! too many resets ... bye bye
        goto 100  
     end if
     dclambda = cmplx(lambda, dzero, dp)  
     !
     !        move to new position
     !
     call mzaxpy(len * m, dclambda, h(1, 1), 1, x(1), 1)  
     call mzaxpy(len * m, dclambda, g(1, 1), 1, y(1, 1), 1)
     if (iand(pw_params%output(1), 8) == 8) call get_timing(time(4),t3)
     !cjp write(9, 991) n, m, rho, ham%shift  
     if (n .eq. maxit .and. rho .gt. .3*rho_init  .and. nover_cnt .lt. pw_params%itdiag_add)then
         nover_cnt=nover_cnt+1
         go to 100
     end if
  end do                              ! ----------- end of big iteration loop

991  format('iter',I5,' m=',I5,' rho=',G14.8,' shift',F14.8)  
199  continue  

  iconv = .false.  

200  continue  

  !
  !  Done with the iterative part.
  !
  !  Now first orthogonalize the vectors. and rotate y also.
  !
  if (iand(pw_params%output(1), 8) == 8) call get_timing(time(4),t3)
  if (.not. iconv) &
    write(9,'(a,i3,a,g12.6)') ' after ',n+nover_cnt-1,' iter., residual error = ',rho
  call cheap_overlap_matrix(s,x,m,len);
  if (iand(pw_params%output(1), 8) == 8) call get_timing(time(6),t3)

  call occsp_sqrtinv_scal('U','n',myproc,s,m,ham%gspace%nproc)

  if (iand(pw_params%output(1), 8) == 8) call get_timing(time(3),t3)
  call mztrmm('R','U','N','N',len,m,zone,s(1),m,x(1),len)     
  call mztrmm('R','U','N','N',len,m,zone,s(1),m,y(1,1),len)
  if (iand(pw_params%output(1), 8) == 8) call get_timing(time(7),t3)
  !  
  !     the vectors x we have now span the right subspace, but are not yet
  !     eigenvectors. to get those, we diagonalize the epsilon matrix.
  !     the eigenvalues are the eigenvalues we are looking for.
  !     the eigenvectors tell us how to superimpose the vectors x to get
  !     the true eigenvalues.
  !
  call overlap_matrix(eps(1),x(1),y(1,1),m,len) ! eps = x^T*y

  if (iand(pw_params%output(1), 8) == 8) call get_timing(time(6),t3)

  call  occsp_diag_scal(myproc,eps,m,ham%gspace%nproc,eval,g,len*m,rwork)

  if (iand(pw_params%output(1), 8) == 8) call get_timing(time(3),t3)
  !
  !     combine the eigenvectors into g
  if (iand(pw_params%output(1), 1) == 1) t2 = gimmetime()  
  call mzgemm('N', 'N', len, m, m, zone, x(1), len, eps(1), m, &
       zzero, g(1, 1), len)
  if (iand (pw_params%output(1), 1) == 1) then  
     t1 = t1 + gimmetime() - t2  
     nops = nops + 8.0d0 * real(m * m * len, dp)
  end if
  !
  !     save the eigenvalues
  !
  do i = 1, m  
     eta(i) = cmplx(eval(i) + ham%shift, dzero, dp)  
  end do     
  neig = m                         ! return number of eigenvectors calculated
  !
  !     copy return eigenvectors into x array
  !
  call mzcopy(len * m, g(1, 1), 1, x(1), 1)                        !   x <- g
  !
  !     retrieve eigenvalues
  !
  eval(1:m) = real(eta(1:m), dp)  
  eval(m + 1) = eval(m) + egap                             ! and remember gap
  !
  !     performance printout
  !
  if (iand(pw_params%output(1), 1) == 1) write(9, 920) t1, nops / t1 * 1.0d-6
  !
  !     timing printout
  !
  if (iand (pw_params%output(1), 8) == 8) then  
     write(9, 910) n, gimmetime() - t0  
     call myflush(9)  
  end if
  if (iand(pw_params%output(1), 8) == 8) call get_timing(time(4),t3)
  !
  !     deallocate the work spaces
  !
  deallocate(eps)
  deallocate(s)
  deallocate(xsi) 
  deallocate(eta)
  deallocate(nu)
  deallocate(mu)
  deallocate(rwork)
  if (ham%vnloc%nvecs > 0 .and. .not. ham%pspot%NL_rspace(1) ) deallocate(ews)  
  deallocate(g)
  deallocate(gold)
  deallocate(h)
  deallocate(gp)
  deallocate(y)
  deallocate(iord)

  return
  
910  format(' DIAGONALIZATION TOOK ',i3,' ITERATIONS AND ', &
          &     f12.3, ' SECONDS')
920  format(' TIME FOR MATRIX-MATRIX:',f12.6,', MFLOPS=',f12.6)  
930  format('*** DIAGONALIZE_HAMILTONIAN: ALLOCATE(',a,'(',i6, &
          &        ',',i4,') FAILED. ERROR CODE=',i4)

end subroutine diagonalize_hamiltonian

