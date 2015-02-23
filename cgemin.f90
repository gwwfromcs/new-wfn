!
subroutine cgemin (time,t0,pw_params, maxit, iconv, bands, delta, kp, ham, &
     wavefn, gs, crys, energs, denc, rden, den, vion, vbareion)
  !
  !     1996 Bernd Pfrommer
  !
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
  type (crystal), intent(in) :: crys  
  integer, intent(in) :: maxit                      ! max number of iterations
  type(kpoint), intent(in) :: kp  
  type(hamiltonian), intent(inout) :: ham(kp%nrk,crys%nspin)      ! hamiltonia at kpoints
  type(pw_parameter), intent(in) :: pw_params             ! various parameters
  type(parallel_gspace), intent(in) :: gs  ! the gspace for the charge density

  complex(dp), intent(in) :: denc(gs%length), &      ! the core charge density
       vbareion(gs%length,crys%nspin)      ! the bare ionic potential in gspace
  real(dp), intent(in) :: delta                  ! accuracy of diagonalization
  !
  !     INPUT/OUTPUT:
  !     ------------
  !
  type(complex_gspace_array), intent(inout) :: wavefn  ! wavefunctions for all
                                                              ! kpoints, bands
  type(band), intent(inout) :: bands      ! eigenvalues and occupation numbers
                                                                 ! are touched
  type(energy), intent(inout) :: energs                 ! various energy terms
  complex(dp), intent(inout) :: vion (gs%length,crys%nspin)   ! vhxc on input,
                                               ! screened ionic pot. on output
  real(dp) time(num_time)
  !
  !     OUTPUT:
  !     -------
  !
  logical :: iconv                             ! true if converged, false else
  !
  !     WORK ARRAYS:
  !     -----------
  !
  complex(dp) :: rden(gs%r_size,crys%nspin), den(gs%length,crys%nspin)  
  !
  !     DESCRIPTION:
  !     -----------
  !
  !     Starting from given trial wavefunctions, this subroutine
  !     minimizes the energy functional with a conjugate gradient
  !     method. All bands and kpoints are moved simultaneously.
  !
  !
  !     --------------------- local variables ----------------------
  !
  real(dp) :: t0, t1, t2,t3, &                              ! timing variables
       oldrho, a, b, c, d, demax, de, lambda, eprecond, &
       delrho, a0, b0, c0, d0, &
       rho, &                                                  ! current error
       nops, &                                          ! number of operations
       egap                 ! gap between highest computed and next eigenvalue
  real(dp), pointer :: eval(:,:,:)  
  integer, parameter :: &
       maxrescg = 10             ! num of resets to switch to steepest descent
  integer :: myproc, &
       nreset, &                       ! counter for the number of resets done
       ierrmall, &
       nrk, &                                              ! number of kpoints
       mnew, &
       neig(crys%nspin), &                          ! number of eigenvalues to be computed
       diaginfo, &                   ! return value of diagonalization routine
       m(crys%nspin), &                                        ! current number of vectors
       m_max, &                                              ! max m
       i, j, &                                                       ! dummies
       irk, &                                                 ! kpoint counter
       imax, &                                             ! largest gap index
       n, &                                                  ! itertion number
       ldawfn, &                   ! leading dimension of wave function arrays
       length(kp%nrk),&                                       ! gspace length
       is,nspin
  logical :: ishift, ireset, iconjug, isteepest_desc  
  real(dp) denom,energy1,oldenergy,num,dweight,spinf
  logical dogcg
  integer line_min

  character(len=11) :: wfnname  
  complex(dp), pointer :: x(:,:,:)                  ! alias for wave functions
  complex(dp) :: dclambda, dclambdasq, dcgamma, dcweight, dclambda2  
  complex(dp) :: vin (gs%length,crys%nspin)              ! old potenital
  integer         num_HPsi,&                          ! op count of H*Psi 
         tstart              
  !
  !     work arrays
  !
  complex(dp), allocatable :: eps(:,:,:), s(:,:,:), xsi(:), &     ! of size m x m
       eta(:), nu(:), mu(:), epsdum(:)
  complex(dp), allocatable :: ews(:,:)                     ! of size nanl x m
  complex(dp), allocatable :: g(:,:,:), gold(:,:,:), &      ! of size len x m
       gp(:,:,:), h(:,:,:), y(:,:,:)
  integer, allocatable :: iord(:)                                 ! size of m
  real(dp), allocatable :: rwork(:)  
  !     local routines
  !
  real(dp), external :: trace, contract, average_ekin, gimmetime 
  !
  !     ------------------------------------------------------------------
  !
  num_HPsi=0
  tstart = gimmetime() 
  if (iand(pw_params%output(1), 8) == 8) t3 = tstart
  nrk = kp%nrk                        ! number of kpoints
  nspin=crys%nspin                                            
  if (bands%nspin >= 2) then
     spinf = dhalf
  else
     spinf = done 
  end if

  neig(1:nspin) = bands%min(1:nspin)
  eval => bands%energy         ! alias the energies to a more convenient name
  x => wavefn%data                             ! alias for the wave functions

  if (iand(pw_params%optimize, 32) /= 32) write(9, 140)  
  if (iand(pw_params%optimize, 128) == 128) then  
     isteepest_desc = .true.  
  else  
     isteepest_desc = .false.  
  end if
  nreset = 0  
  m(1:nspin) = neig(1:nspin)  
  
  do irk = 1, nrk  
     length(irk) = ham(irk,1)%gspace%length  
  end do
  ldawfn = size(x, 1)  
  myproc = ham(1,1)%gspace%myproc  
  !
  !     allocate all the work spaces
  !
  m_max=max(m(1),m(nspin))
  allocate(eps(m_max * m_max, nrk,nspin)) ; eps = zzero 
  allocate(s(m_max * m_max, nrk,nspin)) ; s = zzero
  allocate(xsi(m_max * m_max)) ; xsi = zzero
  allocate(eta(m_max * m_max)) ; eta = zzero  
  allocate(nu(m_max * m_max)) ; nu = zzero  
  allocate(mu(m_max * m_max)) ; mu = zzero  
  allocate(epsdum(m_max * m_max)) ; epsdum = zzero  
  allocate(rwork(max(1, 3 * m_max - 2))) ; rwork = dzero  
  if (ham(1,1)%vnloc%nvecs > 0 .and. .not. ham(1,1)%pspot%NL_rspace(1) ) then  
     allocate(ews(ham(1,1)%vnloc%nvecs, m_max))  
     ews = zzero
  end if
  allocate(g(ldawfn, nrk,nspin), stat = ierrmall) 
  if (ierrmall /= 0) then  
     write(9, 930) 'g', ldawfn, nrk, ierrmall  
     call mystop  
  end if
  g = zzero  
  allocate(gold(ldawfn, nrk,nspin), stat = ierrmall)
  if (ierrmall /= 0) then  
     write(9, 930) 'gold', ldawfn, nrk, ierrmall  
     call mystop  
  end if
  gold = zzero  
  allocate(h(ldawfn, nrk,nspin), stat = ierrmall)
  if (ierrmall /= 0) then  
     write(9, 930) 'h', ldawfn, nrk, ierrmall  
     call mystop  
  end if
  h = zzero  
  allocate(y(ldawfn, nrk,nspin), stat = ierrmall)
  if (ierrmall /= 0) then  
     write(9, 930) 'y', ldawfn, nrk, ierrmall  
     call mystop  
  end if
  y = zzero  
  allocate(gp(ldawfn, nrk,nspin), stat = ierrmall)
  if (ierrmall /= 0) then  
     write(9, 930) 'gp', ldawfn, nrk, ierrmall  
     call mystop  
  end if
  gp = zzero  

  allocate (iord(m_max)) ; iord = 0
  ireset = .false.                                   ! generally no resetting
  iconjug = .false.                      ! don't make conjugate in first step
  ishift = .true.               ! want to shift if error has decreased enough
  rho = 1.0d4                                       ! give some initial error
  nops = 0                                          ! reset operation counter
  t1 = dzero                                             ! reset time counter
  write(9, 100)  
  do irk = 1, nrk  
     write(9, 105) irk, eval(1:neig(1) + 1, irk, 1)  
  end do

  egap = minval(eval(neig(1) + 1,:, 1) - eval(neig(1),:, 1))  
  if (egap.lt.0) then                           ! no good gap as input, guess
     write(9, *) ' *** MINIMIZATION: INPUT EIGENVALUES BAD!'  
     call mystop  
  end if
  !
  !     shift hamiltonian just enough to make the hessian matrix
  !     negative definite   around the minimum
  !
  ham(1:nrk,1:nspin)%shift = maxval(eval(neig(1:nspin),:, 1)) + egap*dqtr + &
       pw_params%shiftsafety
  !
  !     --------------------- start iterative loop --------------------
  !
  do is=1,nspin
   do irk = 1, nrk  
    call cheap_overlap_matrix(s(1,irk,is), x(1,irk,is), m(is), length(irk)) 
    call occsp_sqrtinv_scal('U','n',myproc,s(1,irk,is),m(is),&
                                        ham(irk,is)%gspace%nproc)  
    call mztrmm('R','U','N','N',length(irk),m(is), zone,&
                 s(1,irk,is),m(is),x(1,irk,is),length(irk))
   end do
  end do

  oldenergy=dzero
  do n = 1, maxit  
911 continue                                         ! go here for resetting
    energy1=dzero
    if (iand(pw_params%output(1), 8) == 8) &
         write(9,*)  gimmetime()-t0,' CGEMIN TIME  BEFORE ITERATION',n 
    do is=1,nspin
     do irk = 1, nrk  
        ! y = H * x
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(4),t3) 
        call apply_ham_wavefn(time,pw_params%output(1), y(1, irk,is),&
              ham(irk,is), x(1, irk, is), m(is), ews(1,1), pw_params%nbandsfft)
        num_HPsi= num_HPsi + 1
        if (iand(pw_params%output(1), 1) == 1) t2 = gimmetime()  
        if (iand(pw_params%output(1), 8) == 8) t3 = gimmetime()
        ! s = x^T*x
        call cheap_overlap_matrix(s(1, irk,is), x(1, irk,is),m(is),length(irk))
        ! eps=x^T*y
        call overlap_matrix(eps(1,irk,is),x(1,irk,is),y(1,irk,is), &
                   m(is),length(irk))
        !
        !           calculate the gradient for each kpoint with proper weight
        !
        do i=1,m(is)
          energy1=energy1+real( eps((i-1)*m(is)+i,irk,is),dp )
        end do 
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(6),t3)
        call mzgemm('N', 'N', length(irk), m(is), m(is), zone, x(1, irk, is), &
             length(irk), eps(1,irk,is),m(is),zzero, g(1, irk,is), length(irk))
        call mzgemm('N', 'N', length(irk), m(is), m(is), zone, y(1, irk,is), &
             length(irk), s(1, irk,is), m(is), zone, g(1, irk,is), length(irk))
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(7),t3)
        dcweight = dmtwo
        call mzaxpy(length(irk) *m(is),dcweight,y(1,irk,is),1, g(1, irk,is), 1)
        if (iand(pw_params%output(1), 1) == 1) then  
           t1 = t1 + gimmetime() - t2  
           nops = nops + (12.d0 + 16.d0) * real(m(is) * ldawfn, dp)  
        end if
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(4),t3)
     end do
    end do
    !
    !        precondition and compute rho
    !
    rho = dzero
    do is=1,nspin
     do irk = 1, nrk  
        !prepare prec
        call mzcopy(length(irk) * m(is), g(1, irk,is), 1, gp(1, irk,is), 1)  

        if (iand(pw_params%optimize, 32) == 32) then  
           eprecond = 1.35d0 * average_ekin(x(1, irk, is), xsi(1), &
           ham(irk,is)%gspace%ekin(1), m(is), length(irk))
           call precondition(gp(1, irk,is), eprecond, ham(irk,is)%ekinmod(1), &
                ham(irk,is)%gspace%ekin(1), m(is), length(irk))
        end if

        dweight = kp%w(irk)  
        rho = rho + dweight*real(parallel_zdotc(length(irk)*m(is),g(1,irk,is),&
             1, gp(1, irk,is), 1), dp)
     end do
    end do

    write(9, '(1x,a,i3,a,g12.6)') 'iteration', n, &
          ' total energy estimated error [Ryd]=', rho
    call myflush(9)  
    irk = 1  
    !
    !        make conjugate
    !
    if ((n > 1) .and. (.not.ireset) .and. (.not.isteepest_desc)) then
       delrho = dzero
       do is=1,nspin
        do irk = 1, nrk
          dweight = kp%w(irk) 

          delrho = delrho + dweight*real(parallel_zdotc(length(irk)*m(is), & 
             gold(1, irk,is),  1, gp(1, irk, is),1), dp)
        end do
       end do

       dcgamma = cmplx((rho - delrho) / oldrho, dzero, dp)  
 
       do is=1,nspin   
        do irk = 1, nrk  
          call mzaxpy(length(irk)*m(is),dcgamma,h(1,irk,is), 1,gp(1,irk,is),1)
        end do
       end do    

    else  
        ireset = .false.  
    end if

    do is=1,nspin
     do irk = 1, nrk  
        call mzcopy(length(irk) * m(is), gp(1, irk,is), 1, h(1, irk,is), 1)  
        !           save gradient into old gradient
        call mzcopy(length(irk) * m(is), g(1, irk,is), 1, gold(1, irk,is), 1)
     end do
    end do

    oldrho = rho  
    if (rho < delta) then  
       iconv = .true.  
       goto 200  
    end if
    !
    !        compute the matrices and coefficients for the line minimization
    !
    if (iand(pw_params%output(1), 8) == 8) call get_timing(time(4),t3)
    a = dzero ; b = dzero ; c = dzero ; d = dzero  

    do is=1,nspin
     do irk = 1, nrk 
        if (iand(pw_params%output(1), 1) == 1) t2 = gimmetime()  
        call hermit_overlap_matrix(eta(1), h(1, irk,is), x(1, irk, is), &
             m(is), length(irk))
        call cheap_overlap_matrix(nu(1), h(1, irk,is), m(is), length(irk))  
        call hermit_overlap_matrix(mu(1), h(1, irk,is), y(1, irk,is), m(is),&
                                                            length(irk))
        if (iand (pw_params%output(1), 1) == 1) then  
           t1 = t1 + gimmetime() - t2
           nops = nops + 12.d0 * real(m(is) * m(is) * length(irk), dp)
        end if
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(6),t3)
        call apply_ham_wavefn(time,pw_params%output(1), g(1, irk,is),&
             ham(irk,is), h(1, irk,is), m(is), ews(1, 1), pw_params%nbandsfft)
        num_HPsi= num_HPsi + 1
        ! g=H
        if (iand(pw_params%output(1), 1) == 1) t2 = gimmetime()  
        if (iand(pw_params%output(1), 8) == 8) t3 = gimmetime()
        call overlap_matrix(xsi(1),h(1,irk,is),g(1,irk,is), m(is), length(irk))
        if (iand(pw_params%output(1), 1) == 1) then  
           t1 = t1 + gimmetime() - t2  
           nops = nops + 8.0d0 * real(m(is) * m(is) *length(irk), dp)  
        end if
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(6),t3)
        !
        !           calculate the coefficients for the minimization routine
        !
        a0 = kp%w(irk) * (dtwo * trace(mu(1), m(is)) - contract(eta(1), &
             eps(1, irk,is), m(is)) - contract(mu(1), s(1, irk,is), m(is)))
        b0 = kp%w(irk) * (dtwo * trace(xsi(1), m(is)) - contract(mu(1), &
             eta(1), m(is)) - contract(eps(1, irk,is), nu(1), m(is)) - &
             contract(xsi(1), s(1, irk,is), m(is)))
        c0 = kp%w(irk) * (-contract(xsi(1), eta(1), m(is)) - &
             contract(nu(1), mu(1), m(is)))
        d0 = kp%w(irk) * (-contract(xsi(1), nu(1), m(is)))  
        a = a + a0 ; b = b + b0 ; c = c + c0 ; d = d + d0  
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(3),t3)
     end do
    end do
    !
    !        do the line minimization
    !
    call findroot(lambda, dfour * d, dthree * c, dtwo * b, a)  
    if (lambda < dzero) then  
       nreset = nreset + 1  
       if (nreset > maxrescg) then  
          write(9, 700)  
          isteepest_desc = .true.  
       end if
       if (nreset + n > maxit) goto 199        ! too many resets ... bye bye
       ham(1:nrk,1:nspin)%shift = ham(1:nrk,1:nspin)%shift + dhalf  
       write(9, 800) ham(1,1)%shift  
       ireset = .true.  
       goto 911  
    end if
    dclambda = cmplx(lambda, dzero, dp)  
    !
    !        move wave functions (contained in x) to new position
    !                     x = x + lambda * h
    !
    do is=1,nspin
     do irk = 1, nrk  
        call mzaxpy(length(irk) * m(is), dclambda,h(1,irk,is),1,x(1,irk,is), 1)
     end do
    end do
    if (iand(pw_params%output(1), 8) == 8) call get_timing(time(4),t3)
    !
    !        compute new overlap matrix, and (3/2 -1/2*S) |x>
    !
    do is=1,nspin
     do irk = 1, nrk  
        if (iand(pw_params%output(1), 1) == 1) t2 = gimmetime()  
        call cheap_overlap_matrix(s(1, irk,is), x(1,irk,is),m(is),length(irk))
        ! s = x^T*x
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(6),t3)
        call mzgemm('N', 'N', length(irk), m(is), m(is), zmhalf,x(1,irk, is), &
             length(irk), s(1, irk,is),m(is),zzero, gp(1, irk,is), length(irk))
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(7),t3)
        call mzaxpy(length(irk) * m(is), ztrhalf,x(1,irk,is),1,gp(1,irk,is), 1)
        if (iand(pw_params%output(1), 1) == 1) then  
           t1 = t1 + gimmetime() - t2  
           nops = nops + 12.d0 * real(m(is) * m(is) * length(irk), dp)  
       end if       
         if (iand(pw_params%output(1), 8) == 8) call get_timing(time(4),t3)
     end do  
    end do
    !
    !        diagonalize to find eigenvalues
    !
    if (iand(pw_params%output(1), 8) == 8) call get_timing(time(4),t3)
    call update_hamiltonian(kp, ham(1,1), ldawfn, length(1), m, gp(1,1,1),gs, &
          crys, energs, pw_params, denc(1), rden(1,1), den(1,1), vion(1,1), &
          vbareion(1,1),vin,time,t3 )
    oldenergy=energs%total
 
    if (mod(n, pw_params%checkpoint_wfn) == 0) then  
       write(9, 130)  
       do is=1,nspin
         call writegsdat(0, gs, den(1,is), 1, 1, 1, 'CD', 2)  

         do irk = 1, kp%nrk  
           write(wfnname, 120) irk, 1  
           call writegsdat(0,ham(irk,is)%gspace,gp(1, irk,is),m(is),m(is), 1, &
                wfnname, 10)
         end do
       end do
       if (myproc == 0) then  
          call mysystem('rm -f BAND')  
          call write_band('BAND', bands)  
       end if
    end if
    if (iand(pw_params%output(1), 8) == 8) &
        write(9,*)  num_HPsi,' CGEMIN # H*Psi  AFTER ITERATION',n 
     ! ----------- end of big iteration loop
  end do
  
199 iconv = .false.  

200 continue  
  ! converged in first step, need to update vion
  if (n == 1) then  
    do is=1,nspin
       vion(:, is) = vion(:, is) + vbareion(:,is) 
       if (gs%ig0 > 0) vion(gs%ig0,is) = zzero
     end do 
  end if

  if (iand(pw_params%output(1), 8) == 8) call get_timing(time(4),t3)

  do is=1,nspin
   do irk = 1, nrk
     !
     !     Done with the iterative part.
     !
     !     Now first orthogonalize the vectors approximately. This step reduces
     !     linear errors in the eigenvalues, and is important.
     !
     if (iand(pw_params%output(1), 1) == 1) t2 = gimmetime()  
     call mzgemm('N', 'N', length(irk), m(is), m(is), zmhalf, x(1, irk, is), &
          length(irk), s(1, irk,is), m(is), zzero,gold(1, irk,is), length(irk))
     if (iand(pw_params%output(1), 1) == 1) then  
        t1 = t1 + gimmetime() - t2  
        nops = nops + 8.0d0 * real(m(is) * m(is) *length(irk))  
     end if
     if (iand(pw_params%output(1), 8) == 8) call get_timing(time(7),t3)
     call mzaxpy(length(irk) * m(is),ztrhalf,x(1,irk,is), 1, gold(1, irk,is),1)
     call mzcopy(length(irk) * m(is), gold(1, irk,is), 1, x(1, irk, is), 1)
     !        calculate transformed eps=1/4 (9eps-3St*eps-3eps*S+St*eps*S)
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(4),t3)
     call matrix_polish(s(1, irk,is), eps(1, irk,is), eta, nu, m)  
     !     the vectors x we have now span the right subspace, but are not yet
     !     eigenvectors. to get those, we diagonalize the epsilon matrix.
     !     the eigenvalues are the eigenvalues we are looking for.
     !     the eigenvectors tell us how to superimpose the vectors x to get
     !     the true eigenvalues
     !     There are parallel issues here: If the eps matrix is not
     !     exactly identical, the eigenvectors will acquire a different
     !     phase on different processors. Therefore, we diagonalize only
     !     on proc 0, and broadcast the result.
     !
     !
     if (myproc == 0) then 
        call mzheev('V', 'L', m(is), eps(1, irk,is), m(is),eval(1, irk, is), &
             g(1, irk,is), m(is) * length(irk), rwork(1), diaginfo)
        if (diaginfo /= 0) then  
           if (myproc == 0) write(9, *) &
                'diagonalization failed with info=', diaginfo
           call mystop  
        end if
     end if
     call my_broadcast(eps(1, irk,is), m(is) * m(is), 0)  
     call my_broadcast(eval(1, irk, is), m(is), 0)  
     if (iand(pw_params%output(1), 8) == 8) call get_timing(time(3),t3)
     !
     !        combine the eigenvectors into g
     if (iand(pw_params%output(1), 1) == 1) t2 = gimmetime()  
     call mzgemm('N', 'N', length(irk), m(is), m(is), zone, x(1, irk, is), &
          length(irk), eps(1, irk,is), m(is), zzero, g(1, irk,is), length(irk))
     if (iand(pw_params%output(1), 1) == 1) then  
        t1 = t1 + gimmetime() - t2  
        nops = nops + 8.0d0 * real(m(is) * m(is) * length(irk), dp)  
     end if
     if (iand(pw_params%output(1), 8) == 8) call get_timing(time(7),t3)
     !
     !        save the eigenvalues
     !
     eta(1:m(is)) = cmplx(eval(1:m(is), irk, 1)+ham(irk,is)%shift, dzero, dp)  
     !        in case the calculation did not converge, preform
     !        the orthonormalization of the  subspace:
     !
     !        a) Compute S^-1/2 = AT * D^-1/2 * A
     !        b) Compute x_new= S^-1/2 * x
     !
     !        The eigenvalue array is used as temporary storage
     !
     if (.not.iconv) then  
        write(9, '(a,i3,a,g12.6)') ' minimization not conv. after ', &
             n, ' iter., residual error = ', rho

        call cheap_overlap_matrix(s(1, irk,is), g(1, irk,is),m(is),length(irk))
        !           diagonalize S, find eigenvectors A (returned in S)
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(6),t3)
        call mzheev('V', 'L', m(is), s(1, irk,is), m(is), eval(1, irk, is), &
             x(1, irk, is), m(is) * length(irk), rwork(1), diaginfo)
        !        compute At * D^-1/2 * A
        do i = 1, m(is)  
           eval(i, irk, is) = done / sqrt(eval(i, irk, is))  
           do j = 1, m(is)  
              mu(i + m(is) * (j - 1) ) = eval(i, irk, is) * &
                   conjg(s(j + m(is) * (i - 1), irk,is))
           end do
        end do
        !        compute nu = A * D^-1/2* AT
        call mzgemm('N', 'N', m(is), m(is), m(is), zone, s(1, irk,is), m(is),&
                 mu(1), m(is), zzero, nu(1), m(is))
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(3),t3)
        !        compute g = nu * x
        call mzcopy(length(irk) * m(is), g(1, irk,is), 1, x(1, irk, is), 1)
        !   x <- g
        if (iand(pw_params%output(1), 1) == 1) t2 = gimmetime()  
        call mzgemm('N', 'N', length(irk), m(is), m(is), zone, x(1, irk, is), &
             length(irk), nu(1), m(is), zzero, g(1, irk,is), length(irk))
        if (iand(pw_params%output(1), 1) == 1) then  
           t1 = t1 + gimmetime() - t2  
           nops = nops + 8.0d0 * real(m(is) * m(is) * length(irk))  
        end if
     end if
     !
     !     copy return eigenvectors into x array
     !
     !   x <- g
     call mzcopy(length(irk) * m(is), g(1, irk,is), 1, x(1, irk, is), 1)  
     !
     !     retrieve eigenvalues
     !
     eval(1:m(is), irk, is) = real(eta(1:m(is)), dp)  
     ! and remember gap
     eval(m(is) + 1, irk, is) = eval(m(is), irk, is) + egap  
     ! end of loop over kpoints
     if (iand(pw_params%output(1), 8) == 8) call get_timing(time(7),t3)

   end do
  end do
  !
  !     performance printout
  !
  if (iand(pw_params%output(1), 1) == 1) write(9, 920) t1, nops / t1 * 1.0d-6
  !
  !     timing printout
  !
  if (iand(pw_params%output(1), 8) == 8) then  
     write(9, 910) n, gimmetime() - tstart  
     call myflush(9)  
  end if
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
  deallocate(epsdum)
  if (ham(1,1)%vnloc%nvecs>0 .and. .not. ham(1,1)%pspot%NL_rspace(1)) &
                                         deallocate(ews)
  deallocate(g)
  deallocate(gold)
  deallocate(h)
  deallocate(gp)
  deallocate(y)
  deallocate(iord)

  if (iand(pw_params%output(1), 8) == 8) call get_timing(time(4),t3)

  return

100 format(/' EIGENVALUES FROM SUBMATRIX DIAGONALIZATION:', &
       &     /1000(5f12.6/))
105 format(/' KPOINT:',i4,/1000(5f12.6/))  
120 format('WFN',i5.5,'.',i1)  
130 format(' <<<<<< checkpointing wave functions and charge >>>>>>')  

140 format(/' DIRECT ENERGY MINIMZATION WITHOUT PRECONDITIONING!', &
       &     /' THIS CAN BE A REAL PERFORMANCE HIT.')
700 format(' *** cgemin: switching over to steepest descend')  
800 format(' *** cgemin: resetting and shifting to: ', f12.6)  
910 format(/' MINIMIZATION TOOK ',i3,' ITERATIONS AND ', &
       &     f12.3, ' SECONDS')
920 format(' TIME FOR MATRIX-MATRIX:',f12.6,', MFLOPS=',f12.6)  
930 format('*** DIAGONALIZE_HAMILTONIAN: ALLOCATE(',a,'(',i6, &
       &        ',',i4,') FAILED. ERROR CODE=',i4)

end subroutine cgemin
