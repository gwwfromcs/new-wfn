!
subroutine cgemin_metal_gcg (time,t0,pw_params, maxit, iconv, bands, delta, kp, ham, &
     wavefn, gs, crys, energs, denc, rden, den, vion, vbareion)
  !
  !     2001 David Raczkowski
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
  integer, intent(in) :: maxit                      ! max number of iterations
  type(kpoint), intent(in) :: kp  
  type (crystal), intent(in) :: crys 
  type(hamiltonian), intent(inout) ::&
                             ham(kp%nrk,crys%nspin) ! hamiltonian at kpoints
  type(pw_parameter), intent(in) :: pw_params             ! various parameters
  type(parallel_gspace), intent(in) :: gs  ! the gspace for the charge density
  complex(dp), intent(in) :: denc(gs%length), &      ! the core charge density
       vbareion(gs%length,crys%nspin )     ! the bare ionic potential in gspace
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
  complex(dp), intent(inout) :: vion (gs%length, crys%nspin ) ! vhxc on input,
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
  real(dp) :: t0, t1, t2, t3,&                            ! timing variables
       oldrho,   lambda, eprecond, &
       delrho,  &
       rho, &                                                  ! current error
       nops, &                                          ! number of operations
       egap,&               ! gap between highest computed and next eigenvalue
       ee,eesqh,sq2h, FI,oldenergy_conv,dsmear
  real(dp), pointer :: eval(:,:,:)  
  integer :: myproc, &
       nreset, &                       ! counter for the number of resets done
       ierrmall, &
       nrk, &                                              ! number of kpoints
       mnew, &
       neig, &                          ! number of eigenvalues to be computed
       diaginfo, &                   ! return value of diagonalization routine
       m, &                                        ! current number of vectors
       m_squared, &                                              ! square of m
       i, j, &                                                       ! dummies
       irk, &                                                 ! kpoint counter
       imax, &                                             ! largest gap index
       n, &                                                  ! itertion number
       ldawfn, &                   ! leading dimension of wave function arrays
       length(kp%nrk)                                          ! gspace length
  logical :: ishift, ireset, iconjug, isteepest_desc, iconv1, F_conv  
  real(dp) denom,energy1,oldenergy,num,numF,dweight, tstart,lambdaF
  integer line_min

  character(len=11) :: wfnname  
  complex(dp), pointer :: x(:,:,:)                  ! alias for wave functions
  complex(dp) :: dclambda, dcgamma, dclambda2,dclambdaF,dclambdaF2
  complex(dp) :: vin (gs%length,crys%nspin)              ! old potenital
  integer         num_HPsi,&                          ! op count of H*Psi 
       ispn,is              
  !
  !     work arrays
  !
  complex(dp), allocatable :: eps(:,:,:), s(:,:,:), xsi(:), &  ! of size m x m
       eta(:), nu(:), mu(:), epsdum(:)
  complex(dp), allocatable :: ews(:,:)                     ! of size nanl x m
  complex(dp), allocatable :: g(:,:,:), gold(:,:,:), &    ! of size len x m
       gp(:,:,:), h(:,:,:), y(:,:,:)
  integer, allocatable :: iord(:)                                 ! size of m
  real(dp), allocatable :: rwork(:)  

  complex(dp), allocatable :: gc(:,:,:),gpc(:,:,:),gcold(:,:,:),&
            s_invhalf(:,:,:),F(:,:,:),del_F(:,:,:), grad_F(:,:,:)
  !     local routines
  !
  real(dp), external :: trace, contract, average_ekin, gimmetime
  integer outer_line,ipr
  logical do_update
  !
  !     ------------------------------------------------------------------
  !
  ee=exp(done)
  eesqh=sqrt(ee)*0.5
  sq2h=sqrt(dtwo)*0.5
  iconv1=.false.
  F_conv=.false.
  dsmear=pw_params%smearing/ryd
  rho=1d4
  do_update=.false.

  ispn=crys%nspin
  num_HPsi=0
  tstart = gimmetime()  
  if (iand(pw_params%output(1), 8) == 8) t3 = tstart
  nrk = kp%nrk                                            ! number of kpoints
  neig = bands%min(1)
  eval => bands%energy         ! alias the energies to a more convenient name
  x => wavefn%data                             ! alias for the wave functions

  if (iand(pw_params%optimize, 32) /= 32) write(9, 140)  
  if (iand(pw_params%optimize, 128) == 128) then  
     isteepest_desc = .true.  
  else  
     isteepest_desc = .false.  
  end if
  nreset = 0  
  m = neig  
  m_squared = m * m
  
  do irk = 1, nrk  
     length(irk) = ham(irk,1)%gspace%length  
  end do
  ldawfn = size(x, 1)  
  myproc = ham(1,1)%gspace%myproc  
  !
  !     allocate all the work spaces
  !
  allocate(eps(m * m, nrk,ispn)) ; eps = zzero 
  allocate(s(m * m, nrk,ispn )) ; s = zzero
  allocate(xsi(m * m)) ; xsi = zzero
  allocate(eta(m * m)) ; eta = zzero  
  allocate(nu(m * m)) ; nu = zzero  
  allocate(mu(m * m)) ; mu = zzero  
  allocate(epsdum(m * m)) ; epsdum = zzero  
  allocate(rwork(max(1, 3 * m - 2))) ; rwork = dzero  
  if (ham(1,1)%vnloc%nvecs > 0 .and. .not. ham(1,1)%pspot%NL_rspace(1)) then   
     allocate(ews(ham(1,1)%vnloc%nvecs, m))  
     ews = zzero
  end if
  allocate(g(ldawfn, nrk,ispn ), stat = ierrmall) 
  if (ierrmall /= 0) then  
     write(9, 930) 'g', ldawfn, nrk, ierrmall  
     call mystop  
  end if
  g = zzero  
  allocate(gold(ldawfn, nrk,ispn ), stat = ierrmall)
  if (ierrmall /= 0) then  
     write(9, 930) 'gold', ldawfn, nrk, ierrmall  
     call mystop  
  end if
  gold = zzero  
  allocate(h(ldawfn, nrk,ispn ), stat = ierrmall)
  if (ierrmall /= 0) then  
     write(9, 930) 'h', ldawfn, nrk, ierrmall  
     call mystop  
  end if
  h = zzero  
  allocate(y(ldawfn, nrk,ispn ), stat = ierrmall)
  if (ierrmall /= 0) then  
     write(9, 930) 'y', ldawfn, nrk, ierrmall  
     call mystop  
  end if
  y = zzero  
  allocate(gp(ldawfn, nrk,ispn ), stat = ierrmall)
  if (ierrmall /= 0) then  
     write(9, 930) 'gp', ldawfn, nrk, ierrmall  
     call mystop  
  end if
  gp = zzero  

  allocate(gc(ldawfn, nrk,ispn ), stat = ierrmall)
  if (ierrmall /= 0) then  
     write(9, 930) 'gc', ldawfn, nrk, ierrmall  
     call mystop  
  end if
  gc = zzero  

  allocate(gcold(ldawfn, nrk,ispn ), stat = ierrmall)
  if (ierrmall /= 0) then  
     write(9, 930) 'gcold', ldawfn, nrk, ierrmall  
     call mystop  
  end if
  gcold = zzero  

  allocate(gpc(ldawfn, nrk,ispn ), stat = ierrmall)
  if (ierrmall /= 0) then  
     write(9, 930) 'gpc', ldawfn, nrk, ierrmall  
     call mystop  
  end if
  gpc = zzero 

  allocate(s_invhalf(m*m, nrk,ispn ), stat = ierrmall)
  if (ierrmall /= 0) then  
     write(9, 930) 's_invhalf', m*m, nrk, ierrmall  
     call mystop  
  end if
  s_invhalf = zzero 

  allocate(F(m*m, nrk,ispn ), stat = ierrmall)
  if (ierrmall /= 0) then  
     write(9, 930) 'F', m*m, nrk, ierrmall  
     call mystop  
  end if
  F = zzero 

  allocate(grad_F(m*m, nrk,ispn ), stat = ierrmall)
  if (ierrmall /= 0) then  
     write(9, 930) 'grad_F', m*m, nrk, ierrmall  
     call mystop  
  end if
  grad_F = zzero 

  allocate(del_F(m*m, nrk,ispn ), stat = ierrmall)
  if (ierrmall /= 0) then  
     write(9, 930) 'del_F', m*m, nrk, ierrmall  
     call mystop  
  end if
  del_F = zzero

  allocate (iord(m)) ; iord = 0
  ireset = .false.                                   ! generally no resetting
  iconjug = .false.                      ! don't make conjugate in first step
  nops = 0                                          ! reset operation counter
  t1 = dzero                                             ! reset time counter
  write(9, 100)  
  !
  !     --------------------- start iterative loop --------------------
  !
  oldenergy_conv=dzero
  oldenergy=dzero
  !
  ! set initial occupation matrix using eigenvalues from submatrix diag.
  !
!  call flevel(1, 2, pw_params, energs, crys, bands, kp)  
!  do is=1,ispn
!   do irk = 1, nrk   
!    do i=1,m
!     F(i+(i-1)*m,irk,is)=cmplx(bands%occup(i,irk,is)/(kp%w(irk)*dtwo),dzero,dp)
!     end do
!    end do
!  end do 
   vin=vion
!        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(4),t3)
!      call update_hamiltonian_metal(2,kp, ham(1,1), ldawfn, length(1), m, &
!         x(1, 1,1), gs,bands, crys, energs, pw_params, denc(1), rden(1,1),&
!        den(1,1), vion(1,1),vbareion(1,1),1,vin(1,1),.true.,oldenergy,time,t3)
!      oldenergy=energs%total

  ham(1:nrk,1:crys%nspin)%shift =dzero
  do n = 1, maxit  
911 continue                                         ! go here for resetting

    energy1=dzero
    write(9,*)  gimmetime()-t0,' CGEMIN TIME  BEFORE ITERATION',n 

    if (.not. do_update) then
      call flevel(1, 2, pw_params, energs, crys, bands, kp)  
      F=zzero
      do is=1,ispn
        do irk = 1, nrk   
          do i=1,m
            F(i+(i-1)*m,irk,is)= &
                  cmplx(bands%occup(i,irk,is)/(kp%w(irk)*dtwo),dzero,dp)
          end do
        end do
      end do 
    end if


    do is=1,ispn
      do irk = 1, nrk  
        dweight = kp%w(irk)
        ! y = H * x
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(4),t3)
        if(n .eq. 1 .or. do_update) then
        call apply_ham_wavefn(time,pw_params%output(1), y(1, irk,is), &
              ham(irk,is), x(1, irk, is), m, ews(1, 1), pw_params%nbandsfft)
        num_HPsi= num_HPsi + 1
        if (iand(pw_params%output(1), 1) == 1) t2 = gimmetime()  
        if (iand(pw_params%output(1), 8) == 8) t3 = gimmetime()
        ! eps=x^T*y
        call overlap_matrix(eps(1, irk,is), x(1, irk, is), y(1, irk,is),m,&
                                   length(irk))
!        else
!        call mzaxpy(length(irk)*m, -dclambda, g(1,irk,is), 1, y(1,irk,is), 1)
        end if

        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(6),t3)
        !
        !           calculate the gradient for each kpoint with proper weight
        !
        if(iand(pw_params%output(1),1).eq.1) t2=gimmetime()
        call mzgemm('N','N',length(irk),m,m, zmone, x(1,irk,is),length(irk),&
             eps(1,irk,is), m,  zzero,g(1,irk,is),length(irk)) 
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(7),t3)
        call mzaxpy(length(irk)*m, zone, y(1,irk,is), 1, g(1,irk,is), 1)
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(4),t3)
      end do
    end do
    !
    !        precondition and compute rho
    !
    rho = dzero
    do is=1,ispn
      do irk = 1, nrk  
        !prepare prec
        !
        call mzcopy(length(irk) * m, g(1, irk,is), 1, gp(1, irk,is), 1)  
        if (iand(pw_params%optimize, 32) == 32) then  
          eprecond = 1.35d0 * average_ekin(x(1, irk, is), xsi(1), &
          ham(irk,is)%gspace%ekin(1), m, length(irk))
          call precondition(gp(1, irk,is), eprecond, ham(irk,is)%ekinmod(1), &
                ham(irk,is)%gspace%ekin(1), m, length(irk))
        end if
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(4),t3) 
        !
        ! PROJECT gradient back onto tangent plane
        !
        call overlap_matrix(eta(1),x(1,irk,is),gp(1,irk,is),m,length(irk))
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(6),t3)
        call mzgemm('N','N',length(irk),m,m, zmone,x(1,irk,is),length(irk),&
             eta(1), m, zone,gp(1,irk,is),length(irk))  
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(7),t3)
        !
        ! mult. occup. metric & gradient gives gpc: used in innner products
        !
        call mzgemm('N','N',length(irk),m,m, zone, gp(1,irk,is),length(irk),&
             F(1,irk,is), m, zzero,gpc(1,irk,is),length(irk)) 
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(7),t3)
        call precondition_inv(gpc(1,irk,is),eprecond,ham(irk,is)%ekinmod(1), &
             ham(irk,is)%gspace%ekin(1), m,length(irk))

        dweight = kp%w(irk)  
        rho =rho+ dweight*real((parallel_zdotc(length(irk)*m,gp(1,irk,is),&
               1,gpc(1,irk,is),1)),dp)  
      end do
    end do
    write(9, '(1x,a,i3,a,g12.6)') 'iteration', n, &
          ' total energy estimated error [Ryd]=', rho
    call myflush(9)  
    !
    !        make conjugate
    !
!    if (n <4 ) ireset=.true.   ! steepest descent for first 4 steps
    if ((n > 1) .and. (.not.ireset) .and. (.not.isteepest_desc)) then
      delrho = dzero
      do is=1,ispn
        do irk = 1, nrk
          dweight = kp%w(irk)  
          delrho = delrho +  dweight*real(parallel_zdotc(length(irk)*m, & 
                gold(1, irk,is), 1, gpc(1, irk,is), 1), dp)
          call mzcopy(length(irk) * m, gp(1, irk,is), 1, gold(1, irk,is), 1) 
        end do
      end do 
      dcgamma = cmplx((rho - delrho) / oldrho, dzero, dp)  
    
      do is=1,ispn
        do irk = 1, nrk  
          call mzaxpy(length(irk)*m,dcgamma,h(1, irk,is), 1, gp(1, irk,is), 1)
        end do
      end do

    else  
      write(9,*) 'steepest descent'
       ireset = .false.  
    end if
    !
    !  copy gp to h - search direction 
    do is=1,ispn
      do irk = 1, nrk  
        call mzcopy(length(irk) * m, gp(1, irk,is), 1, h(1, irk,is), 1)  
      end do
    end do
    oldrho = rho
    !
    !  check for convergence
    if (rho < delta) then  
      iconv = .true.  
      goto 200  
    end if
   if (rho .lt. pw_params%init_dir_econv) do_update =.true.

    if (iand(pw_params%output(1), 8) == 8) call get_timing(time(4),t3)
    !
    ! do line minimization
    num=dzero; denom=dzero;
    do is=1,ispn
      do irk = 1, nrk 
        dweight = kp%w(irk) 

        call apply_ham_wavefn(time,pw_params%output(1), g(1,irk,is),&
            ham(irk,is), h(1,irk,is),m,ews(1,1),pw_params%nbandsfft)  ! g=H * h
        num_HPsi= num_HPsi + 1
        if (iand(pw_params%output(1), 8) == 8) t3 = gimmetime()
        !
        ! get (h.gpc) numerator in LM equation
        num = num+dweight*real(parallel_zdotc(length(irk)*m,h(1,irk,is),1,&
                gpc(1,irk,is),1),dp)
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(4),t3)
        !
        !  get denominator in LM equation
        call mzgemm('N','N',length(irk),m,m, zmone, h(1,irk,is),length(irk),&
             eps(1,irk,is), m, zzero,gp(1,irk,is),length(irk))          
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(7),t3)

        call mzaxpy(length(irk)*m, zone, g(1,irk,is), 1, gp(1,irk,is), 1)
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(4),t3)
        call mzgemm('N','N',length(irk),m,m, zone, gp(1,irk,is),length(irk),&
           F(1,irk,is),  m,  zzero,gpc(1,irk,is),length(irk)) 
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(7),t3)
        denom=denom + dweight*real(parallel_zdotc(length(irk)*m,&
                h(1,irk,is),1,gpc(1,irk,is),1),dp)
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(4),t3)
      end do
    end do

    lambda=abs(num)/abs(denom)
    dclambda=cmplx(lambda,dzero,dp)      
    !
    !        move wave functions (contained in x) to new position
    !                     x = x + lambda * h
    !        orthogonzalize, and get wavefunctions in same space as F
    gpc=x
    do is=1,ispn
      do irk = 1, nrk

        call mzaxpy(length(irk)*m, -dclambda , h(1,irk,is), 1, x(1,irk,is), 1)
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(4),t3)
        call cheap_overlap_matrix(s(1,irk,is),x(1,irk,is),m,length(irk));
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(6),t3)
!         diagonalize S, find eigenvectors A (returned in S) 
        call  occsp_diag_scal(myproc,s(1,irk,is),m,ham(irk,is)%gspace%nproc,&
            eval(1,irk,is),gp(1,irk,is),length(irk)*m,rwork)

!        compute At * D^-1/2 * A  
        do i=1,m
          eval(i,irk,is)=1.0/sqrt(eval(i,irk,is))
          do j=1,m
            mu(i+m*(j-1))=eval(i,irk,is)*conjg(s(j+m*(i-1),irk,is))
          end do
        end do
!        compute nu = A * D^-1/2* AT 
        call mzgemm('N','N', m,m,m, zone, s(1,irk,is),m,mu(1),m,zzero, &
             nu(1),m) 
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(3),t3)

!        compute g = nu * x  
        if(iand(pw_params%output(1),1).eq.1) t2=gimmetime()
        call mzgemm('N','N',length(irk),m,m, zone,x(1,irk,is),length(irk),&
             nu(1),m,zzero, gc(1,irk,is),length(irk)) 
          if(iand(pw_params%output(1),1).eq.1) then
            t1=t1+gimmetime()-t2
            nops = nops + 8*dble(m*m)*dble(length(irk))
          endif

        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(7),t3)
        call mzcopy(m*m,nu(1),1,s_invhalf(1,irk,is),1)
        call mzaxpy(length(irk)*m,-dclambda, g(1,irk,is), 1, y(1,irk,is), 1)

        if (.not. do_update) then
          call mzcopy(length(irk)*m,gc(1,irk,is),1,x(1,irk,is),1)
          call mzgemm('N','N',length(irk),m,m, zone,y(1,irk,is),length(irk),&
             s_invhalf(1,irk,is),m,zzero, gc(1,irk,is),length(irk)) 
          call mzcopy(length(irk)*m,gc(1,irk,is),1,y(1,irk,is),1)
        call overlap_matrix(eps(1, irk,is), x(1, irk, is), y(1, irk,is),&
                                                             m,length(irk))
        call  occsp_diag_scal(myproc,eps(1,irk,is),m,ham(irk,is)%gspace%nproc,&
            eval(1,irk,is),gp(1,irk,is),length(irk)*m,rwork)
   
       
          call mzgemm('N','N',length(irk),m,m, zone,x(1,irk,is),length(irk),&
             eps(1,irk,is),m,zzero, gc(1,irk,is),length(irk)) 
          call mzcopy(length(irk)*m,gc(1,irk,is),1,x(1,irk,is),1)
          call mzgemm('N','N',length(irk),m,m, zone,y(1,irk,is),length(irk),&
             eps(1,irk,is),m,zzero, gc(1,irk,is),length(irk))
          call mzcopy(length(irk)*m,gc(1,irk,is),1,y(1,irk,is),1)

        call mzgemm('N','N',length(irk),m,m, zone,h(1,irk,is),length(irk),&
             s_invhalf(1,irk,is),m,zzero, gc(1,irk,is),length(irk))
        call mzgemm('N','N',length(irk),m,m, zone,gc(1,irk,is),length(irk),&
             eps(1,irk,is),m,zzero, h(1,irk,is),length(irk))

        call mzgemm('N','N',length(irk),m,m, zone,gold(1,irk,is),length(irk),&
             s_invhalf(1,irk,is),m,zzero, gc(1,irk,is),length(irk))
       call mzgemm('N','N',length(irk),m,m, zone,gc(1,irk,is),length(irk),&
             eps(1,irk,is),m,zzero, gold(1,irk,is),length(irk))

        eps(1:m*m,irk,is)=zzero
        do i=1,m
          eps(i+(i-1)*m,irk,is) = cmplx(eval(i, irk, is),dzero,dp)
        end do

        else

        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(4),t3)
        call overlap_matrix(eps(1, irk,is), x(1, irk, is), y(1, irk,is),&
                                                             m,length(irk))
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(6),t3)

        call mzgemm('C','N',m,m,m, zone,s_invhalf(1,irk,is),m, eps(1,irk,is),&
                                                      m,zzero, xsi(1),m) 
        call mzgemm('N','N',m,m,m, zone,xsi(1),m, s_invhalf(1,irk,is),m,zzero,&
             epsdum(1),m) 

        call mzcopy(m*m,F(1,irk,is),1,mu(1),1)
        call  occsp_diag_scal(myproc,mu,m,ham(irk,is)%gspace%nproc,&
            eval(1,irk,is),gp(1,irk,is),length(irk)*m,rwork)

        bands%occup(1:m,irk,is)=eval(1:m,irk,is)*kp%w(irk)*dtwo 

        call mzgemm('C','N',m,m,m, zone,mu(1),m, epsdum(1),m,zzero, &
             xsi(1),m) 
        call mzgemm('N','N',m,m,m, zone,xsi(1),m, mu(1),m,zzero, nu(1),m) 

        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(3),t3)

        call mzgemm('N','N',length(irk),m,m, zone, gc(1,irk,is),length(irk),&
            mu(1),m, zzero, gp(1,irk,is), length(irk)) 
       
        do i=1,m
          eval(i, irk, is) = nu(i+(i-1)*m)
        end do
        eval(1:m, irk, is) = eval(1:m, irk, is) + ham(irk,is)%shift 
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(7),t3)
        end if

      end do
    end do

    dclambda2=dclambda

    if (F_conv) then
     ipr=2
    else
     ipr=0
    end if

    if (do_update) then
    call update_hamiltonian_metal(ipr,kp, ham(1,1), ldawfn, length(1), m,&
         gp(1, 1,1), gs,bands, crys, energs, pw_params, denc(1), rden(1,1),&
         den(1,1),vion(1,1), vbareion(1,1),n,vin(1,1),.true.,oldenergy,time,t3)
    else
     goto 25
    end if

     if (energs%total .gt. oldenergy) then
       dclambda2=(abs(num)/crys%nspin)/(.5*( (energs%total-oldenergy)/lambda +2*2*abs(num)/crys%nspin)/lambda)
            do is=1,ispn
      do irk = 1, nrk

  call mzaxpy(length(irk)*m, dclambda-dclambda2, h(1,irk,1), 1, x(1,irk,1), 1)   
     call mzaxpy(length(irk)*m,dclambda-dclambda2, g(1,irk,is), 1, &
           y(1,irk,is), 1)
        call overlap_matrix(eps(1, irk,is), x(1, irk, is), y(1, irk,is),m,length(irk))
        call cheap_overlap_matrix(s(1,irk,is),x(1,irk,is),m,length(irk));

!         diagonalize S, find eigenvectors A (returned in S) 
        if(myproc.eq.0) then
         call mzheev('V','L',m, s(1,irk,is), m, eval(1,irk,is),gp(1,irk,is),&
            m*length(irk),  rwork(1), diaginfo);  
        end if
        call my_broadcast(s(1,irk,is),m*m,0)
        call my_broadcast(eval(1,irk,is),m,0)

!        compute At * D^-1/2 * A  
          do i=1,m
            eval(i,irk,is)=1.0/sqrt(eval(i,irk,is))
            do j=1,m
               mu(i+m*(j-1))=eval(i,irk,is)*conjg(s(j+m*(i-1),irk,is))
           end do
          end do

!        compute nu = A * D^-1/2* AT 
          call mzgemm('N','N', m,m,m, zone, s(1,irk,is),m,mu(1),m,zzero, &
             nu(1),m) 

!        compute g = nu * x  
          if(iand(pw_params%output(1),1).eq.1) t2=gimmetime()
          call mzgemm('N','N',length(irk),m,m, zone,x(1,irk,is),length(irk),&
             nu(1),m,zzero, gc(1,irk,is),length(irk)) 
          if(iand(pw_params%output(1),1).eq.1) then
            t1=t1+gimmetime()-t2
            nops = nops + 8*dble(m*m)*dble(length(irk))
          endif

        call mzcopy(m*m,nu(1),1,s_invhalf(1,irk,is),1)

         call mzgemm('C','N',m,m,m, zone,s_invhalf(1,irk,is),m, eps(1,irk,is),m,zzero, &
             xsi(1),m) 

        call mzgemm('N','N',m,m,m, zone,xsi(1),m, s_invhalf(1,irk,is),m,zzero, &
             epsdum(1),m) 

       mu(1:m*m)=F(1:m*m,irk,is)
        if(myproc.eq.0) then
          call mzheev('V','L',m,mu(1),m, eval(1,irk,is),gp(1,irk,is), &
               m*length(irk), rwork(1),diaginfo) 
          if(diaginfo.ne.0) then
            if(myproc.eq.0)&
                write(9,*) 'diagonalization failed with info='&
                ,diaginfo
            call mystop
          endif
       endif
        call my_broadcast(mu(1),m*m,0)
        call my_broadcast(eval(1,irk,is),m,0)

       bands%occup(1:m,irk,is)=eval(1:m,irk,is)*kp%w(irk)*dtwo 

        call mzgemm('C','N',m,m,m, zone,mu(1),m, epsdum(1),m,zzero, &
             xsi(1),m) 
        call mzgemm('N','N',m,m,m, zone,xsi(1),m, mu(1),m,zzero, &
             nu(1),m) 

        call mzgemm('N','N',length(irk),m,m, zone, gc(1,irk,is),length(irk),&
            mu(1),m, zzero, gp(1,irk,is), length(irk)) 
       
       do i=1,m

         eval(i, irk, is) = nu(i+(i-1)*m)
       end do

       eval(1:m, irk, is) = eval(1:m, irk, is) + ham(irk,is)%shift 

      end do
    end do

!      else 
!        dclambda2=dclambda
!      end if
         write(9,*) dclambda,dclambda2,num,denom,' lambdaX'
         write(9,*) 'lambdaX', energs%total,oldenergy 

    call update_hamiltonian_metal(ipr,kp, ham(1,1), ldawfn, length(1), m,&
         gp(1, 1,1), gs,bands, crys, energs, pw_params, denc(1), rden(1,1),&
         den(1,1),vion(1,1), vbareion(1,1),n,vin(1,1),.true.,oldenergy,time,t3)
 
    end if


    oldenergy=energs%total
    write(9,*) energs%total
    x=gc
    if (F_conv) go to 15
    !
    ! do 2 line minization w.r.t F
    !
!    if (rho/m .lt. 2.5d-3) then
      line_min=2
!    else
!      line_min=1
!    end if      

    do outer_line=1,line_min

      numF=dzero; denom=dzero;  
      !
      ! calculate dA/df = eps-TdS/df
      !  
      grad_F=zzero;
      do is=1,ispn
        do irk = 1, nrk  
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(4),t3)
          call apply_ham_wavefn(time,pw_params%output(1), y(1, irk,is), &
           ham(irk,is), x(1, irk, is), m, ews(1, 1), pw_params%nbandsfft)
          num_HPsi= num_HPsi + 1
          if (iand(pw_params%output(1), 8) == 8) t3 = gimmetime()

          call overlap_matrix(eps(1, irk,is), x(1, irk, is), y(1, irk,is),&
                     m,length(irk))
          if (iand(pw_params%output(1), 8) == 8) call get_timing(time(6),t3)

          call mzcopy(m*m,F(1,irk,is),1,mu(1),1)
          call  occsp_diag_scal(myproc,mu,m,ham(irk,is)%gspace%nproc,&
            eval(1,irk,is),gp(1,irk,is),length(irk)*m,rwork)

          !
          ! get dS/df of smearing method in diag. rep.
          do i=1,m
            FI=real(eval(i,irk,is),dp)
            if (pw_params%smearing_method .eq. 2) then
              if(abs(FI) .gt. 1.d-06 .and. abs(FI-done) .gt. 1.d-06) then
                grad_F(i+(i-1)*m,irk,is)= log(FI) - log(done-FI) 
              end if
            else
              if(abs(FI) .gt. 1.d-6 .and. abs(FI-done) .gt. 1.d-6) then
                if (FI .lt. 0.5d0) then
                  grad_F(i+(i-1)*m,irk,is)=(sqrt(log(eesqh/FI))-sq2h)
                else
                  grad_F(i+(i-1)*m,irk,is)=&
                      (sq2h-sqrt(log(eesqh/(done-FI))))
                end if
              end if
            end if

          end do        
          !
          !  transform bacl to original rep.
          call mzgemm('N','N',m,m,m, zone,mu(1),m, grad_F(1,irk,is),m,&
             zzero, xsi(1),m) 
          call mzgemm('N','C',m,m,m, zone,xsi(1),m, mu(1),m,zzero, &
             grad_F(1,irk,is),m) 
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(3),t3)
          ! put in smearing(T) factor
          call mzdscal(m*m, dsmear,grad_F(1,irk,is),1)
          ! now make grad_F = dA/df
          call mzaxpy(m*m,zone,eps(1,irk,is),1,grad_F(1,irk,is),1)
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(4),t3)
          !
          ! diag eps: used later
          call mzcopy(m*m,eps(1,irk,is),1,s(1,irk,is),1)
          call  occsp_diag_scal(myproc,s(1,irk,is),m,ham(irk,is)%gspace%nproc,&
            eval(1,irk,is),gp(1,irk,is),length(irk)*m,rwork)

          eval(1:m, irk, is) = eval(1:m, irk, is) + ham(irk,is)%shift
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(3),t3)
        end do
      end do
      !     
      ! calculate F(1) in a diag. rep.
      !
      if (outer_line .eq. line_min) then
       ipr=2
      else
       ipr=0
      end if
      call flevel(1, ipr, pw_params, energs, crys, bands, kp) 

      lambdaF=done
      dclambdaF=zone   !cmplx(lambdaF,dzero,dp);
      del_F=zzero; 

      do is=1,ispn
       do irk = 1, nrk
         !     
         ! calculate  F(1)-F(0) which is the search direction
         !  
         do i=1,m
          del_F(i+(i-1)*m,irk,is)=cmplx(bands%occup(i,irk,is)/ &
            (kp%w(irk)*dtwo),dzero,dp) 
         end do
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(4),t3)
         ! rotate F(1) into same rep as F(0)
         call mzgemm('N','N',m,m,m, zone,s(1,irk,is),m, del_F(1,irk,is),m,&
             zzero, xsi(1),m) 
         call mzgemm('N','C',m,m,m, zone,xsi(1),m, s(1,irk,is),m,zzero, &
             del_F(1,irk,is),m) 
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(3),t3)
         ! del_f=F(1)-F(0)
         call mzaxpy(m*m,zmone,F(1,irk,is),1,del_F(1,irk,is),1) 
         ! necessary innner products for line minization formula
         denom=denom+kp%w(irk)*mzdotc(m*m,del_F(1,irk,is),1,del_F(1,irk,is),1)
         numF=numF+kp%w(irk)*mzdotc(m*m,grad_F(1,irk,is),1,del_F(1,irk,is),1)
         ! move F
         call mzaxpy(m*m,dclambdaF,del_F(1,irk,is),1,F(1,irk,is),1)
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(4),t3)
         !
         ! diagonalize F
         call mzcopy(m*m,F(1,irk,is),1,mu(1),1)
         call  occsp_diag_scal(myproc,mu,m,ham(irk,is)%gspace%nproc,&
            eval(1,irk,is),gp(1,irk,is),length(irk)*m,rwork)
         bands%occup(1:m,irk,is)=eval(1:m,irk,is) *kp%w(irk)*dtwo 
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(3),t3)
         !
         ! align eps and x to be in same rep. as F
         call mzgemm('N','N',length(irk),m,m, zone,x(1,irk,is),length(irk),&
             mu(1),m,zzero, gp(1,irk,is),length(irk)) 
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(7),t3)
         call mzgemm('C','N',m,m,m, zone,mu(1),m, eps(1,irk,is),m,zzero, &
             xsi(1),m) 
         call mzgemm('N','N',m,m,m, zone,xsi(1),m, mu(1),m,zzero, &
             nu(1),m) 
         ! get trace = eband
         do i=1,m
          eval(i, irk, is) = nu(i+(i-1)*m)
         end do
         eval(1:m, irk, is) = eval(1:m, irk, is) + ham(irk,is)%shift
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(3),t3)
       end do
      end do
      !
      !  calc. A(1)
      !
      call update_hamiltonian_metal(0,kp, ham(1,1), ldawfn, length(1), m,&
         gp(1, 1,1), gs,  bands, crys, energs, pw_params, denc(1),&
          rden(1,1), den(1,1), vion(1,1), &
          vbareion(1,1),n,vin(1,1),.false.,oldenergy,time,t3)
      !
      ! calc. step size in F(1)-F(0) direction.
      !
      dclambdaF2=cmplx((abs(numF)/crys%nspin)/& 
        ((energs%total-oldenergy)/lambdaF+2*abs(numF)/crys%nspin)/lambdaF &
                                        ,dzero,dp)

      if ( real(dclambdaF2,dp) .gt. done ) dclambdaF2=dclambdaF
      if( real(dclambdaF2,dp) .lt. dzero ) then
        if( oldenergy-energs%total .gt. dzero ) then
          dclambdaF2=dclambdaF
        else 
         dclambdaF2=zzero
        end if
      end if
      write(9,*) numF,denom,dclambdaF2,' lambdaF',oldenergy-energs%total
      !
      ! obtain F,eps, and x in proper rep. for calc. of A 
      !
      do is=1,ispn
       do irk = 1, nrk

         call mzaxpy(m*m,dclambdaF2-dclambdaF,del_F(1,irk,is),1,F(1,irk,is),1)
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(4),t3)
         call mzcopy(m*m,F(1,irk,is),1,mu(1),1)
         call  occsp_diag_scal(myproc,mu,m,ham(irk,is)%gspace%nproc,&
            eval(1,irk,is),gp(1,irk,is),length(irk)*m,rwork)
         bands%occup(1:m,irk,is)=eval(1:m,irk,is)*kp%w(irk)*dtwo  
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(3),t3)
         !
         ! align eps and x to be in same rep. as F
         call mzgemm('N','N',length(irk),m,m, zone,x(1,irk,is),length(irk),&
             mu(1),m,zzero, gp(1,irk,is),length(irk)) 
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(7),t3)
         call mzgemm('C','N',m,m,m, zone,mu(1),m, eps(1,irk,is),m,zzero, &
             xsi(1),m) 
         call mzgemm('N','N',m,m,m, zone,xsi(1),m, mu(1),m,zzero, &
             nu(1),m) 
         ! get trace = eband        
         do i=1,m
           eval(i, irk, is) = nu(i+(i-1)*m)
         end do
         eval(1:m, irk, is) = eval(1:m, irk, is) + ham(irk,is)%shift
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(3),t3)

       end do
      end do

      if (outer_line .eq. line_min) then
       ipr=2
      else
       ipr=0
      end if
      ! calc. A
      call update_hamiltonian_metal(ipr,kp, ham(1,1), ldawfn, length(1), m,&
          gp(1, 1,1), gs, bands, crys, energs, pw_params, denc(1), rden(1,1),&
          den(1,1),vion(1,1),vbareion(1,1),n,vin(1,1),.true.,oldenergy,time,t3)

!      write(9,*) 2*numf,(energs%total-oldenergy)/0.01,' deriv_test'

!      call mystop
!      if (n .gt. 1 .and. oldenergy .lt. energs%total) then
!        do is=1,ispn
!          do irk = 1, nrk
!            call mzaxpy(m*m,-dclambdaF2,del_F(1,irk,is),1,F(1,irk,is),1)
!          end do
!        end do
!        write(9,*) oldenergy,energs%total,' bad energy'
!        goto 15
!      end if
      oldenergy=energs%total 

     if (denom/m .lt. 1d-16) then
        F_conv = .true. 
        go to 15
      end if

    end do  ! end line minization w.r.t. F

15  continue

    if(n .gt. 1) then
      if(abs(oldenergy_conv-energs%total) .lt. pw_params%epsce) then
        if(iconv1) then
          goto 200
        else
          iconv1 = .true.
        end if
      else
        iconv1 = .false.
      end if
    end if            

    oldenergy_conv=energs%total

    do is=1,ispn
!       do irk = 1, nrk
!! get gc to update h
!         call cheap_overlap_matrix(nu(1),h(1,irk,is),m,length(irk))
!         call mzgemm('N','N',length(irk),m,m, zone, gpc(1,irk,is),&
!                length(irk),nu(1), m, zzero,gc(1,irk,is),length(irk))
!! update gold
!         call overlap_matrix(eta(1),h(1,irk,is),gold(1,irk,is),m,length(irk))
!         call mzgemm('N','N',length(irk),m,m, zone*dclambda2, gpc(1,irk,is),&
!              length(irk), eta(1), m, zone,gold(1,irk,is),length(irk))
!         call mzaxpy(length(irk)*m, dclambda2, gc(1,irk,is), 1, h(1,irk,is),1)
!       end do
 
      if (mod(n, pw_params%checkpoint_wfn) == 0) then  
        write(9, 130)  
        call writegsdat(0, gs, den(1,1), 1, 1, 1, 'CD', 2)  
        do irk = 1, kp%nrk  
           write(wfnname, 120) irk, 1  
           call writegsdat(0, ham (irk,is)%gspace, gp(1, irk,is), m, m, 1, &
                wfnname, 10)
        end do
        if (myproc == 0) then  
           call mysystem('rm -f BAND')  
           call write_band('BAND', bands)  
        end if
      end if
      if (iand(pw_params%output(1), 8) == 8) call get_timing(time(4),t3)

    end do

    write(9,*)  num_HPsi,' CGEMIN # H*Psi  AFTER ITERATION',n

    do is=1,ispn
      do irk = 1, nrk
        call mzgemm('N','N',length(irk),m,m, zone,h(1,irk,is),length(irk),&
             s_invhalf(1,irk,is),m,zzero, gc(1,irk,is),length(irk))
      end do
    end do
    h=gc

    do is=1,ispn
      do irk = 1, nrk
        call mzgemm('N','N',length(irk),m,m, zone,gold(1,irk,is),length(irk),&
             s_invhalf(1,irk,is),m,zzero, gc(1,irk,is),length(irk))
      end do
    end do
    gold=gc

    if (iand(pw_params%output(1), 8) == 8) call get_timing(time(7),t3)

25 continue
  end do      ! ----------- end of big iteration loop
  
199 iconv = .false.  

200 continue  
  ! converged in first step, need to update vion
  if (n == 1) then  
     do is=1,ispn
       vion(:, is) = vion(:, is) + vbareion(:,is) 
       if (gs%ig0 > 0) vion(gs%ig0,is) = zzero
     end do
  end if

!     in case the calculation did not converge, perform
!     the orthonormalization of the  subspace: 
!     
!     a) Compute S^-1/2 = AT * D^-1/2 * A   
!     b) Compute x_new= S^-1/2 * x
!     
!     The eigenvalue array is used as temporary storage
!
  do is=1,ispn
    do irk = 1, nrk  

      call cheap_overlap_matrix(s(1,irk,is),x(1,irk,is),m,length(irk));
      if (iand(pw_params%output(1), 8) == 8) call get_timing(time(6),t3)
      !         diagonalize S, find eigenvectors A (returned in S) 
      call  occsp_diag_scal(myproc,s(1,irk,is),m,ham(irk,is)%gspace%nproc,&
            eval(1,irk,is),g(1,irk,is),length(irk)*m,rwork)
      !        compute At * D^-1/2 * A  
      do i=1,m
        eval(i,irk,is)=1.0/sqrt(eval(i,irk,is))
        do j=1,m
          mu(i+m*(j-1))=eval(i,irk,is)*conjg(s(j+m*(i-1),irk,is))
        end do
      end do
      !        compute nu = A * D^-1/2* AT 
      call mzgemm('N','N', m,m,m, zone, s(1,irk,is),m,mu(1),m,zzero, nu(1),m) 
      if (iand(pw_params%output(1), 8) == 8) call get_timing(time(3),t3)
      !        compute g = nu * x  
      if(iand(pw_params%output(1),1).eq.1) t2=gimmetime()
      call mzgemm('N','N',length(irk),m,m, zone,x(1,irk,is),length(irk),&
             nu(1),m,zzero, g(1,irk,is),length(irk)) 
      if(iand(pw_params%output(1),1).eq.1) then
        t1=t1+gimmetime()-t2
        nops = nops + 8*dble(m*m)*dble(length(irk))
      endif
      if (iand(pw_params%output(1), 8) == 8) call get_timing(time(7),t3)
      call apply_ham_wavefn(time,pw_params%output(1),y(1, irk,is),ham(irk,is),&
              g(1, irk,is), m, ews(1, 1), pw_params%nbandsfft)
      if (iand(pw_params%output(1), 8) == 8) t3 = gimmetime()
      call overlap_matrix(epsdum, g(1, irk,is), y(1, irk,is),m,length(irk))
      if (iand(pw_params%output(1), 8) == 8) call get_timing(time(6),t3)
      call  occsp_diag_scal(myproc,epsdum,m,ham(irk,is)%gspace%nproc,&
            eval(1,irk,is),gp(1,irk,is),length(irk)*m,rwork)
      if (iand(pw_params%output(1), 8) == 8) call get_timing(time(3),t3)
      if(iand(pw_params%output(1),1).eq.1) t2=gimmetime()
      call mzgemm('N','N',length(irk),m,m, zone, g(1,irk,is),length(irk),&
             epsdum(1),m, zzero, gp(1,irk,is), length(irk)) 
      if(iand(pw_params%output(1),1).eq.1) then
        t1=t1+gimmetime()-t2
        nops = nops + 8*dble(m*m)*dble(length(irk))
      endif
      if (iand(pw_params%output(1), 8) == 8) call get_timing(time(7),t3)
      !
      !     copy return eigenvectors into x array
      !
      !   x <- g
      call mzcopy(length(irk) * m, gp(1, irk,is), 1, x(1, irk, is), 1)  
      if (iand(pw_params%output(1), 8) == 8) call get_timing(time(4),t3)
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
  if(ham(1,1)%vnloc%nvecs>0 .and. &
           .not.ham(1,1)%pspot%NL_rspace(1)) deallocate(ews)
  deallocate(g)
  deallocate(gold)
  deallocate(h)
  deallocate(gp)
  deallocate(y)
  deallocate(iord)

  deallocate(gc)
  deallocate(gpc)
  deallocate(gcold)
  deallocate(s_invhalf)
  deallocate(F)
  deallocate(del_F)
  deallocate(grad_F)

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

end subroutine cgemin_metal_gcg
