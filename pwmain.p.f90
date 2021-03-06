! m4undef.m4
!
! resets various m4 commands to have a prefix m4_
! this can be achieved using option -P but only on some architectures
! any file to be preprocessed with m4 should use this set of macro definitions
!
! David Prendergast, June 6, 2006



! fft_macros.m4
!
!     fft_aux_space(aux,naux,aux2,naux2,dec1,dec2)
!     fft_local_init(N,A,LDA,LOT,DIREC)
!     fft_perm_init(N,AUX)
!
!     fft_multiple_backward(N,A,LOTSTRIDE,LOT,LOOPDUMMY,VECOUT,OFFSET)
!     fft_multiple_forward (N,A,LOTSTRIDE,LOT,LOOPDUMMY,VECOUT,OFFSET)
!
!     fft_backward(N,A,LOTSTRIDE)
!     fft_forward (N,A,LOTSTRIDE)
!
!     fft_convol(N,A,B,POS)
!
!     fft_multiple_scale(N,A,LOTSTRIDE,LOT,LOOPDUMMY,SCALE)
!
!     timeget (t0)
!     timediff(t1,t0)
!
!
!     fft_fcblock     defines blocking factor for fast convolute
!     fft_local_free  release memory
!

!-*-F90-*- 
!                                Parallel 
!                     Kleinman-Bylander Pseudopotential 
!                      Total Energy Plane Wave Program
!
!     1995/96 Bernd Pfrommer, Department of Physics,     UC Berkeley
!
!     based on a code by J.L. Martins (1990)
!
!     -----------------------------------------------------------------------
!
subroutine pwmain(t0,iwfnreuse, enew, ssum, fsum, &  
                  pw_params,altkpoints,bands,altbands,syms,energs,crys)
  !
  !     ---------------------------------------------------------------
  !  
  use all_to_all_module
  include 'use.h'
  implicit none             ! never uncomment this line. 
  include 'mpif.h'
  include 'interface.h'
  include 'all_to_all.h'
  include 'flibcalls.ph'
  include 'param.i'               ! the parameters for dimensioning the arrays
  !
  !     INPUT/OUTPUT:
  !     ------------
  !
  type(kpoint), intent(inout) :: altkpoints ! alternative BZ integration kts
  type(band), intent(inout) :: bands        ! eigvalues, occupation numbers etc
  type(band), intent(inout) :: altbands     ! another set of eigenvalues etc..
  type(symmetry), intent(inout) :: syms     ! symmetry operations
  type(energy), intent(inout) :: energs     ! energies and other information
  type(crystal), intent(inout) :: crys      ! crystal structure
  type(pw_parameter), intent(inout) :: pw_params ! plane wave parameters
  !
  !     INPUT:
  !     ------
  !
  integer, intent(in) :: &
       iwfnreuse                ! noreuse=0, firstcall=1, reuse=2
  real(dp), intent(in) :: &
       t0                      ! starting time
  !
  !     OUTPUT:
  !     -------
  !
  real(dp), intent(out) :: &
       ssum(6), &                 ! metric stress tensor Ryd
       fsum(3, crys%mxdatm, crys%ntype), &  ! forces in avec*latt. coord
       enew                       ! new energy
  !
  !     -------------------------- local variables -----------------------
  !
  type(symmetry) :: id_sym        ! symmetry structure for identity only
  type(fft_struc) :: ffts         ! information for the fast fourier transform
  type(pseudo_potential) :: pspot
  type (blochl_operator) :: blop  ! blochl operator
  type(force) :: ewaldf           ! ewald force and stress
  type(kpoint),save :: kpoints    ! BZ integration kpoints
  type(complex_gspace_array) :: &
       vion, &              ! screened ionic potential
       denc, &              ! core charge density
       denv, &              ! valence charge density
       dena, &              ! charge density, alternate
       vin,  &              ! input potential
       vout, &              ! output potential
       ddc,  &              ! derivative of core charge 
       dvql, &              ! derivative of local ionic potential
       vnl_ave              ! average nonlocal in real space for TF mixing  
  type(complex_gspace_array), save :: &
         den_del1,den_del2,den_del3  !DBR
  type(complex_gspace_array), save :: &
       den                  ! total charge density
  type(complex_gspace_array) :: &
       wavefnbs             ! wave functions for band structure
  type(complex_gspace_array), save :: &
       wavefn               ! wave functions for all kpoints, bands
  type(complex_gspace_array), save :: & !DBR
       wavefn_old               ! wave functions for all kpoints, bands 
  type(complex_gspace_array), save :: & !DBR
       wavefn_del1               ! wave functions for all kpoints, bands
  type(complex_gspace_array), save :: & !DBR
       wavefn_del2               ! wave functions for all kpoints, bands
  type(double_gspace_array) :: &
       vql, &               ! local pseudo pot in gspace for each type
       dnc                  ! the core charge in gspace for each type
  type(complex_rspace_array) :: &
       chdr                 ! work array: charge density in realspace
  type(parallel_gspace) :: &
       pot_gspace
  type(parallel_gspace), allocatable :: k_gspace(:)
  complex(dp), allocatable :: &
       chi_trace(:), e_field(:,:), &
       j_induced(:,:), &    ! induced current in gspace
       chi_g(:,:,:)         ! rho_magnetic in gspace (only for pw_nmr_plot)
  complex(dp), allocatable :: &    
       U(:), A_t(:),temp(:), rwork(:), U_temp(:)  !DBR
  real(dp), allocatable :: eval(:)
  integer len,m    
  !
  integer :: &
       i, j, k, n, is, irk, info, &
       nproc, &             ! number of processors started up
       myproc, &            ! number of this processor: 0,1,...nproc-1
       nprow, npcol, &      ! number of processor row or columns
       id9, idspin, ierr, &
       mtrxdim, &           ! total size of hamilton matrix
       maxlength, &         ! maximum number of plane waves at any k-point
       ipr, &               ! print flag
       ijob                 ! job counter
  !
  integer :: &
       nfn                  ! number of functions to FFT simultaneously
  integer, parameter :: &
       nsafety = 5          ! number of extra states to start diag with
  !
  integer :: iscr2  ! hjchoi
  real(dp) :: &
       ddum, &              ! double dummy
       factor,tstart
  real(dp), allocatable :: &
       p_weights(:,:,:,:,:,:) ! weights for angular momentum DOS
  character(len=11) :: wfnname
  character(len=255) :: strbuf1     ! dummy string buffer 
  character(len=9) :: cdname      ! file name for cdfile
  !
  !     maps job number to intelligible strings
  !
! davegp
!  character(len=22), parameter :: jobstring(6) = (/ &
!       'SELF-CONSISTENT FIELD ', 'NMR SHIFT             ', &
!       'NMR PLOT              ', 'BAND STRUCTURE        ', &
!       'PLOT POTENTIAL        ', 'TDDFT                 ' /)
  character(len=22), parameter :: jobstring(7) = (/ &
       'SELF-CONSISTENT FIELD ', 'NMR SHIFT             ', &
       'NMR PLOT              ', 'BAND STRUCTURE        ', &
       'PLOT POTENTIAL        ', 'TDDFT                 ', &
       'NON-SELF-CONSISTENT   ' /)
! davegp
  logical :: logidum, &
       lkwfn, &             ! whether wave function should be kept
       lrwfn, &             ! whether wave function should be reused
       did_nmr              ! did we calculate the magnetic susceptability
!
!  Extra variables for wannier=1
!
  !
  !     ################## external functions and subroutines #############
  !
!-------------------------------------
! wannier orbitals
   integer::n_band,n_orb,n_kp,nband0_wan,nband1_wan,ncenter
!-------------------------------------

  real(dp), external :: &
       gimmetime            ! function which returns elapsed CPU time
  !     ###  added by bshih ###
  character(len=11) :: wanname
  !
  !     --------- start up parallel stuff if not serial----------
  !
  
  call mpi_comm_rank(MPI_COMM_WORLD, myproc, info)
  call mpi_comm_size(MPI_COMM_WORLD, nproc, info)
  
  call layout_scalapack(100, info, nproc, nprow, npcol)
  call blacs_get(-1, 0, blacs_context) ! get the default blacs context
  call blacs_gridinit(blacs_context, 'R', nprow, npcol) ! commit to  grid
  
  
  call allocate_shmem_work(nproc, 1)
  !
  !     ------------------------------------
  !           FAILSAFE OUTPUT
  !
  enew = dzero
  ssum = dzero
  fsum = dzero
  !
  !     ----------------------------------------------------------------
  !           TRANSFER DATA FROM THE PARAMETER LIST TO STRUCTURES 
  !
  lrwfn = .false.
  lkwfn = .false.
  if (iwfnreuse > 1) lrwfn = .true. ! can reuse wave function
  if (iwfnreuse > 0) lkwfn = .true. ! should keep wave function

  fsum = dzero         ! to avoid crashes if forces are not computed
 
  pw_params%epsmag = pw_params%epsdiag      ! set same accuracy as in diag.
  pw_params%maxitcgsol= pw_params%maxitdiag ! and same number of iterations
  !
  !----------------------------------------------------------------------------
  !
  !                 OPEN OUTPUT FILES
  !
  !
  !     --- open output files -------
  !
  if (myproc == 0) then
    if (iand(pw_params%output(1), 268435456 ) == 268435456 .and. &
          pw_params%ilog > 1) then
      open(9, file = 'PW_LOG', position = 'append', &
      status = 'unknown', form = 'formatted')
    else     
      open(9, file = 'OUT', status = 'unknown', form = 'formatted')
    end if
!    open(34, file = 'PW_LOG', position = 'append', &
!          status = 'unknown', form = 'formatted')
  else
     open(9, file = '/dev/null', form = 'formatted')
  end if
  !
  !
  !
  call tpage(nproc)              ! prints title page
  !
  !     look at crystal structure. This has to be done always.
  !     
!  call create_crystal(crys, mxdatm, ntype, coord, &
!       natom, nameat, avec, szet, mtsphere, dipole) 
  pspot%mxdlqp = mxdlqp
  !
  
  call mpi_barrier(MPI_COMM_WORLD, info)
  
  !     
  !     gets crystal data and determines symmetry
  !
  crys%ilsdau=pw_params%lsdau
  crys%igwout=0
  if ((iand(pw_params%output(1), 524288) == 524288).or. &
     (iand(pw_params%output(1), 1048576) == 1048576)) crys%igwout=1

  call crstl(pw_params%output(1), crys, syms, pw_params%miscflag,pw_params%ilog)
  !      call symmetrize_crystal(crys, syms)


  ! davegp
  write(9,*) ' out of crstl'
  call myflush(9)

  if (.not. lrwfn) then
     kpoints%nrk = abs( altkpoints%nrk )
     kpoints%grid = altkpoints%grid
     kpoints%shift = altkpoints%shift
  end if

  kpoints%generated = .false. 
  if ( altkpoints%nrk <= 0) kpoints%generated = .true.
  kpoints%reduced = .false.
  if ( altkpoints%nrk   >= 0) kpoints%reduced = .true.
  altkpoints%generated = kpoints%generated
  altkpoints%reduced = kpoints%reduced

  altkpoints%nrk = abs(altkpoints%nrk) 

  call print_pwparam(1, pw_params)

  call ball_and_stick(pw_params%output(1), crys, myproc)
  call neighbor(crys, pw_params)
  ! davegp
  write(9,*) ' out of ball_and_stick and neighbor'
  call myflush(9)

  nfn = pw_params%nbandsfft
  !
  !     the main loop. Implement and document new subroutines from here
  !
  !
  ! davegp
  write(9,*) ' about to start SCF loop'
  call myflush(9)

  do ijob = 1, pw_params%njobs
     write(9, 30) trim(jobstring(pw_params%joblist(ijob) + 1))
! davegp
! non-self-consistent field calculations piggy-back on the self-consistent
! field algorithm without updating the potential
     pw_params%nscf = .false.
     if( pw_params%joblist(ijob) == 6 ) then
       ! set nscf flag to true
       pw_params%nscf = .true.
       ! masquerade as SCF
       pw_params%joblist(ijob) = 0
     endif
! davegp
     tstart = gimmetime()
     select case(pw_params%joblist(ijob))
        !
        !        ------------ just the scf loop  ----------------------------
     case(0)
        !
        !           generate large gspace for the potential
        !
        call generate_potential_gspace(pot_gspace, crys, syms, pw_params, &
             myproc, nproc)
       
  ! davegp
  write(9,*) ' about to start SCF loop'
  call myflush(9)

        if (pw_params%nbandsfft .eq. -1) then
          pw_params%nbandsfft = min(max(int((pot_gspace%nproc**2*2000)/pot_gspace%totlength),1),max(bands%min(1),bands%min(2)) + nsafety )
          nfn = pw_params%nbandsfft
          write(9,*) 'NUMBER of COMPUTED bands for FFT',  pw_params%nbandsfft
        end if 
        !  
        !           allocate work arrays for fft
        !
        call create_fftstruc(ffts, pot_gspace, nfn)
  ! davegp
  write(9,*) ' create_fftstruc'
  call myflush(9)

        !
        !           read pseudopotential and set up valence/core charge
        !           (the ddc,dvql, and dnc arrays are for force/stress calc)
        !
        call read_pseudo(1, crys%nspin, pw_params, pot_gspace, energs, &
             pspot, crys, vion, denc, denv, ddc, dvql, vql,dnc,vnl_ave,ffts)
  ! davegp
  write(9,*) ' read_pseudo'
  call myflush(9)

        !  
        !           calculate ewald energy and forces 
        !
        call ewald(pw_params%output(1), crys, pot_gspace, energs%ewald, ewaldf)
        !
        !           screen the ionic potential with the valence and core
        !           charge density
        !           also creates arrays vin, den, chdr.
        !
        if (0<=pw_params%extrapolation .and. pw_params%extrapolation<10) then
          if (.not. lrwfn ) then ! &
            call create_gspace_array(den_del1,pot_gspace%length,1,crys%nspin)
            call create_gspace_array(den_del2,pot_gspace%length,1,crys%nspin)
            call create_gspace_array(den_del3,pot_gspace%length,1,crys%nspin)
            den_del1%data=zzero;  den_del2%data=zzero;  den_del3%data=zzero;
          else
            if (pw_params%ilog .ge. 2) then
              call mzcopy(pot_gspace%length* crys%nspin,& 
                denv%data(1,1,1),1, den%data(1,1,1),1)
              call mzaxpy(pot_gspace%length* crys%nspin,zone,& 
                  den_del1%data(1,1,1),1, den%data(1,1,1),1)
              call mzaxpy(pot_gspace%length* crys%nspin,&
                      -cmplx(pw_params%MD_alpha,dzero,dp),& 
                 den_del2%data(1,1,1),1, den%data(1,1,1),1)
              call mzaxpy(pot_gspace%length* crys%nspin,&
                  -cmplx(pw_params%MD_beta,dzero,dp),&
                   den_del3%data(1,1,1),1, den%data(1,1,1),1)
            else if (pw_params%ilog .eq. 1) then
              call mzcopy(pot_gspace%length* crys%nspin,& 
                  denv%data(1,1,1),1, den%data(1,1,1),1)
              call mzaxpy(pot_gspace%length* crys%nspin,zone,& 
                  den_del1%data(1,1,1),1, den%data(1,1,1),1)  
            end if     
          end if
        end if
        !
        !          reads or calculates integration k-points, initializes bands.
        !
        if (.not. lrwfn) &
             call generate_kpoints(1, kpoints, crys, syms, bands, &
             pw_params%output(1), pw_params%q, nsafety)
        !
        !           now set up all the small gspaces at the various kpoints
        !

! output kmap info

        if (myproc.eq.0) then
           open(144,file="kmap.dat")
           write(144,*) kpoints%grid(:)
           write(144,*) kpoints%nrk,kpoints%nrk_fbz
           do i=1,kpoints%nrk_fbz
              write(144,1440) i,kpoints%rk_fbz(:,i),abs(kpoints%kmap(i)),kpoints%rk(:,abs(kpoints%kmap(i)))
           end do
           close(144)
        end if
           
 1440   format(i5,3f12.8,i5,3f12.8)


        allocate(k_gspace(kpoints%nrk))
        call generate_k_gspaces(k_gspace, maxlength, pot_gspace, kpoints, &
             crys, syms, pw_params%emax, pw_params%output(1))


        if (lrwfn .and. pw_params%extrapolation .gt.0) then ! &
          if (pw_params%ilog .ge. 3) then
            call mzaxpy(maxlength * bands%max*kpoints%nrk* crys%nspin,&
                cmplx(pw_params%MD_alpha,dzero,dp),& 
                 wavefn_del1%data(1,1,1),1, wavefn%data(1,1,1),1)
           call mzaxpy(maxlength * bands%max*kpoints%nrk* crys%nspin,&
                  cmplx(pw_params%MD_beta,dzero,dp),& 
                 wavefn_del2%data(1,1,1),1, wavefn%data(1,1,1),1)
          else if (pw_params%ilog .eq. 2) then
            call mzaxpy(maxlength * bands%max*kpoints%nrk* crys%nspin,zone,& 
                 wavefn_del1%data(1,1,1),1, wavefn%data(1,1,1),1) 
          end if

          if (pw_params%extrapolation .ge. 10) then
            allocate(U(bands%max*bands%max)); U=zzero
            do is=1,crys%nspin
              do irk = 1, kpoints%nrk  
                m= bands%nband(irk,is)
                call cheap_overlap_matrix(U(1), wavefn%data(1,irk,is), m, &
                      k_gspace(irk)%length)
                call occsp_sqrtinv_scal('U','n',pot_gspace%myproc,U(1),m,&
                      k_gspace(irk)%nproc)  
                call mztrmm('R','U','N','N',k_gspace(irk)%length,m, zone,&
                 U(1),m,wavefn%data(1,irk,is),k_gspace(irk)%length)
              end do
            end do

            call charge(ffts, pw_params, pot_gspace, k_gspace, bands, &
               kpoints, wavefn, den%data(1, 1, 1), energs,crys)
            deallocate(U)
          end if
        end if

        call screen(pw_params%iscr, lrwfn, pot_gspace, crys, ffts, energs,&
            pw_params,   vion, denv, denc, vin, den, chdr)
        !
        call create_gspace_array(vout, pot_gspace%length, 1, crys%nspin)
        !
        !       allocate space for wave functions, and do selfconsistent loop
        if (.not. lrwfn) then ! &
          call create_gspace_array(wavefn, maxlength * bands%max, &
             kpoints%nrk, crys%nspin)
          if (0<= pw_params%extrapolation) then
             call create_gspace_array(wavefn_old, maxlength * bands%max, &
             kpoints%nrk, crys%nspin)         
             call create_gspace_array(wavefn_del1, maxlength * bands%max, &
             kpoints%nrk, crys%nspin)  
             call create_gspace_array(wavefn_del2, maxlength * bands%max, &
             kpoints%nrk, crys%nspin)
             wavefn%data=zzero;
             wavefn_old%data=zzero;wavefn_del1%data=zzero; wavefn_del2%data=zzero;
          end if
        end if

        if (iand(pw_params%optimize,8) == 8 .or. &
             iand(pw_params%optimize,64) == 64) then ! direct minimization
           call direct_emin(t0,pot_gspace, k_gspace(1), energs, pw_params, &
                crys, syms, ffts, pspot, kpoints, bands, wavefn, lrwfn, &
                vion%data(1, 1, 1), vin%data(1, 1, 1), vout%data(1, 1, 1), &
                den%data(1, 1, 1), denc%data(1, 1, 1), chdr%data(1, 1, 1))
        else                ! old scf iteration
           call scfloop(t0,1, pot_gspace, k_gspace(1), energs, pw_params, &
                crys, syms, ffts, pspot, kpoints, bands, wavefn,lrwfn, &
                vion, vin, vout, den, denc, chdr,vnl_ave)
        end if
        energs%inputfermi = energs%efermi(1)
        if (0<= pw_params%extrapolation .and. pw_params%extrapolation <10) then
          if (pw_params%ilog .ge. 2) then
            call mzcopy(pot_gspace%length* crys%nspin,& 
                 den_del2%data(1,1,1),1, den_del3%data(1,1,1),1)
          else
            den_del3%data=zzero
          end if
          if (pw_params%ilog .ge. 1) then
            call mzcopy(pot_gspace%length* crys%nspin,& 
                 den_del1%data(1,1,1),1, den_del2%data(1,1,1),1)
            call mzcopy(pot_gspace%length* crys%nspin,& 
                 den%data(1,1,1),1, den_del1%data(1,1,1),1)
            call mzaxpy(pot_gspace%length* crys%nspin,zmone,& 
                 denv%data(1,1,1),1, den_del1%data(1,1,1),1)
            call mzaxpy(pot_gspace%length* crys%nspin,zmone,& 
                 den_del1%data(1,1,1),1, den_del2%data(1,1,1),1) 
          end if
          if (pw_params%ilog .eq. 0) then
            call mzcopy(pot_gspace%length* crys%nspin,& 
                 den%data(1,1,1),1, den_del1%data(1,1,1),1)
            call mzaxpy(pot_gspace%length* crys%nspin,zmone,& 
                 denv%data(1,1,1),1, den_del1%data(1,1,1),1)
          end if
        end if


!-------------------------------------------------------------------
! calculate wavefunction overlap, for constructing wannier functions

     if((pw_params%wannier.eq.1).or.(pw_params%wannier.eq.2)) then

        call wannier(k_gspace,energs,crys,syms,kpoints,bands,wavefn,nband0_wan,nband1_wan)

     if(pw_params%wannier.eq.2) then
!----------------------------------------------
! output wavefunction

       do irk=1, kpoints%nrk
         do is =1, crys%nspin
            write(wanname,1200) irk, is
            if (myproc == 0) then
              open(unit =18 , file=wanname , status= 'unknown', form='unformatted')
              write(18) ffts%fftsize, irk, nband1_wan-nband0_wan
!              close(18)
            end if

           do n = (1+nband0_wan), nband1_wan

             call dxplot_wannier(4, pw_params%output(1), &
                      wavefn%data(k_gspace(irk)%length *(n - 1) + 1, &
                      irk, is), k_gspace(irk), ffts, crys, irk, is)

           end do
           if (myproc == 0 )  close(18)
         end do
      end do
1200 format ('UNK',i5.5,'.',i1)
!1201 format (5i5)   ! for (formatted UNK files)
! END 
!----------------------------------------------

     end if
     end if

!     if(pw_params%wannier.eq.2) then
!        open(111,file="Umn.dat")
!        read(111,*) n_orb,n_kp,n_band
!     if(n_kp.ne.kpoints%nrk) then
!        write(9,*) "k points do not match"
!        stop
!     end if
!     end if

!-------------------------------------------------------------------

        !
        ! print out gw screening
        !
        if (iand(pw_params%output(1),2097152) == 2097152) then
           factor = 8.0d0 * pi / crys%vcell
           !
           ! in GW calculations, the exchange-correlation potential
           ! should not include partial core correction. Assign denc%data = 0.
           !

! BUG here, P. Zhang
!
!           denc%data(1:pot_gspace%length, 1, 1:crys%nspin) = zzero
           denc%data(1:pot_gspace%length, 1, 1) = zzero
           call velect(1, pot_gspace, crys, ffts, energs, pw_params, &
                denc%data(1, 1, 1), vout%data(1, 1, 1), den%data(1, 1, 1), &
                chdr%data(1, 1, 1))
           ! subtract out core charge
           ! subtract out hartree potential
           do i = 1,pot_gspace%length
              if (pot_gspace%ekin(i) > dzero) then
                 if (crys%nspin == 1) then
                    vout%data(i, 1, 1) = vout%data(i, 1, 1) - &
                         factor / pot_gspace%ekin(i) * den%data(i, 1, 1)
                 end if
                 !     when nspin=2, we should extract spin part separately and
                 !     den%data(i,1,1) is now the sum of charge density.
                 if (crys%nspin == 2) then
                    vout%data(i, 1, 1) = vout%data(i, 1, 1) - &
                         factor / pot_gspace%ekin(i) * &
                         (den%data(i, 1, 1) + den%data(i, 1, 2))
                    vout%data(i, 1, 2) = vout%data(i, 1, 2) - &
                         factor / pot_gspace%ekin(i) * &
                         (den%data(i, 1, 1) + den%data(i, 1, 2))
                 end if
              end if
           end do
           call writegsdat2(0, pot_gspace, vout%data(1, 1, 1), &
                1, 1, crys%nspin, 'VXC', 3)

! PZ output charge density without semicore electrons
!           call charge2(ffts, pw_params, pot_gspace, k_gspace, bands, &
!               kpoints, wavefn, den%data(1, 1, 1), energs,crys)

           call writegsdat2(0, pot_gspace, den%data(1, 1, 1), &
                1, 1, crys%nspin, 'CD95', 4)

        end if
        !
        if (iand(pw_params%output(1), 524288) == 524288 .and. pw_params%gw_sum .eq. 0) then
             write(9,*) "before gwout"
             call gwout('GWR',crys, syms, bands, pot_gspace, k_gspace, &
             wavefn, kpoints, 1)
             write(9,*) "after gwout"
             call myflush(9)
             call myflush(6)
        end if
        if (iand(pw_params%output(1), 524288) == 524288 .and. pw_params%gw_sum .ne. 0) then
             write(9,*) "before gwout, use new format wavefunctions"
             call gwoutnewfmt('GWR',crys, syms, bands, pot_gspace, k_gspace, &
             wavefn, kpoints, 1, pw_params)
             write(9,*) "after gwout"
             call myflush(9)
             call myflush(6)
        end if

        if (iand(pw_params%output(1), 1048576) == 1048576 .and. pw_params%gw_sum .eq. 0) then
             write(9,*) "before gwout"
             call gwout('GWC', crys, syms, bands, pot_gspace, k_gspace, &
             wavefn, kpoints,2)
             write(9,*) "after gwout"
             call myflush(9)
             call myflush(6)
        end if
        if (iand(pw_params%output(1), 1048576) == 1048576 .and. pw_params%gw_sum .ne. 0) then
             write(9,*) "before gwout"
             call gwoutnewfmt('GWC', crys, syms, bands, pot_gspace, k_gspace, &
             wavefn, kpoints,2, pw_params)
             write(9,*) "after gwout"
             call myflush(9)
             call myflush(6)
        end if
        !     
        !     output data to Balazs VQMC code
        !
        if (iand(pw_params%output(1), 2) == 2) &
             call vqmc('VQMC', crys, syms, bands, pot_gspace, k_gspace, &
             wavefn, kpoints, energs%efermi(1), pw_params%smearing / ryd)
        !
        !           compute forces and stress
        !

! If output GW screening, the potentials are modified, cannot calculate forces and stress

        if ((.not. iand(pw_params%output(4), 8) == 8).and.(.not.(iand(pw_params%output(1),2097152) == 2097152))) then


        write(9,*) "before force_stress"
        write(6,*) "before force_stress"
        call force_stress(crys, syms, bands, kpoints, pspot, energs, &
             pw_params, pot_gspace, k_gspace, vion%data(1, 1, 1), &
             vin%data(1, 1, 1), vout%data(1, 1, 1), den%data(1, 1, 1), &
             denc%data(1, 1, 1), dvql%data(1, 1, 1), ddc%data(1, 1, 1), &
             vql%data(1, 1, 1), dnc%data(1, 1, 1), wavefn, ewaldf, ssum, &
             fsum, ffts)

        write(6,*) "after force_stress"
        write(9,*) "after force_stress"
        call myflush(9)
        call myflush(6)
        end if

        call logfile(ipr, myproc, nproc, pw_params%ilog, gimmetime() - tstart,&
                 fsum, ssum, crys, pw_params, syms, kpoints, energs)
        !
        !           plot charge density 
        !
        if (iand(pw_params%output(1), 16) == 16) then
           if (crys%nspin == 1) then
              call dxplot(1, pw_params%output(1), den%data(1, 1, 1), &
                   pot_gspace, ffts,  crys, 'CHARGE___',9)
              call lineplot(1, crys, done / crys%vcell, den%data(1, 1, 1), &
                   pot_gspace, pw_params%nlineplot, pw_params%lineplot, &
                   pw_params%bra, pw_params%ket, 1,'CHARGE___',9)
           else
              call dxplot(1, pw_params%output(1), den%data(1, 1, 1), &
                   pot_gspace, ffts, crys, 'CHARGE_UP',9)
              call dxplot(1, pw_params%output(1), den%data(1, 1, 2), &
                   pot_gspace, ffts, crys, 'CHARGE_DN',9)
              call lineplot(1, crys, done, den%data(1, 1, 1), &
                   pot_gspace, pw_params%nlineplot, pw_params%lineplot, &
                   pw_params%bra, pw_params%ket, 1, 'CHARGE_UP',9)
              call lineplot(1, crys, done, den%data(1, 1, 2), &
                   pot_gspace, pw_params%nlineplot, pw_params%lineplot, &
                   pw_params%bra, pw_params%ket, 1,  'CHARGE_DN',9)
           end if

           do i = 1, pw_params%nenergywindow
              call resolved_charge(pw_params%energywindow(1+ (i-1)*2), &
                   ffts, pw_params, pot_gspace, k_gspace(1), &
                   bands, kpoints, wavefn, den%data(1, 1, 1), energs)
              do is = 1, crys%nspin
                 write(cdname, 20) is, i
                 call dxplot(1, pw_params%output(1), &
                   den%data(1, 1, is), pot_gspace, ffts, crys, cdname,9)
                 call lineplot(1, crys, done / crys%vcell, &
                   den%data(1, 1, is), pot_gspace, pw_params%nlineplot, &
                   pw_params%lineplot, pw_params%bra, pw_params%ket, 1, cdname,9)
              end do
           end do
        end if
        !
        !           compute and plot electric field
        !
        if (iand(pw_params%output(1),131072) == 131072) then
           call electric_field(pw_params, crys, den%data(1, 1, 1), &
                vion%data(1, 1, 1), vin%data(1, 1, 1), pot_gspace, &
                pw_params%nlineplot, pw_params%lineplot, pw_params%bra, &
                pw_params%ket)
        end if
        !
        !           plot wave functions
        !
        if (iand(pw_params%output(1),256) == 256) then
           do is = 1, crys%nspin
              do j = 1, pw_params%nplotwv
                 irk = pw_params%plotwv(1+ (j-1)*2) ! get kpoint number
                 if (irk > kpoints%nrk) then
                    write(9, 530) 'kpoint', irk, &
                         'number of kpoints', kpoints%nrk
                    irk = kpoints%nrk
                 end if
                 n = pw_params%plotwv(2+ (j-1)*2) ! get band number
                 if (n > bands%nband(irk, is)) then
                    write (9,530) 'band', n, 'number of bands', &
                         bands%nband(irk, is)
                    n = bands%nband(irk, is)
                 end if
                 ! Write the square of modulus of the wavefunctions
                 write(strbuf1, '(a3,i2.2,a1,i3.3,a1,i3.3)') 'WF.',is,'.',irk, '.', n
                 call dxplot(3, pw_params%output(1), &
                      wavefn%data(k_gspace(irk)%length * &
                      (n - 1) + 1, irk, is), k_gspace(irk), &
                      ffts, crys, strbuf1,13)
                 ! Write the real part of the wavefunctions
                 write(strbuf1, '(a3,i2.2,a1,i3.3,a1,i3.3)') 'WR.',is,'.', irk, '.', n
                 call dxplot(1,pw_params%output(1), &
                      wavefn%data(k_gspace(irk)%length * &
                      (n - 1) + 1, irk, is), k_gspace(irk), &
                      ffts, crys, strbuf1,13)
                 ! Write the imaginary part of the wavefunctions
                 write(strbuf1, '(a3,i2.2,a1,i3.3,a1,i3.3)') 'WI.',is,'.', irk, '.', n
                 call dxplot(2, pw_params%output(1), &
                      wavefn%data(k_gspace(irk)%length * &
                      (n - 1) + 1, irk, is), k_gspace(irk), &
                      ffts, crys, strbuf1,13)
                 call lineplot(3, crys, done, &
                      wavefn%data(k_gspace(irk)%length * (n - 1) + 1, &
                      irk, is), k_gspace(irk), pw_params%nlineplot, &
                      pw_params%lineplot, pw_params%bra, pw_params%ket, 1,&
                      strbuf1,13)
              end do
           end do
        end if
        !
        !           do lineplot and printout for the data explorer
        !
        if (iand(pw_params%output(1), 262144) == 262144) then
           if (crys%nspin == 1) then
             call lineplot(1, crys, done, vion%data(1, 1, 1), &
                pot_gspace, pw_params%nlineplot, pw_params%lineplot, &
                 pw_params%bra, pw_params%ket, 1, 'VEFFAVG',7)
             call dxplot(1, pw_params%output(1), vion%data(1, 1, 1), &
                pot_gspace, ffts, crys, 'VEFFAVG',7)
           else
             call lineplot(1, crys, done, vion%data(1, 1, 1), &
                pot_gspace, pw_params%nlineplot, pw_params%lineplot, &
                 pw_params%bra, pw_params%ket, 1, 'V_UP',4)
             call dxplot(1, pw_params%output(1), vion%data(1, 1, 1), &
                pot_gspace, ffts, crys, 'V_UP',4)
             call lineplot(1, crys, done, vion%data(1, 1, 2), &
                pot_gspace, pw_params%nlineplot, pw_params%lineplot, &
                 pw_params%bra, pw_params%ket, 1, 'V_DOWN',6)
             call dxplot(1, pw_params%output(1), vion%data(1, 1, 2), &
                pot_gspace, ffts, crys, 'V_DOWN',6)             
           end if

        end if
        !
        !                 compute the weights for the angular momentum DOS
        !
        allocate(p_weights(llmax, crys%mxdatm, crys%ntype, bands%max, &
             bands%nrk, crys%nspin)) 
!        if (iand(pw_params%output(1),1024) == 1024) then
!           call angular(p_weights, k_gspace, ffts, crys, bands, &
!                kpoints, wavefn, syms, pw_params, pw_params%nrg) 
!        end if
!        if (iand(pw_params%output(1),128) == 128 .or. &
!             iand(pw_params%output(1),1024) == 1024) then
           !
           !              compute and print density of states
           !

!--------------------------------------------------------------------
!        if (iand(pw_params%output(1), 4194304) == 4194304) &
!
!             call project(myproc, pw_params%emax, bands%min(1), crys, &
!             wavefn, kpoints%nrk, k_gspace(1), 'project.dat')
!
! PZ

        if (iand(pw_params%output(1), 4194304) == 4194304) then
           call project_pz(myproc,nproc, crys, syms,wavefn, kpoints, k_gspace,bands,energs,1)


        end if
        if (iand(pw_params%output(1),128) == 128) then
           call dos(myproc, nproc, p_weights, crys, bands, energs, &
                kpoints, pw_params)
        end if
        if (iand(pw_params%output(1),2048) == 2048) then
           call write_eigvals(myproc, pw_params, crys, &
                energs, kpoints, bands, p_weights) 
        end if
         if (iand(pw_params%output(4),4) == 4) then
           call fermi_surface(myproc, pw_params, crys, &
                energs, kpoints, bands) 
        end if
        deallocate(p_weights)
        !
        ! align subspace
        !
        if (0<= pw_params%extrapolation) then
          if (pw_params%ilog .ge. 2) then 
            allocate(U(bands%max*bands%max), A_t(bands%max*bands%max),&
                    temp(bands%max*maxlength), rwork(max(1, 3*bands%max-2)),&
                     U_temp(bands%max*bands%max), eval(bands%max)    )

            do is=1,crys%nspin
             do irk=1,kpoints%nrk

              m=bands%nband(irk, is)
              len=k_gspace(irk)%length 
              call overlap_matrix(U(1), wavefn%data(1,irk,is), & 
                wavefn_old%data(1,irk, is),m,len) 
              call mzgemm('C','N',m,m,m,zone,U(1),m, &
                   U(1),m,zzero,A_t(1),m)
              call  occsp_diag_scal(pot_gspace%myproc,A_t,m,pot_gspace%nproc,&
                  eval, U_temp,m*m,rwork)  
              call mzgemm('N','N',m,m,m,zone,U(1),m, &
                   A_t(1),m,zzero,U_temp(1),m)
              do i=1,m
                call mzdscal(m,done/sqrt(eval(i)),U_temp(1+(i-1)*m),1)
              end do
              call mzgemm('N','C',m,m,m,zone,U_temp(1),m, &
                   A_t(1),m,zzero,U(1),m)         

              call mzgemm('N','C',len,m,m,zone,wavefn_del1%data(1,irk,is),len,&
                   U(1),m,zzero,temp(1),len)
              call mzcopy(len*m,temp(1),1, wavefn_del1%data(1,irk,is),1)

              call mzgemm('N','C',len,m,m,zone,wavefn_old%data(1,irk,is),len, &
                   U(1),m,zzero,temp(1),len)
              call mzcopy(len*m,temp(1),1, wavefn_old%data(1,irk,is),1)
             end do
            end do
            call mzcopy(maxlength * bands%max*kpoints%nrk* crys%nspin,& 
                 wavefn_del1%data(1,1,1),1, wavefn_del2%data(1,1,1),1)
            call mzcopy(maxlength * bands%max*kpoints%nrk* crys%nspin,& 
                 wavefn%data(1,1,1),1, wavefn_del1%data(1,1,1),1)
            call mzaxpy(maxlength * bands%max*kpoints%nrk* crys%nspin,zmone,& 
                 wavefn_old%data(1,1,1),1, wavefn_del1%data(1,1,1),1)
            call mzcopy(maxlength * bands%max*kpoints%nrk* crys%nspin,& 
                 wavefn%data(1,1,1),1, wavefn_old%data(1,1,1),1)
            deallocate(U, A_t,temp, rwork,U_temp,eval)
          end if
          if (pw_params%ilog .eq. 1) then 
           allocate(U(bands%max*bands%max), A_t(bands%max*bands%max),&
                    temp(bands%max*maxlength), rwork(max(1, 3*bands%max-2)),&
                     U_temp(bands%max*bands%max), eval(bands%max)    )

           do is=1,crys%nspin
            do irk=1,kpoints%nrk
             m=bands%nband(irk, is)
             len=k_gspace(irk)%length 
             call overlap_matrix(U(1), wavefn%data(1,irk,is), & 
             wavefn_old%data(1,irk, is),m,len) 
             call mzgemm('C','N',m,m,m,zone,U(1),m, &
                   U(1),m,zzero,A_t(1),m)
             call  occsp_diag_scal(pot_gspace%myproc,A_t,m,pot_gspace%nproc,&
                  eval, U_temp,m*m,rwork)  
             call mzgemm('N','N',m,m,m,zone,U(1),m, &
                   A_t(1),m,zzero,U_temp(1),m)
             do i=1,m
               call mzdscal(m,done/sqrt(eval(i)),U_temp(1+(i-1)*m),1)
             end do
             call mzgemm('N','C',m,m,m,zone,U_temp(1),m, &
                   A_t(1),m,zzero,U(1),m)         

             call mzgemm('N','C',len,m,m,zone,wavefn_old%data(1,irk,is),len, &
                   U(1),m,zzero,temp(1),len)
             call mzcopy(len*m,temp(1),1, wavefn_old%data(1,irk,is),1)
            end do
           end do
           call mzcopy(maxlength * bands%max*kpoints%nrk* crys%nspin,& 
                 wavefn%data(1,1,1),1, wavefn_del1%data(1,1,1),1)
           call mzaxpy(maxlength * bands%max*kpoints%nrk* crys%nspin,zmone,& 
                 wavefn_old%data(1,1,1),1, wavefn_del1%data(1,1,1),1)
           call mzcopy(maxlength * bands%max*kpoints%nrk* crys%nspin,& 
                 wavefn%data(1,1,1),1, wavefn_old%data(1,1,1),1)

           deallocate(U, A_t,temp, rwork,U_temp,eval)
          end if

          if (pw_params%ilog .eq. 0) then 
             call mzcopy(maxlength * bands%max*kpoints%nrk* crys%nspin,& 
                 wavefn%data(1,1,1),1, wavefn_old%data(1,1,1),1)
          end if
        end if
        !
        !           deallocate all the arrays which have dynamical allocation
        !
        call destroy(ffts)
        call destroy(ewaldf)
        if (.not. lkwfn) call destroy(bands)
        call destroy(pspot)
        call destroy(pot_gspace)
        do irk = 1, kpoints%nrk
           call destroy(k_gspace(irk))
        end do
        deallocate(k_gspace)
        if (.not. lkwfn) call destroy(kpoints)
        !     
        !           deallocate all the gspace arrays
        !
        if (.not. lkwfn) call destroy(wavefn)
        call destroy(vion)
        call destroy(denc)
        call destroy(denv)
        call destroy(ddc)
        call destroy(dvql)
        call destroy(vql)
        call destroy(dnc)
        call destroy(vin)
        if (.not. lkwfn) call destroy(den)
        call destroy(chdr)
        call destroy(vout)
        ! beginNMR  - don't touch this line!!
        !    **************************************************     
        !    *** Compute  t h e   N M R   s h i f t  ***
        !    **************************************************
!pub case(1); call mystop( ' NMR option not availible' )
     case(1)  
        !
        !    Generate large gspace for the potential
        !
        call generate_potential_gspace(pot_gspace, crys, syms,pw_params, &
             myproc, nproc)
        !
        !    Allocate work arrays for fft
        !
        call create_fftstruc(ffts, pot_gspace, nfn)
        !
        !    Read pseudopotential and set up valence/core charge
        !      (the ddc,dvql, and dnc arrays are for force/stress calc)
        !
        call read_pseudo(1, crys%nspin, pw_params, pot_gspace, energs, &
             pspot, crys, vion, denc, denv, ddc, dvql, vql,dnc,vnl_ave,ffts)
        ! 
        !    Read in the wavefunctions (ae/ps) and create projectors
        !      and matrix elements for the Blochl operator
        !
        call read_blochl_operator(pot_gspace,blop,pspot,crys) 
        !
        !   Screen the ionic potential with the valence and core charge 
        !     density, also creates arrays vin, den, chdr.
        !
        call screen(4, .false., pot_gspace, crys, ffts, energs, pw_params, &
             vion, denv, denc, vin, dena, chdr)
        !
        ! Estimate memory usage
        !
        write(9,*)
        write(9,*) 'NMR Memory Estimate'
        write(9,150) 22.75*pot_gspace%totlength*16/(1024*1024)
        write(9,151) ((5.0+pspot%ntype)*pot_gspace%totlength+ &
            (ffts%r_size*nfn)+(ffts%sendbufsize*nfn)+ &
            (ffts%recbufsize*nfn)+(ffts%t1bufsize * nfn)+&
            (ffts%unfoldbufsize*nfn))*16/(1024*1024)
        write(9,152) 0.875*bands%min(1)*pot_gspace%totlength*16/(1024*1024)
        write(9,153) 1.625*pspot%nanl*pot_gspace%totlength*16/(1024*1024)
        write(9,154) 0.645*bands%min(1)*pot_gspace%totlength*16/(1024*1024)
        write(9,155)
        write(9,156) ((1.625*pspot%nanl+0.875*bands%min(1)+0.645*bands%min(1) &
             +22.75+7.0+pspot%ntype) *pot_gspace%totlength +&
             (ffts%r_size*nfn)+(ffts%sendbufsize*nfn)+ &
             (ffts%recbufsize*nfn)+(ffts%t1bufsize * nfn)+&
             (ffts%unfoldbufsize*nfn))*16/(1024*1024)


        !
        !   Deallocate arrays that are no longer needed
        ! 
        !
        call destroy(denc)
        call destroy(denv)
        call destroy(ddc)
        call destroy(dvql)
        call destroy(vql)
        call destroy(dnc)
        call destroy(vin)
        call destroy(dena)
        call destroy(chdr)
        !
        !   Reads or calculates integration k-points, initializes bands.
        !
        call generate_kpoints(1, altkpoints, crys, syms, altbands, &
             pw_params%output(1), pw_params%q, nsafety)
        !

        did_nmr = .false.

        !   Calculate the magnetic susceptability 
        !     using Equation 3 of Gregor et al
        !     - fastest for molecules

        if (iand(pw_params%output(1),8388608) == 8388608) then       
           call magnetic_suscept_eqn3(pw_params%chi_q,pw_params%g0mask, &
                pw_params,ffts,altbands,&
                pot_gspace,vion,crys,syms,pspot,blop,altkpoints)
           did_nmr  = .true.
        end if

        !   Calculate the magnetic susceptability 
        !     using Equation 8 of Gregor et al
        !     - second fastest for molecules

        if (iand(pw_params%output(1),16777216) == 16777216) then
           call magnetic_suscept_eqn8(pw_params%chi_q,pw_params%g0mask, &
                pw_params,ffts,altbands,&
                pot_gspace,vion,crys,syms,pspot,blop,altkpoints)
           did_nmr = .true.
        end if

        !   Calculate the magnetic susceptability 
        !     using the Mauri approach for periodic systems
        !     and pseudopotential corrections
        !     - slowest for molecules       

        if (iand(pw_params%output(1),33554432) == 33554432) then
           call magnetic_suscept_crys(pw_params%chi_q,pw_params%g0mask(1), &
                pw_params,ffts,altbands,&
                pot_gspace,vion,crys,syms,pspot,blop,altkpoints)
           did_nmr = .true.
        end if

        !   Calculate the magnetic susceptability 
        !     using the Mauri approach for periodic systems
        !     This is the original incorrect version
        !     - do not use!

        if (iand(pw_params%output(1),67108864) == 67108864) then
           call magnetic_suscept_orig(pw_params%chi_q,pw_params%g0mask(1), &
                pw_params,ffts,altbands,&
                pot_gspace,vion,crys,syms,pspot,altkpoints)
           did_nmr = .true.
        end if

        if(.not.did_nmr) then
           call magnetic_suscept_crys(pw_params%chi_q,pw_params%g0mask(1), &
                pw_params,ffts,altbands,&
                pot_gspace,vion,crys,syms,pspot,blop,altkpoints)
        end if

        !
        !  Deallocate all the arrays which are dynamically allocated
        !
        call destroy(ffts)
        call destroy(altkpoints)
        call destroy(altbands)
        call destroy(pot_gspace)
        call destroy(pspot)
        call destroy(blop)
        !     
        !  Deallocate all the gspace arrays
        !
        call destroy(vion)! ------------ endNMR   don't touch this line!!
        !
        ! beginplotNMR   !  don't touch this line!!          
        ! ------------ read and plotNMRshift  ----------------------------
        !
        !
!pub case(2);  call mystop( ' plot NMR shift option not availible' )
     case(2)  
        !
        !           generate large gspace for the potential
        !
        call generate_potential_gspace(pot_gspace, crys, syms, pw_params, &
             myproc, nproc)
        !
        !           allocate work arrays for fft
        !
        call create_fftstruc(ffts, pot_gspace, nfn)
        allocate(chi_g(pot_gspace%length, 3, 3))
        idspin = 1
        id9 = 9
        ierr = 0

        call readgsdat(1, ierr, pot_gspace, chi_g, 9, id9, idspin, 'CHI', 3)
        if (id9 /= 9 .or. idspin /= 1 .or. ierr /= 0) then
           write(9, 40) 
        else
           call nmr_shift(pw_params, crys, pot_gspace, chi_g(1, 1, 1))


           !
           ! calculate the induced current from CHI, plot out
           ! onto lines and/or planes
           !
           allocate(j_induced(pot_gspace%length, 3))
           call induced_current(j_induced(1, 1), chi_g(1, 1, 1), &
                pw_params%ket, pot_gspace, crys)
!           call dxplot_field(2, pw_params%output(1), j_induced(1, 1), &
!                pot_gspace, ffts, crys, 'J_INDUC')
    !
    ! do any line plots
    !
           if (pw_params%nlineplot >= 1) then
              if (iand(pw_params%output(2),2) == 2) then
       !
       !project the current in the direction of bra
       !
                 call lineplot(6,crys,1.0,j_induced(1,1),pot_gspace, &
                      pw_params%nlineplot,pw_params%lineplot, &
                      pw_params%bra,pw_params%ket,3,'clineplot',9)
              else
       !
       !plot the current as a vector
       !
                 call lineplot(2,crys,1.0,j_induced(1,1),pot_gspace, &
                      pw_params%nlineplot,pw_params%lineplot, &
                      pw_params%bra,pw_params%ket,3,'clineplot',9)
              endif
           endif

   ! 
   ! do any slice plots
   !
           if (pw_params%nsliceplot >= 1) then
              if (iand(pw_params%output(2),2) == 2 ) then
       !
       !project the current in the direction of bra
       !
                 call sliceplot(6,crys,1.0,j_induced(1,1),pot_gspace, &
                      pw_params%nsliceplot,pw_params%sliceplot, &
                      pw_params%bra,pw_params%ket,3,'clineplot')
              else
       !
       !plot the current as a vector
       !
                 
                 call sliceplot(2,crys,1.0,j_induced(1,1),pot_gspace, &
                      pw_params%nsliceplot,pw_params%sliceplot, &
                      pw_params%bra,pw_params%ket,3,'clineplot')
              endif
           endif
           deallocate(j_induced)
                                            
!
!   this reads the direct current and plots it jry
!
           if (iand(pw_params%output(2),1) == 1) then
              Call direct_current_process(crys, ffts,pot_gspace,pw_params,&
                   pw_params%nlineplot,pw_params%lineplot, &
                   pw_params%nsliceplot,pw_params%sliceplot, &
                   pw_params%bra,pw_params%ket,pw_params%chi_q)
           endif
 
           !various routines to plot properties of CHI
           !uncomment out what you need               
           
           !
           !   ... and plot if desired
           !
           if (iand(pw_params%output(1), 4) == 4) then
           allocate(chi_trace(pot_gspace%length))
           ddum = dmone / (dthree * 137.036**2) * 1.0d6
           chi_trace(:) = ddum * (chi_g(:, 1, 1) + chi_g(:, 2, 2) + &
                chi_g(:, 3, 3))
           call dxplot(1, pw_params%output(1), chi_trace, pot_gspace, ffts,&
                crys, 'MCHI_TR',7)
           call lineplot(1, crys, done, chi_trace(1), pot_gspace, &
                pw_params%nlineplot, pw_params%lineplot, pw_params%bra, &
                pw_params%ket, 1, 'MCHI_TR',7)
           deallocate(chi_trace)
        end if
        
        !
        !              plot certain component of CHI along a line or plane
        !
        !        call lineplot(1, crys, dmone / (137.036**2) * 1.0d6, &
        !        chi_g(1, 1, 1), pot_gspace, nlineplot, lineplot, bra, ket, &
        !                9, 'MCHI___')
        !
        !                                    
        !
        !     $              (chi_g(:,1,1)+chi_g(:,2,2)+chi_g(:,3,3))/3.d0
        !               call lineplot(1, crys, 1.d0/(137.036**2.d0)*1.d6,
        !     $              chi_trace(1), pot_gspace,  nlineplot,
        !     $              lineplot, bra,ket, 1, 'MCHI_TR')
        !               call dxplot_field(1,pw_params%output(1),chi_g(1,1,1),
        !     $              pot_gspace,ffts,   crys, 'B_INDUC')
        end if
        !                                   
        !
        !           deallocate all the arrays which have dynamical allocation
        !
        deallocate(chi_g)
        call destroy(ffts)
        call destroy(pot_gspace)!endplotNMR   - don't touch this line!!
        !            
        !  ------------ compute band structure  ----------------------------
        !
     case(3)  
        call generate_potential_gspace(pot_gspace, &
             crys, syms,pw_params, myproc, nproc)
        call create_fftstruc(ffts, pot_gspace, nfn)
        call read_pseudo(1, crys%nspin, pw_params, pot_gspace, &
        energs, pspot, crys, vion, denc, denv, ddc, dvql, vql,dnc,vnl_ave,ffts)
!       iscr2 was 4 for bandstructure,  hjchoi
        iscr2 = pw_params%iscr

        if (iscr2 .eq. 5) iscr2 = 6   ! hjchoi
        call screen(iscr2, .false., pot_gspace, crys, ffts, energs, pw_params,&
             vion, denv, denc, vin, dena, chdr)
        call create_gspace_array(vout, pot_gspace%length, 1, crys%nspin)

        logidum = altkpoints%generated
        altkpoints%generated = .false.
        call generate_kpoints(1, altkpoints, crys, syms, &
             altbands, pw_params%output(1), pw_params%q, nsafety)
        altkpoints%generated = logidum
       
        allocate(k_gspace(altkpoints%nrk))
        call generate_k_gspaces(k_gspace, maxlength, pot_gspace, &
             altkpoints,  crys, syms, pw_params%emax, pw_params%output(1))
        !
        !  Distinguish between the case when the wavefunctions
        !  are needed to do further calculations like in the 
        !  momemtum distribution, or when we do not need them. 
        !  In this later case, only allocate the memory of the 
        !  wavefunction at one k-point and reuse it.
        !

!. P. Zhang, I need wavefunctions for LSDA+U

!        if((.not.(iand(pw_params%output(1),4096).eq.4096)) & ! momentum density
!             .and.(.not.(iand(pw_params%output(1),4194304) &
!             .eq.4194304))) then ! projectors

        if((.not.(iand(pw_params%output(1),4096).eq.4096)) & ! momentum density
             .and.(.not.(iand(pw_params%output(1),4194304) &
             .eq.4194304)).and.(.not.(pw_params%lsdau.eq.1))) then ! projectors

           call create_gspace_array(wavefnbs, maxlength * altbands%max, &
                1,crys%nspin)

           i=pw_params%maxitscf
           pw_params%maxitscf=1

           ! call with imode=2
           call scfloop(t0,2, pot_gspace, k_gspace(1), energs, pw_params, & 
                crys, syms, ffts, pspot, altkpoints, altbands, wavefnbs, &
                .false.,vion, vin, vout, dena, denc, chdr,vnl_ave)
           pw_params%maxitscf = i

        else

           call create_gspace_array(wavefnbs, maxlength * altbands%max, &
                altkpoints%nrk, crys%nspin)

           i = pw_params%maxitscf
           pw_params%maxitscf = 1
           call scfloop(t0,0, pot_gspace, k_gspace(1), energs, pw_params,  &
             crys, syms, ffts, pspot, altkpoints, altbands, wavefnbs,  &
             .false., vion, vin, vout, dena, denc, chdr, vnl_ave)
           pw_params%maxitscf = i

        end if
       
        if (myproc == 0) call plot_bandstruc(altkpoints, altbands, &
                                     crys, energs, pw_params%bslabels)

        ! output wavefunction in G space in case needed

        if (iand(pw_params%output(1), 524288) == 524288) &
             call gwout('GWR',crys, syms, altbands, pot_gspace, k_gspace, &
             wavefnbs, altkpoints, 1)


        if (iand(pw_params%output(1), 1048576) == 1048576) &
             call gwout('GWC', crys, syms, altbands, pot_gspace, k_gspace, &
             wavefnbs, altkpoints,2)
        !
        if (iand(pw_params%output(1), 4096) == 4096) then
           call momentum_density(crys,altkpoints,altbands,k_gspace,wavefnbs)
        end if

        !
        !           plot wave functions
        !


        if (iand(pw_params%output(1),256) == 256) then
           do is = 1, crys%nspin
              do j = 1, pw_params%nplotwv
                 irk = pw_params%plotwv(1+ (j-1)*2) ! get altkpoint number
                 if (irk > altkpoints%nrk) then
                    write(9, 530) 'altkpoint', irk, &
                         'number of kpoints', altkpoints%nrk
                    irk = altkpoints%nrk
                 end if
                 n = pw_params%plotwv(2+ (j-1)*2) ! get band number
                 if (n > altbands%nband(irk, is)) then
                    write (9,530) 'band', n, 'number of bands', &
                         altbands%nband(irk, is)
                    n = altbands%nband(irk, is)
                 end if
                 ! Write the square of modulus of the wavefunctions
                 write(strbuf1, '(a3,i2.2,a1,i3.3,a1,i3.3)') 'WF.',is,'.', irk, '.', n
                 call dxplot(3, pw_params%output(1), &
                      wavefnbs%data(k_gspace(irk)%length * &
                      (n - 1) + 1, irk, is), k_gspace(irk), &
                      ffts, crys, strbuf1,13)
                 ! Write the real part of the wavefunctions
                 write(strbuf1, '(a3,i2.2,a1,i3.3,a1,i3.3)') 'WR.',is,'.', irk, '.', n
                 call dxplot(1,pw_params%output(1), &
                      wavefnbs%data(k_gspace(irk)%length * &
                      (n - 1) + 1, irk, is), k_gspace(irk), &
                      ffts, crys, strbuf1,13)
                 ! Write the imaginary part of the wavefunctions
                 write(strbuf1, '(a3,i2.2,a1,i3.3,a1,i3.3)') 'WI.',is,'.', irk, '.', n
                 call dxplot(2, pw_params%output(1), &
                      wavefnbs%data(k_gspace(irk)%length * &
                      (n - 1) + 1, irk, is), k_gspace(irk), &
                      ffts, crys, strbuf1,13)
                 call lineplot(1, crys, done, &
                      wavefnbs%data(k_gspace(irk)%length * (n - 1) + 1, &
                      irk, is), k_gspace(irk), pw_params%nlineplot, &
                      pw_params%lineplot, pw_params%bra, pw_params%ket, 1,&
                      strbuf1,13)
              end do
           end do
        end if

        if (iand(pw_params%output(1), 2048) == 2048) then
           call write_eigvals(myproc, pw_params, crys, &
                energs, altkpoints, altbands, p_weights) 
        end if

!----------------------------------------------------------------------
! Peihong's version of angular momentum projection

!        if (iand(pw_params%output(1), 4194304) == 4194304) &
!             call project(myproc, pw_params%emax, altbands%min, crys, &
!             wavefnbs, altkpoints%nrk, k_gspace(1), 'project.bs.dat')


        if (iand(pw_params%output(1), 4194304) == 4194304) then
           call project_pz(myproc,nproc, crys, syms,wavefnbs, altkpoints, k_gspace,altbands,energs,0)
        end if
!----------------------------------------------------------------------


        do irk = 1, altkpoints%nrk
           call destroy(k_gspace(irk))
        end do
        deallocate(k_gspace)

        call destroy(ffts)
        call destroy(altkpoints)
        call destroy(altbands)
        call destroy(pot_gspace)
        call destroy(vion)
        call destroy(denc)
        call destroy(denv)
        call destroy(ddc)
        call destroy(dvql)
        call destroy(vql)
        call destroy(dnc)
        call destroy(vin)
        call destroy(dena)
        call destroy(chdr)
        call destroy(vout)
        call destroy(wavefnbs)
        call destroy(pspot)
        !            
        ! ------------ plot the potential and the charge ---------------
        !
        !
        !
     case(4)  
        !
        !           generate large gspace for the potential
        !
        call generate_potential_gspace(pot_gspace, crys, syms, pw_params, &
             myproc, nproc)
        !
        !           allocate work arrays for fft
        !
        call create_fftstruc(ffts, pot_gspace, nfn)
        !
        !           read pseudopotential and set up valence/core charge
        !           (the ddc,dvql, and dnc arrays are for force/stress calc)
        !
        call read_pseudo(1, crys%nspin, pw_params, pot_gspace, energs, &
             pspot, crys, vion, denc, denv, ddc, dvql, vql, dnc,vnl_ave,ffts)
        !
        !   screen the ionic potential with the valence and core charge density
        !   also creates arrays vin, den, chdr.
        call screen(4, .false., pot_gspace, crys, ffts, energs, pw_params, &
             vion, denv, denc, vin, dena, chdr)
        !
        !           do lineplot and printout for the data explorer
        !
        if (iand(pw_params%output(1), 262144) == 262144) then
           if (crys%nspin == 1) then
             call lineplot(1, crys, done, vion%data(1, 1, 1), &
                pot_gspace, pw_params%nlineplot, pw_params%lineplot, &
                 pw_params%bra, pw_params%ket, 1, 'VEFFAVG',7)
             call dxplot(1, pw_params%output(1), vion%data(1, 1, 1), &
                pot_gspace, ffts, crys, 'VEFFAVG',7)
           else
             call lineplot(1, crys, done, vion%data(1, 1, 1), &
                pot_gspace, pw_params%nlineplot, pw_params%lineplot, &
                 pw_params%bra, pw_params%ket, 1, 'V_UP',4)
             call dxplot(1, pw_params%output(1), vion%data(1, 1, 1), &
                pot_gspace, ffts, crys, 'V_UP',4)
             call lineplot(1, crys, done, vion%data(1, 1, 2), &
                pot_gspace, pw_params%nlineplot, pw_params%lineplot, &
                 pw_params%bra, pw_params%ket, 1, 'V_DOWN',6)
             call dxplot(1, pw_params%output(1), vion%data(1, 1, 2), &
                pot_gspace, ffts, crys, 'V_DOWN',6)             
           end if

        end if
        !
        !           plot the charge density if the flag is set
        !
        if (iand(pw_params%output(1), 16) == 16) then
           if (crys%nspin == 1) then
              call dxplot(1, pw_params%output(1), dena%data(1, 1, 1), &
                   pot_gspace, ffts, crys, 'CHARGE_',7)
              call lineplot(1, crys, done / crys%vcell, dena%data(1, 1, 1), &
                   pot_gspace, pw_params%nlineplot, pw_params%lineplot, &
                   pw_params%bra, pw_params%ket, 1,'CHARGE_',7)
           else
              call dxplot(1, pw_params%output(1), dena%data(1, 1, 1), &
                   pot_gspace, ffts, crys, 'CHAR_UP',7)
              call dxplot(1, pw_params%output(1), dena%data(1, 1, 2), &
                   pot_gspace, ffts, crys, 'CHAR_DN',7)
              call lineplot(1, crys, done, dena%data(1, 1, 1), & 
                   pot_gspace, pw_params%nlineplot, pw_params%lineplot,&
                   pw_params%bra, pw_params%ket, 1,&
                   'CHAR_UP',7)
              call lineplot(1, crys, done, dena%data(1, 1, 2), & 
                   pot_gspace, pw_params%nlineplot, pw_params%lineplot,&
                   pw_params%bra, pw_params%ket, 1, 'CHAR_DN',7)
           end if
        end if
        !
        !           compute and plot electric field
        !
        if (iand(pw_params%output(1), 131072) == 131072) then
           call electric_field(pw_params, crys, dena%data(1, 1, 1), &
                vion%data(1, 1, 1), vin%data(1, 1, 1), pot_gspace, &
                pw_params%nlineplot, pw_params%lineplot, pw_params%bra, &
                pw_params%ket)
        end if

        call destroy(ffts)
        call destroy(pot_gspace)
        call destroy(pspot)
        !
        !  begintddft  - don't touch this line!!          
        !     ------ timedependent density-functional theory (TDDFT) -----
        !
        !
        !
!pub case(5);  call mystop( ' TDDFT option not availible' )

     case(5)
        !
        !    first calculate new bands needed (similar to band structure calc.)
        !
        call generate_potential_gspace(pot_gspace, crys, syms, pw_params, &
             myproc, nproc)
        call create_fftstruc(ffts, pot_gspace, nfn)
        call read_pseudo(1, crys%nspin, pw_params, pot_gspace, &
       energs, pspot, crys, vion, denc, denv, ddc, dvql, vql, dnc,vnl_ave,ffts)
        call screen(4, .false., pot_gspace, crys, ffts, energs, pw_params, &
             vion, denv, denc, vin, dena, chdr)
        call create_gspace_array(vout, pot_gspace%length, 1, crys%nspin)
        !
        !    generate symmetry structure corresponding to identity only
        !

       id_sym%tnp = dzero
        id_sym%rsymmat = dzero
        id_sym%mtrx = 0
        id_sym%rmtrx = 0
        do i = 1, 3
           id_sym%rsymmat(i, i, 1) = done
           id_sym%mtrx(i, i, 1) = 1
           id_sym%rmtrx(i, i, 1) = 1
        end do
        id_sym%ntrans = 1
        !
        !    generate k-points, reduced by time-inversion symmetry only
        !
        altbands%min(1) = pw_params%nband_tddft
        call generate_kpoints(1, kpoints, crys, id_sym, &
             altbands, pw_params%output(1), pw_params%q, nsafety)

        allocate(k_gspace(kpoints%nrk))
        call generate_k_gspaces(k_gspace, maxlength, &
      pot_gspace, kpoints,  crys, syms, pw_params%emax, pw_params%output(1))
        call create_gspace_array(wavefn, maxlength * altbands%max, &
             kpoints%nrk,crys%nspin)

        i = pw_params%maxitscf
        pw_params%maxitscf = 1
        energs%inputfermi = energs%efermi(1)
        call scfloop(t0,0, pot_gspace, k_gspace(1), energs, pw_params, crys, &
             syms, ffts, pspot, kpoints, altbands, wavefn, .false., &
             vion, vin, vout, dena, denc, chdr, vnl_ave)
        pw_params%maxitscf = i
        !
        !        write eigenvectors to disk 
        !
        if (iand(pw_params%output(1), 512) == 512) then
           write(9, 130)
           if (pot_gspace%myproc == 0) then
              call mysystem('rm -f BAND')
              call write_band('BAND', altbands)
           end if
          do is = 1, crys%nspin
              do irk = 1, kpoints%nrk
                 write(wfnname, 120) irk, is
                 call writegsdat(0, k_gspace(irk), wavefn%data(1, irk, is), &
                      altbands%nband(irk, is), altbands%nband(irk, is), 1, &
                      wfnname, 10)
              end do
           end do
        end if

        if (iand(pw_params%output(1), 2048) == 2048) then
           call write_eigvals(myproc, pw_params, crys, energs, kpoints, &
                altbands, p_weights) 
        end if



        !
        !   now do TDDFT
        !
        call tddft(1, pot_gspace, k_gspace(1), pw_params, crys, &
             syms, ffts, kpoints, altbands, wavefn, denc, dena)
        !
        !   annihilate everything
        !
        do irk = 1, altkpoints%nrk
           call destroy(k_gspace(irk))
        end do
        deallocate(k_gspace)
        call destroy(ffts)
        call destroy(altkpoints)
        call destroy(altbands)
        call destroy(pot_gspace)
        call destroy(vion)
        call destroy(denc)
        call destroy(denv)
        call destroy(ddc)
        call destroy(dvql)
        call destroy(vql)
        call destroy(dnc)
        call destroy(vin)
        call destroy(dena)
        call destroy(chdr)
        call destroy(vout)
        call destroy(wavefn)
        call destroy(pspot)   !endtddft  - don't touch this line !!
        !
        !
        ! ------------------------------------------------------------------
        !
     case default
        write(0, *) 'ERROR: job ', pw_params%joblist(ijob), ' not implemented yet!'
        call mystop
     end select
     write(9, 300) gimmetime() - tstart
     write(9,301)   gimmetime() - t0
     call myflush(9)
  end do                    ! end of gigantic loop 

  enew = energs%total

!  call destroy(crys)

  call deallocate_shmem_work()
  
  call blacs_gridexit(blacs_context)
  
  close(9)
!  close(34)

  return
  !

20 format('CDW',i1,i5.5)
30 format(//' ================ JOB: ',a,' ====================='//)

40 format(//' *** MAGNETIC SUSCEPTIBILITY FILE IS CORRUPT!')
110 format(/' FIXING FERROMAGNETIC MOMENT TO BE ZERO')

120 format('WFN',i5.5,'.',i1)
130 format(/' <<<<<<<< WRITING WAVEFUNCTIONS TO DISK >>>>>>>>>>')


150  format(/'Static   : ',F10.3,' Mb')
151  format( 'Gs-fft   : ',F10.3,' Mb')
152  format( 'Bands    : ',F10.3,' Mb')
153  format( 'Pseudo   : ',F10.3,' Mb')
154  format( 'CG solve : ',F10.3,' Mb')
155  format( '========================')
156  format( 'Total    : ',F10.3,' Mb')

200 format(/,'  COMPUTING TIME FOR ITERATION ',I5,3X,F10.2)
250 format(/,'  COMPUTING TIME FOR STARTING ',3X,F12.6,/)
260 format(/,'  TOTAL TIME s ',3X,F12.6,/)

300 format(/'  TOTAL COMPUTING TIME FOR JOB: ',F10.2)
301 format(/'  TOTAL ELAPSED COMPUTING TIME: ',F10.2)
400 format('  E= ',F12.6)
500 format(/,20x,'START OF SELFCONSISTENT CYCLE',/1x,78('-'))
510 format(' SPIN(',i1,'):')
520 format(' COMPUTED K-POINT: ', 3f12.6,' TIMEs:',f12.6) 
525 format(' diagonalization did not converge in ',i3, &
       ' iterations:')

530 format(' *** WARNING: ',a6,' NUMBER=',i3, &
       ' IN PLOTWV EXCEEDS ',a,' = ',i3)
540 format(' *** REDUCED ',a6,' NUMBER TO ',a6,' = ',i5)

600 format(' *** WARNING: ERROR WHEN READING OLD EIGENVECTORS',&
       /' ***          GENERATING STARTING GUESS')   
605 format(' *** WARNING: ERROR WHEN READING NCALC OR ESHIFT FILE',&
       /' ***          GENERATING STARTING GUESS')   
610 format(' NUMBER OF PROCESSORS: ',i4)

620 format(' Spawning information:',&
       /' -------------------- ',&
       /' Program to spawn: ',a,&
       /' Scratch directory: ',a,&
       /' Group name:',a,&
       /' Number of processors:',i4)

end subroutine pwmain
