!-*-Fortran-*-
!     @process extchk
!
subroutine direct_emin(t0,pot_gspace, k_gspace, energs, pw_params, &
     crys, syms, ffts, pspot, kpoints, bands, wavefn, lrwfn, vion, &
     vhxc, vout, den, denc, rden)

  include 'use.h'  
  implicit none  
  include 'interface.h'  
  !
  !     INPUT:
  !     -----
  !
  real(dp)        t0                      ! starting time
  type(symmetry), intent(in) :: syms                    ! symmetry operations
  type(crystal), intent(in) :: crys                       ! crystal structure
  type(fft_struc), intent(in) :: ffts               ! information for the FFT
  type(pw_parameter), intent(inout) :: pw_params      ! plane wave parameters
  type(pseudo_potential), intent(in) :: pspot    ! info about the pseudopot'l
  type(kpoint), intent(in) :: kpoints                ! BZ integration kpoints
  type(parallel_gspace), intent(in) :: pot_gspace, k_gspace(*)  
  complex(dp), intent(in) :: &
       denc(pot_gspace%length)                          ! core charge density
  logical, intent(in) :: lrwfn        ! true if wavefn is to be reused, false
  !
  !     INPUT/OUTPUT:
  !     -------------
  !
  type(energy), intent(inout) :: energs      ! energies and other information
  type(band), intent(inout) :: bands    ! eigenvalues, occupation numbers etc
  complex(dp), intent(inout) :: &
       vion(pot_gspace%length, crys%nspin), &       !screened ionic potential
       vhxc(pot_gspace%length, crys%nspin)                 ! V_Hartree + V_xc
  !
  !     OUTPUT:
  !     ------
  !
  complex(dp), intent(out) :: &
       den(pot_gspace%length, crys%nspin), &         ! valence charge density
       vout(pot_gspace%length, crys%nspin)                 ! V_Hartree + V_xc
  type(complex_gspace_array), intent(inout) :: &
       wavefn                         ! wave functions for all kpoints, bands
  !
  !     WORK:
  !     ----
  !
  complex(dp) :: rden(pot_gspace%r_size)  
  !
  !     DESCRIPTION:
  !     -----------
  !
  !     Does a direct minimization of the energy functional. The system
  !     must be an insulator for this to work. Also, the memory costs
  !     are very high when multiple k-points are used.
  !
  !
  !     --------------- local variables -------------------------
  !
  type(parallel_gspace) :: sub_gspace                ! for the starting guess
  complex(dp), pointer :: vloc(:),vloc2(:)              ! local potential in realspace
  type(hamiltonian) :: ham(kpoints%nrk,crys%nspin)! hamiltonian at all kpoints
  type(dense_hamiltonian) :: densham                   ! subspace hamiltonian
  type(energy) :: oldenergs                           ! as a dummy for etotal
  complex(dp), allocatable :: vbareion(:,:)              ! bare ionic potential
  real(dp) :: dsmear, sqrpi, delta,gimmetime, &
       ddum                                                ! dummy for etotal
  real(dp) time(14),t_temp,qk(3),qcar(3)
  integer :: i, j, n, &
       nsp, nvecs, &
       ierr, &                                                   ! error flag
       is, &                                                   ! spin counter
       irk,&                                                   ! kpoint counter
       iranm
  logical :: iconv                   ! indicates if minimization is converged
  character(len=11) :: wfnname  
  real Etot_corr_old,Eband_diff,Etot_corr,ELSDA_corr

! LDA+U direction minimization not implemented
!
  Etot_corr_old=0.d0
  ELSDA_corr=0.d0
  Eband_diff=0.d0

  !
  !     --------------------------------------------------------------
  !
  t_temp=gimmetime(); time=dzero;
  write(9, 100) ; call myflush(9)  
  dsmear = pw_params%smearing / ryd  
  sqrpi = sqrt(pi)  
  allocate(vloc(ffts%r_size))                  ! local potential in realspace
  if (crys%nspin .eq. 2)  allocate(vloc2(ffts%r_size))
  !
  !     Set up the Hamiltonians at all k-points. Only one spin
  !     component is set up, and later the local potential component of
  !     the Hamiltonians is set to either spin up or spin down local
  !     potential. This is to save memory. After all, the Kleinman
  !     -Bylander projectors are the same for spin up and spin down
  !
!  do irk = 1, kpoints%nrk  
!     call create_hamiltonian(ham(irk), k_gspace(irk), pspot, vloc, ffts)
!     ham(irk)%ekinmod = pw_params%ekinmod  
!     call setup_nonlocal_potential(ham(irk), pspot, crys)  
!  end do
  !
  !     ------ Now read or generate the starting guesses ----------
  !                  also set up the local potential
  !
  do is = 1, crys%nspin  

   if (is .eq. 1) then
     call setup_local_potential(0, ffts, pot_gspace, vion(1, is), vloc(1))
   else
     call setup_local_potential(0, ffts, pot_gspace, vion(1, is), vloc2(1))
   end if

   do irk = 1, kpoints%nrk

    if (is .eq. 1) then
      call create_hamiltonian(ham(irk,is), k_gspace(irk), pspot, vloc, ffts)
    else
      call create_hamiltonian(ham(irk,is), k_gspace(irk), pspot, vloc2, ffts)
    end if
    ham(irk,is)%ekinmod = pw_params%ekinmod  
    !
    !  if no submatrix_diag then set iter(last argument in setup of 
    !  non-local potential) equal to 2 and this will have the G-space
    !  projectors not to be calculated if real space projectors are used
    !
    if (.not. iand(pw_params%miscflag, 4) == 4) then  
      call setup_nonlocal_potential(ham(irk,is), pspot, crys,1)
    else
      call setup_nonlocal_potential(ham(irk,is), pspot, crys,2)
    end if 
    ierr = 1  
    !
    !  if not re-using wavefunctions
    !
    if (.not.lrwfn) then  
      !
      !  read from file
      !
      if (iand(pw_params%input, 1) == 1) then          
        call read_band('BAND', bands)  
        if (bands%nrk.le.0) then                                ! error
          write(9, 155)  
          bands%nrk = kpoints%nrk
          goto 20  
        end if
        write(wfnname, 120) irk, is  
        wavefn%data(:, irk, is) = zzero  
        nsp = 1
        nvecs = bands%nband(irk, is)  
        call readgsdat(1, ierr, k_gspace(irk), &
                   wavefn%data(1, irk, is), bands%nband(irk, is), nvecs, &
                   nsp, wfnname, 10)
        if (ierr /= 0) then  
          write(9, 150) wfnname
          goto 20  
        end if
        write(9, 140) nvecs, wfnname  
        ierr = 0  
20      continue  
      end if
      !
      !  use random wavefunctions with no submatrix diagnalization
      !
      if (ierr == 1  .and.  iand(pw_params%miscflag, 4) == 4 ) then  

       iranm=-4513*(is-1)-1371*irk-5616*(pot_gspace%myproc+1)
       call init_rdm_wf(wavefn%data(1, irk, is) ,iranm, bands%min(is),ffts,&
           k_gspace(irk),crys)          

      else if (ierr == 1) then         ! ------- have to do submatrix diag --
              sub_gspace%gmax = sqrt(pw_params%emaxsub) 
              sub_gspace%rk = ham (irk,is)%gspace%rk    ! shift from big brother
              sub_gspace%nproc = 1            ! scatter across all processors
              sub_gspace%name = 'sub'
              sub_gspace%myproc = 0 
              sub_gspace%fftsize = (/ 0, 0, 0, 0, 0 /)  
              sub_gspace%istar = .false.
              sub_gspace%igvec = .true.  
              sub_gspace%iekin = .true.
              sub_gspace%iinverse = .false.  
              sub_gspace%istruc = .false.  
              sub_gspace%ipackinfo = .false.  
!              call generate_gspace(pw_params%output(1), sub_gspace, crys, syms, 1)
              call generate_gspace(pw_params%output(1), sub_gspace, crys, syms)
              if(is.eq.1.and.irk.eq.1)call print_mem_use(pot_gspace,ffts,crys,&
                    k_gspace,pspot,pw_params,sub_gspace,ham,kpoints,bands)
              call create_dense_hamiltonian(densham, sub_gspace)  
              call start_vectors(1, pw_params, ham(irk,is), densham, &
                   bands%energy(1, irk, is), bands%max, &
                   wavefn%data(1, irk, is), pot_gspace, sub_gspace, &
                   vion(1, is), crys)
              call destroy(densham)  
              call destroy(sub_gspace)  
              if (pspot%NL_rspace(1)) &
                      call destroy_complex_gspace_array(ham(irk,is)%vnloc)  
           end if
        end if
        if (pw_params%bandgap > dzero) then          ! put in bandgap by hand
           bands%energy(bands%min(is) + 1, irk, is) = &
                bands%energy(bands%min(is), irk, is) + pw_params%bandgap
        end if
        bands%nband(irk, is) = bands%min(is)            ! it is an insulator
     end do
  end do
  !
  !     ------- Do the direct minimization ----------------------------
  !
  if (lrwfn) write(9, 160)  
  pw_params%isymmetrize = .false.                        ! need to symmetrize
  if ((kpoints%nrk > 1) .or. (kpoints%shift(1) /= dzero) .or. &
       (kpoints%shift(2) /= dzero) .or. (kpoints%shift(3) /= dzero)) then
     if (syms%ntrans > 1) then  
        pw_params%isymmetrize = .true.  
     end if
  end if
!  if (crys%nspin > 1 .and. iand(pw_params%optimize,8) == 8 ) then  
!     write(9,*) '*** ENERGY MINIMIZATION WORKS ONLY FOR NSPIN=1 if OPTIMIZE METAL'  
!     call mystop  
!  end if
  allocate(vbareion(pot_gspace%length, crys%nspin ) )    ! bare ionic potential
  do is=1, crys%nspin
   vbareion(:,is) = vion(:, is) - vhxc(:, is)  
   if (pot_gspace%ig0 > 0) vbareion(pot_gspace%ig0,is) = zzero
  end do

     write(9,*)  gimmetime()-t0,' CGEMIN TIME' 

    if (iand(pw_params%output(1), 8) == 8) call get_timing(time(8),t_temp)
  if (iand(pw_params%optimize,64) == 64 ) then
    call cgemin_metal_gcg(time,t0,pw_params, pw_params%maxitdiag, iconv, bands, &
      pw_params%epsdiag, kpoints, ham(1,1), wavefn, pot_gspace, crys, &
      energs, denc(1), rden(1), den(1, 1), vhxc(1, 1), vbareion(1,1))
  else

    energs%esmear= dzero

    if(pw_params%diag_method == 'Mauri') then
      call cgemin(time,t0,pw_params, pw_params%maxitdiag, iconv, bands, &
       pw_params%epsdiag, kpoints, ham(1,1), wavefn, pot_gspace, crys, &
       energs, denc(1), rden(1), den(1, 1), vhxc(1, 1), vbareion(1,1))
    else if (pw_params%diag_method == 'Grassmann') then
  
       call cgemin_gcg(time,t0,pw_params, pw_params%maxitdiag, iconv, bands, &
       pw_params%epsdiag, kpoints, ham(1,1), wavefn, pot_gspace, crys, &
       energs, denc(1), rden(1), den(1, 1), vhxc(1, 1), vbareion(1,1))
    else 
     call mystop( ' direct_emin' )
    end if 
  end if
  if (iand(pw_params%output(1), 8) == 8)  t_temp=gimmetime()
  !
  !     ------- compute and print total energy ----------------------
  !
  !
  !     compute kinetic and band energy and set band occupation numbers
  !
  energs%ektot = dzero
  energs%eband = dzero
  energs%esmear= dzero

  if (iand(pw_params%optimize,64) == 64 ) then
    call flevel(1, 2, pw_params, energs, crys, bands, kpoints)  
  else
   do is=1,crys%nspin 
    do irk = 1, kpoints%nrk 
      bands%occup(1:bands%min(is),irk,is) = dtwo*kpoints%w(irk)
    end do
    bands%ifmax(:,is)=bands%min(is)
   end do
  end if


  do is = 1, crys%nspin  
    do irk = 1, kpoints%nrk  
      call kinet(ham(irk,is), wavefn%data(1, irk, is), bands%nband(irk, is), &
             bands%ekn(1, irk, is))
      do n = 1, bands%nband(irk, is)  
         energs%ektot = energs%ektot + bands%occup(n, irk, is) * &
                bands%ekn(n, irk, is)
         energs%eband = energs%eband+bands%occup(n, irk, is) * &
                bands%energy(n, irk, is)
        end do
     end do
  end do
  energs%eband = energs%eband / real(bands%nspin, dp)
  energs%ektot = energs%ektot / real(bands%nspin, dp)

  write(9, * ) 'vxc0=', energs%vxc0(1)  
  write(9, 110)  
  do is = 1, crys%nspin  
   do irk = 1, kpoints%nrk  
     write(9, 115) irk, bands%energy(1:bands%min(is), irk, 1)  
   end do
  end do
  !
  !     prepare potentials for output
  !
  oldenergs = energs  
  do is = 1, crys%nspin  
     vion(:, is) = vhxc(:, is)  
     vhxc(:, is) = vhxc(:, is) - vbareion(:,is)  
     vout(:, is) = vhxc(:, is)  
     if (pot_gspace%ig0 > 0) vout(pot_gspace%ig0, is) = energs%vxc0(is)
  end do

  call etotal(crys, 1, iconv, ddum, energs, oldenergs, pw_params, &
       pot_gspace, vion(1, 1), vhxc(1, 1), vhxc(1, 1), den(1, 1), &
       Etot_corr,ELSDA_corr,Etot_corr_old,Eband_diff)
  !
  !     ------- Write charge density to disk -------------------------
  !
  call writegsdat(0, pot_gspace, den(1, 1), 1, 1, crys%nspin, 'CD', 2)
  !
  !     ------- Write eigenvectors to disk if requested ---------------
  !
  if (iand(pw_params%output(1), 512) == 512) then  
     write(9, 130)  
     do is = 1, crys%nspin  
        do irk = 1, kpoints%nrk  
           write(wfnname, 120) irk, is  
           call writegsdat(0, k_gspace(irk), wavefn%data(1, irk, is), &
                bands%nband(irk, is), bands%nband(irk, is), 1, wfnname, 10)
        end do
     end do
     if (pot_gspace%myproc == 0) then  
        call mysystem ('rm -f BAND')  
        call write_band ('BAND', bands)  
     end if
  end if
  !
  !     ----------- deallocate all memory -----------------------
  !
  do is = 1, crys%nspin 
  do irk = 1, kpoints%nrk  
     call destroy_hamiltonian(ham(irk,is))  
  end do
  end do
  deallocate(vloc)
  if (crys%nspin .eq. 2)  deallocate(vloc2)  
  deallocate(vbareion)

    if (iand(pw_params%output(1), 8) == 8) call get_timing(time(5),t_temp)
    if (iand(pw_params%output(1), 8) == 8)&
     write(9,200) time(1),time(2),time(3),time(4),time(5),time(6),time(7),&
                  time(8), time(1)+time(2)+time(3)+time(4)+time(5)+time(6)&
                  +time(7)+time(8) 
            
  return

100 format(/,17('-'),' DIRECT MINIMIZATION OF ENERGY FUNCTIONAL ', &
       &     17('-'))
110 format(/' FINAL EIGENVALUES (WITHOUT VXC0):')  

115 format(/' KPOINT:',i4,/1000(5f12.6/))  
120 format('WFN',i5.5,'.',i1)  
130 format(/' <<<<<<<< WRITING WAVEFUNCTIONS TO DISK >>>>>>>>>>')  
140 format(' <<<<<<<< READ ',i4,' BANDS FROM FILE ',a,'>>>>>>>>')  
150 format(' *** COULD NOT FIND WAVE FUNCTION ',a,'ON DISK!')  
155 format(' *** COULD NOT FIND BAND FILE ON DISK!')  
160 format(' <<<<<<<< REUSING WAVE FUNCTIONS >>>>>>>>')  
200 format(/' DIAGONALIZATION TIME -  BREAKDOWN:',/1x,27('-'), &
       &            /20x,'ACTIONS',15x,'seconds'&
       &     /,11X,                             &
       &     /,11X,'FFT time               =',2X,F12.4, &
       &     /,11X,'NL time                =',2X,F12.4, &
       &     /,11X,' occ space time        =',2X,F12.4, &
       &     /,11X,' other cg time         =',2X,F12.4, &
       &     /,11X,' other scf time        =',2X,F12.4, &
       &     /,11X,' overlap cg time       =',2X,F12.4, &
       &     /,11X,' mzgemm m*m*len        =',2X,F12.4, &
       &     /,11X,' init diag             =',2X,F12.4, &
       &     /,11X,' total                 =',2X,F12.4)
end subroutine direct_emin
