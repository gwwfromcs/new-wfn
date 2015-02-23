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

! Main program of Paratec. It handles the input from the file "input" and
!  branches to do a relaxation of the atomic positions and of the unit cell 
!  shape or branhes to do a molecular dynamics run.
! 
!  originally relax.c written by B. Pfrommer.
! 
!  Changed over to Fortran90 by Peter Haynes March 2000
! 
!  modified extensively by adding documentation and organizing variables 
!  also made the relaxation of the atoms a subroutine and renamed to main.f90p
!   D. Raczkowski 2001
!
!------------------------------------------------------------------------

program main

  use constants
  use pw_parameter_module
  use kpoint_module
  use band_module
  use symmetry_module
  use energy_module
  use crystal_module
  use esdf
  use all_to_all_module
  use molecular_dynamics_module
  implicit none
  include 'flibcalls.ph' 
  include 'all_to_all.h'
include 'mpif.h'
  !
  !     DEFINED TYPES: - located in structures.f90
  !     -------------
  !
  type(pw_parameter) :: pw_params ! plane wave parameters
  type(kpoint) :: altkpoints      ! alternative BZ integration kpoints
  type(band) :: bands             ! eigenvalues, occupation numbers etc
  type(band) :: altbands          ! another set of eigenvalues etc..
  type(symmetry) :: syms          ! symmetry operations
  type(energy) :: energs          ! energies and other information
  type(crystal) :: crys           ! crystal structure
  type(molecular_dynamics) :: MD  ! molecular dynamics structure
  !
  !     I/O variables:
  !     ------------
  !

  integer, parameter :: stdout = 6

  integer, parameter :: logfile_c = 7
  integer, parameter :: convergefile = 8
  integer, parameter :: pwlogfile = 34
  character(len=5), parameter :: input_file = 'input'
  character(len=6), parameter :: pw_logfile = 'PW_LOG'
  character(len=9), parameter :: struct_logfile = 'STRUC_LOG'
  character(len=14), parameter :: converge_logfile = 'STRUC_CONVERGE'
  integer :: ios, istherecheckpoint
  !
  !     : BASIC VARAIABLES
  !     ---------------------
  ! 
  real(dp), allocatable :: &
      x(:),&   ! 1st 9 l.v. - rest atomic postions in direct units 
      g(:)     ! 1st 9 stress - rest forces in units of lvec*F_cart
  real(dp) :: f1, t0
  integer :: plen,&               ! packed vector length
             pwcount = 0,&        ! # of times total energy (pwmain) is calc.
             iwfnreuse_c          ! whether wavefunctions can be eused 
                                  ! for extrapolation option   
  character(len=llength) &
         jobstring,&              ! choose for types of relaxation or MD
         relaxstring              ! type of relaxation
  !
  !     RELAXATION VARIABLES:
  !     ---------------------
  ! 
  real(dp), allocatable :: &
      H(:), & ! approximation in inverse of dynamical matrix
      H0(:)   ! inital value of H - kept for checkpointing

  real(dp) relaxstress,&          ! flag if force are to be relaxed
           relaxforce,&           ! flag if force are to be relaxed
           relaxepsfac, &         ! init stress part of H  -  step size for 
           relaxcoordfac,&        ! init force part of H   -  1st move
           accuracy,&             ! value of sqrt(g*g) deemed fully relaxed
           avmass, &              ! average mass of atoms
           pressure, &            ! pressure at which relaxation is done 
           lambdalimit,&          ! max value of step size in relaxation 
           lambda_tolerance, &    ! tol from 1d0 that lambda  is refined 
           lambda2_tolerance, &   ! tol from refinement if 2nd refine is done
           eps(9), &              ! inital deformation tensor 
           fixed_vector(3), &     ! fixed vector in unit cell relaxation 
           vector_plane(3),&      ! with prev. vector forms plane to stay fixed
           vinit                  ! initial volume

  integer  relax_method,&         ! BFGS routine used 1- orginal 2  - lbfgs
           adjustpressure,&       ! flag to adjust pressure
           iupdateh,&             ! 0 if steepest descent not BFGS used 
           itmax,&                ! max number of relaxation steps
           decoup_rc              ! decouples coordinates -for molecular solids
  !
  !      VARIABLES:
  !     -------------
  !
  integer :: nproc = -1
  integer :: myproc = -1
  integer :: comm
  integer :: ierr = 0
  !
  !     MISC VARIABLES:
  !     -------------
  !
  integer ::  i, j, k,ip,packlen,natom_tot
  real(dp)   lvecf(3, 3), feps(3, 3), lvect_inv(3,3) 
  real(dp) ::gimmetime, matinvert, dotprod, funeval,temp
  logical :: switchlv= .false.
  logical irelax
  !
  !     -------------------------- BEGIN  -----------------------
  !


  call MPI_Init(ierr) 
  call MPI_Comm_rank(MPI_COMM_WORLD, myproc, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
  call MPI_Comm_dup(MPI_COMM_WORLD, comm, ierr)


  iwfnreuse_c = 1  ! reuse wavefunctions by default
  t0 = gimmetime()
  !
  ! --------- setup i/o ------------------
  !
  if (myproc == 0) then
     ! Move old logfiles out of the way
     call mvfile(pw_logfile, '.old')
     call mvfile(struct_logfile, '.old')
     call mvfile(converge_logfile, '.old')

     call esdf_init(input_file)

     open(34, file = 'PW_LOG',  &
      status = 'unknown', form = 'formatted')

     open(unit = logfile_c, file = struct_logfile, iostat = ios, &
          form = 'formatted', status = 'replace')
     if (ios /= 0) call mystop( 'Cannot open log file!' )

     open(unit = convergefile, file = converge_logfile, iostat = ios, &
          form = 'formatted', status = 'replace')
     if (ios /= 0) call mystop( 'Cannot open converge file!' )

  end if
  !
  ! --------- read crystal structure and initialize variable crys -------------
  !
  call crystal_struc(myproc,comm, convergefile,switchlv,crys)
  feps=crys%lvec
  vinit=matinvert(feps(1,1))
  !
  ! --------- get some other parameters needed ------------------
  !
  ! mass for the plane wave code
  crys%totalmass = dzero
  do i = 1, crys%ntype
   crys%totalmass  =  crys%totalmass + crys%mass(i) * crys%natom(i)
  end do

  ! relaxation parameters for the coordinates
  avmass = dzero  ! compute average mass
  do i = 1, crys%ntype
     avmass = avmass + crys%mass(i)
  end do
  avmass = avmass / real(crys%ntype, dp)
  !
  ! Get some other parameters needed
  !
  call read_param(myproc, comm, eps, fixed_vector, vector_plane, &
       decoup_rc, relax_method, MD,&
       relaxepsfac, relaxcoordfac, accuracy, relaxstring, pressure, &
       adjustpressure, jobstring, lambdalimit, relaxstress, relaxforce, &
       iupdateh,lambda_tolerance, lambda2_tolerance, itmax, avmass, crys%vcell)

  pw_params%extrapolation=MD%extrapolation
  call fix_rotation(eps,fixed_vector, vector_plane)

  if (myproc == 0) then
     if (decoup_rc == 1) then
        write(convergefile, '(A)') &
             'DECOUPLING REALSPACE COORDINATES! COULD BREAK THE SYMMETRY!'
     end if

     write(logfile_c, '(3A)') '--------RELAXATION STYLE: ', &
          trim(relaxstring), ' ---------'
     write(logfile_c, '(A,F10.6)') 'relax_eps_factor = ', relaxepsfac
     write(logfile_c, '(A,F10.6)') 'relax_coord_factor = ', relaxcoordfac
     write(logfile_c, '(A,I1)') 'Hessian update: ', iupdateh

     write(logfile_c, '(3A)') '--------JOB: ', trim(jobstring), ' ---------'

  end if
  !
  ! read all the parameters for the planewave code
  !
  call read_pwparam(myproc, comm,pw_params,altkpoints,bands,&
                      altbands,syms,energs,crys,crys%lvec)
  !
  ! correct k-points if input lattice vectors and thus k-points are not
  ! in a right-handed cooredinate system
  !
  if(switchlv) then
    i=altkpoints%grid(3)  
    altkpoints%grid(3)=altkpoints%grid(2)
    altkpoints%grid(2)=i
  
    temp= altkpoints%shift(3)
    altkpoints%shift(3)  = altkpoints%shift(2)
    altkpoints%shift(2)  =temp
  end if

  if (myproc == 0) then
     if (pw_params%epsdiag > accuracy * accuracy) then
        write(logfile_c, '(A)') &
             '*** WARNING: Matrix diagonalization not accurate enough ***'
        write(logfile_c, '(A,G12.6)') &
             '***          Decrease accuracy_diag to be less than ', &
             accuracy * accuracy
     end if

     if (pw_params%epscv > accuracy) then
        write(logfile_c, '(A)') &
             '*** WARNING: Potential is not computed accurately enough ***'
        write(logfile_c,'(A,G12.6)') '***          Decrease &
             &potential_convergence_criterion to be less than ', &
             accuracy * accuracy
     end if

  end if
  !
  ! --------- allocate the memory needed -------
  !
  plen = packlen(crys%ntype, crys%natom)
  allocate(x(plen), g(plen), H(plen * plen), H0(plen * plen))
  irelax=.true.
  !
  ! depending on jobstring do different inital tasks
  !
  if (trim(jobstring) == 'relax_from_checkpoint') then
     call pack(crys%ntype, crys%mxdatm, plen, crys%natom, x, crys%coord, eps) 
     call readcheckpoint(plen, myproc, comm, x, g, H, H0, f1, crys%lvec)
     pw_params%ilog = 0
  else if (trim(jobstring) == 'relax_recycle_H') then
     call pack(crys%ntype, crys%mxdatm, plen, crys%natom, x, crys%coord, eps)
     call initialize_hessian(plen, x, H, relaxepsfac, relaxcoordfac, -1, &
            decoup_rc, crys%lvec)
     H0 = H
     call readhessian(plen, myproc, comm, H, H0, crys%lvec)
     pw_params%ilog = 0
     ! calculate stress and forces
     f1 = funeval(t0, plen, x, pressure, g, myproc, relaxstress, relaxforce, &
           iwfnreuse_c, adjustpressure,pw_params,altkpoints,&
          bands,altbands,syms,energs,crys,vinit)
     pwcount = pwcount + 1
  else if (trim(jobstring) == 'band_structure') then
     write(*, '(A)') ' *** job band_structure is not allowed any more'
     write(*, '(A)') ' *** use the pw_job variable instead.'
     call mystop
  else if (trim(jobstring) == 'relax') then ! pack initial values in x
     call pack(crys%ntype, crys%mxdatm, plen, crys%natom, x, crys%coord, eps)
     call initialize_hessian(plen, x, H, relaxepsfac, relaxcoordfac, -1, &
          decoup_rc, crys%lvec)
     H0 = H    
     pw_params%ilog = 0
     ! calculate stress and forces
     f1 = funeval(t0, plen, x, pressure, g, myproc, relaxstress, relaxforce, &
           iwfnreuse_c, adjustpressure,pw_params,altkpoints,&
          bands,altbands,syms,energs,crys,vinit)
     pwcount = pwcount + 1
  else if (trim(jobstring) == 'MD') then ! pack initial values in x
     MD%n_deg_freedom=0
     natom_tot=0
     do i=1,crys%ntype
       do j=1,crys%natom(i)
        natom_tot=natom_tot+1
        do k=1,3
          MD%n_deg_freedom=MD%n_deg_freedom+1
        end do
       end do
     end do
     MD%n_deg_freedom=max(MD%n_deg_freedom-3,3)
     lvect_inv = transpose(crys%lvec)
     temp = matinvert(lvect_inv(1, 1))
     irelax=.false.
     allocate(MD%velocity(plen))
     if (myproc == 0) then     
       call initialize_velocity(MD,plen,crys%mxdatm,crys%ntype,crys%natom,&
                          lvect_inv,crys%mass,natom_tot,crys%lvec) 
     end if
     call my_broadcast(MD%velocity,plen,0)

     call pack(crys%ntype, crys%mxdatm, plen, crys%natom, x, crys%coord, eps)
     relaxstress=dzero
     ! calculate stress and forces
     f1 = funeval(t0, plen, x, pressure, g, myproc, relaxstress, relaxforce, &
           iwfnreuse_c, adjustpressure,pw_params,altkpoints,&
          bands,altbands,syms,energs,crys,vinit)
     pwcount = pwcount + 1
  else
     if (myproc == 0) write(logfile_c, '(2A)') &
          'No valid job variable set: ',trim(jobstring)
     call mystop
  end if
  !
  ! depenending on job string either relax crystal or do MD
  !
  if (irelax) then   

     call relax(f1,plen,H,H0,g,x,t0, pressure, myproc,pwcount,accuracy, &
          iupdateh, relaxepsfac, relaxcoordfac,decoup_rc,lambdalimit,&
          relaxstress, relaxforce, iwfnreuse_c, adjustpressure,&
          relax_method,stdout,&
          lambda_tolerance, lambda2_tolerance,itmax,fixed_vector,vector_plane,&
          pw_params,altkpoints,bands,altbands,syms,energs,crys,vinit)     

  else

     call mol_dyn(plen,MD,x,g,pwcount,pw_params,t0, pressure,&
                     myproc, relaxstress, relaxforce, &
                     iwfnreuse_c, adjustpressure,altkpoints,&
                     bands,altbands,syms,energs,crys)

  end if  ! for relax or MD
  !
  ! FINAL OUTPUT
  !
  call apply_eps_lat(x, crys%lvec, lvecf)

  if (myproc == 0) then
     write(logfile_c, '(A)') 'final lattice vectors:'
     do i = 1, 3
        write(logfile_c, '(3F9.6)') lvecf(1,i), lvecf(2,i), lvecf(3,i)
     end do

     write(logfile_c, '(A)') 'H matrix:'

     k = 1
     do i = 1, plen
        do j = 1, plen, 3
           write(logfile_c, '(3(F9.6,X))') H(k:k+2)
           k = k + 3
        end do
        write(logfile_c, *)
     end do

     write(logfile_c, '(A,I3,A)') 'It took ', pwcount, &
          ' call(s) to the pw program!'

  end if


  deallocate(x, g, H, H0)

  if (myproc == 0) then
     close(logfile_c)
     close(convergefile)
  end if


  call MPI_Finalize(ierr)


end program main
!
! ============================ support routines ===================
!
function istherecheckpoint(myproc)

  use constants
  implicit none
  include 'flibcalls.ph'
include 'mpif.h'

  integer :: istherecheckpoint
  integer, parameter :: logfile_c = 7

  integer, intent(in) :: myproc
  character(len=10), parameter :: checkpointfile = 'CHECKPOINT'
  integer :: ios
  integer, parameter :: ckpfile = 10

  if (myproc == 0) then
     open(unit = ckpfile, file = checkpointfile, iostat = ios, &
          form = 'formatted', status = 'old')
     if (ios /= 0) then
        write(logfile_c,'(A)') 'THERE IS NO CHECKPOINT FILE.'
        istherecheckpoint = 0
     else 
        istherecheckpoint = 1
     endif
     close(unit = ckpfile)
  endif
  return
end function istherecheckpoint

subroutine readcheckpoint(plen, myproc, comm, x, g, H, H0, f1, lvec)

  use constants
  implicit none
  include 'flibcalls.ph'
include 'mpif.h'

  integer, parameter :: logfile_c = 7

  integer, intent(in) :: plen, myproc, comm
  real(dp), intent(out) :: x(plen), g(plen), H(plen * plen), H0(plen * plen), &
       f1
  real(dp), intent(in) :: lvec(3, 3)
  real(dp) :: dumdum

  character(len=10), parameter :: checkpointfile = 'CHECKPOINT'
  integer :: i, ios, icount, ierr
  integer, parameter :: ckpfile = 10
  character(len=10) :: cdummy

  if (myproc == 0) then

     write(logfile_c, '(A)') ' READING FROM CHECKPOINT FILE:'
     open(unit = ckpfile, file = checkpointfile, iostat = ios, &
          form = 'formatted', status = 'old')
     if (ios /= 0) call mystop( 'CANNOT OPEN CHECKPOINT FILE -- STOP' )

     read(ckpfile, *)
     do i = 1, plen
        read(ckpfile, *) x(i)
     end do

     read(ckpfile, *)
     do i = 1, plen
        read(ckpfile, *) g(i)
     end do

     read(ckpfile, *)
     do i = 1, plen * plen
        read(ckpfile, *) H(i)
     end do

     read(ckpfile, *)
     do i = 1, plen * plen
        read(ckpfile, *) H0(i)
     end do

     read(ckpfile, '(A10,G30.20)') cdummy, f1

     close(unit = ckpfile)

     write(logfile_c, '(A)') 'checkpoint state:'
     call printpacked(plen, x, lvec)
     dumdum = mdnrm2(plen, g(1), 1)
     write(logfile_c, '(A,F9.6,A)') 'checkpoint gradient has length ', &
          dumdum, ':'
     call printpacked(plen, g, lvec)
     write(logfile_c, '(A,F14.6)') 'checkpoint enthalpy = ', f1

     call myflush(logfile_c)

  end if


  icount = plen
  call MPI_Bcast(x(1), icount, MPI_DOUBLE_PRECISION, 0, comm, ierr)
  call MPI_Bcast(g(1), icount, MPI_DOUBLE_PRECISION, 0, comm, ierr)
  icount = plen * plen
  call MPI_Bcast(H(1), icount, MPI_DOUBLE_PRECISION, 0, comm, ierr)
  call MPI_Bcast(H0(1), icount, MPI_DOUBLE_PRECISION, 0, comm, ierr)
  call MPI_Bcast(f1, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)


end subroutine readcheckpoint

 
 subroutine readhessian(plen, myproc, comm, H, H0, lvec)
 
   use constants
   implicit none
   include 'flibcalls.ph'
 include 'mpif.h'
 
   integer, parameter :: logfile_c = 7
 
   integer, intent(in) :: plen, myproc, comm
   real(dp), intent(out) :: H(plen * plen), H0(plen * plen)
   real(dp), intent(in) :: lvec(3, 3)
   real(dp) :: dumdum
 
   character(len=10), parameter :: checkpointfile = 'CHECKPOINT'
   integer :: i, ios, icount, ierr
   integer, parameter :: ckpfile = 10
 
   if (myproc == 0) then
 
      write(logfile_c, '(A)') ' READING HESSIAN FROM CHECKPOINT FILE:'
      open(unit = ckpfile, file = checkpointfile, iostat = ios, &
           form = 'formatted', status = 'old')
      if (ios /= 0) then
         write(logfile_c, '(A)') '  CANNOT OPEN CHECKPOINT FILE!'
         write(logfile_c, '(A)') &
              '  Changing job to just relax (no H to recycle).'
      else
         read(ckpfile, *)
         do i = 1, plen
            read(ckpfile, *) dumdum
         end do
         read(ckpfile, *)
         do i = 1, plen
            read(ckpfile, *) dumdum
         end do
         
         read(ckpfile, *)
         do i = 1, plen * plen
            read(ckpfile, *) H(i)
         end do
         
         read(ckpfile, *)
         do i = 1, plen * plen
            read(ckpfile, *) H0(i)
         end do
         
         close(unit = ckpfile)
         
         write(logfile_c, '(A)') 'Finished reading hessian.'
      endif
      call myflush(logfile_c)
 
   end if
 
 
   icount = plen
   icount = plen * plen
   call MPI_Bcast(H(1), icount, MPI_DOUBLE_PRECISION, 0, comm, ierr)
   call MPI_Bcast(H0(1), icount, MPI_DOUBLE_PRECISION, 0, comm, ierr)
 
 
 end subroutine readhessian


