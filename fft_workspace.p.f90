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

!
subroutine fft_workspace(nproc, fftsize, naux, ufbs, t1bufsize, mpistatsize)
  !
  !     determines size of fft workspaces for different machines
  !
  !     1996 Bernd Pfrommer
  !
  use constants
  implicit none             ! nope 
include 'mpif.h'
  !
  !     INPUT:
  !     -----
  !
  integer, intent(in) :: &
       nproc, &             ! number of processors
       fftsize(3)           ! size of the 3d fft grid
  !
  !     OUTPUT:
  !     ------
  !
  integer, intent(out) :: &
       ufbs, &              ! unfold buffer size
       t1bufsize, &         ! t1 buffer size
       naux(4), &           ! sizes of the various work arrays needed
       mpistatsize          ! size of mpi status variable
  !
  !     ----------------------- local variables  -----------------------------
  !     
  integer :: i, nmax, ix, iy

  t1bufsize = ((fftsize(1) * fftsize(3)) / nproc + 1) * fftsize(2) ! t1 buffer 
  ufbs = 1

mpistatsize = MPI_STATUS_SIZE

  naux(1:3) = fftsize(1:3) + 25
  naux(4)   = 1
 
 
    

      
      
      
      


end subroutine fft_workspace

!
!     =====================================================================
!
!     initialize the workspace for those machines where it is
!     necessary and possible
!
!
subroutine fft_iworkspace(workspace, naux, dim)

  !     1997 Bernd Pfrommer
  !
  use constants
  implicit none            
  !
  !     INPUT:
  !     -----
  !
  integer, intent(in) :: &
       naux(4), &           ! sizes of the fft workspace arrays
       dim(3)               ! fft sizes in x, y, z
  !
  !     OUTPUT:
  !     ------
  !
  complex(dp), intent(inout) :: workspace(*)  ! the workspace to be initialized
  !
  !     ---------------------- local variables ---------------------
  !
  
      !
      ! edit davegp
      ! I changed this parameter setting so that you can
      ! link to the correct header associated with the library
      ! that you link. Requires an include directory to be
      ! added to conf.mk
      !
      include 'fftw_f77.i'
      !
      ! INTEGER FFTW_FORWARD,FFTW_BACKWARD
      ! PARAMETER (FFTW_FORWARD=-1,FFTW_BACKWARD=1)
      ! INTEGER FFTW_ESTIMATE,FFTW_MEASURE
      ! PARAMETER (FFTW_ESTIMATE=0,FFTW_MEASURE=1)
      ! INTEGER FFTW_IN_PLACE,FFTW_USE_WISDOM
      ! PARAMETER (FFTW_IN_PLACE=8,FFTW_USE_WISDOM=16)

! for 64-bit machines (integer is not enough)
! it should also work on 32-bit machines

      ! the plan tags from create_plan
      INTEGER*8 FFTW_PLAN_FWD, FFTW_PLAN_BWD


  integer :: dwx, dwy, dwz
  real(dp) :: dum

  dwx = 1
  dwy = naux(1) + 1
  dwz = naux(1) + naux(2) + 1

  
  
  

end subroutine fft_iworkspace
