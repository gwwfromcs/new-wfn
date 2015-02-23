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
subroutine gspace_gather(myproc, nproc, gsstruc, &
     nfft, gsvec, ng, istar, nstar, mstar, phase, inds)
  !
  !     packs a parallel gspace into gsstruc
  !
  use constants
  use vqmc_gspace_module
  implicit none             ! implicit? no!
include 'mpif.h'


  !
  !     INPUT:
  !
  integer, intent(in) :: &
       myproc, nproc, &
       ng, &
       gsvec(3,ng), &       ! the local gvectors
       nfft(3), &           ! the fft grid size
       nstar, &             ! number of stars
       inds(ng), &          ! maps gvectors to stars
       mstar(nstar)         ! size of each star
  complex(kind=8), intent(in) :: &
       phase(ng)            ! phase factors
  logical, intent(in) :: &
       istar                ! if star information is present
  !
  !     OUTPUT:
  !     
  type(gspace), intent(out) :: &
       gsstruc              ! a structure containing the whole gspace
  !
  !     --------- local variables ---------------------------
  !
  integer :: i, j, ngtot, info
  !
  !     set the fft grid dimensions
  !
  gsstruc%fftsize(1) = nfft(1)
  gsstruc%fftsize(2) = nfft(2)
  gsstruc%fftsize(3) = nfft(3)
  !
  !     sum up across all processors to get total number of gvecs
  !
  ngtot = ng

  if (nproc > 1) &
       call mpi_allreduce(ng, ngtot, 1, MPI_INTEGER, &
       MPI_SUM, MPI_COMM_WORLD, info)


  gsstruc%length = ngtot
  !
  !     allocate space, and collect all the data
  !
  allocate(gsstruc%gvec(3, gsstruc%length))
  if (nproc > 1) then

     call mpi_allgather(gsvec, 3 * ng, MPI_INTEGER, &
          gsstruc%gvec, 3 * ng, MPI_INTEGER, &
          MPI_COMM_WORLD, info)

  else
     do i = 1, ng
        gsstruc%gvec(1, i) = gsvec(1, i)
        gsstruc%gvec(2, i) = gsvec(2, i)
        gsstruc%gvec(3, i) = gsvec(3, i)
     end do
  end if
  !
  !     deal with the star information
  !
  gsstruc%istar = istar

  if (istar) then
     !
     !        the mstar array is distributed already. No need to collect.
     !         
     gsstruc%nstar = nstar
     allocate(gsstruc%mstar(nstar))
     gsstruc%mstar = mstar
     !
     !        the other ones have to be collected.
     !
     !
     !        deal with the phase factors
     !
     allocate(gsstruc%phase(gsstruc%length))

     if (nproc > 1) then

        call mpi_allgather(phase(1), ng, MPI_DOUBLE_COMPLEX, &
             gsstruc%phase(1), ng, MPI_DOUBLE_COMPLEX, &
             MPI_COMM_WORLD, info)

     else
        gsstruc%phase = phase
     end if
     !
     !        handle the star index array
     !
     allocate(gsstruc%inds(gsstruc%length))
     if (nproc > 1) then

        call mpi_allgather(inds(1), ng, MPI_INTEGER, &
             gsstruc%inds(1), ng, MPI_INTEGER, &
             MPI_COMM_WORLD, info)

     else
        gsstruc%inds = inds
     end if
  else
     gsstruc%nstar = 0
     nullify(gsstruc%inds)
     nullify(gsstruc%mstar)
     nullify(gsstruc%phase)
  end if

  return

end subroutine gspace_gather







