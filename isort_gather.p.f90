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
subroutine isort_gather(ipr, mysort, gsstruct, gs)
  !
  !     creates global isort array
  !


  use vqmc_gspace_module
  include 'use.h'
  implicit none             ! implicit? no!
include 'mpif.h'
  include 'interface.h'
  !
  !     1996 Bernd Pfrommer
  !
  !
  !     INPUT:
  !     -----
  !
  type(gspace), intent(in) :: &
       gsstruct                      ! a structure containing the whole gspace
  type(parallel_gspace), intent(in) :: gs                ! the parallel gspace 
  integer, intent(in) :: &
       ipr                  ! print flag
  !
  !     OUTPUT:
  !     ------
  !
  integer, intent(out) :: &
       mysort(1:gs%length)          ! the global map between data and gsstruct
  !
  !
  !     DESCRIPTION:
  !     -----------
  !
  !     Sets up the isort array from the parallel gspace, which
  !     associates:
  !
  !      data(mysort(i)) <=> gstruct%gvec(:,i)
  !
  !
  !     --------- local variables ---------------------------
  !
  integer :: i, j, k, info
  integer, external :: findvec
  integer, allocatable :: recvcount(:), displs(:), ibuf(:)
  !
  !     allocate buffer space
  !
  allocate(ibuf(gs%length))

  k = 0

  do i = 1, gsstruct%length
     j = findvec(gsstruct%gvec(1, i), gs)
     if (j > 0) then
        k = k + 1
        ibuf(j) = i
     end if
  end do

  if (k /= gs%length) then
     write(0, *) 'gspace does not enclose data completely!'
     call mystop
  end if

  if (gs%nproc > 1) then

     !
     write(9, *) '** WARNGING isort_gather: not debugged!'

     !
     !        set up displs and recvcount arrays
     !
     allocate(recvcount(gs%nproc))
     allocate(displs(gs%nproc))
     recvcount = 0
     displs = 0

     displs(gs%myproc) = k

     call mpi_reduce(displs, recvcount, MPI_INTEGER, MPI_SUM, 0, &
          MPI_COMM_WORLD, info)

     j = 0
     do i = 1, k
        displs(i) = j
        j = j + recvcount(i)
     end do
     !
     !        now gather everything to processor 0
     !            
     call mpi_gatherv(ibuf, k, MPI_INTEGER, &
          mysort, recvcount, displs, MPI_INTEGER, &
          0, MPI_COMM_WORLD, info)
     deallocate(recvcount)
     deallocate(displs)

  else
     mysort(1:gs%length) = ibuf(1:gs%length)
  end if
  !
  !     now processor 0 has the global isort array
  !
  deallocate(ibuf)

  return

end subroutine isort_gather
