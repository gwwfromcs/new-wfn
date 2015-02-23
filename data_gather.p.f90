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
subroutine data_gather(ipr, myproc, nproc, indata, outdata, ng)
  !
  !     collects the array data from all procs to processor 0
  !
  use constants
  implicit none             ! implicit? no!
  !


include 'mpif.h'
  !
  !     INPUT:
  !
  integer, intent(in) :: &
       myproc, nproc, &
       ipr, &                                                     ! print flag
       ng                                  ! number of vectors in LOCAL gspace
  complex(dp), intent(in) :: indata(*)                  ! the local data array 
  !
  !     OUTPUT:
  !     
  complex(dp), intent(out) :: outdata(*)               ! the global data array
  !
  !     --------- local variables ---------------------------
  !
  integer :: i, j, k, info
  integer, allocatable :: recvcount(:), displs(:)
  !
  if (nproc > 1) then

     !
     !        set up displs and recvcount arrays
     !
     allocate(recvcount(nproc))
     allocate(displs(nproc))
     recvcount = 0
     displs= 0

     displs(myproc + 1) = ng
     call mpi_allreduce(displs, recvcount, nproc, MPI_INTEGER, &
          MPI_SUM, MPI_COMM_WORLD, info)

     j = 0
     do i = 1, nproc
        displs(i) = j
        j = j + recvcount(i)
     end do
     !
     !        now gather everything to processor 0
     !            
     call mpi_gatherv(indata, ng, MPI_DOUBLE_COMPLEX, &
          outdata, recvcount, displs, MPI_DOUBLE_COMPLEX, &
          0, MPI_COMM_WORLD, info)
     deallocate(recvcount)
     deallocate(displs)

  else
     do i = 1, ng
        outdata(i) = indata(i)
     end do
  end if
  !
  !     now processor 0 has the global data array
  !
  return

end subroutine data_gather


