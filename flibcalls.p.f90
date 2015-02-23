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

!-*-Fortran-*-
!
module flibcalls_module

  contains
  !
  !     -------------- diagonalization routine ---------------------
  !
  subroutine zdiagh(n, a, lda, w, work, lwork, rwork) 

    use constants
    implicit none

    integer :: lda, lwork, n, info
    real(dp) :: rwork(*)
    complex(dp) :: work






    complex(dp), intent(inout) :: a(*)
    real(dp), intent(out) :: w(*)
!    complex(dp), intent(inout) :: a(lda * n)
!    real(dp), intent(out) :: w(n)





    !
    !     --------- local variables ----------------------------
    !
    integer :: i, j, k
    complex(dp), allocatable :: dmat(:)
    character :: jobz, uplo

    uplo = 'L'
    jobz = 'V'



    call zheev(jobz, uplo, n, a, lda, w, work, lwork, rwork(1), info)



    if(info /= 0) write(9,110) info
    return

110 format('Error in submatrix diagonalization:',i10)

  end subroutine zdiagh

end module flibcalls_module
