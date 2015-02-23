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
!       Provide correct interfaces for the various libraries
!
!	Bernd, 1996	
!
!
!
!     ---------------------- LAPACK -------------------------------
!

interface mdnrm2

  real(kind(1.d0)) function dnrm2(n,x,incx)

    integer :: n,incx
    real(kind(1.d0)) :: x
  end function 

end interface

interface mdsyev

  subroutine dsyev(v, l, m, a1, n, a2, a3, k, i)

    character :: v, l
    real(kind(1.0d0)) :: a1
    real(kind(1.0d0)):: a2, a3
    integer :: m, n, k, i
  end subroutine
end interface

 
interface mdspev

  subroutine dspev(v, l, m, a1, n, a2, a3, k, i)

    character :: v, l
    complex(kind(1.0d0)) :: a1
    real(kind(1.0d0)):: a2, a3
    integer :: m, n, k, i
  end subroutine 



end interface 

interface mzheev

  subroutine zheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)

      character :: jobz, uplo
      integer :: info, lda, lwork, n
      real(kind(1.0d0)) :: rwork, w
      complex(kind(1.0d0)) :: a, work


  end subroutine zheev



end interface

interface mzpotrf

  subroutine ZPOTRF(uplo,n ,a,lda,info )

      character ::  uplo
      integer :: info, lda,  n
      complex(kind(1.0d0)) :: a


  end subroutine zpotrf


end interface



interface mzpotri

  subroutine ZPOTRI(uplo,n ,a,lda,info )

      character ::  uplo
      integer :: info, lda,  n
      complex(kind(1.0d0)) :: a


  end subroutine zpotri


end interface

interface mztrtri

  subroutine ZTRTRI(uplo,diag,n ,a,lda,info )

      character ::  uplo,diag
      integer :: info, lda,  n
      complex(kind(1.0d0)) :: a


  end subroutine ztrtri



end interface

!
!     
!     --------------------------- BLAS 1 ---------------------------
!

interface mzaxpy

  subroutine zaxpy(n, alpha, x, incx, y, incy)

    integer :: n, incx, incy
    complex(kind(1.0d0)) :: x, y, alpha

  end subroutine
end interface

interface mzcopy

  subroutine zcopy(n, x, incx, y, incy)

    integer :: n, incx, incy
    complex(kind(1.0d0)) x, y

  end subroutine
end interface 

interface mzdotc

  complex(kind(1.0d0)) function zdotc(n, x, incx, y, incy)

    integer :: n, incx, incy
    complex(kind(1.0d0)) x, y

  end function
end interface 

interface mzdotu

  complex(kind(1.0d0)) function zdotu(n, x, incx, y, incy)

    integer :: n, incx, incy
    complex(kind(1.0d0)) x, y

  end function
end interface 

interface mzscal

  subroutine zscal(n, alpha, x, incx)

    integer :: n, incx
    complex(kind(1.0d0)) :: x
    complex(kind(1.0d0)) :: alpha

  end subroutine

end interface 

interface mzdscal

  subroutine zdscal(n, alpha, x, incx)

    integer :: n, incx
    complex(kind(1.0d0)) :: x
    real(kind(1.0d0)) :: alpha

  end subroutine

end interface 


interface mdaxpy

  subroutine daxpy(n, alpha, x, incx, y, incy)

    integer :: n, incx, incy
    real(kind(1.0d0)) :: x, y, alpha

  end subroutine
end interface

interface mdcopy

  subroutine dcopy(n, x, incx, y, incy)

    integer :: n, incx, incy
    real(kind(1.0d0)) x, y

  end subroutine
end interface 

interface mddot

  real(kind(1.0d0)) function ddot(n, x, incx, y, incy)

    integer :: n, incx, incy
    real(kind(1.0d0)) x, y

  end function
end interface 
 
interface mdscal

  subroutine dscal(n, alpha, x, incx)

    integer :: n, incx
    real(kind(1.0d0)) :: x
    real(kind(1.0d0)) :: alpha

  end subroutine

end interface 


!
!     
!     --------------------------- BLAS2 ---------------------------
!

interface mzgemv

  subroutine zgemv(ta, m, n, alf, a, lda, x, incx, bet, y, incy)

    integer :: m, n, lda, incx, incy
    complex(kind(1.0d0)) alf, a, x, bet, y
    character :: ta

  end subroutine

end interface 

!
!     
!     --------------------------- BLAS3 ---------------------------
!

interface mzsymm

  subroutine zsymm(ta, tb, m, n, alf, a, lda, b, ldb, bet, c, ldc)

    integer :: m, n, k, lda, ldb, ldc
    complex(kind(1.0d0)) :: alf, a, b, bet, c
    character :: ta, tb

  end subroutine

end interface 

interface mztrmm

  subroutine ztrmm(side, uplo, transa ,diag , m, n, alf, a, lda, b, ldb)

    integer :: m, n, lda, ldb
    complex(kind(1.0d0)) :: alf, a, b, bet
    character :: side, uplo, transa ,diag

  end subroutine

end interface 

interface mdgemm

  subroutine dgemm(ta, tb, m, n, k, alf, a, lda, b, ldb, bet, c, ldc)

    integer :: m, n, k, lda, ldb, ldc
    real(kind(1.0d0)) :: alf, a, b, bet, c
    character :: ta, tb

  end subroutine

end interface 

interface mzgemm

  subroutine zgemm(ta, tb, m, n, k, alf, a, lda, b, ldb, bet, c, ldc)

    integer :: m, n, k, lda, ldb, ldc
    complex(kind(1.0d0)) :: alf, a, b, bet, c
    character :: ta, tb

  end subroutine

end interface 

interface mzher2k

  subroutine zher2k(uplo, trans, n, k, alf, a, lda, b, ldb, bet, c, ldc)

    integer :: n, k,lda, ldb, ldc
    complex(kind(1.0d0)) :: alf, a, b, bet, c
    character :: uplo, trans

  end subroutine

end interface 

interface mzherk

  subroutine zherk(uplo, trans, n, k, alf, a, lda, bet, c, ldc)

    integer :: n, k, lda, ldc
    complex(kind(1.0d0)) :: alf, a, bet, c
    character :: uplo, trans

  end subroutine

end interface 
