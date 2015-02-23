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
subroutine zdiag_occ_scalapack(context, n, nb, nprow, npcol, a, nc, nr, w,comm)
  !
  use all_to_all_module
  implicit none
  include 'all_to_all.h'
  include 'mpif.h'
  !
  !     INPUT:
  !     -----
  !
  integer, intent(in) :: &
       context, &           ! the blacs context 
       n, &                 ! matrix size
       nb, &                ! block size
       nc,nr, &             ! number of rows and columns local to this proc
       nprow, &             ! number of processor grid rows
       npcol, &                ! number of processor grid columns
       comm
  !
  !     INPUT/OUTPUT:
  !     ------------
  !
  complex(dp), intent(inout) :: &
       a(nr,nc)             ! the matrix to be diagonalized (destroyed)
  !
  !     OUTPUT:
  !     ------
  !
  real(dp), intent(inout) :: &
       w(n)                 ! all eigenvalues found
  !
  !     -------------- local variables ----------------------------
  !

  include 'scalapack.h'

  integer :: &
       info, i, j, &
       nn, nnp, np0, mq0, lwork, lrwork, liwork, &
       neig, &              ! number of eigenvalues sought for
       nfound, &            ! number of eigenvalues found
       nzfound, &           ! number of eigenvectors found
       lda, &               ! leading dimension
       iam, &               ! blacs myproc
       nproc                ! hopefully the same as mpi nproc
  integer :: &
       desca(lld_), &       ! descriptor for the matrix
       descz(lld_)          ! descriptor for the vectors
  integer, parameter :: izero = 0, &
       clustersize = 40     ! max size of eigenvalue clusters
  real(dp), parameter :: abstol = 1.0d-8, orfac = abstol * 1.0d2

  complex(dp), allocatable :: work(:), z(:,:)
  real(dp), allocatable :: rwork(:)
  integer, allocatable :: iwork(:)
  !
  !     those can always be non-aligned -> use F90
  !
  real(dp), allocatable :: gap(:)
  integer, allocatable :: ifail(:), iclustr(:)
  complex(dp) :: ssum
  integer mpibuf,ierr

  lda = nr                  ! close packing
  neig = n                  ! want all eigenvalues

!  call all_max_all(lda,comm)
  call descinit(desca, n, n, nb, nb, 0, 0, context, lda, info)
  call descinit(descz, n, n, nb, nb, 0, 0, context, lda, info)
 
  nn = max(n, nb, 2)
  nnp = max(n, nprow * npcol + 1, 4)
  np0 = numroc(nn, nb, 0, 0, nprow)
  mq0 = numroc(max(neig, nb, 2), nb, 0, 0, npcol)

  lwork = n + (np0 + mq0 +nb) * nb ! new documentation, but different case
  lrwork = 4 * n + max(5 * nn, np0 * mq0) + iceil(neig, nprow * npcol) * nn + &
       (clustersize - 1) * n
  liwork = 6 * nnp
 
  allocate(z(nr, nc))
  allocate(work(lwork))
  allocate(rwork(lrwork))
  allocate(iwork(liwork))

  allocate(ifail(n))
  allocate(iclustr(2 * nprow * npcol))
  allocate(gap(nprow * npcol))




  call pzheevx('V', 'A', 'U', n, a(1, 1), 1, 1, desca(1), dzero, dzero, &
       izero, izero, abstol, nfound, nzfound, w(1), orfac, &
       z(1, 1), 1, 1, descz(1), work(1), lwork, rwork(1), lrwork, &
       iwork(1), liwork, ifail(1), iclustr(1), gap(1), info)
  a = z


  deallocate(z)
  deallocate(iwork)
  deallocate(rwork)
  deallocate(work)

  deallocate(gap)
  deallocate(iclustr)
  deallocate(ifail)

  if (info == 2) then
     write(9, 100)
  else if (info /= 0) then
     write(9, *) ' *** zdiag_scalapack warning:', &
          ' *** diag failed with info=', info
  end if

  if (nzfound /= nfound) then
     write(9, *) ' *** zdiag_scalapack warning:', &
          ' *** found ', nzfound,' eigenvalues only !'
  end if

100 format(' ***  WARNING: DEGENERATE EIGENVECTORS NOT',&
       ' ORTHOGONALIZED')


end subroutine zdiag_occ_scalapack
