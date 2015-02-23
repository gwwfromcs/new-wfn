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
subroutine adjustfft(nfft, symms)

  include 'use.h'
  implicit none             ! implicit? Just say no
  include 'interface.h'

  !
  !     1996 Bernd Pfrommer
  !     
  !     

  !
  !     INPUT:
  !     -----

  type(symmetry) :: symms       ! the symmetry operations of the crystal
  !
  !     INPUT/OUTPUT:
  !     ------------

  integer :: nfft(3)            ! fft grid dimensions
  !
  !     DESCRIPTION
  !     -----------
  !
  !     Increases the fft grid in nfft(3) to fit a suitable value, e.g.
  !     a power of two, three ...
  !
  !     Also checks to make sure the resulting grid is compatible with
  !     the nonprimitive translations, e.g. if there is a translation
  !     of 1/3, the grid must be a multiple of 3.
  !
  !
  !     ------------------- local arrays ---------------------------

  integer :: i, j, k, ns
  real(dp) :: tau

  logical :: ifit




















  do i=1,3
    ! redundant check
    ifit = .false.
    do while( .not.(ifit) )
      nfft(i) = fastnum( nfft(i) )
      ! check if nfft(i) is compatible with symmetry
      ifit=.true. ! this check redundant for some reason
      ! if not increase nfft(i)
      if( .not.(ifit) ) nfft(i)=nfft(i)+1
    enddo
  enddo
 

  return

110 format(' *** ERROR IN ADJUSTFFT: ', &
       'FFT GRIDSIZE EXCEEDS IN DIRECTION',i3, &
       /' FFT GRIDSIZE = ', 3i6, &
       /' MAX GRIDSIZE = ', i6, &
       /'     THIS CAN BE DUE TO STRANGE NONPRIMITIVE', &
       /'     TRANSLATIONS IN THE SYMMETRY OPERATIONS', &
       /'     OR JUST BECAUSE THE SYSTEM IS TOO BIG.', &
       /'     ADD LARGER FFT SIZES IN ADJUSTFFT OR', &
       /'     IMPROVE THE COORDINATES/LATTICE VECTORS.')


  contains

  function fastnum( nfft_in )

  integer :: fastnum
  integer,intent(in) :: nfft_in

  integer :: nfft_out
  integer :: n, i
  integer :: prim(6), primpower(6)
  logical :: notfound

!  - best grids are given by : 2^a*3^b*5^c*7^d*11^e*13^f
!        a,b,c,d are arbitary, e,f are 0,1

  ! allowed prime factors
  prim = (/ 2, 3, 5, 7, 11, 13 /)

  nfft_out=nfft_in
  notfound=.true.
  do while( notfound )

    notfound=.false.
    n=nfft_out
    primpower=0

    ! generate prime power expansion
    do i=1,size(prim)

      do while( mod(n,prim(i))==0 )

        n=n/prim(i)
        primpower(i)=primpower(i)+1

      enddo

    enddo

    n=1
    do i=1,size(prim)
      n=n*prim(i)**primpower(i)
    enddo

    if( n/=nfft_out .or. primpower(5) > 1 .or. primpower(6) > 1 ) then
      notfound=.true.
      nfft_out=nfft_out+1
    endif

  enddo

  ! final result
  fastnum = nfft_out

  end function fastnum


end subroutine adjustfft


