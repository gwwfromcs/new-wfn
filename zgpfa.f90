!        SUBROUTINE 'ZGPFA'
!        SELF-SORTING IN-PLACE GENERALIZED PRIME FACTOR (COMPLEX) FFT
!
!        *** THIS IS THE ALL-FORTRAN VERSION ***
!            -------------------------------
!
!        CALL ZGPFA(Z,TRIGS,LDZ,N,NTRANS,ISIGN,IER)
!
!        Z IS FIRST COMPLEX INPUT/OUTPUT VECTOR
!        TRIGS IS A TABLE OF TWIDDLE FACTORS, PRECALCULATED
!              BY CALLING SUBROUTINE 'SETGPFA'
!        INC IS THE INCREMENT WITHIN EACH DATA VECTOR
!        JUMP IS THE INCREMENT BETWEEN DATA VECTORS
!        N IS THE LENGTH OF THE TRANSFORMS:
!        ier is 0 on successful completion
!          -----------------------------------
!            N = (2**IP) * (3**IQ) * (5**IR)
!          -----------------------------------
!        LOT IS THE NUMBER OF TRANSFORMS
!        ISIGN = +1 FOR FORWARD TRANSFORM
!              = -1 FOR INVERSE TRANSFORM
!
!        WRITTEN BY CLIVE TEMPERTON
!        RECHERCHE EN PREVISION NUMERIQUE
!        ATMOSPHERIC ENVIRONMENT SERVICE, CANADA
!
!----------------------------------------------------------------------
!
!        DEFINITION OF TRANSFORM
!        -----------------------
!
!        X(J) = SUM(K=0,...,N-1)(C(K)*EXP(ISIGN*2*I*J*K*PI/N))
!
!---------------------------------------------------------------------
!
!        FOR A MATHEMATICAL DEVELOPMENT OF THE ALGORITHM USED,
!        SEE:
!
!        C TEMPERTON : "A GENERALIZED PRIME FACTOR FFT ALGORITHM
!          FOR ANY N = (2**P)(3**Q)(5**R)",
!          SIAM J. SCI. STAT. COMP., MAY 1992.
!
!----------------------------------------------------------------------

subroutine zgpfa(z, trigs, ldz, n, ntrans, isign, ier)

  use constants
  implicit none

  complex(dp), intent(inout) :: z(0:*)
  complex(dp), intent(in) :: trigs(0:*)
  integer, intent(in) :: ldz
  integer, intent(in) :: n
  integer, intent(in) :: ntrans
  integer, intent(in) :: isign
  integer, intent(out) :: ier
  integer :: nj(3)
  integer :: i
  integer :: nn, ifac
  integer :: ip, iq, ir
  integer :: kk, ll
  !
  !     decompose n into factors 2,3,5
  !     ------------------------------
  !
  nn = n
  ifac = 2

  do ll = 1, 3
     kk = 0
     do while (mod(nn, ifac) == 0)
        kk = kk + 1
        nn = nn / ifac
     end do
     nj(ll) = kk
     ifac = ifac + ll
  end do

  if (nn /= 1) then
     write(*, '(i10,a)') n, ' is not a legal value of transform length'
     ier = -1
     return
  end if

  ip = nj(1)
  iq = nj(2)
  ir = nj(3)
  !
  !     compute the transform
  !     ---------------------
  !
  i = 0

  if (ip > 0) then
     call zgpfa2f(z, trigs(i), ldz, n, ip, ntrans, isign)
     i = i + 2**ip
  end if

  if (iq > 0) then
     call zgpfa3f(z, trigs(i), ldz, n, iq, ntrans, isign)
     i = i + 3**iq
  end if

  if (ir > 0) then
     call zgpfa5f(z, trigs(i), ldz, n, ir, ntrans, isign)
  end if

  ier = 0

  return

end subroutine zgpfa
