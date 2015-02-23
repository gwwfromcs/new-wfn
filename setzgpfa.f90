!
subroutine setzgpfa(trigs, n, ndim, ier)
  !
  use constants
  implicit none
  !
  integer, intent(in) :: ndim, n
  integer, intent(out) :: ier
  complex(dp), intent(out) :: trigs(ndim)
  integer :: nj(3), nk(3), i, k, kk, ll, nn, ifac, kink, ni, irot
  real(dp) :: angle, del
  !
  !     decompose n into factors 2,3,5
  !     ------------------------------
  nn = n
  ifac = 2

  do ll = 1, 3
     kk = 0
     do while(mod(nn, ifac) == 0)
        kk = kk + 1
        nn = nn / ifac
     end do
     nj(ll) = kk
     ifac = ifac + ll
  end do

  if (nn /= 1) then
     write(*, '(i10,a)') n ,' is not a legal value of transform length'
     ier = -1
     return 
  end if
  !
  !     compute list of rotated twiddle factors
  !     ---------------------------------------
  nk(1) = 2**nj(1)
  nk(2) = 3**nj(2)
  nk(3) = 5**nj(3)

  i = 1
  do ll = 1, 3
     ni = nk(ll)
     if (ni > 1) then
        del = pi2 / real(ni, dp)
        irot = n / ni
        kink = mod(irot, ni)
        kk = 0
        do k = 1, ni
           angle = real(kk, dp) * del
           if (i < ndim) then
              trigs(i) = exp(cmplx(dzero, angle, dp))
           else
              ier = -2
           end if
           i = i + 1
           kk = kk + kink
           if (kk > ni) kk = kk - ni
        end do
     end if
  end do
  ier = 0

  return

end subroutine setzgpfa
