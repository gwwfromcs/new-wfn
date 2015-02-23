!     @process extchk
!
subroutine reg_grid_spline(y, n, y2)  
  !
  !     computes the second derivative array y2 for a function
  !     on a regular grid spaced 1
  !
  !     second derivatives at boundaries set to zero
  !
  !
  !     (1996) F. Mauri B. Pfrommer
  !
  use constants
  implicit none    ! never
  integer, intent(in) :: n     ! number of array points
  real(dp), intent(in) :: y(n) ! function values at 1,2,3....n
  !
  !     OUTPUT:
  !     ------
  !
  real(dp), intent(out) :: y2(n)    ! second derivatives
  !
  !     ------------- local variables --------------------------
  !
  real(dp) :: p, qn, sig, un  
  integer :: i, k  
  real(dp), allocatable :: u(:)  
  allocate(u(n))  
  !     ..
  y2(1) = dzero  
  u(1) = dzero
  sig = dhalf
  do i = 2, n - 1  
     p = sig * y2(i - 1) + dtwo
     y2(i) = (sig - done) / p  
     u(i) = (dthree * ((y(i + 1) - y(i)) - (y(i) - y(i - 1))) - &
          sig * u(i - 1)) / p
  end do
  qn = dzero
  un = dzero 

  y2(n) = (un - qn * u(n - 1)) / (qn * y2(n - 1) + done)  
  do k = n - 1, 1, -1  
     y2(k) = y2(k) * y2(k + 1) + u(k)  
  end do

  deallocate(u)  

  return

end subroutine reg_grid_spline
