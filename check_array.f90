!
subroutine check_array(data, datasiz, rmax, rmin, cmax, cmin, &
     eps, ierr)
  !
  !     examines real and imaginary part of an array.
  !     SETS IMAGINARY PART TO ZERO!
  !
  !
  !     1995 Bernd Pfrommer
  !
  use all_to_all_module  
  implicit none                 ! never do implicit!
  include 'all_to_all.h'  
  !
  !     INPUT:
  !     -----
  !
  integer, intent(in) :: datasiz              ! number of data points in array
  complex(dp), intent(inout) :: data(datasiz)                ! the data itself
  real(dp) :: eps                   ! maximum tolerated absolute value of imag
  !
  !     OUTPUT
  !
  real(dp), intent(out) :: &
       rmax, rmin, &                  ! max and min value of real part of data
       cmax, cmin                     ! max and min value of imag part of data
  integer, intent(out) :: ierr           ! number of data points with imag>eps
  !
  !     -------- local
  !
  integer :: i  
  !
  rmax = real(data(1), dp)  
  rmin = real(data(1), dp)  
  cmax = aimag(data(1))  
  cmin = aimag(data(1))  

  ierr = 0  
  do i = 1, datasiz
     if (real(data(i), dp) > rmax) rmax = real(data(i), dp)  
     if (real(data(i), dp) < rmin) rmin = real(data(i), dp)  
     if (abs(aimag(data(i))) > eps) ierr = ierr + 1  
     if (aimag(data(i)) > cmax) cmax = aimag(data(i))  
     if (aimag(data(i)) < cmin) cmin = aimag(data(i))  

     data(i) = cmplx(real(data(i), dp), dzero, dp)   ! clear complex component
  end do

  !     find max across all processors
  call all_max_all(rmax)  
  call all_max_all(cmax)  
  call all_min_all(rmin)  
  call all_min_all(cmin)  
  call all_sum_all(ierr)  

  return  

end subroutine check_array
