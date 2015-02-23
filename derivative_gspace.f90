!
subroutine derivative_gspace (dir, data, gs, crys, grd)  
  !
  !     computes derivative wrt spatial direction: dir
  !     we compute data(g)*(k+G), so set the k-point of the gspace
  !     to zero if you just want the gradient!
  !
  !     Strictly speaking, it does NOT compute the gradient, but
  !     -i * gradient.
  !
  !
  !     (1996) Bernd Pfrommer
  !
  include 'use.h'  
  implicit none        ! never
  include 'interface.h'  
  !
  !     INPUT:
  !     ------
  !
  integer, intent(in) :: dir                   ! the spatial direction, 1,2,3
  type(crystal), intent(in) :: crys                      ! for the bvec array
  type(parallel_gspace), intent(in) :: gs                        ! the gspace
  complex(dp), intent(in) :: data(gs%length)      ! the data to be gradiented
  !
  !     OUTPUT:
  !     ------
  !
  complex(dp), intent(out) :: grd(gs%length)    ! the derivative of this data
  !
  !     ----------   local variables  --------------------------
  !
  !
  !     ----- variables for the gspace loop ---------
  !
  integer :: igv(4), fftn(4), ffth(4), igv3, irod, iord, igs  
  real(dp) :: gv(3)                                             ! if required
  !     ----------------------------------------------
  if ((dir < 1) .or. (dir > 3)) then  
     write(9, *) ' derivative_gspace: illegal direction:', dir  
     call mystop  
  end if
  fftn(1:3) = gs%fftsize(1:3)  
  fftn(4) = gs%fftsize(3)  
  ffth(:) = fftn(:) / 2  
  igs = 0  
  do iord = 1, gs%lorder                            ! loop through x/y gspace
     irod = gs%order(1, iord)  
     igv(1) = irod / gs%fftsize(2)  
     igv(2) = mod(irod, gs%fftsize(2))  
     igv(3) = gs%order(2, iord)  
     igv(4) = gs%order(3, iord)  
     igv(:) = mod(igv(:) + ffth(:), fftn(:)) - ffth(:)  
     gv(1:2) = real(igv(1:2), dp) + gs%rk(1:2)  
     do igv3 = igv(3), igv(4)                              ! loop over z axis
        gv(3) = real(igv3, dp) + gs%rk(3)  
        igs = igs + 1  
        grd(igs) = data(igs) * (crys%bvec(dir, 1) * gv(1) + &
             crys%bvec(dir, 2) * gv(2) + crys%bvec(dir, 3) * gv(3))
     end do
  end do

  return  

end subroutine derivative_gspace
