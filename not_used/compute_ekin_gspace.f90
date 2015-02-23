!     @process extchk
!
subroutine compute_ekin_gspace(gs, crys)
 
  include 'use.h'  
  implicit none                ! implicit? Just say no!
  include 'interface.h'  
  !
  !      INPUT:
  !      ------
  !
  type(crystal), intent(in) :: crys  
  !
  !     recomputes the array ekin of the gspace for the value of gs%rk.
  !     assumes the ekin array is attached to the gspace, and has correct
  !     1996 Bernd Pfrommer, UC Berkeley
  !
  !     OUTPUT:     gs%ekin
  !     ------
  !
  type(parallel_gspace), intent(inout) :: gs  
  !
  !     ------------- local variables ------------------------------------
  !
  real(dp) :: qk(3)  
  !
  !     ----- variables for the gspace loop ---------
  !
  integer :: igv(4), fftn(4), ffth(4), igv3, irod, iord, igs  
  real(dp) :: gv(3)                                              ! if required
  !
  !     ----------------------------------------------
  !
  fftn(1:3) = gs%fftsize(1:3)  
  fftn(4) = gs%fftsize(3)  
  ffth(:) = fftn(:) / 2  
  igs = 0  
  do iord = 1, gs%lorder                             ! loop through x/y gspace
     irod = gs%order(1, iord)  
     igv(1) = irod / gs%fftsize(2)  
     igv(2) = mod(irod, gs%fftsize(2))  
     igv(3) = gs%order(2, iord)  
     igv(4) = gs%order(3, iord)  
     igv(:) = mod(igv(:) + ffth(:), fftn(:)) - ffth(:)  
     gv(1:2) = real(igv(1:2), dp)  
     do igv3 = igv(3), igv(4)                               ! loop over z axis
        gv(3) = real(igv3, dp)  
        igs = igs + 1  
        !
        !           compute kinetic energy here
        !
        qk(:) = gv(:) + gs%rk(:)  
        gs%ekin(igs) = (qk(1) * crys%bdot(1, 1) + &
             qk(2) * crys%bdot(2, 1) + qk(3) * crys%bdot(3, 1)) * qk(1) + &
             (qk(1) * crys%bdot(1, 2) + qk(2) * crys%bdot(2, 2) + &
             qk(3) * crys%bdot(3, 2)) * qk(2) + &
             (qk(1) * crys%bdot(1, 3) + qk(2) * crys%bdot(2, 3) + &
             qk(3) * crys%bdot(3, 3)) * qk(3)
     end do
  end do

  return  

end subroutine compute_ekin_gspace
