!
subroutine layout_scalapack(matsize, b, nproc, p, q)  
  !
  use constants
  implicit none
  !
  !     INPUT:
  !     -----
  !
  integer, intent(in) :: nproc, &        ! number of processors
       matsize                           ! size of matrix
  !
  !     OUTPUT:
  !     ------
  !
  integer, intent(out) :: b, &           ! block size
       p, &                              ! processor grid row
       q                                 ! processor grid column
  !
  !     DESCRIPTION:
  !     -----------
  !
  !     Figures out a  p by q processor grid layout for the scalapack
  !     library. This p*q grid is used to partition the matrix with
  !     a block size b. For more, see the scalapack documentation
  !
  !     The goal is to get a processor grid which is as close to
  !     "square" as possible.
  !
  !     1997 Andrew Canning, Bernd Pfrommer
  !
  !     ----------------------- local variables ----------------------
  !
  integer :: i  
  !
  !     Find processor grid
  !
  p = int(sqrt(real(nproc, dp) + 1.0d-6))  
  do i = p, 1, -1  
     if (mod(nproc, i) == 0) exit  
  end do

  p = i  
  q = nproc / p  
  !
  !     now for the block size
  !
  b = min(32, matsize / (max(p, q)))  
  !
  !     ensure nonzero
  !
  b = max(b, 1)  

  return  

end subroutine layout_scalapack
