!     @process extchk
!
subroutine kinet(ham, xvec, neig, ekn)

  use all_to_all_module  
  include 'use.h'  
  implicit none             ! implicit? Just say no!
  include 'interface.h'  
  include 'all_to_all.h'  
  !
  !     calculates the kinetic energy
  !     written february 18 1990. jlm
  !
  !     parallel fortran 90 version 1996 Bernd Pfrommer, UC Berkeley
  !
  !
  !     input:
  !     ham         hamiltonian information
  !     neig        number of eigenvectors
  !     xvec(j,i)  real part of the component j of eigenvector i
  !
  !     output:
  !     ekn(i)      kinetic energy of eigenvector i. (rydberg)
  !
  !     INPUT:
  !     -----
  !
  integer, intent(in) :: neig                ! number of eigenvectors computed
  type(hamiltonian), intent(in) :: ham ! to access the gspace and the ekinmods
  complex(dp), intent(in) :: &
       xvec(ham%dim, neig)                 ! eigenvecs as they come from diag.
  !
  !     OUTPUT:
  !     ------
  !
  real(dp), intent(out) :: ekn(neig)             ! kinetic energy of the bands
  !
  !     ----------- local variables --------------------------------------
  !
  real(dp) :: t 
  real(dp), external :: myderf  
  integer :: i, j, len  
  !
  len = ham%gspace%length  
  do i = 1, neig  
     ekn(i) = dzero  
     do j = 1, len  
        t = ham%gspace%ekin(j) + ham%ekinmod(1) * (1 + &
             myderf((ham%gspace%ekin(j) - ham%ekinmod(2)) / ham%ekinmod(3)))
        ekn(i) = ekn(i) + t * xvec(j, i) * conjg(xvec(j, i))  
     end do
  end do

  call all_sum_all(ekn, neig)  

  return

100 format(' KINETIC ENERGY OF BANDS (RYD): ', &
       &      /' -----------------------------')

110 format(7f10.4)  

end subroutine kinet
