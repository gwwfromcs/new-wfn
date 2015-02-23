!@process extchk
!
subroutine add_to_j_bare(ffts, gs, psi1r, psi2, weight, rho)  

  include 'use.h'  
  implicit none                                   ! implicit? Just say no!
  include 'interface.h' 
  !
  !     INPUT:
  !     -----
  !
  type(fft_struc) :: ffts  

  type(parallel_gspace) :: gs         ! gspaces where wavefn were computed

  real(dp) :: weight              ! the weight of the kpoint, with right s

  complex(dp) :: psi1r(ffts%r_size), &       ! first function in realspace
       psi2(gs%length)                             ! psi2 in fourier space
  !
  !     INPUT/OUTPUT:
  !     -------------
  !
  complex(dp) :: rho(ffts%r_size)
  !  
  ! perform rho(r) = rho(r) + conjg(psi1(r)) * p
  !
  !
  !     1996 Bernd Pfrommer, UCB
  !
  !
  !     ---------------- local variables ---------------------------------

  integer :: i  

  complex(dp) :: cdsum, &
       psi2r(ffts%r_size)                  ! the fourier transform of psi2


  call fourier_transform(-1, ffts, gs, psi2, psi2r, 1)  
  !
  !     now add to charge with right weight
  !
  !      cdsum = zzero
  do i = 1, ffts%r_size  
     !         cdsum = cdsum + conjg(psi1r(i)) * psi1r(i)
     rho(i) = rho(i) + conjg(psi1r(i)) * psi2r(i) * weight  

  enddo
  !      write(9,*) 'length:', cdsum

  return  

end subroutine add_to_j_bare
