!-*-Fortran-*-
!*
subroutine setup_local_potential(ipr, ffts, gs, viondata, vlocfft)
  !
  include 'use.h'  
  implicit none                ! implicit? Just say no!
  include 'interface.h'  
  !
  !     INPUT:
  !     -----
  !
  integer, intent(in) :: ipr                 ! print flag
  type(fft_struc), intent(in) :: ffts        ! work arrays for the FFT
  type(parallel_gspace), intent(in) :: &
       gs                     ! the gspace for which the potential is set up
  complex(dp), intent(in) :: &
       viondata(*)            ! the potential in gspace viondata(gs%length)
  !
  !     OUTPUT:
  !     ------
  !
  complex(dp), intent(out) :: &
       vlocfft(*)             ! local pot in realspace vlocfft(ffts%r_size)
  !
  !
  !     Prepares the local part of the hamiltonian in real space: vlocfft
  !
  !
  !     --------- local variables ----------------------------------------
  !
  integer :: ierr, i  
  real(dp) :: dmax, dmin, cmax, cmin, warnmax  
  real(dp), parameter :: small = 1.0d-12 
  !
  !     ------------------------------------------------------------------
  !
  !     fourier transform the local potential. vlocfft contains local part
  !     of potential.
  !
  call fourier_transform(-1, ffts, gs, viondata, vlocfft, 1)  
  call check_array(vlocfft, ffts%r_size, dmax, dmin, cmax, cmin, &
       small, ierr)

  warnmax = max(abs(cmax), abs(cmin))  
  if (ipr > 0) then  
     if (ierr /= 0) call warn (271, warnmax, ierr, 'setup_local_potential')
     write(9, 104) dmax, dmin  
     write(9, 106) cmax, cmin  
     call myflush(9)  
  end if
  !
  return

104 FORMAT(/,'  MAX AND MIN OF POTENTIAL ',2F10.4)  
106 FORMAT('  MAX AND MIN COMPLEX COMPONENT   ',2g10.2/)  

end subroutine setup_local_potential
