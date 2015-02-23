!
subroutine nmr_shift(pw_params, crys, gs, chi_g)  
  !
  use all_to_all_module  
  include 'use.h'  
  implicit none           ! implicit? Just say no!
  include 'interface.h'  
  include 'all_to_all.h'  
  !
  !     computes NMR shift at the nuclei
  !
  !
  !     INPUT:
  !     -----
  !
  type(crystal), intent(in) :: crys            ! for the lattice vectors
  type(parallel_gspace), intent(in) :: gs      ! potential gspace
  type(pw_parameter), intent(in) :: pw_params  ! plane-wave parameters
  complex(dp), intent(in) :: &
       chi_g(gs%length, 3, 3)                  ! chi in gspace
  !
  !
  !     1996 Bernd Pfrommer, UCB
  !
  !
  !     ---------------- local variables ---------------------------------
  !
  complex(dp) :: ft(3, 3)  
  integer :: i, iat, nt, n  
  real(dp) :: r(3), t0
  real(dp), allocatable :: core_chi(:,:,:), ftbare(:,:,:)  
  real(dp), external :: gimmetime  
  !
  !     ----- variables for the gspace loop ---------
  !
  integer :: igv(4), fftn(4), ffth(4), igv3, irod, iord, igs  
  real(dp) :: gv(3)    ! if required

  t0 = gimmetime()  
  !
  !     first compute the correction due to the core dipoles
  !
  allocate(core_chi(3, 3, crys%mxdatm * crys%ntype))  
  allocate(ftbare(3, 3, crys%mxdatm * crys%ntype))  

  write(9, 50)  
  fftn(1:3) = gs%fftsize(1:3)  
  fftn(4) = gs%fftsize(3)  

  ffth(:) = fftn(:) / 2  
  n = 0  
  do nt = 1, crys%ntype  
     do iat = 1, crys%natom(nt)  
        n = n + 1  
        ft = 0  
        r = crys%rat(:, iat, nt)  
        igs = 0  
        do iord = 1, gs%lorder                      ! loop through x/y gspace
           irod = gs%order(1, iord)  
           igv(1) = irod / gs%fftsize(2)  
           igv(2) = mod(irod, gs%fftsize(2))  
           igv(3) = gs%order(2, iord)  
           igv(4) = gs%order(3, iord)  
           igv(:) = mod(igv(:) + ffth(:), fftn(:)) - ffth(:)  
           gv(1:2) = real(igv(1:2), dp) + gs%rk(1:2)  
           do igv3 = igv(3), igv(4)                 ! loop over z axis
              gv(3) = real(igv3, dp) + gs%rk(3)  
              igs = igs + 1  
              ft = ft + chi_g(igs, :,:) * exp(cmplx(dzero, &
                   r(1) * gv(1) + r(2) * gv(2) + r(3) * gv(3), dp))
           end do
        end do
        call all_sum_all(ft, 9)          ! 3x3 matrix sum across all processor
        ftbare(:,:, n) = ft  
        ft = -transpose(ft) / (137.036**2) * 1.0d6  
        write(9, 100) crys%nameat(nt), iat, crys%rat(:, iat, nt) / pi2, &
             real(ft(1, 1) + ft(2, 2) + ft(3, 3), dp) * dthird, real(ft, dp)
     end do
  end do

  if (iand(pw_params%output(1), 8) == 8) write(9, 940) gimmetime() - t0
  call ewald_dipole(pw_params%output(1), crys, gs, core_chi)  

  write(9, 60)  
  n = 0  
  do nt = 1, crys%ntype  
     do iat = 1, crys%natom(nt)  
        n = n + 1  
        ft = -core_chi(:,:, n) / (137.036**2) * 1.0d6  
        write(9, 100) crys%nameat(nt), iat, crys%rat(:, iat, nt) / pi2, &
             real(ft(1, 1) + ft(2, 2) + ft(3, 3), dp) * dthird, real(ft, dp)
     end do
  end do
  write(9, 55)  
  n = 0  
  do nt = 1, crys%ntype  
     do iat = 1, crys%natom(nt)  
        n = n + 1  
        ft = (-transpose(ftbare (:,:, n)) - core_chi(:,:, n)) / &
             (137.036**2) * 1.0d6
        write(9, 100) crys%nameat(nt), iat, crys%rat(:, iat, nt) / pi2, &
             real(ft(1, 1) + ft(2, 2) + ft(3, 3), dp) * dthird, real(ft, dp)
     end do
  end do

  call myflush(9)  
  deallocate(core_chi)  
  deallocate(ftbare)  

  return  

50 format(/' -------- BARE MAGNETIC SHIFTS -----------'//, &
       &     ' ATOM  NR.           POSITION',21x,'SHIFT [ppm] '/)
55 format(/' --- CORE-CORRECTED MAGNETIC SHIFTS ------'//, &
       &     ' ATOM  NR.           POSITION',21x,'SHIFT [ppm] '/)
60 format(/' ----------- CORE-CORRECTION -------------'//, &
       &     ' ATOM  NR.           POSITION',21x,'SHIFT [ppm] '/)

100 format(1x,a2,3x,i3,3f12.6,3x,f12.4,//3(9x,3f12.4/))  

940 format(' TIME FOR NMR FOURIER SUM:',f12.3)  

end subroutine nmr_shift
