!     @process extchk
!
subroutine resolved_charge(energy_window, ffts, pw_params, cd_gs, &
     gs, bands, kpoints, wavefn, den, energs)
  !
  use all_to_all_module  
  include 'use.h'  
  implicit none              ! implicit? Just say no!
  include 'interface.h'  
  include 'all_to_all.h'  
  !
  !     INPUT:
  !     -----
  !
  type(fft_struc), intent(in) :: ffts  
  type(parallel_gspace), intent(in) :: &
       cd_gs, &                    ! the gspace for the charge density and
       gs(*)                       ! gspaces where wavefn were computed
  type(pw_parameter), intent(in) :: pw_params    ! for the smearing parameter
  type(band), intent(in) :: bands  ! for nrk, occupation numbers ...
  type(kpoint), intent(in) :: kpoints    ! for the weights kpoints%w
  type(complex_gspace_array), intent(in) :: wavefn  
  type(energy), intent(in) :: energs     ! needed for the vxc0 shift
  real(dp), intent(in) :: &
       energy_window(2)      ! energy window as start=(1), end=(2) in [eV]
  !
  !     OUTPUT:
  !     ------
  !
  complex(dp), intent(out) :: &
       den(cd_gs%length, bands%nspin)    ! charge density
  !
  !     DESCRIPTION:
  !     -----------
  !
  !     computes the energy-resolved charge density which is generated
  !     by wavefunctions with eigenvalues inside the energy window
  !
  !     the weights are computed as:
  !
  !                              eps_upper
  !     f(eps) = 1/sqrt(pi) *   /             exp(-((eps-e)/sigma)^2) de
  !                            eps_lower
  !
  !     with the definition t=eps/sigma, this becomes:
  !
  !     f(t) = 1 - 0.5*erfc(t-t_lower) - 0.5* erfc(t_upper-t)
  !
  !
  !
  !     1996 Bernd Pfrommer, UCB
  !     ---------------- local variables ---------------------------------
  !
  real(dp) :: t0, erfc1c, erfc2c, occp, desm, &
       dsmear, &                 ! gaussian smearing parameter [Ryd]
       dmin, dmax, cmin, cmax, &
       sqrpi, spinf, delta
  real(dp), external :: gimmetime, myderfc  
  complex(dp) :: cdsum, cdsumtot, cdsumtotr  
  real(dp), parameter :: small = 1.0d-12, range = dfour
  integer :: is, &    ! spin counter
       ns, &          ! number of stars
       ntot, &        ! total number of points in realspace block
       irk, &         ! kpoint counter
       neig, ndim, iel, i, j, k, ierr, idum, ntottot
  complex(dp), allocatable :: &
       dhd(:), &  ! realspace array
       rden(:)    ! accu for rspace dens
  !
  !     ------------------------------------------------------------------
  !
  !     transfer data from structures to local variables
  !
  t0 = gimmetime()  
  ns = cd_gs%nstar  
  ntot = ffts%r_size              ! total number of points in realspace block
  dsmear = pw_params%smearing / ryd  
  sqrpi = sqrt(pi)  
  !
  !     get dimensions of realspace block
  !
  ntottot = ntot  
  call all_sum_all(ntottot)  
  !
  !     allocate various charge densities
  !
  allocate(rden(ntot))  
  allocate(dhd(ntot))  
  cdsumtot = zzero  
  cdsumtotr = zzero  

  write(9, 100) energy_window(1), energy_window(2)  
  if (bands%nspin >= 2) then
     spinf = dhalf
  else
     spinf = done  
  end if

  do is = 1, bands%nspin  ! loop over spins
     rden = zzero
     do irk = 1, bands%nrk  ! ------------   loop over k points ----
        neig = bands%nband(irk, is)
        ndim = gs(irk)%length
        do j = 1, neig
           !
           !              compute occupation number
           !
           !              first treat lower boundary
           !
           desm = (bands%energy(j, irk, is) + energs%evalshift(is) - &
                energy_window(1) / ryd) / max(dsmear, 1.0d-6)
           if (desm >= range) then  ! well above lower margin
              erfc1c = dzero
           else if (desm < -range) then  ! well below lower margin
              erfc1c = done
           else
              erfc1c = dhalf * myderfc(desm)
           end if
           !
           !              now treat upper bound
           !
           desm = (energy_window(2) / ryd - bands%energy(j, irk, is) - &
                energs%evalshift(is)) / max(dsmear, 1.0d-6)
           if (desm >= range) then  ! well above lower margin
              erfc2c = dzero
           else if (desm < -range) then  ! well below lower margin
              erfc2c = done
           else
              erfc2c = dhalf * myderfc(desm)
           end if
           occp = done - erfc1c - erfc2c
           occp = occp * kpoints%w(irk)
           !
           !              fourier transform to real space
           !
           call fourier_transform(-1, ffts, &
                gs(irk), wavefn%data((j - 1) * ndim + 1, irk, is), dhd, 1)
           !
           !              square and add
           !
           do i = 1, ntot
              rden(i) = rden(i) + dhd(i) * conjg(dhd(i)) * occp
           end do
        end do
     end do  ! ------------ end of loop over k points
     !
     !        fourier transform back to momentum space
     !
     cdsum = zzero
     do i = 1, ntot
        cdsum = cdsum + rden(i)
     end do

     call all_sum_all(cdsum)
     cdsum = cdsum / ntottot

     call fourier_transform(1, ffts, cd_gs, den(1, is), rden, 1)
     !
     !        symmetrize the charge density
     !        only if cd_gs%istar .eq. true
     !        Francesco Mauri
     !
     call symmetrize_scalar_global(den(1, is), cd_gs)
  end do  ! end of loop over spins

  deallocate(dhd)
  deallocate(rden)

  if (iand(pw_params%output(1),8) == 8) write(9, 920) gimmetime() - t0

  return

100 format(/' COMPUTING CHARGE DENSITY WITH ENERGY WINDOW:', &
       f10.2, ' TO ', f10.2, 'eV')

920 format(' TIME FOR ENERGY-RESOLVED CHARGE', &
       ' DENSITY CALCULATION:',f12.3)

end subroutine resolved_charge
