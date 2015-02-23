!     @process extchk
!
subroutine charge2(ffts, pw_params, cd_gs, gs, bands, kpoints, &
     wavefn, den, energs,crys)

  use all_to_all_module  
  include 'use.h'  
  implicit none           ! implicit? Just say no!
  include 'interface.h'  
  include 'all_to_all.h'  
  !
  !     INPUT:
  !     -----
  !
  type(fft_struc) :: ffts  
  type(parallel_gspace) :: cd_gs, &    ! the gspace for the charge density and
       gs(*)                              ! gspaces where wavefn were computed
  type(pw_parameter) :: pw_params                 ! for the smearing parameter
  type(band) :: bands                        ! for nrk, occupation numbers ...
  type(kpoint) :: kpoints                          ! for the weights kpoints%w
  type(complex_gspace_array) :: wavefn  
  type(crystal), intent(in) :: crys            ! for vcell, nspin, ztot
  !
  !     OUTPUT:
  !     ------
  !
  type(energy) :: energs        ! efermi is input. compute eband,ektot
  complex(dp) :: den(cd_gs%length, bands%nspin)       ! charge density
  !
  !     computes the charge density and eband, ektot, from the eigenvector
  !     and eigenvalues
  !     1996 Bernd Pfrommer, UCB
  !
  !     ---------------- local variables ---------------------------------
  !
  real(dp), parameter :: small = 1.0d-12
  real(dp) :: t0, &
       dsmear, &                           ! gaussian smearing parameter [Ryd]
       dmin, dmax, cmin, cmax, sqrpi, spinf, delta,s,dasum
  real(dp), external :: gimmetime  
  complex(dp) :: cdsum, cdsumtot, cdsumtotr  
  integer :: is, &                                              ! spin counter
       ns, &                                                 ! number of stars
       ntot, &                     ! total number of points in realspace block
       irk, &                                                 ! kpoint counter
       neig, ndim, iel, &
       jmin, jmax, &
       i, j, k, &
       ierr, idum, ntottot,iwfn,nfn
  complex(dp), allocatable :: dhd(:), &                      ! realspace array
       rden(:)                                          ! accu for rspace dens
  !
  !     ------------------------------------------------------------------
  !
  !     transfer data from structures to local variables
  !

  t0 = gimmetime ()  
  ns = cd_gs%nstar  
  ntot = ffts%r_size               ! total number of points in realspace block
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
  allocate(dhd(ntot * pw_params%nbandsfft))  
  cdsumtot = zzero  
  cdsumtotr = zzero  
  if (bands%nspin >= 2) then
     spinf = dhalf
  else
     spinf = done 
  end if

  energs%eband = dzero
  energs%ektot = dzero

  do is = 1, bands%nspin        ! ------------   loop over spins ------------
     rden = zzero
     do irk=1, bands%nrk        ! ------------   loop over k points ---------

        neig = bands%nband(irk,is)
        ndim = gs(irk)%length
        !
        !        adds to sum of occupied eigenvalues and kinetic energy
        !
        do j = 1, neig
           energs%eband = energs%eband + &
                bands%occup(j, irk, is) * bands%energy(j, irk, is)
           energs%ektot = energs%ektot + bands%occup(j, irk, is) * &
                bands%ekn(j, irk, is)
        end do
        !
        !     -----------       find min and max band with nonzero occupation
        jmin = 0
        jmax = 0
        do i = 1, neig
           if (abs(bands%occup(i, irk, is)) >= small .and. jmin == 0) jmin = i
           if (abs(bands%occup(i, irk, is)) >= small .and. jmin /= 0) jmax = i
        end do

        if (jmin /= 0) then
           do j = jmin, jmax,pw_params%nbandsfft
              nfn=min(pw_params%nbandsfft, jmax - j+1)
              !
              !                 fourier transform to real space
              !
              call fourier_transform(-1, ffts, gs(irk), &
                   wavefn%data((j-1)*ndim+1, irk, is), dhd, nfn)
              !
              !                 square and add
              !
              do iwfn=0,nfn-1

! excluding Zn 3s and 3p
              do i = 5, ntot
                 rden(i) = rden(i) +dhd(i+iwfn*ntot)*conjg(dhd(i+iwfn*ntot))*&
                      bands%occup(j+iwfn, irk, is) * spinf
              end do
              end do
           end do
        end if
     end do                          ! ------------ end of loop over k points
     !
     !        fourier transform back to momentum space
     !
!     cdsum = zzero
!     do i = 1, ntot
!        cdsum = cdsum + rden(i)
!     end do
!
!     call all_sum_all(cdsum)
!     cdsum = cdsum / ntottot

     call fourier_transform(1, ffts, cd_gs, den(1, is), rden, 1)
     !
     !     symmetrize only if the star information are available
     !     Mauri Francesco
     !
! why is this commented out?
!     call symmetrize_scalar_global(den(1, is), cd_gs)

     call symmetrize_scalar_global(den(1, is), cd_gs)

  end do                                              ! end of loop over spins

  energs%eband = energs%eband / real(bands%nspin, dp)
  energs%ektot = energs%ektot / real(bands%nspin, dp)

  deallocate(dhd)
  deallocate(rden)

  if (iand(pw_params%output(1),8) == 8) write(9, 920) gimmetime() - t0

  return

920 format(' TIME FOR CHARGE DENSITY CALCULATION:',f12.3)

end subroutine charge2
