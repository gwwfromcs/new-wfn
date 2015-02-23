!
subroutine momentum_density(crys, kp, bands, k_gs, wfn)  

  use all_to_all_module  
  include 'use.h'  
  implicit none              ! implicit? Just say no!
  include 'interface.h'  
  include 'all_to_all.h'  
  !
  !     INPUT:
  !     -----
  !
  type(crystal), intent(in) :: crys  
  type(kpoint), intent(in) :: kp  
  type(band), intent(in) :: bands  
  type(parallel_gspace), intent(in) :: k_gs(kp%nrk)  
  type(complex_gspace_array), intent(in) :: wfn  
  !
  !     1997 Bernd Pfrommer. Based on a routine by Balazs Kralik.
  !
  !
  !
  !     DESCRIPTION:
  !     -----------
  !
  !     Computes the momentum density along special lines in the BZ.
  !     Those are given by the kpoint mesh for which the wave functions
  !     are computed. Make sure the occupation numbers are correct, i.e.
  !     between 0 and 1 with spin, or 0 and 2 without spin.
  !
  !
  !
  !
  !     ------------------- local variables ----------------------
  !
  real(dp) :: rk(3), rkc(3), rkm(3), klensq, dotp, ecut, ekin, &
       gvec(3), occp
  real(dp), allocatable :: momdens(:)  
  integer :: ndim, n, irk, is  
  !
  !     variables for the gspace loop
  !
  integer :: igv(4), fftn(4), ffth(4), igv3, irod, iord, igs  
  !
  !     -----------------------------------------------------------
  !
  ! does not make sense, don't use

  ecut = 2.0d10  
  !
  !     do sanity checks
  !
  if (k_gs(1)%myproc > 1) then  
     write(9, *) 'momentum_density works only SERIAL!'  
     call mystop  
  end if

  allocate(momdens(bands%max))  
  do is = 1, crys%nspin  
     write(9, *) '========= SPIN:', is, '=========='  
     do irk = 1, kp%nrk  
        rk = kp%rk(:, irk)  
        rkc = matmul(crys%bvec, rk)  
        rkm = matmul(crys%bdot, rk)  
        klensq = rkc(1) * rkc(1) + rkc(2) * rkc(2) + rkc(3) * rkc(3)  
        !
        !           loop over gspace to get the desired components
        !
        fftn(1:3) = k_gs(irk)%fftsize(1:3)  
        fftn(4) = k_gs(irk)%fftsize(3)  

        ffth(:) = fftn(:) / 2  
        igs = 0  
        ndim = k_gs(irk)%length  
        do iord = 1, k_gs(irk)%lorder          ! loop through x/y gspace
           irod = k_gs(irk)%order(1, iord)  
           igv(1) = irod / k_gs(irk)%fftsize(2)  
           igv(2) = mod(irod, k_gs(irk)%fftsize(2))  
           igv(3) = k_gs(irk)%order(2, iord)  
           igv(4) = k_gs(irk)%order(3, iord)  
           igv(:) = mod(igv(:) + ffth(:), fftn(:)) - ffth(:)  
           do igv3 = igv(3), igv(4)             ! loop over z axis
              igs = igs + 1  
              ekin = k_gs(irk)%ekin(igs)  
              gvec(1) = real(igv(1), dp) + rk(1)  
              gvec(2) = real(igv(2), dp) + rk(2)  
              gvec(3) = real(igv3, dp) + rk(3)  
              !
              !                 check if this gvector is parallel to kpoint
              !
              dotp = gvec(1) * rkm(1) + gvec(2) * rkm(2) + gvec(3) * rkm(3)  
              if ((klensq /= 0 .and. dotp * dotp + 1.0d-8 > ekin * klensq &
                   .and. ekin <= ecut) .or. (klensq == 0 .and. ekin == dzero)) &
                   then
                 !
                 !                    sum across all bands:
                 !
                 momdens = dzero
                 do n = 1, bands%nband(irk, is)  
                    occp = bands%occup(n, irk, is) / (kp%w(irk) * &
                         real(crys%nspin, dp))
                    momdens(n) = momdens(n) + occp * &
                         abs(wfn%data(ndim * (n - 1) + igs, irk, is))**2
                 end do
                 write(9, '(3f12.6,200g18.6)') gvec, &
                      sum(momdens(1:bands%nband(irk, is))), &
                      momdens(1:bands%nband(irk, is))
              end if
           end do
        end do
     end do
  end do

  deallocate(momdens)  

end subroutine momentum_density
