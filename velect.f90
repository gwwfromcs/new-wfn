!*
subroutine velect(ipr, gs, crys, ffts, energs, pwp, denc, vout, den, chdr)
  !
  !     1995/96 Bernd Pfrommer
  !
  use all_to_all_module  
  include 'use.h'  
  implicit none             ! implicit? Just say no!
  include 'interface.h'  
  include 'all_to_all.h'  
  include 'flibcalls.ph'  
  !     ---------------- arguments ---------------------------------------
  !
  !     INPUT:
  !     -----
  !      pwp%icorr       control parameter for correlation. it is
  !                  used only if iscr = 3 or 4. for explanation
  !                  see subroutine excorr
  !
  integer, intent(in) :: &
       ipr                     ! print flag
  type(parallel_gspace), intent(in) ::  gs  
  type(crystal), intent(in) :: crys  
  type(fft_struc), intent(in) :: ffts  
  type(pw_parameter), intent(in) :: pwp  
  complex(dp), intent(in) :: &
       denc(gs%length)         !  core charge density
  !
  !     OUTPUT:
  !     ------
  !
  type(energy), intent(out) :: energs  
  complex(dp), intent(out) :: &
       vout(gs%length, crys%nspin), &  ! the new V_hartree+V_xc on output
       den(gs%length, crys%nspin)      ! the valence charge density on output
  !
  !     WORK:
  !     ----
  !
  complex(dp) :: chdr(gs%r_size, crys%nspin)  
  !
  !     DESCRIPTION:
  !     -----------
  !
  !     * computes the xc energy  energs%xc
  !     * computes the V_xc(G=0) component for both spins
  !     * computes the spin-energy due to a magnetic field
  !     * computes the antiferromagnetic moment energs%afmagmom
  !
  !     ------------- local variables
  !
  complex(dp) :: td  
  real(dp) :: fac, tden(2), t0, t1, dum, &
       afmagm                 ! total antiferromagnetic moment
  integer :: i, j, k, &
       imagvec, &             ! index of  G of the magnetic field
       imagveci, &            ! index of -G
       is, n123, ntottot, my_ipr, isymmet, &
       invgmagnet(4)          ! g-vector for magnetic wave
  complex(dp) :: densmag(2, 2)  ! densmag(i,j) = spin i component of n(G),n(-G)
  real(dp), external :: gimmetime
  integer, external :: findvec  
  logical :: isp  
  !
  !     ------------------------------------------------------------------
  !
  t0 = gimmetime()
  isp = .false.  
  if (crys%nspin == 2) isp = .true.  

  ntottot = gs%r_size  
  call all_sum_all(ntottot)  
  !     ---- add core charge density ------------
  if (isp) then  
     vout(1:gs%length, 1:crys%nspin) = zzero
     call mzaxpy(gs%length, zhalf, denc(1), 1, vout(1, 1), 1)
                                                      ! take HALF of core cd
     call mzaxpy(gs%length, zhalf, denc(1), 1, vout(1, 2), 1)
     call mzaxpy(2 * gs%length, zone, den(1, 1), 1, vout(1, 1), 1)
                                                      ! add valence charge
  else  
     call mzcopy(gs%length, denc(1), 1, vout(1, 1), 1)  
     call mzaxpy(gs%length, zone, den(1, 1), 1, vout(1, 1), 1)
  end if
  !
  !     find charge components corresponding to staggered magnetic field
  !
  if (isp .and. pwp%mfield%gvec(4) == 1) then  
     imagvec = findvec(pwp%mfield%gvec, gs)  
     if (imagvec >= 1) then  
        densmag(1, 1) = vout(imagvec, 1)  
        densmag(2, 1) = vout(imagvec, 2)  
     end if

     invgmagnet(1:3) = -pwp%mfield%gvec(1:3)  
     !        cos(Gr) = 0.5*(exp(iGr) + exp(-iGr)), so have to take inverse a
     imagveci = findvec(invgmagnet, gs)  
     if (imagveci >= 1) then  
        densmag(1, 2) = vout(imagveci, 1)  
        densmag(2, 2) = vout(imagveci, 2)  
     end if
  end if
  !
  !      compute exchange and correlation potential/energy
  !
  !
  my_ipr = ipr  
  if (pwp%icorr /= 'pw' .and. pwp%icorr /= 'pb') then  
     !
     !        inverse fourier transform the charge density to real space
     !
     do is = 1, crys%nspin  
        call fourier_transform(-1, ffts, gs, vout(1, is), chdr(1, is), 1)
     end do

     call excorr(my_ipr, pwp%icorr, crys%nspin, gs%r_size, &
          gs%fftsize(1) * gs%fftsize(2) * gs%fftsize(3), chdr(1, 1), &
          crys%vcell, energs%xc, afmagm)
     do is = 1, crys%nspin  
        call fourier_transform(1, ffts, gs, vout(1, is), chdr(1, is), 1)
     end do
  else  
     if (pwp%icorr == 'pw') then  
        call pw91_excorr(my_ipr, ffts, crys, gs, vout(1, 1), energs%xc, afmagm)
     else
        call pbe_excorr(my_ipr, ffts, crys, gs, vout(1, 1), energs%xc, afmagm)
     end if
  end if
  !
  !      compute and add v hartree to the exchange correlation potential
  !
  fac = 8.0d0 * pi / crys%vcell  

  energs%vxc0 = dzero
  if (isp) then  
     do i = 1, gs%length  
        if (gs%ekin(i) > dzero) then  
           td = fac * (den(i, 1) + den(i, 2)) / gs%ekin(i)  
           vout(i, 1) = vout(i, 1) + td  
           vout(i, 2) = vout(i, 2) + td  
        else  
           ! save G=0 of xc pot
           energs%vxc0(1) = real(vout(i, 1), dp)  
           energs%vxc0(2) = real(vout(i, 2), dp)  
        end if
     end do
  else  
     do i = 1, gs%length  
        if (gs%ekin(i) > dzero) then  
           vout(i, 1) = vout(i, 1) + fac * den(i, 1) / gs%ekin(i)  
        else  
           ! save G=0 of xc pot
           energs%vxc0(1) = real(vout(i, 1), dp)  
        end if
     end do
  end if
  !
  !     add staggered magnetic field if necessary
  !
  if (isp .and. pwp%mfield%gvec(4) == 1) then  
     if (imagvec >= 1) then  
        vout(imagvec, 1) = vout(imagvec, 1) - dhalf * &
             cmplx(pwp%mfield%h, dzero, dp)               ! v=-mu * H
        vout(imagvec, 2) = vout(imagvec, 2) + dhalf * &
             cmplx(pwp%mfield%h, dzero, dp)               ! v=+mu * H
     end if
     !        cos(Gr) = exp(iGr) + exp(-iGr), so have to take inverse as well
     if (imagveci >= 1) then  
        vout(imagveci, 1) = vout(imagveci, 1) - dhalf * &
             cmplx(pwp%mfield%h, dzero, dp)               ! v=-mu * H
        vout(imagveci, 2) = vout(imagveci, 2) + dhalf * &
             cmplx(pwp%mfield%h, dzero, dp)               ! v=+mu * H
     end if
     energs%magnetic = dmhalf * pwp%mfield%h * abs((densmag(1, 1) - &
          densmag(2, 1) + densmag(1, 2) - densmag(2, 2)))
     !
     !        must symmetrize afterwards
     !
     call symmetrize_scalar_global(vout(1, 1), gs)  
     call symmetrize_scalar_global(vout(1, 2), gs)  
  else  
     energs%magnetic = dzero
  end if

  call all_sum_all(energs%vxc0, crys%nspin)  
  energs%afmagmom = afmagm  

  if (ipr >= 2 .and. isp) write(9, 300) afmagm  

  if (ipr > 0 .and. iand(pwp%output(1), 8) == 8) write(9, 930) gimmetime() - t0 

  return

930 format(' TIME FOR VELECT:',f12.3)  
300 format(' TOTAL ANTIFERROMAGNETIC MOMENT:',f12.6)  

end subroutine velect
