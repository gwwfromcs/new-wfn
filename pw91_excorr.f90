!     -*-Fortran-*-
!
subroutine pw91_excorr(ipr, ffts, crys, gs, inputcd, exc, afmagm)

  use all_to_all_module  
  include 'use.h'  
  implicit none  
  include 'interface.h'  
  include 'all_to_all.h'  
  !
  !     INPUT:
  !     ------
  !
  type(fft_struc), intent(in) :: ffts  ! for the fourier transform
  type(crystal), intent(in) :: crys    ! for the bvec array to take the gradient
  type(parallel_gspace), intent(in) :: &
       gs                              ! the gspace for the gradient computation
  integer, intent(in) :: ipr           ! printout flag
  !
  !     INPUT/OUTPUT:
  !     ------------
  !
  !     the charge density in fourier space on input, the exchange
  !     -correlation potential (Rydbergs) in fourier space on output
  !
  complex(dp), intent(inout) :: inputcd(gs%length, crys%nspin)  
  !
  !
  !     OUTPUT:
  !     ------
  !
  real(dp), intent(out) :: afmagm, &  ! the antiferromagnetic moment in mu_Bohr
       exc                            ! the exchange-correlation energy in Ryd
  !
  !     DESCRIPTION:
  !     -----------
  !
  !     Computes the exchange correlation potential a la Perdew-Wang 91
  !     as given in:
  !     Y.-M. Juan and E. Kaxiras, PRB 48, 14944 (1993),
  !     (in Rydberg) given the charge density (in electrons per unit cell)
  !     The potential replaces the density in the array inputcd.
  !     This version supports spin polarization.
  !
  !     The implementation follows:
  !     J.A. White and D.M. Bird PRB 50, 4954 (1994).
  !     i.e.:
  !
  !     Exc=V_cell/N_r Sum_r F_xc(n(r),|grad n(r)|)
  !
  !     Vxc(r)= dF_xc/dn(r) -  1/N_r Sum_{r',G} dF_xc/d|grad n(r')|
  !             iG*grad n(r')/|grad n(r')|*exp[iG*(r-r')]
  !
  !     To evaluate eq. (2) 8 FFTs are involved.
  !
  !
  !     1998 by Young-Gui Yoon. Based on an earlier version by Francesco
  !     Mauri, Bernd Pfrommer, and Alfredo Pasquarello.
  !
  !     --------------------- local variables -------------------------
  real(dp) :: alfr, alfden, alflog, dalfdrs, &
       zeta, zeta3, zeta4, fzeta, fzetap, fzetazeta4, &
       phi, phi2, phi3, phi4, invphi, invphi3, phip, &
       vcell, invcell, spscale, invspscale, roe, &
       grx, gry, grz, roeth, rs, gkf, sd, s, s2, s3, s4, ysr, ys, ysl, &
       sysl, fxexp, fxnum, fxden, fx, rs12, rs32, rs2, rs3, invrs12, td, &
       t, t2, t3, t4, aexp, abig, abig2, h0num, h0den, h0arg, h0, ccrnum, &
       ccrden, ccr, h1, dys, dfs, decd, ecr_dif, dh0da, dadec, dadz, &
       decdz, dh0d, dh0dt, dh0dg, dh0dz, dcdrs, dh1drs, dh1d, dh1dt, &
       dh1ds, dh1dg, dh1dz, byagr, dexcdg, normfac
  !
  !     _________________________________________________________________
  !
  integer :: i, k, ir, nspin, rsize, is, isign, norm  
  real(dp), dimension(1:3) :: aroe, agr, ecr, decdrs  
  real(dp), dimension(1:2) :: xchge, ecden, eclog, dfxd, dfxdg, dexcd
  real(dp), allocatable :: rhoe(:,:), &  ! charge density in rspace
       gradr(:,:,:)                      ! gradient of charge density in rspace
  complex(dp), allocatable :: &
       gswork(:), &   ! work array in gspace
       gswork1(:), &  ! work array in gspace
       gswork2(:), &  ! work array in gspace
       rswork(:)      ! work array in rspace
  !
  !     ------------------------------------------------------------------
  !
  !     pw91 parameters
  !
  !     The following three parameters are overwritten by exact expression
  !     Reference: J. Perdew et. al., PRB 46, 6671 (1992)
  !
  !     beta=(16/pi)*[(3*pi**2)**(1/3)]*0.004235 instead of 0.0667263212
  !     cc0=(16/pi)*[(3*pi**2)**(1/3)] instead of 15.75592
  !     cc1=0.004235-(3/7)*0.001667 instead of 0.003521
  !
  real(dp), parameter :: ax = -0.738558766382022405884230032680836d0
  real(dp), parameter :: gam = 0.5198420997897463295344212145565d0
  real(dp), parameter :: invgam = done / gam, invfzz = 9.0d0 * gam / 8.0d0
  real(dp), parameter :: eta = 1.0d-12, delt = 1.0d-12
  real(dp), parameter :: bb1 = 0.19645d0, bb2 = 0.27430d0, bb3 = -0.15084d0, &
       bb4 = 0.004d0, bb5 = 7.7956d0, alfa = 0.09d0, &
       beta = 0.06672632268006d0, cc0 = 15.75592034948314d0, &
       cc1 = 0.00352057142857d0, c1 = 0.001667d0, c2 = 0.002568d0, &
       c3 = 0.023266d0, c4 = 7.389d-6, c5 = 8.723d0, c6 = 0.472d0, &
       c7 = 7.389d-2, a = 0.0621814d0, alfa1 = 0.2137d0, bt1 = 7.5957d0, &
       bt2 = 3.5876d0, bt3 = 1.6382d0, bt4 = 0.49294d0, a_p = 0.0310907d0, &
       alfa1_p = 0.20548d0, bt1_p = 14.1189d0, bt2_p = 6.1977d0, &
       bt3_p = 3.3662d0, bt4_p = 0.62517d0, a_alf = 0.033774d0, &
       alfa1_alf = 0.11125d0, bt1_alf = 10.357d0, bt2_alf = 3.6231d0, &
       bt3_alf = 0.88026d0, bt4_alf = 0.49671d0
  !
  !     ------------------------------------------------------------------
  !
  ! spin polarization
  logical :: isp  
  !
  !     _________________________________________________________________
  real(dp), parameter :: x13 = dthird, x16m = dmhalf * dthird, &
       x23 = dtwo * dthird, x43 = dftrd, x76 = 7.0d0 / 6.0d0
  !     _________________________________________________________________
  !     derived parameters from pi
  !
  real(dp), parameter :: pisq = pi * pi
  real(dp) :: pider1, pider2, pider3
  !     _________________________________________________________________
  !     derived parameters from alfa and beta
  !
  real(dp), parameter :: abder1 = beta * beta / (dtwo * alfa), &  
       abder2 = done / abder1, abder3 = dtwo * alfa / beta  
  !
  pider1 = (0.75d0 / pi)**x13
  pider2 = (dthree * pisq)**x13
  pider3 = (dthree * pisq / 16.0d0)**x13

  norm = gs%fftsize(1) * gs%fftsize(2) * gs%fftsize(3)  
  normfac = crys%vcell / norm  
  vcell = crys%vcell  
  invcell = done / crys%vcell  
  !
  nspin = crys%nspin  
  spscale = real(nspin, dp)  
  invspscale = done / spscale  
  isp = .false.  
  if (nspin == 2) isp = .true.  
  !
  rsize = gs%r_size  
  exc = dzero
  !     _________________________________________________________________
  !     compute charge density in real space (rhoe) and gradient of the
  !     charge density in real space (gradr).
  !     rhoe(:,is) grad(:,:,is) coresponds to spin is; note that
  !     gradr(:,:,3) corresponds to gradient of total charge.
  !
  phi = done 
  invphi = done
  phi2 = done
  phi3 = done  
  invphi3 = done
  phi4 = done
  !
  allocate(rswork(rsize))  
  allocate(gswork(gs%length))  
  allocate(gswork1(gs%length))  
  allocate(gswork2(gs%length))  
  !
  allocate(rhoe(rsize, nspin))  
  allocate(gradr(rsize, 3, 2 * nspin - 1))  
  !
  do is = 1, nspin  
     gswork(:) = inputcd(:, is) * invcell  
     call fourier_transform(-1, ffts, gs, gswork(1), rswork(1), 1)
     rhoe(:, is) = real(rswork(:), dp)  
     do k = 1, 3  
        call derivative_gspace(k, inputcd(1, is), gs, crys, gswork(1))  
        gswork(:) = gswork(:) * invcell * zi
        call fourier_transform(-1, ffts, gs, gswork(1), rswork(1), 1)
        gradr(:, k, is) = real(rswork(:), dp)  
     end do
  end do
  !     _________________________________________________________________
  !     compute antiferromagnetic moment
  !
  afmagm = dzero
  if (isp) then  
     do i = 1, rsize  
        afmagm = afmagm + abs(rhoe(i, 1) - rhoe(i, 2))  
     end do
     call all_sum_all(afmagm)  
     if (norm > 0) then  
        afmagm = vcell * afmagm / real(norm, dp)          ! normalize
     else  
        print *, 'pw91_excorr: norm=0. Something is wrong.'  
        call mystop  
     end if
     !
     gradr(:,:, 3) = gradr(:,:, 1) + gradr(:,:, 2)  
     !
  end if
  !
  !     _________________________________________________________________
  !     main loop
  !
  do ir = 1, rsize  
     do is = 1, nspin  
        roe = rhoe(ir, is) * spscale  
        if ((.not. isp) .and. roe <= dzero) then  
           gradr(ir,:, 1) = dzero
           goto 100  
        end if
        if (isp .and. roe < dzero) then  
           gradr(ir, :,:) = dzero
           goto 100  
        end if
        aroe(is) = abs(roe)  
        grx = gradr(ir, 1, is) * spscale  
        gry = gradr(ir, 2, is) * spscale  
        grz = gradr(ir, 3, is) * spscale  
        agr(is) = sqrt(grx * grx + gry * gry + grz * grz)  
        roeth = aroe(is)**x13  
        if (isp .and. roe == dzero) then  
           xchge(is) = dzero
           dfxd(is) = dzero
           dfxdg(is) = dzero
        else  
           rs = pider1 / roeth  
           gkf = pider2 * roeth  
           sd = done / (dtwo * gkf * aroe(is))  
           s = agr(is) * sd  
           s2 = s * s  
           s3 = s * s2  
           s4 = s2 * s2  
           !  _________________________________________________________________
           !  exchange
           !
           ysr = sqrt(done + bb5 * bb5 * s2)  
           ys = bb5 * s + ysr  
           ysl = log(ys) * bb1  
           sysl = s * ysl  
           fxexp = exp(-1.0d2 * s2)  
           fxnum = done + sysl + (bb2 + bb3 * fxexp) * s2  
           fxden = done / (done + sysl + bb4 * s4)  
           fx = fxnum * fxden  
           xchge(is) = ax * fx * roeth  
           !  _________________________________________________________________
           !  first part xc-potential from exchange
           !
           dys = bb5 * (done + bb5 * s / ysr) / ys  
           dfs = -fxnum * (ysl + bb1 * s * dys + dfour * bb4 * s3) * &
                fxden * fxden + (ysl + bb1 * s * dys + dtwo * s * &
                (bb2 + bb3 * fxexp) - 2.0d2 * s3 * bb3 * fxexp) * fxden
           dfxd(is) = (ax * roeth * x43) * (fx - dfs * s)  
           dfxdg(is) = ax * roeth * dfs * sd  
        end if
     end do
     !  _________________________________________________________________
     !  back to total charge density
     !
     if (isp) then  
        roe = rhoe(ir, 1) + rhoe(ir, 2)  
        if (roe == dzero) then  
           rhoe(ir, :) = dzero
           gradr(ir, :,:) = dzero
           goto 100  
        end if
        aroe(3) = abs(roe)  
        grx = gradr(ir, 1, 3)  
        gry = gradr(ir, 2, 3)  
        grz = gradr(ir, 3, 3)  
        agr(3) = sqrt(grx * grx + gry * gry + grz * grz)  
        roeth = aroe(3)**x13  
        rs = pider1 / roeth  
        gkf = pider2 * roeth  
        sd = done / (dtwo * gkf * aroe(3))  
        s = agr(3) * sd  
        s2 = s * s  
     end if
     !  _________________________________________________________________
     !  correlation ecr=ec(rho,zeta)
     !
     rs12 = sqrt(rs)  
     invrs12 = done / rs12  
     rs32 = rs12 * rs  
     rs2 = rs * rs  
     rs3 = rs * rs2  
     ecden(1) = a * (bt1 * rs12 + bt2 * rs + bt3 * rs32 + bt4 * rs2)
     eclog(1) = log(done + (done / ecden(1)))  
     ecr(1) = -a * (done + alfa1 * rs) * eclog(1)  
     if (isp) then  
        ecden(2) = a_p * (bt1_p * rs12 + bt2_p * rs + bt3_p * rs32 + &
             bt4_p * rs2)
        eclog(2) = log(done + (done / ecden(2)))  
        ecr(2) = -a_p * (done + alfa1_p * rs) * eclog(2)  
        !
        alfden = a_alf * (bt1_alf * rs12 + bt2_alf * rs + bt3_alf * &
             rs32 + bt4_alf * rs2)
        alflog = log(done + (done / alfden))  
        alfr = a_alf * (done + alfa1_alf * rs) * alflog  
        !
        zeta = (rhoe(ir, 1) - rhoe(ir, 2)) / roe  
        if (abs(zeta) > done) zeta = real(sign(done, zeta), dp)  
        zeta3 = zeta * zeta * zeta  
        zeta4 = zeta3 * zeta  
        !           invgam=1.0/(2.d0**(3.d0/4.d0)-2.d0)
        fzeta = ((done + zeta)**x43 + (done - zeta)**x43 - dtwo) * invgam
        fzetap = x43 * ((done + zeta)**x13 - (done - zeta)**x13) * invgam
        fzetazeta4 = fzeta * zeta4  
        phi = ((done + zeta)**x23 + (done - zeta)**x23) * dhalf
        phip = (((done + zeta)**2 + eta)**x16m - &
             ((done - zeta)**2 + eta)**x16m) * x13
        invphi = done / phi  
        phi2 = phi * phi  
        phi3 = phi2 * phi  
        invphi3 = done / phi3  
        phi4 = phi2 * phi2  
        !
        ecr_dif = ecr(2) - ecr(1)  
        ecr(3) = ecr(1) + alfr * fzeta * invfzz * (done - zeta4) + &
             ecr_dif * fzetazeta4
        decdz = dfour * zeta3 * fzeta * (ecr_dif - alfr * invfzz) + &
             fzetap * (zeta4 * ecr_dif + (done - zeta4) * alfr * invfzz)
     end if
     !  _________________________________________________________________
     !  correlation h0(t,ecr,zeta)
     !
     td = pider3 * sd * invrs12 * invphi  
     t = agr(2 * nspin - 1) * td  
     t2 = t * t  
     t3 = t * t2  
     t4 = t2 * t2  
     aexp = exp(-abder2 * ecr(2 * nspin - 1) * invphi3) - done
     abig = abder3 / aexp  
     abig2 = abig * abig  
     h0num = t2 + abig * t4  
     h0den = done / (done + abig * t2 + abig2 * t4)  
     h0arg = done + abder3 * h0num * h0den  
     h0 = abder1 * log(h0arg) * phi3  
     !  _________________________________________________________________
     !  correlation h1(t,s,aroe,zeta)
     !
     ccrnum = c2 + c3 * rs + c4 * rs2  
     ccrden = done / (done + c5 * rs + c6 * rs2 + c7 * rs3)  
     ccr = c1 + ccrnum * ccrden  
     fxexp = exp (-1.0d2 * s2 * phi2)  
     h1 = phi3 * cc0 * (ccr - cc1) * t2 * fxexp  
     !  _________________________________________________________________
     !  updating of xc-energy
     !
     exc = exc + (ecr(2 * nspin - 1) + h0 + h1) * aroe(2 * nspin - 1)
     do is = 1, nspin  
        exc = exc + xchge(is) * aroe(is) * invspscale  
     end do
     !  _________________________________________________________________
     !  first part xc-potential from ecr
     !
     decdrs(1) = -a * alfa1 * eclog(1) * rs + a * (done + alfa1 * rs) * a * &
          (dhalf * bt1 * rs12 + bt2 * rs + dtrhalf * bt3 * rs32 + &
          dtwo * bt4 * rs2) / (ecden(1) * ecden(1) + ecden(1))
     if (isp) then  
        decdrs(2) = -a_p * alfa1_p * eclog(2) * rs + a_p * &
             (done + alfa1_p * rs) * a_p * (dhalf * bt1_p * rs12 + bt2_p * &
             rs + dtrhalf * bt3_p * rs32 + dtwo * bt4_p * rs2) / &
             (ecden(2) * ecden(2) + ecden(2))
        dalfdrs = a_alf * alfa1_alf * alflog * rs - a_alf * (done + &
             alfa1_alf * rs) * a_alf * (dhalf * bt1_alf * rs12 + bt2_alf * &
             rs + dtrhalf * bt3_alf * rs32 + dtwo * bt4_alf * rs2) / &
             (alfden * alfden + alfden)
        decdrs(3) = decdrs(1) * (done - fzetazeta4) + decdrs(2) * &
             fzetazeta4 + dalfdrs * fzeta * invfzz * (done - zeta4)
     end if
     decd = -x13 * decdrs(2 * nspin - 1)  
     !  _________________________________________________________________
     !  first part xc-potential from h0
     !
     dh0da = phi3 * abder1 / h0arg * abder3 * h0den * (t4 - h0num * &
          h0den * (t2 + dtwo * abig * t4))
     dadec = abder3 * abder2 * (aexp + done) / (aexp * aexp) * invphi3
     if (isp) dadz = dadec * (decdz - ecr(2 * nspin - 1) * dthree * &
          invphi * phip)
     dh0d = dh0da * dadec * decd  
     dh0dt = phi3 * abder1 / h0arg * abder3 * h0den * (dtwo * t + &
          dfour * abig * t3 - h0num * h0den * (dtwo * abig * t + dfour * &
          abig2 * t3))
     dh0d = dh0d - x76 * t * dh0dt  
     dh0dg = dh0dt * td  
     if (isp) dh0dz = dh0da * dadz + (dthree * h0 - dh0dt * t) * phip * invphi
     !  _________________________________________________________________
     !  first part xc-potential from h1
     !
     dcdrs = (c3 + dtwo * c4 * rs - ccrnum * ccrden * (c5 + dtwo * &
          c6 * rs + dthree * c7 * rs2)) * ccrden
     dh1drs = phi3 * cc0 * t2 * fxexp * dcdrs  
     dh1d = - x13 * rs * dh1drs  
     dh1dt = phi3 * dtwo * t * cc0 * (ccr - cc1) * fxexp  
     dh1d = dh1d - x76 * t * dh1dt  
     dh1ds = -2.0d2 * phi2 * s * cc0 * (ccr - cc1) * phi3 * t2 * fxexp
     dh1d = dh1d - x43 * s * dh1ds  
     dh1dg = dh1dt * td + dh1ds * sd  
     if (isp) dh1dz = (dthree * phip * invphi - 2.0d2 * phi * phip * s2) * h1
     !  _________________________________________________________________
     !  first part xc-potential
     !
     do is = 1, nspin  
        dexcd(is) = dfxd(is) + decd + dh0d + dh1d + ecr(2 * nspin - 1) &
             + h0 + h1
     end do
     isign = sign(done, agr(2 * nspin - 1) - delt)  
     byagr = dhalf * (1 + isign) / (agr(2 * nspin - 1) + (1 - isign) * delt)
     !
     if (isp) then  
        dexcdg = (dh0dg + dh1dg) * aroe(2 * nspin - 1) * byagr  
        do is = 1, nspin  
           rhoe(ir, is) = dexcd(is) - (zeta - dtwo * (dtrhalf - is)) * &
                (decdz + dh0dz + dh1dz)
           isign = sign(done, agr(is) - delt)  
           byagr = dhalf * (1 + isign) / (agr(is) + (1 - isign) * delt)
           gradr(ir, :, is) = spscale * gradr(ir, :, is) * dfxdg(is) * &
                aroe(is) * byagr
        end do
     else  
        rhoe(ir, 1) = dexcd(1)  
        dexcdg = (dfxdg(1) + dh0dg + dh1dg) * aroe(2 * nspin - 1) * byagr
     end if
     !
     gradr(ir, :, 2 * nspin - 1) = gradr(ir, :, 2 * nspin - 1) * dexcdg
100 end do
  !     _________________________________________________________________
  !     sum exc over all processors
  call all_sum_all(exc)  
  !     _________________________________________________________________
  !     au to ry
  !
  exc = exc * dtwo * normfac  
  !  _________________________________________________________________
  !  second part xc-potential: 3 forward ffts
  !
  !  gradr(:,k) now contains df_xc/d|grad n(r')| grad n(r')/|grad n(r')
  !
  do is = 1, nspin  
     gswork2(:) = zzero
     do k = 1, 3  
        if (isp) then  
           rswork(:) = gradr(:, k, 3) + gradr(:, k, is)  
        else  
           rswork(:) = gradr(:, k, 1)  
        end if
        call fourier_transform(1, ffts, gs, gswork(1), rswork(1), 1)
        call derivative_gspace(k, gswork(1), gs, crys, gswork1(1))  
        gswork2(:) = gswork2(:) + gswork1(:) * zi
     end do
     !
     !  _________________________________________________________________
     !  second part xc-potential: 1 inverse fft
     !
     call fourier_transform(-1, ffts, gs, gswork2(1), rswork(1), 1)
     rswork(:) = rhoe(:, is) - rswork(:)  
     !
     !     _________________________________________________________________
     !
     call fourier_transform(1, ffts, gs, inputcd(1, is), rswork(1), 1)
  end do
  !     au to ry
  !
  inputcd(:, :) = inputcd(:, :) * dtwo
  !     _________________________________________________________________
  deallocate(rswork)  
  deallocate(gswork)  
  deallocate(gswork1)  
  deallocate(gswork2)  
  !
  deallocate(rhoe)  
  deallocate(gradr)  
  !     _________________________________________________________________
  return

end subroutine pw91_excorr
!
!
!
!     _________________________________________________________________
!     -*-Fortran-*-
!
subroutine pbe_excorr(ipr, ffts, crys, gs, inputcd, exc, afmagm)

  use all_to_all_module  
  include 'use.h'  
  implicit none  
  include 'interface.h'  
  include 'all_to_all.h'  
  !
  !     INPUT:
  !     ------
  !
  type(fft_struc), intent(in) :: ffts  ! for the fourier transform
  type(crystal), intent(in) :: crys    ! for the bvec array to take the gradient
  type(parallel_gspace), intent(in) :: &
       gs                              ! the gspace for the gradient computation
  integer, intent(in) :: ipr           ! printout flag
  !
  !     INPUT/OUTPUT:
  !     ------------
  !
  !     the charge density in fourier space on input, the exchange
  !     -correlation potential (Rydbergs) in fourier space on output
  !
  complex(dp), intent(inout) :: inputcd(gs%length, crys%nspin)  
  !
  !
  !     OUTPUT:
  !     ------
  !
  real(dp), intent(out) :: afmagm, &  ! the antiferromagnetic moment in mu_Bohr
       exc                            ! the exchange-correlation energy in Ryd
  !
  !     DESCRIPTION:
  !     -----------
  !
  !     Computes the exchange correlation potential of  PBE  as given in
  !     J. Perdew, K. Burke, and M. Ernzerhorf, PRL 77, 3865 (1996) in
  !     Rydberg  given  the  charge  density  in electrons per unit cell.
  !     The potential  replaces the density in the array inputcd. This
  !     supports spin polarization.
  !
  !     The implementation follows:
  !     J.A. White and D.M. Bird PRB 50, 4954 (1994).
  !     i.e.:
  !
  !     Exc=V_cell/N_r Sum_r F_xc(n(r),|grad n(r)|)
  !
  !     Vxc(r)= dF_xc/dn(r) -  1/N_r Sum_{r',G} dF_xc/d|grad n(r')|
  !             iG*grad n(r')/|grad n(r')|*exp[iG*(r-r')]
  !
  !     To evaluate eq. (2) 8 FFTs are involved.
  !
  !
  !
  !     1998 by Young-Gui Yoon. Based on an earlier version of pw91_excorr
  !
  !     --------------------- local variables -------------------------
  real(dp) :: alfr, alfden, alflog, dalfdrs, ecr_dif, zeta, zeta3, zeta4, &
       fzeta, fzetap, fzetazeta4, phi, phi3, invphi, invphi3, phip, &
       vcell, invcell, spscale, invspscale, roe, grx, &
       gry, grz, roeth, rs, gkf, sd, s, s2, fxden, fx, rs12, rs32, rs2, &
       rs3, dfs, decd, dexcdg, byagr, td, t, t2, t3, t4, aexp, abig, &
       abig2, hnum, hden, harg, h, dhda, dadec, dadz, decdz, dhd, dhdt, &
       dhdg, dhdz, normfac
  !
  !     _________________________________________________________________
  !
  integer :: k, ir, nspin, rsize, is, isign, norm  
  real(dp), dimension(1:3) :: aroe, agr, ecr, decdrs  
  real(dp), dimension(1:2) ::xchge, ecden, eclog, dfxd, dfxdg, dexcd
  real(dp), allocatable :: rhoe(:,:), &  ! charge density in rspace
       gradr(:,:,:)                      ! gradient of charge density in rspace
  complex(dp), allocatable :: &
       gswork(:), &   ! work array in gspace
       gswork1(:), &  ! work array in gspace
       gswork2(:), &  ! work array in gspace
       rswork(:)      ! work array in rspace
  !
  !     ------------------------------------------------------------------
  !
  !     pbe parameters
  !
  !     ------------------------------------------------------------------
  !     Formulas:
  !        gam= 2^(4/3)-2
  !        fzz=f''(0)= 8/(9*gam)
  !        gamma=(1-log(2))/pi^2
  !        bet=coefficient in gradient expansion for correlation.
  !        eta=small number to stop d phi/ dzeta from blowing up.
  !        e_x[unif]=ax*rho^(4/3)  [LDA]
  !        ax = -0.75*(3/pi)^(1/3)
  !        e_x[PBE]=e_x[unif]*FxPBE(s)
  !        FxPBE(s)=1+uk-uk/(1+(um/uk)*s*s)
  !     ------------------------------------------------------------------
  !
  real(dp), parameter :: ax = -0.738558766382022405884230032680836d0
  real(dp), parameter :: um = 0.2195149727645171d0, uk = 0.804d0
  real(dp), parameter :: gam = 0.5198420997897463295344212145565d0
  real(dp), parameter :: invgam = done / gam, invfzz = 9.0d0 * gam / 8.0d0
  real(dp), parameter :: gamma = 0.03109069086965489503494086371273d0
  real(dp), parameter :: bet = 0.06672455060314922d0
  real(dp), parameter :: eta = 1.0d-12, delt = 1.0d-12
  real(dp), parameter :: a = 0.0621814d0, alfa1 = 0.2137d0, bt1 = 7.5957d0, &
       bt2 = 3.5876d0, bt3 = 1.6382d0, bt4 = 0.49294d0, a_p = 0.0310907d0, &
       alfa1_p = 0.20548d0, bt1_p = 14.1189d0, bt2_p = 6.1977d0, &
       bt3_p = 3.3662d0, bt4_p = 0.62517d0, a_alf = 0.033774d0, &
       alfa1_alf = 0.11125d0, bt1_alf = 10.357d0, bt2_alf = 3.6231d0, &
       bt3_alf = 0.88026d0, bt4_alf = 0.49671d0
  !
  logical :: isp             ! spin polarization
  !
  !     ------------------------------------------------------------------
  !
  real(dp), parameter :: x13 = dthird, x16m = dmhalf * dthird, &  
       x23 = 2.0d0 / 3.0d0, x43 = 4.0d0 / 3.0d0, x76 = 7.0d0 / 6.0d0  
  !     _________________________________________________________________
  !     derived parameters from pi
  !
  real(dp), parameter :: pisq = pi * pi
  real(dp) :: pider1, pider2, pider3
  !     _________________________________________________________________
  !     derived parameters from beta, gamma, and kappa
  !
  real(dp), parameter :: abder1 = gamma, abder2 = done / abder1, &  
       abder3 = bet * abder2, uk2 = uk * uk  
  !
  pider1 = (0.75d0 / pi)**x13
  pider2 = (3.0d0 * pisq)**x13
  pider3 = (3.0d0 * pisq / 16.0d0)**x13  

  norm = gs%fftsize(1) * gs%fftsize(2) * gs%fftsize(3)  
  normfac = crys%vcell / norm  
  vcell = crys%vcell  
  invcell = done / crys%vcell  
  !
  nspin = crys%nspin  
  spscale = real(nspin, dp)  
  invspscale = done / spscale  
  isp = .false.  
  if (nspin == 2) isp = .true.  
  !
  rsize = gs%r_size  
  exc = dzero
  !     _________________________________________________________________
  !     compute charge density in real space (rhoe) and gradient of the
  !     charge density in real space (gradr).
  !     rhoe(:,is) grad(:,:,is) coresponds to spin is; note that
  !     gradr(:,:,3) corresponds to gradient of total charge.
  !
  phi = done 
  invphi = done
  phi3 = done
  invphi3 = done
  !
  allocate(rswork(rsize))  
  allocate(gswork(gs%length))  
  allocate(gswork1(gs%length))  
  allocate(gswork2(gs%length))  
  !
  allocate(rhoe(rsize, nspin))  
  allocate(gradr(rsize, 3, 2 * nspin - 1))  
  !
  do is = 1, nspin  
     gswork(:) = inputcd(:, is) * invcell  
     call fourier_transform(-1, ffts, gs, gswork(1), rswork(1), 1)
     rhoe(:, is) = real(rswork(:), dp)  
     do k = 1, 3  
        call derivative_gspace(k, inputcd(1, is), gs, crys, gswork(1))  
        gswork(:) = gswork(:) * invcell * zi
        call fourier_transform(-1, ffts, gs, gswork(1), rswork(1), 1)
        gradr(:, k, is) = real(rswork(:), dp)  
     end do
  end do
  !     _________________________________________________________________
  !     compute antiferromagnetic moment and gradient of total charge
  !
  afmagm = dzero
  ! calculate antiferromagnetic moment
  if (isp) then  
     do ir = 1, rsize  
        afmagm = afmagm + abs(rhoe(ir, 1) - rhoe(ir, 2))  
     end do
     call all_sum_all(afmagm)  
     if (norm > 0) then  
        afmagm = normfac * afmagm          ! normalize
     else  
        print *, 'pw91_excorr: norm=0. Something is wrong.'  
        call mystop  
     end if
     !
     gradr(:,:, 3) = gradr(:,:, 1) + gradr(:,:, 2)  
     !
  end if
  !
  !     _________________________________________________________________
  !     main loop
  !
  do ir = 1, rsize  
     do is = 1, nspin  
        roe = rhoe(ir, is) * spscale          ! spin-scaled for exchange
        if ((.not. isp) .and. roe <= dzero) then  
           gradr(ir, :, 1) = dzero
           goto 100  
        end if
        if (isp .and. roe < dzero) then  
           gradr(ir, :,:) = dzero
           goto 100  
        end if
        aroe(is) = abs(roe)  
        grx = gradr(ir, 1, is) * spscale  
        gry = gradr(ir, 2, is) * spscale  
        grz = gradr(ir, 3, is) * spscale  
        agr(is) = sqrt(grx * grx + gry * gry + grz * grz)  
        roeth = aroe(is)**x13  
        if (isp .and. roe == dzero) then  
           xchge(is) = dzero
           dfxd(is) = dzero
           dfxdg(is) = dzero
        else  
           rs = pider1 / roeth  
           gkf = pider2 * roeth  
           sd = done / (dtwo * gkf * aroe(is))  
           s = agr(is) * sd  
           s2 = s * s  
           !  _________________________________________________________________
           !  exchange and first part xc-potential from exchange
           !           d[]d means (d[]/dn) for dfxd
           !           d[]d means n(d[]/dn) otherwise
           !           d[]dg means d[]/d|grad(n)|/n for dfxdg
           !           d[]dg means d[]/d|grad(n)| otherwise
           !           Be careful! fx for dfxd and dfxdg gives exchange energy
           !           when integrated; others are exchange enhancement factor.
           !
           fxden = done / (uk + um * s2)  
           fx = done + uk - uk2 * fxden  
           xchge(is) = ax * fx * roeth  
           !
           dfs = dtwo * uk2 * fxden * fxden * um * s  
           dfxd(is) = (ax * roeth * x43) * (fx - dfs * s)  
           dfxdg(is) = ax * roeth * dfs * sd  
        end if
     end do
     !  _________________________________________________________________
     !  back to total charge density
     !     Be careful! Total charge density is not spin-scaled!
     !
     if (isp) then  
        roe = rhoe(ir, 1) + rhoe(ir, 2)  
        if (roe == dzero) then  
           rhoe(ir, :) = dzero
           gradr(ir, :,:) = dzero
           goto 100  
        end if
        aroe(3) = abs(roe)  
        grx = gradr(ir, 1, 3)  
        gry = gradr(ir, 2, 3)  
        grz = gradr(ir, 3, 3)  
        agr(3) = sqrt(grx * grx + gry * gry + grz * grz)  
        roeth = aroe(3)**x13  
        rs = pider1 / roeth  
        gkf = pider2 * roeth  
        sd = done / (dtwo * gkf * aroe(3))  
        s = agr(3) * sd  
        s2 = s * s  
     end if
     !  _________________________________________________________________
     !  correlation ecr=ec(rho,zeta)
     !
     rs12 = sqrt(rs)  
     rs32 = rs12 * rs  
     rs2 = rs * rs  
     rs3 = rs * rs2  
     ecden(1) = a * (bt1 * rs12 + bt2 * rs + bt3 * rs32 + bt4 * rs2)
     eclog(1) = log(done + (done / ecden(1)))  
     ecr(1) = -a * (done + alfa1 * rs) * eclog(1)  
     if (isp) then  
        ecden(2) = a_p * (bt1_p * rs12 + bt2_p * rs + bt3_p * rs32 + &
             bt4_p * rs2)
        eclog(2) = log(done + (done / ecden(2)))  
        ecr(2) = -a_p * (done + alfa1_p * rs) * eclog(2)  
        !
        alfden = a_alf * (bt1_alf * rs12 + bt2_alf * rs + bt3_alf * &
             rs32 + bt4_alf * rs2)
        alflog = log(done + (done / alfden))  
        !yy note the sign change for alfr
        alfr = a_alf * (done + alfa1_alf * rs) * alflog  
        !
        zeta = (rhoe(ir, 1) - rhoe(ir, 2)) / roe  
        if (abs(zeta) > done) zeta = real(sign(done, zeta), dp)  
        zeta3 = zeta * zeta * zeta  
        zeta4 = zeta3 * zeta  
        !           invgam=1.0/(2.d0**(3.d0/4.d0)-2.d0)
        fzeta = ((done + zeta)**x43 + (done - zeta)**x43 - dtwo) * invgam
        fzetap = x43 * ((done + zeta)**x13 - (done - zeta)**x13) * invgam
        fzetazeta4 = fzeta * zeta4  
        phi = ((done + zeta)**x23 + (done - zeta)**x23) * dhalf
        phip = (((done + zeta)**2 + eta)**x16m - ((done - zeta)**2 + &
             eta)**x16m) * x13
        invphi = done / phi  
        phi3 = phi * phi * phi  
        invphi3 = done / phi3  
        !
        ecr_dif = ecr(2) - ecr(1)  
        ecr(3) = ecr(1) + alfr * fzeta * invfzz * (done - zeta4) + &
             ecr_dif * fzetazeta4
        decdz = dfour * zeta3 * fzeta * (ecr_dif - alfr * invfzz) + & 
             fzetap * (zeta4 * ecr_dif + (done - zeta4) * alfr * invfzz)
     end if
     !  _________________________________________________________________
     !  correlation h(t,ecr,zeta)
     !
     td = pider3 * sd / rs12 * invphi  
     t = agr(2 * nspin - 1) * td  
     t2 = t * t  
     t3 = t * t2  
     t4 = t2 * t2  
     aexp = exp(-abder2 * ecr(2 * nspin - 1) * invphi3) - done
     abig = abder3 / aexp  
     abig2 = abig * abig  
     hnum = t2 + abig * t4  
     hden = done / (done + abig * t2 + abig2 * t4)  
     harg = done + abder3 * hnum * hden  
     h = abder1 * log(harg) * phi3  
     !  _________________________________________________________________
     !  updating of xc-energy
     !        Be cafeful! Exchange energy has 1/2 factor in front;
     !        however, exchange potential do not have 1/2 factor if we
     !        just use the form of unpolarized form with twice the n(+,-).
     !
     exc = exc + (ecr(2 * nspin - 1) + h) * aroe(2 * nspin - 1)  
     do is = 1, nspin  
        exc = exc + xchge(is) * aroe(is) * invspscale  
     end do
     !  _________________________________________________________________
     !  first part xc-potential from ecr
     !
     decdrs(1) = -a * alfa1 * eclog(1) * rs + a * (done + alfa1 * &
          rs) * a * (dhalf * bt1 * rs12 + bt2 * rs + dtrhalf * bt3 * rs32 + &
          dtwo * bt4 * rs2) / (ecden(1) * ecden(1) + ecden(1))
     if (isp) then  
        decdrs(2) = -a_p * alfa1_p * eclog(2) * rs + a_p * &
             (done + alfa1_p * rs) * a_p * (dhalf * bt1_p * rs12 + bt2_p * &
             rs + dtrhalf * bt3_p * rs32 + dtwo * bt4_p * rs2) / &
             (ecden(2) * ecden(2) + ecden(2))
        dalfdrs = a_alf * alfa1_alf * alflog * rs - a_alf * (done + &
             alfa1_alf * rs) * a_alf * (dhalf * bt1_alf * rs12 + bt2_alf * &
             rs + dtrhalf * bt3_alf * rs32 + dtwo * bt4_alf * rs2) / &
             (alfden * alfden + alfden)
        decdrs(3) = decdrs(1) * (done - fzetazeta4) + decdrs(2) * &
             fzetazeta4 + dalfdrs * fzeta * invfzz * (done - zeta4)
     end if
     decd = -x13 * decdrs(2 * nspin - 1)  
     !  _________________________________________________________________
     !  first part xc-potential from h
     !
     dhda = phi3 * abder1 / harg * abder3 * hden * (t4 - hnum * &
          hden * (t2 + dtwo * abig * t4) )
     dadec = abder3 * abder2 * (aexp + done) / (aexp * aexp) * invphi3
     if (isp) dadz = dadec * (decdz - ecr(2 * nspin - 1) * dthree * &
          invphi * phip)
     dhd = dhda * dadec * decd  
     dhdt = phi3 * abder1 / harg * abder3 * hden * (dtwo * t + &
          dfour * abig * t3 - hnum * hden * (dtwo * abig * t + dfour * &
          abig2 * t3))
     dhd = dhd - x76 * t * dhdt  
     dhdg = dhdt * td  
     if (isp) dhdz = dhda * dadz + (dthree * h - dhdt * t) * phip * &
          invphi
     !  _________________________________________________________________
     !  first part xc-potential
     !        Note the following for the spin-polarized calculation.
     !        1. dexcd include dfxd: Exchange potential from the derivative
     !           w.r.t spin-scale*n(+,-) is included in the xc-potential.
     !           This is stored at rhoe(is) temporarily.
     !        2. dexcdg ( defined by d exc d |grad n|) at this point do not
     !           include the contribution from exchange: The dependence is
     !           through |grad n(+,-)| for exchange whereas the dependence
     !           is through |grad n_tot| with zeta for correlation of PBE.
     !        3. Therefore, contribution from the exchange is prepared at
     !           this stage and not included. Compare with unpolarized case.
     !
     do is = 1, nspin  
        dexcd(is) = dfxd(is) + decd + dhd + ecr(2 * nspin - 1) + h  
     end do
     !        if agr is less than delt byagr is set to zero.
     isign = sign(done, agr(2 * nspin - 1) - delt)  
     byagr = dhalf * (1 + isign) / (agr(2 * nspin - 1) + (1 - isign) * delt)
     !
     if (isp) then  
        dexcdg = dhdg * aroe(2 * nspin - 1) * byagr  
        do is = 1, nspin  
           rhoe(ir, is) = dexcd(is) - (zeta - dtwo * (dtrhalf - is)) * &
                (decdz + dhdz)
           isign = sign(done, agr(is) - delt)  
           byagr = dhalf * (1 + isign) / (agr(is) + (1 - isign) * delt)
           gradr(ir, :, is) = spscale * gradr(ir, :, is) * dfxdg(is) * &
                aroe(is) * byagr
        end do
     else  
        rhoe(ir, 1) = dexcd(1)  
        dexcdg = (dfxdg(1) + dhdg) * aroe(2 * nspin - 1) * byagr  
     end if
     !
     gradr(ir, :, 2 * nspin - 1) = gradr(ir, :, 2 * nspin - 1) * dexcdg
100 end do
  !     _________________________________________________________________
  !     sum exc over all processors
  call all_sum_all(exc)  
  !     _________________________________________________________________
  !     au to ry
  !
  exc = exc * dtwo * normfac  
  !     _________________________________________________________________
  !     second part xc-potential: 3 forward ffts
  !
  !     gradr(:,k) now contains df_xc/d|grad n(r')| grad n(r')/|grad n(r')
  !     for unpolarized case. For polarized case exchange is separate.
  !
  do is = 1, nspin  
     gswork2(:) = zzero
     do k = 1, 3  
        if (isp) then  
           rswork(:) = gradr(:, k, 3) + gradr(:, k, is)  
        else  
           rswork(:) = gradr(:, k, 1)  
        end if
        call fourier_transform(1, ffts, gs, gswork(1), rswork(1), 1)
        call derivative_gspace(k, gswork(1), gs, crys, gswork1(1))  
        gswork2(:) = gswork2(:) + gswork1(:) * zi
     end do
     !
     !     _________________________________________________________________
     !     second part xc-potential: 1 inverse fft
     !
     call fourier_transform( - 1, ffts, gs, gswork2(1), rswork(1), 1)
     rswork(:) = rhoe(:, is) - rswork(:)  
     !
     !     _________________________________________________________________
     !
     call fourier_transform(1, ffts, gs, inputcd(1, is), rswork(1), 1)
  end do
  !
  !     _________________________________________________________________
  !     au to ry
  !
  inputcd(:,:) = inputcd(:,:) * dtwo
  !     _________________________________________________________________
  deallocate(rswork)  
  deallocate(gswork)  
  deallocate(gswork1)  
  deallocate(gswork2)  
  !
  deallocate(rhoe)  
  deallocate(gradr)  
  !     _________________________________________________________________
  return  

end subroutine pbe_excorr
