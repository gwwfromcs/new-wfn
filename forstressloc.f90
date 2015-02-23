!
subroutine forstressloc(ipr, ispn, fsum, ssum, aadot, exc, &
     ealpha, ng, gsvec, ek, vion, vin, vout, den, denc, dvql, ddc, vql, &
     dnc, ntype, natom, rat, bdot, vcell, mxdatm)
  !
  !     computes the local potential contributions to the Hellman-Feynman
  !     Forces and the stress. Adapted from a code by J.L. Martins 1988
  !
  !     parallel version 1995 Bernd Pfrommer
  !
  !
  use constants
  implicit none             ! never do implicit !
  !
  !     INPUT:
  !     -----
  !
  integer, intent(in) :: &
       ipr, &                 ! print flag
       ispn, &                ! number of spin components: 1 or 2
       ng, &                  ! number of gvectors
       ntype, &               ! number of atomic types
       gsvec(3, *), &         ! gvectors
       natom(*), &            ! number of atoms of each type
       mxdatm                 ! dimensions for arrays
  real(dp), intent(in) :: &
       vcell, &               ! unit cell volume
       bdot(3, 3), &          ! reciprocal lattice metric
       vql(ng, ntype), &      ! local pseudopotential of different atoms
       dnc(ng, ntype), &      ! core charge density of atoms
       rat(3, mxdatm, *), &   ! atomic coordinates in lattice *2*pi
       exc, &                 ! the exchange correlation energy
       ealpha, &              ! the divergent term to the energy contribution
       ek(ng)                 ! kinetic energy of gvectors
  complex(dp), intent(in) :: &
       vion(ng, ispn), &      ! the screened ionic potential
       vin(ng, ispn), &       ! electronic V_hxc input potential
       vout(ng, ispn), &      ! electronic V_hxc output potential
       den(ng, ispn), &       ! valence charge density
       denc(ng), &            ! core charge density
       dvql(ng), &            ! derivative of the local pseudopotential
       ddc(ng)                ! derivative of core charge density
  !
  !     OUTPUT:
  !     ------
  !
  real(dp), intent(out) :: &
       fsum(3, mxdatm, *), &  ! unsymmetrized forces[lattice coordinates]
                              ! from local part on each atom
       aadot(6), &            ! metric in realspace:
                              ! adot(11,22,33,12,23,31)/(2*pi)**2
       ssum(6)                ! unsymmetrized, metric stress tensor [Ryd]
  !
  !     ---------------- local variables ---------------------------------
  !
  integer :: i, j, k, l, is  
  real(dp) :: fi, fac, vld, dvld, &
       vxd, dvxd, exp1, exp2, vhd, denr, deni
  complex(dp) :: vdg, vd, vxc, densum, dentot  
  logical :: isp  

  isp = .false.  
  if (ispn == 2) isp = .true.  
  !
  !     initialize fsum
  !
  do i = 1, ntype  
     do j = 1, natom(i)  
        do k = 1, 3  
           fsum(k, j, i) = dzero  
        end do
     end do
  end do

  fac = dfour * pi2 / vcell  
  !     -----------------------------------------------------------------
  !              first the local contribution to the forces
  !

  do i = 1, ntype                                ! loop over atomic types
     do j = 1, ng                       ! loop over gvectors, exclude G=0
        if (ek(j) > dzero) then  
           !
           !           vql * cc(den) + cc(vxc) * denc
           !           vxc = vout - vhartree
           !
           !    spinpolarized:  (vxc_up-vhartree + vxc_down-vhartree) * dnc/2
           !
           !
           if (isp) then  
              dentot = den(j, 1) + den(j, 2)  
              vdg = vql(j, i) * conjg(dentot) + &
                   conjg((vout(j, 1) + vout(j, 2)) * dhalf - &
                   fac * dentot / ek(j)) * dnc (j, i)
           else  
              vdg = vql(j, i) * conjg(den(j, 1)) + &
                   conjg(vout(j, 1) - fac * den(j, 1) / ek(j)) * dnc(j, i)
           end if
           do l = 1, natom (i)             ! loop over atoms of same type
              fi = real(gsvec(1, j), dp) * rat(1, l, i) + &
                   real(gsvec(2, j), dp) * rat(2, l, i) + &
                   real(gsvec(3, j), dp) * rat(3, l, i)
              exp1 = sin(fi) * real(vdg, dp) - cos(fi) * aimag(vdg)  
              !
              !                 add to forces
              !
              fsum(1, l, i) = fsum(1, l, i) + real(gsvec(1, j), dp) * exp1  
              fsum(2, l, i) = fsum(2, l, i) + real(gsvec(2, j), dp) * exp1  
              fsum(3, l, i) = fsum(3, l, i) + real(gsvec(3, j), dp) * exp1  
           end do
        end if
     end do
  end do
  !     ------------------------------------------------------------------
  !                            local part of the stress
  !
  !     initialize ssum array
  !
  do i = 1, 6  
     ssum(i) = dzero  
  end do
  !
  !     compute metric in real space
  !     aadot(1-6) = adot(11, 22, 33, 12, 23, 31) / (2*pi)**2
  !
  fac = vcell / (pi2 * pi2 * pi2)  
  fac = fac * fac  
  aadot(1) = fac * (bdot(2, 2) * bdot(3, 3) - bdot(2, 3) * bdot(2, 3))
  aadot(2) = fac * (bdot(3, 3) * bdot(1, 1) - bdot(3, 1) * bdot(3, 1))
  aadot(3) = fac * (bdot(1, 1) * bdot(2, 2) - bdot(1, 2) * bdot(1, 2))
  aadot(4) = fac * (bdot(1, 3) * bdot(3, 2) - bdot(1, 2) * bdot(3, 3))
  aadot(5) = fac * (bdot(2, 1) * bdot(1, 3) - bdot(2, 3) * bdot(1, 1))
  aadot(6) = fac * (bdot(1, 2) * bdot(2, 3) - bdot(1, 3) * bdot(2, 2))
  !     -------------------------------------------------------------
  !                 local part of stress contribution
  !
  !

  fac = dfour * pi2 / vcell  

  do i = 1, ng                                     ! loop over gvectors
     if (ek(i) > dzero) then  
        !
        !           vlocal * den
        !
        if (isp) then  
           dentot = den(i, 1) + den(i, 2)  
        else  
           dentot = den(i, 1)  
        end if
        denr = real(dentot, dp)  
        deni = aimag(dentot)  
        ! ionic potential, same for both spins
        vd = vion(i, 1) - vin(i, 1)  
        vld = real(vd, dp) * denr + aimag(vd) * deni  
        !
        !           derivative of v local * den
        !
        dvld = real(dvql(i), dp) * denr + aimag(dvql(i)) * deni  
        !
        !           vhartree * den / 2
        !
        vhd = fac * (denr * denr + deni * deni) / (dtwo * ek(i))  
        !
        !           vxc * den and vxc * d(den)/d(g**2)  (for g .ne. 0)
        !           (vxc (g=0) and exc done below)
        !
        !
        !           in the spin polarized case, take  denc/2.0, and ddc/2.0,
        !           but sum over both spins
        !
        vxd = dzero
        dvxd = dzero  
        do is = 1, ispn  
           vxc = vout(i, is) - fac * dentot / ek(i)  
           densum = den(i, is) + denc(i) / real(ispn, dp)  
           vxd = vxd + real(vxc, dp) * real(densum, dp) + &
                aimag(vxc) * aimag(densum)
           dvxd = dvxd + (real(vxc, dp) * real(ddc(i), dp) + &
                aimag(vxc) * aimag(ddc(i))) / real(ispn, dp)
        end do
        !
        !           diagonal (exp1) and nondiag (exp2) contributions
        !
        exp1 =  (vld + vhd + vxd)  
        exp2 = dtwo * (dvld - vhd / ek(i) + dvxd)  
        !
        !           add diagonal terms to stress tensor
        !
        do j = 1, 6  
           ssum(j) = ssum(j) + exp1 * aadot(j)  
        end do
        !
        !           and add nondiagonal terms
        !
        ssum(1) = ssum(1) + exp2 * real(gsvec(1, i) * gsvec(1, i), dp)
        ssum(2) = ssum(2) + exp2 * real(gsvec(2, i) * gsvec(2, i), dp)
        ssum(3) = ssum(3) + exp2 * real(gsvec(3, i) * gsvec(3, i), dp)
        ssum(4) = ssum(4) + exp2 * real(gsvec(1, i) * gsvec(2, i), dp)
        ssum(5) = ssum(5) + exp2 * real(gsvec(2, i) * gsvec(3, i), dp)
        ssum(6) = ssum(6) + exp2 * real(gsvec(3, i) * gsvec(1, i), dp)
     else                                            ! ha! here we have G=0
        !
        !           vxc (g=0), exc, and ealpha contribution
        !
        if (isp) then  
           exp1 = real(vout(i, 1), dp) * real(den(i, 1) + &
                denc(i) * dhalf, dp) + real(vout(i, 2), dp) * &
                real(den(i, 2) + denc(i) * dhalf) - exc + ealpha
        else  
           exp1 = real(vout(i, 1), dp) * real(den(i, 1) + denc(i))  - &
                exc + ealpha
        end if
        do j = 1, 6  
           ssum(j) = ssum(j) + exp1* aadot(j)  
        end do
     end if
  end do                                               ! loop over gvectors

  return  

end subroutine forstressloc
!
!
!
!
!    _________________________________________________________________
!
subroutine stresspw91(ffts, crys, gs, den, denc, ssum)

  use all_to_all_module  
  include 'use.h'  
  implicit none  
  include 'interface.h'  
  include 'all_to_all.h'  
  !
  !     INPUT:
  !     ------
  !
  type(fft_struc), intent(in) :: ffts              ! for the fourier transform
  type(crystal), intent(in) :: crys  ! for the bvec array to take the gradient
  type(parallel_gspace), intent(in) :: &
       gs                            ! the gspace for the gradient computation
  !
  !     INPUT/OUTPUT:
  !     ------------
  !
  !     the charge density in fourier space on input, the exchange
  !     -correlation potential (Rydbergs) in fourier space on output
  !
  complex(dp), intent(inout) :: den(gs%length, crys%nspin), denc(gs%length)
  real(dp), intent(inout) :: ssum(6)
  !
  !     --------------------- local variables -------------------------
  real(dp) :: ggss(6), exp3, exp4, alfr, &
       alfden, alflog, dalfdrs, zeta, zeta3, zeta4, fzeta, fzetap, &
       fzetazeta4, phi, phi2, phi3, invphi, invphi3, &
       phip, pider1, pider2, &
       pider3, abder1, abder2, abder3, vcell, invcell, spscale, &
       invspscale, roe, grx, gry, grz, roeth, rs, gkf, sd, s, s2, s3, s4, &
       ysr, ys, ysl, sysl, fxexp, fxnum, fxden, rs12, rs32, rs2, rs3, &
       invrs12, td, t, t2, t3, t4, aexp, abig, abig2, h0num, h0den, &
       h0arg, h0, ccrnum, ccrden, ccr, h1, dys, bys, dfs, dfs2, ecr_dif, &
       dh0dt2, dh1dt2, dh1ds2, normfac
  !
  !     _________________________________________________________________
  !
  integer :: k, ir, nspin, rsize, is, isign, norm, i, j  
  real(dp), dimension(1:3) :: aroe, agr, ecr  
  real(dp), dimension(1:2) :: ecden, eclog, dfxd, dfxdg, dexcd
  real(dp), dimension(3, 3) :: avec  
  real(dp), allocatable :: rhoe(:,:), &            ! charge density in rspace
       gradr(:,:,:)                    ! gradient of charge density in rspace
  complex(dp), allocatable :: gswork(:), &             ! work array in gspace
       gswork1(:), &                                   ! work array in gspace
       rswork(:)                                       ! work array in rspace
  !
  real(dp), parameter :: ax = -0.738558766382022405884230032680836d0, &
       gam = 0.5198420997897463295344212145565d0, invgam = done / gam, &
       invfzz = 9.0d0 * gam / 8.0d0, eta = 1.0d-12, delt = 1.0d-12, &
       bb1 = 0.19645d0, bb2 = 0.27430d0, bb3 = -0.15084d0, &
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
  logical :: isp                                          ! spin polarization
  !
  !     ------------------------------------------------------------------
  !
  real(dp), parameter :: x13 = dthird, x16m = dmone / dsix, &
       x23 = dtwo / dthree, x43 = dfour / dthree, x76 = 7.0d0 / dsix  
  !     _________________________________________________________________
  !     derived parameters from pi
  !
  real(dp), parameter :: pisq = pi * pi
  
  pider1 = (0.75d0 / pi)**x13  
  pider2 = (dthree * pisq)**x13  
  pider3 = (dthree * pisq / 16.0d0)**x13  
  !     _________________________________________________________________
  !     derived parameters from alfa and beta
  !
  abder1 = beta * beta / (dtwo * alfa)  
  abder2 = done / abder1  
  abder3 = dtwo * alfa / beta  
  !
  norm = gs%fftsize(1) * gs%fftsize(2) * gs%fftsize(3)  
  normfac = crys%vcell / real(norm, dp)  
  vcell = crys%vcell  
  invcell = done / crys%vcell  
  avec(:,:) = crys%avec(:,:)  
  !
  nspin = crys%nspin  
  spscale = real(nspin, dp)  
  invspscale = done / spscale  
  isp = .false.  
  if (nspin == 2) isp = .true.  
  !
  rsize = gs%r_size  
  phi = done
  invphi = done 
  phi2 = done
  phi3 = done
  invphi3 = done
  !
  allocate(rswork(rsize))  
  allocate(gswork(gs%length))  
  allocate(gswork1(gs%length))  
  !
  allocate(rhoe(rsize, nspin))  
  allocate(gradr(rsize, 3, 2 * nspin - 1))  
  !
  do i = 1, 6  
     ggss(i) = dzero  
  end do
  do is = 1, nspin  
     gswork(:) = (den(:, is) + denc(:) * invspscale) * invcell  
     call fourier_transform(-1, ffts, gs, gswork(1), rswork(1), 1)
     rhoe(:, is) = real(rswork(:), dp)  
     do k = 1, 3  
        call derivative_gspace(k, gswork(1), gs, crys, gswork1(1))  
        gswork1(:) = gswork1(:) * zi
        call fourier_transform(-1, ffts, gs, gswork1(1), rswork(1), 1)
        gradr(:, k, is) = real(rswork(:), dp)  
     end do
  end do
  if (isp) gradr(:,:, 3) = gradr(:,:, 1) + gradr(:,:, 2)  
  !
  !     _________________________________________________________________
  !     main loop
  !
  do 100 ir = 1, rsize  
     do is = 1, nspin  
        roe = rhoe(ir, is) * spscale  
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
           exp3 = dzero
        else  
           rs = pider1 / roeth  
           gkf = pider2 * roeth  
           sd = done / (dtwo * gkf * aroe(is))  
           s = agr(is) * sd  
           s2 = s * s  
           s3 = s * s2  
           s4 = s2 * s2  
           ysr = sqrt(done + bb5 * bb5 * s2)  
           ys = bb5 * s + ysr  
           ysl = log(ys) * bb1  
           sysl = s * ysl  
           fxexp = exp(-1.0d2 * s2)  
           fxnum = done + sysl + (bb2 + bb3 * fxexp) * s2  
           fxden = done / (done + sysl + bb4 * s4)  
           dys = bb5 * (done + bb5 * s / ysr) / ys  
           dfs = -fxnum * (ysl + bb1 * s * dys + dfour * bb4 * s3) * &
                fxden * fxden + (ysl + bb1 * s * dys + dtwo * s * &
                (bb2 + bb3 * fxexp) - 2.0d2 * s3 * bb3 * fxexp) * fxden
           isign = sign(done, s - delt)  
           bys = dhalf * real(1+isign, dp) / (s + real(1-isign, dp) * delt)  
           dfs2 = dhalf * dfs * bys  
           exp3 = spscale * ax * dfs2 * roeth * normfac / (gkf * gkf * &
                aroe(is) * dfour * pisq)
           do i = 1, 3  
              do j = 1, 3  
                 ggss(1) = ggss(1) + exp3 * avec(i, 1) * gradr(ir, i, is) * &
                      gradr(ir, j, is) * avec(j, 1)
                 ggss(2) = ggss(2) + exp3 * avec(i, 2) * gradr(ir, i, is) * &
                      gradr(ir, j, is) * avec(j, 2)
                 ggss(3) = ggss(3) + exp3 * avec(i, 3) * gradr(ir, i, is) * &
                      gradr(ir, j, is) * avec(j, 3)
                 ggss(4) = ggss(4) + exp3 * avec(i, 1) * gradr(ir, i, is) * &
                      gradr(ir, j, is) * avec(j, 2)
                 ggss(5) = ggss(5) + exp3 * avec(i, 2) * gradr(ir, i, is) * &
                      gradr(ir, j, is) * avec(j, 3)
                 ggss(6) = ggss(6) + exp3 * avec(i, 3) * gradr(ir, i, is) * &
                      gradr(ir, j, is) * avec(j, 1)
              end do
           end do
        end if
     end do
     !     _________________________________________________________________
     !     back to total charge density
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
     !     _________________________________________________________________
     !     correlation ecr=ec(rho,zeta)
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
        if (abs(zeta) > done) zeta = sign(done, zeta)
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
        !
        ecr_dif = ecr(2) - ecr(1)  
        ecr(3) = ecr(1) + alfr * fzeta * invfzz * (done - zeta4) + &
             ecr_dif * fzetazeta4
     end if
     !     _________________________________________________________________
     !     correlation h0(t,ecr,zeta)
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
     !     _________________________________________________________________
     !     correlation h1(t,s,aroe,zeta)
     !
     ccrnum = c2 + c3 * rs + c4 * rs2  
     ccrden = done / (done + c5 * rs + c6 * rs2 + c7 * rs3)  
     ccr = c1 + ccrnum * ccrden  
     fxexp = exp(-1.0d2 * s2 * phi2)  
     h1 = phi3 * cc0 * (ccr - cc1) * t2 * fxexp  
     dh0dt2 = phi3 * abder1 / h0arg * abder3 * h0den * (done + &
          dtwo * abig * t2 - h0num * h0den * (done * abig + dtwo * &
          abig2 * t2))
     dh1dt2 = phi3 * cc0 * (ccr - cc1) * fxexp  
     dh1ds2 = -1.0d2 * phi2 * cc0 * (ccr - cc1) * phi3 * t2 * fxexp
     is = 2 * nspin - 1  
     exp4 = (dh0dt2 + dh1dt2) / (1.6d1 * phi2 * gkf * aroe(is) * pi) * &
          normfac + dh1ds2 * normfac / (gkf * gkf * aroe(is) * dfour * pisq)
     do i = 1, 3
        do j = 1, 3  
           ggss(1) = ggss(1) + exp4 * avec(i, 1) * gradr(ir, i, is) * &
                gradr(ir, j, is) * avec(j, 1)
           ggss(2) = ggss(2) + exp4 * avec(i, 2) * gradr(ir, i, is) * &
                gradr(ir, j, is) * avec(j, 2)
           ggss(3) = ggss(3) + exp4 * avec(i, 3) * gradr(ir, i, is) * &
                gradr(ir, j, is) * avec(j, 3)
           ggss(4) = ggss(4) + exp4 * avec(i, 1) * gradr(ir, i, is) * &
                gradr(ir, j, is) * avec(j, 2)
           ggss(5) = ggss(5) + exp4 * avec(i, 2) * gradr(ir, i, is) * &
                gradr(ir, j, is) * avec(j, 3)
           ggss(6) = ggss(6) + exp4 * avec(i, 3) * gradr(ir, i, is) * &
                gradr(ir, j, is) * avec(j, 1)
        end do
     end do
100 end do
  call all_sum_all(ggss(1))  
  call all_sum_all(ggss(2))  
  call all_sum_all(ggss(3))  
  call all_sum_all(ggss(4))  
  call all_sum_all(ggss(5))  
  call all_sum_all(ggss(6))  
  ssum(:) = ssum(:) + ggss(:)  
  !     _________________________________________________________________
  deallocate(rswork)  
  deallocate(gswork)  
  deallocate(gswork1)  
  !
  deallocate(rhoe)  
  deallocate(gradr)  
  !     _________________________________________________________________
  return

end subroutine stresspw91
!     _________________________________________________________________
!
!
!
!
!
subroutine stresspbe(ffts, crys, gs, den, denc, ssum)

  use all_to_all_module  
  include 'use.h'  
  implicit none  
  include 'interface.h'  
  include 'all_to_all.h'  
  !
  !     INPUT:
  !     ------
  !
  type(fft_struc), intent(in) :: ffts              ! for the fourier transform
  type(crystal), intent(in) :: crys  ! for the bvec array to take the gradient
  type(parallel_gspace), intent(in) :: &
       gs                            ! the gspace for the gradient computation
  !
  !     INPUT/OUTPUT:
  !     ------------
  !
  !     stress and the charge density in fourier space on input,
  !     gga contribution of stress is added to ssum.
  !
  complex(dp), intent(inout) :: den(gs%length, crys%nspin), denc(gs%length)
  real(dp), intent(inout) :: ssum(6)
  !
  !     --------------------- local variables -------------------------
  real(dp) :: ggss(6), &
       uk2, alfr, alfden, alflog, dalfdrs, ecr_dif, zeta, &
       zeta3, zeta4, fzeta, fzetap, fzetazeta4, phi, phi2, phi3, &
       invphi, invphi3, pider1, &
       pider2, pider3, abder1, abder2, abder3, vcell, invcell, spscale, &
       invspscale, roe, grx, gry, grz, roeth, rs, gkf, sd, s, s2, fxden, &
       fxden2, fx, rs12, rs32, rs2, rs3, dfs, td, t, t2, t3, t4, aexp, &
       abig, abig2, hnum, hden, harg, h, dhdt2, normfac, exp3, exp4
  !
  !     _________________________________________________________________
  !
  integer :: k, ir, nspin, rsize, is, isign, norm, i, j  
  real(dp), dimension(1:3) :: aroe, agr, ecr  
  real(dp), dimension(1:2) :: ecden, eclog  
  real(dp), dimension(3, 3) :: avec  
  real(dp), allocatable :: rhoe(:,:), &             ! charge density in rspace
       gradr(:,:,:)                     ! gradient of charge density in rspace
  complex(dp), allocatable :: gswork(:), &              ! work array in gspace
       gswork1(:), &                                    ! work array in gspace
       rswork(:)                                        ! work array in rspace
  !
  real(dp), parameter :: ax = -0.738558766382022405884230032680836d0, &
       um = 0.2195149727645171d0, uk = 0.804d0, & 
       gam = 0.5198420997897463295344212145565d0, invgam = done / gam, &
       invfzz = 9.0d0 * gam / 8.0d0, &
       gamma = 0.03109069086965489503494086371273d0, & 
       bet = 0.06672455060314922d0, &
       eta = 1.0d-12, delt = 1.0d-12, & 
       a = 0.0621814d0, alfa1 = 0.2137d0, bt1 = 7.5957d0, bt2 = 3.5876d0, &
       bt3 = 1.6382d0, bt4 = 0.49294d0, a_p = 0.0310907d0, &
       alfa1_p = 0.20548d0, bt1_p = 14.1189d0, bt2_p = 6.1977d0, &
       bt3_p = 3.3662d0, bt4_p = 0.62517d0, a_alf = 0.033774d0, &
       alfa1_alf = 0.11125d0, bt1_alf = 10.357d0, bt2_alf = 3.6231d0, &
       bt3_alf = 0.88026d0, bt4_alf = 0.49671d0
  !
  logical :: isp                                           ! spin polarization
  !
  !     ------------------------------------------------------------------
  !
  real(dp), parameter :: x13 = dthird, &  
       x16m = -done / dsix, &  
       x23 = dtwo / dthree, & 
       x43 = dfour / dthree, &  
       x76 = 7.0d0 / 6.0d0  
  !     _________________________________________________________________
  !     derived parameters from pi
  !
  real(dp), parameter :: pisq = pi * pi

  pider1 = (0.75d0 / pi)**x13  
  pider2 = (dthree * pisq)**x13  
  pider3 = (dthree * pisq / 1.6d1)**x13  
  !     _________________________________________________________________
  !     derived parameters from beta, gamma, and kappa
  !
  abder1 = gamma  
  abder2 = done / abder1  
  abder3 = bet * abder2  
  uk2 = uk * uk  
  !
  norm = gs%fftsize(1) * gs%fftsize(2) * gs%fftsize(3)  
  normfac = crys%vcell / real(norm, dp)  
  vcell = crys%vcell  
  invcell = done / crys%vcell  
  avec(:,:) = crys%avec(:,:)  
  !
  nspin = crys%nspin  
  spscale = real(nspin, dp)  
  invspscale = done / spscale  
  isp = .false.  
  if (nspin == 2) isp = .true.  
  !
  rsize = gs%r_size  
  phi = done
  invphi = done
  phi2 = done
  phi3 = done
  invphi3 = done
  !
  allocate(rswork(rsize))  
  allocate(gswork(gs%length))  
  allocate(gswork1(gs%length))  
  !
  allocate(rhoe(rsize, nspin))  
  allocate(gradr(rsize, 3, 2 * nspin - 1))  
  do i = 1, 6  
     ggss (i) = dzero
  end do
  !
  do is = 1, nspin  
     gswork(:) = (den(:, is) + denc(:) * invspscale) * invcell  
     call fourier_transform(-1, ffts, gs, gswork(1), rswork(1), 1)
     rhoe(:, is) = real(rswork(:), dp)  
     do k = 1, 3  
        call derivative_gspace(k, gswork(1), gs, crys, gswork1(1))  
        gswork1(:) = gswork1(:) * zi 
        call fourier_transform(-1, ffts, gs, gswork1(1), rswork(1), 1)
        gradr(:, k, is) = real(rswork(:), dp)  
     end do
  end do
  if (isp) gradr(:,:, 3) = gradr(:,:, 1) + gradr(:,:, 2)  
  !
  !     _________________________________________________________________
  !     main loop
  !
  do 100 ir = 1, rsize  
     do is = 1, nspin  
        roe = rhoe(ir, is) * spscale                ! spin-scaled for exchange
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
           exp3 = dzero
        else  
           rs = pider1 / roeth  
           gkf = pider2 * roeth  
           sd = done / (dtwo * gkf * aroe(is))  
           s = agr(is) * sd  
           s2 = s * s  
           fxden = done / (uk + um * s2)  
           fxden2 = fxden * fxden  
           fx = done + uk - uk2 * fxden  
           ! au to ry     exp3=spscale*ax*um*uk2*fxden2*roeth
           !                   /(2.d0*gkf*gkf*aroe(is))*normfac*2.d0/(4.d0*pi2)
           exp3 = spscale * ax * um * uk2 * fxden2 * roeth * normfac / &
                (gkf * gkf * aroe(is) * dfour * pisq)
           do i = 1, 3  
              do j = 1, 3  
                 ggss(1) = ggss(1) + exp3 * avec(i, 1) * gradr(ir, i, is) * &
                      gradr(ir, j, is) * avec(j, 1)
                 ggss(2) = ggss(2) + exp3 * avec(i, 2) * gradr(ir, i, is) * &
                      gradr(ir, j, is) * avec(j, 2)
                 ggss(3) = ggss(3) + exp3 * avec(i, 3) * gradr(ir, i, is) * &
                      gradr(ir, j, is) * avec(j, 3)
                 ggss(4) = ggss(4) + exp3 * avec(i, 1) * gradr(ir, i, is) * &
                      gradr(ir, j, is) * avec(j, 2)
                 ggss(5) = ggss(5) + exp3 * avec(i, 2) * gradr(ir, i, is) * &
                      gradr(ir, j, is) * avec(j, 3)
                 ggss(6) = ggss(6) + exp3 * avec(i, 3) * gradr(ir, i, is) * &
                      gradr(ir, j, is) * avec(j, 1)
              end do
           end do
        end if
     end do
     !     _________________________________________________________________
     !     back to total charge density
     !        Be careful! Total charge density is not spin-scaled!
     !
     if (isp) then  
        roe = rhoe(ir, 1) + rhoe(ir, 2)  
        if (roe == dzero) goto 100  
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
     !     _________________________________________________________________
     !     correlation ecr=ec(rho,zeta)
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
        alfden = a_alf * (bt1_alf * rs12 + bt2_alf * rs + bt3_alf * &
             rs32 + bt4_alf * rs2)
        alflog = log(done + (done / alfden))  
        alfr = a_alf * (done + alfa1_alf * rs) * alflog  
        zeta = (rhoe(ir, 1) - rhoe(ir, 2)) / roe  
        if (abs(zeta) > done) zeta = sign(done, zeta)  
        zeta3 = zeta * zeta * zeta  
        zeta4 = zeta3 * zeta  
        fzeta = ((done + zeta)**x43 + (done - zeta)**x43 - dtwo) * invgam
        fzetap = x43 * ((done + zeta)**x13 - (done - zeta)**x13) * invgam
        fzetazeta4 = fzeta * zeta4  
        phi = ((done + zeta)**x23 + (done - zeta)**x23) * dhalf  
        invphi = done / phi  
        phi2 = phi * phi  
        phi3 = phi * phi * phi  
        invphi3 = done / phi3  
        ecr_dif = ecr(2) - ecr(1)  
        ecr (3) = ecr(1) + alfr * fzeta * invfzz * (done - zeta4) + &
             ecr_dif * fzetazeta4
     end if
     !     _________________________________________________________________
     !     correlation h(t,ecr,zeta)
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
     dhdt2 = phi3 * abder1 / harg * abder3 * hden * (done + dtwo * &
          abig * t2 - hnum * hden * (done * abig + dtwo * abig2 * t2))
     is = 2 * nspin - 1  
     !        exp4=dhdt2/(2.d0*phi2*4.0d0*gkf/pi*aroe(is))
     !             2.d0*normfac/(4.d0*pi2)
     exp4 = dhdt2 / (1.6d1 * phi2 * gkf * aroe(is) * pi) * normfac  
     do i = 1, 3  
        do j = 1, 3  
           ggss(1) = ggss(1) + exp4 * avec(i, 1) * gradr(ir, i, is) * &
                gradr(ir, j, is) * avec(j, 1)
           ggss(2) = ggss(2) + exp4 * avec(i, 2) * gradr(ir, i, is) * &
                gradr(ir, j, is) * avec(j, 2)
           ggss(3) = ggss(3) + exp4 * avec(i, 3) * gradr(ir, i, is) * &
                gradr(ir, j, is) * avec(j, 3)
           ggss(4) = ggss(4) + exp4 * avec(i, 1) * gradr(ir, i, is) * &
                gradr(ir, j, is) * avec(j, 2)
           ggss(5) = ggss(5) + exp4 * avec(i, 2) * gradr(ir, i, is) * &
                gradr(ir, j, is) * avec(j, 3)
           ggss(6) = ggss(6) + exp4 * avec(i, 3) * gradr(ir, i, is) * &
                gradr(ir, j, is) * avec(j, 1)
        end do
     end do
100 end do

  call all_sum_all(ggss(1))  
  call all_sum_all(ggss(2))  
  call all_sum_all(ggss(3))  
  call all_sum_all(ggss(4))  
  call all_sum_all(ggss(5))  
  call all_sum_all(ggss(6))  
  ssum(:) = ssum(:) + ggss(:)  
  !     _________________________________________________________________
  deallocate(rswork)  
  deallocate(gswork)  
  deallocate(gswork1)  
  !
  deallocate(rhoe)  
  deallocate(gradr)  
  !     _________________________________________________________________
  return

end subroutine stresspbe
