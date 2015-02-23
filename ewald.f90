!*
subroutine ewald(ioutput, crys, gs, eewald, ew)  
  !
  use all_to_all_module  
  include 'use.h'  
  implicit none                    ! implicit? Just say no!
  include 'interface.h'
  include 'all_to_all.h'
  !
  !     INPUT:
  !     -----
  !
  type(crystal), intent(in) :: crys  
  type(parallel_gspace), intent(in) :: gs  
  integer, intent(in) :: ioutput                                  ! print flag
  !
  !     OUTPUT:
  !     ------
  !
  real(dp), intent(out) :: eewald                         ! ewald energy (Ryd)
  type(force), intent(out) :: ew             ! structure containing the forces
  !
  !     DESCRIPTION:
  !     -----------
  !
  !     computes the ewald contribution to the
  !     total energy, forces and stresses.
  !     adapted from Sverre Froyen plane wave program
  !     written december 16 1987. JLM
  !
  !     parallel version:          1996 Bernd Pfrommer, UC Berkeley
  !
  !     better scaling of the gspace loop 1997 Bernd Pfrommer
  !
  !     ---------------- local variables ---------------------------------
  !
  real(dp), parameter :: small = 1.0d-12
  real(dp) :: tol, t0, &
       eps, &                                    ! corresponds to Kittel's eta
       seps, sepi, &
       rmax, &
       zz2, &                                        ! sum_{ij}    zv(i)*zv(j)
       adot (3, 3), &                                   ! metric in real space
       factor, rmod, arg, &
       esumg, esumr, esub, enorm, &
       ssumg(6), &                                      ! stress sum in gspace
       ssumr(6), &                                  ! stress sum in real space
       ssum0(6), &                            ! stress sum for the case of a=b
       tau(3), &                                           ! difference vector
       esum0, exp1, exp2, exp3, &
       emax, &                                                 ! energy cutoff
       ssub(6), fsub(3), rp(3), gdt
  real(dp), external :: gimmetime  
  !
  !     variables for the gspace loop:
  !
  integer :: igs, iord, irod, igv(4), igv3, ffth(4), fftn(4)  
  real(dp) :: gv(3)
  !  
  !     WORK ARRAYS:
  !     -----------
  !
  real(dp), allocatable :: zz(:), &          ! size zz(crys%ntype*crys%mxdatm)
       rc(:,:), &                          ! size rc(3,crys%ntype*crys%mxdatm)
       fsumr(:,:), fsumg(:,:), &         ! size fsum(3,crys%ntype*crys%mxdatm)
       expg(:)                                          ! size expg(gs%length)
  complex(dp) :: sfac  
  complex(dp), allocatable :: cfsumg(:,:)  
  real(dp) :: cutfac  
  integer :: imx, jmx, kmx, &
       natot, &                 ! total number of atoms, no matter which kind
       i, j, k, n, im, l, m, ii, jj, kk, &
       ir(3), im2, jm2, km2
  real(dp), external :: myderfc  

  !     ------------------------------------------------------------------

  t0 = gimmetime()  

  !     alloc mem
  allocate(expg(gs%length))  
  allocate(zz(crys%ntype * crys%mxdatm))  
  allocate(rc(3, crys%ntype * crys%mxdatm))  
  allocate(fsumr(3, crys%ntype * crys%mxdatm))  
  allocate(fsumg(3, crys%ntype * crys%mxdatm))  
  allocate(ew%force(3, crys%mxdatm, crys%ntype))  
  !
  !     compute metric in real space
  !
  factor = (crys%vcell / (pi2 * pi2)) * (crys%vcell / (pi2 * pi2))
  adot(1, 1) = factor * (crys%bdot(2, 2) * crys%bdot(3, 3) - &
       crys%bdot(2, 3) * crys%bdot(2, 3))
  adot(2, 2) = factor * (crys%bdot(3, 3) * crys%bdot(1, 1) - &
       crys%bdot(3, 1) * crys%bdot(3, 1))
  adot(3, 3) = factor * (crys%bdot(1, 1) * crys%bdot(2, 2) - &
       crys%bdot(1, 2) * crys%bdot(1, 2))
  adot(1, 2) = factor * (crys%bdot(1, 3) * crys%bdot(3, 2) - &
       crys%bdot(1, 2) * crys%bdot(3, 3))
  adot(1, 3) = factor * (crys%bdot(1, 2) * crys%bdot(2, 3) - &
       crys%bdot(1, 3) * crys%bdot(2, 2))
  adot(2, 3) = factor * (crys%bdot(2, 1) * crys%bdot(1, 3) - &
       crys%bdot(2, 3) * crys%bdot(1, 1))
  adot(2, 1) = adot(1, 2)  
  adot(3, 1) = adot(1, 3)  
  adot(3, 2) = adot(2, 3)  
  !
  !     store data for atoms into new arrays
  !
  natot = 0  
  do i = 1, crys%ntype  
     do j = 1, crys%natom(i)  
        natot = natot + 1  
        zz(natot) = crys%zv(i)  
        rc(1, natot) = crys%rat(1, j, i) / pi2  
        rc(2, natot) = crys%rat(2, j, i) / pi2  
        rc(3, natot) = crys%rat(3, j, i) / pi2  
        rc(1, natot) = rc(1, natot) - real(int(rc(1, natot)), dp)  
        if (rc(1, natot) < dzero) rc(1, natot) = rc(1, natot) + done  
        rc(2, natot) = rc(2, natot) - real(int(rc(2, natot)), dp)  
        if (rc(2, natot) < dzero) rc(2, natot) = rc(2, natot) + done  
        rc(3, natot) = rc(3, natot) - real(int(rc(3, natot)), dp)  
        if (rc(3, natot) < dzero) rc(3, natot) = rc(3, natot) + done 
     end do
  end do
  !
  zz2 = dzero  
  do i = 1, natot  
     do j = 1, natot  
        zz2 = zz2 + zz(i) * zz(j)  
     end do
  end do
  !
  !     compute various parameters
  !
  tol = -log(small)  

  eps = dtwo * gs%gmax / (dfour * tol)  
  seps = sqrt(eps)  
  sepi = dtwo * seps / sqrt(pi)  
  ! increase precision with number of atoms:
  rmax = sqrt((tol + log(real(natot, dp))) / eps)  
  !
  !     Changed the following to have +2 instead of +1
  !     to avoid problems with large gspace cutoffs. The max()
  !     expression provides a safeguard against the case where the
  !     atoms are way outside of the unit cell. Then, the rspace sum
  !     must be extended appropriately.
  !
  imx = int(rmax * sqrt(crys%bdot(1, 1)) / pi2) + 2  
  jmx = int(rmax * sqrt(crys%bdot(2, 2)) / pi2) + 2  
  kmx = int(rmax * sqrt(crys%bdot(3, 3)) / pi2) + 2  
  !
  !      initialize sums : esum(g,r)   -   energy
  !                        fsum(g,r)   -   force
  !                        ssum(g,r)   -   stress
  !
  !      note that the sums are in units of basis vectors
  !      (g) in units of bi's, and (r) in units of ai's
  !
  esumg = dzero  
  esumr = dzero  
  do i = 1, 3  
     ssumg(i + 3) = dzero  
     ssumr(i + 3) = dzero  
     ssumg(i) = dzero  
     ssumr(i) = dzero  
     do j = 1, natot  
        fsumg(i, j) = dzero  
        fsumr(i, j) = dzero  
     end do
  end do
  !     ------------------------------------------------------------------
  !                          DO SUM IN GSPACE
  !
  !     (this is done in parallel)
  !
  !     compute forces
  !
  allocate(cfsumg(3, natot))  
  cfsumg = zzero 

  esumg = dzero
  fftn(1:3) = gs%fftsize(1:3)  
  fftn(4) = gs%fftsize(3)  

  ffth(:) = fftn(:) / 2  
  igs = 0  
  do iord = 1, gs%lorder                             ! loop through x/y gspace
     irod = gs%order(1, iord)  
     igv(1) = irod / gs%fftsize(2)  
     igv(2) = mod(irod, gs%fftsize(2))  
     igv(3) = gs%order(2, iord)  
     igv(4) = gs%order(3, iord)  
     igv(:) = mod(igv(:) + ffth(:), fftn(:)) - ffth(:)  
     gv (1:2) = real(igv(1:2), dp)  
     do igv3 = igv(3), igv(4)                               ! loop over z axis
        gv(3) = real(igv3, dp)  
        igs = igs + 1  
        if (gs%ekin(igs) > dzero) then  
           cutfac = exp(-gs%ekin(igs) / (dfour * eps)) / gs%ekin(igs)  
           !
           !              first precompute the structure factor
           !
           sfac = zzero
           do i = 1, natot  
              gdt = pi2 * (gv(1) * rc(1, i) + gv(2) * rc(2, i) + &
                   gv(3) * rc(3, i))
              sfac = sfac + zz(i) * cmplx(cos(gdt), -sin(gdt), dp)  
           end do
           !
           !              contribution to the energy ...
           !
           esumg = esumg + real(sfac * conjg(sfac)) * cutfac  
           !
           !              .... to the stress ...
           !
           exp3 = -(done / (dtwo * eps) + dtwo / gs%ekin(igs))  
           exp3 = exp3 * real(sfac * conjg(sfac), dp) * cutfac  
           ssumg(1) = ssumg(1) + exp3 * gv(1) * gv(1)  
           ssumg(2) = ssumg(2) + exp3 * gv(2) * gv(2)  
           ssumg(3) = ssumg(3) + exp3 * gv(3) * gv(3)  
           ssumg(4) = ssumg(4) + exp3 * gv(1) * gv(2)  
           ssumg(5) = ssumg(5) + exp3 * gv(2) * gv(3)  
           ssumg(6) = ssumg(6) + exp3 * gv(3) * gv(1)  
           !
           !              and to the forces.
           !
           sfac = sfac * cutfac
           do i = 1, natot  
              gdt = pi2 * (gv(1) * rc(1, i) + gv(2) * rc(2, i) + &
                   gv(3) * rc(3, i))
              cfsumg(:, i) = cfsumg(:, i) + sfac * &
                   cmplx(cos(gdt), sin(gdt), dp) * gv(:)
           end do
        end if
     end do
  end do
  !
  !     put in the prefactors etc
  !
  do j = 1, natot  
     cfsumg(1:3, j) = 8.0d0 * pi * zz(j) * cfsumg(1:3, j) / crys%vcell
  end do
  esumg = pi4 * esumg / crys%vcell  
  ssumg(:) = pi4 * ssumg(:) / crys%vcell  
  !
  !     sum across the processors
  !
  call all_sum_all(cfsumg, 3 * natot)  
  call all_sum_all(esumg)  
  call all_sum_all(ssumg, 6)  

  fsumg(1:3, 1:natot) = aimag(cfsumg(1:3, 1:natot))  

  deallocate(cfsumg)  
  !
  !     ------------------------------------------------------------
  !
  !                          START SUM IN R SPACE
  !     (no parallization done. all processors compute the same thing)
  im2 = 2 * imx + 1  
  jm2 = 2 * jmx + 1  
  km2 = 2 * kmx + 1  
  esum0 = dzero  

  ssum0 = dzero  
  do i = 1, im2  
     ir(1) = i - imx - 1  
     do j = 1, jm2  
        ir(2) = j - jmx - 1  
        do k = 1, km2  
           ir(3) = k - kmx - 1  
           rmod = dzero  
           do l = 1, 3  
              do m = 1, 3  
                 rmod = rmod + real(ir(l), dp) * adot(l, m) * real(ir(m), dp)  
              end do
           end do
           if (rmod /= dzero) then  
              rmod = sqrt(rmod)  
              arg = seps * rmod  
              if (arg < 2.5d1) then  
                 exp1 = myderfc(arg) / rmod  
                 exp2 = (exp1 + sepi * exp(-arg * arg)) / (rmod * rmod)  
                 esum0 = esum0 + exp1  
                 ssum0(1) = ssum0(1) + exp2 * real(ir(1) * ir(1), dp)  
                 ssum0(2) = ssum0(2) + exp2 * real(ir(2) * ir(2), dp)  
                 ssum0(3) = ssum0(3) + exp2 * real(ir(3) * ir(3), dp)  
                 ssum0(4) = ssum0(4) + exp2 * real(ir(1) * ir(2), dp)  
                 ssum0(5) = ssum0(5) + exp2 * real(ir(2) * ir(3), dp)  
                 ssum0(6) = ssum0(6) + exp2 * real(ir(3) * ir(1), dp)  
              end if
           end if
        end do
     end do
  end do

  esum0 = esum0 - pi / (eps * crys%vcell) - sepi  
  !
  !      start loop over atoms in cell
  !
  do i = 1, natot  
     !
     !        term with a=b
     !
     esumr = esumr + zz(i) * zz(i) * esum0  
     do j = 1, 6  
        ssumr(j) = ssumr(j) + zz(i) * zz(i) * ssum0(j)  
     end do
     im = i - 1  
     if (im /= 0) then  
        !
        !           terms with a#b
        !
        do j = 1, im                      ! loop over other atoms in unit cell
           esub = dzero  
           do k = 1, 3  
              fsub(k) = dzero  
              ssub(k + 3) = dzero  
              ssub(k) = dzero  
           end do
           do ii = 1, im2                ! ---------  loop over lattice points
              ir(1) = ii - imx - 1  
              do jj = 1, jm2  
                 ir(2) = jj - jmx - 1  
                 do kk = 1, km2  
                    ir(3) = kk - kmx - 1  
                    rp(1) = real(ir(1), dp) + rc(1, i) - rc(1, j)  
                    rp(2) = real(ir(2), dp) + rc(2, i) - rc(2, j)  
                    rp(3) = real(ir(3), dp) + rc(3, i) - rc(3, j)  
                    rmod = dzero  
                    do l = 1, 3  
                       do m = 1, 3  
                          rmod = rmod + rp(l) * adot(l, m) * rp(m)  
                       end do
                    end do
                    rmod = sqrt(rmod)  
                    arg = seps * rmod  
                    if (arg < 2.5d1) then  
                       exp1 = myderfc(arg) / rmod  
                       exp2 = (exp1 + sepi * exp(-arg * arg)) / (rmod * rmod)
                       esub = esub + exp1  
                       fsub(1) = fsub(1) + rp(1) * exp2  
                       fsub(2) = fsub(2) + rp(2) * exp2  
                       fsub(3) = fsub(3) + rp(3) * exp2  
                       ssub(1) = ssub(1) + rp(1) * exp2 * rp(1)  
                       ssub(2) = ssub(2) + rp(2) * exp2 * rp(2)  
                       ssub(3) = ssub(3) + rp(3) * exp2 * rp(3)  
                       ssub(4) = ssub(4) + rp(1) * exp2 * rp(2)  
                       ssub(5) = ssub(5) + rp(2) * exp2 * rp(3)  
                       ssub(6) = ssub(6) + rp(3) * exp2 * rp(1)  
                    end if
                 end do
              end do
           end do
           esub = esub - pi / (eps * crys%vcell)  
           esumr = esumr + dtwo * zz(i) * zz(j) * esub  
           do k = 1, 6  
              ssumr(k) = ssumr(k) + dtwo * zz(i) * zz(j) * ssub(k)  
           end do
           do k = 1, 3  
              fsumr(k, i) = fsumr(k, i) + dtwo * zz(i) * zz(j) * fsub(k)
              fsumr(k, j) = fsumr(k, j) - dtwo * zz(i) * zz(j) * fsub(k)
           end do
        end do                     ! end of loop over other atoms in unit cell
     end if
  end do
  !
  !      end r sum
  !
  eewald = esumg + esumr  
  !
  !      force
  !      note - returned force is in units of real space lattice vectors
  !             printed force is in cartesian coordinates
  !
  n = 0  
  do i = 1, crys%ntype  
     do j = 1, crys%natom(i)  
        n = n + 1  
        do k = 1, 3  
           ew%force(k, j, i) = fsumr(k, n)  
           do l = 1, 3  
              ew%force(k, j, i) = ew%force(k, j, i) + crys%bdot(k, l) * &
                   fsumg(l, n) / pi2
           end do
        end do
     end do
  end do
  !
  !     stress
  !      note - both returned and printed stress are
  !             in cartesian coordinates
  !
  do i = 1, 6  
     j = i  
     k = i  
     if (i > 3) j = i - 3  
     if (i > 3) k = j + 1  
     if (k > 3) k = 1  
     ew%stress(i) = crys%bvec(j, 1) * ssumg(1) * crys%bvec(k, 1) + &
          crys%avec(j, 1) * ssumr(1) * crys%avec(k, 1) + &
          crys%bvec(j, 2) * ssumg(2) * crys%bvec(k, 2) + &
          crys%avec(j, 2) * ssumr(2) * crys%avec(k, 2) + &
          crys%bvec(j, 3) * ssumg(3) * crys%bvec(k, 3) + &
          crys%avec(j, 3) * ssumr(3) * crys%avec(k, 3) + &
          crys%bvec(j, 1) * ssumg(4) * crys%bvec(k, 2) + &
          crys%avec(j, 1) * ssumr(4) * crys%avec(k, 2) + &
          crys%bvec(j, 2) * ssumg(5) * crys%bvec(k, 3) + &
          crys%avec(j, 2) * ssumr(5) * crys%avec(k, 3) + &
          crys%bvec(j, 3) * ssumg(6) * crys%bvec(k, 1) + &
          crys%avec(j, 3) * ssumr(6) * crys%avec(k, 1) + &
          crys%bvec(k, 1) * ssumg(4) * crys%bvec(j, 2) + &
          crys%avec(k, 1) * ssumr(4) * crys%avec(j, 2) + &
          crys%bvec(k, 2) * ssumg(5) * crys%bvec(j, 3) + &
          crys%avec(k, 2) * ssumr(5) * crys%avec(j, 3) + &
          crys%bvec(k, 3) * ssumg(6) * crys%bvec(j, 1) + &
          crys%avec(k, 3) * ssumr(6) * crys%avec(j, 1)
     if (i <= 3) then  
        ew%stress(i) = ew%stress(i) + esumg - zz2 * pi / (eps * crys%vcell)
     end if
  end do
  !
  !      printout
  !
  enorm = eewald * crys%vcell**dthird  

  write(9, 100) enorm  
  write(9, 101) eewald, enorm  
  n = 0  
  do i = 1, crys%ntype  
     do j = 1, crys%natom(i)  
        n = n + 1  
        do k = 1, 3  
           tau(k) = crys%avec(k, 1) * ew%force(1, j, i) + &
                crys%avec(k, 2) * ew%force(2, j, i) + &
                crys%avec(k, 3) * ew%force(3, j, i)
        end do
        write(9, 102) n, (crys%rat(k, j, i) / pi2, k = 1, 3), &
             (tau(k), k = 1, 3)
     end do
  end do

  write(9, 103) (ew%stress(i), i = 1, 6)  
  deallocate(expg)  
  deallocate(zz)  
  deallocate(rc)  
  deallocate(fsumr)  
  deallocate(fsumg)  

  if (iand (ioutput, 8) == 8) write(9, 940) gimmetime() - t0  
  call myflush(9)  

  return  

100 FORMAT (/,' NORMALIZED EWALD ENERGY =',F16.8)  
101 FORMAT (/,' EWALD ANALYSIS',/,1X,14('*'),/,50X,'1/3',/, &
       & ' ENERGY :',12X,'ENERGY (RY)',16X, &
       & '*V    (RY*A.U.)',/,16X,F16.10,14X,F16.9,// &
       & ' FORCES(CARTESIAN) :',3X,'N',11X,'COORD',19X,'FORCE (RY/A.U.)',/ &
       & 19X,'A1',4X,'A2',4X,'A3',11X,'-X-',7X,'-Y-',7X,'-Z-')
102 FORMAT (10X,I3,3X,3F6.3,5X,3F10.6)  
103 FORMAT (/,' STRESS(CARTESIAN) :',25X,'SIGMA * V (RY)',/, &
       & 14X,'-XX-',6X,'-YY-',6X,'-ZZ-',6X,'-XY-',6X,'-YZ-',6X,'-ZX-',/, &
       & 9X,6F10.6)

940 format(' TIME FOR EWALD CALCULATION:',f12.3)  

end subroutine ewald
