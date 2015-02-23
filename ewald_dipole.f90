!*
subroutine ewald_dipole(ioutput, crys, gs, ewald_tens)  
  !
  use all_to_all_module  
  include 'use.h'  
  implicit none                       ! implicit? Just say no!
  include 'interface.h'  
  include 'all_to_all.h'  
  !
  !     INPUT:
  !     -----
  !
  type(crystal), intent(in) :: crys         ! the crystal structure parameters
  type(parallel_gspace), intent(in) :: &
       gs                          ! the gspace structure and everything in it
  integer, intent(in) :: ioutput                                  ! print flag

  !
  !     OUTPUT:
  !     ------
  !
  real(dp), intent(out) :: ewald_tens(3, 3, crys%mxdatm * crys%ntype)  
  !
  !     DESCRIPTION:
  !     -----------
  !
  !     computes the ewald field on each atomd due to the presence of
  !     a dipole array. At each atom, there is a dipole attached of
  !     strength crys%dipole(ntype), and all of the dipoles are
  !     oriented in the same, but arbitrary direction (this direction
  !     is the first tensor index). The dipole array gives rise to a
  !     field, which is in the direction given by the second tensor
  !     index, i.e.:
  !
  !     Field(i) = Sum_j ewald_tens(j,i) * dipole_orientation(j)
  !
  !     In case of the magnetic field, the dipole orientation is along
  !     the external applied magnetic field, and the individual dipole
  !     strengths are given by crys%dipole(ntype) from an atomic
  !     calculation.
  !
  !     The summation used here is essentially from Akio Honma,
  !     "Science of Light", vol25, No. 2 p 41 (1976). It is equation (2
  !     .37), but without the factor of 4pi/3 P up front, and
  !     generalized to give the dipole field at an arbitrary dipole
  !     site in the unit cell.
  !
  !     To get the physical field at the dipole site, one has to add:
  !
  !     4pi/3 * polarization  + macroscopic field    (electric field)
  !
  !     -8pi/3 * magnetization + macroscopic field   (magnetic field)
  !
  !     So in the sense of Jackson, we have computed E_near, without
  !     the Lorentz field.
  !
  !     ---------------- local variables ---------------------------------
  !
  integer :: alp, bet, &                                 ! atom basis indices
       imx, jmx, kmx, ieta, &
       natot, &                 ! total number of atoms, no matter which kind
       i, j, k, l, m, ir(3), im2, jm2, km2, itot, idum
  real(dp), parameter :: small = 1.0d-12, &
       etalist(6) = (/ 2d-2, dhalf, done, dtwo, dfive, 1.0d1 /)
  real(dp) :: tol, t0, &
       eta, &                                    ! corresponds to Kittels eta
       seta, &                                                    ! sqrt(eta)
       aconst, kfundivx2, kfunthird, pi4invcell, inv4eta, inv2pi, &
       rmax, xsi, kfun, invg2, xseta, xseta2, &
       x, x2, x5, invx5, ep, er, f, g, ralp(3), &
       prefac, g_times_r_alp, gc(3)
  real(dp), external :: gimmetime, myderfc
  complex(dp) :: sumalp, dum  
  complex(dp), allocatable :: stfac(:), tens(:,:,:)  

  !     work arrays
  !
  real(dp), allocatable :: dipole(:), &  ! size dipole(crys%ntype*crys%mxdatm)
       rc(:,:)                             ! size rc(3,crys%ntype*crys%mxdatm)
  !
  !     variables for the gspace loop:
  !
  integer :: igs, iord, irod, igv(4), igv3, ffth(4), fftn(4)  
  real(dp) :: gv(3)  
  !     ------------------------------------------------------------------

  t0 = gimmetime()  
  allocate(rc(3, crys%ntype * crys%mxdatm)) ; rc = dzero
  allocate(dipole(crys%ntype * crys%mxdatm)) ; dipole = dzero  
  !
  !     compute various parameters
  !
  tol = -log(small)  

  eta = dtwo * gs%gmax / (dfour * tol)  
  seta = sqrt(eta)  

  xsi = seta  

  aconst = dfour * seta**3 / sqrt(pi)  
  rmax = sqrt(tol / eta)  
  imx = int(rmax * sqrt(crys%bdot(1, 1)) / pi2) + 1  
  jmx = int(rmax * sqrt(crys%bdot(2, 2)) / pi2) + 1  
  kmx = int(rmax * sqrt(crys%bdot(3, 3)) / pi2) + 1  
  pi4invcell = pi4 / crys%vcell  
  inv4eta = done / (dfour * eta)              ! 1/(4 eta) for the gspace decay
  inv2pi = done / pi2  
  !
  !     ----------------- store data for atoms into new arrays -----------
  !
  natot = 0  
  do i = 1, crys%ntype  
     do j = 1, crys%natom(i)  
        natot = natot + 1  
        rc(:, natot) = crys%rat(:, j, i)  
        dipole(natot) = crys%dipole(i)  
     end do
  end do
  allocate(stfac(natot))  
  allocate(tens(3, 3, natot))  
  !     ------------------------------------------------------------------
  !                          DO SUM IN GSPACE
  !
  !     loop over gspace
  fftn(1:3) = gs%fftsize(1:3)  
  fftn(4) = gs%fftsize(3)  
  ffth(:) = fftn(:) / 2  
  tens = zzero
  igs = 0  
  do iord = 1, gs%lorder                             ! loop through x/y gspace
     irod = gs%order(1, iord)  
     igv(1) = irod / gs%fftsize(2)  
     igv(2) = mod(irod, gs%fftsize(2))  
     igv(3) = gs%order(2, iord)  
     igv(4) = gs%order(3, iord)  
     igv(:) = mod(igv(:) + ffth(:), fftn(:)) - ffth(:)  
     gv(1:2) = real(igv(1:2), dp)  
     do igv3 = igv(3), igv(4)                               ! loop over z axis
        gv(3) = real(igv3, dp)  
        igs = igs + 1  
        if (gs%ekin(igs) > dzero) then  
           !
           !              compute exp(-G^2/(4 eta)) G/|G|^2
           !
           prefac = exp(-gs%ekin(igs) * inv4eta)  
           !
           !              precompute the "structure factors" stfac. Here,
           !              rc are the relative coordinates MULTIPLIED BY TWO PI!
           !
           do alp = 1, natot  
              g_times_r_alp = gv(1) * rc(1, alp) + gv(2) * rc(2, alp) + &
                   gv(3) * rc(3, alp)                    ! rc has fac of twopi
              stfac(alp) = cmplx(cos(g_times_r_alp), sin(g_times_r_alp), dp)
           end do
           !
           !              now run the inner double loop over the pairs of atoms
           !
           !              tens(external_field, induced_field,atom)
           sumalp = zzero
           do alp = 1, natot  
              sumalp = sumalp + conjg(stfac(alp)) * dipole(alp)  
           end do
           sumalp = sumalp * prefac  
           gc = matmul(crys%bvec, gv)  
           invg2 = done / gs%ekin(igs)  
           do bet = 1, natot  
              dum = sumalp * stfac(bet)  
              do i = 1, 3  
                 tens(i,:, bet) = tens(i,:, bet) - gc(i) * gc(:) * dum * invg2
                 tens(i, i, bet) = tens(i, i, bet) + dthird * dum  
              end do
           end do
        end if
     end do
  end do
  !
  !     put in the factor of -4pi/cell_volume
  !
  tens = pi4invcell * tens  
  !
  !     --------------------- now the realspace sum ----------------------
  !
  !     Notice that the lattice sum does not scale with the size of the
  !     system, i.e. the number of R vectors summed up is constant.
  !     What does scale is the innermost double loop over the atoms,
  !     and hence the routine is of complexity N**2.
  ! remove factor of two pi

  rc(:, 1:natot) = rc(:, 1:natot) * inv2pi  
  im2 = 2 * imx + 1  
  jm2 = 2 * jmx + 1  
  km2 = 2 * kmx + 1  
  idum = 0
  do itot = gs%myproc, im2 * jm2 - 1, gs%nproc  
     i = itot / jm2 + 1  
     j = mod(itot, jm2) + 1  
     ir(1) = i - imx - 1  
     ir(2) = j - jmx - 1  
     do k = 1, km2  
        ir(3) = k - kmx - 1  
        if ((ir(1) /= 0) .or. (ir(2) /= 0) .or. (ir(3) /= 0)) then
           idum = idum + 1  
           do bet = 1, natot  
              do alp = 1, natot  
                 !                    rc is the relative lattice coordinate
                 ralp(:) = rc(:, alp) - rc(:, bet) + real(ir(:), dp)  
                 ralp = matmul(crys%avec, ralp)  
                 x2 = ralp(1) * ralp(1) + ralp(2) * ralp(2) + ralp(3) * ralp (3)
                 x = sqrt(x2)  
                 ep = exp(-eta * x2)  
                 er = myderfc(seta * x)  
                 xseta = x * seta  
                 xseta2 = x2 * eta  
                 kfun = (done + dtrhalf / xseta2) * ep + &
                      dthree * sqrt(pi) / (dfour * xseta * xseta2) * er
                 !
                 !                       add to the tensor:
                 !
                 !
                 kfun = kfun * aconst * dipole(alp)  
                 kfundivx2 = kfun / x2  
                 kfunthird = kfun * dthird  
                 do l = 1, 3  
                    tens(l,:, bet) = tens(l,:, bet) + ralp(l) * ralp(:) * &
                         kfundivx2
                    tens(l, l, bet) = tens(l, l, bet) - kfunthird  
                 end do
              end do
           end do
        end if
     end do
  end do

  ! sum across all processors
  call all_sum_all(tens, natot * 9)  
  !
  !     treat the atoms inside the unit cell, i.e. R=0
  !
  do bet = 1, natot  
     do alp = 1, natot  
        if (alp /= bet) then  
           ralp(:) = rc(:, alp) - rc(:, bet)  
           ralp = matmul(crys%avec, ralp)  
           x2 = ralp(1) * ralp(1) + ralp(2) * ralp(2) + &
                ralp(3) * ralp(3)
           x = sqrt(x2)  
           ep = exp (-eta * x2)  
           er = myderfc(seta * x)  
           xseta = x * seta  
           xseta2 = x2 * eta  
           kfun = (done + dtrhalf / xseta2) * ep + dthree * sqrt(pi) / &
                (dfour * xseta * xseta2) * er
           kfun = kfun * aconst * dipole(alp)  
           kfundivx2 = kfun / x2  
           kfunthird = kfun * dthird  
           do l = 1, 3  
              tens(l,:, bet) = tens(l,:, bet) + ralp(l) * ralp(:) * kfundivx2
              tens(l, l, bet) = tens(l, l, bet) - kfunthird  
           end do
        end if
     end do
  end do
  ewald_tens = dzero  
  ewald_tens(:,:, 1:natot) = real(tens(:,:, 1:natot), dp)  
  write(9, 110)  
  do i = 1, crys%ntype  
     write(9, 120) crys%nameat(i), crys%dipole(i)  
  end do

  if (iand(ioutput, 8) == 8) write(9, 940) gimmetime() - t0  
  call myflush(9)

  deallocate(tens)  
  deallocate(stfac)  
  deallocate(dipole)  
  deallocate(rc)  

  return  

110 format(/' ATOM TYPE       MAGNETIC POLARIZABILITY IN A.U.'/)  

120 format(' ',a2,25x,f12.7)  
910 format(' ------------- EWALD DIPOLE TENSORS --------------')  
920 format(' ATOM NUMBER:', i5)  

930 format(3(10x,3f20.10/))  

940 format(/' TIME FOR EWALD DIPOLE CALCULATION:',f12.3)  

end subroutine ewald_dipole
