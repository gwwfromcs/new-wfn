!
subroutine electric_field(pw_params, crys, den, vion, vin, gs, &
     nlineplot, line_plot, bra, ket)
  !
  use all_to_all_module  
  include 'use.h'  
  implicit none                     ! implicit? Just say no!
  include 'interface.h'  
  include 'all_to_all.h'  
  !
  !     computes and plots electric field
  !
  !
  !     INPUT:
  !     -----
  !
  integer, intent(in) :: nlineplot               ! number of line plots to do
  real(dp), intent(in) :: bra(3), &             ! for tensor data left vector
       ket(3), &                               ! for tensor data right vector
       line_plot(7, *)                  ! line plot ((start,end,bins), lines)
  type(crystal), intent(in) :: crys                 ! for the lattice vectors
  type(parallel_gspace), intent(in) :: gs                  ! potential gspace
  type(pw_parameter), intent(in) :: pw_params         ! plane-wave parameters
  complex(dp), intent(in) :: den(gs%length, crys%nspin), &   ! charge density
       vion(gs%length, crys%nspin), &              ! screened ionic potential
       vin(gs%length, crys%nspin)                              ! vxc+vhartree
  !
  !
  !     1997 Bernd Pfrommer, UCB
  !
  !     computes the electric field due to the ionic and hartree
  !     potential. The electric field is
  !
  !     E = grad V,
  !
  !     where V is V_hartree + V_ion. There is no minus sign, because
  !     the electronic potential is the negative of the electric
  !     potential.
  !
  !     The spin polarized version is programmed, but NOT DEBUGGED!
  !
  !
  !     ---------------- local variables ---------------------------------
  !
  integer :: i, iat, nt, idir, is  
  real(dp) :: r(3), t0, fac, g, gr, gvc(3)  
  complex(dp) :: td, gvfac(3)  
  complex(dp), allocatable :: v(:), egs(:,:), efield(:,:,:)
  real(dp), external :: gimmetime  
  !
  !     ----- variables for the gspace loop ---------
  !
  integer :: igv(4), fftn(4), ffth(4), igv3, irod, iord, igs  
  real(dp) :: gv(3)                                             ! if required
  !
  !     ----------------------------------------------------------
  !
  fac = dfour * pi2 / crys%vcell  
  if (crys%nspin > 1) then  
     write(9, *) 'SPIN POLARIZED VERSION OF ELECTRIC_FIELD'  
     write(9, *) 'IS NOT DEBUGGED!'  
  end if
  t0 = gimmetime()  
  !
  !     first compute the bare ionic potential
  !
  allocate(v(gs%length))  
  v(1:gs%length) = vion(1:gs%length, 1) - vin(1:gs%length, 1)  
  !
  !     now add the hartree potential
  !
  if (crys%nspin > 1) then  
     do i = 1, gs%length  
        if (gs%ekin(i) > dzero) then  
           td = fac * (den(i, 1) + den(i, 2)) / gs%ekin(i)  
           v(i) = v(i) + td  
        end if
     end do
  else  
     do i = 1, gs%length  
        if (gs%ekin(i) > dzero) then  
           v(i) = v(i) + fac * den(i, 1) / gs%ekin(i)  
        end if
     end do
  end if
  !
  !     take the derivative in gspace
  !
  !     E = - gradient (V_hartree + V_ion)
  !
  allocate(egs(gs%length, 3))  
  do i = 1, 3  
     call derivative_gspace(i, v(1), gs, crys, egs(1, i))  
  end do
  !
  !     to get the values at the atoms, run the fourier sum
  !     explicitly rather than FFT
  !
  allocate(efield(3, crys%mxdatm, crys%ntype))  

  efield = zzero
  fftn(1:3) = gs%fftsize(1:3)  
  fftn(4) = gs%fftsize(3)  

  ffth(:) = fftn(:) / 2  
  do idir = 1, 3  
     do nt = 1, crys%ntype  
        do iat = 1, crys%natom(nt)  
           r = crys%rat(:, iat, nt)  
           igs = 0  
           do iord = 1, gs%lorder                   ! loop through x/y gspace
              irod = gs%order(1, iord)  
              igv(1) = irod / gs%fftsize(2)  
              igv(2) = mod(irod, gs%fftsize(2))  
              igv(3) = gs%order(2, iord)  
              igv(4) = gs%order(3, iord)  
              igv(:) = mod(igv(:) + ffth(:), fftn(:)) - ffth(:)  
              gv(1:2) = real(igv(1:2), dp) + gs%rk(1:2)  
              do igv3 = igv (3), igv (4)                   ! loop over z axis
                 gv(3) = real(igv3, dp) + gs%rk(3)  
                 igs = igs + 1  
                 efield(idir, iat, nt) = efield(idir, iat, nt) + &
                      egs (igs, idir) * exp(cmplx(dzero, &
                      (r(1) * gv(1) + r(2) * gv(2) + r(3) * gv(3)), dp))
              end do
           end do
        end do
     end do
  end do

  call all_sum_all(efield, 3 * crys%mxdatm * crys%ntype)  
  write(9, 100)  
  do nt = 1, crys%ntype  
     do iat = 1, crys%natom(nt)  
        write(9, 1254) crys%nameat(nt), iat, dzero, -aimag(efield(:, iat, nt))
     end do
  end do

1254 format(1x,a2,i3,2x,f9.3,3x,3f12.6)  

  call lineplot(2, crys, dmone, egs(1, 1), gs, nlineplot, &
       line_plot, bra, ket, 3, 'E_FIELD')

  call myflush(9)  
  !     do an intergration of the MT sphere
  !

  efield = zzero
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
     gv(1:2) = real(igv(1:2), dp) + gs%rk(1:2)  
     do igv3 = igv(3), igv(4)                               ! loop over z axis
        igs = igs + 1  
        g = sqrt(gs%ekin(igs))  
        if (g > dzero) then  
           gv(3) = real(igv3, dp) + gs%rk(3)  
           gvc = matmul(crys%bvec, gv)  
           do nt = 1, crys%ntype  
              gr = g * crys%mtsphere(nt)  
              gvfac = gvc * v(igs) * (sin(gr) - gr * cos(gr)) / gr**3  
              do iat = 1, crys%natom(nt)  
                 r = crys%rat(:, iat, nt)  
                 efield(:, iat, nt) = efield(:, iat, nt) + &
                      gvfac * exp(cmplx(dzero, &
                      (r(1) * gv(1) + r(2) * gv(2) + r(3) * gv(3)), dp))
              end do
           end do
        end if
     end do
  end do
  efield = cmplx(dzero, dthree, dp) * efield  

  call all_sum_all(efield, 3 * crys%mxdatm * crys%ntype)  
  write(9, 110)  
  do nt = 1, crys%ntype  
     do iat = 1, crys%natom(nt)  
        write(9, 1254) crys%nameat(nt), iat, crys%mtsphere(nt), &
             real(efield(:, iat, nt), dp)
     end do
  end do

  write(9, 120)  
  deallocate(v)  
  deallocate(egs)  
  deallocate(efield)  

  if (iand(pw_params%output(1), 8) == 8) write(9, 940) gimmetime() - t0

  return

100 format(/' -------- ELECTRIC FIELD ----------'/, &
       &     1x,' Atom     Radius',9x,'E x',9x,'E y',9x,'E z')
110 format(/' ---- AVERAGED ELECTRIC FIELD -----'/, &
       &     1x,' Atom     Radius',9x,'E x',9x,'E y',9x,'E z')
120 format(/' Electric field must be multiplied by 8.5763*10^6', &
       &     ' to convert to statvolt/cm')

940 format(/' TIME FOR ELECTRIC FIELD:',f12.3)  

end subroutine electric_field
