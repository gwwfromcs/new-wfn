!
subroutine lukman(crys, gs, pot)

  include 'use.h'  
  implicit none              ! never uncomment this line.
  include 'interface.h'
  !
  !     1997 Bernd Pfrommer, Paul Delaney
  !
  !
  !
  !     INPUT:
  !     -----
  !
  type(crystal), intent(in) :: crys  
  type(parallel_gspace), intent(in) :: gs  
  !
  !     OUTPUT:
  !     ------
  !
  complex(dp), intent(out) :: &
       pot(gs%length)                      ! the effective potential in gspace
  !
  !     DESCRIPTION:
  !     -----------
  !
  !     sets up a purely local potential for model calculations
  !
  !
  !
  !     ---------------- local variables ----------------------
  !
  integer :: igv(3), nparam, igs, i  
  integer, allocatable :: ipotgv(:,:)  

  integer :: atomnum, NUMATOMS  
  real(dp) :: a, b, c, d, temp, mu, V, sf, sfr, sfi, q, x, y, z, &
       GdotTau, sinqona
  real(dp) :: depth, width  
  real(dp) :: volratio  
  real(dp), allocatable :: coeff(:)  
  real(dp) :: taux(1000), tauy(1000), tauz(1000)  
  !
  !     ---------------------------------------------------------------
  !
  if (.not. gs%igvec) then  
     write(9, *) 'lukman: gvectors not set up!'  
     call mystop  
  end if

  write(0, *) 'gmax is ', gs%gmax  
  write(0, *) 'rk is ', gs%rk(1) , gs%rk(2) , gs%rk(3)  
  write(0, *) 'length is ', gs%length  
  write(0, *) 'totlength is ', gs%totlength  
  write(0, *) 'nstar is ', gs%nstar  
  write(0, *) 'ig0 is ', gs%ig0  
  write(0, *) 'igvec is ', gs%igvec  

  write(0, *) 'iekin is ', gs%iekin  

  pot = zzero    ! clear potential
  open(17, file = 'input.lukman', status = 'old')  
  rewind(17)  
  read(17, *) nparam  
  if (nparam > 0) then  
     allocate(coeff(nparam))  
     allocate(ipotgv(3, nparam))  
     do i = 1, nparam  
        read(17, *) ipotgv(:, i), coeff(i)  
     end do
  else if (nparam == 0 .or. nparam == -1) then  
     read(17, *) a, b, c  
     !
     !     set up the coefficients according to a,b,c
     !     nparam = 0 means to the simple cubic: <0 means bcc.
     write(9, *) ' FOUND A,B,C, SETTING UP PARAMETERS'  

     write(9, *) 'a=', a, 'b=', b, 'c=', c  

     allocate(coeff(19))
     allocate(ipotgv(3, 19))  
     if (nparam == 0) then  
        !     This stuff is for the larger cell: simple cubic.
        write(9, *) 'lukman: doing the setup for the larger cell:', &
             'simple cubic'
        coeff(1) = (a * a + b * b + c * c) * dhalf 

        ipotgv(:, 1) = 0
        coeff(2:5) = a * b * dhalf
        ipotgv(1, 2) = 1
        ipotgv(2, 2) = 1
        ipotgv(3, 2) = 0  
        ipotgv(1, 3) = -1
        ipotgv(2, 3) = 1
        ipotgv(3, 3) = 0  
        ipotgv(1, 4) = 1
        ipotgv(2, 4) = -1
        ipotgv(3, 4) = 0  
        ipotgv(1, 5) = -1
        ipotgv(2, 5) = -1
        ipotgv(3, 5) = 0
        coeff(6:9) = a * c * dhalf
        ipotgv(1, 6) = 1
        ipotgv(2, 6) = 0
        ipotgv(3, 6) = 1  
        ipotgv(1, 7) = -1
        ipotgv(2, 7) = 0
        ipotgv(3, 7) = 1  
        ipotgv(1, 8) = 1
        ipotgv(2, 8) = 0
        ipotgv(3, 8) = -1  
        ipotgv(1, 9) = -1
        ipotgv(2, 9) = 0
        ipotgv (3, 9) = -1
        coeff(10:13) = b * c * dhalf
        ipotgv(1, 10) = 0
        ipotgv(2, 10) = 1
        ipotgv(3, 10) = 1
        ipotgv(1, 11) = 0
        ipotgv(2, 11) = -1
        ipotgv(3, 11) = 1
        ipotgv(1, 12) = 0
        ipotgv(2, 12) = 1
        ipotgv(3, 12) = -1
        ipotgv(1, 13) = 0
        ipotgv(2, 13) = -1
        ipotgv(3, 13) = -1
        coeff(14:19) = a * a * dqtr
        ipotgv(1, 14) = 2
        ipotgv(2, 14) = 0
        ipotgv(3, 14) = 0
        ipotgv(1, 15) = -2
        ipotgv(2, 15) = 0
        ipotgv(3, 15) = 0
        ipotgv(1, 16) = 0
        ipotgv(2, 16) = 2
        ipotgv(3, 16) = 0
        ipotgv(1, 17) = 0
        ipotgv(2, 17) = -2
        ipotgv(3, 17) = 0
        ipotgv(1, 18) = 0
        ipotgv(2, 18) = 0
        ipotgv(3, 18) = 2
        ipotgv(1, 19) = 0
        ipotgv(2, 19) = 0
        ipotgv(3, 19) = -2

        nparam = 19  

     else if (nparam == -1) then  
        !     This stuff is for the smaller cell:bcc lattice.
        !     Invert the potential from the one lukman previously gave me (repla
        !     V by -V) to see if we get more localisation in real-space, as the
        !     old potential had minima which were very soft.

        write(9, *) 'lukman: doing the setup for the bcc lattice...'  
        coeff(1) = - (a * a + b * b + c * c) * dhalf
        ipotgv(:, 1) = 0
        coeff(2:5) = -a * b * dhalf
        ipotgv(1, 2) = 0
        ipotgv(2, 2) = 0
        ipotgv(3, 2) = 1  
        ipotgv(1, 3) = 0
        ipotgv(2, 3) = 0
        ipotgv(3, 3) = -1  
        ipotgv(1, 4) = -1
        ipotgv(2, 4) = 1
        ipotgv(3, 4) = 0  
        ipotgv(1, 5) = 1
        ipotgv(2, 5) = -1
        ipotgv(3, 5) = 0  
        coeff(6:9) = -a * c * dhalf
        ipotgv(1, 6) = 0
        ipotgv(2, 6) = 1
        ipotgv(3, 6) = 0  
        ipotgv(1, 7) = 0
        ipotgv(2, 7) = -1
        ipotgv(3, 7) = 0  
        ipotgv(1, 8) = 1
        ipotgv(2, 8) = 0
        ipotgv(3, 8) = -1  
        ipotgv(1, 9) = -1
        ipotgv(2, 9) = 0
        ipotgv(3, 9) = 1  
        coeff (10:13) = -b * c * dhalf
        ipotgv(1, 10) = 1
        ipotgv(2, 10) = 0
        ipotgv(3, 10) = 0
        ipotgv(1, 11) = -1
        ipotgv(2, 11) = 0
        ipotgv(3, 11) = 0
        ipotgv(1, 12) = 0
        ipotgv(2, 12) = 1
        ipotgv(3, 12) = -1
        ipotgv(1, 13) = 0
        ipotgv(2, 13) = -1
        ipotgv(3, 13) = 1
        coeff(14:15) = -a * a * dqtr
        ipotgv(1, 14) = -1
        ipotgv(2, 14) = 1
        ipotgv(3, 14) = 1
        ipotgv(1, 15) = 1
        ipotgv(2, 15) = -1
        ipotgv(3, 15) = -1
        coeff(16:17) = -b * b * dqtr
        ipotgv(1, 16) = 1
        ipotgv(2, 16) = -1
        ipotgv(3, 16) = 1
        ipotgv(1, 17) = -1
        ipotgv(2, 17) = 1
        ipotgv(3, 17) = -1
        coeff(18:19) = -c * c * dqtr
        ipotgv(1, 18) = 1
        ipotgv(2, 18) = 1
        ipotgv(3, 18) = -1
        ipotgv(1, 19) = -1
        ipotgv(2, 19) = -1
        ipotgv(3, 19) = 1

        nparam = 19  
     end if

  end if
  if (nparam /= -2) then  
     write(9, '('' NUMBER OF LUKMAN GVECTORS SET:'', i4,/)') nparam
     do i = 1, nparam  
        write(9, '(3i4,f12.6)') ipotgv(:, i), coeff(i)  
     end do
     do igs = 1, gs%length  
        igv = gs%gvec(:, igs)  
        do i = 1, nparam  
           if (igv(1) == ipotgv(1, i) .and. igv(2) == ipotgv(2, i) .and. &
                igv(3) == ipotgv(3, i)) then
              pot(igs) = cmplx(coeff(i), dzero, dp)  
              write(9, *) igs, igv, pot(igs)  
           end if
        end do
     end do
     deallocate(coeff)  
     deallocate(ipotgv)  
  end if
  !     Now see if we're doing tube stuff
  if (nparam ==  -2) then  
     read(17, *) a  
     read(17, *) b  
     read(17, *) temp  
     read(17, *) d  
     read(17, *) mu  
     read(17, *) depth  
     read(17, *) width  
     write(0, *) 'a,b,temp,d,mu=', a, b, temp, d, mu  
     write(0, *) 'depth,width=', depth, width  
     write(0, *) 'have ', crys%ntype, ' types of atoms'  
     write(0, *) 'have ', crys%natom (1) , ' atoms of the first type'
     do atomnum = 1, crys%natom(1)  
        write(0, *) 'atom ', atomnum, ' has lattice coords', &
             crys%rat(1, atomnum, 1) / pi2, crys%rat(2, atomnum, 1) / pi2, &
             crys%rat(3, atomnum, 1) / pi2
     end do

     volratio = 76.8383d0 / crys%vcell  
     if (.not. gs%iekin) then  
        write(0, *) 'Oops: lukman: kinetic energies not available...'  
        write(0, *) 'Stopping...'  
        write(9, *) 'Oops: lukman: kinetic energies not available...'  
        write(9, *) 'Stopping...'  
        call mystop  
     end if
     if (gs%rk(1) /= dzero .or. gs%rk(2) /= dzero .or. &
          gs%rk(3) /= dzero) then
        write(0, *) 'Cannot use ekin because rk is ', gs%rk(1), gs%rk(2), &
             gs%rk(3)
        write(0, *) 'stopping...'  
        write(9, *) 'Cannot use ekin because rk is ', gs%rk(1), gs%rk(2), &
             gs%rk(3)
        write(9, *) 'stopping...'  
        call mystop  
     end if
     do igs = 1, gs%length  
        igv = gs%gvec(:, igs)  
        q = gs%ekin(igs)  
        q = sqrt(q)  
        sinqona = sin(q / a)
        V = (d * sinqona**3 + b) / (exp((q - mu) / temp) + done)
        V = V + depth * exp(-q * q / (width * width))  
        sfr = dzero
        sfi = 0.0d0  
        do atomnum = 1, crys%natom(1)  
           GdotTau = (real(igv(1), dp) * crys%rat(1, atomnum, 1) + &
                real(igv(2), dp) * crys%rat(2, atomnum, 1) + &
                real(igv(3), dp) * crys%rat(3, atomnum, 1))
           sfr = sfr + cos(-GdotTau)  
           sfi = sfi + sin(-GdotTau)  
        end do
        sfr = sfr / crys%natom(1)  
        sfi = sfi / crys%natom(1)  
        if (abs(sfi) >= 0.0001d0) then  
           write(9, *) ' Oops - trouble in calculating the ', &
                'structure factors'
           write(9, *) 'for G-vector (', igv(1), igv(2), igv(3), &
                ') ', 'the real/imag parts are ', sfr, sfi
           call mystop  
        end if
        sf = sfr  
        pot(igs) = cmplx(V * sf * volratio * crys%natom(1) * dhalf, &
             dzero, dp)
     end do
  end if

  close(17) 
 
  return  

end subroutine lukman
