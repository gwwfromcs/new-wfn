!     @process extchk
!
subroutine plot_bandstruc(kpoints, bands, crys, energs, bslabels)  
  !
  include 'use.h'  
  implicit none  
  include 'interface.h'  
  !
  !     INPUT:
  !     -----
  !
  type(kpoint), intent(in) :: kpoints  
  type(band), intent(in) :: bands  
  type(crystal), intent(in) :: crys  
  type(energy), intent(in) :: energs  
  character(len=100), intent(in) :: bslabels(*)  ! names of the k-points
  !
  !     DESCRIPTION:
  !     -----------
  !
  !     plots band structure, with correct metric, into a file
  !     The file can be read for instance with xmgr
  !
  !     note that the label field is used to avoid jumps when
  !     two different BZ lines follow each other. this works
  !     because kpoints%label(k)=0.0 in that case.
  !
  !
  !
  !     ------------------ local variables ------------------------
  !
  integer :: is, k, n, nregion, i, ii, l, itick  
  character(len=100) :: dummy, labs(3)  
  real(dp) :: efermi(2), &
       plabs(3), &                    ! positions of the three labels
       s, d, dk(3), slab, smid, emin, emax, xmax
  real(dp), parameter :: &
       tmarg = done, bmarg = done     ! top and bottom margins
  real(dp), parameter :: vxmin = 0.15d0, vymin = 0.15d0, &
       vxmax = 0.85d0, vymax = 0.85d0  

  open(11, file = 'BANDSTRUC', form = 'formatted', status = 'unknown')
  open(12, file = 'band.dat', form = 'formatted', status = 'unknown')
  !     efermi=energs%efermi
  efermi(1) = energs%inputfermi  
  efermi(2) = energs%inputfermi  
  !
  !     figure out max and minimum band energy
  !
  emin = 1.0d6  
  emax = -1.0d6  
  do is = 1, bands%nspin  
     do k = 1, kpoints%nrk  
        if (maxval(bands%energy(1:bands%nband(k, is), k, is)) + &
!             energs%evalshift(is) - efermi(is) > emax) emax = &
            energs%vxc0(is) - efermi(is) > emax) emax = &
             maxval(bands%energy(1:bands%nband(k, is), k, is)) + &
!             energs%evalshift(is) - efermi(is)
            energs%vxc0(is) - efermi(is)
        if (minval(bands%energy(1:bands%nband(k, is), k, is)) + &
!             energs%evalshift(is) - efermi(is) < emin) emin = &
            energs%vxc0(is) - efermi(is) < emin) emin = &
             minval(bands%energy(1:bands%nband(k, is), k, is)) + &
!             energs%evalshift(is) - efermi(is)
           energs%vxc0(is) - efermi(is)
     end do
  end do
  emin = emin * ryd
  emax = emax * ryd
  !
  !     find maximum abscissa value
  !
  do is = 1, bands%nspin  
     s = dzero
     dk = dzero
     do k = 1, kpoints%nrk  
        if (k < kpoints%nrk) dk = kpoints%rk(:, k + 1) - kpoints%rk(:, k)
        d = kpoints%label(k) * sqrt(dot_product(dk, matmul(crys%bdot, dk)))
        s = s + d
     end do
  end do
  xmax = s  
  !
  !     write first part of the header
  !
  write(11, 5)  
  if (bands%nspin > 1) then
     write(11, 6) efermi(1:bands%nspin) * ryd
  else
     write(11, 66) efermi(1) * ryd
  endif
!  write(11, 7) energs%evalshift(1:bands%nspin) * ryd
  write(11, 7) energs%vxc0(1:bands%nspin) * ryd 
  write(11, 10) xmax, real(int(emin - bmarg), dp), &
       real(int(emax + tmarg), dp), vxmin, vxmax, vymin, vymax
  !
  !     write the labels
  !
  call myflush(9)  
  dk = dzero
  nregion = 1  
  slab = dzero
  s = dzero
  itick = 0  
  do k = 1, kpoints%nrk  
     if (k < kpoints%nrk) dk = kpoints%rk(:, k + 1) - kpoints%rk(:, k)
     d = kpoints%label(k) * sqrt(dot_product(dk, matmul(crys%bdot, dk)))
     if (d == dzero) then       ! switched to new region
        nregion = nregion + 1  
        !
        !           decompose into left/middle/right
        !
        dummy = trim(adjustl(bslabels(nregion - 1)))  
        l = len(dummy)  
        if (l <= 0) exit          ! no label written, continue with do loop
        is = 0  
        do i = 1, 3  
           ii = index(dummy(is + 1:l), '_') + is  
           if (ii <= is) then  
              labs(i) = dummy(is + 1:index(dummy, ' '))  
           else  
              labs(i) = dummy(is + 1:ii - 1)  
           end if
           if (ii <= is) exit  
           is = ii             ! start new search after '_'
        end do
        if (i == 1) then  
           plabs(1) = slab + (s - slab) * dhalf
        else if (i == 2) then  
           plabs(1) = slab
           plabs(2) = s  
        else  
           plabs(1) = slab
           plabs(2) = slab + (s - slab) * dhalf
           plabs(3) = s
        end if
        do l = 1, i  
           write(11, 30) itick, plabs(l), itick, trim(labs(l))  
           itick = itick + 1  
        end do
        slab = s  
     end if
     s = s + d
  end do
  write(11, 40) itick  
  !
  !     output the lines
  !
  !      write(11,19)
  dk = dzero
  nregion = 1  
  slab = dzero
  s = dzero
  do k = 1, kpoints%nrk  
     if (k < kpoints%nrk) dk = kpoints%rk(:, k + 1) - kpoints%rk(:, k)
     d = kpoints%label(k) * sqrt(dot_product(dk, matmul(crys%bdot, dk)))
     if (d == dzero) then       ! switched to new region
        nregion = nregion + 1  
        !           draw the vertical line
        !            write(11,20) s,dble(int(emin-bmarg)),s,
        !     $           dble(int(emax+tmarg))
        write(11, 20) s / xmax * (vxmax - vxmin) + vxmin, vymin, &
             s / xmax * (vxmax - vxmin) + vxmin, vymax
        slab = s  
     end if
     s = s + d
  end do
  !
  !     write the bands one by one
  !
  do is = 1, bands%nspin  
     do n = 1, bands%min(is)
        write(11, '(''@TYPE xy'')')  
        dk = dzero
        ii = 0  
        s = dzero
        do k = 1, kpoints%nrk  
           if (k < kpoints%nrk) dk = kpoints%rk(:, k + 1) - kpoints%rk(:, k)
           d = kpoints%label(k) * sqrt(dot_product(dk, matmul(crys%bdot, dk)))
           if (bands%nband(k, is) >= n) then  
              write(11, 100) s, (bands%energy(n, k, is) + &
!                  energs%evalshift(is) - efermi(is)) * ryd
                  energs%vxc0(is) - efermi(is)) * ryd
              write(12, 101) s/0.529177d0, bands%energy(n, k, is)*ryd
              ii = ii + 1  
           end if
           s = s + d
        end do
        write(11, '(''&'')')  
        write(12, *)
     end do
  end do

  close(11)  
  close(12)  

  write(9, 200)  

  return

5 format(3('#',/),'#      BANDSTRUCTURE FOR USE WITH XMGR', &
       &     3('#',/),'#    1996 Bernd Pfrommer')
6 format(3('#',/),'# THE EIGENVALUES ARE SHIFTED BY VXC(0), AND', &
       &     /'# ARE GIVEN WITH RESPECT TO THE INPUT FERMI LEVEL:', &
       &     2f20.8,' eV',/'# MAKE SURE THE INPUT FERMI IS CORRECT!')
66 format(3('#',/),'# THE EIGENVALUES ARE SHIFTED BY VXC(0), AND', &
       &     /'# ARE GIVEN WITH RESPECT TO THE INPUT FERMI LEVEL:', &
       &     1f20.8,' eV',/'# MAKE SURE THE INPUT FERMI IS CORRECT!')

7 format('# VXC0-SHIFT [eV]:',2f12.6)  

10 format('@with g0'/, &
       &     '@    world xmin 0'/ &
       &     '@    world xmax ',f12.6/, &
       &     '@    world ymin ',f12.6/, &
       &     '@    world ymax ',f12.6/, &
       &     '@    view xmin ',f12.6/, &
       &     '@    view xmax ',f12.6/, &
       &     '@    view ymin ',f12.6/, &
       &     '@    view ymax ',f12.6/, &
       &     '@    yaxis  tick on'/, &
       &     '@    yaxis  tick major 5'/, &
       &     '@    yaxis  tick minor 2.5'/, &
       &     '@    xaxis  tick major off'/, &
       &     '@    xaxis  tick minor off'/, &
       &     '@    xaxis  ticklabel on'/, &
       &     '@    xaxis  ticklabel type spec'/, &
       &     '@    xaxis  tick type spec')
19 format('@    line loctype world')  
20 format('@with line'/, &
       &     '@    line',f12.6,',',f12.6,',',f12.6,',',f12.6/, &
       &     '@line def')
30 format('@    xaxis  tick ',i2,', ',f12.6/, &
       &     '@    xaxis  ticklabel ',i2,', "',a,'"')

40 format('@    xaxis  tick spec ',i3)  
100 format(f12.6, 500(f13.6))  
101 format(f16.8, 500(f16.8))  

200 format(' WROTE BANDSTRUCTURE TO FILE BANDSTRUC!')

end subroutine plot_bandstruc
