!
subroutine gwreal(bands, kgs, wfn, kpoints)  
  !
  use all_to_all_module  
  include 'use.h'
  implicit none  
  include 'interface.h'  
  include 'all_to_all.h' 
  ! 
  !     include 'mpif.h'
  !
  !     INPUT:
  !     -----
  !
  type(complex_gspace_array), intent(inout) :: wfn    ! the wave functions
  type(band), intent(in) :: bands  
  type(parallel_gspace), intent(in) :: &
       kgs(bands%nrk)                   ! the gspaces at the different kpoints
  type(kpoint), intent(in) :: kpoints   ! all the kpoints
  !     ---------------- local variables ---------------------------------
  integer :: dimension_span, reduced_span, i, idum, ig, inc, inum, &
       irk, is, k, kprime, n, nprime, ng, wfninc, deg, mdeg
  integer, allocatable :: imap(:), inull(:), null_map(:,:)  
  real(dp) :: ii, rr, norm  
  complex(dp), allocatable :: psi(:,:), phi(:,:), scr1(:), scr2(:)
  !
  ! the code may crash after printing out " bands output for GW = "
  ! in file OUT if mdeg is set to a small number like this:
  ! integer, parameter :: mdeg = 10
  !
  ! let us determine the maximum degeneracy and store it in mdeg
  !
  mdeg = 1  
  do is = 1, wfn%nspin  
     do irk = 1, kpoints%nrk  
        do n = 1, bands%nband(irk, is) - 1  
           deg = 1  
           do nprime = n + 1, bands%nband(irk, is)  
              if (abs(bands%energy(n, irk, is) - &
                 bands%energy(nprime, irk, is)) < 1.0d-6 + 1.0d-10) &
                 deg = deg + 1  
           enddo  
           if (deg > mdeg) mdeg = deg  
        enddo  
     enddo  
  enddo  
  !
  ! then multiply mdeg by two to account for the
  ! real and imaginary parts of the wavefunction
  !
  mdeg = 2 * mdeg  

  allocate(imap(bands%max))  
  allocate(inull(bands%max))  
  allocate(null_map(mdeg, bands%max))
 
  do is = 1, wfn%nspin  
     do irk = 1, kpoints%nrk  
        ng = kgs(irk)%length  
        !
        ! figure out degeneracies
        ! store degeneracy information in imap(inum)
        !
        inum = 1  
        do n = 1, bands%nband(irk, is)  
           if (n == bands%nband(irk, is)) then  
              imap(inum) = n  
              inum = inum + 1  
              goto 10  
           end if
           if (abs(bands%energy(n, irk, is) - &
                bands%energy(n + 1, irk, is)) > 1.0d-6) then
              imap(inum) = n  
              inum = inum + 1  
           end if
10         continue  
        end do

        inum = inum - 1  
        !
        ! construct null vector map
        !
        inc = 1  
        allocate(scr1(ng))  
        allocate(scr2(ng))  
        do i = 1, inum  
           inull(i) = 1  
           do n = inc, imap(i)  
              do ig = 1, ng  
                 scr1(ig) = real(wfn%data(ig + (n - 1) * ng, irk, is), dp)  
                 scr2(ig) = aimag(wfn%data(ig + (n - 1) * ng, irk, is))  
              end do
              rr = parallel_zdotc(ng, scr1, 1, scr1, 1)  
              if (sqrt(rr) < 1.0d-2) null_map(inull(i), i) = 0  
              if (sqrt(rr) > 1.0d-2) null_map(inull(i), i) = 1  
              inull(i) = inull(i) + 1  
              !
              ii = parallel_zdotc(ng, scr2, 1, scr2, 1)  
              if (sqrt(ii) < 1.0d-2) null_map(inull(i), i) = 0  
              if (sqrt(ii) > 1.0d-2) null_map(inull(i), i) = 1  
              inull(i) = inull(i) + 1  
           end do
           inull(i) = inull(i) - 1  
           inc = imap(i) + 1  
        end do
        deallocate(scr1)  
        deallocate(scr2)  
        !
        ! i loops over each degenerate multiplet
        ! k loops over the index of the wavefunction in each
        !   multiplet counting real and imaginary parts
        !
        inc = 1  
        wfninc = 1  
        allocate(psi(ng, mdeg))  
        allocate(phi(ng, mdeg))  
        do i = 1, inum  
           kprime = 1  
           do k = 1, 2 * (imap(i) - inc) + 1, 2  
              if (null_map(k, i) == 1) then  
                 do ig = 1, ng  
                    phi(ig, kprime) = cmplx(real(wfn%data(ig + &
                         (wfninc - 1) * ng, irk, is), dp), dzero, dp)
                 end do
                 kprime = kprime + 1  
              end if
              if (null_map(k + 1, i) == 1) then  
                 do ig = 1, ng  
                    phi(ig, kprime) = cmplx(aimag(wfn%data(ig + &
                         (wfninc - 1) * ng, irk, is)), dzero, dp)
                 end do
                 kprime = kprime + 1  
              end if
              wfninc = wfninc + 1  
           end do
           dimension_span = kprime - 1  
           ! normalize basis of spanning functions
           do n = 1, dimension_span  
              norm = parallel_zdotc(ng, phi(:, n), 1, phi(:, n), 1)  
              do ig = 1, ng  
                 phi(ig, n) = phi(ig, n) / sqrt(norm)  
              end do
           end do
           ! perform gramm-schimdt process
           call gramm_schmidt(psi, reduced_span, phi, dimension_span, ng)  
           if (reduced_span < (imap(i) - inc + 1)) then  
              write(9, *) 'you really messed up n1 /= n2 ', reduced_span, &
                   (imap(i) - inc + 1)
              call mystop  
           end if
           !
           ! here change the wavefunctions to their real versions
           !
           do n = 1, (imap(i) - inc + 1)  
              do ig = 1, ng  
                 wfn%data(ig + (inc + n - 1 - 1) * ng, irk, is) = psi(ig, n)  
              end do
           end do
           !
           inc = imap(i) + 1  
        end do
        deallocate(psi)  
        deallocate(phi)  
     end do            ! end of k-point loop
  end do               ! end of spin loop
  deallocate(imap)  
  deallocate(inull)  
  deallocate(null_map)  
  ! end of gramm-schmidt process

end subroutine gwreal
!-----------------------------------------------------------------------
subroutine gramm_schmidt(psi, n_psi, phi, n_phi, npts)  
  !
  use constants
  use all_to_all_module  
  implicit none
  include 'all_to_all.h'  
  !
  integer, intent(in) :: n_phi, npts  
  integer, intent(out) :: n_psi
  complex(dp), intent(out) :: psi(npts, n_phi)
  complex(dp), intent(in) :: phi(npts, n_phi)  
  complex(dp), allocatable :: vec(:)  
  real(dp) :: norm  
  real(dp) :: scalar  
  integer :: i, l, ii 
  !
  allocate(vec(npts))  
  norm = dzero
  n_psi = 1  
  ! intialize
  norm = parallel_zdotc(npts, phi(:, 1), 1, phi(:, 1), 1)  
  norm = sqrt(norm)  
  do l = 1, npts  
     psi(l, 1) = phi(l, 1) / norm  
  end do
  ! begin gramm-scmidt algorithm
  do i = 1, n_phi - 1  
     ! project out
     do l = 1, npts  
        vec(l) = phi(l, i + 1)  
     end do
     do ii = 1, n_psi  
        scalar = parallel_zdotc(npts, phi(:, i + 1), 1, psi(:, ii), 1)
        do l = 1, npts  
           vec(l) = vec(l) - scalar * psi(l, ii)  
        end do
     end do
     norm = parallel_zdotc(npts, vec, 1, vec, 1)  
     norm = sqrt(norm)  
     ! if we have a non-null projection, set new psi
     if (norm > 1.0d-6) then  
        n_psi = n_psi + 1  
        do l = 1, npts  
           psi(l, n_psi) = vec(l) / norm  
        end do
     end if
  end do
  deallocate(vec)
  
end subroutine gramm_schmidt
