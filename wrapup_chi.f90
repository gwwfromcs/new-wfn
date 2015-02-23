!     @process extchk
!
subroutine wrapup_chi(crys, ffts, gs, rho, qmag, g0mask, chi_g, chi,writeJ)
  !
  include 'use.h'  
  implicit none        ! implicit? Just say no!
  include 'interface.h'  
  !
  !     INPUT:
  !     -----
  !
  type(crystal), intent(in) :: crys           ! for the lattice vectors
  type(fft_struc), intent(in) :: ffts  
  type(parallel_gspace), intent(in) :: gs     ! potential gspace
  complex(dp), intent(in) :: &
       chi(3, 3), &                           ! G=0 component
       rho(ffts%r_size, 3, 3)                 ! raw f(r, p0, p1)
  real(dp), intent(in) :: &
       g0mask(3, 3), &                        ! mask for g=0 component
       qmag                                   ! magnitude of q
  integer,intent(in) :: writeJ                ! >0 write out current
  !
  !     OUTPUT:
  !     ------
  !
  complex(dp), intent(out) :: &
       chi_g(gs%length, 3, 3)                 ! chi in gspace
  !
  !
  !     1996 Bernd Pfrommer, UCB
  !
  !
  !     ---------------- local variables ---------------------------------
  !
  integer :: i, j, ig, ia, ib, itwoq, iq, iubxq, p0, p1, info  
  real(dp) :: fac, uaxg(3), ua(3), ubxqsign, pi4vcell, &
       met(3, 3, 3)       ! metric : met(i, j, alpha) = (u_alpha x b_j)_i
  complex(dp), allocatable :: r_g(:,:,:)         ! rho_magnetic in gspace
  !
  !
  !
  allocate(r_g(gs%length, 3, 3), stat = info)  
  if (info /= 0) then  
     write(9, *) '*** ERROR: ALLOCATION FAILED IN WRAPUP_CHI!'  
     write(9, *) '*** CURE:  REDUCE MEMORY CONSUMPTION:'  
     call mystop  
  end if

  pi4vcell = dfour * pi / crys%vcell  
  !
  !     transform to fourier space
  !
  do p0 = 1, 3  
     do p1 = 1, 3  
        call fourier_transform(1, ffts, gs, r_g(1, p0, p1), rho(1, p0, p1), 1)
     end do
  end do
  if (writeJ > 0) then  !save the direct current
     Call writegsdat(1, gs, r_g(1,1,1), 9,9,1,'JGS',3)
  endif
  !
  !     compute metric tensor met(i,j,alpha) = (u_alpha x b_j)_i
  !
  do ia = 1, 3  
     ua = dzero
     ua(ia) = done
     do j = 1, 3  
        call myvecprod(met(1, j, ia), ua, crys%bvec(1, j))  
     end do
  end do
  !
  !     pick up contributions to chi in gspace.
  !
  do ig = 1, gs%length                      ! loop over gspace
     chi_g(ig, :,:) = zzero  
     if (gs%ekin(ig) > dzero) then  
        !
        !        conversion factor fac:
        !
        !        2   for spin up/down
        !        2   for Hartree-> Rydberg
        !        1/2 for averaging over the two q directions perpendicular to b0
        !        1/2q  because derivative is ( f(q)-f(-q) )/2q
        !        4pi because B = 4pi M
        !        1/vcell for the definition of the Fourier transform
        fac = pi4vcell / (qmag * gs%ekin(ig))  
        do ia = 1, 3          ! loop over alpha index
           uaxg = matmul(met(:,:, ia), real(gs%gvec(:, ig), dp))  ! uaxg = u_a
           do ib = 1, 3       ! loop over beta index
              do j = 1, 3                ! handle the u_alpha x G
                 chi_g(ig, ia, ib) = chi_g(ig, ia, ib) + uaxg(j) * &
                      r_g(ig, ib, j)
              end do
              !
              !                 normalize
              !
              chi_g(ig, ia, ib) = chi_g(ig, ia, ib) * fac  
           end do
        end do
     else  
        !
        !           put in the G=0 component. See Jackson for details.
        !           by default, assume spherical sample, i.e 2/3 of full screeni
        do ia = 1, 3  
           do ib = 1, 3  
              chi_g(ig, ia, ib) = chi(ia, ib) * pi4vcell * g0mask(ia, ib)  
           end do
        end do
        !            write(9,*) 'G=0 component:',chi_g(ig,:,:)
     end if
  end do
  !
  !
  !
  deallocate(r_g)
  
  return  

end subroutine wrapup_chi
!
!     ------- local subroutines ------------------------
!
subroutine myvecprod(a, b, c)

  use constants
  implicit none  
  !
  !     computes a = b x c
  !
  real(dp), intent(in) :: b(3), c(3)
  real(dp), intent(out) :: a(3)

  a(1) = b(2) * c(3) - b(3) * c(2)  
  a(2) = b(3) * c(1) - b(1) * c(3)  
  a(3) = b(1) * c(2) - b(2) * c(1)  

end subroutine myvecprod
