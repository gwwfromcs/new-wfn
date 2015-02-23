!
subroutine qualitykp(crys, syms, kp)

  include 'use.h'  
  implicit none           ! implicit? Just say no!
  include 'interface.h'  
  !
  !     INPUT:
  !     -----
  !
  type(crystal), intent(in) :: crys  
  type(symmetry), intent(in) :: syms  
  type(kpoint), intent(in) :: kp  
  !
  !     DESCRIPTION
  !     -----------
  !
  !     Computes R_min, which is a measure of the quality of the kpoint
  !     grid.
  !
  !     1998 F. Mauri, B. Pfrommer, based on a code by A. de Vita.
  !
  !
  !     ------------ local variables ------------------------------
  !
  integer :: mxr(3), &  ! maximum grid limit
       nkpts            ! total number of symmetrized kpoints
  real(dp), allocatable :: &
       symk(:,:), &     ! symmetrized kpoints
       ws(:)            ! weights
  real(dp) :: rmax, &   ! maximum radius for which to check
       sum, phase, rr, rmin, rvec(3), vec(3)
  integer :: i, ix, iy, iz, k, ks  
  real(dp), parameter :: eps = 1.0d-4  ! criterion for which to decide on rmin
  !
  !     unfold symmetrize kpoints into full grid by using the symmetry
  !     operations.
  !
  nkpts = syms%ntrans * kp%nrk  
  allocate(symk(3, nkpts))  
  allocate(ws(nkpts))  
  do i = 1, syms%ntrans  
     do k = 1, kp%nrk  
        ks = k + (i - 1) * kp%nrk  
        ws(ks) = kp%w(k)  
        symk(:, ks) = matmul(syms%mtrx(:,:, i), kp%rk(:, k))  
     end do
  end do
  rmax = 10.d0 * (dthree * crys%vcell / pi4)**dthird
  !
  !     find grid boundaries that match nicely
  !
  call findbound(rmax, crys%avec, mxr(1), mxr(2), mxr(3))  
  !
  !     now find rmin
  !
  rmin = rmax  
  do ix = -mxr(1), mxr(1)  
     do iy = -mxr(2), mxr(2)  
        do iz = -mxr(3), mxr(3)  
           sum = dzero
           do k = 1, nkpts  
              phase = pi2 * (symk(1, k) * ix + symk(2, k) * iy + &
                   symk(3, k) * iz)
              sum = sum + cos(phase) * ws(k)  
           end do
           sum = abs(sum)  
           if (sum > eps) then  
              vec(1) = real(ix, dp)  
              vec(2) = real(iy, dp)  
              vec(3) = real(iz, dp)  
              rvec = matmul(crys%avec, vec)  
              rr = sqrt(rvec(1)**2 + rvec(2)**2 + rvec(3)**2)  
              if (rr > eps**2) rmin = min(rr, rmin)  
           end if
        end do
     end do
  end do
  if (rmin < rmax) then  
     write(9, 1010) rmin  
  else  
     write(9, 1020) rmin  
  end if
  deallocate(symk)  
  deallocate(ws)

1010 format(/'    K POINT QUALITY ESTIMATE: R_MIN NONZERO =',f12.3)  

1020 format(/'    K POINT QUALITY ESTIMATE: R_MIN NONZERO >',f12.3)  

end subroutine qualitykp
