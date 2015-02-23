!
subroutine symmetrize_crystal(crys, syms)

  include 'use.h'  
  implicit none  
  include 'interface.h'  
  !
  !     INPUT
  !     -----
  !
  type(symmetry), intent(in) :: syms  
  !
  !     INPUT/OUTPUT
  !     ------------
  !
  type(crystal), intent(inout) :: crys  
  !
  !     DESCRIPTION:
  !     -----------
  !
  !     symmetrizes the atomic coordinates and the lattice vectors (not
  !     yet) according to the symmetry  operations.
  !
  !     1998 Bernd Pfrommer Bugs Inc.
  !
  !
  !     ------------------------------------------------------------------
  !
  integer :: i, j, k, l, m, kk, idif  
  real(dp) :: cro(3), xdif, drc(3)  
  real(dp), parameter :: eps = 1.0d-6
  real(dp), allocatable :: rtmp(:,:,:)  

  allocate(rtmp(3, crys%mxdatm, crys%ntype))  
  do j = 1, crys%ntype  
     do k = 1, crys%natom(j)  
        rtmp(:, k, j) = crys%rat(:, k, j) / pi2  
     end do
  end do

  crys%rat = dzero
  do j = 1, crys%ntype  
     do k = 1, crys%natom(j)  
        do i = 1, syms%ntrans  
           !
           !                     -1
           !            find mtrx    * (rat - tnp)
           !
           do l = 1, 3  
              cro(l) = dzero
              do m = 1, 3  
                 cro(l) = cro(l) + real(syms%mtrx(m, l, i), dp) * &
                      (rtmp(m, k, j) - syms%tnp(m, i) / pi2)
              end do
           end do
           !
           !     look for matching atom
           !
           do l = 1, crys%natom(j)  
              do m = 1, 3  
                 xdif = abs(cro(m) - rtmp(m, l, j))  
                 idif = xdif + eps  
                 if (abs(xdif - real(idif, dp)) > eps) goto 16  
              end do
              kk = l  
              goto 17  
16            continue  
           end do
           write(9, *) 'subroutine symmetrize_coordinates'  
           write(9, *) 'unable to find equivalent atom for k=', k  
           write(9, *) 'please symmetrize coordinates and retry.'  
           call mystop  
17         continue  
           !
           !              add up the transformed coordinates
           !
           drc = cro(:) - rtmp(:, kk, j)  
           drc(1) = drc(1) - real(nint(drc(1)), dp)  
           drc(2) = drc(2) - real(nint(drc(2)), dp)  
           drc(3) = drc(3) - real(nint(drc(3)), dp)  
           cro = rtmp(:, kk, j) + drc  
           crys%rat(:, kk, j) = crys%rat(:, kk, j) + cro(:)  
        end do
     end do
  end do

  crys%rat = crys%rat / real(syms%ntrans, dp) * pi2  
  
  deallocate(rtmp)  
  !
  !     symmetrize lattice vectors
  !
end subroutine symmetrize_crystal
