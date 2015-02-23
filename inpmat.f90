!
subroutine inpmat (crys, syms)  
  !
  include 'use.h'  
  implicit none  
  include 'interface.h'
  !
  type(crystal), intent(in) :: crys  
  type(symmetry), intent(out) :: syms  
  !
  !     ------------------------------------------------------
  !
  integer :: lmtrx(48, 3, 3), tnp(48, 3), nops, i, k, j, ierr  
  real(dp) :: dummy(3, 3)  
  !
  write(9, 9000)  
9000 format(/' read transformation matrices',/1x, &
       &     '(nonprimitive translations) from file:'/)
  !
  open(unit = 7, file = 'SYMMETRYOPS', status = 'unknown', form = &
       'formatted')
  i = 1  
  nops = 0  
  do while (i > 0)  
     nops = nops + 1  
     read(7, *, end = 911) i, ((syms%rmtrx(j, k, nops), k = 1, 3), j = 1, 3), &
          (syms%tnp(k, nops), k = 1, 3)
     write(9, 9010) nops, ((syms%rmtrx(j, k, nops), k = 1, 3), j = 1, 3), &
          (syms%tnp(k, nops), k = 1, 3)
9010 format   (i5,3(3x,3i3),4x,3f9.5)  
  end do
911 continue  

  close(7)  

  nops = nops - 1  
  if (nops == 0) then  
     write(9, *) '*** SYMMETRYOPS has no symmetry operations!'  
     call mystop  
  end if

  syms%ntrans = nops  
  !
  !     one is just the transpose of the other
  !
  do k = 1, 3  
     do j = 1, 3  
       syms%mtrx(j, k, :)  =  syms%rmtrx(k, j, :) 
     end do
  end do
  !      do i=1, syms%ntrans
  !         write(9,9013) i, syms%rmtrx(:,:,i), syms%tnp(:,i)
  !      end do
  ! 9013 format(i5,3(3x,3i3),4x,3f8.5)
  !
  !     compute the realspace symmetry operations
  !
  do k = 1, 48  
     dummy = matmul(syms%rmtrx(:, :, k), transpose(crys%bvec))  
     syms%rsymmat(:, :, k) = matmul(crys%avec, dummy) / pi2  
  end do
  !
  !     transform to make symchk happy
  !
  do k = 1, 48  
     do i = 1, 3  
        lmtrx(k, i, :) = syms%rmtrx(i, :, k)  
     end do
     tnp(k, :) = syms%tnp(:, k)  
  end do

  syms%tnp = pi2 * syms%tnp  
  !      call symchk(0, ierr,nops,lmtrx,tnp,0)
  !
  return  
  !
end subroutine inpmat
