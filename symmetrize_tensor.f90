!*
subroutine symmetrize_tensor(tens, crys, syms)  
  !
  !     symmetrizes the tensor tens
  !
  use flibcalls_module  
  include 'use.h'  
  implicit none                  ! implicit? Just say no!
  include 'interface.h'  
  !
  !     INPUT:
  !     -----
  !
  type(crystal), intent(in) :: crys  
  type(symmetry), intent(in) :: syms  
  !
  !     INPUT/ OUTPUT:
  !     -------------
  !
  real(dp), intent(inout) :: tens(3, 3)  
  !
  !     -------- local variables ---------
  !
  integer :: nops, n  
  real(dp) :: symtens(3, 3), d, matinvert,alf(3, 3), alfm1(3, 3)  

  symtens = dzero
  nops = 0  
  do n = 1, syms%ntrans  
     alf = syms%rsymmat(:, :, n)  
     alfm1 = alf  
     d=matinvert(alfm1)
     if (abs(d) > dzero) then  
        nops = nops + 1  
        symtens = symtens + matmul(alf, matmul(tens, alfm1))  
     else
        write(9, *) 'warning: symmetry operation', n, 'is singular'  
     end if

  end do

  tens = symtens / real(nops, dp)  

end subroutine symmetrize_tensor
