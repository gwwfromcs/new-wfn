!*
subroutine symmetrize_tensor_j(tens, crys, syms)  
  !
  !     symmetrizes the tensor tens which depends of atom coordinates
  !
  use flibcalls_module  
  use constants
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
  real(dp), intent(inout) :: tens(3,3,crys%ntype,crys%mxdatm)  
  !
  !     -------- local variables ---------
  !
  integer :: n, j, natomj, k, l, m, kk, idif  
  real(dp) :: d, alf(3, 3), alfm1(3, 3), matinvert, cro(3), xdif
  real(dp), allocatable :: symtens(:,:,:,:)
  real(dp), parameter :: eps =1.0d-6

  allocate(symtens(3,3,crys%ntype,crys%mxdatm))
  symtens = dzero
  do n = 1, syms%ntrans  
     alf = syms%rsymmat(:, :, n)  
     alfm1 = alf 
     d=matinvert(alfm1)   
     if (abs(d) > dzero) then  
        do j=1,crys%ntype
           natomj=crys%natom(j)
           do k = 1, natomj  
              !
              !                     -1
              !            find mtrx    * (rat - tnp)
              !
              do l = 1, 3  
                 cro(l) = dzero  
                 do m = 1, 3  
                    cro(l) = cro(l) + real(syms%mtrx(m, l, n), dp) * &
                         (crys%rat(m, k, j) - syms%tnp(m, n))
                 end do
              end do
              do l = 1, natomj  
                 do m = 1, 3  
                    xdif = abs(cro(m) - crys%rat(m, l, j)) / pi2  
                    idif = int(xdif + eps)  
                    if (abs(xdif - real(idif, dp)) > eps) goto 16  
                 end do
                 kk = l  
                 goto 17  
16               continue  
              end do
              write(9, *) 'subroutine symmetrize_tensor_j'  
              write(9, *) 'unable to find equivalent atom for k=', k  
              call mystop  
17            continue  
                !
                !              rotate force and add
                !
              symtens(:,:,j,k) = symtens(:,:,j,k) + &
                   matmul(alf, matmul(tens(:,:,j,kk), alfm1))  
           enddo
        enddo
     else
        write(9, *) 'warning: symmetry operation', n, 'is singular'  
     end if
  end do

  tens = symtens / real(syms%ntrans, dp)  
  
  deallocate(symtens)

end subroutine symmetrize_tensor_j
