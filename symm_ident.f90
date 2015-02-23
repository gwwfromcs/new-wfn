!
subroutine symm_ident(ntrans, mtrx, tnp, id)  
  !
  use constants
  implicit none  
  !
  integer, intent(in) :: ntrans  
  integer, intent(in) :: mtrx(48, 3, 3)  
  real(dp), intent(in) :: tnp(48, 3)  
  character(len=10), intent(out) :: id(48)
  !
  integer :: i, j, oper, ierr, ipt  
  logical :: proper, nonsym(48)  
  real(dp) :: det, trace, te, tt  
  integer :: iv1(3)  
  real(dp) :: a (3, 3), z (3, 3), wr (3), wi (3), fv1 (3)  
  real(dp) :: t (3), e (3), canon_t (3)  
  real(dp) :: rdot  
  character(len=2) :: axes(-2:2)  
  !
  data axes / 'C2', 'C3', 'C4', 'C6', 'E ' /  

  do oper = 1, ntrans  
     !
     !       Copy the matrix to a double precision format for
     !       processing with EISPACK. Also, copy the non-primitive
     !       translation.
     !
     do i = 1, 3  
        do j = 1, 3  
           a(i, j) = mtrx(oper, i, j)  
           t(j) = tnp(48, j)  
        end do
     end do
     !
     !       Compute determinant  and trace
     !
     det = a(1, 1) * a(2, 2) * a(3, 3) + a(2, 1) * a(3, 2) * a(1, 3) + &
          a(1, 2) * a(2, 3) * a(3, 1) - a(3, 1) * a(2, 2) * a(1, 3) - &
          a(2, 1) * a(1, 2) * a(3, 3) - a(1, 1) * a(3, 2) * a(2, 3)
     !
     proper = (nint(det) == 1)  
     !
     trace = a(1, 1) + a(2, 2) + a(3, 3)  
     !
     if (proper) then  
        !
        !           Proper operation
        !
        id(oper) = axes(nint(trace - 1))  
        !
     else  
        !
        !           R = IS , where S is proper
        !
        call sscal(9, dmone, a, 1)  
        id(oper) = 'I'//axes(nint(-trace - 1))  
        !
     end if
  end do
  !
  return  
  !
end subroutine symm_ident



