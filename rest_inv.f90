!
logical function rest_inv(ipr, n, ntype, natom, coorat, nameat, &
     tnp, mxdatm)
  !
  use constants
  implicit none  
  !
  !     Restores the inversion symmetry if possible.
  !     Alberto Garcia, 1990-1991
  !
  !     n:   Operation number associated with the inversion
  !
  integer, intent(in) :: n, ntype, mxdatm  
  integer, intent(in) :: natom(ntype)  
  real(dp), intent(inout) :: coorat(ntype, mxdatm, 3)  
  character(len=2), intent(in) :: nameat(ntype)
  real(dp), intent(in) :: tnp(48, 3)
  integer, intent(in) :: ipr
  !
  !
  integer :: i, j, k
  !
  !
  rest_inv = .false.  
  !
  !
  !
  if (n <= 0) return  
  !
  !     Check whether the inversion has a non-symmorphic translation
  !
  if (max(abs(tnp(n, 1)), abs(tnp(n, 2)), abs(tnp(n, 3))) > 1.d-6) then
     !
     rest_inv = .true.  
     !
     if (ipr == 1) then  
        write(9, '(/,a,//,a,//)') ' Restoring inversion symmetry...', &
             ' New coordinates:'
     end if
     !
     !
     do i = 1, ntype  
        do j = 1, natom(i)  
           do k = 1, 3  
              coorat(i, j, k) = coorat(i, j, k) - dhalf * tnp(n, k)  
           end do
           if (ipr == 1) then  
              write(9, 9030) nameat(i), j, (coorat(i, j, k), k = 1, 3)
           end if
9030       format(1x,a2,' #',i4,':',3f18.12)  
        end do
     end do
     !
  end if
  !
  if (ipr == 1) then  
     write(9, '(/)')  
  end if
  !
  return  
  !
end function rest_inv
