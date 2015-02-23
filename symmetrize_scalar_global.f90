!
subroutine symmetrize_scalar_global(scal, gs)  
  !
  !
  use all_to_all_module  
  include 'use.h'  
  implicit none              ! implicit? Just say no!
  include 'interface.h'  
  include 'all_to_all.h' 
  ! 
  !     INPUT:
  !     -----
  !
  type(parallel_gspace), intent(in) :: &
       gs                          ! the gspace corresponding to the data
  !
  !     INPUT/OUTPUT:
  !     ------------
  !
  complex(dp), intent(inout) :: &
       scal(gs%length)             ! the quantity to be symmetrized
  !
  !     DESCRIPTION:
  !     -----------
  !
  !     Symmetrizes a scalar field across multiple processors.
  !
  !     1996 Bernd Pfrommer
  !     1997 changes by Francesco Mauri
  !     ---------------------- local variables ---------------------------
  !
  complex(dp), allocatable :: symscal(:)  
  integer :: i, ii  

  if (.not. gs%istar) return             ! Francesco Mauri's shortcut
  allocate(symscal(gs%nstar))  
  symscal = zzero
  do i = 1, gs%length  
     ii = gs%inds(i)  
     symscal(ii) = symscal(ii) + scal(i) * gs%phase(i)  
  end do

  symscal(1:gs%nstar) = symscal(1:gs%nstar) / real(gs%mstar(1:gs%nstar))

  call all_sum_all(symscal, gs%nstar)  
  do i = 1, gs%length  
     scal(i) = symscal(gs%inds(i)) * conjg(gs%phase(i))  
  end do

  deallocate(symscal)

end subroutine symmetrize_scalar_global
