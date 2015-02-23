!     @process extchk
!
subroutine induced_current(curr, chi, bfield, gs, crys) 
  ! 
  include 'use.h'  
  implicit none              ! implicit? Just say no!
  include 'interface.h'  
  !
  !     1996 Bernd Pfrommer
  !
  !     INPUT
  !     -----
  !

  type(parallel_gspace), intent(in) :: &
       gs                               ! gspace on which the data is defined
  type(crystal), intent(in) :: crys    
  real(dp), intent(in) :: &
       bfield(3)                  ! the unnormalized direction of the B-field
  complex(dp), intent(in) :: &
       chi(gs%length, 3, 3)       ! susceptibility tensor field
  !
  !     OUTPUT:
  !     ------
  !
  complex(dp), intent(out) :: &
       curr(gs%length, 3)      ! induced current as computed from chi, bfield
  !
  !     DESCRIPTION:
  !     -----------
  !
  !     Given chi and the direction of the B field, this subroutine
  !     computes the induced current in for a B-field of 1 Gauss.
  !     The units are gauss/seconds = esu/(cm^2 s)
  !
  !      J = c/4pi * nabla x B_induced
  !
  !
  !     Only the imaginary part of J should be used:
  !
  !     J_true(r) = Imag(FFT(J_returned))(r)
  !
  !
  !     ------------- local variables
  !
  real(dp) :: bflen, gv(3)  
  complex(dp) :: em(3), jm(3), bf(3)  
  integer :: ig, i

  bflen = dot_product(bfield, bfield)  
  if (bflen > 1.0d-10) bf = bfield / sqrt(bflen)  
  !
  !     convert to gauss/seconds.
  !
  bf = bf / (137.036**2 * pi4) * 2.99792458d10  
  do ig = 1, gs%length  
     do i = 1, 3  
        em(i) = chi(ig, i, 1) * bf(1) + chi(ig, i, 2) * bf(2) + &
             chi(ig, i, 3) * bf(3)
     end do
     gv = matmul(crys%bvec, real(gs%gvec(:, ig), dp))  
     call mycmplxvecprod(jm(1), gv(1), em(1))           ! jm = gv x em
     !
     !        and stick it into final array
     !
     curr(ig, :) = jm(:)  
  end do

  return  

end subroutine induced_current
!
!     --------- computes complex vector cross-product ------------
!
subroutine mycmplxvecprod(a, b, c)  
  !
  use constants
  implicit none  
  !
  !     computes a = b x c
  !
  !     where b is real, c is complex, and a is complex
  !
  complex(dp), intent(in) :: c(3)  
  complex(dp), intent(out) :: a(3)
  real(dp), intent(in) :: b(3)

  a(1) = b(2) * c(3) - b(3) * c(2)  
  a(2) = b(3) * c(1) - b(1) * c(3)  
  a(3) = b(1) * c(2) - b(2) * c(1)
  
end subroutine mycmplxvecprod

