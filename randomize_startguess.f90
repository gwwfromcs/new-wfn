!
subroutine randomize_startguess(crys, gs, gsub, random_startvec, neig, xvec)
  !
  include 'use.h'  
  use all_to_all_module  
  implicit none              ! implicit? Just say no!
  include 'interface.h'  
  include 'all_to_all.h' 
  include 'flibcalls.ph'  
  !
  !
  !     INPUT:
  !     -----
  !
  type(crystal), intent(in) :: crys  
  type(parallel_gspace), intent(in) :: gs  
  real(dp), intent(in) :: &
       random_startvec, &  ! percent of random components
       gsub                ! sub gspace gmax
  integer, intent(in) :: &
       neig                ! number of eigenvectors
  !
  !     INPUT/OUTPUT:
  !     ------------
  !
  complex(dp), intent(inout) :: &
       xvec(gs%length * neig)    ! the wave functions without/with random
  !
  !     1997 Michel Cote, Bernd Pfrommer
  !
  !
  !     DESCRIPTION:
  !     -----------
  !
  !     Throws in random components into the starting guess. This is
  !     only done for those gvectors that are outsided the starting
  !     guess energy cutoff, and would be zero otherwise.
  !
  !
  !
  !     ---------------------- local variables ------------------------
  !
  real(dp) :: gm, gm2, a, a2, temp_norm, temp, temp2, prefac, dsqtemp, gsub2
  complex(dp) :: random  
  real,allocatable:: tempt(:)
  integer :: l, j  

  gm = gs%gmax
  gm2 = gm * gm
  a = done / done       ! one over the decay length in a.u.
  a2 = a * a
  gsub2 = gsub * gsub

  temp_norm = a2 * crys%vcell / (dtwo * pi**2) * (-gm / (dsix * &
       (a2 + gm2)**3) + gsub / (dsix * (a2 + gsub2)**3) + &
       gm / (24.0d0 * a2 * (a2 + gm2)**2) - gsub / (24.0d0 * a2 * &
       (a2 + gsub2)**2) + gm / (16.0d0 * a2 * a2 * (a2 + gm2)) - &
       gsub / (16.0d0 * a2 * a2 * (a2 + gsub2)) + (atan(gm / a) - &
       atan(gsub / a)) / (16.0d0 * a * a2 * a2))

  prefac = sqrt(dsix * random_startvec / (1.0d2 * temp_norm))  
  allocate (tempt(2*gs%length))  !much quicker on IBM to allocate problems  
                                 ! when too large though
 
  do l = 1, neig  
     call random_number(tempt)          !and calculate a single big array
     do j = 1, gs%length  
        random = tempt(2*j-1) * prefac * a / ((a2 + gs%ekin(j))**2) * &
             exp(cmplx(dzero, pi2, dp) * tempt(2*j))
        if (gs%ekin(j) > gsub2) then  
           if (.not. crys%icomplex) then  
              xvec((l - 1) * gs%length + j) = xvec((l - 1) * gs%length + j) + &
                   real(random, dp)
           else  
              xvec((l - 1) * gs%length + j) = xvec((l - 1) * gs%length + j) + &
                   random
           end if
        end if
     end do
     temp = dzero
     do j = 1, gs%length  
        temp = temp + abs(xvec((l - 1) * gs%length + j) )**2  
     end do
     call all_sum_all(temp)  
     dsqtemp = 1 / sqrt(temp)  
     do j = 1, gs%length  
        xvec((l - 1) * gs%length + j) = xvec((l - 1) * gs%length + j) * &
             dsqtemp
     end do
  end do

  return
  
end subroutine randomize_startguess
