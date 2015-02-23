!
!     This file contains a set of routines which are used for the
!     diagonalization and the direct energy minimization
!
!     ==================================================================
!
subroutine overlap_matrix(ovlmat, lvec, rvec, m, len)  
  !
  !     FLOPS: 8*m*m*len
  !
  use constants
  use all_to_all_module  
  implicit none  
  include 'all_to_all.h'  
  include 'flibcalls.ph'
 
  complex(dp), intent(out) :: ovlmat(*)
  complex(dp), intent(in) :: lvec(*), rvec(*)  
  integer, intent(in) :: m, len  

  call mzgemm('C', 'N', m, m, len, zone, lvec(1), len, rvec(1), len, &
       zzero, ovlmat(1), m)
  call fast_all_sum_all_alloc_dc(m * m)  
  call fast_all_sum_all_complex(ovlmat(1), m * m)
  
  return  

end subroutine overlap_matrix
!
!     ==================================================================
!
subroutine cheap_overlap_matrix(ovlmat, vec, m, len)  
  !
  !     uses zherk, and assumes lvec = rvec
  !
  !
  !     FLOPS: 4*m*m*len
  !
  use all_to_all_module  
  implicit none  
  include 'all_to_all.h'  
  include 'flibcalls.ph'

  complex(dp), intent(out) :: ovlmat(*)
  complex(dp), intent(in) :: vec(*)  
  integer, intent(in) :: m, len  

  call mzherk('L', 'C', m, len, zone, vec(1), len, zzero, ovlmat(1), m)
  call fast_all_sum_all_alloc_dc(m * m)  
  call fast_all_sum_all_complex(ovlmat(1), m * m)  
  call matrix_complete(ovlmat, m)                ! mirror to get all elements

  return  

end subroutine cheap_overlap_matrix
!
!     ==================================================================
!
subroutine hermit_overlap_matrix(ovlmat, lvec, rvec, m, len)  
  !
  !     FLOPS: 4*m*m*len
  !
  use all_to_all_module  
  implicit none  
  include 'all_to_all.h'  
  include 'flibcalls.ph'

  complex(dp), intent(out) :: ovlmat(*)
  complex(dp), intent(in) :: lvec(*), rvec(*)  

  integer :: m, len  

  call mzher2k('L', 'C', m, len, zone, lvec(1), len, rvec(1), &
       len, zzero, ovlmat(1), m)
  call fast_all_sum_all_alloc_dc(m * m)  
  call fast_all_sum_all_complex(ovlmat(1), m * m)

  return  

end subroutine hermit_overlap_matrix
!
!     ==================================================================
!
subroutine printmatrix(mat, msize, nsize)

  use constants
  implicit none

  integer, intent(in) :: msize, nsize  
  integer :: i, j  
  complex(dp), intent(in) :: mat(msize, nsize)
  
  do i = 1, msize  
     !         write(9,'(100(''('',f15.12,'','',f15.12,'')''))')
     !     $        (mat(i,j),j=1,nsize)
     write(9, '(100(''('',f23.20,'','',f23.20,'')''))') (mat(i, j), &
          j = 1, nsize)
  end do

end subroutine printmatrix
!
!     ==================================================================
!
subroutine printdmatrix(mat, msize, nsize)

  use constants
  implicit none  

  integer, intent(in) :: msize, nsize  
  integer :: i, j  
  real(dp), intent(in) :: mat(msize, nsize)
  
  do i = 1, msize  
     write(9, '(100f15.6)') (mat(i, j), j = 1, nsize)  
  end do

end subroutine printdmatrix
!
!     ==================================================================
!
subroutine matrix_complete(mat, size)  
  !
  use constants
  implicit none  

  integer, intent(in) :: size
  integer :: i, j  
  complex(dp), intent(inout) :: mat(size, size)

  do i = 1, size  
     do j = 1, i - 1  
        mat(j, i) = conjg(mat(i, j))  
     end do
  end do

end subroutine matrix_complete
!
!     ==================================================================
!
subroutine matrix_polish (s, mat, w1, w2, m)  
  !
  use constants
  implicit none  
  include 'flibcalls.ph'  
  !
  !     removes second order error from matrices
  !
  !
  !     INPUT
  !
  complex(dp), intent(in) :: s(*)                        ! s-matrix see below
  integer, intent(in) :: m                                   ! size of matrix
  !
  !     INPUT/OUTPUT:
  !
  !     mat = (9*mat -3*mat*s - 3*s*mat + s^H * mat *s)/4
  !
  complex(dp), intent(inout) :: mat(*)  
  !
  !     WORK:
  !
  complex(dp) :: w1(*), w2(*)  
  !
  !     ------ local variables ---------------------
  !
  integer :: i  

  call mzgemm('N', 'N', m, m, m, zone, mat(1), m, s(1), m, zzero, w1(1), m)
  call mzgemm('C', 'N', m, m, m, zone, s(1), m, w1(1), m, zzero, w2(1), m)
  call mzgemm('N', 'N', m, m, m, zmthree, mat(1), m, s(1), m, zone, w2(1), m)
  call mzgemm('C', 'N', m, m, m, zmthree, s(1), m, mat(1), m, zone, w2(1), m)
  call mzaxpy(m * m, znine, mat(1), 1, w2(1), 1)  
  do i = 1, m * m
     mat(i) = w2(i) * dqtr  
  end do

end subroutine matrix_polish
!
!     ==================================================================
!
subroutine wavefn_rearrange(rows, cols, svec, dvec, iord, mnew)  
  !
  !     rearranges wave functions according to order array
  !
  use constants
  implicit none  
  include 'flibcalls.ph'
  
  complex(dp), intent(in) :: svec(*)
  complex(dp), intent(out) :: dvec(*)  
  integer, intent(in) :: rows, cols, iord(*), mnew  
  !
  !
  !
  integer :: i  
  do i = 1, mnew  
     call mzcopy(rows, svec(rows * (iord(i) - 1) + 1), 1, &
          dvec(rows * (i - 1) + 1), 1)
  end do

end subroutine wavefn_rearrange
!
!     ================== compute average kinetic energy ================
!
function average_ekin(x, tkinav, ekin, m, len)  

  use constants
  use all_to_all_module  
  implicit none  
  include 'all_to_all.h'  
  !
  !     INPUT:
  !     -----
  !
  integer :: m, &                          ! number of vectors to precondition
       len                                                 ! length of vectors
  complex(dp), intent(in) :: x(len, m) ! wave functions used to compute tkinav
  real(dp), intent(in) :: ekin(len)     ! kinetic energy of the gvectors |k+G|
  !
  !     OUTPUT:
  !     ------
  !
  real(dp) :: average_ekin
  !
  real(dp) :: tkinav(m)                                           ! work array
  !
  !     DESCRIPTION:
  !     computes and returns average kinetic energy
  !     ---------- local variables -----------------
  !
  integer :: i, j  
  real(dp) :: eprecond, xsq

  !     calculate maximum average kinetic energy of all trial vectors
  do i = 1, m  
     tkinav(i) = dzero  
     do j = 1, len  
        xsq = real(x(j, i) * conjg(x(j, i)), dp)  
        tkinav(i) = tkinav(i) + xsq * ekin(j)  
     end do
  end do
  call all_sum_all(tkinav, m)  
  average_ekin = maxval(tkinav)  

end function average_ekin
!
!     ==================================================================
!
!     @process extchk
!
subroutine precondition(v, eprecond, ekinmod, ekin, m, len)

  use constants
  use all_to_all_module  
  implicit none  
  include 'all_to_all.h'  
  !
  !     INPUT:
  !     -----
  !
  integer, intent(in) :: m, &              ! number of vectors to precondition
       len                                                 ! length of vectors
  real(dp), intent(in) :: eprecond, &      ! the preconditioning energy cutoff
       ekin(len), &                     ! kinetic energy of the gvectors |k+G|
       ekinmod(3)                     ! kinetic energy modification parameters
  !
  !     INPUT/OUTPUT:
  !     ------------
  !
  complex(dp), intent(inout) :: v(len, m)       ! vectors to be preconditioned
  !
  !     DESCRIPTION:
  !     preconditions vectors v(len,m) with the preconditioner
  !     of Teter, Payne, and Allen
  !
  !     ---------- local variables -----------------
  !
  integer :: i, j  
  complex(dp) :: gt  
  real(dp) :: pcfn, ooekinmod2, ooeprecond
  real(dp), external :: precfn, myderf

  ooekinmod2 = done / ekinmod(2)
  ooeprecond = done / eprecond

  do i = 1, m  
     do j = 1, len  
        pcfn = precfn((ekin(j) + ekinmod(1) * (done + myderf((ekin(j) - &
             ekinmod(1)) * ooekinmod2))) * ooeprecond)
        v(j, i) = v(j, i) * pcfn  
     end do
  end do

  return  

end subroutine precondition
!
!     ==================================================================
!
function precfn(x)  

  use constants
  implicit none  

  real(dp) :: precfn
  real(dp), intent(in) :: x  
  !
  !     teter payne allan preconditioning
  !
  !
  !      double precision x1,x2,x3,x4
  !      x2=x*x
  !      x3=x2*x
  !      x4=x3*x
  !
  !      x1=27.d0+18.d0*x+12.d0*x2+8.d0*x3
  !      precfn=x1/(x1+16.0*x4)
  !
  !      preconditioning a la f. mauri
  !
  precfn = done / max(done, x)
  
  return

end function precfn
!
!     ==================================================================
!
function trace (mat, m)  

  use constants
  implicit none  

  real(dp) :: trace
  integer, intent(in) :: m  
  complex(dp), intent(in) :: mat(m, m)  
  integer :: i 
 
  trace = dzero
  do i = 1, m  
     trace = trace + real(mat(i, i), dp)  
  end do

  return  

end function trace
!
!     ==================================================================
!
function contract(mat1, mat2, m)  

  use constants
  implicit none

  real(dp) :: contract  
  integer, intent(in) :: m  
  complex(dp), intent(in) :: mat1(m, m), mat2(m, m)  
  integer :: i, j  
  real(dp) :: sum, odsum
  
  sum = dzero
  odsum = dzero  
  do i = 1, m  
     sum = sum + real(mat1(i, i), dp) * real(mat2(i, i), dp)  
     do j = i + 1, m  
        odsum = odsum + real(mat1(j, i) * conjg(mat2(j, i)), dp)  
     end do
  end do
  contract = sum + dtwo * odsum
  
  return  

end function contract
