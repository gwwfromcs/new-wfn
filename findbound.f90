!
!=======================================================================
!
function scp(v1, v2)

  use constants
  implicit none

  real(dp) :: scp                    !p    scalar product without metric
  real(dp), intent(in) :: v1(3), v2(3)  

  scp = v1(1) * v2(1) + v1(2) * v2(2) + v1(3) * v2(3)
  
  return  

end function scp
!
!=======================================================================
!

subroutine findbound(qmax, b, maxi, maxj, maxk)  
  !p
  !p    adapted from subroutine ban 3.2.93, Bernd Pfrommer
  !p
  !p    ! bi,bj,bk should form a right-handed system !
  !p
  !p    subroutine findbound finds boundaries of the space built up by
  !p    the vectors bi,bj,bk if a sphere of
  !p    radius qmax should be placed inside that space. In maxi, maxj, max
  !p    integers values are returned representing the coordinates along bi
  !p    that are necessary to involve the qmax-sphere entirely.
  !p
  !     -------------------- documentation of variables ------------------
  !p    qmax      radius of involved sphere
  !p    b(i,k)    ith component of recip. lattice vector k.
  !p    maxi,maxj,maxk integer variables used to return required minimum
  !p              boundaries.
  !p    volp      volume produkt bi*(bj x bk)
  !p    proi,proj,prok  see documentation in source
  !p    bj_times_bk, bk_times_bj, bi_times_bk  absolute value of cross product
  !p

  use constants
  implicit none  
  real(dp), intent(in) :: qmax                ! radius of involved sphere
  integer, intent(out) :: maxi, maxj, maxk    ! return values: boundaries
  real(dp) :: scp                             ! standard scalar product
  real(dp), intent(in) :: b(3, 3)             ! reciprocal lattice vectors
  !
  !p    ------- local variables -------------------
  !
  real(dp) :: bi(3), bj(3), bk(3)      ! rec lattice vectors (or any vector)
  integer :: i  
  real(dp) :: volp                     ! volume product
  real(dp) :: proi, proj, prok         ! see below

  real(dp) :: bj_times_bk, bk_times_bi, bi_times_bj

! abs. value of cross prod.  
  do i = 1, 3  
     bi(i) = b(i, 1)  
     bj(i) = b(i, 2)  
     bk(i) = b(i, 3)  
  end do

  !p    volume product bi*(bj x bk)

  volp = bi(1) * (bj(2) * bk(3) - bj(3) * bk(2)) + &
       bi(2) * (bj(3) * bk(1) - bj(1) * bk(3)) + &
       bi(3) * (bj(1) * bk(2) - bj(2) * bk(1))
  !
  !p    project bi on normal vector to the plane (bk,bj) to determine
  !p    the minimum component a vector has to have to be surely outside
  !p    the sphere regardless of the component along bk and bj.
  !p    proi=bi*(bj x bk)/|bj x bk|
  !
  bj_times_bk = sqrt(scp(bj, bj) * scp(bk, bk) - scp(bj, bk)**2)
  if (bj_times_bk <= 1.0d-6) then  
     write (9, *) 'error in findbound: bj,bk parallel'  
     call mystop  
  else  
     proi = abs(volp / bj_times_bk)  
     maxi = int(qmax / proi)  
  end if
  !
  !p    project bj on normal vector to the plane (bk,bi) to determine
  !p    the minimum component a vector has to have to be surely outside
  !p    the sphere regardless of the component along bk and bi.
  !p    proj=bj*(bk x bi)/|bk x bi|
  !
  bk_times_bi = sqrt(scp(bk, bk) * scp(bi, bi) - scp(bk, bi)**2)
  if (bk_times_bi <= 1.0d-6) then  
     write(9, *) 'error in findbound: bk,bi parallel'  
     call mystop  
  else  
     proj = abs(volp / bk_times_bi)  
     maxj = int(qmax / proj)  
  end if
  !
  !p    project bk on normal vector to the plane (bi,bj) to determine
  !p    the minimum component a vector has to have to be surely outside
  !p    the sphere regardless of the component along bi and bj.
  !p    prok=bk*(bi x bj)/|bi x bj|
  !
  bi_times_bj = sqrt(scp(bi, bi) * scp(bj, bj) - scp(bi, bj)**2)
  if (bi_times_bj <= 1.0d-6) then  
     write(9, *) 'error in findbound: bi,bj parallel'  
     call mystop  
  else  
     prok = abs(volp / bi_times_bj)  
     maxk = int(qmax / prok)  
  end if

end subroutine findbound
