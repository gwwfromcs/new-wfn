!
subroutine symgen(ipr, a, coorat, natom, ntype, invers_no, &
    nops,  ntrans, gmtrx, tnp, mxdatm, rmtrx, trans,ind_rot)
  !
  use constants
  use flibcalls_module  
  implicit none  
  !
  integer, intent(in) :: &
       ipr, &        ! ipr=1 print
       mxdatm        ! array dimensioning
  !
  !     OUTPUT:
  !     ------
  !
  integer, intent(out) :: ntrans, gmtrx(48, 3, 3),ind_rot(48,5000)  
  real(dp), intent(out) :: tnp(48, 3), trans(3, 3, 48)
  integer, intent(out) :: &
       rmtrx(48, 3, 3)    ! rotation matrices in cartesian coordinates
  !
  !     This is the driver routine to generate the symmetry
  !     operations.
  !
  !     Adapted by J.L. Martins from the program GROUP
  !     written in 1974 by Warren and Worlton
  !     Computer Physics Communications, vol 8, 71-74 (1974)
  !     incorporated by Eckold et al into UNISOFT.
  !
  !     Modified by Alberto Garcia (1990)
  !
  !
  !     Let's see if we can get the rotations both in real and
  !     in reciprocal space:
  !
  !     gmtrx : g-space representation
  !     rmtrx: r-space representation
  !
  !     Input:
  !
  !     a(i,j) is the i-th cartesian component of the j-th primitive
  !     translation vector of the direct lattice and thus it is the
  !     transpose of the matrix A defined by Jones in 'The Theory
  !     of Brillouin Zones and Electronic States in Crystals'
  !
  !     coorat(i,j,k) is the k-th component (lattice coordinates) of
  !     the position of the j-th atom of type i.
  !
  !     natom(i) is the number of atoms of type i.
  !
  !     ntype is the total number of types of atoms.
  !
  !
  !     nops indicates whether or not to count nonzero fractional
  !     translations
  !
  !     internal:
  !
  !     b contains the reciprocal lattice vectors in the
  !     crystallographic usage, that is, WITHOUT the 2pi factor.
  !     This matrix IS Jones'.
  !
  !     na is the number of atoms.
  !
  !     ity(i) is an integer distinguishing atoms of different type
  !     i.e. different atomic species.
  !
  !     x(j,i) is the j-th cartesian component of the position vector for
  !     the i-th atom in the unit cell.
  !
  !     output:
  !
  !     ntrans is the number of point group operations.
  !
  !     gmtrx is the matrix of rotations in g-lattice coordinates.
  !     rmtrx is the matrix of rotations in r-lattice coordinates.
  !
  !     tnp(oper,ilat) is the fractional translation in latt. coor.
  !
  !     invers_no is the operation number of the inversion (0 if not
  !     present). It is used to restore the inversion symmetry by
  !     a simple change of origin when possible
  !
  !     .. Scalar Arguments ..
  integer, intent(in) :: ntype, nops  
  integer, intent(out) :: invers_no  
  !     ..
  !     .. Array Arguments ..
  real(dp), intent(in) :: a(3, 3), coorat(ntype, mxdatm, 3)  
  integer, intent(in) :: natom(ntype)  
  !     ..
  !     .. Local Scalars ..
  real(dp) :: xdum,matinvert  
  integer :: i, ihg, ind, ipm, j, k, l, li, m, na, nat, ierr,&
     include_fractional
  !     ..
  !     .. Local Arrays ..
  real(dp) :: b(3, 3), r(49, 3, 3), r1(3, 3), rlat(48, 3, 3)  
  real(dp), allocatable :: x(:,:)  
  integer, allocatable :: ity(:)  
  integer :: ib(48)  
  character(len=10) :: id(48)
  !     ..
  !     .. External Subroutines ..
  external atftmt, pgl, symm_ident, symchk  
  !     ..
  ipm = 0  
  include_fractional = nops

  allocate(x(3, ntype * mxdatm))  
  allocate(ity(ntype * mxdatm))  
  !
  !     Calculate cartesian coordinates, atom types, and number of atoms.
  !
  na = 0  
  do i = 1, ntype  
     nat = natom(i)  
     do j = 1, nat  
        ind = na + j  
        do k = 1, 3  
           x(k, ind) = dzero
           do l = 1, 3  
              x(k, ind) = x(k, ind) + a(k, l) * coorat(i, j, l)  
           end do
        end do
        ity(ind) = i  
     end do
     na = na + natom(i)  
  end do
  !
  !     Determine reciprocal lattice basis vectors.
  !     We know that
  !
  !               T                             T -1         -1
  !              B A  = 2pi (1), so   B  = 2pi(A )   = 2pi(a )
  !
  !     and we can use the linpack (SCILIB version) routines
  !     sgefa and sgedi to
  !     compute the determinant (celvol) and the inverse.
  !     But note that the 2pi factor is NOT needed here.
  !
  do i = 1, 3  
     do j = 1, 3  
        b(i, j) = a(i, j)  
     end do
  end do
  !
   xdum=matinvert(b)
  !
  call pgl(a, b, r, ntrans, ib, ihg)  
  !
  !     Subroutine pgl determines the point group of the lattice and the
  !     crystal system. The array ib contains the locations of the group
  !     operations and ntrans is the order of the group.
  !
  call atftmt(ipr, a, b, x, r, tnp, trans, ity, na, ib, ihg, &
      include_fractional, li, ntrans, invers_no, ntype, mxdatm,ind_rot)
  !
  !     Subroutine atftmt determines the point group of the crystal,
  !     the atom transformation table f0, the fractional translations
  !     tnp associated with each rotation and the multiplication
  !     table mt for the point group of the crystal. The array ib now
  !     contains operations in the point group of the crystal and ntrans
  !     is the order of this group.
  !
  !     if(li.gt.0) write(9,180)
  ! 180 format (5x,'the point group of the crystal contains the',
  !    1' inversion, therefore,',/,5x,'time reversal invariance will',
  !    2' be invoked for all wave vectors.')
  !
  !     We have the rotations in cartesian coordinates.
  !     Transform into lattice coordinates (r-space and g-space)
  !
  do l = 1, ntrans  
     !
     !           In terms of the real-space basis vectors:
     !                      T
     !              y' = ( B R a ) y     ( y are the r-lattice coord.)
     !
     !                              T
     !              watch out: b = B
     !
     !           Trans * a ...
     !
     do j = 1, 3  
        do k = 1, 3  
           r1(j, k) = dzero
           do m = 1, 3  
              r1(j, k) = r1(j, k) + trans(m, j, l) * a(m, k)  
           end do
        end do
     end do
     !
     !           B * Trans * a
     !
     do j = 1, 3  
        do k = 1, 3  
           rlat(l, j, k) = dzero
           do m = 1, 3  
              rlat(l, j, k) = rlat(l, j, k) + b(j, m) * r1(m, k)  
           end do
           rmtrx(l, j, k) = nint(rlat(l, j, k))  
        end do
     end do
  end do
  !
  !        Identify the symmetry operations
  !
  call symm_ident(ntrans, rmtrx, tnp, id)  
  !
  !
  do l = 1, ntrans  
     !
     !           In terms of the g-space basis vectors:
     !                      T  T
     !              z' = ( a  R  b ) z    ( z are the g-lattice coord.)
     !
     do k = 1, 3  
        gmtrx(l, :, k) = rmtrx(l, k, :)  
     end do
  end do
  !
  !       write the matrices and fractional translations
  !
  if (ipr == 1) then  
     write(9, 9000)  
9000 format(///4x,'Rotation matrices (r-lattice) and fractional', &
          ' translations (r-lattice)',//)
     !
     do i = 1, ntrans  
        write(9, 9010) i, ((rmtrx(i, j, k), k = 1, 3), j = 1, 3), &
             (tnp(i, k), k = 1, 3), id(i)
        !
     end do
9010 format(i5,3(3x,3i3),4x,3f9.5,1x,a5)  
     !
     write(9, *) ' skipping symmetry check!'  
  end if
  !      call symchk(ipr,ierr,ntrans,rmtrx,tnp,1)
  !      call symchk(ipr,ierr,ntrans,gmtrx,tnp,0)
  !
  deallocate(x)  
  deallocate(ity)
  return  
  !
end subroutine symgen
