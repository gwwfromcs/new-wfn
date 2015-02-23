!
subroutine atftmt(ipr, a, ai, x, r, tnp, trans, ity, na, ib, ihg, &
    include_fractional, li, nc, invers_no, mxdtyp, mxdatm,ind_rot)
  !
  use constants
  implicit none  
  !
  !     INPUT:
  !     -----
  !
  integer, intent(in) :: &
       ipr, &       ! print flag
       ihg, &       ! holohedral group number
       na, &        ! total number of atoms (of all kinds)
       mxdtyp, &    ! max number of atom types
       mxdatm, &      ! max number of atoms of each type
       include_fractional  ! -2 means don't include frac translations

  real(dp), intent(in) :: &
       x(3, mxdtyp * mxdatm), &  ! compact list of all coordinates (cartesian)
       a(3, 3), &                ! realspace lattice vectors
       ai(3, 3)                  ! inverse of lattice vectors (no 2pi)
  !
  !     INPUT/OUTPUT:
  !
  integer, intent(inout) :: &
       ity(mxdtyp * mxdatm), &   ! compact list of all the types of all atoms
       ib(48), &                 ! index map for symmetry operations
       nc                        ! number of symm-ops without/with basis
  real(dp), intent(inout) :: &
       r(49, 3, 3)               ! the rotation matrices as they come
                                 ! out of the pgl subroutine
  !
  !     OUTPUT:
  !     ------
  !
  integer, intent(out) :: &
       li, &                     ! something to do with inversion?
       invers_no,ind_rot(48,500)                 ! which operation is the inversion
  real(dp), intent(out) :: &
       trans(3, 3, 48), &        ! the cartesian realspace transf. matrices
       tnp(48, 3)                ! the nonprimitive translations
  !
  !     DESCRIPTION:
  !     -----------
  !
  !
  !     subroutine atftmt determines the point group of the crystal, the
  !     atom transformation table,f0, the fractional translations,tnp,
  !     associated with each rotation and finally the multiplication table
  !     mt, for the point group of the crystal. ib now contains
  !     operations in the p.g. of the crystal and ntrans is the order of
  !     this group.
  !
  !
  !     1997 Bernd Pfrommer, based on a routine from Alberto Garcia's
  !     code. I cleaned up some more, put in dynamic memory allocation,
  !     explicit typing, fortran 90 style, and more comments. Also,
  !     some shortcuts were put in to speed up the case where the
  !     symmetry is low. The rdiff array precomputes the rotated
  !     vectors to speed up things
  !
  !     ------------------------ local variables -----------------------
  !
  !     ..
  real(dp) :: v(3, 48), vr(3), vt(3), xb(3)  
  integer :: ia(48), ic(48), mt(48, 48), &
       ipm                 ! print flag for multiplication table

  real(dp) :: da, dif, ts, vs, difmax  
  real(dp), parameter :: eps = 1.0d-8
  integer :: i, il, is, isy, iu, j, k, k1, k2, k3, k4, ks, l, m, n, &
       n1, n2, n3, nca, ni
  real(dp), allocatable :: rx(:,:), rdiff(:,:,:)  
  integer, allocatable :: if0(:,:)  
  character(len=12), parameter :: cst(7) = &
       (/ 'triclinic   ', 'monoclinic  ', 'orthorhombic', 'tetragonal  ', &
       'cubic       ', 'trigonal    ', 'hexagonal   ' /)

  invers_no = 0  
  ipm = 0   

  allocate(rx(3, mxdtyp * mxdatm))  
  allocate(if0(48, mxdtyp * mxdatm))  
  allocate(rdiff(3, na, na))  
  !
  !     eps should be slightly larger than computer precision
  !
  nca = 0  
  ni = 13  
  if (ihg < 6) ni = 25  
  li = 0  
  do n = 1, nc                     ! loop over all lattice symmetry operations
     l = ib(n)                     ! get index of symmetry operation
     ic(n) = ib(n)  
     !
     !        operate on all atoms with symmetry operation l, and store in
     !        list rx
     !
     do k = 1, na  
        do i = 1, 3  
           rx(i, k) = dzero  
           do j = 1, 3
              rx(i, k) = rx(i, k) + r(l, i, j) * x(j, k)  
           end do
        end do
     end do
     !
     !
     !     This piece of code is pretty naive, and scales like order(n**3).
     !     It basically checks if for each
     !
     !       R*x_1 - x_2   there is a matching R*x_3-x4
     !
     !     excluding the trivial case of x_1 .eq. x_3  and  x_2 .eq. x_4
     !
     !
     ! precompute the rdiffs first
     do k1 = 1, na  
        do k2 = 1, na  
           xb(:) = rx(:, k1) - x(:, k2)  
           ! rdiff = R*x_1 -x_2
           call rlv(ai, xb, rdiff(1, k2, k1), il)  
           !     subroutine rlv removes a direct lattice vector from xb
           !     leaving the remainder in rdiff. if a nonzero lattice vector was
           !     removed, il is made nonzero.
        end do
     end do
     !
     !
     !
     difmax = dzero
     do k1 = 1, na                       ! double loop over compact atom list
        do k2 = 1, na  
           if (ity(k1) == ity(k2)) then                    ! same type atoms?
              vr = rdiff(:, k2, k1)                !vr stands for v-reference.
              ks = 0  
              do k3 = 1, na  
                 do k4 = 1, na  
                    if (ity(k3) == ity(k4)) then  
                       vt = rdiff(:, k4, k3)       ! vt stands for v-test
                       dif = dzero
                       do i = 1, 3  
                          da = abs(vr(i) - vt(i)) + eps  
                          dif = dif + mod(da, done)  
                       end do
                       if (dif <= 10.0d0 * eps) then  
                          if0(l, k3) = k4  

                          ! if0 is the function defined in maradudin and 
                          ! vosko by eq.(2.35). it defines the atom
                          ! transformation table
                          ks = ks + k4  
!           if (ks == na * (na + 1) / 2) goto 110       ! found
                          if (ks .eq. na*(na+1)/2) then


!     Reject all symops with non-zero fractional translations.  It
!     should be helpful when you have weird fractional translations that
!     confuse adjustfft.fp.  -- David Roundy
                            if (include_fractional.eq.-1.or. &
                               vr(1)*vr(1)+vr(2)*vr(2)+  &
                               vr(3)*vr(3).eq. 0.0d0)  &
                               go to 110 ! found all
                          endif
                          goto 80  
                       end if
                    end if
                 end do
                 exit  
80               continue  
              end do
              !
              !  BP put in this shortcut (check carefully). If there
              !  is a single R*x_1- x_2 without match, then give up.
              !  this is not a symmetry operation
              if (ks == 0) goto 130  
           end if
        end do
     end do
     !
     ! this was not a symmetry operation
     goto 130  
     !
     ! we found a symmetry operation
110  continue  

     nca = nca + 1  
     !
     !        v(i,l) is the i-th cartesian component of the fractional
     !        translation associated with the rotation r(l).
     !
     v(1:3, l) = vr(1:3)  
     !
     ib(nca) = l  
     if (l == ni) then  
        li = l  
        invers_no = nca  
     end if

130 end do
  deallocate(rdiff)  
  !
  !     -------------- there are no gotos across this line ---------------
  !
  if (ipr == 1) then  
     if ((ihg == 7 .and. nca == 24) .or. (ihg == 5 .and. nca == 48)) then
        write(9, 9010) cst(ihg)  
9010    format      (/' The point group of the crystal is the full ' &
             &           ,a12,'group')
     else  
        write(9, 9000) cst(ihg), (ic(i), i = 1, nc)  
9000    format      (/' The crystal system is ',a12 &
             &           ,' with operations: ',/5x,24i3,/5x,24i3,/)
     end if
  end if
  vs = dzero  
  nc = nca  
  do n = 1, nc  
     l = ib(n)  
     vs = sum(abs(v(1:3, l)))  
  end do
  if (vs > eps) then  
     if (ipr == 1) write(9, 9030)  
9030 format   (/' The space group is non-symmorphic',/, &
          &        ' (Or a non-standard origin of coordinates is used)',/)
     isy = 0  
     is = 0  
  else  
     if (ipr == 1) write(9, 9020)  
9020 format   (/' the space group of the crystal is symmorphic',/)  
     isy = 1  
     is = 1  
  end if
  !
  !     construct the multiplication table
  !
  do n1 = 1, nc  
     do n2 = 1, nc  
        l = ib(n1)  
        m = ib(n2)  
        do i = 1, 3  
           do j = 1, 3  
              r(49, i, j) = dzero  
              do k = 1, 3  
                 r(49, i, j) = r(49, i, j) + r(l, i, k) * r(m, k, j)  
              end do
           end do
        end do
        do n3 = 1, nc  
           n = ib(n3)  
           ts = dzero  
           do i = 1, 3  
              do j = 1, 3  
                 ts = ts + abs(r(49, i, j) - r(n, i, j))  
              end do
           end do
           if (ts > 1.0d2 * eps) cycle  
           mt(l, m) = n  
           exit  
        end do
     end do
  end do
  !
  il = 1  
  iu = nc  
  if (iu > 24) iu = 24  
280 continue  
  if (ipr == 1) then  
     write(9, 9040) (ib(i), i = il, iu)  
  end if
9040 format(' Operation number  ',24i3)  

  do i = 1, na  
     do j = 1, nc  
        l = ib(j)  
        ia(j) = if0(l, i)  
        ind_rot(j,ia(j))=i
     end do
  end do

! 

!  write(9,*) " Sym OP , Atom_i , Atom_J "
!  do j = 1, nc  
!  do i = 1, na  
!     write(9,222) j,i,ind_rot(j,i)
!  end do
!  end do

222 format(3i7)

  if (nc - iu) 320, 320, 310  
310 continue  
  if (ipr == 1) then  
     write(9, 9050)  
  end if
9050 format(//)  
  il = 25  
  iu = nc  
  !
  goto 280  
  !
  !     Print multiplication table and fractional translations.
  !
320 continue  
  if (ipm == 0) goto 410  
  il = 1  
  iu = nc  
  if (nc > 24) iu = 24  
  if (is) 330, 330, 340  
330 continue  
  if (ipr == 1) then  
     write(9, 9060)  
9060 format   ('0',57x,'Multiplication table',30x, &
          &      'Fractional translations')
     write(9, 9070) (ib(i), i = il, iu)  
9070 format   ('0',4x,24i4)  
     write(9, 9080)  
9080 format   ('+',107x,'v(1)      v(2)      v(3)')  
  end if
  !
  goto 360  
  !
340 continue  
  if (ipr == 1) then  
     write(9, 9090)  
  end if
9090 format('0',57x,'Multiplication table')  
350 continue  
  if (ipr == 1) then  
     write(9, 9100) (ib(i), i = il, iu)  
  end if
9100 format('0',4x,24i4)  

360 continue  
  do j = 1, nc  
     l = ib(j)  
     do i = il, iu  
        n = ib(i)  
        ia(i) = mt(l, n)  
     end do
     if (is) 380, 380, 390  
380  continue  
     if (ipr == 1) then  
        write(9, 9120) ib(j), (ia(i), i = il, iu)  
        write(9, 9110) (v(i, l), i = 1, 3)  
9110    format      ('+',102x,3f10.4)  
     end if
     !
     goto 400  
     !
390  continue  
!     if (ipr == 1) then  
        write(9, 9120) ib(j), (ia(i), i = il, iu)  
!     end if
400 end do
9120 format(i5,24i4)  
  if (iu == nc) goto 410  
  il = 25  
  iu = nc  
  is = 1  
  !
  goto 350  
  !

410 continue  
  !
  do i = 1, nc  
     l = ib(i)  
     do j = 1, 3  
        tnp(i, j) = -v(j, l)  
        do k = 1, 3  
           trans(j, k, i) = r(l, j, k)  
        end do
     end do
  end do
  !
  deallocate(rx)  
  deallocate(if0)
  
  return  
  !
end subroutine atftmt
!
!     ------------------------------------------------------------------
!
!
subroutine pgl(a, b, r, nc, ib, ihg)  
  !
  use constants
  implicit none
  !
  !     .. Scalar Arguments ..
  integer, intent(out) :: ihg, nc  
  !     ..
  !     .. Array Arguments ..
!  real(dp), intent(in) :: a(3, 3), b(3, 3), r(49, 3, 3)  
  real(dp), intent(in) :: a(3, 3), b(3, 3)
  real(dp), intent(out) :: r(49, 3, 3)  
  integer, intent(out) :: ib(48)  
  !     ..
  !     .. Local Scalars ..
  real(dp) :: tr  
  !
  !     eps should be slightly larger than computer precision
  !
  real(dp), parameter :: eps = 1.0d-8 
  integer :: i, ihc, j, k, lx, n, nr  
  !     ..
  !     .. Local Arrays ..
  real(dp) :: vr(3), xa(3)  
  !     ..
  ihc = 0  
  !
  !     ihc is 0 for hexagonal groups and 1 for cubic groups.
  !
  nr = 24  
10 continue  
  nc = 0  
  call rot(r, nr)  

  do n = 1, nr  
     ib(n) = 0  
     tr = dzero  
     do k = 1, 3  
        do i = 1, 3  
           xa(i) = dzero  
           do j = 1, 3  
              xa(i) = xa(i) + r(n, i, j) * a(j, k)  
           end do
        end do
        call rlv(b, xa, vr, lx)  
        do i = 1, 3  
           tr = tr + abs(vr(i))  
        end do
     end do
     if (tr <= 10.0d0 * eps) then  
! PZ
! ifort has some problem.
! not sure what was the problem
! if I remove this stupid write(*,*)
! The code cannot find the correct symmetry of 2-atom silicon
! 
! it seems that the problem is resolved, r(49,3,3) should be defined as
! intent(out) 
! 
!        write(*,*)
        nc = nc + 1  
        ib(nc) = n  
     end if
  end do
  if (ihc == 0) then  
     if (nc == 12) then  
        ihg = 6  
        !
        return  
        !
     end if
     if (nc > 12) then  
        ihg = 7  
        !
        return  
        !
     end if
     if (nc < 12) then  
        nr = 48  
        ihc = 1  
        !
        goto 10  
        !
     end if
  else  
     if (nc == 16) then  
        ihg = 4  
        !
        return  
        !
     end if
     if (nc > 16) then  
        ihg = 5  
        !
        return  
        !
     end if
     if (nc < 16) then  
        if (nc == 4) then  
           ihg = 2  
           !
           return  
           !
        end if
        if (nc > 4) then  
           ihg = 3  
           !
           return  
           !
        end if
        if (nc < 4) then  
           ihg = 1  
           !
           return  
           !
        end if
     end if
  end if
  !
  !     ihg stands for holohedral group number.
  !
end subroutine pgl
!
subroutine rot(r, nr)  
  !
  use constants
  implicit none
  !
  !     .. Scalar Arguments ..
  integer, intent(in) :: nr  
  !     ..
  !     .. Array Arguments ..
  real(dp), intent(out) :: r(49, 3, 3)  
  !     ..
  !     .. Local Scalars ..
  real(dp) :: f  
  integer :: i, j, k, n, nv  
  !     ..
  do n = 1, nr  
     do i = 1, 3  
        do j = 1, 3  
           r(n, i, j) = dzero
        end do
     end do
  end do
  !
  if (nr <= 24) then  
     !
     !        define the generators for the rotation matrices
     !                                 --hexagonal group
     !
     f = drt3 * dhalf
     r(2, 1, 1) = dhalf  
     r(2, 1, 2) = -f  
     r(2, 2, 1) = f  
     r(2, 2, 2) = dhalf
     r(7, 1, 1) = dmhalf
     r(7, 1, 2) = -f  
     r(7, 2, 1) = -f  
     r(7, 2, 2) = dhalf  
     do n = 1, 6  
        r(n, 3, 3) = done
        r(n + 18, 3, 3) = done
        r(n + 6, 3, 3) = dmone
        r(n + 12, 3, 3) = dmone
     end do
     !
     !     generate the rest of the rotation matrices
     !
     do i = 1, 2  
        r(1, i, i) = done
        do j = 1, 2  
           r(6, i, j) = r(2, j, i)  
           do k = 1, 2  
              r(3, i, j) = r(3, i, j) + r(2, i, k) * r(2, k, j)  
              r(8, i, j) = r(8, i, j) + r(2, i, k) * r(7, k, j)  
              r(12, i, j) = r(12, i, j) + r(7, i, k) * r(2, k, j)
           end do
        end do
     end do
     do i = 1, 2  
        do j = 1, 2  
           r(5, i, j) = r(3, j, i)  
           do k = 1, 2  
              r(4, i, j) = r(4, i, j) + r(2, i, k) * r(3, k, j)  
              r(9, i, j) = r(9, i, j) + r(2, i, k) * r(8, k, j)  
              r(10, i, j) = r(10, i, j) + r(12, i, k) * r(3, k, j)
              r(11, i, j) = r(11, i, j) + r(12, i, k) * r(2, k, j)
           end do
        end do
     end do
     !
     do n = 1, 12  
        nv = n + 12  
        do i = 1, 2  
           do j = 1, 2  
              r(nv, i, j) = -r(n, i, j)  
           end do
        end do
     end do
  else  
     !
     !        define the generators for the rotation matrices
     !                                          --cubic group
     !
     r(9, 1, 3) = done
     r(9, 2, 1) = done
     r(9, 3, 2) = done
     r(19, 1, 1) = done
     r(19, 2, 3) = dmone
     r(19, 3, 2) = done
     do i = 1, 3  
        r(1, i, i) = done
        do j = 1, 3  
           r(20, i, j) = r(19, j, i)  
           r(5, i, j) = r(9, j, i)  
           do k = 1, 3  
              r(2, i, j) = r(2, i, j) + r(19, i, k) * r(19, k, j)
              r(16, i, j) = r(16, i, j) + r(9, i, k) * r(19, k, j)
              r(23, i, j) = r(23, i, j) + r(19, i, k) * r(9, k, j)
           end do
        end do
     end do
     do i = 1, 3  
        do j = 1, 3  
           do k = 1, 3  
              r(6, i, j) = r(6, i, j) + r(2, i, k) * r(5, k, j)  
              r(7, i, j) = r(7, i, j) + r(16, i, k) * r(23, k, j)
              r(8, i, j) = r(8, i, j) + r(5, i, k) * r(2, k, j)  
              r(10, i, j) = r(10, i, j) + r(2, i, k) * r(9, k, j)
              r(11, i, j) = r(11, i, j) + r(9, i, k) * r(2, k, j)
              r(12, i, j) = r(12, i, j) + r(23, i, k) * r(16, k, j)
              r(14, i, j) = r(14, i, j) + r(16, i, k) * r(2, k, j)
              r(15, i, j) = r(15, i, j) + r(2, i, k) * r(16, k, j)
              r(22, i, j) = r(22, i, j) + r(23, i, k) * r(2, k, j)
              r(24, i, j) = r(24, i, j) + r(2, i, k) * r(23, k, j)
           end do
        end do
     end do
     do i = 1, 3  
        do j = 1, 3  
           do k = 1, 3  
              r(3, i, j) = r(3, i, j) + r(5, i, k) * r(12, k, j)  
              r(4, i, j) = r(4, i, j) + r(5, i, k) * r(10, k, j)  
              r(13, i, j) = r(13, i, j) + r(23, i, k) * r(11, k, j)
              r(17, i, j) = r(17, i, j) + r(16, i, k) * r(12, k, j)
              r(18, i, j) = r(18, i, j) + r(16, i, k) * r(10, k, j)
              r(21, i, j) = r(21, i, j) + r(12, i, k) * r(15, k, j)
           end do
        end do
     end do
     do n = 1, 24  
        nv = n + 24  
        do i = 1, 3  
           do j = 1, 3  
              r(nv, i, j) = -r(n, i, j)  
           end do
        end do
     end do
  end if
  !
end subroutine rot
!     ----------------------------------------------

subroutine rlv(p, g, y, l)  
  !
  use constants
  implicit none  
  !
  !     INPUT:
  !     -----
  !
  real(dp), intent(in) :: g(3), &     ! vector to be multiplied
       p(3, 3)              ! multiplication matrix, e.g. lattice vectors
  !
  !
  !     OUTPUT:
  !     ------
  !  
  real(dp), intent(out) :: y(3)  ! mod(multiplied vector,lattice vector)
  integer, intent(out) :: l      ! is nonzero if a lattice vector was removed
  !
  !     subroutine rlv removes a direct lattice vector from g by
  !     operation p, leaving the remainder in y. If a nonzero lattice
  !     vector was removed, l is made nonzero.
  !
  y(1:3) = matmul(p(1:3, 1:3), g(1:3))  
  l = sum(nint(abs(y(1:3))))  
  y(1:3) = y(1:3) - done * nint(y(1:3))  

  return  
  !
end subroutine rlv
