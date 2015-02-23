!
subroutine symchk(ipr, ierr, ntrans, mtrx, tnp, flag)  
  !
  !     symchk checks if the symmetry operations defined by
  !     mtrx and tnp really forms a group.
  !
  !     Jan 7, 1990: AG/ 2pi factor removed ! /
  !
  use constants
  implicit none
  !
  !
  ! Out  ierr                returns 0 if no error, otherwise it returns
  !                          the number of the suspected operation.
  !
  ! Inp  flag                0 for g-space check, 1 for r-space check
  !
  !     .. Scalar Arguments ..
  integer, intent(in) :: ntrans, flag, ipr
  integer, intent(out) :: ierr
  !     ..
  !     .. Array Arguments ..
  real(dp), intent(in) :: tnp(48, 3)  
  integer, intent(in) :: mtrx(48, 3, 3)  
  !     ..
  !     .. Arrays in Common ..
  integer :: mult(48, 48), nerr(48)  
  !     ..
  !     .. Local Scalars ..
  real(dp) :: ttest, xdum  
  integer :: i, im, itest, j, k, l, m, maxerr  
  character(len=80) :: not_a_grp  
  !     ..
  !     .. Local Arrays ..
  integer :: mtest(3, 3)  
  !     ..
  !     .. External Subroutines ..
  external fatal  
  !     ..
  !     .. Common blocks ..
  common mult, nerr  
  !     ..
  !
  not_a_grp = 'The symmetry operations'// &
       ' do not form a group. Check operation no%i4$\n'
  !
  !      check for duplicate operations
  !
  if (ntrans == 1) return  
  do i = 2, ntrans  
     im = i - 1  
     do j = 1, im  
        do k = 1, 3  
           if (abs(tnp(i, k) - tnp(j, k)) > 1.0d-8) goto 30  
           do l = 1, 3  
              if (mtrx(i, k, l) /= mtrx(j, k, l)) goto 30  
           end do
        end do
        if (ipr == 1) then  
           write(9, 9000) j, i  
        end if
9000    format(/' symmetry operations',i3,' and',i3,' are equal')  
        ierr = i  
        call fatal('symchk', not_a_grp, xdum, ierr)  
        !
        return  
        !
30   end do
  end do
  !
  !      construct muliplication table
  !
  do i = 1, ntrans  
     nerr(i) = 0  
     do j = 1, ntrans  
        mult(i, j) = 0  
        !
        !      mulitiply i and j
        !
        do k = 1, 3  
           do l = 1, 3  
              mtest(k, l) = 0  
              do m = 1, 3  
                 mtest(k, l) = mtest(k, l) + mtrx(i, k, m) * mtrx(j, m, l)
              end do
           end do
        end do
        !
        !      check for match
        !
        do k = 1, ntrans  
           do l = 1, 3  
              do m = 1, 3  
                 if (mtest(l, m) /= mtrx(k, l, m)) goto 100  
              end do
           end do
           mult(i, j) = k  
100     end do
     end do
  end do
  !
  !      if translations not correct set mult(i,j) to -1
  !
  if (flag == 0) then  
     !
     if (ipr == 1) then  
        write(9, *) 'Checking in g-space'  
     end if
     do i = 1, ntrans  
        do j = 1, ntrans  
           k = mult(i, j)  
           if (k == 0) cycle
           do l = 1, 3  
              ttest = tnp(j, l)  
              do m = 1, 3  
                 ttest = ttest + real(mtrx(i, m, l), dp) * &
                      (tnp(i, m) - tnp(k, m))
              end do
              ttest = abs(ttest)  
              itest = nint(ttest)  
              if (abs(ttest - real(itest, dp)) < 1.0d-4) cycle
              !     if (abs(ttest-dble(itest)) .lt. 1.d-6) cycle
              write(9, *) i, j, l, abs(ttest - real(itest, dp))  
              mult(i, j) = -1  
              !
              goto 150  
              !
           end do
150     end do
     end do
     !
  else  
     !
     if (ipr == 1) then  
        write(9, *) 'Checking in r-space'  
     end if
     do i = 1, ntrans  
        do j = 1, ntrans  
           k = mult(i, j)  
           if (k == 0) cycle
           do l = 1, 3  
              ttest = tnp(i, l) - tnp(k, l)  
              do m = 1, 3  
                 ttest = ttest + mtrx(i, l, m) * tnp(j, m)  
              end do
              itest = nint(ttest)  
              if (abs(ttest - real(itest, dp)) < 1.0d-4) cycle
              !               if (abs(ttest-dble(itest)) .lt. 1.d-6) cycle
              write(9, *) i, j, k, l, abs(ttest - real(itest, dp))  
              mult(i, j) = -1  
              !
              goto 550  
              !
           end do
550     end do
     end do
     !
  end if
  !
  !      check multiplication table
  !
  do i = 1, ntrans  
     do j = 1, ntrans  
        if (mult(i, j) > 0) cycle
        nerr(i) = nerr(i) + 1  
        nerr(j) = nerr(j) + 1  
     end do
  end do
  !
  !      find element with max error
  !
  ierr = 0  
  maxerr = 0  
  do i = 1, ntrans  
     if (nerr(i) <= maxerr) cycle
     maxerr = nerr(i)  
     ierr = i  
  end do
  if (ierr == 0) return  
  if (ipr == 1) then  
     write(9, 9010)  
9010 format('1Multiplication table',/)  
     do i = 1, ntrans  
        write(9, 9020) (mult(i, j), j = 1, ntrans)  
9020    format(1x,48i2)  
     end do
  end if
  call fatal('symchk', not_a_grp, xdum, ierr)  
  !
  return  
  !
end subroutine symchk
