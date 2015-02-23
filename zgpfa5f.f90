!
subroutine zgpfa5f(z, trigs, jump, n, mm, lot, isign)

  use constants
  implicit none

  complex(dp), intent(inout) :: z(0:*)
  complex(dp), intent(in) :: trigs(0:*)
  integer, intent(in) :: jump, n
  integer, intent(in) :: mm, lot, isign

  integer :: n5, inq, jstepx, ipass
  integer :: mu, m, mh, la, error
  complex(dp), allocatable :: work(:,:)

  allocate(work(10, 0:lot - 1), stat = error)
  if (error /= 0) then
     write(*, *) 'Error allocating array `work'''
     return
  end if

  n5 = 5**mm
  inq = n / n5
  jstepx = n5 - n
  mu = mod(inq, 5)
  if (isign == -1) mu = 5 - mu

  m = mm
  mh = (m + 1) / 2

  la = 1

  !  loop on type I radix-5 passes
  !  -----------------------------
  do ipass = 1, mh

     call rad5I(z, la, n, jump, lot, inq, jstepx, isign)

     !     finished if n5 = 5
     !     ------------------
     if (n5 == 5) exit

     call rad5Itwid(z, trigs, la, n, jump, lot, inq, jstepx, isign)
     la = 5 * la

  end do

  !  loop on type II radix-5 passes
  !  ------------------------------
  do ipass = mh + 1, m

     call rad5II(z, la, n, jump, lot, inq, jstepx, isign, work)

     !     finished if last pass
     !     ---------------------
     if (ipass == mm) exit

     call rad5IItwid(z, trigs, la, n, jump, lot, inq, jstepx, isign, work)
     la = 5 * la

  end do

  deallocate(work, stat = error)
  if (error /= 0) then
     write(*, *) 'Error deallocating array `work'''
     return
  end if

  return

contains

  subroutine rad5I(z, la, n, jump, lot, inq, jstepx, isign)

    implicit none

    complex(dp), intent(inout) :: z(0:*)
    integer, intent(in) :: la, n, jump, lot
    integer, intent(in) :: inq, jstepx, isign

    integer :: ja, jb, jc, jd, je
    integer :: j, jjj, mu, nu
    integer :: jstep, jstepl
    real(dp) :: c1, c2, c3, cswap
    complex(dp) :: zja, zjb, zjc, zjd, zje
    complex(dp) :: z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11

    jstep = n / (5 * la)
    jstepl = jstep - n

    mu = mod(inq, 5)
    if (isign == -1) mu = 5 - mu

    c1 = drt5 * dqtr
    c2 = sin(0.4d0 * pi)
    c3 = sin(0.2d0 * pi)
    if (mu == 2 .or. mu == 3) then
       c1 = -c1
       cswap = c2
       c2 = c3
       c3 = cswap
    end if
    if (mu == 3 .or. mu == 4) c2 = -c2
    if (mu == 2 .or. mu == 4) c3 = -c3

    !  k = 0 (no twiddle factors)
    !  --------------------------

    do jjj = 0, n - 1, 5 * jstep
       ja = jjj

       !     "transverse" loop
       !     -----------------
       do nu = 1, inq
          jb = ja + jstepl
          if (jb < 0) jb = jb + n
          jc = jb + jstepl
          if (jc < 0) jc = jc + n
          jd = jc + jstepl
          if (jd < 0) jd = jd + n
          je = jd + jstepl
          if (je < 0) je = je + n

          !  loop across transforms
          !  ----------------------

!voption indep(z)
          do j = 0, (lot - 1) * jump, jump
             zja = z(ja + j)
             zjb = z(jb + j)
             zjc = z(jc + j)
             zjd = z(jd + j)
             zje = z(je + j)

             z1 = zjb + zje
             z2 = zjc + zjd
             z3 = zjb - zje
             z4 = zjc - zjd
             z5 = z1 + z2
             z6 = c1 * (z1 - z2)
             z7 = zja - dqtr * z5
             z8 = z7 + z6
             z9 = z7 - z6
             z10 = c3 * z3 - c2 * z4
             z11 = c2 * z3 + c3 * z4

             z(ja + j) = zja + z5
             z(jb + j) = z8 + zi * z11
             z(jc + j) = z9 + zi * z10
             z(jd + j) = z9 - zi * z10
             z(je + j) = z8 - zi * z11
          end do
          !-----( end of loop across transforms )

          ja = ja + jstepx
          if (ja < 0) ja = ja + n
       end do
    end do
    !-----( end of loop along transforms )

  end subroutine rad5I

  subroutine rad5Itwid(z, trigs, la, n, jump, lot, inq, jstepx, isign)

    implicit none

    complex(dp), intent(inout) :: z(0:*)
    complex(dp), intent(in) :: trigs(0:*)
    integer, intent(in) :: la, n, jump, lot
    integer, intent(in) :: inq, jstepx, isign

    integer :: ja, jb, jc, jd, je
    integer :: j, k, kk, jjj, mu, nu
    integer :: jstep, jstepl
    real(dp) :: c1, c2, c3, cswap
    complex(dp) :: zja, zjb, zjc, zjd, zje
    complex(dp) :: z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11

    jstep = n / (5 * la)
    jstepl = jstep - n
    mu = mod(inq, 5)
    if (isign == -1) mu = 5 - mu

    c1 = drt5 * dqtr
    c2 = sin(0.4d0 * pi)
    c3 = sin(0.2d0 * pi)
    if (mu == 2 .or. mu == 3) then
       c1 = -c1
       cswap = c2
       c2 = c3
       c3 = cswap
    end if
    if (mu == 3 .or. mu == 4) c2 = -c2
    if (mu == 2 .or. mu == 4) c3 = -c3

    kk = la

    !  loop on nonzero k
    !  -----------------
    do k = inq, jstep - inq, inq

       !  loop along transform
       !  --------------------
       do jjj = k, n - 1, 5 * jstep
          ja = jjj

          !     "transverse" loop
          !     -----------------
          do nu = 1, inq
             jb = ja + jstepl
             if (jb < 0) jb = jb + n
             jc = jb + jstepl
             if (jc < 0) jc = jc + n
             jd = jc + jstepl
             if (jd < 0) jd = jd + n
             je = jd + jstepl
             if (je < 0) je = je + n

             !  loop across transforms
             !  ----------------------
!voption indep(z)
             do j = 0, (lot - 1) * jump, jump
                zja = z(ja + j)
                zjb = z(jb + j)
                zjc = z(jc + j)
                zjd = z(jd + j)
                zje = z(je + j)

                z1 = zjb + zje
                z2 = zjc + zjd
                z3 = zjb - zje
                z4 = zjc - zjd
                z5 = z1 + z2
                z6 = c1 * (z1 - z2)
                z7 = zja - dqtr * z5
                z8 = z7 + z6
                z9 = z7 - z6
                z10 = c3 * z3 - c2 * z4
                z11 = c2 * z3 + c3 * z4

                z(ja + j) = zja + z5
                z(jb + j) = z8 + zi * z11
                z(jc + j) = z9 + zi * z10
                z(jd + j) = z9 - zi * z10
                z(je + j) = z8 - zi * z11
             end do

             call twiddle5(z, lot, jump, jb, jc, jd, je, isign, kk, trigs)
             !-----( end of loop across transforms )

             ja = ja + jstepx
             if (ja < 0) ja = ja + n
          end do
       end do
       !-----( end of loop along transforms )
       kk = kk + la
    end do

  end subroutine rad5Itwid

  subroutine rad5II(z, la, n, jump, lot, inq, jstepx, isign, work)

    implicit none

    complex(dp), intent(inout) :: z(0:*)
    integer, intent(in) :: la, n, jump
    integer, intent(in) :: lot, inq, jstepx, isign
    complex(dp) :: work(10, 0:lot - 1)

    integer :: ja, jb, jc, jd, je, jf, jg, jh, ji
    integer :: jj, jk, jl, jm, jn, jo, jp, jq
    integer :: jr, js, jt, ju, jv, jw, jx, jy
    integer :: j, j1, jjj, ll, mu, nu
    integer :: jstep, jstepl, laincl
    real(dp) :: c1, c2, c3, cswap
    complex(dp) :: zja, zjb, zjc, zjd, zje
    complex(dp) :: z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11

    jstep = n / (5 * la)
    jstepl = jstep - n
    laincl = la * inq - n
    mu = mod(inq, 5)
    if (isign == -1) mu = 5 - mu

    c1 = drt5 * dqtr
    c2 = sin(0.4d0 * pi)
    c3 = sin(0.2d0 * pi)
    if (mu == 2 .or. mu == 3) then
       c1 = -c1
       cswap = c2
       c2 = c3
       c3 = cswap
    end if
    if (mu == 3 .or. mu == 4) c2 = -c2
    if (mu == 2 .or. mu == 4) c3 = -c3

    !  k = 0 (no twiddle factors)
    !  --------------------------

    do ll = 0, (la - 1) * inq, 5 * jstep

       do jjj = ll, n - 1, 5 * la * inq
          ja = jjj

          !     "transverse" loop
          !     -----------------
          do nu = 1, inq
             jb = ja + jstepl
             if (jb < 0) jb = jb + n
             jc = jb + jstepl
             if (jc < 0) jc = jc + n
             jd = jc + jstepl
             if (jd < 0) jd = jd + n
             je = jd + jstepl
             if (je < 0) je = je + n
             jf = ja + laincl
             if (jf < 0) jf = jf + n
             jg = jf + jstepl
             if (jg < 0) jg = jg + n
             jh = jg + jstepl
             if (jh < 0) jh = jh + n
             ji = jh + jstepl
             if (ji < 0) ji = ji + n
             jj = ji + jstepl
             if (jj < 0) jj = jj + n
             jk = jf + laincl
             if (jk < 0) jk = jk + n
             jl = jk + jstepl
             if (jl < 0) jl = jl + n
             jm = jl + jstepl
             if (jm < 0) jm = jm + n
             jn = jm + jstepl
             if (jn < 0) jn = jn + n
             jo = jn + jstepl
             if (jo < 0) jo = jo + n
             jp = jk + laincl
             if (jp < 0) jp = jp + n
             jq = jp + jstepl
             if (jq < 0) jq = jq + n
             jr = jq + jstepl
             if (jr < 0) jr = jr + n
             js = jr + jstepl
             if (js < 0) js = js + n
             jt = js + jstepl
             if (jt < 0) jt = jt + n
             ju = jp + laincl
             if (ju < 0) ju = ju + n
             jv = ju + jstepl
             if (jv < 0) jv = jv + n
             jw = jv + jstepl
             if (jw < 0) jw = jw + n
             jx = jw + jstepl
             if (jx < 0) jx = jx + n
             jy = jx + jstepl
             if (jy < 0) jy = jy + n

             !  loop across transforms
             !  ----------------------

             j1 = 0
!voption indep(z)
             do j = 0, (lot-1)*jump, jump
                work(1, j1) = z(jf + j)
                work(2, j1) = z(jk + j)
                work(3, j1) = z(jl + j)
                work(4, j1) = z(jp + j)
                work(5, j1) = z(jq + j)
                work(6, j1) = z(jr + j)
                work(7, j1) = z(ju + j)
                work(8, j1) = z(jv + j)
                work(9, j1) = z(jw + j)
                work(10, j1) = z(jx + j)
                j1 = j1 + 1
             end do

!voption indep(z)
             do j = 0, (lot - 1) * jump, jump
                zja = z(ja + j)
                zjb = z(jb + j)
                zjc = z(jc + j)
                zjd = z(jd + j)
                zje = z(je + j)

                z1 = zjb + zje
                z2 = zjc + zjd
                z3 = zjb - zje
                z4 = zjc - zjd
                z5 = z1 + z2
                z6 = c1 * (z1 - z2)
                z7 = zja - dqtr * z5
                z8 = z7 + z6
                z9 = z7 - z6
                z10 = c3 * z3 - c2 * z4
                z11 = c2 * z3 + c3 * z4

                z(ja + j) = zja + z5
                z(jf + j) = z8 + zi * z11
                z(jk + j) = z9 + zi * z10
                z(jp + j) = z9 - zi * z10
                z(ju + j) = z8 - zi * z11
             end do

             j1 = 0
!voption indep(z)
             do j = 0, (lot - 1) * jump, jump
                zja = work(1, j1)
                zjb = z(jg + j)
                zjc = z(jh + j)
                zjd = z(ji + j)
                zje = z(jj + j)
                j1 = j1 + 1

                z1 = zjb + zje
                z2 = zjc + zjd
                z3 = zjb - zje
                z4 = zjc - zjd
                z5 = z1 + z2
                z6 = c1 * (z1 - z2)
                z7 = zja - dqtr * z5
                z8 = z7 + z6
                z9 = z7 - z6
                z10 = c3 * z3 - c2 * z4
                z11 = c2 * z3 + c3 * z4

                z(jb + j) = zja + z5
                z(jg + j) = z8 + zi * z11
                z(jl + j) = z9 + zi * z10
                z(jq + j) = z9 - zi * z10
                z(jv + j) = z8 - zi * z11
             end do

             j1 = 0
!voption indep(z)
             do j = 0, (lot - 1) * jump, jump
                zja = work(2, j1)
                zjb = work(3, j1)
                zjc = z(jm + j)
                zjd = z(jn + j)
                zje = z(jo + j)
                j1 = j1 + 1

                z1 = zjb + zje
                z2 = zjc + zjd
                z3 = zjb - zje
                z4 = zjc - zjd
                z5 = z1 + z2
                z6 = c1 * (z1 - z2)
                z7 = zja - dqtr * z5
                z8 = z7 + z6
                z9 = z7 - z6
                z10 = c3 * z3 - c2 * z4
                z11 = c2 * z3 + c3 * z4

                z(jc + j) = zja + z5
                z(jh + j) = z8 + zi * z11
                z(jm + j) = z9 + zi * z10
                z(jr + j) = z9 - zi * z10
                z(jw + j) = z8 - zi * z11
             end do

             j1 = 0
!voption indep(z)
             do j = 0, (lot - 1) * jump, jump
                zja = work(4, j1)
                zjb = work(5, j1)
                zjc = work(6, j1)
                zjd = z(js + j)
                zje = z(jt + j)
                j1 = j1 + 1

                z1 = zjb + zje
                z2 = zjc + zjd
                z3 = zjb - zje
                z4 = zjc - zjd
                z5 = z1 + z2
                z6 = c1 * (z1 - z2)
                z7 = zja - dqtr * z5
                z8 = z7 + z6
                z9 = z7 - z6
                z10 = c3 * z3 - c2 * z4
                z11 = c2 * z3 + c3 * z4

                z(jd + j) = zja + z5
                z(ji + j) = z8 + zi * z11
                z(jn + j) = z9 + zi * z10
                z(js + j) = z9 - zi * z10
                z(jx + j) = z8 - zi * z11
             end do

             j1 = 0
!voption indep(z)
             do j = 0, (lot - 1) * jump, jump
                zja = work(7, j1)
                zjb = work(8, j1)
                zjc = work(9, j1)
                zjd = work(10, j1)
                zje = z(jy + j)
                j1 = j1 + 1

                z1 = zjb + zje
                z2 = zjc + zjd
                z3 = zjb - zje
                z4 = zjc - zjd
                z5 = z1 + z2
                z6 = c1 * (z1 - z2)
                z7 = zja - dqtr * z5
                z8 = z7 + z6
                z9 = z7 - z6
                z10 = c3 * z3 - c2 * z4
                z11 = c2 * z3 + c3 * z4

                z(je + j) = zja + z5
                z(jj + j) = z8 + zi * z11
                z(jo + j) = z9 + zi * z10
                z(jt + j) = z9 - zi * z10
                z(jy + j) = z8 - zi * z11
             end do

             !-----(end of loop across transforms)

             ja = ja + jstepx
             if (ja < 0) ja = ja + n
          end do

       end do
    end do

  end subroutine rad5II

  subroutine rad5IItwid(z, trigs, la, n, jump, lot, inq, jstepx, isign, work)

    implicit none

    complex(dp), intent(inout) :: z(0:*)
    complex(dp), intent(in) :: trigs(0:*)
    integer, intent(in) :: la, n, jump
    integer, intent(in) :: lot, inq, jstepx, isign
    complex(dp) :: work(10, 0:lot - 1)

    integer :: ja, jb, jc, jd, je, jf, jg, jh, ji
    integer :: jj, jk, jl, jm, jn, jo, jp, jq
    integer :: jr, js, jt, ju, jv, jw, jx, jy
    integer :: j, j1, k, jjj, kk, ll, mu, nu
    integer :: jstep, jstepl, laincl
    real(dp) :: c1, c2, c3, cswap
    complex(dp) :: zja, zjb, zjc, zjd, zje
    complex(dp) :: z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11

    jstep = n / (5 * la)
    jstepl = jstep - n
    laincl = la * inq - n
    kk = la
    mu = mod(inq, 5)
    if (isign == -1) mu = 5 - mu

    c1 = drt5 * dqtr
    c2 = sin(0.4d0 * pi)
    c3 = sin(0.2d0 * pi)
    if (mu == 2 .or. mu == 3) then
       c1 = -c1
       cswap = c2
       c2 = c3
       c3 = cswap
    end if
    if (mu == 3 .or. mu == 4) c2 = -c2
    if (mu == 2 .or. mu == 4) c3 = -c3

    !  loop on nonzero k
    !  -----------------
    do k = inq, jstep - inq, inq

       !  double loop along first transform in block
       !  ------------------------------------------
       do ll = k, (la - 1) * inq, 5 * jstep

          do jjj = ll, n - 1, 5 * la * inq
             ja = jjj

             !     "transverse" loop
             !     -----------------
             do nu = 1, inq
                jb = ja + jstepl
                if (jb < 0) jb = jb + n
                jc = jb + jstepl
                if (jc < 0) jc = jc + n
                jd = jc + jstepl
                if (jd < 0) jd = jd + n
                je = jd + jstepl
                if (je < 0) je = je + n
                jf = ja + laincl
                if (jf < 0) jf = jf + n
                jg = jf + jstepl
                if (jg < 0) jg = jg + n
                jh = jg + jstepl
                if (jh < 0) jh = jh + n
                ji = jh + jstepl
                if (ji < 0) ji = ji + n
                jj = ji + jstepl
                if (jj < 0) jj = jj + n
                jk = jf + laincl
                if (jk < 0) jk = jk + n
                jl = jk + jstepl
                if (jl < 0) jl = jl + n
                jm = jl + jstepl
                if (jm < 0) jm = jm + n
                jn = jm + jstepl
                if (jn < 0) jn = jn + n
                jo = jn + jstepl
                if (jo < 0) jo = jo + n
                jp = jk + laincl
                if (jp < 0) jp = jp + n
                jq = jp + jstepl
                if (jq < 0) jq = jq + n
                jr = jq + jstepl
                if (jr < 0) jr = jr + n
                js = jr + jstepl
                if (js < 0) js = js + n
                jt = js + jstepl
                if (jt < 0) jt = jt + n
                ju = jp + laincl
                if (ju < 0) ju = ju + n
                jv = ju + jstepl
                if (jv < 0) jv = jv + n
                jw = jv + jstepl
                if (jw < 0) jw = jw + n
                jx = jw + jstepl
                if (jx < 0) jx = jx + n
                jy = jx + jstepl
                if (jy < 0) jy = jy + n

                !  loop across transforms
                !  ----------------------

                j1 = 0
!voption indep(z)
                do j = 0, (lot - 1) * jump, jump
                   work(1, j1) = z(jf + j)
                   work(2, j1) = z(jk + j)
                   work(3, j1) = z(jl + j)
                   work(4, j1) = z(jp + j)
                   work(5, j1) = z(jq + j)
                   work(6, j1) = z(jr + j)
                   work(7, j1) = z(ju + j)
                   work(8, j1) = z(jv + j)
                   work(9, j1) = z(jw + j)
                   work(10, j1) = z(jx + j)
                   j1 = j1 + 1
                end do

!voption indep(z)
                do j = 0, (lot - 1) * jump, jump
                   zja = z(ja + j)
                   zjb = z(jb + j)
                   zjc = z(jc + j)
                   zjd = z(jd + j)
                   zje = z(je + j)

                   z1 = zjb + zje
                   z2 = zjc + zjd
                   z3 = zjb - zje
                   z4 = zjc - zjd
                   z5 = z1 + z2
                   z6 = c1 * (z1 - z2)
                   z7 = zja - dqtr * z5
                   z8 = z7 + z6
                   z9 = z7 - z6
                   z10 = c3 * z3 - c2 * z4
                   z11 = c2 * z3 + c3 * z4

                   z(ja + j) = zja + z5
                   z(jf + j) = z8 + zi * z11
                   z(jk + j) = z9 + zi * z10
                   z(jp + j) = z9 - zi * z10
                   z(ju + j) = z8 - zi * z11
                end do

                call twiddle5(z, lot, jump, jf, jk, jp, ju, isign, kk, trigs)

                j1 = 0
!voption indep(z)
                do j = 0, (lot - 1) * jump, jump
                   zja = work(1, j1)
                   zjb = z(jg + j)
                   zjc = z(jh + j)
                   zjd = z(ji + j)
                   zje = z(jj + j)
                   j1 = j1 + 1

                   z1 = zjb + zje
                   z2 = zjc + zjd
                   z3 = zjb - zje
                   z4 = zjc - zjd
                   z5 = z1 + z2
                   z6 = c1 * (z1 - z2)
                   z7 = zja - dqtr * z5
                   z8 = z7 + z6
                   z9 = z7 - z6
                   z10 = c3 * z3 - c2 * z4
                   z11 = c2 * z3 + c3 * z4

                   z(jb + j) = zja + z5
                   z(jg + j) = z8 + zi * z11
                   z(jl + j) = z9 + zi * z10
                   z(jq + j) = z9 - zi * z10
                   z(jv + j) = z8 - zi * z11

                end do

                call twiddle5(z, lot, jump, jg, jl, jq, jv, isign, kk, trigs)

                j1 = 0
!voption indep(z)
                do j = 0, (lot - 1) * jump, jump
                   zja = work(2, j1)
                   zjb = work(3, j1)
                   zjc = z(jm + j)
                   zjd = z(jn + j)
                   zje = z(jo + j)
                   j1 = j1 + 1

                   z1 = zjb + zje
                   z2 = zjc + zjd
                   z3 = zjb - zje
                   z4 = zjc - zjd
                   z5 = z1 + z2
                   z6 = c1 * (z1 - z2)
                   z7 = zja - dqtr * z5
                   z8 = z7 + z6
                   z9 = z7 - z6
                   z10 = c3 * z3 - c2 * z4
                   z11 = c2 * z3 + c3 * z4

                   z(jc + j) = zja + z5
                   z(jh + j) = z8 + zi * z11
                   z(jm + j) = z9 + zi * z10
                   z(jr + j) = z9 - zi * z10
                   z(jw + j) = z8 - zi * z11

                end do

                call twiddle5(z, lot, jump, jh, jm, jr, jw, isign, kk, trigs)

                j1 = 0
!voption indep(z)
                do j = 0, (lot - 1) * jump, jump
                   zja = work(4, j1)
                   zjb = work(5, j1)
                   zjc = work(6, j1)
                   zjd = z(js + j)
                   zje = z(jt + j)
                   j1 = j1 + 1

                   z1 = zjb + zje
                   z2 = zjc + zjd
                   z3 = zjb - zje
                   z4 = zjc - zjd
                   z5 = z1 + z2
                   z6 = c1 * (z1 - z2)
                   z7 = zja - dqtr * z5
                   z8 = z7 + z6
                   z9 = z7 - z6
                   z10 = c3 * z3 - c2 * z4
                   z11 = c2 * z3 + c3 * z4

                   z(jd + j) = zja + z5
                   z(ji + j) = z8 + zi * z11
                   z(jn + j) = z9 + zi * z10
                   z(js + j) = z9 - zi * z10
                   z(jx + j) = z8 - zi * z11

                end do

                call twiddle5(z, lot, jump, ji, jn, js, jx, isign, kk, trigs)

                j1 = 0
!voption indep(z)
                do j = 0, (lot - 1) * jump, jump
                   zja = work(7, j1)
                   zjb = work(8, j1)
                   zjc = work(9, j1)
                   zjd = work(10, j1)
                   zje = z(jy + j)
                   j1 = j1 + 1

                   z1 = zjb + zje
                   z2 = zjc + zjd
                   z3 = zjb - zje
                   z4 = zjc - zjd
                   z5 = z1 + z2
                   z6 = c1 * (z1 - z2)
                   z7 = zja - dqtr * z5
                   z8 = z7 + z6
                   z9 = z7 - z6
                   z10 = c3 * z3 - c2 * z4
                   z11 = c2 * z3 + c3 * z4

                   z(je + j) = zja + z5
                   z(jj + j) = z8 + zi * z11
                   z(jo + j) = z9 + zi * z10
                   z(jt + j) = z9 - zi * z10
                   z(jy + j) = z8 - zi * z11
                end do

                call twiddle5(z, lot, jump, jj, jo, jt, jy, isign, kk, trigs)

                !-----(end of loop across transforms)

                ja = ja + jstepx
                if (ja < 0) ja = ja + n
             end do

          end do
       end do

       !-----( end of double loop for this k )
       kk = kk + la
    end do

  end subroutine rad5IItwid

  subroutine twiddle5(z, lot, jump, jb, jc, jd, je, isign, kk, trigs)

    implicit none

    complex(dp), intent(inout) :: z(0:*)
    complex(dp),  intent(in)    :: trigs(0:*)
    integer, intent(in) :: jump, lot, isign, kk
    integer, intent(in) :: jb, jc, jd, je

    integer :: j
    complex(dp) :: exp1, exp2, exp3, exp4

    if (isign == 1) then
       exp1 = trigs(kk)
       exp2 = trigs(2 * kk)
       exp3 = trigs(3 * kk)
       exp4 = trigs(4 * kk)
    else
       exp1 = conjg(trigs(kk))
       exp2 = conjg(trigs(2 * kk))
       exp3 = conjg(trigs(3 * kk))
       exp4 = conjg(trigs(4 * kk))
    end if

!voption indep(z)
    do j = 0, (lot - 1) * jump, jump
       z(jb + j) = exp1 * z(jb + j)
       z(jc + j) = exp2 * z(jc + j)
       z(jd + j) = exp3 * z(jd + j)
       z(je + j) = exp4 * z(je + j)
    end do

  end subroutine twiddle5

end subroutine zgpfa5f
