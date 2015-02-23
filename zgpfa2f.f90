!
subroutine zgpfa2f(z, trigs, jump, n, mm, lot, isign)

  use constants
  implicit none

  complex(dp), intent(inout) :: z(0:*)
  complex(dp), intent(in) :: trigs(0:*)
  integer, intent(in) :: jump, n, mm, lot, isign

  integer :: jstepx, m, la, mh, ipass, inq, n2
  logical :: m2, m8

  n2 = 2**mm
  inq = n / n2
  jstepx = n2 - n

  m2 = .false.
  m8 = .false.
  if (mod(mm, 2) == 0) then
     m = mm / 2
  else if (mod(mm, 4) == 1) then
     m = (mm - 1) / 2
     m2 = .true.
  else if (mod(mm, 4) == 3) then
     m = (mm - 3) / 2
     m8 = .true.
  end if
  mh = (m + 1) / 2

  la = 1

  !  loop on type I radix-4 passes
  !  -----------------------------

  do ipass = 1, mh

     call rad4I(z, la, n, jump, lot, inq, jstepx, isign)

     !  finished if n2=4
     !  ----------------
     if (n2 == 4) return

     call rad4Itwid(z, trigs, la, n, jump, lot, inq, jstepx, isign)
     la = 4 * la

  end do

  !  central radix-2 pass
  !  --------------------
  if (m2) then

     call rad2(z, la, n, jump, lot, inq, jstepx)

     !  finished if n2=2
     !  ----------------
     if (n2 == 2) return

     call rad2twid(z, trigs, la, n, jump, lot, inq, jstepx, isign, n2)
     la = 2 * la

  end if

  !  central radix-8 pass
  !  --------------------
  if (m8) then

     call rad8I(z, la, n, jump, lot, inq, jstepx, isign)
     call rad8II(z, la, n, jump, lot, inq, jstepx, isign)

     !  finished if n2=8
     !  ----------------
     if (n2 == 8) return

     call rad8IItwid(z, trigs, la, n, jump, lot, inq, jstepx, isign)
     la = 8 * la

  end if

  !  loop on type II radix-4 passes
  !  ------------------------------

  do ipass = mh + 1, m

     call rad4II(z, la, n, jump, lot, inq, jstepx, isign)

     !  finished if last pass
     !  ---------------------
     if (ipass == m) return

     call rad4IItwid(z, trigs, la, n, jump, lot, inq, jstepx, isign)
     la = 4 * la

  end do

  return

contains

  subroutine rad4I(z, la, n, jump, lot, inq, jstepx, isign)

    implicit none

    complex(dp), intent(inout) :: z(0:*)
    integer, intent(in) :: la, n, jump
    integer, intent(in) :: lot, inq
    integer, intent(in) :: jstepx, isign

    integer :: ja, jb, jc, jd
    integer :: jjj, j, nu, mu
    integer :: jstep, jstepl
    complex(dp) :: zja, zjb, zjc, zjd
    complex(dp) :: z0, z1, z2, z3
    real(dp) :: ss

    jstep = n / (4 * la)
    jstepl = jstep - n
    mu = mod(inq, 4)
    if (isign == -1) mu = 4 - mu
    ss = done
    if (mu == 3) ss = dmone

    !  k = 0 loop (no twiddle factors)
    !  -------------------------------
    do jjj = 0, n - 1, 4 * jstep
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

          !  loop across transforms
          !  ----------------------
!voption indep(z)
          do j = 0, (lot - 1) * jump, jump
             zja = z(ja + j)
             zjb = z(jb + j)
             zjc = z(jc + j)
             zjd = z(jd + j)

             z0 = zja + zjc
             z1 = zjb + zjd
             z2 = zja - zjc
             z3 = ss * (zjb - zjd)

             z(ja + j) = z0 + z1
             z(jb + j) = z2 + zi * z3
             z(jc + j) = z0 - z1
             z(jd + j) = z2 - zi * z3
          end do
          ja = ja + jstepx
          if (ja < 0) ja = ja + n
       end do
    end do

    return

  end subroutine rad4I

  subroutine rad4Itwid(z, trigs, la, n, jump, lot, inq, jstepx, isign)

    implicit none

    complex(dp), intent(inout) :: z(0:*)
    complex(dp), intent(in) :: trigs(0:*)
    integer, intent(in) :: la, n, jump
    integer, intent(in) :: lot, inq
    integer, intent(in) :: jstepx, isign

    integer :: ja, jb, jc, jd
    integer :: jjj, j, nu, mu, k, kk
    integer :: jstep, jstepl

    complex(dp) :: zja, zjb, zjc, zjd
    complex(dp) :: z0, z1, z2, z3
    complex(dp) :: exp1, exp2, exp3
    real(dp) :: ss

    jstep = n / (4 * la)
    jstepl = jstep - n
    mu = mod(inq, 4)
    if (isign == -1) mu = 4 - mu
    ss = done
    if (mu == 3) ss = dmone

    kk = la

    !  loop on nonzero k
    !  -----------------
    do k = inq, jstep - inq, inq
       if( isign == 1 )then
          exp1 = trigs(kk)
          exp2 = trigs(2 * kk)
          exp3 = trigs(3 * kk)
       else
          exp1 = conjg(trigs(kk))
          exp2 = conjg(trigs(2 * kk))
          exp3 = conjg(trigs(3 * kk))
       end if

       !  loop along transform
       !  --------------------
       do jjj = k, n - 1, 4 * jstep
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

             !  loop across transforms
             !  ----------------------
!voption indep(z)
             do j = 0, (lot - 1) * jump, jump
                zja = z(ja + j)
                zjb = z(jb + j)
                zjc = z(jc + j)
                zjd = z(jd + j)

                z0 = zja + zjc
                z1 = zjb + zjd
                z2 = zja - zjc
                z3 = ss * (zjb - zjd)

                z(ja + j) = z0 + z1
                z(jb + j) = exp1 * (z2 + zi * z3)
                z(jc + j) = exp2 * (z0 - z1 )
                z(jd + j) = exp3 * (z2 - zi * z3)
             end do
             !-----( end of loop across transforms )
             ja = ja + jstepx
             if (ja < 0) ja = ja + n
          end do
       end do
       !-----( end of loop along transforms )
       kk = kk + la
    end do
    !-----( end of loop on nonzero k )

    return

  end subroutine rad4Itwid

  subroutine rad4II(z, la, n, jump, lot, inq, jstepx, isign)

    implicit none

    complex(dp), intent(inout) :: z(0:*)
    integer, intent(in) :: la, n, jump
    integer, intent(in) :: lot, inq
    integer, intent(in) :: jstepx, isign

    integer :: ja, jb, jc, jd, je, jf, jg, jh
    integer :: ji, jj, jk, jl, jm, jn, jo, jp
    integer :: ll, jjj, j, nu, mu
    integer :: jstep, jstepl, laincl
    complex(dp) :: zja, zjb, zjc, zjd, zje, zjf, zjg, zjh
    complex(dp) :: zji, zjj, zjk, zjl, zjm, zjn, zjo, zjp
    complex(dp) :: z0, z1, z2, z3
    real(dp) :: ss

    jstep = n / (4 * la)
    jstepl = jstep - n
    mu = mod(inq, 4)
    if (isign == -1) mu = 4 - mu
    ss = done
    if (mu == 3) ss = dmone
    laincl = la * inq - n

    !  k=0 loop (no twiddle factors)
    !  -----------------------------
    do ll = 0, (la - 1) * inq, 4 * jstep

       do jjj = ll, n - 1, 4 * la * inq
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
             je = ja + laincl
             if (je < 0) je = je + n
             jf = je + jstepl
             if (jf < 0) jf = jf + n
             jg = jf + jstepl
             if (jg < 0) jg = jg + n
             jh = jg + jstepl
             if (jh < 0) jh = jh + n
             ji = je + laincl
             if (ji < 0) ji = ji + n
             jj = ji + jstepl
             if (jj < 0) jj = jj + n
             jk = jj + jstepl
             if (jk < 0) jk = jk + n
             jl = jk + jstepl
             if (jl < 0) jl = jl + n
             jm = ji + laincl
             if (jm < 0) jm = jm + n
             jn = jm + jstepl
             if (jn < 0) jn = jn + n
             jo = jn + jstepl
             if (jo < 0) jo = jo + n
             jp = jo + jstepl
             if (jp < 0) jp = jp + n

             !  loop across transforms
             !  ----------------------
!voption indep(z)
             do j = 0, (lot - 1) * jump, jump
                zja = z(ja + j)
                zjb = z(jb + j)
                zjc = z(jc + j)
                zjd = z(jd + j)
                zje = z(je + j)
                zjf = z(jf + j)
                zjg = z(jg + j)
                zjh = z(jh + j)
                zji = z(ji + j)
                zjj = z(jj + j)
                zjk = z(jk + j)
                zjl = z(jl + j)
                zjm = z(jm + j)
                zjn = z(jn + j)
                zjo = z(jo + j)
                zjp = z(jp + j)

                z0 = zja + zjc
                z1 = zjb + zjd
                z2 = zja - zjc
                z3 = ss * (zjb - zjd)

                z(ja + j) = z0 + z1
                z(je + j) = z2 + zi * z3
                z(ji + j) = z0 - z1
                z(jm + j) = z2 - zi * z3
                !----------------------
                z0 = zje + zjg
                z1 = zjf + zjh
                z2 = zje - zjg
                z3 = ss * (zjf - zjh)

                z(jb + j) = z0 + z1
                z(jf + j) = z2 + zi * z3
                z(jj + j) = z0 - z1
                z(jn + j) = z2 - zi * z3
                !----------------------
                z0 = zji + zjk
                z1 = zjj + zjl
                z2 = zji - zjk
                z3 = ss * (zjj - zjl)

                z(jc + j) = z0 + z1
                z(jg + j) = z2 + zi * z3
                z(jk + j) = z0 - z1
                z(jo + j) = z2 - zi * z3
                !----------------------
                z0 = zjm + zjo
                z1 = zjn + zjp
                z2 = zjm - zjo
                z3 = ss * (zjn - zjp)

                z(jd + j) = z0 + z1
                z(jh + j) = z2 + zi * z3
                z(jl + j) = z0 - z1
                z(jp + j) = z2 - zi * z3
             end do
             !-----( end of loop across transforms )
             ja = ja + jstepx
             if (ja < 0) ja = ja + n
          end do
       end do
    end do
    !-----( end of double loop for k=0 )

    return

  end subroutine rad4II

  subroutine rad4IItwid(z, trigs, la, n, jump, lot, inq, jstepx, isign)

    implicit none

    complex(dp), intent(inout) :: z(0:*)
    complex(dp), intent(in) :: trigs(0:*)
    integer, intent(in) :: la, n, jump
    integer, intent(in) :: lot, inq
    integer, intent(in) :: jstepx, isign

    integer :: ja, jb, jc, jd, je, jf, jg, jh
    integer :: ji, jj, jk, jl, jm, jn, jo, jp
    integer :: ll, jjj, j, nu, mu, k, kk
    integer :: jstep, jstepl, laincl

    complex(dp) :: zja, zjb, zjc, zjd, zje, zjf, zjg, zjh
    complex(dp) :: zji, zjj, zjk, zjl, zjm, zjn, zjo, zjp
    complex(dp) :: z0, z1, z2, z3
    complex(dp) :: exp1, exp2, exp3
    real(dp) :: ss

    jstep = n / (4 * la)
    jstepl = jstep - n
    mu = mod(inq, 4)
    if (isign == -1) mu = 4 - mu
    ss = done
    if (mu == 3) ss = dmone
    laincl = la * inq - n
    kk = la

    !     loop on nonzero k
    !     -----------------
    do k = inq, jstep - inq, inq
       if( isign  ==  1 )then
          exp1 = trigs(kk)
          exp2 = trigs(2 * kk)
          exp3 = trigs(3 * kk)
       else
          exp1 = conjg(trigs(kk))
          exp2 = conjg(trigs(2 * kk))
          exp3 = conjg(trigs(3 * kk))
       end if

       !  double loop along first transform in block
       !  ------------------------------------------
       do ll = k, (la - 1) * inq, 4 * jstep
          !
          do jjj = ll, n - 1, 4 * la * inq
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
                je = ja + laincl
                if (je < 0) je = je + n
                jf = je + jstepl
                if (jf < 0) jf = jf + n
                jg = jf + jstepl
                if (jg < 0) jg = jg + n
                jh = jg + jstepl
                if (jh < 0) jh = jh + n
                ji = je + laincl
                if (ji < 0) ji = ji + n
                jj = ji + jstepl
                if (jj < 0) jj = jj + n
                jk = jj + jstepl
                if (jk < 0) jk = jk + n
                jl = jk + jstepl
                if (jl < 0) jl = jl + n
                jm = ji + laincl
                if (jm < 0) jm = jm + n
                jn = jm + jstepl
                if (jn < 0) jn = jn + n
                jo = jn + jstepl
                if (jo < 0) jo = jo + n
                jp = jo + jstepl
                if (jp < 0) jp = jp + n

                !  loop across transforms
                !  ----------------------
!voption indep(z)
                do j = 0, (lot - 1) * jump, jump
                   zja = z(ja + j)
                   zjb = z(jb + j)
                   zjc = z(jc + j)
                   zjd = z(jd + j)
                   zje = z(je + j)
                   zjf = z(jf + j)
                   zjg = z(jg + j)
                   zjh = z(jh + j)
                   zji = z(ji + j)
                   zjj = z(jj + j)
                   zjk = z(jk + j)
                   zjl = z(jl + j)
                   zjm = z(jm + j)
                   zjn = z(jn + j)
                   zjo = z(jo + j)
                   zjp = z(jp + j)

                   z0 = zja + zjc
                   z1 = zjb + zjd
                   z2 = zja - zjc
                   z3 = ss * (zjb - zjd)

                   z(ja + j) = z0 + z1
                   z(je + j) = exp1 * (z2 + zi * z3)
                   z(ji + j) = exp2 * (z0 - z1)
                   z(jm + j) = exp3 * (z2 - zi * z3)
                   !----------------------------------------
                   z0 = zje + zjg
                   z1 = zjf + zjh
                   z2 = zje - zjg
                   z3 = ss * (zjf - zjh)

                   z(jb + j) = z0 + z1
                   z(jf + j) = exp1 * (z2 + zi * z3)
                   z(jj + j) = exp2 * (z0 - z1)
                   z(jn + j) = exp3 * (z2 - zi * z3)
                   !----------------------------------------
                   z0 = zji + zjk
                   z1 = zjj + zjl
                   z2 = zji - zjk
                   z3 = ss * (zjj - zjl)

                   z(jc + j) = z0 + z1
                   z(jg + j) = exp1 * (z2 + zi * z3)
                   z(jk + j) = exp2 * (z0 - z1)
                   z(jo + j) = exp3 * (z2 - zi * z3)
                   !----------------------------------------
                   z0 = zjm + zjo
                   z1 = zjn + zjp
                   z2 = zjm - zjo
                   z3 = ss * (zjn - zjp)

                   z(jd + j) = z0 + z1
                   z(jh + j) = exp1 * (z2 + zi * z3)
                   z(jl + j) = exp2 * (z0 - z1)
                   z(jp + j) = exp3 * (z2 - zi * z3)
                end do
                !-----(end of loop across transforms)
                ja = ja + jstepx
                if (ja < 0) ja = ja + n
             end do
          end do
       end do
       !-----( end of double loop for this k )
       kk = kk + la
    end do
    !-----( end of loop over values of k )

    return

  end subroutine rad4IItwid

  subroutine rad2(z, la, n, jump, lot, inq, jstepx)

    implicit none

    complex(dp), intent(inout) :: z(0:*)
    integer, intent(in) :: la, n, jump
    integer, intent(in) :: lot, inq, jstepx

    integer :: ja, jb
    integer :: jjj, j, nu
    integer :: jstep, jstepl
    complex(dp) :: zja, zjb

    jstep = n / (2 * la)
    jstepl = jstep - n

    !  k=0 loop (no twiddle factors)
    !  -----------------------------
    do jjj = 0, n - 1, 2 * jstep
       ja = jjj

       !     "transverse" loop
       !     -----------------
       do nu = 1, inq
          jb = ja + jstepl
          if (jb < 0) jb = jb + n

          !  loop across transforms
          !  ----------------------
!voption indep(z)
          do j = 0, (lot - 1) * jump, jump
             zja = z(ja + j)
             zjb = z(jb + j)
             z(ja + j) = zja + zjb
             z(jb + j) = zja - zjb
          end do
          !-----(end of loop across transforms)
          ja = ja + jstepx
          if (ja < 0) ja = ja + n
       end do
    end do

    return

  end subroutine rad2

  subroutine rad2twid(z, trigs, la, n, jump, lot, inq, jstepx, isign, n2)

    implicit none

    complex(dp), intent(inout) :: z(0:*)
    complex(dp), intent(in) :: trigs(0:*)
    integer, intent(in) :: la, n, jump
    integer, intent(in) :: lot, inq
    integer, intent(in) :: jstepx, isign, n2

    integer :: ja, jb
    integer :: jjj, j, nu, mu, kk, k
    integer :: jstep, jstepl
    complex(dp) :: zja, zjb, exp1
    real(dp) :: ss

    jstep = n / (2 * la)
    jstepl = jstep - n
    mu = mod(inq, 4)
    if (isign == -1) mu = 4 - mu
    ss = done
    if (mu == 3) ss = dmone
    kk = la

    !  loop on nonzero k
    !  -----------------
    do k = inq, jstep - inq, inq
       if (isign  ==  1) then
          exp1 = trigs(kk)
       else
          exp1 = conjg(trigs(kk))
       end if

       !  loop along transforms
       !  ---------------------
       do jjj = k, n - 1, 2 * jstep
          ja = jjj

          !     "transverse" loop
          !     -----------------
          do nu = 1, inq
             jb = ja + jstepl
             if (jb < 0) jb = jb + n

             !  loop across transforms
             !  ----------------------
             if (kk == n2 / 4) then
!voption indep(z)
                do j = 0, (lot - 1) * jump, jump
                   zja = z(ja + j)
                   zjb = z(jb + j)
                   z(ja + j) = zja + zjb
                   z(jb + j) = cmplx(dzero, ss, dp) * (zja - zjb)
                end do
             else
!voption indep(z)
                do j = 0, (lot - 1) * jump, jump
                   zja = z(ja + j)
                   zjb = z(jb + j)
                   z(ja + j) = zja + zjb
                   z(jb + j) = exp1 * (zja - zjb)
                end do
             end if

             !-----(end of loop across transforms)
             ja = ja + jstepx
             if (ja < 0) ja = ja + n
          end do
       end do
       !-----(end of loop along transforms)
       kk = kk + la
    end do
    !-----(end of loop on nonzero k)

    return

  end subroutine rad2twid

  subroutine rad8I(z, la, n, jump, lot, inq, jstepx, isign)

    implicit none

    complex(dp), intent(inout) :: z(0:*)
    integer, intent(in) :: la, n, jump, lot
    integer, intent(in) :: inq, isign, jstepx

    integer :: ja, jb, jc, jd, je, jf, jg, jh
    complex(dp) :: zja, zjb, zjc, zjd, zje, zjf, zjg, zjh
    complex(dp) :: z2, z3
    real(dp) :: c1, c2, c3
    integer :: jstep, jstepl, mu, jjj, j, nu, k

    jstep = n / (8 * la)
    jstepl = jstep - n
    mu = mod(inq, 8)
    if (isign == -1) mu = 8 - mu
    c1 = done
    if (mu == 3 .or. mu == 7) c1 = dmone
    c2 = oort2
    if (mu == 3 .or. mu == 5) c2 = -c2
    c3 = c1 * c2

    !  stage 1
    !  -------
    do k = 0, jstep - inq, inq
       do jjj = k, n - 1, 8 * jstep
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
             jf = je + jstepl
             if (jf < 0) jf = jf + n
             jg = jf + jstepl
             if (jg < 0) jg = jg + n
             jh = jg + jstepl
             if (jh < 0) jh = jh + n

!voption indep(z)
             do j = 0, (lot - 1) * jump, jump
                zja = z(ja + j)
                zjb = z(jb + j)
                zjc = z(jc + j)
                zjd = z(jd + j)
                zje = z(je + j)
                zjf = z(jf + j)
                zjg = z(jg + j)
                zjh = z(jh + j)

                z2 = zjb - zjf
                z3 = zjd - zjh

                z(ja + j) = zja + zje
                z(jb + j) = zja - zje
                z(jc + j) = zjb + zjf
                z(jd + j) = c2 * (z2 - z3)
                z(je + j) = zjc + zjg
                z(jf + j) = c1 * (zjc - zjg)
                z(jg + j) = zjd + zjh
                z(jh + j) = c3 * (z2 + z3)
             end do

             ja = ja + jstepx
             if (ja < 0) ja = ja + n
          end do
       end do
    end do

    return

  end subroutine rad8I

  subroutine rad8II(z, la, n, jump, lot, inq, jstepx, isign)

    implicit none

    complex(dp), intent(inout) :: z(0:*)
    integer, intent(in) :: la, n, jump, lot
    integer, intent(in) :: inq, isign, jstepx

    integer :: ja, jb, jc, jd, je, jf, jg, jh
    complex(dp) :: zja, zjb, zjc, zjd, zje, zjf, zjg, zjh
    complex(dp) :: z0, z1, z2, z3
    real(dp) :: c1, c2
    integer :: jstep, jstepl, mu, jjj, j, nu

    jstep = n / (8 * la)
    jstepl = jstep - n
    mu = mod(inq, 8)
    if (isign == -1) mu = 8 - mu
    c1 = done
    if (mu == 3 .or. mu == 7) c1 = dmone
    c2 = oort2
    if (mu == 3 .or. mu == 5) c2 = -c2

    !  stage 2
    !  -------

    !  k=0 (no twiddle factors)
    !  ------------------------
    do jjj = 0, n - 1, 8 * jstep
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
          jf = je + jstepl
          if (jf < 0) jf = jf + n
          jg = jf + jstepl
          if (jg < 0) jg = jg + n
          jh = jg + jstepl
          if (jh < 0) jh = jh + n

!voption indep(z)
          do j = 0, (lot - 1) * jump, jump
             zja = z(ja + j)
             zjc = z(jc + j)
             zje = z(je + j)
             zjg = z(jg + j)

             z0 = zja + zje
             z2 = zja - zje
             z1 = zjc + zjg
             z3 = c1 * (zjc - zjg)

             z(ja + j) = z0 + z1
             z(jc + j) = z2 + zi * z3
             z(je + j) = z0 - z1
             z(jg + j) = z2 - zi * z3
          end do

!voption indep(z)
          do j = 0, (lot - 1) * jump, jump
             zjb = z(jb + j)
             zjd = z(jd + j)
             zjf = z(jf + j)
             zjh = z(jh + j)

             z0 = zjb + zjd
             z2 = zjb - zjd
             z1 = zjf - zjh
             z3 = zjf + zjh

             z(jb + j) = z0 + zi * z3
             z(jd + j) = z2 - zi * z1
             z(jf + j) = z2 + zi * z1
             z(jh + j) = z0 - zi * z3
          end do

          ja = ja + jstepx
          if (ja < 0) ja = ja + n
       end do
    end do

    return

  end subroutine rad8II

  subroutine rad8IItwid(z, trigs, la, n, jump, lot, inq, jstepx, isign)

    implicit none

    complex(dp), intent(inout) :: z(0:*)
    complex(dp), intent(in) :: trigs(0:*)
    integer, intent(in) :: la, n, jump
    integer, intent(in) :: lot, inq
    integer, intent(in) :: jstepx, isign

    integer :: ja, jb, jc, jd, je, jf, jg, jh
    integer :: jjj, j, nu, k, kk, mu
    integer :: jstep, jstepl

    complex(dp) :: zja, zjb, zjc, zjd, zje, zjf, zjg, zjh
    complex(dp) :: z0, z1, z2, z3
    complex(dp) :: exp1, exp2, exp3, exp4, exp5, exp6, exp7
    real(dp) :: c1, c2

    jstep = n / (8 * la)
    jstepl = jstep - n
    mu = mod(inq, 8)
    if (isign == -1) mu = 8 - mu
    c1 = done
    if (mu == 3 .or. mu == 7) c1 = dmone
    c2 = oort2
    if (mu == 3 .or. mu == 5) c2 = -c2

    !  loop on nonzero k
    !  -----------------
    kk = la

    do k = inq, jstep - inq, inq

       if (isign  ==  1) then
          exp1 = trigs(kk)
          exp2 = trigs(2 * kk)
          exp3 = trigs(3 * kk)
          exp4 = trigs(4 * kk)
          exp5 = trigs(5 * kk)
          exp6 = trigs(6 * kk)
          exp7 = trigs(7 * kk)
       else
          exp1 = conjg(trigs(kk))
          exp2 = conjg(trigs(2 * kk))
          exp3 = conjg(trigs(3 * kk))
          exp4 = conjg(trigs(4 * kk))
          exp5 = conjg(trigs(5 * kk))
          exp6 = conjg(trigs(6 * kk))
          exp7 = conjg(trigs(7 * kk))
       end if

       do jjj = k, n - 1, 8 * jstep
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
             jf = je + jstepl
             if (jf < 0) jf = jf + n
             jg = jf + jstepl
             if (jg < 0) jg = jg + n
             jh = jg + jstepl
             if (jh < 0) jh = jh + n

!voption indep(z)
             do j = 0, (lot - 1) * jump, jump
                zja = z(ja + j)
                zjc = z(jc + j)
                zje = z(je + j)
                zjg = z(jg + j)

                z0 = zja + zje
                z2 = zja - zje
                z1 = zjc + zjg
                z3 = c1 * (zjc - zjg)

                z(ja + j) = z0 + z1
                z(jc + j) = exp2 * (z2 + zi * z3)
                z(je + j) = exp4 * (z0 - z1 )
                z(jg + j) = exp6 * (z2 - zi * z3)
             end do

!voption indep(z)
             do j = 0, (lot - 1) * jump, jump
                zjb = z(jb + j)
                zjd = z(jd + j)
                zjf = z(jf + j)
                zjh = z(jh + j)

                z0 = zjb + zjd
                z2 = zjb - zjd
                z1 = zjf - zjh
                z3 = zjf + zjh

                z(jb + j) = exp1 * (z0 + zi * z3)
                z(jd + j) = exp3 * (z2 - zi * z1)
                z(jh + j) = exp7 * (z0 - zi * z3)
                z(jf + j) = exp5 * (z2 + zi * z1)
             end do
             ja = ja + jstepx
             if (ja < 0) ja = ja + n
          end do
       end do
       kk = kk + la
    end do

    return

  end subroutine rad8IItwid

end subroutine zgpfa2f
