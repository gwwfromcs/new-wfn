!
subroutine zgpfa3f(z, trigs, jump, n, mm, lot, isign)

  use constants
  implicit none

  complex(dp), intent(inout) :: z(0:*)
  complex(dp), intent(in) :: trigs(0:*)
  integer, intent(in) :: jump, n, mm, lot, isign

  integer :: jstepx, inq, la, mh, n3, ipass

  n3 = 3**mm
  inq = n / n3
  jstepx = n3 - n
  mh = (mm + 1) / 2
  la = 1

  !     loop on type I radix-3 passes
  !     -----------------------------
  do ipass = 1, mh

     call rad3I(z, la, n, jump, lot, inq, jstepx, isign)

     !     finished if n3 = 3
     !     ------------------
     if (n3 == 3) return

     call rad3Itwid(z, trigs, la, n, jump, lot, inq, jstepx, isign)
     la = 3 * la

  end do

  !     loop on type II radix-3 passes
  !     ------------------------------

  do ipass = mh + 1, mm

     call rad3II(z, la, n, jump, lot, inq, jstepx, isign)

     !     finished if last pass
     !     ---------------------
     if (ipass == mm) return

     call rad3IItwid(z, trigs, la, n, jump, lot, inq, jstepx, isign)
     la = 3 * la

  end do

  return

contains

  subroutine rad3I(z, la, n, jump, lot, inq, jstepx, isign)

    implicit none

    complex(dp), intent(inout) :: z(0:*)
    integer, intent(in) :: la, n, jump
    integer, intent(in) :: lot, inq
    integer, intent(in) :: jstepx, isign

    integer :: ja, jb, jc
    integer :: jjj, j, nu, mu
    integer :: jstep, jstepl
    real(dp) :: c1
    complex(dp) :: zja, zjb, zjc
    complex(dp) :: z1, z2, z3

    jstep = n / (3 * la)
    jstepl = jstep - n
    mu = mod(inq, 3)
    if (isign == -1) mu = 3 - mu
    c1 = drt3 * dhalf
    if (mu == 2) c1 = -c1

    !     k = 0 loop (no twiddle factors)
    !     -------------------------------
    do jjj = 0, n - 1, 3 * jstep
       ja = jjj

       !     "transverse" loop
       !     -----------------
       do nu = 1, inq
          jb = ja + jstepl
          if (jb < 0) jb = jb + n
          jc = jb + jstepl
          if (jc < 0) jc = jc + n

          !     loop across transforms
          !     ----------------------

!voption indep(z)
          do j = 0, (lot - 1) * jump, jump
             zja = z(ja + j)
             zjb = z(jb + j)
             zjc = z(jc + j)
             z1 = zjb + zjc
             z2 = zja - dhalf * z1
             z3 = c1 * (zjb - zjc)
             z(ja + j) = zja + z1
             z(jb + j) = z2 + zi * z3
             z(jc + j) = z2 - zi * z3
          end do

          ja = ja + jstepx
          if (ja < 0) ja = ja + n
       end do
    end do

    return

  end subroutine rad3I

  subroutine rad3Itwid(z, trigs, la, n, jump, lot, inq, jstepx, isign)

    implicit none

    complex(dp), intent(inout) :: z(0:*)
    complex(dp), intent(in) :: trigs(0:*)
    integer, intent(in) :: la, n, jump
    integer, intent(in) :: lot, inq
    integer, intent(in) :: jstepx, isign

    integer :: ja, jb, jc
    integer :: jjj, j, k, kk, nu, mu
    integer :: jstep, jstepl
    real(dp) :: c1
    complex(dp) :: zja, zjb, zjc
    complex(dp) :: z1, z2, z3
    complex(dp) :: exp1, exp2

    jstep = n / (3 * la)
    jstepl = jstep - n
    mu = mod(inq, 3)
    if (isign == -1) mu = 3 - mu
    c1 = drt3 * dhalf
    if (mu == 2) c1 = -c1

    kk = la

    !     loop on nonzero k
    !     -----------------
    do k = inq, jstep - inq, inq

       if (isign  ==  1) then
          exp1 = trigs(kk)
          exp2 = trigs(2 * kk)
       else
          exp1 = conjg(trigs(kk))
          exp2 = conjg(trigs(2 * kk))
       end if

       !     loop along transform
       !     --------------------
       do jjj = k, n - 1, 3 * jstep
          ja = jjj

          !     "transverse" loop
          !     -----------------
          do nu = 1, inq
             jb = ja + jstepl
             if (jb < 0) jb = jb + n
             jc = jb + jstepl
             if (jc < 0) jc = jc + n

             !     loop across transforms
             !     ----------------------
!voption indep(z)
             do j = 0, (lot - 1) * jump, jump
                zja = z(ja + j)
                zjb = z(jb + j)
                zjc = z(jc + j)
                z1 = zjb + zjc
                z2 = zja - dhalf * z1
                z3 = c1 * (zjb - zjc)
                z(ja + j) = zja + z1
                z(jb + j) = exp1 * (z2 + zi * z3)
                z(jc + j) = exp2 * (z2 - zi * z3)
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

  end subroutine rad3Itwid

  subroutine rad3II(z, la, n, jump, lot, inq, jstepx, isign)

    implicit none

    complex(dp), intent(inout) :: z(0:*)
    integer, intent(in) :: la, n, jump
    integer, intent(in) :: lot, inq
    integer, intent(in) :: jstepx, isign

    integer :: ja, jb, jc, jd, je, jf, jg, jh, ji
    integer :: ll, jjj, j, nu, mu
    integer :: jstep, jstepl, laincl
    real(dp) :: c1
    complex(dp) :: zja, zjb, zjc, zjd, zje, zjf, zjg, zjh, zji
    complex(dp) :: z1, z2, z3

    jstep = n / (3 * la)
    jstepl = jstep - n
    laincl = la * inq - n
    mu = mod(inq, 3)
    if (isign == -1) mu = 3 - mu
    c1 = drt3 * dhalf
    if (mu == 2) c1 = -c1

    !     k=0 loop (no twiddle factors)
    !     -----------------------------
    do ll = 0, (la - 1) * inq, 3 * jstep

       do jjj = ll, n - 1, 3 * la * inq
          ja = jjj

          !     "transverse" loop
          !     -----------------
          do nu = 1, inq
             jb = ja + jstepl
             if (jb < 0) jb = jb + n
             jc = jb + jstepl
             if (jc < 0) jc = jc + n
             jd = ja + laincl
             if (jd < 0) jd = jd + n
             je = jd + jstepl
             if (je < 0) je = je + n
             jf = je + jstepl
             if (jf < 0) jf = jf + n
             jg = jd + laincl
             if (jg < 0) jg = jg + n
             jh = jg + jstepl
             if (jh < 0) jh = jh + n
             ji = jh + jstepl
             if (ji < 0) ji = ji + n

             !     loop across transforms
             !     ----------------------
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

                z1 = zjb + zjc
                z2 = zja - dhalf * z1
                z3 = c1 * (zjb - zjc)
                z(ja + j) = zja + z1
                z(jd + j) = z2 + zi * z3
                z(jg + j) = z2 - zi * z3
                !----------------------
                z1 = zje + zjf
                z2 = zjd - dhalf * z1
                z3 = c1 * (zje - zjf)
                z(jb + j) = zjd + z1
                z(je + j) = z2 + zi * z3
                z(jh + j) = z2 - zi * z3
                !----------------------
                z1 = zjh + zji
                z2 = zjg - dhalf * z1
                z3 = c1 * (zjh - zji)
                z(jc + j) = zjg + z1
                z(jf + j) = z2 + zi * z3
                z(ji + j) = z2 - zi * z3
             end do
             !-----( end of loop across transforms )
             ja = ja + jstepx
             if (ja < 0) ja = ja + n
          end do
       end do
    end do
    !-----( end of double loop for k=0 )

    return

  end subroutine rad3II

  subroutine rad3IItwid(z, trigs, la, n, jump, lot, inq, jstepx, isign)

    implicit none

    complex(dp), intent(inout) :: z(0:*)
    complex(dp), intent(in) :: trigs(0:*)
    integer, intent(in) :: la, n, jump
    integer, intent(in) :: lot, inq
    integer, intent(in) :: jstepx, isign

    integer :: ja, jb, jc, jd, je, jf, jg, jh, ji
    integer :: k, kk, ll, jjj, j, nu, mu
    integer :: jstep, jstepl, laincl
    real(dp) :: c1
    complex(dp) :: zja, zjb, zjc, zjd, zje, zjf, zjg, zjh, zji
    complex(dp) :: z1, z2, z3
    complex(dp) :: exp1, exp2

    jstep = n / (3 * la)
    jstepl = jstep - n
    laincl = la * inq - n
    mu = mod(inq, 3)
    if (isign == -1) mu = 3 - mu
    c1 = drt3 * dhalf
    if (mu == 2) c1 = -c1

    kk = la

    !     loop on nonzero k
    !     -----------------
    do k = inq, jstep - inq, inq

       if (isign == 1) then
          exp1 = trigs(kk)
          exp2 = trigs(2 * kk)
       else
          exp1 = conjg(trigs(kk))
          exp2 = conjg(trigs(2 * kk))
       end if

       !     double loop along first transform in block
       !     ------------------------------------------
       do ll = k, (la - 1) * inq, 3 * jstep

          do jjj = ll, n - 1, 3 * la * inq
             ja = jjj

             !     "transverse" loop
             !     -----------------
             do nu = 1, inq
                jb = ja + jstepl
                if (jb < 0) jb = jb + n
                jc = jb + jstepl
                if (jc < 0) jc = jc + n
                jd = ja + laincl
                if (jd < 0) jd = jd + n
                je = jd + jstepl
                if (je < 0) je = je + n
                jf = je + jstepl
                if (jf < 0) jf = jf + n
                jg = jd + laincl
                if (jg < 0) jg = jg + n
                jh = jg + jstepl
                if (jh < 0) jh = jh + n
                ji = jh + jstepl
                if (ji < 0) ji = ji + n

                !     loop across transforms
                !     ----------------------
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

                   z1 = zjb + zjc
                   z2 = zja - dhalf * z1
                   z3 = c1 * (zjb - zjc)
                   z(ja + j) = zja + z1
                   z(jd + j) = exp1 * (z2 + zi * z3)
                   z(jg + j) = exp2 * (z2 - zi * z3)
                   !----------------------
                   z1 = zje + zjf
                   z2 = zjd - dhalf * z1
                   z3 = c1 * (zje - zjf)
                   z(jb + j) = zjd + z1
                   z(je + j) = exp1 * (z2 + zi * z3)
                   z(jh + j) = exp2 * (z2 - zi * z3)
                   !----------------------
                   z1 = zjh + zji
                   z2 = zjg - dhalf * z1
                   z3 = c1 * (zjh - zji)
                   z(jc + j) = zjg + z1
                   z(jf + j) = exp1 * (z2 + zi * z3)
                   z(ji + j) = exp2 * (z2 - zi * z3)
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

    return

  end subroutine rad3IItwid

end subroutine zgpfa3f
