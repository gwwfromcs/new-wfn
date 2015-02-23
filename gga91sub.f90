!
!     THE PERDEW-WANG SUBROUTINES AS RECEIVED FROM THE AUTHORS
!
!     variable declaration by Alberto.
!     removed common blocks 1/97 Bernd
!
!
subroutine exch(d, s, u, v, ex, vx)  
  !
  !     gga91 exchange for a spin-unpolarized electronic system
  !     input d : density
  !     input s:  abs(grad d)/(2*kf*d)
  !     input u:  (grad d)*grad(abs(grad d))/(d**2 * (2*kf)**3)
  !     input v: (laplacian d)/(d*(2*kf)**2)
  !     output:  exchange energy per electron (ex) and potential (vx)
  !
  !
  use constants
  implicit none
  !
  !     .. Scalar Arguments ..
  real(dp), intent(in) :: d, s, u, v
  real(dp), intent(out) :: ex, vx
  !     ..
  !     .. Local Scalars ..
  real(dp) :: f, fac, fs, fss, p0, &
       p1, p10, p11, p2, p3, p4, p5, p6, p7, p8, p9, s2, s3, s4
  !     ..
  !     .. Parameters ..
  real(dp), parameter :: a1 = 0.19645d0, a2 = 0.27430d0, a3 = 0.15084d0, &
       a4 = 1.0d2, ax = -0.7385588d0, a = 7.7956d0, b1 = 0.004d0
  !     ..
  fac = ax * d**dthird
  s2 = s * s  
  s3 = s2 * s  
  s4 = s2 * s2  
  p0 = done / sqrt(done + a * a * s2)  
  p1 = log(a * s + done / p0)  
  p2 = exp(-a4 * s2)  
  p3 = done / (done + a1 * s * p1 + b1 * s4)  
  p4 = done + a1 * s * p1 + (a2 - a3 * p2) * s2  
  f = p3 * p4  

  ex = fac * f  
  !     LOCAL EXCHANGE OPTION
  !      EX = FAC
  !     ENERGY DONE. NOW THE POTENTIAL:
  p5 = b1 * s2 - (a2 - a3 * p2)  
  p6 = a1 * s * (p1 + a * s * p0)  
  p7 = dtwo * (a2 - a3 * p2) + dtwo * a3 * a4 * s2 * p2 - dfour * &
       b1 * s2 * f
  fs = p3 * (p3 * p5 * p6 + p7)  
  p8 = dtwo * s * (b1 - a3 * a4 * p2)  
  p9 = a1 * p1 + a * a1 * s * p0 * (dthree - a * a * s2 * p0 * p0)  
  p10 = dfour * a3 * a4 * s * p2 * (dtwo - a4 * s2) - 8.0d0 * b1 * &
       s * f - dfour * b1 * s3 * fs
  p11 = -p3 * p3 * (a1 * p1 + a * a1 * s * p0 + dfour * b1 * s3)  
  fss = p3 * p3 * (p5 * p9 + p6 * p8) + dtwo * p3 * p5 * p6 * p11 + &
       p3 * p10 + p7 * p11

  vx = fac * (dftrd * f - (u - dftrd * s3) * fss - v * fs)  
  !     LOCAL EXCHANGE OPTION:
  !      VX = FAC*THRD4
  return  
  !
end subroutine exch
!
!     =================================================================
!
subroutine corlsd(rs, zet, ec, vcup, vcdn, ecrs, eczet, alfc)  
  !     UNIFORM-GAS CORRELATION OF PERDEW AND WANG 1991
  !     INPUT: SEITZ RADIUS (RS), RELATIVE SPIN POLARIZATION (ZET)
  !     OUTPUT: CORRELATION ENERGY PER ELECTRON (EC), UP- AND DOWN-SPIN
  !     POTENTIALS (VCUP,VCDN), DERIVATIVES OF EC WRT RS (ECRS) & ZET
  !     (ECZET)
  !     OUTPUT: CORRELATION CONTRIBUTION (ALFC) TO THE SPIN STIFFNESS
  !
  use constants
  implicit none
  !
  !     .. Scalar Arguments ..
  real(dp), intent(in) :: rs, zet
  real(dp), intent(out) :: ec, vcdn, vcup, ecrs, eczet, alfc
  !     ..
  !     .. Local Scalars ..
  real(dp) :: alfm, alfrsm, comm, ep, eprs, eu, eurs, f, fz, z4
  !     ..
  !     .. Parameters ..
  real(dp), parameter :: gam =  0.5198421d0, fzz = 1.709921d0
  !     ..
  f = ((done + zet)**dftrd + (done - zet)**dftrd - dtwo) / gam
  call gcor(0.0310907d0, 0.21370d0, 7.5957d0, 3.5876d0, 1.6382d0, &
       0.49294d0, done, rs, eu, eurs)
  call gcor(0.01554535d0, 0.20548d0, 14.1189d0, 6.1977d0, 3.3662d0, &
       0.62517d0, done, rs, ep, eprs)
  call gcor(0.0168869d0, 0.11125d0, 10.357d0, 3.6231d0, 0.88026d0, &
       0.49671d0, done, rs, alfm, alfrsm)
  !     ALFM IS MINUS THE SPIN STIFFNESS ALFC
  alfc = -alfm  
  z4 = zet**4  
  ec = eu * (done - f * z4) + ep * f * z4 - alfm * f * (done - z4) / fzz
  !     ENERGY DONE. NOW THE POTENTIAL:
  ecrs = eurs * (done - f * z4) + eprs * f * z4 - alfrsm * f * &
       (done - z4) / fzz
  fz = dftrd * ((done + zet)**dthird - (done - zet)**dthird) / gam
  eczet = dfour * zet**3 * f * (ep - eu + alfm / fzz) + fz * &
       (z4 * ep - z4 * eu - (done - z4) * alfm / fzz)
  comm = ec - rs * ecrs * dthird - zet * eczet  
  vcup = comm + eczet  
  vcdn = comm - eczet  
  !
  return  
  !
end subroutine corlsd
!
!     ================================================================
!
subroutine corgga(rs, zet, t, uu, vv, ww, h, dvcup, dvcdn, ec, &
     ecrs, eczet, fk, g, sk)
  !
  !     GGA91 CORRELATION
  !     INPUT RS: SEITZ RADIUS
  !     INPUT ZET: RELATIVE SPIN POLARIZATION
  !     INPUT T: ABS(GRAD D)/(D*2.*KS*G)
  !     INPUT UU: (GRAD D)*GRAD(ABS(GRAD D))/(D**2 * (2*KS*G)**3)
  !     INPUT VV: (LAPLACIAN D)/(D * (2*KS*G)**2)
  !     INPUT WW:  (GRAD D)*(GRAD ZET)/(D * (2*KS*G)**2)
  !     OUTPUT H: NONLOCAL PART OF CORRELATION ENERGY PER ELECTRON
  !     OUTPUT DVCUP,DVCDN:  NONLOCAL PARTS OF CORRELATION POTENTIALS
  !
  use constants
  implicit none
  !
  !     .. Scalar Arguments ..
  real(dp), intent(in) :: rs, zet, t, uu, vv, ww, ec, ecrs, eczet, fk, g, sk
  real(dp), intent(out) :: h, dvcdn, dvcup
  real(dp) :: skc5, cc, ccrs, coeff, comm, delt, fac, &
       fact0, fact1, fact2, fact3, fact4, fact5, g3, g4, gz, h0, h0b, &
       h0bt, h0rs, h0rst, h0t, h0tt, h0z, h0zt, h1, h1rs, h1rst, h1t, &
       h1tt, h1z, h1zt, hrs, hrst, ht, htt, hz, hzt, pon, pref, q4, q5, &
       q6, q7, q8, q82, q9, r0, r1, r2, r3, r4, rs2, rs3, rsthrd, t2, t4, t6
  !     ..
  !     ------------- local variables ----------------------------
  !
  real(dp) :: b, b2, bec, bet, bg
  !
  !     .. Parameters ..
  !
  real(dp), parameter :: xnu = 15.75592d0, cc0 = 0.004235d0, &
       cx = -0.001667212d0, alf = 0.09d0, c1 = 0.002568d0, &
       c2 = 0.023266d0, c3 = 7.389d-6, c4 = 8.723d0, c5 = 0.472d0, &
       c6 = 7.389d-2, a4 = 1.0d2
  !     ..
  bet = xnu * cc0  
  delt = dtwo * alf / bet  
  g3 = g**3  
  g4 = g3 * g  
  pon = -delt * ec / (g3 * bet)  
  b = delt / (exp(pon) - done)  
  b2 = b * b  
  t2 = t * t  
  t4 = t2 * t2  
  t6 = t4 * t2  
  rs2 = rs * rs  
  rs3 = rs2 * rs  
  q4 = done + b * t2  
  q5 = done + b * t2 + b2 * t4  
  q6 = c1 + c2 * rs + c3 * rs2  
  q7 = done + c4 * rs + c5 * rs2 + c6 * rs3  
  cc = -cx + q6 / q7  
  r0 = (sk / fk)**2  
  r1 = a4 * r0 * g4  
  coeff = cc - cc0 - dthree * cx / 7.0D0  
  r2 = xnu * coeff * g3  
  r3 = exp (-r1 * t2)  
  h0 = g3 * (bet / delt) * log(done + delt * q4 * t2 / q5)  
  h1 = r3 * r2 * t2  
  h = h0 + h1  
  !  LOCAL CORRELATION OPTION:
  !     H = 0.0D0
  !  ENERGY DONE. NOW THE POTENTIAL:
  ccrs = (c2 + dtwo * c3 * rs) / q7 - q6 * (c4 + dtwo * c5 * rs + &
       dthree * c6 * rs2) / q7**2
  rsthrd = rs * dthird  
  r4 = rsthrd * ccrs / coeff  
  gz = ((done + zet)**dmthird - (done - zet)**dmthird) * dthird  
  fac = delt / b + done
  bg = -dthree * b2 * ec * fac / (bet * g4)  
  bec = b2 * fac / (bet * g3)  
  q8 = q5 * q5 + delt * q4 * q5 * t2  
  q9 = done + dtwo * b * t2  
  q82 = q8 * q8  
  h0b = -bet * g3 * b * t6 * (dtwo + b * t2) / q8  
  h0rs = -rsthrd * h0b * bec * ecrs  
  fact0 = dtwo * delt - dsix * b  
  fact1 = q5 * q9 + q4 * q9 * q9  
  h0bt = dtwo * bet * g3 * t4 * (q4 * q5 * fact0 - delt * fact1) / q82
  h0rst = rsthrd * t2 * h0bt * bec * ecrs  
  h0z = dthree * gz * h0 / g + h0b * (bg * gz + bec * eczet)  
  h0t = dtwo * bet * g3 * q9 / q8  
  h0zt = dthree * gz * h0t / g + h0bt * (bg * gz + bec * eczet)  
  fact2 = q4 * q5 + b * t2 * (q4 * q9 + q5)  
  fact3 = dtwo * b * q5 * q9 + delt * fact2  
  h0tt = dfour * bet * g3 * t * (dtwo * b / q8 - q9 * fact3 / q82)  
  h1rs = r3 * r2 * t2 * (-r4 + r1 * t2 * dthird)  
  fact4 = dtwo - r1 * t2  
  h1rst = r3 * r2 * t2 * (dtwo * r4 * (done - r1 * t2) - dttrd * &
       r1 * t2 * fact4)
  h1z = gz * r3 * r2 * t2 * (dthree - dfour * r1 * t2) / g  
  h1t = dtwo * r3 * r2 * (done - r1 * t2)  
  h1zt = dtwo * gz * r3 * r2 * (dthree - 11.0d0 * r1 * t2 + dfour * &
       r1 * r1 * t4) / g
  h1tt = dfour * r3 * r2 * r1 * t * (dmtwo + r1 * t2)  
  hrs = h0rs + h1rs  
  hrst = h0rst + h1rst  
  ht = h0t + h1t  
  htt = h0tt + h1tt  
  hz = h0z + h1z  
  hzt = h0zt + h1zt  
  comm = h + hrs + hrst + t2 * ht * dsixth + 7.0d0 * t2 * t * htt * dsixth
  pref = hz - gz * t2 * ht / g  
  fact5 = gz * (dtwo * ht + t * htt) / g  
  comm = comm - pref * zet - uu * htt - vv * ht - ww * (hzt - fact5)  
  dvcup = comm + pref  
  dvcdn = comm - pref  
  !  LOCAL CORRELATION OPTION:
  !     DVCUP = 0.0D0
  !     DVCDN = 0.0D0
  !
  return  
  !
end subroutine corgga
!
!     ==============================================================
!
subroutine gcor(a, a1, b1, b2, b3, b4, p, rs, gg, ggrs)
  !
  !     CALLED BY SUBROUTINE CORLSD
  !
  use constants
  implicit none
  !
  !     .. Scalar Arguments ..
  real(dp), intent(in) :: a, a1, b1, b2, b3, b4, p, rs  
  real(dp), intent(out) :: gg, ggrs
  !     ..
  !     .. Local Scalars ..
  real(dp) :: p1, q0, q1, q2, q3, rs12, rs32, rsp  
  !     ..
  p1 = p + done
  q0 = dmtwo * a * (done + a1 * rs)  
  rs12 = sqrt(rs)  
  rs32 = rs12**3  
  rsp = rs**p  
  q1 = dtwo * a * (b1 * rs12 + b2 * rs + b3 * rs32 + b4 * rs * rsp)  
  q2 = log(done + done / q1)  
  gg = q0 * q2  
  q3 = a * (b1 / rs12 + dtwo * b2 + dthree * b3 * rs12 + dtwo * b4 * &
       p1 * rsp)
  ggrs = dmtwo * a * a1 * q2 - q0 * q3 / (q1**2 + q1)  
  !
  return  
  !
end subroutine gcor
