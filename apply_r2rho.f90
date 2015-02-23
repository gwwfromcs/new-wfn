Subroutine apply_r2rho(dir,ffts,rinr,rho,rrho) 

  Use all_to_all_module

  Include 'use.h'
  Implicit None             ! implicit? Just say no!
  Include 'interface.h'
  Include 'all_to_all.h'
  Include 'flibcalls.ph'

  !
  !     INPUT:
  !     -----

  Type (fft_struc)       :: ffts ! FFT structure

  Integer                :: dir  ! direction or r to apply

  Complex(dp) :: rho(ffts%r_size)  ! rho to which r is applied
  Complex(dp) :: rinr(ffts%r_size,3) ! Position operator

  !
  !     OUTPUT:
  !     ------
  !
  Complex(dp) :: rrho(ffts%r_size) ! r*rho
  !
  !     DESCRIPTION:
  !     -----------
  !
  !  *  Applies the position operator to rho in realspace
  !
  !     (1999) Chris Pickard
  !
  !
  !     ----------- local variables ---------------------
  !

  Real(dp) :: sgn

  If(Abs(dir).Le.0) Then
     rrho=Cmplx(0.d0,0.d0)
     Return
  Endif

  sgn = dir/Abs(dir)
  dir = Abs(dir)

  If((dir.Lt.1).Or.(dir.Gt.3)) Then
     Write(9,*) ' apply_r2rho: illegal direction: ',dir
     Stop
  Endif

  rrho = sgn*rho*rinr(:,dir)

  Return

End Subroutine apply_r2rho

