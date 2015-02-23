Subroutine apply_r2psi(dir,gs,ffts,rinr,psi,rpsi) 

  Use all_to_all_module

  Include 'use.h'
  Implicit None             ! implicit? Just say no!
  Include 'interface.h'
  Include 'all_to_all.h'
  Include 'flibcalls.ph'

  !
  !     INPUT:
  !     -----

  Type (parallel_gspace) :: gs   ! G-space structure
  Type (fft_struc)       :: ffts ! FFT structure

  Integer                :: dir  ! direction or r to apply

  Complex(dp) :: psi(gs%length)      ! the wavefn to which r is applied
  Complex(dp) :: rinr(ffts%r_size,3) ! Position operator
  !
  !     OUTPUT:
  !     ------
  !
  Complex(dp) :: rpsi(gs%length) ! r|psi>
  !
  !     DESCRIPTION:
  !     -----------
  !
  !  *  Applies the position operator (shifted by d) to a wavefunction
  !
  !     (1999) Chris Pickard
  !
  !
  !     ----------- local variables -----------------------------------
  !
  Complex(dp), Allocatable :: psi_r(:)


  If((dir.Lt.1).Or.(dir.Gt.3)) Then
     Write(9,*) ' apply_r2psi: illegal direction: ',dir
     Stop
  Endif

  !     ----- zero the output array -----

  rpsi=cmplx(0.d0,0.d0)

  !
  !     ----- Psi to real space -----
  !  

  Allocate(psi_r(ffts%r_size))

  Call fourier_transform(-1,ffts,gs,psi(1),psi_r(1),1)

  !
  !     ----- Apply r to psi in real space -----
  !

  psi_r = psi_r * rinr(:,dir)

  !
  !     ----- Psi x r to reciprocal space -----
  !

  Call fourier_transform(+1,ffts,gs,rpsi(1),psi_r(1),1)

  ! 
  !     ----- Finish off -----
  !

  Deallocate(psi_r)

  Return

End Subroutine apply_r2psi
