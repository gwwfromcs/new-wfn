Subroutine position_in_rspace(ffts,gs,rinr,crys) 

  Use all_to_all_module

  Include 'use.h'
  Implicit None             ! implicit? Just say no!
  Include 'interface.h'
  Include 'all_to_all.h'
  Include 'flibcalls.ph'

  !
  !     INPUT:
  !     -----

  Type (fft_struc)       :: ffts    ! The FFT structure
  Type (parallel_gspace) :: gs      ! G-space structure
  Type (crystal)         :: crys    ! Crystal structure

  !
  !     OUTPUT:
  !     ------

  Complex(dp) :: rinr(ffts%r_size,3)

  !
  !     LOCAL:
  !     ------

  Integer :: igs,iv(1),dir

  Complex(dp), Allocatable :: ring(:,:)

  !
  !     ========================================================

  Allocate(ring(gs%length,3))

  !     ---- First build in reciprocal space -----

  !     Mostly zero ...

  ring = Cmplx(0.d0,0.d0)

  Do igs=1,gs%length

     !     If contribution not zero

     If (Count(Abs(gs%gvec(:,igs)).Gt.0).Eq.1) Then

        iv          = Maxloc(Abs(gs%gvec(:,igs)))                        
        ring(igs,:) = real((-1)**gs%gvec(iv(1),igs),dp)*crys%avec(:,iv(1))/&
             real(gs%gvec(iv(1),igs),dp)

     Endif

  End Do

  ring = ring*Cmplx(0.d0,1.d0)/pi2

  !     ----- Now transform to real space ----

  Do dir=1,3
     Call fourier_transform(-1,ffts,gs,ring(1,dir),rinr(1,dir),1) 
  End Do

  Deallocate(ring)

  Return

End Subroutine position_in_rspace



