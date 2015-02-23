Subroutine apply_mom(dir,gs,rk,psi,dpsi,crys) 
 
  Use all_to_all_module

  Include 'use.h'
  Implicit None             ! implicit? Just say no!
  Include 'interface.h'
  Include 'all_to_all.h'
  Include 'flibcalls.ph'

  !
  !     INPUT:
  !     -----

  Type(parallel_gspace) :: gs 
  Type(crystal)         :: crys
  Integer               :: dir  ! Direction into which to take the derivative

  Real(dp)              :: rk(3) ! The k-point
  Complex(dp)           :: psi(gs%length) ! The wavefn

  !
  !     OUTPUT:
  !     ------
  !
  Complex(dp)           :: dpsi(gs%length) ! The P |psi>

  !
  !     DESCRIPTION:
  !     -----------
  !
  !     * Computes derivative wrt spatial direction: dir
  !       We compute data(g)*(k+G), so set the k-point of the gspace
  !       to zero if you just want the gradient!
  !
  !     (1996) Bernd Pfrommer
  !
  !     (1999) Chris Pickard

  !
  !     ----------- Local variables -----------------------------------
  !

  !
  !     ----- Variables for the gspace loop ---------
  !
  Integer :: igv(4),fftn(4),ffth(4),igv3,irod,iord,igs

  Real(dp) :: gv(3)

  If((dir.Lt.1).Or.(dir.Gt.3)) Then
     Write(9,*) ' apply_mom : illegal direction: ',dir
     Stop
  Endif

  dpsi = cmplx(0.d0,0.d0)

  !
  !     --------- Take the gradient ----------
  !
  fftn(1:3) = gs%fftsize(1:3)
  fftn(4)   = gs%fftsize(3)
  ffth(:)   = fftn(:)/2
  igs = 0 
  Do iord=1, gs%lorder    ! Loop through x/y gspace 
     irod=gs%order(1,iord)
     igv(1)= irod/gs%fftsize(2)
     igv(2)= Mod(irod,gs%fftsize(2))
     igv(3)= gs%order(2,iord)
     igv(4)= gs%order(3,iord)
     igv(:)= Mod(igv(:)+ffth(:),fftn(:))-ffth(:) 
     gv(1:2)=Dble(igv(1:2))+rk(1:2)
     Do igv3=igv(3),igv(4)     ! Loop over z axis
        gv(3)=Dble(igv3)+rk(3)
        igs = igs + 1
        dpsi(igs)=dpsi(igs)+psi(igs)*(crys%bvec(dir,1)*gv(1)+&
             crys%bvec(dir,2)*gv(2)+crys%bvec(dir,3)*gv(3))
     End Do
  End Do
  Return

End Subroutine apply_mom 
