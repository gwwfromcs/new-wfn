Subroutine j2chi_mol(crys,ffts,gs,J_r,chi_g,chi_0,g0mask)

  Include 'use.h'
  Implicit None             ! implicit? Just say no!
  Include 'interface.h'
  !
  !     INPUT:
  !     -----     

  Type (crystal)         :: crys          ! for the lattice vectors
  Type (fft_struc)       :: ffts
  Type (parallel_gspace) :: gs            ! potential gspace

  Complex(dp) :: J_r(ffts%r_size,3,3) ! The current in realspace
  Real(dp)    :: chi_0(3,3),g0mask(3,3)

  !
  !     OUTPUT:
  !     ------

  Complex(dp) :: chi_g(gs%length,3,3) ! chi in gspace
  !
  !
  !     1996 Bernd Pfrommer, UCB
  !
  !     1999 Chris Pickard, Modified for molecular case
  !          Paris/Kiel
  !
  !     ---------------- local variables ----------------------------------

  Integer :: i,j,ig,ia,ib,itwoq,iq, iubxq,p0,p1, info

  Real(dp) :: fac, pi4vcell, g(3)

  Complex(dp), Allocatable :: J_g(:,:,:) ! J_r_magnetic in gspace
  Complex(dp) :: ct1(3),ct2(3)
  Allocate(J_g(gs%length,3,3))

  pi4vcell=16.d0*datan(1.d0)/crys%vcell

  !     
  !     Transform to fourier space
  !

  Do ib=1,3
     Do ia=1,3
        Call fourier_transform(1,ffts,gs,J_g(1,ib,ia),J_r(1,ib,ia),1)
     End Do
  End Do

  J_g = J_g*Cmplx(0.d0,1.d0)  ! Factor of i required
  !
  !     Pick up contributions to chi in gspace.
  !

  Do ig=1, gs%length ! Loop over gspace

     chi_g(ig,:,:)=0     

     If(gs%ekin(ig).Gt.0.d-8) Then 

        !
        !     Conversion factor fac:
        !     
        !     4pi     because B = 4pi M
        !     1/vcell for the definition of the Fourier transform
        !
        g = Matmul(crys%bvec(:,:),Dble(gs%gvec(:,ig)))
        fac= pi4vcell/(gs%ekin(ig))

        Do ib=1,3           ! Loop over beta index                 
           ct2(:) = J_g(ig,ib,:)
           Call vecprod_mol(ct1,g,ct2)
           chi_g(ig,:,ib) = ct1(:)*fac
        End Do

     Else
        !
        ! Put in the G=0 component. See Jackson for details.
        !  by default, assume spherical sample, i.e 2/3 of full screening
        !
        chi_g(ig,:,:) = chi_0(:,:)*pi4vcell*g0mask(:,:)

     Endif
  End Do ! End loop over gspace
  !     
  !     
  !
  Deallocate(J_g)

  Return
End Subroutine j2chi_mol

Subroutine vecprod_mol(a,b,c)

  Include 'use.h'
  Implicit None
  !
  !     Computes a = b x c
  !     Type:   "c" "r" "c"

  Complex(dp)  :: a(3),c(3)
  Real(dp) :: b(3)

  a(1)=b(2)*c(3)-b(3)*c(2)
  a(2)=b(3)*c(1)-b(1)*c(3)
  a(3)=b(1)*c(2)-b(2)*c(1)

End Subroutine vecprod_mol

