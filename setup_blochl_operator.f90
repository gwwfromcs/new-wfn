!-*-F90-*- 
Subroutine setup_blochl_operator(ham,blop,crys,blop_proj,&
     atom_type,atom_num,lm,l)

  Include 'use.h'
  Implicit None             ! implicit? Just say no!
  Include 'interface.h'

  Type (hamiltonian)     :: ham
  Type (blochl_operator) :: blop
  Type (crystal)         :: crys
  integer, intent(in) 	 :: lm,l  
  Complex(dp) :: blop_proj(*)  !ham%gspace%length) o200 can't have elements of
                               ! structures as array arguments 
  !    
  !     Set up the nonlocal part of the Hamiltonian for a given k-point.
  !     Kleinman-Bylander Pseudopotential.
  !
  !     all it does is setting up the arrays ham%vnloc and ham%xnorm
  !     at the kpoint specified by ham%gspace%rk
  !     

  !     --------------------- local variables ------------------------------

  !
  !     ----- variables for the gspace loop ---------
  !
  Integer :: igv(4),fftn(4),ffth(4),igv3,irod,iord,igs

  Real(dp) :: gv(3)

  Integer :: i,atom_type,atom_num,ni,ni1,lmmin,lmmax 

  Real(dp) :: qi,qk(3),qcar(3),qinv,fi,xni,sum,fpi
  Real(dp) :: zero,one,two,three,half,eps,vq,xi,xdum,xa,xb


  Parameter ( zero = 0.0d0, one = 1.0d0, two = 2.0d0)
  Parameter ( three = 3.0d0, half = 0.5d0)
  Parameter ( eps = 1.0d-8)

  Complex(dp) :: flm(9)
  Complex(dp), Allocatable :: st(:) ! work array for structure factors

  fpi = 16.d0 * Atan(1.d0)
  !
  !     -----------------------------------------------------------------------

  Allocate(st(crys%mxdatm))

  fftn(1:3) = ham%gspace%fftsize(1:3)
  fftn(4)   = ham%gspace%fftsize(3)
  ffth(:)   = fftn(:)/2
  !
  !      starts loop over g-vectors in small gspace 
  !

  igs = 0 
  Do iord=1, ham%gspace%lorder    ! loop through x/y gspace 
     irod=ham%gspace%order(1,iord)
     igv(1)= irod/ ham%gspace%fftsize(2)
     igv(2)= Mod(irod, ham%gspace%fftsize(2))
     igv(3)= ham%gspace%order(2,iord)
     igv(4)= ham%gspace%order(3,iord)
     igv(:)= Mod(igv(:)+ffth(:),fftn(:))-ffth(:) 
     gv(1:2)=Dble(igv(1:2))
     Do igv3=igv(3),igv(4)     ! loop over z axis
        gv(3)=Dble(igv3)
        igs = igs + 1
        qi    = Sqrt(ham%gspace%ekin(igs)) ! get kinetic energy from array
        qk(:) = ham%gspace%rk(:) +gv(:)
        !
        !           cartesian of k+G
        !
        qcar(1) = crys%bvec(1,1)*qk(1)+crys%bvec(1,2)*qk(2)+&
             crys%bvec(1,3)*qk(3)
        qcar(2) = crys%bvec(2,1)*qk(1)+crys%bvec(2,2)*qk(2)+&
             crys%bvec(2,3)*qk(3)
        qcar(3) = crys%bvec(3,1)*qk(1)+crys%bvec(3,2)*qk(2)+&
             crys%bvec(3,3)*qk(3)
        !
        !           ----- compute angular functions
        !     
        !     ** NOTE The projectors are the c.c. of the standard defns.
        !     
        flm(1) = one
        If(qi .Gt. eps) Then
           qinv = one/qi
           flm(2) = +Sqrt(1.5d0)*Cmplx(+qcar(1),+qcar(2))*qinv
           flm(3) = +Sqrt(3.0d0)*Cmplx(+qcar(3),+0.000d0)*qinv
           flm(4) = -Sqrt(1.5d0)*Cmplx(+qcar(1),-qcar(2))*qinv
           qinv   = qinv*qinv
           flm(5) = +Sqrt(15.d0/8.d0)*Cmplx(+qcar(1),+qcar(2))**2*qinv
           flm(6) = +Sqrt(15.d0/2.d0)*Cmplx(qcar(1),+qcar(2))*qcar(3)*qinv
           flm(7) = +Sqrt(5.d0/4.d0)*(3.d0*qcar(3)**2*qinv-1.d0)
           flm(8) = -Sqrt(15.d0/2.d0)*Cmplx(qcar(1),-qcar(2))*qcar(3)*qinv
           flm(9) = +Sqrt(15.d0/8.d0)*Cmplx(+qcar(1),-qcar(2))**2*qinv
        Else
           Do i=2,9
              flm(i) = zero
           End Do
        Endif
        !          
        ! starts loop over second index
        !



        fi = gv(1)*crys%rat(1,atom_num,atom_type) + &
             gv(2)*crys%rat(2,atom_num,atom_type) + &
             gv(3)*crys%rat(3,atom_num,atom_type)
        st(atom_num) = Exp(cmplx(0.d0, fi,dp))
        
        
        ! loop over lo is replaced by loop over l
           ! and l quantum number is given  explicitely by referencing blop%lo
        
        xni = qi/blop%delqnl(atom_type) + two ! interpolate potential
        !
                 ! cubic spline interpolation
                 !
        vq=zero
        ni = xni 
        If(ni .Le. 2) ni = 3
        ni1 = ni+1
        
        If(ni .Lt. blop%nqnl(atom_type)) Then
           xa = Dble(ni1)-xni
           xb = xni-Dble(ni)
           vq = xa*blop%vkb(ni,l,atom_type) + xb*blop%vkb(ni1,l,atom_type) + &
                ((xa**3-xa)*blop%d2vkbdq2(ni-2,l,atom_type) + &
                (xb**3-xb)*blop%d2vkbdq2(ni1-2,l,atom_type))/6.d0
        Endif

        blop_proj(igs) = &
             Conjg(st(atom_num)*flm(lm)*vq)*Sqrt(fpi/crys%vcell)

     End Do                 ! end of loop over 3rd dimension

  End Do                    ! end of loop over small gspace

  Deallocate(st)

  Return

End Subroutine setup_blochl_operator

