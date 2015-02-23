Subroutine apply_rixdvnl(dir,gs,psi,nanl,nloc_deriv,dpsi,crys,pspot) 

  Use all_to_all_module

  Include 'use.h'
  Implicit None             ! implicit? Just say no!
  Include 'interface.h'
  Include 'all_to_all.h'
  Include 'flibcalls.ph'

  !
  !     INPUT:
  !     -----

  Type (parallel_gspace)  :: gs 
  Type (pseudo_potential) :: pspot
  Type (crystal)          :: crys

  Integer :: nanl ! Number of nonlocal projectors
  Integer :: dir  ! Direction into which to take the derivative

  Complex(dp) :: psi(gs%length) ! The wavefn to which DVnl is applied
  Complex(dp) :: nloc_deriv(gs%length+nanl,nanl,2,3) ! Deriv of nl pot at k
  !
  !     OUTPUT:
  !     ------
  !
  Complex(dp) :: dpsi(gs%length) ! R_I x DVnl  * |psi>

  !
  !     DESCRIPTION:
  !     -----------
  !
  !     * Applies the nonlocal potential derivative operator to the wave
  !       function 
  !
  !     * Computes derivative wrt spatial direction: dir
  !     
  !     1999 Chris Pickard
  !
  !
  !     ----------- Local variables -----------------------------------
  !

  Complex(dp), Allocatable :: work(:),phi(:)
  Real(dp) :: Ri(3)
  Integer  :: i,j,k,indx(-1:+1)
  !
  !
  Select Case(dir)
  Case(1)
     indx(+1)=2 ; indx(-1)=3
  Case(2)
     indx(+1)=3 ; indx(-1)=1
  Case(3)
     indx(+1)=1 ; indx(-1)=2
  Case default
     Write(9,*) ' apply_dvnl: illegal direction: ',dir
     Stop
  End Select

  !
  !     --------- Apply the nonlocal part of the hamiltonian ----------
  !

  dpsi=zzero
 
  If(nanl.Gt.0) Then

     Allocate(work(nanl),phi(gs%length))

     Do i=1,-1,-2

        phi = zzero

        !
        !        the +q part
        !
        Call mzgemm('C', 'N', nanl,1,gs%length,zone,&
             nloc_deriv(1,1,2,indx(-i)),gs%length+nanl,&
             psi(1),gs%length,zzero,work(1),nanl)

        Do k=1,nanl

           Ri=proj2rion(k,pspot,crys,gs%fftsize(1:3))

           work(k)=work(k)*Real(nloc_deriv(gs%length+k,1,1,1),dp)*&
                Real(i,dp)*Ri(indx(i))

        End Do

        Call all_sum_all(work,nanl)


        !     
        !     linearly combine the projectors to get +q part
        !     
        Call mzgemm('N','N',gs%length,1,nanl,zone,&
             nloc_deriv(1,1,2,indx(-i)),gs%length+nanl,&
             work(1),nanl,zzero,phi(1),gs%length)
        !     
        !     the -q part
        !     
        Call mzgemm('C', 'N', nanl,1,gs%length,zone,&
             nloc_deriv(1,1,1,indx(-i)), gs%length+nanl, &
             psi(1),gs%length,zzero,work(1),nanl)

        Do k=1,nanl
           Ri=proj2rion(k,pspot,crys,gs%fftsize(1:3))

           work(k)=work(k)*Real(nloc_deriv(gs%length+k,1,1,1),dp)*&
                Real(i,dp)*Ri(indx(i))

        End Do

        Call all_sum_all(work,nanl)


        !     
        !     linearly combine the projectors to get -q part
        !     
        Call mzgemm('N','N',gs%length,1,nanl,zmone,&
             nloc_deriv(1,1,1,indx(-i)),gs%length+nanl,work(1),nanl,&
             zone,phi(1),gs%length) 

        dpsi = dpsi + phi

     End Do

     Deallocate(work,phi)

  Endif

  Return

Contains

  !
  !     Function, to obtain the position of the atom for a given projector
  !     (R_I) in a manner consistent with the r-operator elsewhere
  !     i.e. applied as a sawtooth function
  !     

  Function proj2rion(proj,pspot,crys,nn)

    Implicit None             ! implicit? Just say no!

    !
    !     INPUT:
    !     -----  

    Integer                 :: proj,nn(3)

    Type (pseudo_potential) :: pspot
    Type (crystal)          :: crys

    !
    !     OUTPUT:
    !     ------

    Real(dp) :: proj2rion(3)

    !
    !     ----------- Local variables -----------------------------------
    !
    Integer :: ind,k,kk,l,lm,lmmin,lmmax
    Real(dp) :: n(3),no(3),R(3)

    !     ----- Prepare the sawtooth parameters -----

    no = dhalf

    !     ------ Find the position to which the projector is attached ----

    ind = 0
    loop : Do k=1,crys%ntype  
       Do l=1,5 ! Loop over angular momenta l
          If(pspot%nkb(l,k) .Ne. 0) Then
             lmmin = pspot%lo(l,k)*pspot%lo(l,k)+1
             lmmax = (pspot%lo(l,k)+1)*(pspot%lo(l,k)+1)
             Do lm=lmmin,lmmax ! Loop over m quantum number
                Do kk=1,crys%natom(k)
                   ind = ind + 1
                   If(ind.Eq.proj) Then
                      n=crys%rat(:,kk,k)/pi2
                      Exit loop
                   Endif
                End Do
             End Do
          Endif
       End Do
    End Do loop

    !     ----- Generate the "sawtooth" function -----

    Do k=1,3

       If(n(k).Le.no(k)) Then
          R(k) = n(k) 
       Else
          R(k) = n(k)-done
       Endif

    End Do

    !     ----- Convert to cartesian co-ordinates -----

    proj2rion=Matmul(crys%avec,R)      

  End Function proj2rion

End Subroutine apply_rixdvnl
