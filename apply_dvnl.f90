Subroutine apply_dvnl(dir,gs,psi,nanl,nloc_deriv,dpsi,crys) 

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

  Integer :: nanl   ! Number of nonlocal projectors
  Integer :: dir    ! Direction into which to take the derivative

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
  Complex(dp) :: dcone, dczero, dcmone
  Parameter(dcone=(1.d0,0.d0),dczero=(0.d0,0.d0),dcmone=(-1.d0,0.d0))
  Complex(dp), Allocatable :: work(:)
  Integer :: k

  If((dir.Lt.1).Or.(dir.Gt.3)) Then
     Write(9,*) ' apply_dvnl: illegal direction: ',dir
     Stop
  End If
  !
  !     --------- Apply the nonlocal part of the hamiltonian ----------
  !

  dpsi=dczero

  If(nanl.Gt.0) Then

     Allocate(work(nanl))   

     !
     !        the +q part
     !
     Call mzgemm('C', 'N', nanl,1,gs%length,dcone,nloc_deriv(1,1,2,dir), &
          gs%length+nanl,psi(1),gs%length,dczero,work(1),nanl)

     Do k=1,nanl       
        work(k)=work(k)*Dble(nloc_deriv(gs%length+k,1,1,1))      
     End Do

     Call all_sum_all(work,nanl)
     !     
     !     linearly combine the projectors to get +q part
     !     
     Call mzgemm('N','N',gs%length,1,nanl,dcone,nloc_deriv(1,1,2,dir),&
          gs%length+nanl,work(1),nanl,dczero,dpsi(1),gs%length)
     !     
     !     the -q part
     !     
     Call mzgemm('C', 'N', nanl,1,gs%length,dcone,nloc_deriv(1,1,1,dir), &
          gs%length+nanl, psi(1),gs%length,dczero,work(1),nanl)

     Do k=1,nanl
        work(k)=work(k)*Dble(nloc_deriv(gs%length+k,1,1,1))
     End Do


     Call all_sum_all(work,nanl)
     !     
     !     linearly combine the projectors to get -q part
     !     
     Call mzgemm('N','N',gs%length,1,nanl,dcmone,nloc_deriv(1,1,1,dir),&
          gs%length+nanl,work(1),nanl,dcone,dpsi(1),gs%length) 

     Deallocate(work)

  Endif

  Return
End Subroutine apply_dvnl
