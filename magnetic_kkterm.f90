!     @process extchk
!
Subroutine magnetic_kkterm(ham, nloc_deriv_k, qmag, rk, pspot, u_k, &
     len, e_k, noccup, pw_params, crys, weight, uktil, outvv, outhh)
  !
  !
  !     1996 Bernd Pfrommer, UCB
  !
  !
  Use all_to_all_module  
  Include 'use.h'  
  Implicit None            ! implicit? Just say no!
  Include 'interface.h'  
  Include 'all_to_all.h'  
  !
  !     INPUT:
  !     -----
  !
  Integer, Intent(in) :: &
       len, &             ! length of gspace. SGI compiler needs that
       noccup             ! number of occupied states
  Real(dp), Intent(in) :: &
       rk(3), &           ! the kpoint
       qmag, &            ! magnitude of q
       weight, &          ! kpoint weight
       e_k(noccup)        ! eigenvalues of the occupied states

  Type(hamiltonian), Intent(inout) :: ham  
  Type(pw_parameter), Intent(in) :: pw_params  
  Type(pseudo_potential), Intent(in) :: pspot  
  Type(crystal), Intent(in) :: crys  
  Complex(dp), Intent(in) :: &
       u_k(len, noccup)   ! the wavefunctions for that kpoint
  Complex(dp), Intent(in) :: &
       nloc_deriv_k(ham%dim + ham%vnloc%nvecs, &  ! see note below
       ham%vnloc%nvecs, 2, 3)     !  derivative of nonlocal potential at k
  !
  !     OUTPUT:
  !     ------
  !
  Complex(dp), Intent(out) :: &
       uktil(len, noccup, *), &   ! if cgsolveguess, will contain starting
       outvv(3, 3), &             ! with   < dH/dk ><dH/dk>
       outhh(3, 3)                ! with   < dH/dk >< p >
  !
  !     DESCRIPTION:
  !     -----------
  !
  !     computes  Sum_(j=empty) <u_k,i|p1|u_k,j><u_k,j|p0|u_k,i>/(e_k,i-e_
  !     for all occupied states i
  !
  !     To reduce memory usage, the nonlocal derivative array may be
  !     deallocated and recomputed afterwards. So if pw%optimize is set
  !     for memory optimizing, nloc_deriv_k is ignored.
  !
  !
  !
  !     ---------------- local variables ---------------------------------
  !
  Integer :: i, &       ! band index
       pp, &            ! temp operator index
       p1, p0           ! operator index
  Real(dp) :: t0, &
       oldshift         ! original shift of the hamiltonian
  Real(dp), External :: gimmetime  
  Complex(dp), Pointer :: &      ! fortran 90 automatic array
       nloc_deri(:,:,:,:), &
       phi(:,:)                  ! the rhs for cgsolve

  If (Iand(pw_params%output(1), 8) == 8) t0 = gimmetime()  

  Allocate(phi(ham%gspace%length, noccup))  

  oldshift = ham%shift  
  ham%shift = dzero     ! no shift here

  outvv = zzero
  outhh = zzero

  If (Iand(pw_params%optimize, 1) == 1) Then  
     Allocate(nloc_deri(ham%gspace%length + pspot%nanl, pspot%nanl, 2, 3))
     Call setup_nonlocal_derivative(qmag, rk(1), pspot, ham%gspace, &
          nloc_deri(1, 1, 1, 1), crys)
  End If

  Do p0 = 1, 3    ! loop over operator p0
     pp = 1  
     If (Iand(pw_params%optimize, 4) == 4) pp = p0  
     !
     !        compute gradient of uk for rhs of cgsolve
     !
     If (Iand(pw_params%optimize, 1) == 1) Then  
        Do i = 1, noccup  
           Call take_nonloc_deriv(p0, ham%gspace, rk(1), u_k(1, i), &
                ham%vnloc%nvecs, nloc_deri(1, 1, 1, 1), phi(1, i), crys)
        End Do
        Deallocate(nloc_deri)  
     Else  
        Do i = 1, noccup  
           Call take_nonloc_deriv(p0, ham%gspace, rk(1), u_k(1, i),   &
                ham%vnloc%nvecs, nloc_deriv_k(1, 1, 1, 1), phi(1, i), &
                crys)
        End Do
     End If
     !
     !        now solve for uktil
     !

     Call cg_blocksolve(pw_params%output(1), crys, p0, ham,                   &
          pw_params%maxitcgsol, pw_params%epsmag,ham%gspace%rk(1),noccup,  &
          e_k(1),e_k(1),u_k(1, 1), u_k(1, 1), uktil(1, 1, pp),1, phi(1, 1),&
          pw_params%nbandsfft)
    
     If (Iand(pw_params%optimize, 1) == 1) Then  
        Allocate(nloc_deri(ham%gspace%length + pspot%nanl, pspot%nanl, 2, 3))
        Call setup_nonlocal_derivative(qmag, rk(1), &
             pspot, ham%gspace, nloc_deri(1, 1, 1, 1), crys)
     End If


     Do i = 1, noccup       ! loop over all occupied bands
        !
        !           apply operator p1 to u_k~
        !
        Do p1 = 1, 3  
           !
           !              -- vv term: take dH/dk left and right --
           !
           If (Iand(pw_params%optimize, 1) == 1) Then  
              Call take_nonloc_deriv(p1, ham%gspace, ham%gspace%rk(1),      &
                   uktil(1, i, pp), ham%vnloc%nvecs, nloc_deri(1, 1, 1, 1), &
                   phi(1, 1), crys)
           Else  
              Call take_nonloc_deriv(p1, ham%gspace, ham%gspace%rk(1),        &
                   uktil(1, i, pp), ham%vnloc%nvecs, nloc_deriv_k(1, 1, 1, 1),&
                   phi(1, 1), crys)
           End If

           outvv(p0, p1) = outvv(p0, p1) + parallel_zdotc2(ham%gspace%length, &
                u_k(1, i), 1, phi(1, 1), 1) * weight / Real(noccup, dp)
           !
           !              -- hh term: take p  and dH/dk  ----
           !
           Call take_nonloc_deriv(p1, ham%gspace, ham%gspace%rk(1), &
                uktil(1, i, pp), 0, nloc_deriv_k(1, 1, 1, 1), phi(1, 1), crys)

           outhh(p0, p1) = outhh(p0, p1) +                       &
                parallel_zdotc2(ham%gspace%length, u_k(1, i), 1, &
                phi(1, 1), 1) * weight / Real(noccup, dp)

        End Do
     End Do     
  End Do

  If (Iand(pw_params%optimize, 1) == 1) Then  
     Deallocate(nloc_deri)  
  End If

  Deallocate(phi)  

  ham%shift = oldshift  
  If (Iand(pw_params%output(1), 8) == 8) Then  
     Write(9, 100) gimmetime() - t0  
     Call myflush(9)  
  End If

  Return  

100 Format(' TIME FOR KK-TERM:', f12.3)  

End Subroutine magnetic_kkterm
