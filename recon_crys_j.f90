Subroutine recon_crys_j(j_corr,ham,crys,blop,psi,psj,nwfn)

  Use all_to_all_module
  Include 'use.h'
  Implicit None             ! implicit? Just say no!

  Include 'interface.h'
  Include 'all_to_all.h'
  Include 'flibcalls.ph' 
  !
  !     INPUT:
  !     -----
  Type (hamiltonian)     :: ham           ! hamiltonian at a given k-point
  Type (crystal)         :: crys          ! The crystal structure
  Type (blochl_operator) :: blop          ! blochl op. at a given k-point
  Complex(dp)            :: psi(*),psj(*) ! input wave functions
  Integer                :: nwfn          ! number of wavefunctions
  !
  !     OUTPUT:
  !     ------

  Complex(dp) :: j_corr(crys%ntype,crys%mxdatm,3)

  !
  !     WORK:
  !     -----

  Complex(dp), Allocatable :: worki(:),workj(:)

  !
  !     ------------- local variables --------------------------------
  !
  Integer :: len,nanl,indx,k,l,m,mp,kk,nw,np
  Integer :: Index(nwfn,blop%nlprjmx,-4:4,crys%ntype,crys%mxdatm)
  Integer :: indxi,indxj,indxjp,indxjm,i,j,npi,npj
  Real(dp):: fac
  Complex(dp),Parameter :: dcone=(1.d0,0.d0),dczero=(0.d0,0.d0)
  Complex(dp),Parameter :: dci  =(0.d0,1.d0)

  nanl=blop%nanl

  len=ham%gspace%length    

  !     Allocate the work space

  Allocate(worki(nanl*nwfn),workj(nanl*nwfn))

  j_corr = (0.d0,0.d0)

  If(nanl.Gt.0) Then

     !     Do matrix-matrix multiply - calculate the matrix elements

!     Call mzgemm('C','N',nanl,nwfn,len,dcone,blop%opnloc%data(1,1,1),&
!          len,psi(1),len,dczero,worki(1),nanl)
!     Call mzgemm('C','N',nanl,nwfn,len,dcone,blop%opnloc_q%data(1,1,1),&
!          len,psj(1),len,dczero,workj(1),nanl)

     Call fast_all_sum_all_alloc_dc(nanl*nwfn)
     Call fast_all_sum_all_complex(worki(1),nanl*nwfn)
     Call fast_all_sum_all_complex(workj(1),nanl*nwfn)

     Index=0 !we need this for the SR2201 otherwise when we take
             !uninitialed values from index (but don't use them!)
             !the hitachi complains! jry
     !     Set up the indexing
     indx=0
     Do nw=1,nwfn
        Do k=1,crys%ntype
           Do np=1,blop%nlprjmx 
              If(blop%nkb(np,k) .Ne. 0) Then
                 Do m=-blop%lo(np,k),+blop%lo(np,k)
                    Do kk=1,crys%natom(k)
                       indx = indx + 1
                       Index(nw,np,m,k,kk)=indx
                    End Do
                 End Do
              Endif
           End Do
        End Do
     End Do

     !     Evaluate the current-like correction term

     !     A factor of -i is required --- cf. recon_j

     Do nw=1,nwfn
        Do k=1,crys%ntype
           Do kk=1,crys%natom(k)
              Do l=0,3
                 Do i=1,blop%numl(l,k)
                    npi=blop%ln2np(l,i,k)
                    Do j=1,blop%numl(l,k)
                       npj=blop%ln2np(l,j,k)
                       Do m=-l,+l
                          Do mp=-l,+l


                             indxi =Index(nw,npi,m,k,kk)
                             indxj =Index(nw,npj,m,k,kk)
                             indxjp=Index(nw,npj,m+1,k,kk)
                             indxjm=Index(nw,npj,m-1,k,kk)
                             !
                             !     Lx
                             !
                             If(m.Eq.mp+1) Then
                                fac = Sqrt(Real(l*(l+1)-mp*(mp+1)))/2.d0
                                j_corr(k,kk,1)=j_corr(k,kk,1)-    &
!                                     blop%deltaj(l,i,j,k)*fac*dci*&
                                     Conjg(worki(indxi))*workj(indxjm)
                             End If

                             If(m.Eq.mp-1) Then
                                fac = Sqrt(Real(l*(l+1)-mp*(mp-1)))/2.d0
                                j_corr(k,kk,1)=j_corr(k,kk,1)-    &
!                                     blop%deltaj(l,i,j,k)*fac*dci*&
                                     Conjg(worki(indxi))*workj(indxjp)
                             End If
                             !     
                             !     Ly
                             !     
                             If(m.Eq.mp+1) Then
                                fac = Sqrt(Real(l*(l+1)-mp*(mp+1)))/2.d0
                                j_corr(k,kk,2)=j_corr(k,kk,2)-&
!                                     blop%deltaj(l,i,j,k)*fac*&
                                     Conjg(worki(indxi))*workj(indxjm)
                             End If

                             If(m.Eq.mp-1) Then
                                fac =-Sqrt(Real(l*(l+1)-mp*(mp-1)))/2.d0
                                j_corr(k,kk,2)=j_corr(k,kk,2)-&
!                                     blop%deltaj(l,i,j,k)*fac*&
                                     Conjg(worki(indxi))*workj(indxjp)
                             End If
                             !     
                             !     Lz
                             !     
                             If(m.Eq.mp) Then
                                fac = Real(m)
                                j_corr(k,kk,3)=j_corr(k,kk,3)-&
!                                     blop%deltaj(l,i,j,k)*fac*dci*&
                                     Conjg(worki(indxi))*workj(indxj)

                             End If


                          End Do
                       End Do
                    End Do
                 End Do
              End Do
           End Do
        End Do
     End Do

  Endif

  Deallocate(worki,workj)

  Return

End Subroutine recon_crys_j




