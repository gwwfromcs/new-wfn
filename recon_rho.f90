Subroutine recon_rho(rho_corr,ham,crys,blop,psi,nwfn)

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
  Complex(kind=8)        :: psi(*)        ! input wave functions
  Integer                :: nwfn          ! number of wavefunctions
  !
  !     OUTPUT:
  !     ------

  Real(dp) :: rho_corr(crys%ntype,crys%mxdatm)

  !
  !     WORK:
  !     -----

  Complex(dp), Allocatable :: work(:)

  !
  !     ------------- local variables --------------------------------
  !
  Integer :: len,nanl,indx,k,l,m,kk,nw,np,i,j,npi,npj
  Integer :: Index(nwfn,blop%nlprjmx,-3:3,crys%ntype,crys%mxdatm)
  Complex(dp),Parameter :: dcone=(1.d0,0.d0),dczero=(0.d0,0.d0)

  nanl=blop%nanl

  len=ham%gspace%length    

  !     Allocate the work space

  Allocate(work(nanl*nwfn))

  rho_corr = 0.d0

  If(nanl.Gt.0) Then

     !     Do matrix-matrix multiply - calculate the matrix elements

!     Call mzgemm('C','N',nanl,nwfn,len,dcone,blop%opnloc%data(1,1,1),&
!          len,psi(1),len,dczero,work(1),nanl)

     Call fast_all_sum_all_alloc_dc(nanl*nwfn)
     Call fast_all_sum_all_complex(work(1),nanl*nwfn)

     ! *** Evaluate the charge-like correction term

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

     !     Sum the term


     Do nw=1,nwfn
        Do k=1,crys%ntype
           Do l=0,3
              Do i=1,blop%numl(l,k)
                 npi=blop%ln2np(l,i,k)
                 Do j=1,blop%numl(l,k)
                    npj=blop%ln2np(l,j,k)
                    Do m=-l,+l
                       Do kk=1,crys%natom(k)
                          rho_corr(k,kk)=rho_corr(k,kk)+              &
!                               blop%deltarho(l,i,j,k)*                &
                               Real(Conjg(work(Index(nw,npi,m,k,kk)))*&
                               work(Index(nw,npj,m,k,kk)))
                       End Do
                    End Do
                 End Do
              End Do
           End Do
        End Do
     End Do

  Endif

  Deallocate(work)

  Return

End Subroutine recon_rho




