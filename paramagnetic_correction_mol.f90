Subroutine  paramagnetic_correction_mol (para_corr_shift,ham,crys, &
     blop,psi,psj,nwfn,cfact)

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
  Complex(dp),intent(in) ::cfact          ! 1 or i if we need real or imaginary
                                          !part of para_corr_shift
  !
  !     OUTPUT:
  !     ------
  integer,parameter ::  current=3

  Real(dp) :: para_corr_shift(current,crys%ntype,crys%mxdatm)

  !
  !     WORK:
  !     -----

  Complex(dp), Allocatable :: worki(:),workj(:)
  Complex(dp), Allocatable :: worki2(:,:),workj2(:,:)
  Complex(dp), allocatable :: blop_proj(:)

  !
  !     ------------- local variables --------------------------------
  !
  Integer  :: len,nanl,indx,atom_type,m,mp,atom_num,nw,np
  Integer  :: table(nwfn,blop%nlprjmx,-4:4,crys%ntype,crys%mxdatm)
  Integer  :: indxi,indxj,indxjp,indxjm,i,j,npi,npj
  Real(dp) :: fac
  Integer :: ind      ! counting index for non local operator
  Integer :: lmmin,lmmax,lm,l 
  Complex(dp),Parameter :: dcone=(1.d0,0.d0),dczero=(0.d0,0.d0)
  Complex(dp),Parameter :: dci  =(0.d0,1.d0)

  nanl=blop%nanl

  len=ham%gspace%length    

  !     Allocate the work space

  Allocate(worki(nanl*nwfn),workj(nanl*nwfn))
  Allocate(worki2(nwfn,nanl),workj2(nwfn,nanl))
  Allocate(blop_proj(len))

  para_corr_shift = dzero

  If(nanl.Gt.0) Then

     !     Do matrix-matrix multiply - calculate the matrix elements
     ind=0
     Do atom_type=1,crys%ntype
        Do l=1,blop%nlprjmx              ! loop over angular momenta l
           If(blop%nkb(l,atom_type) .Ne. 0) Then ! found pot of that ang. mom
 
              lmmin = blop%lo(l,atom_type)*blop%lo(l,atom_type) + 1
              lmmax = (blop%lo(l,atom_type)+1)*(blop%lo(l,atom_type)+1)
                 !     
              Do lm=lmmin,lmmax ! loop over m quantum number
                 Do atom_num=1,crys%natom(atom_type)
                    ind = ind + 1
     call setup_blochl_operator(ham,blop,crys,blop_proj,&
          atom_type,atom_num,lm,l)
 

     Call mzgemm('C','N',1,nwfn,len,dcone,blop_proj(1),&
          len,psi(1),len,dczero,worki2(1,ind),1)
     Call mzgemm('C','N',1,nwfn,len,dcone,blop_proj(1),&
          len,psj(1),len,dczero,workj2(1,ind),1)
                 enddo
               enddo
             endif
          enddo
     enddo
     do i=1,nwfn
        do j=1,nanl
           m=(i-1)*nanl+j
           worki(m)=worki2(i,j)
           workj(m)=workj2(i,j)
        enddo
     enddo

     Call fast_all_sum_all_alloc_dc(nanl*nwfn)
     Call fast_all_sum_all_complex(worki(1),nanl*nwfn)
     Call fast_all_sum_all_complex(workj(1),nanl*nwfn)

   table=0 !we need this for the SR2201 otherwise when we take
             !uninitialed values from index (but don't use them!)
             !the hitachi complains! jry
    !     Set up the indexing

     indx=0
     Do nw=1,nwfn
        Do atom_type=1,crys%ntype
           Do np=1,blop%nlprjmx 
              If(blop%nkb(np,atom_type) .Ne. 0) Then
                 Do m=-blop%lo(np,atom_type),+blop%lo(np,atom_type)
                    Do atom_num=1,crys%natom(atom_type)
                       indx = indx + 1
                       table(nw,np,m,atom_type,atom_num)=indx
                    End Do
                 End Do
              Endif
           End Do
        End Do
     End Do

     !     Evaluate the current-like correction term

     Do nw=1,nwfn
        Do atom_type=1,crys%ntype
           Do atom_num=1,crys%natom(atom_type)
              Do l=0,3
                 Do i=1,blop%numl(l,atom_type)
                    npi=blop%ln2np(l,i,atom_type)
                    Do j=1,blop%numl(l,atom_type)
                       npj=blop%ln2np(l,j,atom_type)
                       Do m=-l,+l
                          Do mp=-l,+l
                             
                             indxi =table(nw,npi,m,atom_type,atom_num)
                             indxj =table(nw,npj,m,atom_type,atom_num)
                             indxjp=table(nw,npj,m+1,atom_type,atom_num)
                             indxjm=table(nw,npj,m-1,atom_type,atom_num)

                             !
                             !     Lx
                             !
                             

                             If(m.Eq.mp+1) Then
                                fac = Sqrt(Real(l*(l+1)-mp*(mp+1)))/2.d0
                                para_corr_shift(1,atom_type,atom_num)=&
                                     para_corr_shift(1,atom_type,atom_num)+&
                                     blop%delta_para(l,i,j,atom_type)*&
                                     fac*real(cfact*&
                                     Conjg(worki(indxi))*workj(indxjm),dp)
                             End If


                             If(m.Eq.mp-1) Then
                                fac = Sqrt(Real(l*(l+1)-mp*(mp-1)))/2.d0
                                para_corr_shift(1,atom_type,atom_num)=&
                                     para_corr_shift(1,atom_type,atom_num)+&
                                     blop%delta_para(l,i,j,atom_type)*&
                                     fac*real(cfact*&
                                     Conjg(worki(indxi))*workj(indxjp),dp)
                             End If

                             !     
                             !     Ly
                             !     
                             
                             If(m.Eq.mp+1) Then
                                fac = Sqrt(Real(l*(l+1)-mp*(mp+1)))/2.d0
                                para_corr_shift(2,atom_type,atom_num)=&
                                     para_corr_shift(2,atom_type,atom_num)-&
                                     blop%delta_para(l,i,j,atom_type)*&
                                     fac*real(cfact*dci*&
                                     Conjg(worki(indxi))*workj(indxjm),dp)
                             End If

                             If(m.Eq.mp-1) Then
                                fac =-Sqrt(Real(l*(l+1)-mp*(mp-1)))/2.d0
                                para_corr_shift(2,atom_type,atom_num)=&
                                     para_corr_shift(2,atom_type,atom_num)-&
                                     blop%delta_para(l,i,j,atom_type)*&
                                     fac*real(cfact*dci*&
                                     Conjg(worki(indxi))*workj(indxjp),dp)
                             End If

                             !     
                             !     Lz
                             !     
                             
                             If(m.Eq.mp) Then
                                fac = Real(m)
                                para_corr_shift(3,atom_type,atom_num)=&
                                     para_corr_shift(3,atom_type,atom_num)+&
                                     blop%delta_para(l,i,j,atom_type)*&
                                     fac*real(cfact*&
                                     Conjg(worki(indxi))*workj(indxj),dp)

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
  Deallocate(worki2,workj2)
  Deallocate(blop_proj)
 
  Return

End Subroutine paramagnetic_correction_mol




