Subroutine diamagnetic_correction(dia_corr_shift,ham,crys,blop,psi,nwfn)

  Use all_to_all_module
  Include 'use.h'
  Implicit None             ! implicit? Just say no!

  Include 'interface.h'
  Include 'all_to_all.h'
  Include 'flibcalls.ph' 
  !
  !     INPUT:
  !     -----
  Type (hamiltonian),intent(in):: ham        ! hamiltonian at a given k-point
  Type (crystal),intent(in)    :: crys       ! The crystal structure
  Type (blochl_operator)       :: blop       ! blochl op. at a given k-point
  Complex(dp),intent(in)       :: psi(*)     ! input wave functions
  Integer,intent(in)           :: nwfn       ! number of wavefunctions
  !
  !     OUTPUT:
  !     ------

  Real(dp) :: dia_corr_shift(crys%ntype,crys%mxdatm)

  !
  !     WORK:
  !     -----

  Complex(dp), Allocatable :: work(:),work2(:,:)
  Complex(dp), allocatable :: blop_proj(:)
  !
  !     ------------- local variables --------------------------------
  !
  Integer :: len,nanl,indx,atom_type,m,atom_num,nw,np,i,j,npi,npj
  Integer :: Index(nwfn,blop%nlprjmx,-3:3,crys%ntype,crys%mxdatm)
  Complex(dp),Parameter :: dcone=(1.d0,0.d0),dczero=(0.d0,0.d0)
  Integer :: ind      ! counting index for non local operator
  Integer :: lmmin,lmmax,lm,l 

  nanl=blop%nanl

  len=ham%gspace%length    

  !     Allocate the work space

  Allocate(work(nanl*nwfn))
  Allocate(work2(nwfn,nanl))
  Allocate(blop_proj(len))
  dia_corr_shift = 0.d0

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
          len,psi(1),len,dczero,work2(1,ind),1)
     


                 enddo
               enddo
             endif
          enddo
     enddo
     do i=1,nwfn
        do j=1,nanl
           m=(i-1)*nanl+j
           work(m)=work2(i,j)
        enddo
     enddo



     Call fast_all_sum_all_alloc_dc(nanl*nwfn)
     Call fast_all_sum_all_complex(work(1),nanl*nwfn)

     ! *** Evaluate the charge-like correction term

     !     Set up the indexing
     Index=0
     indx=0
     Do nw=1,nwfn
        Do atom_type=1,crys%ntype
           Do np=1,blop%nlprjmx 
              If(blop%nkb(np,atom_type) .Ne. 0) Then
                 Do m=-blop%lo(np,atom_type),+blop%lo(np,atom_type)
                    Do atom_num=1,crys%natom(atom_type)
                       indx = indx + 1
                       Index(nw,np,m,atom_type,atom_num)=indx
                    End Do
                 End Do
              Endif
           End Do
        End Do
     End Do

     !     Sum the term

     Do nw=1,nwfn
        Do atom_type=1,crys%ntype
           Do l=0,3
              Do i=1,blop%numl(l,atom_type)
                 npi=blop%ln2np(l,i,atom_type)
                 Do j=1,blop%numl(l,atom_type)
                    npj=blop%ln2np(l,j,atom_type)
                    Do m=-l,+l
                       Do atom_num=1,crys%natom(atom_type)
                          dia_corr_shift(atom_type,atom_num)=&
                               dia_corr_shift(atom_type,atom_num)+ & 
                               blop%delta_dia(l,i,j,atom_type)*Real(Conjg&
                               (work(Index(nw,npi,m,atom_type,atom_num)))*&
                               work(Index(nw,npj,m,atom_type,atom_num)))
                       End Do
                    End Do
                 End Do
              End Do
           End Do
        End Do
     End Do

  Endif

  Deallocate(work,work2)
  Deallocate(blop_proj)
  Return

End Subroutine diamagnetic_correction




