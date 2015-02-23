    subroutine lsdau(imode,iter,k_gspace, energs,pw_params, crys,  &
                     syms, kpoints, bands, wavefn,iscfconv,maxiter)

! LSDA+U method by P. Zhang, 2003-2004 while at UC Berkeley

  use lsdau_shared

  include 'use.h'
  implicit none               ! never remove this line.

  !
  !     INPUT:
  !     -----
  !
  integer, intent(in) :: imode            ! operation mode
  integer, intent(in) :: iter,maxiter     ! SCF iteration number
  logical, intent(in) :: iscfconv
  type(symmetry), intent(in) :: syms      ! symmetry operations
  type(crystal) :: crys       ! crystal structure
  type(pw_parameter), intent(inout) :: &
       pw_params                          ! plane wave parameters
  type(kpoint), intent(in) :: kpoints     ! BZ integration kpoints
                                          ! false otherwise
  type(parallel_gspace), intent(in) :: k_gspace(*)
  !
  !     INPUT/OUTPUT:
  !     -------------
  !
  type(energy), intent(inout) :: energs   ! energies and other information
  type(band), intent(inout) :: bands      ! eigenvalues, occupation numbers etc
  type(complex_gspace_array), intent(inout) ::  wavefn 
                            ! wave functions for all kpoints,


! local variables
  type(energy) :: LSDAenergs              ! energies of the previous iteration
  integer::i,j,k,ii,jj,kk,ll,ia,im,jm,ib,jb,is,irk,idr,nkk

  integer,allocatable::lproj_tem2(:,:)

 
  double precision::rc,tau(3),kkx,kky,kkz

  double complex::cnk(5),drcnk(5,3)
  double complex::cnk_f(5,48)


interface

     subroutine construct_dm0(syms,cnk,cnk_f,dmat,Nsite,lproj, &
                              nlmax,ndim,nband0,nrk,bands,myproc,ipr,iAFM,iNM,iFM,reduced)


     include 'use.h'
     implicit none                       ! implicit? No!
     include 'interface.h'

     type(band), intent(in) :: bands
     type(symmetry), intent(in) :: syms

     integer,intent(in)::Nsite,ndim,nrk,nband0,lproj(Nsite),nlmax,myproc,ipr,iAFM,iNM,iFM
     logical::reduced
     double complex,intent(in)::cnk(nlmax,ndim,nrk,2,Nsite)
     double complex,intent(in)::cnk_f(nlmax,ndim,nrk,2,Nsite,syms%nfrac)

     double complex,intent(out)::dmat(nlmax,nlmax,2,Nsite)
     end subroutine

     subroutine sub_diag(bands,wavefn,k_gspace,Hmat,deltaE,deltaEU,Eband_diff,nrk,ndim, &
                         nband0,myproc,imix)


     include 'use.h'
     implicit none                       ! implicit? No!
     include 'interface.h'

     integer::nband0,ndim,myproc,nrk,imix

     type(band) :: bands
     type(complex_gspace_array), intent(inout) :: wavefn
     type(parallel_gspace), intent(in) :: k_gspace(*)

     double precision::deltaE(ndim,nrk,2),deltaEU(ndim,ndim,nrk,2),Eband_diff
     double complex::Hmat(ndim, ndim,nrk,2)

     end subroutine

     subroutine construct_cnk_LSDAU(cnk_LSDAU,cnk_LSDA,cnk_LSDAU_f, &
                    cnk_LSDA_f,Hmat,occup_new,ndim,nrk,Nsite,lproj,nlmax,syms,reduced,nspin)


     include 'use.h'

     implicit none                       ! implicit? No!

     type(symmetry), intent(in) :: syms

     integer,intent(in)::nrk,Nsite,ndim,lproj(Nsite),nlmax,nspin
     logical::reduced
     double complex,intent(out):: cnk_LSDAU(nlmax,ndim,nrk,2,Nsite)
     double complex,intent(in) :: cnk_LSDA(nlmax,ndim,nrk,2,Nsite)
     double complex,intent(out):: cnk_LSDAU_f(nlmax,ndim,nrk,2,Nsite,syms%nfrac)
     double complex,intent(in) :: cnk_LSDA_f(nlmax,ndim,nrk,2,Nsite,syms%nfrac)

     double complex,intent(in) :: Hmat(ndim,ndim,nrk,2)
     double precision::occup_new(ndim,nrk,2,Nsite)
     end subroutine


end interface


  idr=0
  if(iscfconv.or.(iter.eq.maxiter)) idr=1
  if(imode.eq.0) idr=0

  bands%nband1=nband1
  bands%nband0=nband0


!     if(iAFM.eq.1) then

!        if(myproc.eq.0) then 
!          write(9,*) "------------------"
!          write(9,*) "AFM symmetrization"
!          write(9,*) "------------------"
!        end if

! AFM wavefunction symmetrization
! not working 
!        tau(:)=crys%rat(:,1,iup)-crys%rat(:,1,idown)
!        call AFM_symm(bands,wavefn,k_gspace,tau)

! eigenvalue
      
!        do irk=1,kpoints%nrk
!           do ib=1,nband1
!              bands%energy(ib,irk,1)=(bands%energy(ib,irk,1)+bands%energy(ib,irk,2))/2.d0
!              bands%energy(ib,irk,2)=bands%energy(ib,irk,1)
!           end do        
!        end do
!        energs%vxc0(1)=(energs%vxc0(2)+energs%vxc0(1))/2.d0
!        energs%vxc0(2)=energs%vxc0(1)
!
!
!     end if
!-----------------------------------

     if (((iter.ge.miter).and.(imode.eq.1)).or.(imode.eq.0)) then

         if (myproc.eq.0) then
            write(240,*) 
            write(240,*) "========================================"
            write(240,*) "        iteration",iter
            write(240,*) "========================================"
            call myflush(240)
         end if


        write(9,*) "idr",idr
        call myflush(9)

        call angular_wfnpz_dr(wavefn,k_gspace,crys,kpoints,idr)


 if(imode.eq.1) then

     write(9,*) "----------------"
     write(9,*) "LSDA eigenvalues"
     write(9,*) "----------------"

     if(pw_params%occupy.ne.5) then

     call flevel(imode, 2, pw_params, energs, crys, bands, kpoints)

! save LSDA band occupation
     if(iLSDAU.eq.0) then 
        bandoccupLSDA(1:nband1,:,:)=bands%occup(1:nband1,:,:)
        LSDAenergs=energs
     end if


     call myflush(9)

     if (myproc.eq.0) then
        write(89,*) "iteration",iter
        write(89,211) "VXC0",energs%vxc0(1),energs%vxc0(2)
     end if

 211 format(a4,2f14.8)

     if((iter.eq.miter).and.(idmat.eq.0)) then
        call construct_dm0(syms,cnk_LSDA,cnk_LSDA_f,dmat_LSDAU,&
             Nsite,lproj,nlmax,ndim,nband0,kpoints%nrk,bands,myproc,0,iAFM,iNM,iFM,kpoints%reduced)
        dmat_old=dmat_LSDAU
        dmat_old2=dmat_LSDAU
     end if
     end if


! construct dmat_LSDA (only for printing out info)

     call construct_dm0(syms,cnk_LSDA,cnk_LSDA_f,dmat_LSDA, &
             Nsite,lproj,nlmax,ndim,nband0,kpoints%nrk,bands,myproc,0,iAFM,iNM,iFM,kpoints%reduced)

     idmat=1
     Etot_corr_old=Etot_corr

!----------------------------------

! mix old and new dmat

     if(iter.gt.niter_m) then
       if (pw_params%nscf) then 
       dmat_LSDAU=dmat_old
       else
       if(myproc.eq.0) write(240,*)  "mixing new and old density matrix"
!       dmat_LSDAU=mixing*dmat_LSDAU+0.6d0*(1.d0-mixing)*dmat_old+0.4d0*(1.d0-mixing)*dmat_old2
       dmat_LSDAU=mixing*dmat_LSDAU+(1.d0-mixing)*dmat_old
       call symm_dmat(syms,dmat_LSDAU, Nsite,lproj,nlmax,kpoints%nrk,crys%nspin,kpoints%reduced)
       end if
     end if

     if(isph.eq.0) then
        call construct_Hm(cnk_LSDA,dmat_LSDAU,Hmat,Etot_corr,ELSDA_corr,kpoints%nrk,ndim,nband0,HU,SJ, &
                       idc,alpha_all,lambda,Nsite,lproj,nlmax,bands,myproc,iperturb)
     else

        call construct_Hm2(cnk_LSDA,dmat_LSDAU,Vmat,Hmat,Etot_corr,ELSDA_corr,kpoints%nrk,ndim,nband0,HU,SJ, &
                       idc,alpha_all,lambda,Nsite,lproj,nlmax,bands,myproc,iperturb)
     end if

! update the U term may improve the convergency 

     if((iter.le.niter_m).and.(N_updateU.ge.1)) then
     do ii=1,N_updateU

        if(myproc.eq.0) write(240,*) "update U term",ii

        call sub_diag(bands,wavefn,k_gspace,Hmat,deltaE,deltaEU,Eband_diff,kpoints%nrk,ndim, &
                      nband0,myproc,0)

        call construct_cnk_LSDAU(cnk_LSDAU,cnk_LSDA,cnk_LSDAU_f,cnk_LSDA_f, &
            Hmat,occup_new,ndim,kpoints%nrk,Nsite,lproj,nlmax,syms,kpoints%reduced,crys%nspin)

        do is=1,crys%nspin
        do irk=1,kpoints%nrk
        do ib=1,ndim
           bands%energy(ib+nband0,irk,is)=bands%energy(ib+nband0,irk,is)+deltaE(ib,irk,is)
        end do
        end do
        end do

        call flevel(imode, 0, pw_params, energs, crys, bands, kpoints)

        call construct_dm0(syms,cnk_LSDAU,cnk_LSDAU_f,dmat_LSDAU, &
             Nsite,lproj,nlmax,ndim,nband0,kpoints%nrk,bands,myproc,1,iAFM,iNM,iFM,kpoints%reduced)

        do is=1,crys%nspin
        do irk=1,kpoints%nrk
           do ib=1,ndim
              bands%energy(ib+nband0,irk,is)=bands%energy(ib+nband0,irk,is)-deltaE(ib,irk,is)
           end do
        end do
        end do

        if(isph.eq.0) then
           call construct_Hm(cnk_LSDA,dmat_LSDAU,Hmat,Etot_corr,ELSDA_corr,kpoints%nrk,ndim,nband0,HU,SJ, &
                       idc,alpha_all,lambda,Nsite,lproj,nlmax,bands,myproc,iperturb)
        else

           call construct_Hm2(cnk_LSDA,dmat_LSDAU,Vmat,Hmat,Etot_corr,ELSDA_corr,kpoints%nrk,ndim,nband0,HU,SJ, &
                       idc,alpha_all,lambda,Nsite,lproj,nlmax,bands,myproc,iperturb)
        end if

     end do
     end if

     call sub_diag(bands,wavefn,k_gspace,Hmat,deltaE,deltaEU,Eband_diff,kpoints%nrk,ndim, &
                      nband0,myproc,iLSDAU)

     call construct_cnk_LSDAU(cnk_LSDAU,cnk_LSDA,cnk_LSDAU_f,cnk_LSDA_f, &
             Hmat,occup_new,ndim,kpoints%nrk,Nsite,lproj,nlmax,syms,kpoints%reduced,crys%nspin)

     do is=1,crys%nspin
     do irk=1,kpoints%nrk
     do ib=1,ndim
        bands%energy(ib+nband0,irk,is)=bands%energy(ib+nband0,irk,is)+deltaE(ib,irk,is)
     end do
     end do
     end do

!----------------------------------

     write(9,*)
     write(9,*) "------------------"
     write(9,*) "LSDA+U eigenvalues"
     write(9,*) "------------------"

!------------------------------------------------
! AFM symmetrization

!     if(iAFM.eq.1) then
!     write(9,*) "-----------------------------"
!     write(9,*) "AFM eigenvalue symmetrization"
!     write(9,*) "-----------------------------"
!
!     do irk=1,kpoints%nrk
!        do ib=1,nband1
!           bands%energy(ib,irk,1)=(bands%energy(ib,irk,1)+bands%energy(ib,irk,2))/2.d0
!           bands%energy(ib,irk,2)=bands%energy(ib,irk,1)
!        end do        
!     end do

!     do irk=1,kpoints%nrk
!     do ib=1,ndim
!        deltaE(ib,irk,1)=(deltaE(ib,irk,1)+deltaE(ib,irk,2))/2.d0
!        deltaE(ib,irk,2)=deltaE(ib,irk,1)
!     end do
!     end do
!
!     end if

!------------------------------------------------

     if(pw_params%occupy.eq.5) write(9,*) "occupy levels: from_file"
     if(pw_params%occupy.ne.5) then

     call flevel(imode, 2, pw_params, energs, crys, bands, kpoints)

     if (iand(pw_params%output(1),1073741824 ) == 1073741824) then

     if(myproc == 0) then
     open(109,file="occ.dat.old")
     write(109,*) MIN(bands%min(1),bands%min(2)),bands%nspin,bands%nrk

! write-out occupation for all k-points

     do is = 1,bands%nspin
     do irk=1,kpoints%grid(1)*kpoints%grid(2)*kpoints%grid(3)
           nkk=abs(kpoints%kmap(irk))
           do ib = 1, nband1
           write(109,*) 0.5d0*bands%occup(ib,nkk,is)/kpoints%w(nkk)
          end do
        end do
     end do
     close(109)
     end if
     end if
     end if

     call myflush(9)

!----------------------------------------



!     if(myproc.eq.0) then
!        rewind(244)
!        rewind(245)
!        write(244,211) "VXC0",energs%vxc0(1),energs%vxc0(2)
!        write(245,211) "VXC0",energs%vxc0(1),energs%vxc0(2)
!     end if

!----------------------------------------

     if (pw_params%nscf) then
     write(9,*) "Non-scf calculation, no update of the density matrix."
     else
     call construct_dm0(syms,cnk_LSDAU,cnk_LSDAU_f,dmat_LSDAU, &
             Nsite,lproj,nlmax,ndim,nband0,kpoints%nrk,bands,myproc,1,iAFM,iNM,iFM,kpoints%reduced)
     end if

     if(idr.eq.1) then

! construct dervivatives of the density matrix for force evaluation

     call myflush(9)

     call construct_cnk_LSDAU_dr(drcnk_LSDAU,drcnk_LSDA,&
                    Hmat,ndim,kpoints%nrk,Nsite,nlmax,lproj,crys%nspin)
     call myflush(9)

     call construct_dm0_dr(syms,cnk_LSDAU,drcnk_LSDAU,drdmat_LSDAU, &
             Nsite,lproj,nlmax,ndim,nband0,kpoints%nrk,bands,kpoints%reduced)

     call myflush(9)

     call uforce(HU,SJ,force_U,drdmat_LSDAU,dmat_LSDAU,Vmat,Nsite,lproj,nlmax,crys%nspin)


! print out projected dos
     if (proj_dos.eq.1) then
     allocate(lproj_tem2(1,crys%ntype))
     lproj_tem2(1,:)=lproj_tem(:)
     cnk_LSDAUd2=REAL(DCONJG(cnk_LSDAUd)*cnk_LSDAUd)


     call project_dos(myproc, nproc, cnk_LSDAUd2, syms, crys, bands, energs, kpoints,  &
             nband0, ndim,Nsite,lproj,mproj,lproj_tem2,1,0,nproj, &
             ind_site_atom,ind_site_isite,ind_isite_site)

     deallocate(lproj_tem2)
     end if

     end if

     call print_out(myproc,bands,syms,Hmat,dmat_LSDA,dmat_LSDAU, &
                    dmat_LSDAd,dmat_LSDAUd,deltaE,deltaEU,occup_old,occup_new, &
                    cnk_LSDAU,cnk_LSDAUd,cnk_LSDAU_f,cnk_LSDA, &
                    cnk_LSDAd,ndim,nband0,kpoints,Nsite,lproj,nlmax)

!----------------------------------------


     occup_old=occup_new

     dmat_old2=dmat_old
     dmat_old=dmat_LSDAU

     if(myproc.eq.0) then
       call myflush(89)
       call myflush(240)
!       call myflush(242)
!       call myflush(243)
!       call myflush(244)
!       call myflush(245)
     end if

! if don't do LSDA+U, restore band energies
     if(iLSDAU.eq.0) then


      do is=1,crys%nspin
      do irk=1,kpoints%nrk
         do ib=1,ndim
            bands%energy(ib+nband0,irk,is)=bands%energy(ib+nband0,irk,is)-deltaE(ib,irk,is)
         end do
      end do
      end do

! 

      Etot_corr=0.d0
      Etot_corr_old=0.d0
      Eband_diff=0.d0
      ELSDA_corr=0.d0
      energs=LSDAenergs

! restore LSDA band occupation
     bands%occup(1:nband1,:,:)=bandoccupLSDA(1:nband1,:,:)
     end if
     call myflush(9)


 else ! (if band structure calculation)

!------------band structure calculation begins--------------------

! occupation pertubation not allowed for band structure calculation 
     lambda=0.d0

     if (iLSDAU.eq.1) then

     if(idmat.eq.0) then
        if(myproc.eq.0) write(240,*) "need density matrix for LSDA+U band structure calculation."
        stop
     end if

!     call construct_cm(cnk_LSDA,cmat,Nsite,lproj,nlmax,ndim,kpoints%nrk)
        
     if(isph.eq.0) then
     call construct_Hm(cnk_LSDA,dmat_LSDAU,Hmat,Etot_corr,ELSDA_corr,kpoints%nrk,ndim,nband0,HU,SJ, &
                       idc,alpha_all,lambda,Nsite,lproj,nlmax,bands,myproc,0)
     else

     call construct_Hm2(cnk_LSDA,dmat_LSDAU,Vmat,Hmat,Etot_corr,ELSDA_corr,kpoints%nrk,ndim,nband0,HU,SJ, &
                       idc,alpha_all,lambda,Nsite,lproj,nlmax,bands,myproc,0)
     end if

     call sub_diag(bands,wavefn,k_gspace,Hmat,deltaE,deltaEU,Eband_diff,kpoints%nrk,ndim, &
                   nband0,myproc,1)

! LSDA+U eigenvalues

     do is=1,crys%nspin
     do irk=1,kpoints%nrk
     do ib=1,ndim
        bands%energy(ib+nband0,irk,is)=bands%energy(ib+nband0,irk,is)+deltaE(ib,irk,is)
     end do
     end do
     end do
     end if  !if iLSDAU.eq.1

!------------band structure calculation ends--------------------
     end if  ! if imode.eq.1

     end if  ! if iter .ge. miter

     return
     end subroutine lsdau


!------------------------------------------------

     subroutine lsdau_deallocate()
 
     use lsdau_shared

     if (allocated(ylm2)) deallocate(ylm2)
     if (allocated(ylm1)) deallocate(ylm1)
     if (allocated(sb2)) deallocate(sb2)
     if (allocated(sb1)) deallocate(sb1)
     if (allocated(sb0)) deallocate(sb0)

  deallocate(phase_ikgr0)
  deallocate(ikg)
  deallocate(ntp)
  deallocate(ngrid)
  deallocate(ngrid1)
  deallocate(lambda)


  deallocate(HU)
  deallocate(SJ)
  deallocate(HU_tem)
  deallocate(SJ_tem)
  deallocate(lproj)

! do not deallocate here, used in forstresssym.f90
!  deallocate(lproj_tem)

  deallocate(Rnl_intp)
  deallocate(dmat_LSDAU)
  deallocate(dmat_LSDA)
  deallocate(dmat_LSDAUd)
  deallocate(dmat_LSDAd)
  deallocate(dmat_old)
  deallocate(dmat_old2)
  deallocate(cnk_LSDA)
  deallocate(cnk_LSDA_f)
  deallocate(cnk_LSDAU_f)
  deallocate(cnk_LSDAd)
  deallocate(cnk_LSDAU)
  deallocate(cnk_LSDAUd)
  deallocate(drcnk_LSDAU)
  deallocate(drcnk_LSDA)
  deallocate(drdmat_LSDAU)
  deallocate(cnk_LSDAUd2)
  deallocate(Hmat)
  deallocate(occup_new)
  deallocate(occup_old)
  deallocate(deltaE)
  deallocate(deltaEU)
  deallocate(Vmat)

  deallocate(mproj)
  deallocate(nproj)


  return
  end


!------------------------------------------------------------------------------

     subroutine restore_LSDA(kpoints,k_gspace, bands, wavefn)

  use lsdau_shared
  include 'use.h'
  implicit none               ! never remove this line.


  !     INPUT:
  type(kpoint), intent(in) :: kpoints     ! BZ integration kpoints
  type(parallel_gspace), intent(in) :: k_gspace(*)
  !
  !     INPUT/OUTPUT:
  !     -------------
  !
  type(band), intent(inout) :: bands      ! eigenvalues, occupation numbers etc
  type(complex_gspace_array), intent(inout) ::  wavefn 
                            ! wave functions for all kpoints,

! local variables
  integer::is, irk,ib,ll,jb
  double complex,allocatable::wavefn_tem(:,:)

!
! restore LSDA eigenvalues for next scf iteration
!
       if(iLSDAU.eq.1) then

        write(9,*) "restore LSDA eigenvalues for next iteration"
        do is=1,bands%nspin
           do irk=1,kpoints%nrk
              do ib=nband0+1,nband1
                  bands%energy(ib,irk,is)=bands%energy(ib,irk,is)- deltaE(ib-nband0,irk,is)
              end do
           end do
        end do


     end if


! restore LDA wavefunctions

     if((iwfn.eq.1).and.(iLSDAU.eq.1)) then
        write(9,*) "restore LSDA wavefunctions for next iteration"
        do is=1,bands%nspin
        do irk=1, kpoints%nrk
           ll=k_gspace(irk)%length
           allocate(wavefn_tem(ll,ndim))
           wavefn_tem=(0.d0,0.d0)
 
         do ib=1,ndim
            do jb=1,ndim
            wavefn_tem(:,ib)=wavefn_tem(:,ib)+DCONJG(Hmat(ib,jb,irk,is))* &
                             wavefn%data((jb+nband0-1)*ll+1:(jb+nband0)*ll,irk,is)
            end do
         end do
 
        do jb=1,ndim
           wavefn%data((jb+nband0-1)*ll+1:(jb+nband0)*ll,irk,is)=wavefn_tem(:,jb)
        end do
 
        deallocate(wavefn_tem)
        end do
        end do
     end if

     return
     end subroutine restore_LSDA

!--------------------------------------------------------------------------

  subroutine construct_Hm(cnk_LSDA,dmat,Hmat,Etot_corr,ELSDA_corr,nrk,ndim,nband0, &
                          HU,SJ,idc,alpha_all,lambda,Nsite,lproj,nlmax,bands,myproc,iperturb)

  include 'use.h'  
  implicit none                       ! implicit? No!
  include 'interface.h'  

  type(band), intent(in) :: bands  
  integer::nrk,nband0,ndim,Nsite,myproc,idc,lproj(Nsite),nlmax,iperturb

  double precision::HU(Nsite),SJ(Nsite),Etot_corr,ELSDA_corr,alpha_all   ! correction to total energy 
  double precision::lambda(nlmax,2,Nsite)
  double complex::dmat(nlmax,nlmax,2,Nsite),cnk_LSDA(nlmax,ndim,nrk,2,Nsite)
  double complex::Hmat(ndim,ndim,nrk,2)

  double precision,allocatable::V(:)
  double complex,allocatable::corr(:,:,:,:,:)
  double precision,allocatable::alpha(:)
  double complex,allocatable::Vnm(:,:,:,:)

  double precision,allocatable::avg_occup(:,:),ecorr(:)
  double complex,allocatable::delta_dmat(:,:,:,:),dmat2(:,:,:,:),delta_dmat2(:,:,:,:)

  integer::im,jm,is,irk,ib,ia,jb,nl

  allocate(corr(ndim,ndim,nrk,2,Nsite))
  allocate(V(Nsite))
  allocate(Vnm(nlmax,nlmax,2,Nsite))
  allocate(avg_occup(2,Nsite))
  allocate(delta_dmat(nlmax,nlmax,2,Nsite))
  allocate(delta_dmat2(nlmax,nlmax,2,Nsite))
  allocate(dmat2(nlmax,nlmax,2,Nsite))
  allocate(ecorr(Nsite))

!--------------------------------------------------------------------------
! Petukhov, et al PRB 67, 153106, 2003
  
  allocate(alpha(Nsite))

  do ia=1,Nsite
     V(ia)=-(HU(ia)-SJ(ia))
  end do

!--------------------------------------------------------------------------
    avg_occup=0.d0

    do ia=1,Nsite     ! loop over atomic sites
       nl=2*lproj(ia)+1
       do is=1,bands%nspin
          do im=1,nl
             avg_occup(is,ia)=avg_occup(is,ia)+REAL(dmat(im,im,is,ia))
          end do
       end do   ! spin
       avg_occup(:,ia)=avg_occup(:,ia)/nl
    end do   !over atomic sites
!--------------------------------------------------------------------------
    delta_dmat=dmat

    do ia=1,Nsite
       nl=2*lproj(ia)+1
    do is=1,bands%nspin
    do im=1,nl
       delta_dmat(im,im,is,ia)=delta_dmat(im,im,is,ia)-avg_occup(is,ia)
    end do
    end do
    end do

    do ia=1,Nsite
    nl=2*lproj(ia)+1
    do is=1,bands%nspin
       delta_dmat2(1:nl,1:nl,is,ia)=MATMUL(delta_dmat(1:nl,1:nl,is,ia),DCONJG(TRANSPOSE(delta_dmat(1:nl,1:nl,is,ia))))
    end do
    end do

    if(idc.eq.0) then
      alpha(:)=alpha_all
    else

    alpha=0.0d0

    do ia=1,Nsite
       nl=2*lproj(ia)+1
       do is=1,bands%nspin
          do im=1,nl
          alpha(ia)=alpha(ia)+REAL(delta_dmat2(im,im,is,ia))
          end do
       end do
    end do

   if (bands%nspin.eq.1) avg_occup(2,:)=avg_occup(1,:)
          
    do ia=1,Nsite
       nl=2*lproj(ia)+1
       alpha(ia)=alpha(ia)/(nl*(avg_occup(1,ia)+avg_occup(2,ia) &
                 -avg_occup(1,ia)*avg_occup(1,ia)-avg_occup(2,ia)*avg_occup(2,ia)))
       if(myproc.eq.0) write(240,*) "atom, alpha",ia,alpha(ia)
    end do
    end if

!----------------------------------------------------------------------
! The total energy correction due to eigenvalues changes
! contains double counting which must be taken out
!
  do ia=1,Nsite
     nl=2*lproj(ia)+1
  do is=1,bands%nspin
     dmat2(1:nl,1:nl,is,ia)=MATMUL(dmat(1:nl,1:nl,is,ia),DCONJG(TRANSPOSE(dmat(1:nl,1:nl,is,ia))))
  end do
  end do


  ecorr=0.d0
  do ia=1,Nsite
     nl=2*lproj(ia)+1
     do is=1,bands%nspin
        do im=1,nl
        ecorr(ia)=ecorr(ia)+alpha(ia)*REAL(dmat2(im,im,is,ia))+ &
                     (1.d0-alpha(ia))*REAL(delta_dmat2(im,im,is,ia))

        end do
     end do
  end do

  do ia=1,Nsite
     ecorr(ia)=ecorr(ia)*V(ia)/2.d0
  end do

  ELSDA_corr=0.d0
  do ia=1,Nsite
     nl=2*lproj(ia)+1
     ELSDA_corr=ELSDA_corr-ecorr(ia)+V(ia)*alpha(ia)*(avg_occup(1,ia)+avg_occup(2,ia))*nl/2.d0
  end do

  Etot_corr=0.d0
  do ia=1,Nsite
     Etot_corr=Etot_corr-ecorr(ia)
  end do

!--------------------------------------------------------------------------------------

!U potential matrix
 
  do ia=1,Nsite 
     nl=2*lproj(ia)+1
  do is = 1, bands%nspin  
     vnm(1:nl,1:nl,is,ia)=dmat(1:nl,1:nl,is,ia)

     do im=1,nl
        Vnm(im,im,is,ia)=Vnm(im,im,is,ia)-0.5d0*alpha(ia)-(1.d0-alpha(ia))*avg_occup(is,ia)
     end do

  end do  ! spin
  end do  ! atomic site

 
  corr=0.d0
  do ia=1,Nsite 
     nl=2*lproj(ia)+1
  do is = 1, bands%nspin
     do irk = 1, nrk  
        do ib = 1,ndim

        do jb = 1,ndim

           do im=1,nl
              do jm=1,nl
             
                 corr(ib,jb,irk,is,ia)=corr(ib,jb,irk,is,ia)+Vnm(im,jm,is,ia)* &
                  DCONJG(cnk_LSDA(im,ib,irk,is,ia))*cnk_LSDA(jm,jb,irk,is,ia)
              end do
           end do

        end do
        end do
     end do
  end do  ! spin
  end do  ! atomic site


  do ia=1,Nsite
     corr(:,:,:,:,ia)=corr(:,:,:,:,ia)*V(ia)
  end do

!----------------------
! now apply pertubation 
!----------------------

  if(iperturb.eq.1) then
  do ia=1,Nsite 
     nl=2*lproj(ia)+1
  do is = 1, bands%nspin
     do irk = 1, nrk  
        do ib = 1,ndim
        do jb = 1,ndim
           do im=1,nl
              corr(ib,jb,irk,is,ia)=corr(ib,jb,irk,is,ia)+lambda(im,is,ia)* &
                              DCONJG(cnk_LSDA(im,ib,irk,is,ia))*cnk_LSDA(im,jb,irk,is,ia)
           end do
        end do
        end do
     end do
  end do  ! spin
  end do  ! atomic site
  end if

  Hmat(:,:,:,:)=(0.d0,0.d0)
  do ia=1,Nsite
     Hmat(:,:,:,:)=Hmat(:,:,:,:)+corr(:,:,:,:,ia)
  end do

 120 format(40f9.5)

  do is = 1, bands%nspin
     do irk = 1, nrk
       do ib = 1,ndim
           Hmat(ib,ib,irk,is)=Hmat(ib,ib,irk,is)+bands%energy(ib+nband0, irk, is)
        end do
     end do
  end do
   

 deallocate(corr)
 deallocate(V)
 deallocate(Vnm)
 deallocate(alpha)
 deallocate(avg_occup)
 deallocate(delta_dmat)
 deallocate(delta_dmat2)
 deallocate(dmat2)
 deallocate(ecorr)

 200  format(i5,10f10.5)
 210  format(a5,9a10)

  return
end subroutine construct_Hm


!---------------------------------------------------------------------------

     subroutine sub_diag(bands,wavefn,k_gspace,Hmat,deltaE,deltaEU,Eband_diff,nrk,ndim, &
                         nband0,myproc,imix)


     include 'use.h'
     implicit none                       ! implicit? No!
     include 'interface.h'

!      include 'mpif.h'


     integer::nband0,ndim,myproc,nrk,imix

     type(band) :: bands
     type(complex_gspace_array), intent(inout) :: wavefn 
     type(parallel_gspace), intent(in) :: k_gspace(*)

     double precision::deltaE(ndim,nrk,2),deltaEU(ndim,ndim,nrk,2),Eband_diff
     double complex::Hmat(ndim, ndim,nrk,2)


! local variables  

     integer::LWORK,irk,is,ll,ib,jb,kb,ier,ia
     double precision,allocatable::RWORK(:),Eig(:)
     double complex,allocatable::WORK(:),wavefn_tem(:,:)

! temp variables
     double complex::norm0,norm1
 
     LWORK=5*ndim
     allocate(WORK(LWORK))
     allocate(RWORK(max(3*ndim-2,1)))
     allocate(Eig(ndim))

     do is=1,bands%nspin
     do irk = 1, nrk

        call ZHEEV("V", "U", ndim, Hmat(:,:,irk,is), ndim, Eig, WORK, LWORK, RWORK, ier )

        do ib=1,ndim
           deltaE(ib,irk,is)=Eig(ib)-bands%energy(ib+nband0,irk,is)
        end do

     end do
     end do

! correction to the band energy

! <nk_LSDAU|H_LSDAU-H_LSDA|nk_LSDAU>
! diagonal part

     deltaEU=0.d0

     do is=1,bands%nspin
     do irk = 1, nrk
        do ib=1,ndim
           deltaEU(ib,ib,irk,is)=deltaE(ib,irk,is)+ bands%energy(ib+nband0,irk,is)
        end do
     end do
     end do

     do is=1,bands%nspin
     do irk = 1, nrk
        do ib=1,ndim
        do jb=1,ndim
           deltaEU(ib,ib,irk,is)=deltaEU(ib,ib,irk,is)-REAL(CONJG(Hmat(jb,ib,irk,is))*Hmat(jb,ib,irk,is))* &
                      bands%energy(jb+nband0,irk,is)
        end do
        end do
     end do
     end do

! off-diag

     do is=1,bands%nspin
     do irk = 1, nrk
        do ib=1,ndim
        do jb=1,ndim

           if(ib.ne.jb) then
           do kb=1,ndim
           deltaEU(ib,jb,irk,is)=deltaEU(ib,jb,irk,is)-REAL(CONJG(Hmat(kb,ib,irk,is))*Hmat(kb,jb,irk,is))* &
                      bands%energy(kb+nband0,irk,is)
           end do
           end if

        end do
        end do
     end do
     end do

!-----------------
! LSDA band energy
!-----------------
! There is a definition problem here.
! What do we mean by LSDA band energy, do we mean

! \Sum_occ <nk_LSDAU|H_LSDA|nk_LSDAU> , sum over occupied LSDA+U states

! or

! \Sum_occ <nk_LSDA|H_LSDA|nk_LSDA>


     Eband_diff=0.d0
     do is=1,bands%nspin
     do irk = 1, nrk
        do ib=1,ndim
        do jb=1,ndim
           Eband_diff=Eband_diff -REAL(CONJG(Hmat(jb,ib,irk,is))*Hmat(jb,ib,irk,is))* &
                      bands%energy(jb+nband0,irk,is)*bands%occup(ib+nband0,irk,is)/bands%nspin  
        end do
        end do
     end do
     end do


! <nk_LSDAU|H_LSDAU-H_LSDA|nk_LSDAU>

     do is=1,bands%nspin
     do irk = 1, nrk
        do ib=1,ndim
           Eband_diff=Eband_diff+(deltaE(ib,irk,is)+bands%energy(ib+nband0,irk,is))* &
                      bands%occup(ib+nband0,irk,is)/bands%nspin 
        end do
     end do
     end do

 112 format(i5,3f12.6,i5,40f12.6)
 113 format(a5,3a12,a5,a12)

!---------------------------------------------------------------
! mix the wavefunctions

     if (imix.eq.1) then
     do is=1,bands%nspin
     do irk=1, nrk
        ll=k_gspace(irk)%length
      
        allocate(wavefn_tem(ll,ndim))
        wavefn_tem=(0.d0,0.d0)

        do ib=1,ndim

!-------------------------------------------------
! for debugging purpose, should be deleted later
!           wavefn_tem(:,ib)= wavefn%data((ib+nband0-1)*ll+1:(ib+nband0)*ll,irk,is)
!           wavefn_tem(:,ib)= wavefn_tem(:,ib)* DCONJG(wavefn_tem(:,ib))
!           norm0=REAL(SUM(wavefn_tem(:,ib)))
!           call MPI_ALLREDUCE(norm0,norm1,1,MPI_DOUBLE_COMPLEX,  &
!                              MPI_SUM,MPI_COMM_WORLD,ier)
!           if((abs(REAL(norm1)-1.d0).gt.1.d-13).and.(myproc.eq.0))  &
!              write(*,*) "before",irk,is,ib,REAL(norm1)
! 
!           wavefn_tem(:,ib)=(0.d0,0.d0)
!-------------------------------------------------

           do jb=1,ndim
              wavefn_tem(:,ib)=wavefn_tem(:,ib)+Hmat(jb,ib,irk,is)* &
                            wavefn%data((jb+nband0-1)*ll+1:(jb+nband0)*ll,irk,is)
           end do

       end do

       do ib=1,ndim
          wavefn%data((ib+nband0-1)*ll+1:(ib+nband0)*ll,irk,is)=wavefn_tem(:,ib)

!-------------------------------------------------
! for debugging purpose, should be deleted later
!           wavefn_tem(:,ib)= wavefn_tem(:,ib)* DCONJG(wavefn_tem(:,ib))
!
!           norm0=REAL(SUM(wavefn_tem(:,ib)))
!
!           call MPI_ALLREDUCE(norm0,norm1,1,MPI_DOUBLE_COMPLEX,  &
!                              MPI_SUM,MPI_COMM_WORLD,ier)
!           if((abs(REAL(norm1)-1.d0).gt.1.d-13).and.(myproc.eq.0)) &
!              write(*,*) "after",irk,is,ib,REAL(norm1)

!-------------------------------------------------
       end do


       deallocate(wavefn_tem)
    end do
    end do
    end if

     deallocate(WORK)
     deallocate(RWORK)
     deallocate(Eig)   
     
     return
     end


subroutine findtransd(ntrans,mtrx,qtran)

  include 'use.h'
  implicit none

!
!     OUTPUT
!     vtran(k,i,j) vector transformation matrix
!                  for the k-th symmetry operation
!     qtran(k,i,j) 2nd rank tensor transformation matrix
!                  for the k-th symmetry operation
!     avec(i,j)    i-th comp. of j-th primitive vector
!     bvec(i,j)    i-th comp. of j-th reciprocal vector
!
  real(dp) zero,um,six
  parameter (zero=0.0D0, um=1.0D0, six=6.0D0)
!
  real(dp) qtran(5,5,ntrans)
  real(dp) mtrx(3,3,ntrans)
  real(dp) vtran(3,3,48)
  integer ntrans

! local variables
  real(dp) coef(5,3,3)
  real(dp) rt2i,rt6i,cjmn,delta

  integer i,j,k,m,n
!
  rt2i = um/sqrt(2.d0)
  rt6i = um/sqrt(six)
  delta=1.0D-7

!------------------------------------------------------
  vtran(:,:,1:ntrans)=mtrx(:,:,1:ntrans)

!
!      compose the 2nd rank tensor transformation matrix
!
  coef = zero
  qtran = zero

! Y1=xy
  coef(1,1,2) = rt2i
  coef(1,2,1) = rt2i

! Y2=yz
  coef(2,2,3) = rt2i
  coef(2,3,2) = rt2i

! Y3=zx
  coef(3,1,3) = rt2i
  coef(3,3,1) = rt2i

! Y4=x^2-y^2

  coef(4,1,1) = rt2i
  coef(4,2,2) = -rt2i

! Y5=3z^2-1=2z^2-x^2-y^2

  coef(5,1,1) = -rt6i
  coef(5,2,2) = -rt6i
  coef(5,3,3) = 2.d0*rt6i

  do i=1,5
   do j=1,5
    do k=1,ntrans
     do m=1,3
      do n=1,3
           cjmn = vtran(1,m,k)*(coef(j,1,1)*vtran(1,n,k) &
                             + coef(j,1,2)*vtran(2,n,k) &
                             + coef(j,1,3)*vtran(3,n,k)) &
               + vtran(2,m,k)*(coef(j,2,1)*vtran(1,n,k) &
                             + coef(j,2,2)*vtran(2,n,k) &
                             + coef(j,2,3)*vtran(3,n,k)) &
               + vtran(3,m,k)*(coef(j,3,1)*vtran(1,n,k) &
                             + coef(j,3,2)*vtran(2,n,k) &
                             + coef(j,3,3)*vtran(3,n,k)) 
           qtran(i,j,k) = qtran(i,j,k) + coef(i,m,n)*cjmn
      end do
     end do    
    end do
   end do
  end do
!
!  do k=1,ntrans
!     vtran(:,:,k) = abs(vtran(:,:,k))
!     qtran(:,:,k) = abs(qtran(:,:,k))
!  enddo
  return

end subroutine findtransd

!----------------------------------------------------------------
     subroutine construct_cnk_LSDAU(cnk_LSDAU,cnk_LSDA,cnk_LSDAU_f, &
                    cnk_LSDA_f,Hmat,occup_new,ndim,nrk,Nsite,lproj,nlmax,syms,reduced,nspin)


     include 'use.h'

     implicit none                       ! implicit? No!

     type(symmetry), intent(in) :: syms

     integer,intent(in)::nrk,Nsite,ndim,lproj(Nsite),nlmax,nspin
     logical::reduced
     double complex,intent(out):: cnk_LSDAU(nlmax,ndim,nrk,2,Nsite)
     double complex,intent(in) :: cnk_LSDA(nlmax,ndim,nrk,2,Nsite)
     double complex,intent(out):: cnk_LSDAU_f(nlmax,ndim,nrk,2,Nsite,syms%nfrac)
     double complex,intent(in) :: cnk_LSDA_f(nlmax,ndim,nrk,2,Nsite,syms%nfrac)

     double complex,intent(in) :: Hmat(ndim,ndim,nrk,2)
     double precision::occup_new(ndim,nrk,2,Nsite)
 
     integer::ia,is,irk,ib,jb,im,nl,it
     
     cnk_LSDAU=(0.d0,0.d0)
     do ia=1,Nsite
        nl=2*lproj(ia)+1
        do is=1,nspin
           do irk=1,nrk
              do im=1,nl
              do ib=1,ndim
                 do jb=1,ndim
                 cnk_LSDAU(im,ib,irk,is,ia)=cnk_LSDAU(im,ib,irk,is,ia)+ &
                                            cnk_LSDA(im,jb,irk,is,ia)*Hmat(jb,ib,irk,is)
                 end do
              end do
              end do
           end do
        end do
     end do

!     if(syms%nfrac.ne.0) then
!     cnk_LSDAU_f=(0.d0,0.d0)
!     do it=1,syms%nfrac
!     do ia=1,Nsite
!        nl=2*lproj(ia)+1
!        do is=1,nspin
!           do irk=1,nrk
!              do im=1,nl
!              do ib=1,ndim
!                 do jb=1,ndim
!                 cnk_LSDAU_f(im,ib,irk,is,ia,it)=cnk_LSDAU_f(im,ib,irk,is,ia,it)+ &
!                                            cnk_LSDA_f(im,jb,irk,is,ia,it)*Hmat(jb,ib,irk,is)
!                 end do
!              end do
!              end do
!           end do
!        end do
!     end do
!     end do
!     end if


     occup_new=0.d0

     do ia=1,Nsite     ! loop over atomic sites
        nl=2*lproj(ia)+1
        do is=1,nspin
          do irk = 1, nrk
             do ib = 1,ndim
                do im=1,nl
                   occup_new(ib,irk,is,ia)=occup_new(ib,irk,is,ia)+ &
                      DCONJG(cnk_LSDAU(im,ib,irk,is,ia))*cnk_LSDAU(im,ib,irk,is,ia)
                end do
             end do
          end do
       end do   ! spin
     end do      ! atomic sites

     return
     end 


!----------------------------------------------------------------

     subroutine construct_cnk_LSDAU_dr(drcnk_LSDAU,drcnk_LSDA,&
                    Hmat,ndim,nrk,Nsite,nlmax,lproj,nspin)


     include 'use.h'

     implicit none                       ! implicit? No!

     integer,intent(in)::nrk,Nsite,ndim,lproj(Nsite),nlmax,nspin
     logical::reduced
     double complex,intent(out):: drcnk_LSDAU(nlmax,ndim,nrk,2,Nsite,3)
     double complex,intent(in) :: drcnk_LSDA(nlmax,ndim,nrk,2,Nsite,3)

     double complex,intent(in) :: Hmat(ndim,ndim,nrk,2)
 
     integer::ia,is,irk,ib,jb,im,nl,it,ix
     
     drcnk_LSDAU=(0.d0,0.d0)

     do ix=1,3
     do ia=1,Nsite
        nl=2*lproj(ia)+1
        do is=1,nspin
           do irk=1,nrk
              do im=1,nl
              do ib=1,ndim
                 do jb=1,ndim
                 drcnk_LSDAU(im,ib,irk,is,ia,ix)=drcnk_LSDAU(im,ib,irk,is,ia,ix)+ &
                                    drcnk_LSDA(im,jb,irk,is,ia,ix)*Hmat(jb,ib,irk,is)
                 end do
              end do
              end do
           end do
        end do
     end do
     end do


     return
     end 

!---------------------------------------------------------------------------------
     subroutine print_out(myproc,bands,syms,Hmat,dmat_LSDA,dmat_LSDAU, &
                          dmat_LSDAd,dmat_LSDAUd,deltaE,deltaEU,occup_old,occup_new, &
                         cnk_LSDAU,cnk_LSDAUd,cnk_LSDAU_f,cnk_LSDA,cnk_LSDAd,ndim,nband0,kpoints,Nsite,lproj,nlmax)

     include 'use.h'
     implicit none                       ! implicit? No!
     include 'interface.h'

     type(band):: bands
     type(symmetry) syms
     type(kpoint) kpoints
     integer,intent(IN)::myproc,Nsite,ndim,nband0,Lproj(Nsite),nlmax
     double complex::Hmat(ndim,ndim,kpoints%nrk,2)
     double complex::dmat_LSDA(nlmax,nlmax,2,Nsite),dmat_LSDAU(nlmax,nlmax,2,Nsite)
     double complex::dmat_LSDAd(nlmax,nlmax,2,Nsite),dmat_LSDAUd(nlmax,nlmax,2,Nsite)
     double complex::cnk_LSDA(nlmax,ndim,kpoints%nrk,2,Nsite),cnk_LSDAU(nlmax,ndim,kpoints%nrk,2,Nsite)
     double complex::cnk_LSDAd(nlmax,ndim,kpoints%nrk,2,Nsite),cnk_LSDAUd(nlmax,ndim,kpoints%nrk,2,Nsite)
     double complex::cnk_LSDAU_f(nlmax,ndim,kpoints%nrk,2,Nsite,syms%nfrac)
     double precision::occup_new(ndim,kpoints%nrk,2,Nsite),occup_old(ndim,kpoints%nrk,2,Nsite)
     double precision::deltaE(ndim,kpoints%nrk,2)
     double precision::deltaEU(ndim,ndim,kpoints%nrk,2)

! local variables
     integer::ii,ia,is,irk,ib,jb,indd(1),jj,nrk,nl,im,jm
     double complex,allocatable::dmat(:,:,:,:)

     nrk=kpoints%nrk
     allocate(dmat(nlmax,nlmax,2,Nsite))

! print LSDA d-electron information 


    if(myproc.eq.0) then

    write(240,*)
    write(240,*) "            LSDA density matrix               "
    Write(240,*) "----------------------------------------------"
    end if

    call write_dmat(myproc,dmat_LSDA,dmat_LSDAd,nlmax,lproj,Nsite,1,bands%nspin)

! print LSDA+U d-electron information
! bands%occup is calculated using LSDA eigenvalues, so is wrong

    if(myproc.eq.0) then
    write(240,*)
    write(240,*) "           LSDA+U density matrix              "
    Write(240,*) "----------------------------------------------"
    end if

    call write_dmat(myproc,dmat_LSDAU,dmat_LSDAUd,nlmax,lproj,Nsite,1,bands%nspin)
 
    if(myproc.eq.0) then
    write(240,*) 
    write(240,*) "LSDA+U d-electrons, both occupied and unoccupied"
    Write(240,*) "------------------------------------------------"
    end if

    call construct_dm1(syms,cnk_LSDAU,cnk_LSDAU_f,dmat,Nsite,lproj,nlmax,ndim,nband0,kpoints,bands%nspin)

    call write_dmat(myproc,dmat,dmat_LSDAUd,nlmax,lproj,Nsite,0,bands%nspin)


 111  format(10f16.12)
   
!----------------------------------------------------
! print out LSDA+U energy, occupation of each band
    if(myproc.eq.0) then
     do is=1,bands%nspin
     do irk = 1, nrk

        write(240,*)"is irk", is,irk
        write(240,113)"band","  LSDU eng","  LSDA eng"," max cf","  ib", &
                      (ia,ia=1,Nsite)
!                      ("  new    old  ",ia=1,Nsite)

        do ib=1,ndim
           indd=MAXLOC(abs(Hmat(:,ib,irk,is)))
           write(240,112) ib+nband0,bands%energy(ib+nband0,irk,is), &
                          bands%energy(indd(1)+nband0,irk,is)-deltaE(indd(1),irk,is), &
                          abs(Hmat(indd(1),ib,irk,is)), indd(1)+nband0, &
                        (occup_new(ib,irk,is,ia),ia=1,Nsite)
!                        (occup_new(ib,irk,is,ia),occup_old(ib,irk,is,ia),ia=1,min(Nsite,8))
        end do
     end do
     end do
     end if


 112 format(i4,2f10.6,f8.5,i4,80f7.4)
 113 format(a4,2a10,  a8,  a4,40i7)
 118 format(i5,300f10.6)

!---------------------------------------
    if(myproc.eq.0) then
     open(119,file="delta.dat")

     write(119,*) nband0,nband0+ndim
     do is=1,bands%nspin
     do irk = 1, nrk

        write(119,*)"is irk", is,irk
        do ib=1,ndim
           write(119,118) ib,(deltaEU(ib,jb,irk,is),jb=1,ndim)
        end do
     end do
     end do
     close(119)
     end if
!
!---------------------------------------

     occup_old=occup_new


    if(myproc.eq.0) then
!     write(242,*) "Spherical harmonic"
!     write(242,*)
!     write(243,*) "Spherical harmonic"
!     write(243,*)
!     write(244,*) "Eigenvectors of the density matrix"
!     write(244,*)
!     write(245,*) "Eigenvectors of the density matrix"
!     write(245,*)

!     do is=1,2
!     do irk = 1, nrk
!
!        write(242,*)"is irk", is,irk
!        write(242,115)"band","LSDA eng",  "   Ylm(-lmax),  ...,  Ylm(lmax)"
!
!        write(243,*)"is irk", is,irk
!        write(243,115)"band","LSDAU eng", "   Ylm(-lmax),  ...,  Ylm(lmax)"
!
!        do ib=1,ndim
!           indd=MAXLOC(abs(Hmat(:,ib,irk,is)))
!           write(242,114) ib+nband0,bands%energy(indd(1)+nband0,irk,is)-deltaE(indd(1),irk,is), &
!                   ((REAL(CONJG(cnk_LSDA(jj,indd(1),irk,is,ia))*cnk_LSDA(jj,indd(1),irk,is,ia)), &
!                                    jj=1,2*lproj(ia)+1),ia=1,Nsite)
!
!           write(243,114) ib+nband0,bands%energy(ib+nband0,irk,is), &
!                   ((REAL(CONJG(cnk_LSDAU(jj,ib,irk,is,ia))*cnk_LSDAU(jj,ib,irk,is,ia)), &
!                                    jj=1,2*lproj(ia)+1),ia=1,Nsite)
!        end do
!     end do
!     end do
     end if


! construct angular decomposition onto eigenvectors of the LSDA+U density matrix


     cnk_LSDAUd=(0.d0,0.d0)

     do ia=1,Nsite
        nl=2*lproj(ia)+1
        do is=1,bands%nspin
        do irk=1,nrk
           do ib=1,ndim
           do im=1,nl
           do jm=1,nl
              cnk_LSDAUd(im,ib,irk,is,ia)=cnk_LSDAUd(im,ib,irk,is,ia)+cnk_LSDAU(jm,ib,irk,is,ia)*dmat_LSDAUd(jm,im,is,ia)
           end do
           end do
           end do
        end do
        end do
     end do

!    if(myproc.eq.0) then
!
!     do is=1,2
!     do irk = 1, nrk
!
!        write(244,*)"is irk", is,irk
!        write(244,115)"band","LSDAU eng",  "   phi(-lmax),  ...,  Phi(lmax)"
!
!        do ib=1,ndim
!           write(244,114) ib+nband0,bands%energy(ib+nband0,irk,is), &
!                   ((REAL(CONJG(cnk_LSDAUd(jj,ib,irk,is,ia))*cnk_LSDAUd(jj,ib,irk,is,ia)), &
!                                    jj=1,2*lproj(ia)+1),ia=1,Nsite)
!
!        end do
!     end do
!     end do
!     end if
!
     cnk_LSDAd=(0.d0,0.d0)

     do ia=1,Nsite
        nl=2*lproj(ia)+1
        do is=1,bands%nspin
        do irk=1,nrk
           do ib=1,ndim
           do im=1,nl
           do jm=1,nl
              cnk_LSDAd(im,ib,irk,is,ia)=cnk_LSDAd(im,ib,irk,is,ia)+cnk_LSDA(jm,ib,irk,is,ia)*dmat_LSDAd(jm,im,is,ia)
           end do
           end do
           end do
        end do
        end do
     end do
 
!     if(myproc.eq.0) then
!     do is=1,2
!     do irk = 1, nrk
!
!        write(245,*)"is irk", is,irk
!        write(245,115)"band","LSDA eng",  "   phi(-lmax),  ...,  Phi(lmax)"
!
!        do ib=1,ndim
!           indd=MAXLOC(abs(Hmat(:,ib,irk,is)))
!           write(245,114) ib+nband0,bands%energy(indd(1)+nband0,irk,is)-deltaE(indd(1),irk,is), &
!                   ((REAL(CONJG(cnk_LSDAd(jj,indd(1),irk,is,ia))*cnk_LSDAd(jj,indd(1),irk,is,ia)), &
!                                    jj=1,2*lproj(ia)+1),ia=1,Nsite)
!
!        end do
!     end do
!     end do
!
!     end if

     deallocate(dmat)


 115 format(a4,a9,a45)
 114 format(i4,f9.4,200f7.4)
     return

     end 


!---------------------------------------------------------------------------
! density matrix from occupied states
!
     subroutine construct_dm0(syms,cnk,cnk_f,dmat,Nsite,lproj, & 
                              nlmax,ndim,nband0,nrk,bands,myproc,ipr,iAFM,iNM,iFM,reduced)


     include 'use.h'
     implicit none                       ! implicit? No!
     include 'interface.h'

     type(band), intent(in) :: bands   
     type(symmetry), intent(in) :: syms

     integer,intent(in)::Nsite,ndim,nrk,nband0,lproj(Nsite),nlmax,myproc,ipr,iAFM,iNM,iFM
     logical::reduced
     double complex,intent(in)::cnk(nlmax,ndim,nrk,2,Nsite)
     double complex,intent(in)::cnk_f(nlmax,ndim,nrk,2,Nsite,syms%nfrac)

     double complex,intent(out)::dmat(nlmax,nlmax,2,Nsite)

! local variables
     integer::ii,jj,kk,is,ia,irk,ib,jb,nl,it,ntrans,imap
   

     double complex,allocatable::qtranc(:,:,:),vtranc(:,:,:),dmat_tem(:,:,:,:),dmat_f(:,:,:,:,:)

!--------------------------------------------------------------------------


    dmat=(0.d0,0.d0)

    do ia=1,Nsite     ! loop over atomic sites
       nl=2*lproj(ia)+1 
       do is=1,bands%nspin
          do ii=1,nl
             do jj=1,nl

             do irk = 1, nrk
                do ib = 1,ndim
                   dmat(ii,jj,is,ia)=dmat(ii,jj,is,ia)+cnk(ii,ib,irk,is,ia)* &
                   DCONJG(cnk(jj,ib,irk,is,ia))*bands%occup(ib+nband0,irk,is)*0.5d0
                end do
             end do
          end do 
          end do  

       end do   ! spin
    end do      ! atomic sites

!----------------------------------------------
    ntrans=syms%ntrans

!
! sym operations with fraction translation are disabled in this version
!
!    do ii=1,syms%ntrans
!       if ((SUM(abs(syms%tnp(:,ii))).gt.1.d-6).and.(syms%nfrac.eq.0)) ntrans=ntrans-1
!    end do
!----------------------------------------------


    allocate(dmat_tem(nlmax,nlmax,2,Nsite))
    dmat_tem=(0.d0,0.d0)

    if(reduced) then

!    if (syms%nfrac.ne.0) then
!       write(9,*) "not debugged"
!       call myflush(9)
!       stop
!
!       allocate(dmat_f(nlmax,nlmax,2,Nsite,syms%nfrac))
!       dmat_f=(0.d0,0.d0)
!
!       do it=1,syms%nfrac ! loop over fractional translations
!
!       do ia=1,Nsite     ! loop over atomic sites
!       nl=2*lproj(ia)+1
!       do is=1,bands%nspin
!          do ii=1,nl
!             do jj=1,nl
!
!             do irk = 1, nrk
!                do ib = 1,ndim
!                   dmat_f(ii,jj,is,ia,it)=dmat_f(ii,jj,is,ia,it)+cnk_f(ii,ib,irk,is,ia,it)* &
!                   DCONJG(cnk_f(jj,ib,irk,is,ia,it))*bands%occup(ib+nband0,irk,is)*0.5d0
!                end do
!             end do
!          end do
!          end do
!
!       end do   ! spin
!       end do      ! atomic sites
!       end do
!
!    end if

!---------------------------------------------------------
! symmetrize density matrix

    allocate(qtranc(5,5,syms%ntrans))
    allocate(vtranc(3,3,syms%ntrans))

    call findtransc(syms%ntrans,syms%rsymmat,qtranc,vtranc)


    do ia=1,Nsite
       nl=2*lproj(ia)+1
       do is=1,bands%nspin

          if(lproj(ia).eq.2) then

          do ii=1,syms%ntrans
!             if ((SUM(abs(syms%tnp(:,ii))).gt.1.d-6).and.(syms%nfrac.gt.0)) then
!             write(9,*) "not debugged"
!             call myflush(9)
!             stop
!
!             it=it+1
!             imap=syms%ind_rot(ii,ia)
!             dmat_tem(1:nl,1:nl,is,ia)=dmat_tem(1:nl,1:nl,is,ia)+ &
!             MATMUL(qtranc(:,:,ii), MATMUL(dmat_f(1:nl,1:nl,is,imap,it),DCONJG(TRANSPOSE(qtranc(:,:,ii)))))
!             end if
! 
!             if (SUM(abs(syms%tnp(:,ii))).lt.1.d-6) then
      
             imap=syms%ind_rot(ii,ia)

!             dmat_tem(1:nl,1:nl,is,ia)=dmat_tem(1:nl,1:nl,is,ia)+ &
!             MATMUL(qtranc(:,:,ii), MATMUL(dmat(1:nl,1:nl,is,imap),DCONJG(TRANSPOSE(qtranc(:,:,ii)))))

             dmat_tem(1:nl,1:nl,is,imap)=dmat_tem(1:nl,1:nl,is,imap)+ &
             MATMUL(qtranc(:,:,ii), MATMUL(dmat(1:nl,1:nl,is,ia),DCONJG(TRANSPOSE(qtranc(:,:,ii)))))


!             end if
          end do

!          dmat(1:nl,1:nl,is,ia)=dmat_tem(1:nl,1:nl)
          end if

          if(lproj(ia).eq.1) then
          do ii=1,syms%ntrans

!             if ((SUM(abs(syms%tnp(:,ii))).gt.1.d-6).and.(syms%nfrac.gt.0)) then
!             write(9,*) "not debugged"
!             call myflush(9)
!             stop
!             it=it+1
!             imap=syms%ind_rot(ii,ia)
!             dmat_tem(1:nl,1:nl,is,ia)=dmat_tem(1:nl,1:nl,is,ia)+ &
!             MATMUL(vtranc(:,:,ii), MATMUL(dmat_f(1:nl,1:nl,is,imap,it),DCONJG(TRANSPOSE(vtranc(:,:,ii)))))
!             end if

!             if (SUM(abs(syms%tnp(:,ii))).lt.1.d-6) then


             imap=syms%ind_rot(ii,ia)

!             dmat_tem(1:nl,1:nl,is,ia)=dmat_tem(1:nl,1:nl,is,ia)+ &
!             MATMUL(vtranc(:,:,ii), MATMUL(dmat(1:nl,1:nl,is,imap),DCONJG(TRANSPOSE(vtranc(:,:,ii)))))

             dmat_tem(1:nl,1:nl,is,imap)=dmat_tem(1:nl,1:nl,is,imap)+ &
             MATMUL(vtranc(:,:,ii), MATMUL(dmat(1:nl,1:nl,is,ia),DCONJG(TRANSPOSE(vtranc(:,:,ii)))))



!             end if
          end do
!          dmat(1:nl,1:nl,is,ia)=dmat_tem(1:nl,1:nl)

          end if


          if(lproj(ia).eq.0) then
          dmat_tem(1:nl,1:nl,is,ia)=dmat(1:nl,1:nl,is,ia)*ntrans
          end if

       end do
    end do

!    dmat=dmat/syms%ntrans
    dmat=dmat_tem/ntrans


    deallocate(qtranc)
    deallocate(vtranc)
!    if(syms%nfrac.ne.0) deallocate(dmat_f)


! make diagonal terms real
!    do ia=1,Nsite
!       nl=2*lproj(ia)+1
!       do is=1,bands%nspin
!          do ii=1,nl
!             dmat(ii,ii,is,ia)=REAL(dmat(ii,ii,is,ia))
!          end do
!       end do
!    end do

      


    end if   ! if reduced

! make density matrix hermitian

    do ia=1,Nsite
       nl=2*lproj(ia)+1
       do is=1,bands%nspin
          dmat_tem(:,:,is,ia)=dmat(:,:,is,ia)+DCONJG(TRANSPOSE(dmat(:,:,is,ia)))
       end do
    end do

    dmat=dmat_tem/2.d0

    deallocate(dmat_tem)

!    if ((iAFM.eq.2).or.(iAFM.eq.3).or.(iNM.gt.0).or.(iFM.gt.0)) call lsdau_symm2(dmat)

     if(ipr.eq.1) then
     if(myproc.eq.0) then

     open(109,file="dmat.new")

     do ia=1,Nsite
        nl=2*lproj(ia)+1

        do is=1,bands%nspin
           write(89,*) " site ",ia," spin ",is
           write(109,*) " site ",ia," spin ",is
           do ii=1,nl
              write(89,*) (dmat(ii,jj,is,ia),jj=1,nl)
              write(109,*) (dmat(ii,jj,is,ia),jj=1,nl)
           end do
        end do
     end do


     call myflush(89)
     close(109)

     end if
     end if

    return
    end

!---------------------------------------------------------------------------
! symmetrize density matrix
!
     subroutine symm_dmat(syms,dmat,Nsite,lproj,nlmax,nrk,nspin,reduced)


     include 'use.h'
     implicit none                       ! implicit? No!
     include 'interface.h'

     type(symmetry), intent(in) :: syms

     integer,intent(in)::Nsite,nrk,lproj(Nsite),nlmax,nspin
     logical::reduced

     double complex,intent(out)::dmat(nlmax,nlmax,2,Nsite)

! local variables
     integer::ii,jj,kk,is,ia,irk,ib,jb,nl,it,ntrans,imap
   
     double complex,allocatable::qtranc(:,:,:),vtranc(:,:,:),dmat_tem(:,:,:,:),dmat_f(:,:,:,:,:)

!--------------------------------------------------------------------------

    ntrans=syms%ntrans

    allocate(dmat_tem(nlmax,nlmax,2,Nsite))
    dmat_tem=(0.d0,0.d0)

    if(reduced) then

!---------------------------------------------------------
! symmetrize density matrix

    allocate(qtranc(5,5,syms%ntrans))
    allocate(vtranc(3,3,syms%ntrans))

    call findtransc(syms%ntrans,syms%rsymmat,qtranc,vtranc)


    do ia=1,Nsite
       nl=2*lproj(ia)+1
       do is=1,nspin

          if(lproj(ia).eq.2) then

          do ii=1,syms%ntrans
             imap=syms%ind_rot(ii,ia)
             dmat_tem(1:nl,1:nl,is,imap)=dmat_tem(1:nl,1:nl,is,imap)+ &
             MATMUL(qtranc(:,:,ii), MATMUL(dmat(1:nl,1:nl,is,ia),DCONJG(TRANSPOSE(qtranc(:,:,ii)))))
          end do

          end if

          if(lproj(ia).eq.1) then
          do ii=1,syms%ntrans

             imap=syms%ind_rot(ii,ia)

             dmat_tem(1:nl,1:nl,is,imap)=dmat_tem(1:nl,1:nl,is,imap)+ &
             MATMUL(vtranc(:,:,ii), MATMUL(dmat(1:nl,1:nl,is,ia),DCONJG(TRANSPOSE(vtranc(:,:,ii)))))

          end do

          end if


          if(lproj(ia).eq.0) then
          dmat_tem(1:nl,1:nl,is,ia)=dmat(1:nl,1:nl,is,ia)*ntrans
          end if

       end do
    end do

    dmat=dmat_tem/ntrans


    deallocate(qtranc)
    deallocate(vtranc)

    end if   ! if reduced

! make density matrix hermitian

    do ia=1,Nsite
       nl=2*lproj(ia)+1
       do is=1,nspin
          dmat_tem(:,:,is,ia)=dmat(:,:,is,ia)+DCONJG(TRANSPOSE(dmat(:,:,is,ia)))
       end do
    end do

    dmat=dmat_tem/2.d0

    deallocate(dmat_tem)

    return
    end


!--------------------------------------------------------------


     subroutine construct_dm0_dr(syms,cnk,drcnk,drdmat,Nsite,lproj, & 
                              nlmax,ndim,nband0,nrk,bands,reduced)

     include 'use.h'
     implicit none                       ! implicit? No!
     include 'interface.h'

     type(band), intent(in) :: bands   
     type(symmetry), intent(in) :: syms

     logical::reduced
     integer,intent(in)::Nsite,ndim,nrk,nband0,lproj(Nsite),nlmax
     double complex,intent(in)::cnk(nlmax,ndim,nrk,2,Nsite)
     double complex,intent(in)::drcnk(nlmax,ndim,nrk,2,Nsite,3)

     double complex,intent(out)::drdmat(nlmax,nlmax,2,Nsite,3)

! local variables
     integer::ii,jj,kk,is,ia,irk,ib,jb,nl,it,ntrans,ix,imap

     double complex,allocatable::qtranc(:,:,:),vtranc(:,:,:),dmat_tem(:,:,:,:,:)

!--------------------------------------------------------------------------


    ntrans=syms%ntrans
    allocate(dmat_tem(nlmax,nlmax,2,Nsite,3))

!--------------------------------------------------------------------------------------   
! sym operations with fraction translation are disabled in this version
!   
!    do ii=1,syms%ntrans
!       if ((SUM(abs(syms%tnp(:,ii))).gt.1.d-6).and.(syms%nfrac.eq.0)) ntrans=ntrans-1
!    end do
!--------------------------------------------------------------------------------------

    drdmat=(0.d0,0.d0)

    do ix=1,3
    do ia=1,Nsite     ! loop over atomic sites
       nl=2*lproj(ia)+1 
       do is=1,bands%nspin
          do ii=1,nl
             do jj=1,nl

             do irk = 1, nrk
                do ib = 1,ndim
                   drdmat(ii,jj,is,ia,ix)=drdmat(ii,jj,is,ia,ix)+ &
                   (drcnk(ii,ib,irk,is,ia,ix)*DCONJG(cnk(jj,ib,irk,is,ia))+ &
                    cnk(ii,ib,irk,is,ia)*DCONJG(drcnk(jj,ib,irk,is,ia,ix))) &
                    *bands%occup(ib+nband0,irk,is)*0.5d0
               end do
             end do

          end do 
          end do  

       end do   ! spin
    end do      ! atomic sites
    end do      ! direction x,y,z

!---------------------------------------------------------
! symmetrize density matrix

    if(reduced) then

    allocate(qtranc(5,5,syms%ntrans))
    allocate(vtranc(3,3,syms%ntrans))

    call findtransc(syms%ntrans,syms%rsymmat,qtranc,vtranc)

    dmat_tem=(0.d0,0.d0)
    do ix=1,3
    do ia=1,Nsite
       nl=2*lproj(ia)+1
       do is=1,bands%nspin


          if(lproj(ia).eq.2) then
          do ii=1,syms%ntrans

!          if (SUM(abs(syms%tnp(:,ii))).lt.1.d-6) then
             imap=syms%ind_rot(ii,ia)
!             dmat_tem(1:nl,1:nl,is,ia,ix)=dmat_tem(1:nl,1:nl,is,ia,ix)+ &
!             MATMUL(qtranc(:,:,ii), MATMUL(drdmat(1:nl,1:nl,is,imap,ix),DCONJG(TRANSPOSE(qtranc(:,:,ii)))))

             dmat_tem(1:nl,1:nl,is,imap,ix)=dmat_tem(1:nl,1:nl,is,imap,ix)+ &
             MATMUL(qtranc(:,:,ii), MATMUL(drdmat(1:nl,1:nl,is,ia,ix),DCONJG(TRANSPOSE(qtranc(:,:,ii)))))


!          end if
          end do
!            drdmat(1:nl,1:nl,is,ia,ix)=dmat_tem(1:nl,1:nl)
          end if

          if(lproj(ia).eq.1) then
          do ii=1,syms%ntrans

!          if (SUM(abs(syms%tnp(:,ii))).lt.1.d-6) then
             imap=syms%ind_rot(ii,ia)
!             dmat_tem(1:nl,1:nl,is,ia,ix)=dmat_tem(1:nl,1:nl,is,ia,ix)+ &
!             MATMUL(vtranc(:,:,ii), MATMUL(drdmat(1:nl,1:nl,is,imap,ix),DCONJG(TRANSPOSE(vtranc(:,:,ii)))))

             dmat_tem(1:nl,1:nl,is,imap,ix)=dmat_tem(1:nl,1:nl,is,imap,ix)+ &
             MATMUL(vtranc(:,:,ii), MATMUL(drdmat(1:nl,1:nl,is,ia,ix),DCONJG(TRANSPOSE(vtranc(:,:,ii)))))

!          end if
          end do
!          drdmat(1:nl,1:nl,is,ia,ix)=dmat_tem(1:nl,1:nl)

          end if


          if(lproj(ia).eq.0) then
          dmat_tem(1:nl,1:nl,is,ia,ix)=drdmat(1:nl,1:nl,is,ia,ix)*ntrans
          end if

       end do
    end do
    end do

    drdmat=dmat_tem/ntrans

    deallocate(qtranc)
    deallocate(vtranc)
    end if

! make diagonal terms real
    do ix=1,3
    do ia=1,Nsite
       nl=2*lproj(ia)+1
       do is=1,bands%nspin
          do ii=1,nl
             drdmat(ii,ii,is,ia,ix)=REAL(drdmat(ii,ii,is,ia,ix))
          end do
       end do
    end do
    end do


    deallocate(dmat_tem)

    return
    end


subroutine findtransc(ntrans,vtran,qtranc,vtranc)

  include 'use.h'
  implicit none

!
!     OUTPUT
!     vtran(k,i,j) vector transformation matrix
!                  for the k-th symmetry operation
!     qtran(k,i,j) 2nd rank tensor transformation matrix
!                  for the k-th symmetry operation
!     avec(i,j)    i-th comp. of j-th primitive vector
!     bvec(i,j)    i-th comp. of j-th reciprocal vector
!
  real(dp) zero,um,six
  parameter (zero=0.0D0, um=1.0D0, six=6.0D0)
!
  double complex::qtranc(5,5,ntrans),vtranc(3,3,ntrans)
  real(dp) vtran(3,3,ntrans)
  double complex::A(5,5),Ainv(5,5),ione,B(3,3),Binv(3,3)
  real(dp),allocatable:: qtran(:,:,:)

  integer ntrans

! local variables
  real(dp) coef(5,3,3)
  real(dp) rt2i,rt6i,cjmn

  integer i,j,k,m,n
!
  allocate(qtran(5,5,ntrans))

  rt2i = um/sqrt(2.d0)
  rt6i = um/sqrt(six)
  ione=CMPLX(0.d0,1.d0,kind=8)

  A=(0.d0,0.d0)
  Ainv=(0.d0,0.d0)

  A(1,1)=-ione
  A(1,4)=1.d0
  A(2,2)=-ione
  A(2,3)=1.d0
  A(3,5)=dsqrt(2.d0)
  A(4,2)=-ione
  A(4,3)=-1.d0
  A(5,1)=ione
  A(5,4)=1.d0

  A=A/dsqrt(2.d0)

  Ainv(1,1)=ione
  Ainv(1,5)=-ione
  Ainv(2,2)=ione
  Ainv(2,4)=ione
  Ainv(3,2)=1.d0
  Ainv(3,4)=-1.d0
  Ainv(4,1)=1.d0
  Ainv(4,5)=1.d0
  Ainv(5,3)=dsqrt(2.d0)

  Ainv=Ainv/dsqrt(2.d0)

  B=(0.d0,0.d0)
  B(1,1)=1.d0
  B(1,2)=-ione
  B(2,3)=dsqrt(2.d0)
  B(3,1)=-1.d0
  B(3,2)=-ione
  B=B/dsqrt(2.d0)
 
  Binv=(0.d0,0.d0)
  Binv(1,1)=1.d0
  Binv(1,3)=-1.d0
  Binv(2,1)=ione
  Binv(2,3)=ione
  Binv(3,2)=dsqrt(2.d0)

  Binv=Binv/dsqrt(2.d0)
 
!------------------------------------------------------

!
!      compose the 2nd rank tensor transformation matrix
!
  coef = zero
  qtran = zero

! Y1=xy
  coef(1,1,2) = rt2i
  coef(1,2,1) = rt2i

! Y2=yz
  coef(2,2,3) = rt2i
  coef(2,3,2) = rt2i

! Y3=zx
  coef(3,1,3) = rt2i
  coef(3,3,1) = rt2i

! Y4=x^2-y^2

  coef(4,1,1) = rt2i
  coef(4,2,2) = -rt2i

! Y5=3z^2-1=2z^2-x^2-y^2

  coef(5,1,1) = -rt6i
  coef(5,2,2) = -rt6i
  coef(5,3,3) = 2.d0*rt6i

  do i=1,5
   do j=1,5
    do k=1,ntrans
     do m=1,3
      do n=1,3
           cjmn = vtran(1,m,k)*(coef(j,1,1)*vtran(1,n,k) &
                             + coef(j,1,2)*vtran(2,n,k) &
                             + coef(j,1,3)*vtran(3,n,k)) &
               + vtran(2,m,k)*(coef(j,2,1)*vtran(1,n,k) &
                             + coef(j,2,2)*vtran(2,n,k) &
                             + coef(j,2,3)*vtran(3,n,k)) &
               + vtran(3,m,k)*(coef(j,3,1)*vtran(1,n,k) &
                             + coef(j,3,2)*vtran(2,n,k) &
                             + coef(j,3,3)*vtran(3,n,k)) 
           qtran(i,j,k) = qtran(i,j,k) + coef(i,m,n)*cjmn
      end do
     end do    
    end do
   end do
  end do



!  do k=1,ntrans
!     vtran(:,:,k) = abs(vtran(:,:,k))
!     qtran(:,:,k) = abs(qtran(:,:,k))
!  enddo

  do k=1,ntrans
     qtranc(:,:,k)=MATMUL(A,MATMUL(qtran(:,:,k),Ainv))
     vtranc(:,:,k)=MATMUL(B,MATMUL(vtran(:,:,k),Binv))
  end do

  deallocate(qtran)


  return

end subroutine findtransc


      double precision function a_coef( l, k, m_0, m_1, m_2, m_3 )
      
      IMPLICIT none
      INTEGER NDIM,IER
      PARAMETER( NDIM = 40)
      DOUBLE PRECISION THRCOF(NDIM)

      INTEGER l, k, m_0, m_1, m_2, m_3
      !-- The same numbers in double precision format (for function input)
      DOUBLE PRECISION dp_l, dp_k, dp_m_0, dp_m_1, dp_m_2, dp_m_3

!     m    = m_0
!     m'   = m_1
!     m''  = m_2
!     m''' = m_3

      double precision term    !-- Each of the terms in the sum
      DOUBLE PRECISION first_3j, second_3j, third_3j 

      INTEGER q,i
      DOUBLE PRECISION dp_q,m2_max,m2_min,m2
      DOUBLE PRECISION zero

      first_3j=0.d0
      second_3j=0.d0
      third_3j=0.d0

      if((m_0+m_2).ne.(m_1+m_3)) then
         a_coef=0.d0
      else

!      
      zero = 0.d0

      dp_l   = 1.d0 * l
      dp_k   = 1.d0 * k
      dp_m_0 = 1.d0 * m_0
      dp_m_1 = 1.d0 * m_1
      dp_m_2 = 1.d0 * m_2
      dp_m_3 = 1.d0 * m_3

      a_coef = 0.d0 

         q=m_0-m_1
         dp_q = 1.d0 * q 

         THRCOF=0.d0
         CALL DRC3JM( dp_l, dp_k, dp_l, zero, m2_min,m2_max, THRCOF,NDIM,IER)
         m2=m2_min-1.d0
         do i=1,NDIM
            m2=m2+1.d0
            if(abs(m2).lt.1.d-10) then
             first_3j=THRCOF(i)
             goto 111
            end if
         end do

 111     continue

         THRCOF=0.d0
         CALL DRC3JM( dp_l, dp_k, dp_l, -dp_m_0, m2_min,m2_max, THRCOF,NDIM,IER)

         m2=m2_min-1.d0
         do i=1,NDIM
            m2=m2+1.d0
            if(abs(m2-dp_q).lt.1.d-10) then
             second_3j=THRCOF(i)
             goto 112
            end if
         end do

 112     continue

         THRCOF=0.d0
         CALL DRC3JM( dp_l, dp_k, dp_l, -dp_m_2, m2_min,m2_max,  THRCOF,NDIM,IER)

         m2=m2_min-1.d0
         do i=1,NDIM
            m2=m2+1.d0
            if(abs(m2+dp_q).lt.1.d-10) then
             third_3j=THRCOF(i)
             goto 113
            end if
         end do

 113     continue

         a_coef = (2*l+1)*(2*l+1) *( (-1)**(m_0+q+m_2))
         a_coef = a_coef * first_3j * first_3j * second_3j * third_3j
      end if

      RETURN
      END



!     +--------------------------------------------------+
!     | http://www.netlib.org/slatec/src/drc3jm.f        |
!     +--------------------------------------------------+

!*DECK DRC3JM
     SUBROUTINE DRC3JM (L1, L2, L3, M1, M2MIN, M2MAX, THRCOF, NDIM, IER)

!***BEGIN PROLOGUE  DRC3JM
!***PURPOSE  Evaluate the 3j symbol g(M2) = (L1 L2   L3  )
!                                           (M1 M2 -M1-M2)
!            for all allowed values of M2, the other parameters
!            being held fixed.
!***LIBRARY   SLATEC
!***CATEGORY  C19
!***TYPE      DOUBLE PRECISION (RC3JM-S, DRC3JM-D)
!***KEYWORDS  3J COEFFICIENTS, 3J SYMBOLS, CLEBSCH-GORDAN COEFFICIENTS,
!             RACAH COEFFICIENTS, VECTOR ADDITION COEFFICIENTS,
!             WIGNER COEFFICIENTS
!***AUTHOR  Gordon, R. G., Harvard University
!           Schulten, K., Max Planck Institute
!***DESCRIPTION
!
! *Usage:
!
!        DOUBLE PRECISION L1, L2, L3, M1, M2MIN, M2MAX, THRCOF(NDIM)
!        INTEGER NDIM, IER
!
!        CALL DRC3JM (L1, L2, L3, M1, M2MIN, M2MAX, THRCOF, NDIM, IER)
!
! *Arguments:
!
!     L1 :IN      Parameter in 3j symbol.
!
!     L2 :IN      Parameter in 3j symbol.
!
!     L3 :IN      Parameter in 3j symbol.
!
!     M1 :IN      Parameter in 3j symbol.
!
!     M2MIN :OUT  Smallest allowable M2 in 3j symbol.
!
!     M2MAX :OUT  Largest allowable M2 in 3j symbol.
!
!     THRCOF :OUT Set of 3j coefficients generated by evaluating the
!                 3j symbol for all allowed values of M2.  THRCOF(I)
!                 will contain g(M2MIN+I-1), I=1,2,...,M2MAX-M2MIN+1.
!
!     NDIM :IN    Declared length of THRCOF in calling program.
!
!     IER :OUT    Error flag.
!                 IER=0 No errors.
!                 IER=1 Either L1.LT.ABS(M1) or L1+ABS(M1) non-integer.
!                 IER=2 ABS(L1-L2).LE.L3.LE.L1+L2 not satisfied.
!                 IER=3 L1+L2+L3 not an integer.
!                 IER=4 M2MAX-M2MIN not an integer.
!                 IER=5 M2MAX less than M2MIN.
!                 IER=6 NDIM less than M2MAX-M2MIN+1.
!
! *Description:
!
!     Although conventionally the parameters of the vector addition
!  coefficients satisfy certain restrictions, such as being integers
!  or integers plus 1/2, the restrictions imposed on input to this
!  subroutine are somewhat weaker. See, for example, Section 27.9 of
!  Abramowitz and Stegun or Appendix C of Volume II of A. Messiah.
!  The restrictions imposed by this subroutine are
!       1. L1.GE.ABS(M1) and L1+ABS(M1) must be an integer;
!       2. ABS(L1-L2).LE.L3.LE.L1+L2;
!       3. L1+L2+L3 must be an integer;
!       4. M2MAX-M2MIN must be an integer, where
!          M2MAX=MIN(L2,L3-M1) and M2MIN=MAX(-L2,-L3-M1).
!  If the conventional restrictions are satisfied, then these
!  restrictions are met.
!
!     The user should be cautious in using input parameters that do
!  not satisfy the conventional restrictions. For example, the
!  the subroutine produces values of
!       g(M2) = (0.75 1.50   1.75  )
!               (0.25  M2  -0.25-M2)
!  for M2=-1.5,-0.5,0.5,1.5 but none of the symmetry properties of the
!  3j symbol, set forth on page 1056 of Messiah, is satisfied.
!
!     The subroutine generates g(M2MIN), g(M2MIN+1), ..., g(M2MAX)
!  where M2MIN and M2MAX are defined above. The sequence g(M2) is
!  generated by a three-term recurrence algorithm with scaling to
!  control overflow. Both backward and forward recurrence are used to
!  maintain numerical stability. The two recurrence sequences are
!  matched at an interior point and are normalized from the unitary
!  property of 3j coefficients and Wigner's phase convention.
!
!    The algorithm is suited to applications in which large quantum
!  numbers arise, such as in molecular dynamics.
!
!***REFERENCES  1. Abramowitz, M., and Stegun, I. A., Eds., Handbook
!                  of Mathematical Functions with Formulas, Graphs
!                  and Mathematical Tables, NBS Applied Mathematics
!                  Series 55, June 1964 and subsequent printings.
!               2. Messiah, Albert., Quantum Mechanics, Volume II,
!                  North-Holland Publishing Company, 1963.
!               3. Schulten, Klaus and Gordon, Roy G., Exact recursive
!                  evaluation of 3j and 6j coefficients for quantum-
!                  mechanical coupling of angular momenta, J Math
!                  Phys, v 16, no. 10, October 1975, pp. 1961-1970.
!               4. Schulten, Klaus and Gordon, Roy G., Semiclassical
!                  approximations to 3j and 6j coefficients for
!                  quantum-mechanical coupling of angular momenta,
!                  J Math Phys, v 16, no. 10, October 1975,
!                  pp. 1971-1988.
!               5. Schulten, Klaus and Gordon, Roy G., Recursive
!                  evaluation of 3j and 6j coefficients, Computer
!                  Phys Comm, v 11, 1976, pp. 269-278.
!***ROUTINES CALLED  D1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   750101  DATE WRITTEN
!   880515  SLATEC prologue added by G. C. Nielson, NBS; parameters
!           HUGE and TINY revised to depend on D1MACH.
!   891229  Prologue description rewritten; other prologue sections
!           revised; MMATCH (location of match point for recurrences)
!           removed from argument list; argument IER changed to serve
!           only as an error flag (previously, in cases without error,
!           it returned the number of scalings); number of error codes
!           increased to provide more precise error information;
!           program comments revised; SLATEC error handler calls
!           introduced to enable printing of error messages to meet
!           SLATEC standards. These changes were done by D. W. Lozier,
!           M. A. McClain and J. M. Smith of the National Institute
!           of Standards and Technology, formerly NBS.
!   910415  Mixed type expressions eliminated; variable C1 initialized;
!           description of THRCOF expanded. These changes were done by
!           D. W. Lozier.
!***END PROLOGUE  DRC3JM
!
      INTEGER NDIM, IER
      DOUBLE PRECISION L1, L2, L3, M1, M2MIN, M2MAX, THRCOF(NDIM)
!
      INTEGER I, INDEX, LSTEP, N, NFIN, NFINP1, NFINP2, NFINP3, NLIM, NSTEP2
      DOUBLE PRECISION A1, A1S, C1, C1OLD, C2, CNORM, D1MACH, DV, EPS, &
                       HUGE, M2, M3, NEWFAC, OLDFAC, ONE, RATIO, SIGN1, &
                       SIGN2, SRHUGE, SRTINY, SUM1, SUM2, SUMBAC, &
                       SUMFOR, SUMUNI, THRESH, TINY, TWO, X, X1, X2, X3, &
                       Y, Y1, Y2, Y3, ZERO
 
      DATA  ZERO,EPS,ONE,TWO /0.0D0,0.01D0,1.0D0,2.0D0/
!
!***FIRST EXECUTABLE STATEMENT  DRC3JM
      IER=0
!  HUGE is the square root of one twentieth of the largest floating
!  point number, approximately.

!------ Emmanouil Kioupakis: -don't call d1mach, instead use a very large number
!------                      -don't print error messages 


!      HUGE = SQRT(D1MACH(2)/20.0D0)
      HUGE = SQRT(1.0D20/20.0D0)
      SRHUGE = SQRT(HUGE)
      TINY = 1.0D0/HUGE
      SRTINY = 1.0D0/SRHUGE
!
!     MMATCH = ZERO
!
!
!  Check error conditions 1, 2, and 3.
      IF((L1-ABS(M1)+EPS.LT.ZERO).OR.(MOD(L1+ABS(M1)+EPS,ONE).GE.EPS+EPS))THEN
         IER=1
!         CALL XERMSG('SLATEC','DRC3JM','L1-ABS(M1) less than zero or '//
!     +      'L1+ABS(M1) not integer.',IER,1)
         RETURN
      ELSEIF((L1+L2-L3.LT.-EPS).OR.(L1-L2+L3.LT.-EPS).OR.(-L1+L2+L3.LT.-EPS))THEN
         IER=2
!         CALL XERMSG('SLATEC','DRC3JM','L1, L2, L3 do not satisfy '//
!     +      'triangular condition.',IER,1)
         RETURN
      ELSEIF(MOD(L1+L2+L3+EPS,ONE).GE.EPS+EPS)THEN
         IER=3
!         CALL XERMSG('SLATEC','DRC3JM','L1+L2+L3 not integer.',IER,1)
         RETURN
      ENDIF
!
!
!  Limits for M2
      M2MIN = MAX(-L2,-L3-M1)
      M2MAX = MIN(L2,L3-M1)
!
!  Check error condition 4.
      IF(MOD(M2MAX-M2MIN+EPS,ONE).GE.EPS+EPS)THEN
         IER=4
!         CALL XERMSG('SLATEC','DRC3JM','M2MAX-M2MIN not integer.',IER,1)
         RETURN
      ENDIF
      IF(M2MIN.LT.M2MAX-EPS)   GO TO 20
      IF(M2MIN.LT.M2MAX+EPS)   GO TO 10
!
!  Check error condition 5.
      IER=5
!      CALL XERMSG('SLATEC','DRC3JM','M2MIN greater than M2MAX.',IER,1)
      RETURN
!
!
!  This is reached in case that M2 and M3 can take only one value.
   10 CONTINUE
!     MSCALE = 0
      THRCOF(1) = (-ONE) ** INT(ABS(L2-L3-M1)+EPS) /SQRT(L1+L2+L3+ONE)
      RETURN
!
!  This is reached in case that M1 and M2 take more than one value.
   20 CONTINUE
!     MSCALE = 0
      NFIN = INT(M2MAX-M2MIN+ONE+EPS)
      IF(NDIM-NFIN)   21, 23, 23
!
!  Check error condition 6.
   21 IER = 6
!      CALL XERMSG('SLATEC','DRC3JM','Dimension of result array for '// &
!                  '3j coefficients too small.',IER,1)
      RETURN
!
!
!
!  Start of forward recursion from M2 = M2MIN
!
   23 M2 = M2MIN
      THRCOF(1) = SRTINY
      NEWFAC = 0.0D0
      C1 = 0.0D0
      SUM1 = TINY
!
!
      LSTEP = 1
   30 LSTEP = LSTEP + 1
      M2 = M2 + ONE
      M3 = - M1 - M2
!
!
      OLDFAC = NEWFAC
      A1 = (L2-M2+ONE) * (L2+M2) * (L3+M3+ONE) * (L3-M3)
      NEWFAC = SQRT(A1)
!
!
      DV = (L1+L2+L3+ONE)*(L2+L3-L1) - (L2-M2+ONE)*(L3+M3+ONE) &
                                     - (L2+M2-ONE)*(L3-M3-ONE)
!
      IF(LSTEP-2)  32, 32, 31
!
   31 C1OLD = ABS(C1)
   32 C1 = - DV / NEWFAC
!
      IF(LSTEP.GT.2)   GO TO 60
!
!
!  If M2 = M2MIN + 1, the third term in the recursion equation vanishes,
!  hence
!
      X = SRTINY * C1
      THRCOF(2) = X
      SUM1 = SUM1 + TINY * C1*C1
      IF(LSTEP.EQ.NFIN)   GO TO 220
      GO TO 30
!
!
   60 C2 = - OLDFAC / NEWFAC
!
!  Recursion to the next 3j coefficient
      X = C1 * THRCOF(LSTEP-1) + C2 * THRCOF(LSTEP-2)
      THRCOF(LSTEP) = X
      SUMFOR = SUM1
      SUM1 = SUM1 + X*X
      IF(LSTEP.EQ.NFIN)   GO TO 100
!
!  See if last unnormalized 3j coefficient exceeds SRHUGE
!
      IF(ABS(X).LT.SRHUGE)   GO TO 80
!
!  This is reached if last 3j coefficient larger than SRHUGE,
!  so that the recursion series THRCOF(1), ... , THRCOF(LSTEP)
!  has to be rescaled to prevent overflow
!
!     MSCALE = MSCALE + 1
      DO 70 I=1,LSTEP
      IF(ABS(THRCOF(I)).LT.SRTINY)   THRCOF(I) = ZERO
   70 THRCOF(I) = THRCOF(I) / SRHUGE
      SUM1 = SUM1 / HUGE
      SUMFOR = SUMFOR / HUGE
      X = X / SRHUGE
!
!
!  As long as ABS(C1) is decreasing, the recursion proceeds towards
!  increasing 3j values and, hence, is numerically stable.  Once
!  an increase of ABS(C1) is detected, the recursion direction is
!  reversed.
!
   80 IF(C1OLD-ABS(C1))   100, 100, 30
!
!
!  Keep three 3j coefficients around MMATCH for comparison later
!  with backward recursion values.
!
  100 CONTINUE
!     MMATCH = M2 - 1
      NSTEP2 = NFIN - LSTEP + 3
      X1 = X
      X2 = THRCOF(LSTEP-1)
      X3 = THRCOF(LSTEP-2)
!
!  Starting backward recursion from M2MAX taking NSTEP2 steps, so
!  that forwards and backwards recursion overlap at the three points
!  M2 = MMATCH+1, MMATCH, MMATCH-1.
!
      NFINP1 = NFIN + 1
      NFINP2 = NFIN + 2
      NFINP3 = NFIN + 3
      THRCOF(NFIN) = SRTINY
      SUM2 = TINY
!
!
!
      M2 = M2MAX + TWO
      LSTEP = 1
  110 LSTEP = LSTEP + 1
      M2 = M2 - ONE
      M3 = - M1 - M2
      OLDFAC = NEWFAC
      A1S = (L2-M2+TWO) * (L2+M2-ONE) * (L3+M3+TWO) * (L3-M3-ONE)
      NEWFAC = SQRT(A1S)
      DV = (L1+L2+L3+ONE)*(L2+L3-L1) - (L2-M2+ONE)*(L3+M3+ONE) &
                                     - (L2+M2-ONE)*(L3-M3-ONE)
      C1 = - DV / NEWFAC
      IF(LSTEP.GT.2)   GO TO 120
!
!  If M2 = M2MAX + 1 the third term in the recursion equation vanishes
!
      Y = SRTINY * C1
      THRCOF(NFIN-1) = Y
      IF(LSTEP.EQ.NSTEP2)   GO TO 200
      SUMBAC = SUM2
      SUM2 = SUM2 + Y*Y
      GO TO 110
!
  120 C2 = - OLDFAC / NEWFAC
!
!  Recursion to the next 3j coefficient
!
      Y = C1 * THRCOF(NFINP2-LSTEP) + C2 * THRCOF(NFINP3-LSTEP)
!
      IF(LSTEP.EQ.NSTEP2)   GO TO 200
!
      THRCOF(NFINP1-LSTEP) = Y
      SUMBAC = SUM2
      SUM2 = SUM2 + Y*Y
!
!
!  See if last 3j coefficient exceeds SRHUGE
!
      IF(ABS(Y).LT.SRHUGE)   GO TO 110
!
!  This is reached if last 3j coefficient larger than SRHUGE,
!  so that the recursion series THRCOF(NFIN), ... , THRCOF(NFIN-LSTEP+1)
!  has to be rescaled to prevent overflow.
!
!     MSCALE = MSCALE + 1
      DO 111 I=1,LSTEP
      INDEX = NFIN - I + 1
      IF(ABS(THRCOF(INDEX)).LT.SRTINY) THRCOF(INDEX) = ZERO
  111 THRCOF(INDEX) = THRCOF(INDEX) / SRHUGE
      SUM2 = SUM2 / HUGE
      SUMBAC = SUMBAC / HUGE
!
      GO TO 110
!
!
!
!  The forward recursion 3j coefficients X1, X2, X3 are to be matched
!  with the corresponding backward recursion values Y1, Y2, Y3.
!
  200 Y3 = Y
      Y2 = THRCOF(NFINP2-LSTEP)
      Y1 = THRCOF(NFINP3-LSTEP)
!
!
!  Determine now RATIO such that YI = RATIO * XI  (I=1,2,3) holds
!  with minimal error.
!
      RATIO = ( X1*Y1 + X2*Y2 + X3*Y3 ) / ( X1*X1 + X2*X2 + X3*X3 )
      NLIM = NFIN - NSTEP2 + 1
!
      IF(ABS(RATIO).LT.ONE)   GO TO 211
!
      DO 210 N=1,NLIM
  210 THRCOF(N) = RATIO * THRCOF(N)
      SUMUNI = RATIO * RATIO * SUMFOR + SUMBAC
      GO TO 230
!
  211 NLIM = NLIM + 1
      RATIO = ONE / RATIO
      DO 212 N=NLIM,NFIN
  212 THRCOF(N) = RATIO * THRCOF(N)
      SUMUNI = SUMFOR + RATIO*RATIO*SUMBAC
      GO TO 230
!
  220 SUMUNI = SUM1
!
!
!  Normalize 3j coefficients
!
  230 CNORM = ONE / SQRT((L1+L1+ONE) * SUMUNI)
!
!  Sign convention for last 3j coefficient determines overall phase
!
      SIGN1 = SIGN(ONE,THRCOF(NFIN))
      SIGN2 = (-ONE) ** INT(ABS(L2-L3-M1)+EPS)
      IF(SIGN1*SIGN2)  235,235,236
  235 CNORM = - CNORM
!
  236 IF(ABS(CNORM).LT.ONE)   GO TO 250
!
      DO 240 N=1,NFIN
  240 THRCOF(N) = CNORM * THRCOF(N)
      RETURN
!
  250 THRESH = TINY / ABS(CNORM)
      DO 251 N=1,NFIN
      IF(ABS(THRCOF(N)).LT.THRESH)   THRCOF(N) = ZERO
  251 THRCOF(N) = CNORM * THRCOF(N)
!
!
!
      RETURN
      END

        
!-------------------------------------------------------------------------

  subroutine construct_Hm2(cnk_LSDA,dmat,Vmat,Hmat,Etot_corr,ELSDA_corr,nrk,ndim,nband0, &
                          HU,SJ,idc,alpha_all,lambda,Nsite,lproj,nlmax,bands,myproc,iperturb)

  include 'use.h'  
  implicit none                       ! implicit? No!
  include 'interface.h'  

  type(band), intent(in) :: bands  
  integer::nrk,nband0,ndim,Nsite,myproc,idc,lproj(Nsite),nlmax,iperturb

  double precision::HU(Nsite),SJ(Nsite),Etot_corr,ELSDA_corr,alpha_all    ! correction to total energy
  double precision::lambda(nlmax,2,Nsite)
  double precision::Vmat(nlmax,nlmax,nlmax,nlmax,Nsite)
  double complex::dmat(nlmax,nlmax,2,Nsite), cnk_LSDA(nlmax,ndim,nrk,2,Nsite)
  double complex::Hmat(ndim,ndim,nrk,2)

  double complex,allocatable::corr(:,:,:,:,:)
  double precision,allocatable::alpha(:)
  double complex,allocatable::Vnm(:,:,:,:)

  double precision,allocatable::avg_occup(:,:),ecorr(:)

  integer::im,jm,is,irk,ib,ia,jb,m1,m2,m3,m4,nl

  allocate(corr(ndim,ndim,nrk,2,Nsite))
  allocate(Vnm(nlmax,nlmax,2,Nsite))
  allocate(avg_occup(2,Nsite))
  allocate(ecorr(Nsite))

!--------------------------------------------------------------------------

    avg_occup=0.d0

    do ia=1,Nsite     ! loop over atomic sites
       nl=2*lproj(ia)+1
       do is=1,bands%nspin
          do im=1,nl
             avg_occup(is,ia)=avg_occup(is,ia)+REAL(dmat(im,im,is,ia))
          end do
       end do   ! spin
    end do   !over atomic sites

!----------------------------------------------------------------------
! potential matrix

  Vnm=(0.d0,0.d0)
   

  if (bands%nspin.eq.1) then

  do ia=1,Nsite 
     nl=2*lproj(ia)+1

  do is = 1, bands%nspin
     do m1=1,nl
     do m2=1,nl
        do m3=1,nl
        do m4=1,nl

! same spin: Coulomb and exchange
           Vnm(m1,m2,is,ia)=Vnm(m1,m2,is,ia)+(Vmat(m1,m2,m3,m4,ia)-Vmat(m1,m4,m3,m2,ia))*dmat(m3,m4,is,ia)
          

! difference spin: Coulomb only
           Vnm(m1,m2,is,ia)=Vnm(m1,m2,is,ia)+Vmat(m1,m2,m3,m4,ia)*dmat(m3,m4,is,ia)

        end do
        end do
     end do
     end do

  end do  ! spin
  end do  ! atomic site

  else


  do ia=1,Nsite 
     nl=2*lproj(ia)+1

  do is = 1, bands%nspin
     do m1=1,nl
     do m2=1,nl
        do m3=1,nl
        do m4=1,nl

! same spin: Coulomb and exchange
           Vnm(m1,m2,is,ia)=Vnm(m1,m2,is,ia)+(Vmat(m1,m2,m3,m4,ia)-Vmat(m1,m4,m3,m2,ia))*dmat(m3,m4,is,ia)
          

! difference spin: Coulomb only
           if(is.eq.1) then
           Vnm(m1,m2,is,ia)=Vnm(m1,m2,is,ia)+Vmat(m1,m2,m3,m4,ia)*dmat(m3,m4,2,ia)
           else
           Vnm(m1,m2,is,ia)=Vnm(m1,m2,is,ia)+Vmat(m1,m2,m3,m4,ia)*dmat(m3,m4,1,ia)
           end if

        end do
        end do
     end do
     end do

  end do  ! spin
  end do  ! atomic site
  end if



!-----------------------------------------------------------
! The total energy correction due to eigenvalues changes
!

! first term

  if(bands%nspin.eq.1) avg_occup(2,:)=avg_occup(1,:)

  ecorr=0.d0
  do ia=1,Nsite
     nl=2*lproj(ia)+1
  do is=1,bands%nspin
     do im=1,nl
     do jm=1,nl
        ecorr(ia)=ecorr(ia)-0.5d0*REAL(dmat(im,jm,is,ia)*Vnm(im,jm,is,ia))*2.d0/bands%nspin
     end do
     end do
  end do
  end do

! second term

  do ia=1,Nsite
     ecorr(ia)=ecorr(ia)+0.5d0*HU(ia)*(avg_occup(1,ia)+avg_occup(2,ia))* &
                                     (avg_occup(1,ia)+avg_occup(2,ia)-1.d0)
     ecorr(ia)=ecorr(ia)-0.5d0*SJ(ia)*(avg_occup(1,ia)*(avg_occup(1,ia)-1.d0)+ &
                                 avg_occup(2,ia)*(avg_occup(2,ia)-1.d0))
  end do

! sum of the first and second terms gives the energy difference between 
! the LSDA and LSDA+U energy functionals

  ELSDA_corr=0.d0
  do ia=1,Nsite
     ELSDA_corr=ELSDA_corr+ecorr(ia)
  end do

! The third term

  do ia=1,Nsite
     ecorr(ia)=ecorr(ia)+0.5d0*(HU(ia)-SJ(ia))*(avg_occup(1,ia)+avg_occup(2,ia))
  end do

!
! Finally, this is the total energy double counting in the band structure energy
! that need to be taken out 

  Etot_corr=0.d0
  do ia=1,Nsite
     Etot_corr=Etot_corr+ecorr(ia)
  end do
  

!--------------------------------------------------------------------------------------
  do ia=1,Nsite 
  nl=2*lproj(ia)+1
  do is = 1, bands%nspin  
!     write(9,'a10,2i5,2f10.5') "ia,is",ia,is,avg_occup(1,ia),avg_occup(2,ia)
     do m1=1,nl
        Vnm(m1,m1,is,ia)=Vnm(m1,m1,is,ia)-HU(ia)*(avg_occup(1,ia)+avg_occup(2,ia)-0.5d0)+&
                         SJ(ia)*(avg_occup(is,ia)-0.5d0)
!        write(9,*) "m1",m1,Vnm(m1,m1,is,ia)
     end do
  end do
  end do
 
  corr=(0.d0,0.d0)
  do ia=1,Nsite 
     nl=2*lproj(ia)+1
  do is = 1, bands%nspin  
     do irk = 1, nrk  
        do ib = 1,ndim
        do jb = 1,ndim

           do im=1,nl
           do jm=1,nl
              corr(ib,jb,irk,is,ia)=corr(ib,jb,irk,is,ia)+Vnm(im,jm,is,ia)* &
                   DCONJG(cnk_LSDA(im,ib,irk,is,ia))*cnk_LSDA(jm,jb,irk,is,ia)
           end do
           end do

        end do
        end do
     end do
  end do  ! spin
  end do  ! atomic site

  Hmat(:,:,:,:)=(0.d0,0.d0)
  do ia=1,Nsite
     Hmat(:,:,:,:)=Hmat(:,:,:,:)+corr(:,:,:,:,ia)
  end do

 120 format(40f9.5)

  do is = 1, bands%nspin
     do irk = 1, nrk
       do ib = 1,ndim
           Hmat(ib,ib,irk,is)=Hmat(ib,ib,irk,is)+bands%energy(ib+nband0, irk, is)
        end do
     end do
  end do
   
 deallocate(corr)
 deallocate(Vnm)
 deallocate(avg_occup)
 deallocate(ecorr)


 200  format(i5,10f10.5)
 210  format(a5,9a10)

  return
end subroutine construct_Hm2

  subroutine uforce(HU,SJ,force_U,drdmat,dmat,Vmat,Nsite,lproj,nlmax,nspin)

  include 'use.h'  
  implicit none                       ! implicit? No!
  include 'interface.h'  

  integer::Nsite,lproj(Nsite),nlmax,nspin

  double precision::HU(Nsite),SJ(Nsite)
  double precision::force_U(Nsite,3)
  double precision::Vmat(nlmax,nlmax,nlmax,nlmax,Nsite)
  double complex::drdmat(nlmax,nlmax,2,Nsite,3)
  double complex::dmat(nlmax,nlmax,2,Nsite)


  integer::is,ia,m1,m2,m3,m4,nl,ix

  double precision::occup_sum
  double complex::drdmat_sum
  
  double precision,allocatable::occup(:,:)
  double complex,allocatable::Vnm(:,:,:,:)


  allocate(occup(2,Nsite))
  allocate(Vnm(nlmax,nlmax,2,Nsite))

  force_U=0.d0
  Vnm=(0.d0,0.d0)
   
  occup=0.d0

  do ia=1,Nsite     ! loop over atomic sites
     nl=2*lproj(ia)+1
     do is=1,nspin
        do m1=1,nl
           occup(is,ia)=occup(is,ia)+REAL(dmat(m1,m1,is,ia))
        end do
     end do   ! spin
  end do   !over atomic sites


  do ia=1,Nsite 
     nl=2*lproj(ia)+1

  do is = 1, nspin
     do m1=1,nl
     do m2=1,nl
        do m3=1,nl
        do m4=1,nl
  
! same spin: Coulomb and exchange
           Vnm(m1,m2,is,ia)=Vnm(m1,m2,is,ia)+(Vmat(m1,m2,m3,m4,ia)-Vmat(m1,m4,m3,m2,ia))*dmat(m3,m4,is,ia)

! difference spin: Coulomb only
           if(nspin.eq.1) then
             Vnm(m1,m2,is,ia)=Vnm(m1,m2,is,ia)+Vmat(m1,m2,m3,m4,ia)*dmat(m3,m4,is,ia)
           else
             if(is.eq.1) then
             Vnm(m1,m2,is,ia)=Vnm(m1,m2,is,ia)+Vmat(m1,m2,m3,m4,ia)*dmat(m3,m4,2,ia)
             else
             Vnm(m1,m2,is,ia)=Vnm(m1,m2,is,ia)+Vmat(m1,m2,m3,m4,ia)*dmat(m3,m4,1,ia)
             end if
           end if

        end do
        end do
     end do
     end do

  end do  ! spin
  end do  ! atomic site


! force due to E^u

  do ix=1,3
  do ia=1,Nsite 
     nl=2*lproj(ia)+1

  do is = 1, nspin
     do m1=1,nl
     do m2=1,nl

       force_U(ia,ix)=force_U(ia,ix)-REAL(Vnm(m1,m2,is,ia)*drdmat(m1,m2,is,ia,ix))

     end do
     end do

  end do  ! spin
  end do  ! atomic site
  end do

  if(nspin.eq.1) force_U=force_U*2.d0
 

! double counting term  E^DC

! U term

  do ix=1,3
  do ia=1,Nsite

  nl=2*lproj(ia)+1
  if(nspin.eq.2) then
  occup_sum=occup(1,ia)+occup(2,ia)
  else
  occup_sum=occup(1,ia)*2.d0
  end if

  drdmat_sum=0.d0

  do is = 1, nspin
     do m1=1,nl
     drdmat_sum=drdmat_sum+drdmat(m1,m1,is,ia,ix)
     end do
  end do
  if(nspin.eq.1) drdmat_sum=drdmat_sum*2.d0

  force_U(ia,ix)=force_U(ia,ix)+ REAL(HU(ia)*(occup_sum-0.5d0)*drdmat_sum)

  end do
  end do

! J term

  do ix=1,3
  do ia=1,Nsite

  nl=2*lproj(ia)+1

  if(nspin.eq.2) then
  do is = 1, nspin
     drdmat_sum=0.d0
     do m1=1,nl
     drdmat_sum=drdmat_sum+drdmat(m1,m1,is,ia,ix)
     end do
     force_U(ia,ix)=force_U(ia,ix)- REAL(SJ(ia)*(occup(is,ia)-0.5d0)*drdmat_sum)
  end do
  else
     drdmat_sum=0.d0
     do m1=1,nl
     drdmat_sum=drdmat_sum+drdmat(m1,m1,1,ia,ix)
     end do
     force_U(ia,ix)=force_U(ia,ix)- REAL(SJ(ia)*(occup(1,ia)-0.5d0)*drdmat_sum)*2.d0
  end if
  

  end do
  end do

! still need debugging

!  force_U=0.d0

! it seems that a factor of 2 is missing somewhere, check it out please!

!  force_U=force_U/2.d0

  write(9,*) "---------------------------------"
  write(9,*) "Hubbard U contributions to forces"
!
  do ia=1,Nsite 
     write(9,'(6f10.5)')(force_U(ia,ix),ix=1,3) 
  end do
  write(9,*) "---------------------------------"

  deallocate(Vnm)
  deallocate(occup)

  return
end subroutine uforce


!---------------------------------------------------------------------------
! density matrix from occupied and unoccupied states
!

     subroutine construct_dm1(syms,cnk,cnk_f,dmat,Nsite,lproj,nlmax,ndim,nband0,kpoints,nspin)


     include 'use.h'
     implicit none                       ! implicit? No!
     include 'interface.h'

     type(kpoint), intent(in) :: kpoints
     type(symmetry), intent(in) :: syms

     integer,intent(IN)::Nsite,ndim,nband0,lproj(Nsite),nlmax,nspin
     double complex::cnk(nlmax,ndim,kpoints%nrk,2,Nsite)
     double complex,intent(in)::cnk_f(nlmax,ndim,kpoints%nrk,2,Nsite,syms%nfrac)
     double complex::dmat(nlmax,nlmax,2,Nsite)

! local variables
     integer::ii,jj,kk,is,ia,irk,ib,jb,nrk,nl,it,ntrans,imap
      
     double complex,allocatable::qtranc(:,:,:),vtranc(:,:,:)
     double complex,allocatable::dmat_tem(:,:,:,:),dmat_f(:,:,:,:,:)

!--------------------------------------------------------------------------

    nrk=kpoints%nrk

    dmat=(0.d0,0.d0)
    allocate(dmat_tem(nlmax,nlmax,2,Nsite))


    do ia=1,Nsite     ! loop over atomic sites
       nl=2*lproj(ia)+1
       do is=1,nspin
          do ii=1,nl
             do jj=1,nl

             do irk = 1, nrk
                do ib = 1,ndim
                   dmat(ii,jj,is,ia)=dmat(ii,jj,is,ia)+cnk(ii,ib,irk,is,ia)* &
                   DCONJG(cnk(jj,ib,irk,is,ia))*kpoints%w(irk)
                end do
             end do
          end do 
          end do  

       end do   ! spin
    end do      ! atomic sites


    if(kpoints%reduced) then

!    if(syms%nfrac.ne.0) then
!       allocate(dmat_f(nlmax,nlmax,2,Nsite,syms%nfrac))
!       dmat_f=(0.d0,0.d0)
!
!       do it=1,syms%nfrac ! loop over fractional translations
!
!       do ia=1,Nsite     ! loop over atomic sites
!       nl=2*lproj(ia)+1
!       do is=1,nspin
!          do ii=1,nl
!             do jj=1,nl
!
!             do irk = 1, nrk
!                do ib = 1,ndim
!                   dmat_f(ii,jj,is,ia,it)=dmat_f(ii,jj,is,ia,it)+cnk_f(ii,ib,irk,is,ia,it)* &
!                   DCONJG(cnk_f(jj,ib,irk,is,ia,it))*kpoints%w(irk)
!                end do
!             end do
!          end do
!          end do
!
!       end do   ! spin
!       end do      ! atomic sites
!       end do
!
!    end if


!------------------------------------------------------
    ntrans=syms%ntrans
!
!    do ii=1,syms%ntrans
!       if ((SUM(abs(syms%tnp(:,ii))).gt.1.d-6).and.(syms%nfrac.eq.0)) ntrans=ntrans-1
!    end do

!---------------------------------------------------------
! symmetrize density matrix

    allocate(qtranc(5,5,syms%ntrans))
    allocate(vtranc(3,3,syms%ntrans))

    call findtransc(syms%ntrans,syms%rsymmat,qtranc,vtranc)

    dmat_tem=(0.d0,0.d0)
    do ia=1,Nsite
       nl=2*lproj(ia)+1
       do is=1,nspin

          if(lproj(ia).eq.2) then
          it=0
     

          do ii=1,syms%ntrans
!             if ((SUM(abs(syms%tnp(:,ii))).gt.1.d-6).and.(syms%nfrac.gt.0)) then
!             it=it+1
!             dmat_tem(1:nl,1:nl,is,ia)=dmat_tem(1:nl,1:nl,is,ia)+ &
!             MATMUL(qtranc(:,:,ii), MATMUL(dmat_f(1:nl,1:nl,is,ia,it),DCONJG(TRANSPOSE(qtranc(:,:,ii)))))
!             end if

!             if (SUM(abs(syms%tnp(:,ii))).lt.1.d-6) then
             imap=syms%ind_rot(ii,ia)

!             dmat_tem(1:nl,1:nl,is,ia)=dmat_tem(1:nl,1:nl,is,ia)+ &
!             MATMUL(qtranc(:,:,ii), MATMUL(dmat(1:nl,1:nl,is,imap),DCONJG(TRANSPOSE(qtranc(:,:,ii)))))

             dmat_tem(1:nl,1:nl,is,imap)=dmat_tem(1:nl,1:nl,is,imap)+ &
             MATMUL(qtranc(:,:,ii), MATMUL(dmat(1:nl,1:nl,is,ia),DCONJG(TRANSPOSE(qtranc(:,:,ii)))))

!             end if

          end do
!          dmat(1:nl,1:nl,is,ia)=dmat_tem(1:nl,1:nl)
         end if

          if(lproj(ia).eq.1) then
          it=0

          do ii=1,syms%ntrans
!             if ((SUM(abs(syms%tnp(:,ii))).gt.1.d-6).and.(syms%nfrac.gt.0)) then
!             it=it+1
!             dmat_tem(1:nl,1:nl,is,ia)=dmat_tem(1:nl,1:nl,is,ia)+ &
!             MATMUL(vtranc(:,:,ii), MATMUL(dmat_f(1:nl,1:nl,is,ia,it),DCONJG(TRANSPOSE(vtranc(:,:,ii)))))
!             end if

!             if (SUM(abs(syms%tnp(:,ii))).lt.1.d-6) then
             imap=syms%ind_rot(ii,ia)

!             dmat_tem(1:nl,1:nl,is,ia)=dmat_tem(1:nl,1:nl,is,ia)+ &
!             MATMUL(vtranc(:,:,ii), MATMUL(dmat(1:nl,1:nl,is,imap),DCONJG(TRANSPOSE(vtranc(:,:,ii)))))

             dmat_tem(1:nl,1:nl,is,imap)=dmat_tem(1:nl,1:nl,is,imap)+ &
             MATMUL(vtranc(:,:,ii), MATMUL(dmat(1:nl,1:nl,is,ia),DCONJG(TRANSPOSE(vtranc(:,:,ii)))))

!             end if

          end do
!          dmat(1:nl,1:nl,is,ia)=dmat_tem(1:nl,1:nl)
          end if

          if(lproj(ia).eq.0) then
          dmat_tem(1:nl,1:nl,is,ia)=dmat(1:nl,1:nl,is,ia)*ntrans
          end if

       end do
    end do

    dmat=dmat_tem/ntrans

    deallocate(qtranc)
    deallocate(vtranc)
!    if(syms%nfrac.ne.0) deallocate(dmat_f)

    end if

    do ia=1,Nsite
       nl=2*lproj(ia)+1
       do is=1,nspin
          dmat_tem(:,:,is,ia)=dmat(:,:,is,ia)+DCONJG(TRANSPOSE(dmat(:,:,is,ia)))
       end do
    end do

    dmat=dmat_tem/2.d0

! make diagonal terms real
!    do ia=1,Nsite
!       nl=2*lproj(ia)+1
!       do is=1,nspin
!          do ii=1,nl
!             dmat(ii,ii,is,ia)=REAL(dmat(ii,ii,is,ia))
!          end do
!       end do
!    end do

    deallocate(dmat_tem)

    return
    end


    subroutine write_dmat(myproc,dmat,dmat_vec,nlmax,lproj,Nsite,idiag,nspin)


    implicit none
    integer::nlmax,Nsite,idiag,myproc,nspin
    integer::lproj(Nsite)
    double complex::dmat(nlmax,nlmax,2,Nsite),dmat_vec(nlmax,nlmax,2,Nsite)

    integer::ia,is,ii,nl
    double precision,allocatable::avg_occup(:,:)

    integer::LWORK,ier,jj
    double precision,allocatable::RWORK(:),EIG(:) 
    double complex,allocatable::WORK(:),dmat_tem(:,:)

    allocate(avg_occup(2,Nsite))

    avg_occup=0.d0


    do ia=1,Nsite     ! loop over atomic sites
       nl=2*lproj(ia)+1

    if(myproc.eq.0) then
       do is=1,nspin
          write(240,*) "atomic site, spin",ia,is
          do ii=1,nl
             write(240,*) ii,REAL(dmat(ii,ii,is,ia))
             avg_occup(is,ia)=avg_occup(is,ia)+REAL(dmat(ii,ii,is,ia))
          end do
       end do   ! spin

       avg_occup(:,ia)=avg_occup(:,ia)/nl
     
       write(240,*) "------------------------------------------------------"
       write(240,111) "Avg. occ/orb/spin           ",avg_occup(1,ia),avg_occup(2,ia)
       write(240,111) "Average occupation/spin     ",avg_occup(1,ia)*nl,avg_occup(2,ia)*nl
       write(240,111) "Total d(p) electrons        ",(avg_occup(1,ia)+avg_occup(2,ia))*nl
       write(240,111) "Magnetic moment on this site", abs(avg_occup(1,ia)-avg_occup(2,ia))*nl
       write(240,*) "------------------------------------------------------" 
       end if

 111   format(a28,1x,2f14.10)

    
       if(idiag.eq.1) then
          LWORK=5*nl

          allocate(dmat_tem(nl,nl))
          allocate(WORK(LWORK))
          allocate(RWORK(max(3*nl-2,1)))
          allocate(EIG(nl))

          if(myproc.eq.0) then
          write(240,*) 
          write(240,*) "Eigenvalues and eigenvectors of the density matrix"
          write(240,*) "--------------------------------------------------"
          end if

          do is=1,nspin
             dmat_tem=dmat(1:nl,1:nl,is,ia)

             call ZHEEV("V", "U", nl, dmat_tem, nl, EIG, WORK, LWORK, RWORK, ier )

             if(myproc.eq.0) then
             write(240,*) "Spin", is

             do ii=1,nl
                write(240,121) EIG(ii),(dmat_tem(jj,ii),jj=1,nl)
             end do
             end if

             dmat_vec(1:nl,1:nl,is,ia)=dmat_tem(:,:)
          end do

          deallocate(dmat_tem)
          deallocate(WORK)
          deallocate(RWORK)
          deallocate(EIG)
          if(myproc.eq.0) Write(240,*) "--------------------------------------------------"
         end if

    end do   !over atomic sites

    if(myproc.eq.0) call myflush(240)

    deallocate(avg_occup)

 121 format(f10.6,10f8.4)

    return
    end


!---------------------------------------------------------------------------


    subroutine lsdau_init(imyproc,inproc,mb,nrk,imode,crys,syms,kpoints)

! LSDA+U module by P. Zhang, 2003-2004 while at UC Berkeley

  use lsdau_shared

  include 'use.h'
  implicit none               ! never remove this line.

  integer::imyproc,mb,nrk,inproc
  type(crystal) :: crys 
  type(symmetry) :: syms      ! symmetry operations
  type(kpoint), intent(in) :: kpoints     ! BZ integration kpoints



! local variables

  character*1 symbol1,symbol2,symbol_tem
  character (LEN=2):: symbol
  integer:: is,ia,imode,nt,error
  integer:: i,ii,kk,im,jm,ier,ind,ngrid1max,ngrid2,j,k,isite,iatom
  integer::LEFT

  double precision:: delta,step2
  double precision:: r1,r2,r3,r4,r5,r6,r7,r8,norm

  double precision,allocatable:: rr(:),YP(:),YPP(:),YPPP(:),rad_int(:)

!--------------------------------------
  init_flag=999
  myproc=imyproc
  nproc=inproc


!----------------------------------------------
! there was a version in which I believe I fixed
! the problem ...
! for the moment, I will simply get rid of these
! trouble making symmetry operations.
! 
! the generate_kpoint.f90 is modified so that
! when it detects lsda+u calculation, no operations
! with fractional translation will be generated.
! 
! disable rotations with fractional translation
!
!----------------------------------------------

!  syms%nfrac=0
!  if(kpoints%reduced) then
!  do i=1,syms%ntrans
!     if (SUM(abs(syms%tnp(:,i))).gt.1.d-6) then
!        syms%nfrac=syms%nfrac+1
!        syms%indfrac(syms%nfrac)=i
!     end if
!  end do
!  end if


  if(syms%nfrac.ne.0) then
     write(9,*) "Number of fractional translations",syms%nfrac
     write(9,*) (syms%indfrac(ii),ii=1,syms%nfrac)
  end if

! to make program simple, all processors open
! the file, read in the data.
! The better was is only proc0 does the reading
! and broadcast the data. 

  open(89,file="lsdau.in")

  if (myproc .eq. 0) then

! this call causes some problem
!     call mvfile('occup.dat', '.old')
     open(unit=240,file='occup.dat')
!     open(unit=242,file='pdos_LSDA.dat')
!     open(unit=243,file='pdos_LSDAU.dat')  ! projection onto spherical harmonics
!     open(unit=244,file='pdos_LSDAUd.dat') ! projection onto eigenvectors of the density matrix
!     open(unit=245,file='pdos_LSDAd.dat') ! projection onto eigenvectors of the density matrix
  end if

  Nsite=0
  Nssite=0
  Npsite=0
  Ndsite=0
  Nfsite=0

  nt=crys%ntype

  allocate(HU_tem(nt))
  allocate(SJ_tem(nt))
  allocate(lproj_tem(nt))
  allocate(proj_flag(nt))
  allocate(inorm(nt))
  allocate(ngrid1(nt))
  allocate(rad_int(nt))

  lproj_tem=-1

! make it simple:

  inorm(:) = 0
  step = 0.015d0   ! 0.015au
  proj_flag(:)=1

  do ii=1,crys%ntype
     read(89,*) lproj_tem(ii),HU_tem(ii),SJ_tem(ii),rad_int(ii)
     if(lproj_tem(ii).lt.0) then
      ngrid1(ii)=0
     else
     ngrid1(ii)=int(rad_int(ii)/step)
     end if
  end do

  deallocate(rad_int)

!  read(89,*) step

  ngrid1max=MAXVAL(ngrid1)
  nlmax=MAXVAL(lproj_tem)*2+1

  allocate(Rnl_intp(ngrid1max,crys%ntype))

  do ii=1,crys%ntype
     if(lproj_tem(ii).eq.3) Nfsite=Nfsite+crys%natom(ii)
     if(lproj_tem(ii).eq.2) Ndsite=Ndsite+crys%natom(ii)
     if(lproj_tem(ii).eq.1) Npsite=Npsite+crys%natom(ii)
     if(lproj_tem(ii).eq.0) Nssite=Nssite+crys%natom(ii)
  end do
  Nsite=Nfsite+Ndsite+Npsite+Nssite

  allocate(lproj(Nsite))
  allocate(mproj(crys%ntype))
  allocate(nproj(1,crys%ntype))
  allocate(ngrid(Nsite))
  allocate(ntp(Nsite))
  allocate(lambda(nlmax,2,Nsite))


  nproj=1

  mproj=0

  allocate(HU(Nsite))
  allocate(SJ(Nsite))

  kk=1
  do ii=1,crys%ntype
     if(lproj_tem(ii).ge.0) then
       lproj(kk:kk+crys%natom(ii)-1)=lproj_tem(ii)
       mproj(ii)=1 
       HU(kk:kk+crys%natom(ii)-1)=HU_tem(ii)
       SJ(kk:kk+crys%natom(ii)-1)=SJ_tem(ii)
       ngrid(kk:kk+crys%natom(ii)-1)=ngrid1(ii)
       ntp(kk:kk+crys%natom(ii)-1)=ii
 
       kk=kk+crys%natom(ii)
     end if
  end do

  isite=0
  iatom=0
  do ii=1,crys%ntype
     if(mproj(ii).ge.1) then
     do k=1,mproj(ii)
     do j=1,crys%natom(ii)
        isite=isite+1
        ind_site_atom(isite)=iatom+j
        ind_isite_site(k,iatom+j)=isite
        ind_site_isite(isite)=k
     end do
     end do
     end if
     iatom=iatom+crys%natom(ii)
  end do

  interC=0

  miter=1

  read(89,*) nband0,nband1

!  read(89,*) miter
!  if((imode.eq.0).or.(imode.eq.2)) miter=1
!

!  read(89,*) mixing
   mixing=0.66d0

!  read(89,*) niter_m
!  read(89,*) iwfn
  iwfn=1
  niter_m=1

!  read(89,*) N_updateU
!  read(89,*) idc
!  read(89,*) alpha_all
!  read(89,*) isph

  N_updateU=0
  idc=0
  alpha_all=1.d0
  isph=1
  iLSDAU=1
  iperturb=0
  iAFM=0 
  iNM=0
  iFM=0
!  crys%net_charge = 0.d0

  if(iLSDAU.eq.0) allocate(bandoccupLSDA(nband1,kpoints%nrk,2))


  read(89,*) idmat
  read(89,*,end=101) isph
!  read(89,*,end=101) crys%net_charge
!  read(89,*,end=101) proj_dos

!  read(89,*) iLSDAU
!  read(89,*) iperturb

! do wavefunction symmetrization for AFM phase

!  read(89,*,end=101) iAFM
!  if (iAFM.ne.0) then
!      read(89,*) iup,idown
!  else
!      read(89,*,end=101) 
!  end if
!
!-----------------
! dangeous options
!  read(89,*,end=101) iNM
!  if (iNM.ge.1) then
!      read(89,*) (inms(ii),ii=1,iNM)
!  else
!      read(89,*,end=101) 
!  end if
!
!  read(89,*,end=101) iFM,(nfm(ii),ii=1,iFM)
!  if (iFM.ge.1) then
!      do kk=1,iFM
!      read(89,*) (ifms(ii,kk),ii=1,nfm(ii))
!      end do
!  else
!      read(89,*,end=101) 
!  end if
!-----------------

  

 101 continue
  close(89)

  nband1=min(mb,nband1)


  if(myproc.eq.0) then

    write(240,*) "========================================================"
    write(240,*) "          LSDA+U local charge density analysis          "
    write(240,*) "========================================================"
    write(240,*)
    write(240,*) "input parameters"
    write(240,*) "----------------"
    write(240,*)
    write(240,*) "Start U correction from",miter,"iterations."
    write(240,*) "Mixing coeff.",mixing
    write(240,*) "Restore LSDA wavefunctions?",iwfn
    write(240,*) "Update U one more time?",N_updateU
    write(240,*) "Number of low  energy bands not involved in the LSDA+U model",nband0
    write(240,*) "Number of high energy bands not involved in the LSDA+U model",mb-nband1
    if(isph.eq.1) write(240,*) "Fully nonspherical LDA+U scheme"
    if(isph.eq.0) write(240,*) "Spherically averaged LDA+U scheme"
    if(Nssite.gt.0) write(240,*) "Number of s-atomic sites",Nssite
    if(Npsite.gt.0) write(240,*) "Number of p-atomic sites",Npsite
    if(Ndsite.gt.0) write(240,*) "Number of d-atomic sites",Ndsite
    if(Nfsite.gt.0) write(240,*) "Number of f-atomic sites",Nfsite
    write(240,*) "Total number of atomic sites",Nsite
    write(240,*)
    if(.NOT.kpoints%reduced) then
       write(240,*) "Full zone k-point sampling, the density matrix will not be symmetrized."
    else
       write(240,*) "Reduced zone k-point sampling, the density matrix will be symmetrized."
    end if

    if (idmat.eq.0) then
        write(240,*) "No initial density matrix for constructing U potential."
    else
        write(240,*) "Use input density matrix to construct U potential."
    end if
    write(240,*)

    write(240,*) "--------------------------------------------------------"
    write(240,*) "Atom      U        J    radial integration parameters"
    write(240,*) "--------------------------------------------------------"

    do ii=1,crys%ntype
       write(240,311) crys%nameat(ii)(1:2),HU_tem(ii),SJ_tem(ii),ngrid1(ii)
    end do
    write(240,*)

 311  format(1x,a2,2x,2f10.5,5x,i5)

!    write(240,*) "--------------------------------------------------------"
!    write(240,*) "            Occupation perturbations"
!    write(240,*) "--------------------------------------------------------"
!   
!    if(iperturb.eq.0) then
!       lambda=0.d0
!       write(240,*) "No occupation perturbation"
!    end if
!
!    if(iperturb.eq.1) then
!       write(240,*) "Do perturbed SCF"
!    end if
!    if(iperturb.eq.2) then
!       write(240,*) "Do perturbation after SCF"
!    end if
!
!
!    do ii=1,Nsite
!       write(240,411) ii,((lambda(kk,is,ii),kk=1,2*lproj(ii)+1),is=1,2)
!    end do

 411  format(1x,i5,1x,10f10.5)

  end if

  ndim=nband1-nband0
  HU=HU/13.60569172d0
  SJ=SJ/13.60569172d0
  Etot_corr=0.d0
  Etot_corr_old=0.d0

! loop over atomic types

  do kk=1,crys%ntype
     if((lproj_tem(kk).ge.0).and.(proj_flag(kk).eq.1)) then
        symbol=crys%nameat(kk)(1:2)
        if(lproj_tem(kk).eq.0) symbol_tem="s"
        if(lproj_tem(kk).eq.1) symbol_tem="p"
        if(lproj_tem(kk).eq.2) symbol_tem="d"
        if(lproj_tem(kk).eq.3) symbol_tem="f"

        if(symbol(2:2)==' ') then
           open(89,file=symbol(1:1)//"_"//symbol_tem//"_orb.dat")
        else
           open(89,file=symbol(1:1)//symbol(2:2)//"_"//symbol_tem//"_orb.dat")
        end if

        read(89,*) ngrid2
        allocate(rr(ngrid2))
        allocate(Rnl(ngrid2))

        do ii=1,ngrid2
           read(89,*) rr(ii),Rnl(ii)
        end do

! interpolate onto regular grid

        allocate(YP(ngrid2))
        allocate(YPP(ngrid2))
        allocate(YPPP(ngrid2))

        call SPCOEF(ngrid2,rr,Rnl,YP,YPP,YPPP,ier)

        do ii=1,ngrid1(kk)
           ind=LEFT(ngrid2,rr,ii*step,ier)
           delta=ii*step-rr(ind)
           Rnl_intp(ii,kk)=Rnl(ind)+(YP(ind)+(YPP(ind)+YPPP(ind)*delta)*delta)*delta
        end do

        deallocate(rr)
        deallocate(Rnl)
        deallocate(YP)
        deallocate(YPP)
        deallocate(YPPP)
        close(89)

! renormalize atomic orbital
        if(inorm(kk).eq.1) then
          step2=step*step
          norm=0.d0
          allocate(Rnl(ngrid1(kk)))
          do ii=1,ngrid1(kk)
             Rnl(ii)=Rnl_intp(ii,kk)*Rnl_intp(ii,kk)*step2
          end do
 

       do ii=1,ngrid1(kk)-7,7
          r1=ii*ii
          r2=(ii+1)*(ii+1)
          r3=(ii+2)*(ii+2)
          r4=(ii+3)*(ii+3)
          r5=(ii+4)*(ii+4)
          r6=(ii+5)*(ii+5)
          r7=(ii+6)*(ii+6)
          r8=(ii+7)*(ii+7)

          norm=norm+    &
          751*(r1*Rnl(ii)+ r8*Rnl(ii+7))+    &
          3577*(r2*Rnl(ii+1)+r7*Rnl(ii+6))+    &
          1323*(r3*Rnl(ii+2)+r6*Rnl(ii+5))+    &
          2989*(r4*Rnl(ii+3)+r5*Rnl(ii+4))
       end do

       norm=norm*step*(7.d0/17280.d0)
       if(myproc.eq.0) write(240,*) crys%nameat(kk), "atomic orbital renormalization",norm

! this seems to be wrong. No sqrt is needed
!       norm=dsqrt(norm)
    
       Rnl_intp(:,kk)=Rnl_intp(:,kk)/norm
       deallocate(Rnl)
       end if
     end if
  end do


!------------------------
  allocate(deltaE(ndim,nrk,2))
  allocate(deltaEU(ndim,ndim,nrk,2))

  allocate(Hmat(ndim,ndim,nrk,2))
  allocate(cnk_LSDA(nlmax,ndim,nrk,2,Nsite))

  allocate(drcnk_LSDA(nlmax,ndim,nrk,2,Nsite,3))
  allocate(force_U(Nsite,3))
  allocate(drcnk_LSDAU(nlmax,ndim,nrk,2,Nsite,3))

!  force_U=0.d0
  drcnk_LSDA=(0.d0,0.d0)
  drcnk_LSDAU=(0.d0,0.d0)
   
  allocate(cnk_LSDAd(nlmax,ndim,nrk,2,Nsite))
  allocate(cnk_LSDAU(nlmax,ndim,nrk,2,Nsite))
  allocate(cnk_LSDAUd(nlmax,ndim,nrk,2,Nsite))
  allocate(cnk_LSDAUd2(nlmax,ndim,nrk,2,Nsite))
  allocate(dmat_LSDAU(nlmax,nlmax,2,Nsite))
  allocate(drdmat_LSDAU(nlmax,nlmax,2,Nsite,3))
  allocate(dmat_LSDA(nlmax,nlmax,2,Nsite))
  allocate(dmat_LSDAUd(nlmax,nlmax,2,Nsite))
  allocate(dmat_LSDAd(nlmax,nlmax,2,Nsite))
  allocate(dmat_old(nlmax,nlmax,2,Nsite))
  allocate(dmat_old2(nlmax,nlmax,2,Nsite))
  allocate(occup_new(ndim,nrk,2,Nsite))
  allocate(occup_old(ndim,nrk,2,Nsite))

!  if(syms%nfrac.ne.0) then
  allocate(cnk_LSDA_f(nlmax,ndim,nrk,2,Nsite,syms%nfrac))
  allocate(cnk_LSDAU_f(nlmax,ndim,nrk,2,Nsite,syms%nfrac))

!  end if

  occup_old=0.d0

  
  if(myproc.eq.0)  write(240,*) "imode,iLSDAU,idmat",imode,iLSDAU,idmat
  if (iLSDAU.eq.1) then
  if(((imode.eq.0).or.(imode.eq.2)).and.(idmat.eq.0)) then
     if(myproc.eq.0) write(240,*) "Need density matrix for band structure calculation"
     stop
  end if
  end if

  if(idmat.eq.1) then
  open(89,file="dmat.old")

  do ia=1,Nsite
     nl=2*lproj(ia)+1
     do is=1,crys%nspin
        read(89,*)
        do im=1,nl
           read(89,*) (dmat_LSDAU(im,jm,is,ia),jm=1,nl)
        end do
     end do
  end do

  dmat_old=dmat_LSDAU
  dmat_old2=dmat_old

  close(89)
  end if
  if(myproc.eq.0) call myflush(240)

 111  format(10f16.12)

  if(myproc.eq.0) open(89,file="dmat.sav")

!-----------------------------------------------------------------------
! construct coefficients for e-e interaction integrals, for s, p and d electrons
! only. f electron will be implemented later. 
!-----------------------------------------------------------------------
   call construct_a(a20,a22,a24,a10,a12,a00)
   allocate(Vmat(nlmax,nlmax,nlmax,nlmax,Nsite))
   call construct_Vm(Vmat,a20,a22,a24,a10,a12,a00,HU,SJ,Nsite,nlmax,lproj)
!------------------------------------------------------------


   return
   end subroutine lsdau_init

!------------------------------------------------------------
      SUBROUTINE SPCOEF(N,X,F,B,C,D,FLAG)

      IMPLICIT NONE
      integer,intent(IN)::N
      integer,intent(OUT)::FLAG
      DOUBLE PRECISION,intent(IN):: X(N),F(N)
      double precision,intent(OUT):: B(N),C(N),D(N)
!
!  Input parameters:
!    N    = number of data points.
!    X    = vector of values of the independent variable ordered
!           SO THAT  x(i) < x(i+1)  for all I.
!    F    = vector of values of the dependent variable.
!  Output parameters:
!    B    = vector of S'(X(I)) values.
!    C    = vector of S"(X(I))/2 values.
!    D    = vector of S'''(X(I)+)/6 values (I .LT. N).
!    FLAG =  0  normal return;
!         = -1  input N .LE. 1;
!         = -2  X vector is incorrectly ordered.
!
!  The vectors X, F, B, C, D must be dimensioned at least N in the
!  calling program.
!
!  Local variables:
      INTEGER I,K
      DOUBLE PRECISION FP1,FPN,H,P,THREE,TWO,ZERO
      DATA THREE,TWO,ZERO/3.D0,2.D0,0.D0/
!
      IF (N .LE. 1) THEN
         FLAG = -1
         RETURN
      ENDIF
!
!     Calculate coefficients for the tri-diagonal system: store
!     sub-diagonal in B, diagonal in D, difference quotient in C.
!
      B(1) = X(2)-X(1)
      IF (B(1) .LE. ZERO) THEN
         FLAG = -2
         RETURN
      ENDIF
      C(1) = (F(2)-F(1))/B(1)
      IF (N .EQ. 2) THEN
         B(1) = C(1)
         C(1) = ZERO
         D(1) = ZERO
         B(2) = B(1)
         C(2) = ZERO
         FLAG = 0
         RETURN
      ENDIF
      D(1) = TWO*B(1)
      DO 20 I = 2,N-1
         B(I) = X(I+1)-X(I)
         IF (B(I) .LE. ZERO) THEN
            FLAG = -2
            RETURN
         ENDIF
         C(I) = (F(I+1)-F(I))/B(I)
         D(I) = TWO*(B(I)+B(I-1))
   20 CONTINUE
      D(N) = TWO*B(N-1)
!
!     Calculate estimates for the end slopes.  Use polynomials
!     interpolating data nearest the end.
!
      FP1 = C(1)-B(1)*(C(2)-C(1))/(B(1)+B(2))
      IF (N .GT. 3) FP1 = FP1+B(1)*((B(1)+B(2))*(C(3)-C(2)) &
       /(B(2)+B(3))-C(2)+C(1))/(X(4)-X(1))
      FPN = C(N-1)+B(N-1)*(C(N-1)-C(N-2))/(B(N-2)+B(N-1))
      IF (N .GT. 3) FPN = FPN+B(N-1)*(C(N-1)-C(N-2)-(B(N-2) &
       +B(N-1))*(C(N-2)-C(N-3))/(B(N-2)+B(N-3)))/(X(N)-X(N-3))
!
!     Calculate the right-hand-side and store it in C.
!
      C(N) = THREE*(FPN-C(N-1))
      DO 30 K = 2,N-1
         I = N-K+1
         C(I) = THREE*(C(I)-C(I-1))
   30 CONTINUE
      C(1) = THREE*(C(1)-FP1)
!
!     Solve the tridiagonal system.
!
      DO 40 K = 2,N
         P = B(K-1)/D(K-1)
         D(K) = D(K)-P*B(K-1)
         C(K) = C(K)-P*C(K-1)
   40 CONTINUE
      C(N) = C(N)/D(N)
      DO 50 K = N-1,1,-1
         C(K) = (C(K)-B(K)*C(K+1))/D(K)
   50 CONTINUE
!
!     Calculate the coefficients defining the spline.
!
      DO 60 I = 1,N-1
         H = X(I+1)-X(I)
         D(I) = (C(I+1)-C(I))/(THREE*H)
         B(I) = (F(I+1)-F(I))/H-H*(C(I)+H*D(I))
   60 CONTINUE
      B(N) = B(N-1)+H*(TWO*C(N-1)+H*THREE*D(N-1))
      FLAG = 0
      RETURN
      END
!----------------------------------------------------------------
      INTEGER FUNCTION LEFT(LXT,XT,X,MFLAG)
!   LEFT    FINDS INDEX LEFT OF AN ARRAY XT FOR WHICH XT(LEFT)
!           LIES IMMEDIATELY LEFT OF X

!   PURPOSE:
!           FINDS INDEX LEFT OF AN ARRAY XT FOR WHICH XT(LEFT)
!           LIES IMMEDIATELY LEFT OF X

!   INPUT ARGUMENTS:
!     LXT   : NUMBER OF ELEMENTS IN VECTOR XT
!     XT    : VECTOR OF LENGTH LXT STORING THE ABSCISSAE
!     X     : X-VALUE FOR WHICH THE INDEX LEFT IS TO BE FOUND

!   OUTPUT ARGUMENTS:
!     LEFT  : INDEX FOR WHICH XT(LEFT) LIES IMMEDIATELY LEFT OF X
!     MFLAG : FLAG SET IN THE FOLLOWING MANNER
!             LEFT  MFLAG
!              1     -1     IF               X .LT. XT(1)
!              I      0     IF  XT(I)   .LE. X .LT. XT(I+1)
!             LXT     1     IF  XT(LXT) .LE. X

!   METHOD:
!     THAT OF CARL DE BOOR AS DESCRIBED ON PAGE 91 FF. IN:
!     /1/ DE BOOR,C. (1978) A PRACTICAL GUIDE TO SPLINES.
!         APPLIED MATHEMATICAL SCIENCES, VOLUME 27.
!         NEW-YORK-HEIDELBERG-BERLIN: SPRINGER.

!   IMPLEMENTED BY:
!      KRAFT,D., DLR - INSTITUT FUER DYNAMIK DER FLUGSYSTEME
!                D-8031 OBERPFAFFENHOFEN

!   STATUS: 15. JANUARY 1980

!   SUBROUTINES REQUIRED: NONE

      INTEGER LXT,MFLAG,IHI,ILO,ISTEP,MIDDLE
      DOUBLE PRECISION X,XT(LXT)
      SAVE ILO
      DATA ILO/1/

      IHI=ILO+1
      IF(IHI.LT.LXT)                   GOTO  10
      IF(X.GE.XT(LXT))                 GOTO 110
      IF(LXT.LE.1)                     GOTO  90
      ILO=LXT-1
      IHI=LXT
   10 IF(X.GE.XT(IHI))                 GOTO  40
      IF(X.GE.XT(ILO))                 GOTO 100
      ISTEP=1
   20 IHI=ILO
      ILO=IHI-ISTEP
      IF(ILO.LE.1)                     GOTO  30
      IF(X.GE.XT(ILO))                 GOTO  70
      ISTEP=ISTEP+ISTEP
                                       GOTO  20
   30 ILO=1
      IF(X.LT.XT(1))                   GOTO  90
                                       GOTO  70
   40 ISTEP=1
   50 ILO=IHI
      IHI=ILO+ISTEP
      IF(IHI.GE.LXT)                   GOTO  60
      IF(X.LT.XT(IHI))                 GOTO  70
      ISTEP=ISTEP+ISTEP
                                       GOTO  50
   60 IF(X.GE.XT(LXT))                 GOTO 110
      IHI=LXT
   70 MIDDLE=(ILO+IHI)/2
      IF(MIDDLE.EQ.ILO)                GOTO 100
      IF(X.LT.XT(MIDDLE))              GOTO  80
      ILO=MIDDLE
                                       GOTO  70
   80 IHI=MIDDLE
                                       GOTO  70
   90 MFLAG=-1
      LEFT=1
                                       GOTO 120
  100 MFLAG=0
      LEFT=ILO
                                       GOTO 120
  110 MFLAG=1
      LEFT=LXT
  120                                  RETURN

!   END OF LEFT

      END

!--------------------------------------------------------------
      subroutine construct_a(a20,a22,a24,a10,a12,a00)

      implicit none
      double precision::a20(5,5,5,5),a22(5,5,5,5),a24(5,5,5,5)
      double precision::a10(3,3,3,3),a12(3,3,3,3),a00(1,1,1,1)

      double precision::a_coef

      integer::m1,m2,m3,m4
      integer::n1,n2,n3,n4

!      write(9,*) "L=2"
      do m1=1,5
         n1=m1-3
         do m2=1,5
            n2=m2-3
            do m3=1,5
               n3=m3-3
               do m4=1,5
                  n4=m4-3
                  a20(m1,m2,m3,m4)=a_coef( 2, 0, n1, n2, n3, n4 )
                  a22(m1,m2,m3,m4)=a_coef( 2, 2, n1, n2, n3, n4 )
                  a24(m1,m2,m3,m4)=a_coef( 2, 4, n1, n2, n3, n4 )
!                  if(abs(a20(m1,m2,m3,m4)).gt.1.d-5) write(9,111) n1,n2,n3,n4,0, a20(m1,m2,m3,m4)
!                  if(abs(a22(m1,m2,m3,m4)).gt.1.d-5) write(9,111) n1,n2,n3,n4,2, a22(m1,m2,m3,m4)
!                  if(abs(a24(m1,m2,m3,m4)).gt.1.d-5) write(9,111) n1,n2,n3,n4,4, a24(m1,m2,m3,m4)
               end do
            end do
         end do
      end do

!      write(9,*) "L=1"
      do m1=1,3
         n1=m1-2
         do m2=1,3
            n2=m2-2
            do m3=1,3
               n3=m3-2
               do m4=1,3
                  n4=m4-2
                  a10(m1,m2,m3,m4)=a_coef( 1, 0, n1, n2, n3, n4 )
                  a12(m1,m2,m3,m4)=a_coef( 1, 2, n1, n2, n3, n4 )
!                  if(abs(a10(m1,m2,m3,m4)).gt.1.d-5) write(9,111) n1,n2,n3,n4,0, a10(m1,m2,m3,m4)
!                  if(abs(a12(m1,m2,m3,m4)).gt.1.d-5) write(9,111) n1,n2,n3,n4,2, a12(m1,m2,m3,m4)
               end do
            end do
         end do
      end do

      a00(1,1,1,1)=a_coef( 0, 0, 0, 0, 0, 0 )
!      write(9,*) "L=1"
!      write(9,111) 0,0,0,0,0,a00(1,1,1,1)

 111  format(5i3,f10.5)

      return
     end subroutine construct_a


     subroutine construct_Vm(Vmat,a20,a22,a24,a10,a12,a00,HU,SJ,Nsite,nlmax,lproj)

     implicit none
     integer::Nsite,nlmax,lproj(Nsite)

     double precision::Vmat(nlmax,nlmax,nlmax,nlmax,Nsite)
     double precision::a10(3,3,3,3),a12(3,3,3,3),a20(5,5,5,5),a22(5,5,5,5),a24(5,5,5,5),a00(1,1,1,1)
     double precision::HU(Nsite),SJ(Nsite)

     integer::ia,m1,m2,m3,m4
     double precision::F0,F2,F4


! IMPORTANT:
!
! Vmat(m1,m2,m3,m4)=<m1,m3|V|m2,m4>
!                  =a0(m1,m2,m3,m4)*F0+
!                   a2(m1,m2,m3,m4)*F2+
!                   a4(m1,m2,m3,m4)*F4
!
! see A. Liechtenstein, et al, PRB 52,R5467 (1995)

 
     do ia=1,Nsite
        if(lproj(ia).eq.2) then
        F0=HU(ia)
        F2=14.d0*SJ(ia)/1.625d0
        F4=F2*0.625d0

        do m1=1,5
        do m2=1,5
        do m3=1,5
        do m4=1,5
           Vmat(m1,m2,m3,m4,ia)=a20(m1,m2,m3,m4)*F0+a22(m1,m2,m3,m4)*F2+a24(m1,m2,m3,m4)*F4
        end do 
        end do
        end do
        end do
        end if

        if(lproj(ia).eq.1) then
        F0=HU(ia)
        F2=5.d0*SJ(ia)

        do m1=1,3
        do m2=1,3
        do m3=1,3
        do m4=1,3
           Vmat(m1,m2,m3,m4,ia)=a10(m1,m2,m3,m4)*F0+a12(m1,m2,m3,m4)*F2
        end do 
        end do
        end do
        end do
        end if
        if(lproj(ia).eq.0) then
        F0=HU(ia)
        Vmat(1,1,1,1,ia)=a00(1,1,1,1)*F0
        end if


     end do
     return
     end subroutine construct_Vm



!-----------------------------------------------------
    subroutine lsdau_AFM_symm(cd_gs,cd_data,crys)

  use lsdau_shared

  include 'use.h'
  implicit none               ! never remove this line.

  integer::i
  double precision::tau(3)
  double precision::phase

  type(crystal) :: crys
  type(parallel_gspace) :: cd_gs
  complex(dp) :: cd_data(cd_gs%length, 2)

  if((iAFM.eq.1).or.(iAFM.eq.3)) then
    if(myproc.eq.0) write(9,*) "AFM symmetrization of the charge density" 
    tau(:)=crys%rat(:,1,iup)-crys%rat(:,1,idown)

  do i = 1, cd_gs%length
     phase=cd_gs%gvec(1,i)*tau(1)+ &
           cd_gs%gvec(2,i)*tau(2)+ &
           cd_gs%gvec(3,i)*tau(3)
           cd_data(i,1)=cd_data(i,2)*exp(cmplx(0.0d0,phase,kind=8))
  end do
  end if

  if(iAFM.eq.4) then
    if(myproc.eq.0) write(9,*) "NM symmetrization of the charge density" 

  do i = 1, cd_gs%length
     cd_data(i,1)=(cd_data(i,2)+cd_data(i,1))/2.d0
     cd_data(i,2)=cd_data(i,1)

  end do
  end if

  return
  end 

!------------------------------------------
   subroutine lsdau_symm2(dmat)

  use lsdau_shared
  implicit none               ! never remove this line.

  integer::i,j
  double complex::dmat(nlmax,nlmax,2,Nsite)


  if((iAFM.eq.2).or.(iAFM.eq.3)) then
  if(myproc.eq.0) write(9,*) "AFM symmetrization of the density matrix",iup,idown

  dmat(:,:,1,iup)=(dmat(:,:,2,idown)+dmat(:,:,1,iup))/2.d0
  dmat(:,:,2,idown)=dmat(:,:,1,iup)

  dmat(:,:,1,idown)=(dmat(:,:,2,iup)+dmat(:,:,1,idown))/2.d0
  dmat(:,:,2,iup)=dmat(:,:,1,idown)
  end if

  if(iNM.ge.1) then
  if(myproc.eq.0) write(9,*) "symmetrization of the density matrix of nonmag. sites"

  do i=1,iNM
     if(myproc.eq.0) write(9,*) "Site",inms(i)
     dmat(:,:,1,inms(i))=(dmat(:,:,2,inms(i))+dmat(:,:,1,inms(i)))/2.d0
     dmat(:,:,2,inms(i))=dmat(:,:,1,inms(i))
  end do

  end if

  if(iFM.ge.1) then
  if(myproc.eq.0) write(9,*) "symmetrization of the density matrix of ferromag. sites"

  do i=1,iFM
     do j=2,nfm(i)
     if(myproc.eq.0) write(9,*) "Site",ifms(1,i),ifms(j,i)

     dmat(:,:,:,ifms(j,i))=dmat(:,:,:,ifms(1,i))
     end do
  end do

  end if


  return
  end
