    module lsdau_shared


  integer:: myproc,nproc,proj_dos,init_flag
  integer:: iAFM,iup,idown,iNM,inms(200),iFM,ifms(200,200),nfm(200)
  integer:: iwfn,idc,iLSDAU,isph,idmat,iperturb,interC
  integer:: Nsite,Nssite,Npsite,Ndsite,Nfsite,nl,nlmax,nlmax_alt
  integer:: Nsite_alt,Nssite_alt,Npsite_alt,Ndsite_alt,Nfsite_alt
  integer:: miter,N_updateU,nband1,nband0,ndim,ndim_alt,niter_m,nband0_alt,nband1_alt
  integer::ind_site_atom(500), ind_site_isite(500),ind_isite_site(5,500)


  double precision:: mixing
  double precision:: Etot_corr,ELSDA_corr,Etot_corr_old,Eband_diff,alpha_all,step

  double precision::a20(5,5,5,5),a22(5,5,5,5),a24(5,5,5,5)
  double precision::a10(3,3,3,3),a12(3,3,3,3)
  double precision::a00(1,1,1,1),ylm0

  integer,allocatable::lproj(:),lproj_tem(:),ngrid1(:),ngrid(:),proj_flag(:),inorm(:),mproj(:),nproj(:,:)
  integer,allocatable::lproj_alt(:),lproj_tem_alt(:),ngrid1_alt(:),proj_flag_alt(:),inorm_alt(:)
  integer,allocatable::ntp(:)

  double precision,allocatable::HU(:),SJ(:),HU_tem(:),SJ_tem(:)
  double precision,allocatable::HU_alt(:),HU_tem_alt(:),SJ_tem_alt(:),SJ_alt(:)
  double precision,allocatable::Rnl(:),Rnl_intp(:,:),Rnl_intp_alt(:,:),step_alt(:)
  double precision,allocatable::Vmat(:,:,:,:,:)
  double precision,allocatable::lambda(:,:,:)
  double precision,allocatable::deltaE(:,:,:)
  double precision,allocatable::deltaE_alt(:,:,:)
  double precision,allocatable::deltaEU(:,:,:,:)
  double precision,allocatable::bandoccupLSDA(:,:,:)
  double precision,allocatable::occup_new(:,:,:,:) ! nband,nrk,ns,Nsite
  double precision,allocatable::occup_old(:,:,:,:) ! nband,nrk,ns,Nsite
  double precision,allocatable::cnk_LSDAUd2(:,:,:,:,:) ! m,nb,nk,ns,nsite
  double precision,allocatable::sb2(:,:,:),sb1(:,:,:),sb0(:,:,:)

  double complex,allocatable::ylm2(:,:,:),ylm1(:,:,:),phase_ikgr0(:,:,:),ikg(:,:,:)
  double complex,allocatable::Hmat(:,:,:,:) ! nb,nb,nk,ns
  double complex,allocatable::cnk_LSDA(:,:,:,:,:) ! m,nb,nk,ns,nsite
  double complex,allocatable::cnk_LSDAd(:,:,:,:,:) ! m,nb,nk,ns,nsite
  double complex,allocatable::cnk_LSDA_alt(:,:,:,:,:) ! m,nb,nk,ns,nsite

! derivative of projection coefficients
  double precision,allocatable::force_U(:,:) ! Nsite,pol(x,y,z)
  double complex,allocatable::drcnk_LSDA(:,:,:,:,:,:) ! m,nb,nk,ns,nsite,pol(x,y,z)
  double complex,allocatable::drcnk_LSDAU(:,:,:,:,:,:) ! m,nb,nk,ns,nsite,pol(x,y,z)
  double complex,allocatable::drdmat_LSDAU(:,:,:,:,:) ! m,m',ns,nsite,3

! for systems with fractional translations


  double complex,allocatable::cnk_LSDA_f(:,:,:,:,:,:) ! m,nb,nk,ns,nsite,nfrac
  double complex,allocatable::cnk_LSDAU_f(:,:,:,:,:,:) ! m,nb,nk,ns,nsite,nfrac
  double complex,allocatable::cnk_LSDA_alt_f(:,:,:,:,:,:) ! m,nb,nk,ns,nsite


  double complex,allocatable::cnk_LSDAU(:,:,:,:,:) ! m,nb,nk,ns,nsite
  double complex,allocatable::cnk_LSDAUd(:,:,:,:,:) ! m,nb,nk,ns,nsite
  double complex,allocatable::dmat_LSDAU(:,:,:,:) ! m,m',ns,nsite
  double complex,allocatable::dmat_LSDAUd(:,:,:,:) ! m,m',ns,nsite
  double complex,allocatable::dmat_LSDA(:,:,:,:) ! m,m',ns,nsite
  double complex,allocatable::dmat_LSDAd(:,:,:,:) ! m,m',ns,nsite
  double complex,allocatable::dmat_old(:,:,:,:) ! m,m',ns,nsite
  double complex,allocatable::dmat_old2(:,:,:,:) ! m,m',ns,nsite

  end module lsdau_shared
