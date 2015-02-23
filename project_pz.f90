

  subroutine project_pz(myproc,nproc, crys, syms,wavefn, kpoints, k_gspace,bands,energs,iflag)
                                                                                                         
  include 'use.h'
  implicit none

!-------------------------------------------------------------------------
  interface

!        subroutine angular_wfnpz_d2r(wavefn,nkpt,kxyz,tau,rc,cart,  &
!                          kkx,kky,kkz,vcell,myproc,ngrid,step,Rl,cnk,Rlm,proj_flag)
!
!        implicit none
!
!        integer ngrid,myproc,proj_flag,nkpt,kxyz(3,nkpt)
!        double precision tau(3),rc,kkx,kky,kkz,cart(3,3),step,vcell,Rl(ngrid)
!        double complex:: wavefn(nkpt),Rlm(ngrid,5),cnk(5)
!
!        end subroutine angular_wfnpz_d2r
!
!        subroutine angular_wfnpz_p2r(wavefn,nkpt,kxyz,tau,rc,cart,  &
!                          kkx,kky,kkz,vcell,myproc,ngrid,step,Rl,cnk,Rlm,proj_flag)
!
!        implicit none
!
!        integer ngrid,myproc,proj_flag,nkpt,kxyz(3,nkpt)
!        double precision tau(3),rc,kkx,kky,kkz,cart(3,3),step,vcell,Rl(ngrid)
!        double complex:: wavefn(nkpt),Rlm(ngrid,3),cnk(3)
!
!        end subroutine angular_wfnpz_p2r
!
!        subroutine angular_wfnpz_s2(wavefn,nkpt,kxyz,tau,rc,cart,  &
!                          kkx,kky,kkz,vcell,myproc,ngrid,step,Rl,cnk,Rlm,proj_flag)
!
!        implicit none
!
!        integer ngrid,myproc,proj_flag,nkpt,kxyz(3,nkpt)
!        double precision tau(3),rc,kkx,kky,kkz,cart(3,3),step,vcell,Rl(ngrid)
!        double complex:: wavefn(nkpt),Rlm(ngrid,1),cnk(1)
!
!        end subroutine angular_wfnpz_s2


        subroutine project_dos(myproc, nproc, cnk_all, syms, crys, bands, energs, kp,  &
                      nband0,nbandmin,nsites,lproj1,mproj,lproj,max_proj,info,nproj, &
                      ind_site_atom,ind_site_isite,ind_isite_site)



         include 'use.h'
         implicit none

         type(crystal), intent(in) :: crys                   ! the crystal structure
         type(band), intent(in) :: bands
         type(energy), intent(in) :: energs
         type(kpoint), intent(in) :: kp
         type(symmetry), intent(in) :: syms
         integer, intent(in) :: nsites,max_proj,lproj1(nsites),mproj(crys%ntype),lproj(max_proj,crys%ntype)
         integer, intent(in) :: myproc, nproc, nbandmin,nproj(max_proj,crys%ntype),info,nband0
         integer, intent(in) :: ind_site_atom(500),ind_site_isite(500),ind_isite_site(5,500)
         double complex, intent(in) :: cnk_all(5, nbandmin,bands%nrk,crys%nspin,nsites)
         end subroutine project_dos

  end interface
!-------------------------------------------------------------------------


  integer myproc,nproc,iflag,lprojtem(100)
  type(crystal) :: crys       ! crystal structure
  type(symmetry) :: syms
  type(complex_gspace_array), intent(in) ::  wavefn
  type(kpoint), intent(in) :: kpoints     ! BZ integration kpoints
  type(parallel_gspace), intent(in) :: k_gspace(*)
  type(band), intent(in) :: bands 
  type(energy), intent(in) :: energs     ! energies and other information

!--------------------------------------------------------------

  character*1 symbol1,symbol2,symbol_tem
  character (LEN=2):: symbol
  integer::ii,jj,kk,ll,ier,ind,indd(1),ngrid2,is,irk,ib,nbandmin,nsites,isite,jsite,i,j,max_proj
  integer::ind_site_atom(500),ind_site_isite(500),ind_isite_site(5,500),iatom,k
  integer::max_L,max_LM
  integer::LEFT
  integer,allocatable::lproj1(:)

! maximum projection per atom,max L
!
  parameter(max_proj=5,max_L=2,max_LM=max_L*2+1)
  

  double precision:: delta,step2,rc
  double precision:: r1,r2,r3,r4,r5,r6,r7,r8,norm,tau(3),kkx,kky,kkz
  double complex::cnk(max_LM)

  integer,allocatable::ngrid(:,:),mproj(:),nproj(:,:),lproj(:,:),proj_flag(:,:),inorm(:,:)
  double precision,allocatable::step(:,:),Rnl(:),Rnl_intp(:,:,:),anorm(:),dnk_all(:,:,:,:,:)
  double precision,allocatable:: rr(:),YP(:),YPP(:),YPPP(:)
  double complex,allocatable::Rlm(:,:),cnk_all_ibz(:,:,:,:,:)


  open(unit=339,file='project.in')

! loop over atomic types

  allocate(mproj(crys%ntype))
  allocate(lproj(max_proj,crys%ntype))
  allocate(inorm(max_proj,crys%ntype))
  allocate(nproj(max_proj,crys%ntype))
  allocate(ngrid(max_proj,crys%ntype))
  allocate(proj_flag(max_proj,crys%ntype))
  allocate(step(max_proj,crys%ntype))
  lproj=0
  mproj=0
  nproj=0
  ngrid=0
  step=0.0
  proj_flag=1

! mproj: how many projectors on this site
! nproj: principle quantum of the projection

  isite=0
  do ii=1,crys%ntype
     read(339,*) mproj(ii)
     if(mproj(ii).gt.max_proj) then
     write(9,*) "please increase max_proj"
     call myflush(9)
     call mystop()
     end if

     if(mproj(ii).ge.1) then
        do jj=1, mproj(ii)
           read(339,*)nproj(jj,ii),lproj(jj,ii),step(jj,ii),ngrid(jj,ii),proj_flag(jj,ii),inorm(jj,ii)
           lprojtem(isite+1:isite+crys%natom(ii))=lproj(jj,ii)
           isite=isite+crys%natom(ii)
        end do
     end if
  end do


  close(339)

  nbandmin=MIN(bands%min(1),bands%min(2))

  allocate(Rnl_intp(MAXVAL(ngrid),MAXVAL(mproj),crys%ntype))
  allocate(Rlm(MAXVAL(ngrid),MAXVAL(lproj)*2+1))


  nsites=sum(mproj)
  allocate(anorm(nsites))

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
!        write(9,*) "isite,iatom,ind_isite_site",k,iatom+j,isite

     end do
     end do

     end if

     iatom=iatom+crys%natom(ii)

  end do
  nsites=isite

!  do i=1,nsites
!     write(9,*) "site,atom,ind_site_isite",i,ind_site_atom(i),ind_site_isite(i)
!  end do


  allocate(lproj1(nsites))
  lproj1(1:nsites)=lprojtem(1:nsites)

  allocate(cnk_all_ibz(max_LM,nbandmin,kpoints%nrk,crys%nspin,nsites))
  allocate(dnk_all(max_LM,nbandmin,kpoints%nrk,crys%nspin,nsites))

  cnk_all_ibz=(0.d0,0.d0)

  isite=0
  do ii=1,crys%ntype

     symbol=crys%nameat(ii)(1:2)
     if(mproj(ii).ge.1) then
     do jj=1,mproj(ii)
        isite=isite+1
        if(proj_flag(jj,ii).eq.1) then

        if(lproj(jj,ii).eq.0) symbol_tem="s"
        if(lproj(jj,ii).eq.1) symbol_tem="p"
        if(lproj(jj,ii).eq.2) symbol_tem="d"
        if(lproj(jj,ii).eq.3) symbol_tem="f"

        if(symbol(2:2)==' ') then
           open(89,file=symbol(1:1)//"_"//CHAR(nproj(jj,ii)+48)//symbol_tem//"_orb.dat")
        else
           open(89,file=symbol(1:1)//symbol(2:2)//"_"//CHAR(nproj(jj,ii)+48)//symbol_tem//"_orb.dat")
        end if

        read(89,*) ngrid2

        allocate(rr(ngrid2))
        allocate(Rnl(ngrid2))


        do kk=1,ngrid2
           read(89,*) rr(kk),Rnl(kk)
        end do
! interpolate onto regular grid

        allocate(YP(ngrid2))
        allocate(YPP(ngrid2))
        allocate(YPPP(ngrid2))

        call SPCOEF(ngrid2,rr,Rnl,YP,YPP,YPPP,ier)

        do kk=1,ngrid(jj,ii)
           ind=LEFT(ngrid2,rr,kk*step(jj,ii),ier)
           delta=kk*step(jj,ii)-rr(ind)
           Rnl_intp(kk,jj,ii)=Rnl(ind)+(YP(ind)+(YPP(ind)+YPPP(ind)*delta)*delta)*delta
        end do

        deallocate(rr)
        deallocate(Rnl)
        deallocate(YP)
        deallocate(YPP)
        deallocate(YPPP)
        close(89)

! normalize atomic orbital
        step2=step(jj,ii)*step(jj,ii)
        norm=0.d0
        allocate(Rnl(ngrid(jj,ii)))

        do kk=1,ngrid(jj,ii)
           Rnl(kk)=Rnl_intp(kk,jj,ii)*Rnl_intp(kk,jj,ii)*step2
        end do

        do kk=1,ngrid(jj,ii)-7,7
           r1=kk*kk
           r2=(kk+1)*(kk+1)
           r3=(kk+2)*(kk+2)
           r4=(kk+3)*(kk+3)
           r5=(kk+4)*(kk+4)
           r6=(kk+5)*(kk+5)
           r7=(kk+6)*(kk+6)
           r8=(kk+7)*(kk+7)
 
           norm=norm+    &
           751*(r1*Rnl(kk)+ r8*Rnl(kk+7))+    &
           3577*(r2*Rnl(kk+1)+r7*Rnl(kk+6))+    &
           1323*(r3*Rnl(kk+2)+r6*Rnl(kk+5))+    &
           2989*(r4*Rnl(kk+3)+r5*Rnl(kk+4))
        end do

        norm=norm*step(jj,ii)*(7.d0/17280.d0)
        write(9,333) crys%nameat(ii),nproj(jj,ii), symbol_tem," normalization",norm,step(jj,ii),ngrid(jj,ii)
        call myflush(9)

        anorm(isite)=norm

        deallocate(Rnl)
        end if
    end do
    end if
  end do
 333 format(4x,a2,i5,a2,a15,2f10.5,i5)

!------------------------------------------------------------------------
! projection of wavefunctions

   isite=0
   do ii=1, crys%ntype
      if(mproj(ii).ge.1) then
      do jj=1,mproj(ii)

      do kk=1, crys%natom(ii)
         isite=isite+1
         tau(:)=matmul(crys%avec,crys%rat(:,kk,ii))/pi2

         do is = 1, crys%nspin
         do irk = 1, kpoints%nrk
            kkx=kpoints%rk(1,irk)
            kky=kpoints%rk(2,irk)
            kkz=kpoints%rk(3,irk)

            do ib = 1, nbandmin

               if(lproj(jj,ii).eq.2) then
               call angular_wfnpz_d2r(wavefn%data((ib-1)*k_gspace(irk)%length+1,irk,is),  &
                    k_gspace(irk)%length,k_gspace(irk)%gvec, tau,rc,crys%bvec,&
                    kkx,kky,kkz,crys%vcell,myproc,ngrid(jj,ii),step(jj,ii),Rnl_intp(:,jj,ii), &
                    cnk,Rlm,proj_flag(jj,ii))
                    
                    do i=1,5
!                       cnk_all(i,ib,irk,is,isite)=cnk(i)*DCONJG(cnk(i))
                       cnk_all_ibz(i,ib,irk,is,isite)=cnk(i)
                    end do

               end if

               if(lproj(jj,ii).eq.1) then

               call angular_wfnpz_p2r(wavefn%data((ib-1)*k_gspace(irk)%length+1,irk,is),  &
                    k_gspace(irk)%length,k_gspace(irk)%gvec, tau,rc,crys%bvec,&
                    kkx,kky,kkz,crys%vcell,myproc,ngrid(jj,ii),step(jj,ii),Rnl_intp(:,jj,ii), &
                    cnk,Rlm,proj_flag(jj,ii))

                    do i=1,3
!                       cnk_all(i,ib,irk,is,isite)=cnk(i)*DCONJG(cnk(i))
                       cnk_all_ibz(i,ib,irk,is,isite)=cnk(i)
                    end do
               end if

               if(lproj(jj,ii).eq.0) then
               call angular_wfnpz_s2(wavefn%data((ib-1)*k_gspace(irk)%length+1,irk,is),  &
                    k_gspace(irk)%length,k_gspace(irk)%gvec, tau,rc,crys%bvec,&
                    kkx,kky,kkz,crys%vcell,myproc,ngrid(jj,ii),step(jj,ii),Rnl_intp(:,jj,ii), &
                    cnk,Rlm,proj_flag(jj,ii))

!                    cnk_all(1,ib,irk,is,isite)=cnk(1)*DCONJG(cnk(1))
                    cnk_all_ibz(1,ib,irk,is,isite)=cnk(1)

               end if

            end do

         end do
         end do

     end do
     end do
     end if

  end do

!---------------------
! write out occupation

! to be implemented


!---------------------

! normalize if asked

  isite=0
  jsite=0
  do ii=1,crys%ntype

     if(mproj(ii).ge.1) then
     do jj=1,mproj(ii)
        isite=isite+1

        if(inorm(jj,ii).eq.1) then

          norm=anorm(isite)
!          cnk_all(:,:,:,:,jsite+1:jsite+crys%natom(ii))=  &
!          cnk_all(:,:,:,:,jsite+1:jsite+crys%natom(ii))/norm

          cnk_all_ibz(:,:,:,:,jsite+1:jsite+crys%natom(ii))=  &
          cnk_all_ibz(:,:,:,:,jsite+1:jsite+crys%natom(ii))/dsqrt(norm)
        end if
        jsite=jsite+crys%natom(ii)

     end do
     end if

  end do


! print out projection results

  dnk_all=cnk_all_ibz*DCONJG(cnk_all_ibz)

  if(myproc.eq.0) then

  isite=0
  do ii=1,crys%ntype

     symbol=crys%nameat(ii)(1:2)
     if(mproj(ii).ge.1) then
     do jj=1,mproj(ii)
       
        if(lproj(jj,ii).eq.0) symbol_tem="s"
        if(lproj(jj,ii).eq.1) symbol_tem="p"
        if(lproj(jj,ii).eq.2) symbol_tem="d"
        if(lproj(jj,ii).eq.3) symbol_tem="f"

        if(iflag.eq.1) then
        if(symbol(2:2)==' ') then
           open(239,file=symbol(1:1)//CHAR(ii+48)//"_"//CHAR(nproj(jj,ii)+48)//symbol_tem//"_proj.dat")
        else
           open(239,file=symbol(1:1)//symbol(2:2)//CHAR(ii+48)//"_"//CHAR(nproj(jj,ii)+48)//symbol_tem//"_proj.dat")
        end if 
        else
        if(symbol(2:2)==' ') then
           open(239,file=symbol(1:1)//CHAR(ii+48)//"_"//CHAR(nproj(jj,ii)+48)//symbol_tem//"_projbs.dat")
        else
           open(239,file=symbol(1:1)//symbol(2:2)//CHAR(ii+48)//"_"//CHAR(nproj(jj,ii)+48)//symbol_tem//"_projbs.dat")
        end if 
        end if



        do is=1,crys%nspin
        if(crys%nspin.eq.2) write(239,*) "spin",is
  
          do ib = 1, nbandmin
          do irk = 1, kpoints%nrk
             write(239,111) irk,bands%energy(ib,irk,is)+energs%vxc0(is), &
                            ((dnk_all(i,ib,irk,is,j),i=1,lproj(jj,ii)*2+1),j=isite+1,isite+crys%natom(ii))
          end do
          write(239,*) 
          end do
        end do
        close(239)
        isite=isite+crys%natom(ii)
     end do
     end if
  end do

  end if


  if(iflag.eq.1) then
  call project_dos(myproc, nproc, cnk_all_ibz, syms,crys, bands, energs, kpoints,  &
              0,nbandmin,nsites,lproj1,mproj,lproj,max_proj,1,nproj, &
                     ind_site_atom,ind_site_isite,ind_isite_site)
  end if


 111 format(i5,100f10.5)

  deallocate(mproj)
  deallocate(nproj)
  deallocate(inorm)
  deallocate(lproj)
  deallocate(ngrid)
  deallocate(proj_flag)
  deallocate(step)
                                                                                                                       
  deallocate(anorm)
  deallocate(Rnl_intp)
  deallocate(Rlm)
  deallocate(dnk_all)
  deallocate(cnk_all_ibz)

  if(myproc.eq.0) close(339)

  end subroutine project_pz

