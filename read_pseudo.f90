subroutine read_pseudo(ipr, ispn, pw_params, gs, energs, pspot, &
     crys, vion, denc, denv, ddc, dvql, vql, dnc,vnl_ave, ffts)
  !     ------------------------------------------------------------------
  use all_to_all_module  
  include 'use.h'  
  implicit none               ! implicit? Just say no!
  include 'interface.h'  
  include 'all_to_all.h'  
  include 'flibcalls.ph' 
  !
  !     INPUT:
  !     -----
  !
  type(pw_parameter), intent(in) :: pw_params  
  integer, intent(in) :: &
       ipr, &  ! print flag
       ispn    ! number of spins
  type(parallel_gspace), intent(in) :: gs  
  type(fft_struc), intent(in) :: &
         ffts   ! information for the fast fourier transform
  !
  !     OUTPUT:
  !     ------
  !
  type(energy), intent(out) :: energs    ! only alpha term is touched
  type(complex_gspace_array), intent(out) :: &
       vion, &  ! local ionic potential
       denc, &  ! core charge
       denv, &  ! valence charge
       ddc, &   ! derivative of core charge
       dvql, &  ! derivative of local ionic potential
       vnl_ave  ! average nonlocal in real space for TF mixing
  type(double_gspace_array), intent(out) :: &
       vql, &   ! local part of pseudopotl for atom,gvec
       dnc      ! core charge for atom, gvector
  !
  !     INPUT/OUTPUT:
  !     ------------
  !
  type(pseudo_potential), intent(inout) :: pspot
  type(crystal), intent(inout) :: &
       crys    ! only valence charge zv() and ztot are altered
  !
  !     reads the fourier pseudo potentials and the
  !     atomic core and valence charge densities from
  !     files whose name depend on the chemical symbol.
  !     it then computes the same quantities on the prototype
  !     g-vectors.
  !     adapted from sverre froyen plane wave program
  !     written february 2 1990. jlm
  !
  !     Parallel and fortran 90 version 1996 Bernd Pfrommer
  !
  !     Multiple projector version 1997 Young-Gui Yoon and Bernd Pfrommer
  !
  !     Faster parallel i/o version 1997 Bernd Pfrommer.
  !     Now processor 0 only reads the file, and distributes the data
  !     via broadcast. This is much faster than each processor
  !     accessing separately.
  !
  !     added support for LW Wang version of Martins code and FHI98PP
  !     in w_linerq and support for pulay_tf mixing with calc_vnl_ave
  !     DBR May 2001
  !
  !     ------------ local variables ---------------------------------
  !
  character(len=1) :: symb(2)  
  character(len=2) :: namel, icorrt  
  character(len=3) :: irel  
  character(len=4) :: icore  
  character(len=8) :: tfile  
  character(len=10) :: iray(6), ititle(7)  
  character(len=2) :: icorr    ! type of correlation
  character(len=250) :: fnam   ! file name for open unit
  integer :: fnamlen  
  integer, external :: findvec  
  real(dp) :: vql0, fac, glmax, qj, xn, delql, xdum, q2vn, q2vp, &
       q2vm, vqj, dvqj, dcj, ddcj, dvj, fi,sum,dzv,f1,f2,f3,x
  complex(dp) sumc
  real(dp) vda(pspot%mxdlqp),denc_atom(pspot%mxdlqp),rhoq(pspot%mxdlqp),&
             vnl_atom(pspot%mxdlqp) ! just some temporary arrays 
  logical :: isp  
  integer :: i, j, it, nt, izv, nql, norb, lp1, n, idum, ind, ja, &
       is, igvec(3),nrefpts,iq
  integer natom_tot,mrb2,mrb2_matom_node,i00,iatsum,iatsum2,ia,nref,iatom,ierr
  real(dp), allocatable :: wmasktmp(:,:),xyzmaptmp(:)
  integer, allocatable :: indmtmp(:)
  real(dp) s,avenmap,t0
  real(dp), external :: gimmetime  
  !   
  !     ------------------------------------------------------------------
  !
  icorr = pw_params%icorr    ! type of correlation
  t0 = gimmetime() 
  fac = done / sqrt(crys%vcell)  
  !
  !     allocate memory for the pseudo potential structure
  !
  pspot%ntype = crys%ntype  
  allocate(pspot%vkb(pspot%mxdlqp,5,pspot%ntype));        pspot%vkb=dzero; 
  allocate(pspot%vkb_mask(pspot%mxdlqp,5,pspot%ntype));   pspot%vkb_mask=dzero;
  allocate(pspot%d2vkbdq2(pspot%mxdlqp-2,5,pspot%ntype)); pspot%d2vkbdq2=dzero
  allocate(pspot%qi(pspot%mxdlqp, pspot%ntype));          pspot%qi=dzero
  allocate(pspot%nqnl(pspot%ntype));                      pspot%nqnl=0 
  allocate(pspot%delqnl(pspot%ntype));                    pspot%delqnl=dzero
  allocate(pspot%nkb(5, pspot%ntype)) ;                   pspot%nkb = 0
  allocate(pspot%lo(5, pspot%ntype));                     pspot%lo = 0  

  allocate(pspot%nrr(pspot%ntype));                       pspot%nrr=0;
  allocate(pspot%r(2000, pspot%ntype));                   pspot%r=dzero;
  allocate(pspot%vr_pp(2000,5, pspot%ntype));             pspot%vr_pp=dzero; 
  allocate(pspot%wr_pp(2000,5, pspot%ntype));             pspot%wr_pp=dzero;
  !
  !     allocate memory for the gspace arrays
  !
  call create_gspace_array(vion, gs%length, 1, ispn)  
  call create_gspace_array(denv, gs%length, 1, ispn)  
  call create_gspace_array(vql, gs%length, pspot%ntype, 1)  
  call create_gspace_array(dnc, gs%length, pspot%ntype, 1)  
  call create_gspace_array(denc, gs%length, 1, 1)  
  call create_gspace_array(ddc, gs%length, 1, 1)  
  call create_gspace_array(dvql, gs%length, 1, 1)  
  call create_gspace_array(vnl_ave, gs%length, 1, 1)

  isp = .false.  
  if (ispn == 2) isp = .true.  

  if (ipr == 1) write(9, 100)    ! write heading
  !
  !      initialize arrays
  !
  do i = 1, gs%length  
     denc%data(i, 1, 1) = zzero
     dvql%data(i, 1, 1) = zzero
     ddc%data(i, 1, 1) = zzero
     vnl_ave%data(i, 1, 1) = zzero
  end do
  do is = 1, ispn  
     do i = 1, gs%length  
        vion%data(i, 1, is) = zzero
        denv%data(i, 1, is) = zzero
     end do
  end do
  energs%alpha = dzero  
  crys%ztot = dzero  
  !
  !  calculate total # of atoms and allocate arrays for real space NL 
  !
  natom_tot=0
  do nt = 1, crys%ntype 
    natom_tot=natom_tot+crys%natom(nt)
  end do   
  pspot%natom_tot=natom_tot
  pspot%NL_rspace=pw_params%NL_rspace
  iatom=0
  iatsum = 0
  iatsum2=0 
  nrefpts=0
  pspot%mrb2_matom_node=0
  if(pw_params%NL_rspace(1)) then
    ! this formula sometime is not right, if many atoms are located inside 1 PE
    mrb2=min(int(2*(4*pi/3*pw_params%NLPP_rcut(1)**3)/(crys%vcell/&
              (ffts%fftsize(1)*  ffts%fftsize(2) * ffts%fftsize(3) ))),&
               ffts%fftsize(1)*  ffts%fftsize(2) * ffts%fftsize(3))*2      
    mrb2_matom_node=mrb2*natom_tot/gs%nproc
    pspot%mrb2_matom_node=mrb2_matom_node

    allocate(pspot%kb_rsp(9*mrb2_matom_node)); pspot%kb_rsp=dzero 
    allocate(pspot%xyzmap(3*mrb2_matom_node)); pspot%xyzmap=dzero
    allocate(pspot%indm(mrb2_matom_node));     pspot%indm=0  
    allocate(pspot%nmap(natom_tot));           pspot%nmap=0       
    allocate(pspot%numref(natom_tot));         pspot%numref=0

    allocate(wmasktmp(9,mrb2));                wmasktmp=dzero
    allocate(xyzmaptmp(3*mrb2));               xyzmaptmp=dzero
    allocate(indmtmp(mrb2));                   indmtmp=0
   end if    
  !
  !      START LOOP OVER ATOMIC TYPES
  !
  tfile = '_POT.DAT' 
  do nt = 1, crys%ntype  

     it = 10 + nt  
     !
     !        ------ open file ------
     !
     symb(1) = crys%nameat(nt)(1:1)  
     symb(2) = crys%nameat(nt)(2:2)  
     if (symb(1) /= ' ' .and. symb(2) /= ' ') then  
        fnam = crys%nameat(nt)//tfile  
     else if (symb(1) == ' ') then  
        fnam = symb(2)//tfile  
     else  
        fnam = symb(1)//tfile  
     end if

     fnamlen = index(fnam, ' ') - 1
     !
     ! for formats in real space
     !
     if (pw_params%pp_format .eq. 2 .or.  pw_params%pp_format .eq. 3 ) then
       !
       !  set spacing for 1d k-point grid: have -2 not -1 because of shift of
       !    g=0 to (2) later on.
       !
       pspot%delqnl(nt)=sqrt(4*pw_params%emax)*1.2d0/(pspot%mxdlqp-2.d0)  
       delql= pspot%delqnl(nt)
       pspot%nqnl(nt)=pspot%mxdlqp
       nql=pspot%nqnl(nt)
       !
       !  create mask function for real space non-local if needed
       !
       pspot%ri(1)=dzero; pspot%amr(1)=done;
       call maskr(pspot%ri(2),pspot%amr(2))
       !
       !  read in potentials, atomic wavefunctions, and get KB form
       !  w_linerq gets arrays in q-space starting with (2) 
       !  as g=0 is assumed in that position for interploation. 
       !  so pspot%mxdlqp-1 points are used.
       !
       if (gs%myproc == 0) then
         open(unit = it, file = fnam(1:fnamlen) , status = 'old', &
             form = 'formatted')
         vda=dzero; denc_atom=dzero; rhoq=dzero; vnl_atom=dzero;
         call w_linerq(it,nt,crys%zv(nt),vql0,vda,denc_atom,rhoq,&
            vnl_atom,norb,pspot,pw_params)
       end if 
       !
       ! Broadcast values to other procs and set values of some data structures
       !
       call my_broadcast(crys%zv(nt), 0)
       call my_broadcast(vql0, 0)  
       call my_broadcast(pw_params%pp_iloc(nt),0)
       call my_broadcast(pspot%qi(1,nt),pspot%mxdlqp,0)

       energs%alpha = energs%alpha + vql0 * real(crys%natom(nt), dp)  
       crys%ztot = crys%ztot + crys%zv(nt) * real(crys%natom(nt), dp)     
        
       call my_broadcast(norb, 0)  
       call my_broadcast(pspot%nkb(1, nt), min(norb, 5), 0) 
       call my_broadcast(vda(1), nql, 0) 
       call my_broadcast(vnl_atom(1), nql, 0) 

       do lp1 = 1, norb    
         call my_broadcast(pspot%vkb(1, lp1, nt), pspot%nqnl(nt), 0)
         call my_broadcast(pspot%vkb_mask(1, lp1,nt), pspot%nqnl(nt), 0)
       end do
       !
       ! calculate non-local average on unfiorm grid in q-space for TF mixing
       !
       call calc_vnl_ave(vnl_ave,ffts,gs,pspot,crys,vnl_atom,nt)
       !
       !   do setup for real space NON-local
       !
       if(pw_params%NL_rspace(1)) then
         !         
         !   calculate 2nd der. for spline of vkb_mask onto 3d g space
         do j = 1, 5  
          if (pspot%nkb(j, nt) /= 0) then  
           call reg_grid_spline(pspot%vkb_mask(3, j,nt), pspot%nqnl(nt) - 2, &
                pspot%d2vkbdq2(1, j, nt))
          end if
         end do
         !
         ! for each atom getwmask gets pot*mask at each pt within rcut
         !
         do ia=1,crys%natom(nt)
           iatom=iatom+1

           call getwmask(crys%rat(1,ia,nt),pspot%nmap(iatom), indmtmp, &
             pspot%vkb_mask(1,1,nt), pspot%mxdlqp,norb, wmasktmp,xyzmaptmp,&
             mrb2,nref, crys,ffts,gs,pspot%amr,pspot%ri,pw_params%pp_iloc(nt),&
           pw_params%NLPP_rcut(nt),pspot%delqnl(nt),pspot%d2vkbdq2(1,1,nt),&
             pw_params%nbandsfft )

           pspot%numref(iatom)=nref
           if(iatsum+pspot%nmap(iatom).gt.mrb2_matom_node) then
             print*, "iatsum.gt.mbr2_matom_node, stop",&
                  iatsum,pspot%nmap(iatom),mrb2_matom_node,ia,mrb2
             call mystop 
           endif
           !
           !  pack in permanent variables with temporary ones.
           !
           i00=3*iatsum
           do i=1,pspot%nmap(iatom)
             do j=1,nref
               iatsum2=iatsum2+1
               pspot%kb_rsp(iatsum2)=wmasktmp(j,i)
             enddo
           enddo
           do i=1,pspot%nmap(iatom)
             pspot%indm(i+iatsum)=indmtmp(i)
           enddo

           do i=1,3*pspot%nmap(iatom)
             pspot%xyzmap(i+i00)=xyzmaptmp(i)
           enddo
           iatsum = iatsum + pspot%nmap(iatom)     
           nrefpts = nrefpts + pspot%nmap(iatom)*nref

         end do

         write(9,*) nref,' nref  ',pw_params%NLPP_rcut(nt),'  rcut'
         if (nt .eq. crys%ntype) then
           s = real(iatsum,dp)
           pspot%NLPP_rsp_pts=iatsum
           pspot%NLPP_rsp_nrefpts=nrefpts
           call all_sum_all(s)       
           avenmap = s/real(natom_tot,dp)
           if(gs%myproc .eq. 0) write(9,*) "ave nmap=",avenmap
         end if

       endif !NL_rspace
 
     else  ! In the case of the old-fashioned PP format...
     !
     !        ------- read heading ----
     !
       if (gs%myproc == 0) then  
          open(unit = it, file = fnam(1:fnamlen) , status = 'old', &
             form = 'unformatted')
          read(it) namel, icorrt, irel, icore, (iray(j), j = 1, 6), &
             (ititle(j), j = 1, 7)
          read(it) izv, nql, delql, vql0  
       end if
       call my_broadcast(namel, 2, 0)
       call my_broadcast(icorrt, 2, 0)
       call my_broadcast(irel, 3, 0)
       call my_broadcast(icore, 4, 0)  
       do j = 1, 6
         call my_broadcast(iray(j), 10, 0)  
       end do
       do j = 1, 7
         call my_broadcast(ititle(j), 10, 0)  
       end do
       call my_broadcast(izv, 0)
       call my_broadcast(nql, 0)  
       call my_broadcast(delql, 0)
       call my_broadcast(vql0, 0)  
       crys%zv(nt) = real(izv, dp)  
       pspot%delqnl(nt) = delql  
       if (ipr >= 1) write(9, 101) namel, icorrt, irel, icore, &
          (iray(j), j = 1, 6), (ititle(j), j = 1, 7), nql, delql
       if (nql+2 > pspot%mxdlqp) call fatal(130, xdum, pspot%mxdlqp, 'pseukb')
       if (namel /= crys%nameat(nt)) then  
         if (.not. ((namel(1:1) == ' ') .and. (crys%nameat(nt)(2:2) == ' ') &
             .and. (namel(2:2) == crys%nameat(nt)(1:1)) .or. &
             ((namel(2:2) == ' ') .and. (crys%nameat(nt)(1:1) == ' ') .and. &
             (namel(1:1) == crys%nameat(nt)(2:2))))) call warn(132, xdum, &
             idum, 'read_pseudo')
       end if
       if (icorrt /= icorr) then  
         write(9, *) 'icorrt=', icorrt  
         write(9, *) 'icorr=', icorr  
         icorrt = icorr  
         call warn(133, xdum, idum, 'read_pseudo')  
       end if
       !
       !        sum up vql0 for the alpha energy
       !        and zv for the total valence charge
       !
       energs%alpha = energs%alpha + vql0 * real(crys%natom(nt), dp)  
       crys%ztot = crys%ztot + crys%zv(nt) * real(crys%natom(nt), dp)  
       !
       !        read pseudopotentials
       !        note that q=0 is shifted to vda(2) and vkb(2,l,k)
       !
       nql = nql + 2  

       pspot%nqnl(nt) = nql  
       !
       !     --- read normalization and local potential in q-space ------------
       !
       pspot%nkb(:, nt) = 0
       !        maximum allowed norb is increased to 5
       !        norb is actually nprojec in newkb
       !
       if (gs%myproc == 0) read(it) norb, (pspot%nkb(j, nt), &
          j = 1, min(norb, 5))


       call my_broadcast(norb, 0)  
       call my_broadcast(pspot%nkb(1, nt), min(norb, 5), 0)  
       call myflush(9)
       if (norb > 5) then  
         write(9, *) 'TOO MANY NONLOCAL PROJECTORS FOUND IN read_pseudo'  
         write(9, *) 'NUMBER OF PROJECTORS=', norb  
         call mystop  
       end if  
       !
       !     - local part
       !
       if (gs%myproc == 0) read(it) (vda(j), j = 3, nql)  
       call my_broadcast(vda(1), nql, 0)  
       !
       !     --- read nonlocal part of the Pseudopotential -----------
       do lp1 = 1, norb  
         if (gs%myproc == 0) read(it) (pspot%vkb(j, lp1, nt), j = 2, &
             pspot%nqnl(nt))
         call my_broadcast(pspot%vkb(1, lp1, nt), pspot%nqnl(nt), 0)
       end do
 
     end if   ! Now done for all PP_formats
     !
     !        check argument limits and compare with gmax
     !
     glmax = delql * real(nql - 3, dp)  
     if (glmax < gs%gmax) then  
        write(9, *) 'LOCAL POTENTIAL IS SET UP TO GLMAX=', glmax  
        write(9, *) 'BUT GMAX IS LARGER:', gs%gmax  
     end if
     !
     !        COMPUTE IONIC LOCAL POTENTIAL
     !
     do j = 1, gs%length  
        vql%data(j, nt, 1) = dzero  
        !
        !          QUADRATIC INTERPOLATION OF Q**2 * POTENTIAL
        !
        qj = sqrt(gs%ekin(j))  
        if (qj > dzero) then  
           xn = qj / delql + dtwo  
           n = xn + dhalf  
           if (n < nql) then  
              if (n <= 3) n = 4  
              xn = xn - real(n, dp)  
              q2vn = real((n - 2) * (n - 2), dp) * vda(n)  
              q2vp = real((n - 1) * (n - 1), dp) * vda(n + 1)  
              q2vm = real((n - 3) * (n - 3), dp) * vda(n - 1)  

              vqj = q2vn * (done - xn) * (done + xn) + dhalf * (q2vp * &
                   (done + xn) - q2vm * (done - xn)) * xn
              vqj = delql * delql * vqj / (crys%vcell * gs%ekin(j))  
              dvqj = dmtwo * q2vn * xn + q2vp * (dhalf + xn) - q2vm * &
                   (dhalf - xn)
              dvqj = (delql * dvqj / (dtwo * crys%vcell * qj) - vqj) / &
                   gs%ekin(j)
              do is = 1, ispn  
                 vion%data(j, 1, is) = vion%data(j, 1, is) + vqj * &
                      gs%struc(nt, j)
              end do

              !                required for force/stress calculations:
              dvql%data(j, 1, 1) = dvql%data(j, 1, 1) + dvqj * &
                   gs%struc(nt, j)
              ! add all up
              ! keep separately
              vql%data(j, nt, 1) = vqj  
           end if
        end if
     end do
     !     ---------------- deal with core and valence charge ---------------
     !
     !       read core charge
     ! 
     if (gs%myproc == 0) then
       if (pw_params%pp_format .eq. 2 .or.  pw_params%pp_format .eq. 3 ) then 
         call mdcopy(nql,denc_atom(1),1,vda(1),1)
       else    
         read(it) (vda(j), j = 3, nql)  
       end if
     end if

     call my_broadcast(vda (1), nql, 0)

     vda(1) = vda(3)  
     vda(2) = vda(3) - (vda(4) - vda(3) ) * dthird
     !
     !       compute core charge
     !
     do j = 1, gs%length  
        dnc%data(j, nt, 1) = dzero  
        !
        !          interpolate vda
        !
        qj = sqrt(gs%ekin(j))  
        xn = qj / delql + dtwo  
        n = xn + dhalf  
        if (n < nql) then  
           xn = xn - real(n, dp)  
           dcj = vda(n) * (done + xn) * (done - xn) + dhalf * (vda(n + 1) * &
                (done + xn) - vda(n - 1) * (done - xn)) * xn
           dnc%data(j, nt, 1) = dcj  
           denc%data(j, 1, 1) = denc%data(j, 1, 1) + dcj * gs%struc(nt, j)
           if (gs%ekin(j) < 1.0d-16) then  
              ddc%data(j, 1, 1) = ddc%data(j, 1, 1) + gs%struc(nt, j) * &
                   (vda(3) - vda(2)) / (delql * delql)
           else  
              ddcj = dmtwo * vda(n) * xn + vda(n + 1) * (dhalf + xn) - &
                   vda(n - 1) * (dhalf - xn)
              ddcj = ddcj / (dtwo * delql * qj)  
              ddc%data(j, 1, 1) = ddc%data(j, 1, 1) + ddcj * gs%struc(nt, j)
           end if
        end if

     end do 
     !
     !       read valence charge
     !
     if (gs%myproc == 0) then
       if ( pw_params%pp_format .eq. 2 .or.  pw_params%pp_format .eq. 3 )  then
          call mdcopy(nql,rhoq(1),1,vda(1),1)
       else    
          read(it) (vda(j), j = 3, nql)  
       end if
     end if

     call my_broadcast(vda(1), nql, 0)  
     vda(1) = vda(3)  
     vda(2) = crys%zv(nt)  
     !
     !       compute valence charge
     !
     do j = 1, gs%length  
        !
        !          interpolate vda
        !
        xn = sqrt(gs%ekin(j)) / delql + dtwo  
        n = xn + dhalf  

        if (n < nql) then  
           xn = xn - real(n, dp)  
           dvj = vda(n) * (done + xn) * (done - xn) + dhalf * (vda(n + 1) * &
                (done + xn) - vda(n - 1) * (done - xn)) * xn
           !
           !             sum up the charge density
           !
           if (isp) then  
              denv%data(j, 1, 1) = denv%data(j, 1, 1) + dvj * &     ! spin up
                   gs%struc(nt, j) * dhalf * (done + crys%szet(nt))
              denv%data(j, 1, 2) = denv%data(j, 1, 2) + dvj * &     ! spin down
                   gs%struc(nt, j) * dhalf * (done - crys%szet(nt))
           else  
              denv%data(j, 1, 1) = denv%data(j, 1, 1) + dvj * &
                   gs%struc(nt, j)
           end if
        end if
     end do
     !
     !     pspot%lo is read for identification of l quantum number
     !     To ensure compatiblity with old pspots, check if lo info
     !     is available, and initialize correctly if not.
     !
     pspot%lo(:, nt) = 0  
     pspot%lo(1, nt) = 0
     pspot%lo(2, nt) = 1
     pspot%lo(3, nt) = 2
     !     read information about the angular momenta of the nonlocal
     !     pspots
     if (gs%myproc == 0 .and. pw_params%pp_format == 1) then  
        read(it, end = 912) (pspot%lo(j, nt), j = 1, min(norb, 5))
        goto 913  
912     continue  
        write(9, *) ' ... consider upgrading to new kbconv ...'  
913     continue  
     end if
     call my_broadcast(pspot%lo(1, nt), min(norb, 5), 0)  
     !
     !       do some strange extrapolation to 0 i.e. extend vkb to small
     !       negative ql according to parity (-1)**l so that interpolation
     !       works around 0.
     !
     do lp1 = 1, norb  
        pspot%vkb(1, lp1, nt) = -1**pspot%lo(lp1, nt) * pspot%vkb(3, lp1, nt)
     end do
     do lp1 = 1, norb  
        do j = 1, pspot%nqnl(nt)  
           pspot%vkb(j, lp1, nt) = fac * pspot%vkb(j, lp1, nt)  
        end do
     end do
     if (gs%myproc == 0) close(unit = it)  
  end do     ! end of loop over atomic types
  !
  !     compute the alpha energy. this is the contribution
  !     to the total energy coming from g=0.
  !

!------------------------------------------------------------------
! P. Zhang, March, 2004
! For charged system, the alpha term should be computed
! using the total electronic  charge density, not the ionic charge.
!------------------------------------------------------------------
! 
!  energs%alpha = crys%ztot * energs%alpha / crys%vcell  
  energs%alpha = (crys%ztot+crys%net_charge) * energs%alpha / crys%vcell  
  if (pw_params%occupy .eq.  4) then
  energs%alpha = ( crys%nalpha+ crys%nbeta) * energs%alpha / crys%vcell  
  end if

  !
  !     compute the number of components the nonlocal part has
  !
  !     For each atom, we have for a different potential component
  !
  pspot%nanl = 0  
  do i = 1, crys%ntype  
     do j = 1, 5  
        if (pspot%nkb(j, i) /= 0) pspot%nanl = pspot%nanl + &
             (2 * pspot%lo(j, i) + 1) * crys%natom(i)
     end do
  end do
  if (ipr >= 1) write(9, 200) pspot%nanl  
  call myflush(9)  
  !
  !     compute second derivatives for spline interpolation
  !
  do i = 1, crys%ntype  
     do j = 1, 5  
        if (pspot%nkb(j, i) /= 0) then  
           call reg_grid_spline(pspot%vkb(3, j, i), pspot%nqnl(i) - 2, &
                pspot%d2vkbdq2(1, j, i))
        end if
     end do
  end do
  !
  !     put in an external electric potential if desired
  !
  !     the shape is  V_e(r) = eampli * sin(Gr)
  !
  !     but then it is symmetrized with the symmetry operations
  !     of the crystal.
  !
  if (pw_params%epot%gvec(4) > 0) then  
     igvec = pw_params%epot%gvec(1:3)  
     i = findvec(igvec, gs)  

     if (i > 0) vion%data(i, 1, :) = vion%data(i, 1, :) - dhalf * &
          cmplx(dzero, pw_params%epot%h, dp)
     igvec = -pw_params%epot%gvec(1:3)  
     i = findvec(igvec, gs)  
     if (i > 0) vion%data(i, 1, :) = vion%data(i, 1, :) + dhalf * &
          cmplx(dzero, pw_params%epot%h, dp)
     do i = 1, crys%nspin  
        call symmetrize_scalar_global(vion%data(1, 1, is), gs)  
     end do
  end if

  if(pw_params%NL_rspace(1)) then 
    deallocate(wmasktmp)
    deallocate(xyzmaptmp)
    deallocate(indmtmp)
  end if

  write(9,*) gimmetime() - t0,' PP setup time' 

  return  

200 format(/' Nonlocal potential: ',/' ------------------', &
       &     /' Number of projectors: ',i5)
100 FORMAT(//,' POTENTIALS :',/,1X,12('-'))  

101 FORMAT(/,1X,A2,2X,A2,2X,A3,2X,A4,/,1X,6A10,/,1X,7A10, &
       &     /,' NQL=',I4,' DELQL=',F5.3)

end subroutine read_pseudo













