!-*-F90-*-
!     @process extchk
!
subroutine scfloop(t0,imode, pot_gspace, k_gspace, energs, &
     pw_params, crys, syms, ffts, pspot, kpoints, bands, wavefn, lrwfn, &
     vion, vin, vout, den, denc, chdr,vnl_ave)
  !
  !     1996 Bernd Pfrommer
  !
  use esdf
  !
  use all_to_all_module       ! added by hjchoi

!-----------------
! LSDA+U. P. Zhang
!-----------------
  use lsdau_shared         

  include 'use.h'
  implicit none               ! never remove this line.
  include 'interface.h'  
  include 'flibcalls.ph'  
  include 'all_to_all.h'      ! added by hjchoi
  !
  !     INPUT:
  !     -----
  !
  integer, intent(in) :: imode            ! operation mode
  type(symmetry), intent(in) :: syms      ! symmetry operations
  type(crystal), intent(in) :: crys       ! crystal structure
  type(fft_struc), intent(in) :: ffts     ! information for the FFT
  type(pw_parameter), intent(inout) :: &
       pw_params                          ! plane wave parameters
  type(pseudo_potential), intent(in) :: &
       pspot                              ! all info about the pseudopotentials
  type(kpoint), intent(in) :: kpoints     ! BZ integration kpoints
  type(complex_gspace_array), intent(in) :: &
       denc, &                            ! core charge density
       vnl_ave  ! average nonlocal in real space for TF mixing
  logical, intent(in) :: lrwfn            ! true if wavefn is to be reused,
                                          ! false otherwise
  type(parallel_gspace), intent(in) :: &
       pot_gspace, k_gspace(*)
  real(dp) t0  
  !
  !     INPUT/OUTPUT:
  !     -------------
  !
  type(energy), intent(inout) :: energs   ! energies and other information
  type(band), intent(inout) :: bands      ! eigenvalues, occupation numbers etc
  type(complex_gspace_array), intent(inout) :: &
       vion, &                            ! screened ionic potential
       den, &                             ! total charge density
       vin, &                             ! input potential
       vout                               ! output potential
  type(complex_gspace_array), intent(inout) :: &
       wavefn                             ! wave functions for all kpoints,
                                          !  bands
  !
  !     WORK:
  !     ----
  !
  type(complex_rspace_array) :: &
       chdr                               ! work array:
                                          ! charge density in realspace
  !
  !
  !     DESCRIPTION:
  !     -----------
  !
  !     Either computes self-consistent field (if imode =1) or
  !     does just a single iteration (if imode =0, like for a
  !     bandstructure calculation). If imode=2, reuse the memory
  !     of the wavefunction for different k-point. Use this option
  !     when the wavefunction are not needed, it save memory in big
  !     bandstructure calculation.
  !
  !           imode = 0    single iteration, keep all wavefunctions
  !           imode = 1    self-consistent field
  !           imode = 2    single iteration, no wavefunction
  !
  !     above Peter
  !
  !     In case lrwfn (LogicalReuseWaveFuNction) is true, the
  !     subroutine assumes that "bands" and "wavefn" contain a good
  !     starting guess.
  !
  !     3/01 added Pulay-Kerker mixing and Pulay-Thomas-Fermi mixing
  !
  !          also added grassmann and grassmann minimization metal routines
  !
  !     --------------- local variables -------------------------
  !
  type(parallel_gspace) :: &
       sub_gspace, &                     ! for the starting guess
       mix_gspace                        ! gspace for mixing
  complex(dp), pointer :: vloc(:)        ! local potential in realspace
  type(hamiltonian) :: ham               ! hamiltonian at a given k-point
  type(dense_hamiltonian) :: densham     ! subspace hamiltonian
  type(energy) :: oldenergs              ! energies of the previous iteration
  complex(dp), allocatable ::  &
       ajac(:)                            ! inverse of hessian: mixing matrix
  complex(dp),allocatable ::  R0(:),w_in0(:),&
        dw(:,:,:),dR(:,:,:)
  real(dp) ::  &
       errscfconv                        ! largest error in the potential
  integer :: i, j, n,&
       maxiter, &                        ! maximum number of scf iterations
       maxitdiag, &                      ! maximum number of diag iterations
       nsp, nvecs, &
       ierr, &                           ! error flag
       is,is2, &                         ! spin counter
       irk,irk2, &                       ! kpoint counter
       iter,&                            ! iteration number
       num_HPsi, &                          ! op count of H*Psi  
       mypr,mb, nrk
       

 logical :: istartguess, &              ! whether to compute the starting guess
       iscfconv, iscfconv2, &           ! indicates convergence of scf loop
       idiaconv                         ! indicates convergence of diag.
! davegp
  logical :: idiaconv_total
! davegp
  character(len=11) :: wfnname, jacname,&
       hwfname                           ! add by hjchoi
  real(dp), external :: gimmetime  
!  real(dp)  AA(itmaxpul,itmaxpul)
  real(dp) denom,oldenergy,t_temp,occ_tem
  integer::nb,nspin,nnrk,nkk,ib,ipr

  integer length( kpoints%nrk)
  real(dp) time(num_time)!1 fft 2 nl 3occ_spa 4othercg 5otherscf 6 olap_cg 7 rotations
  !  5 init_diag 8 flevel 9 charge 10 velect 11 etotal 12 io 13 mixer 14 setup
  !
!     Following variables added by David Roundy
!
      complex(dp), allocatable :: & ! of size nanl x m
           ews(:,:)
      integer m,len, k          ! Dimensions of my matrices.
      complex(dp), allocatable :: & ! of size m x m
           ham_overlap(:) &                       
           , f(:,:) &                             
           ,wfn_overlap(:) &             ! added by hjchoi
           , f0(:), hf0(:)               ! added by hjchoi    
      integer, external :: findvec       ! added by hjchoi
      real(dp):: save_shift, HG0         ! added by hjchoi
      integer :: idirect,imlen,jgb(3)    ! added by hjchoi
      real(dp) time2(num_time)           ! added by hjchoi
!     End of David's Variables                    
                                                  
  integer itmaxpul  

  !     --------------------------------------------------------------
  !
  iscfconv2=.false.
  iscfconv=.false.
  itmaxpul=20
  mypr=pot_gspace%myproc
  nproc=pot_gspace%nproc

  t_temp=gimmetime()
  time=dzero
  time2=dzero ! hjchoi
  num_HPsi=0
  oldenergs = dzero
  oldenergy=dzero
  ham%ekinmod = pw_params%ekinmod  
 
!-----------------
! LSDA+U. P. Zhang
!-----------------

  Etot_corr=0.d0
  ELSDA_corr=0.d0
  Etot_corr_old=0.d0
  Eband_diff=0.d0

     mb=MIN(bands%min(1),bands%min(2))
     nrk=kpoints%nrk
     if((pw_params%lsdau.eq.1).and.(init_flag.ne.999)) then
         call lsdau_init(mypr,nproc,mb,nrk,imode,crys,syms,kpoints)
         call lsdau_proj_init(k_gspace,kpoints,crys)
     end if

  if(pw_params%mix_method == 'broyden') then  
    !
    !     generate a gspace which is local to all processors
    !     for the mixing. find inversion information, which
    !     is used in the mixing routine.
    !
    mix_gspace%gmax = sqrt(pw_params%emaxmix)  
    mix_gspace%rk = (/ dzero, dzero, dzero /)     ! no shift
    mix_gspace%nproc = 1                          ! the same on all processors
    mix_gspace%name = 'mixing'  
    mix_gspace%myproc = 0
    mix_gspace%fftsize = (/ 0, 0, 0, 0, 0 /)  
    mix_gspace%istar = .true.
    mix_gspace%igvec = .false.
    mix_gspace%iekin = .true.
    mix_gspace%iinverse = .true.  
    mix_gspace%istruc = .false.  
    mix_gspace%ipackinfo = .false. 

!    call generate_gspace(pw_params%output(1), mix_gspace, crys, syms, 1)
    call generate_gspace(pw_params%output(1), mix_gspace, crys, syms)
    allocate(ajac(crys%nspin * crys%nspin * mix_gspace%length * &
       mix_gspace%length))                      ! allocate mix matrix

  else if (pw_params%mix_method == 'pulay_kerker' .or. &
          pw_params%mix_method == 'pulay_tf'  ) then
    allocate(R0(crys%nspin * pot_gspace%length ))
    allocate(w_in0(crys%nspin * pot_gspace%length))
    allocate(dw( pot_gspace%length,crys%nspin,itmaxpul))
    allocate(dR( pot_gspace%length,crys%nspin,itmaxpul))
  end if

  allocate(vloc(ffts%r_size))
 
  if (pw_params%maxitscf > 1) write(9, 100)  
  call myflush(9)

  if (iand(pw_params%input, 1) == 1) call read_band('BAND', bands)
  if (bands%nrk <= 0) write(9, 155)  

  if (iand(pw_params%input, 2) == 2) call lukman(crys, pot_gspace, &
       vion%data(1, 1, 1))
  
  if (imode == 1) then  
     maxiter = pw_params%maxitscf  
  else  
     maxiter = 1  
  end if
  
  if (imode == 2) then  
     maxitdiag = pw_params%maxitdiag_band 
  else  
     maxitdiag = pw_params%maxitdiag
  end if


!------------------------------------
!P.Zhang, fixed occupation from input

! read in occupation number if asked
  if(pw_params%occupy.eq.5) then
     write(9,*) "fixed occupation from input"
     open(109,file="occ.dat",err=234)
     read(109,*) nb,nspin,nnrk
     if(nb.ne.MIN(bands%min(1),bands%min(2))) call mystop
     if(nspin.ne.bands%nspin) call mystop

     bands%occup=0.d0

     nnrk=kpoints%grid(1)*kpoints%grid(2)*kpoints%grid(3)

     do is = 1,nspin
        nkk=0
        do irk=1,nnrk
           if(kpoints%kmap(irk).lt.0) nkk=nkk+1
           do ib = 1, nb
           read(109,*) occ_tem
           if(kpoints%kmap(irk).lt.0) then
              bands%occup(ib,nkk,is)=occ_tem*kpoints%w(nkk)*2.d0
           end if
           end do
        end do
     end do

     close(109)
 234 continue
     write(*,*) "cannot file occ.dat"
     call mystop

   end if

!
! PZ
  if(pw_params%iexact_diag.eq.1) maxiter=1

  do iter = 1, maxiter 
     if(iand(pw_params%output(4),2).eq.2) then
        open(unit=37,file='elpho')
        write(37,*) 'electron phonon matrix element file'
     endif 
     write(9,*)  gimmetime()-t0,' SCFLOOP TIME  BEFORE ITERATION',iter 
     !
     !        compute the wave functions
     !

! davegp
     idiaconv_total = .true.
! davegp

     do is = 1, crys%nspin  
        call setup_local_potential(1, ffts, pot_gspace, vion%data(1, 1, is), &
             vloc(1))
        call myflush(9)  
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(14),t_temp)
        do irk = 1, kpoints%nrk  
           call create_hamiltonian(ham, k_gspace(irk), pspot, vloc, ffts)  
           length(irk)= ham%gspace%length 
           call setup_nonlocal_potential(ham, pspot, crys,iter)  
           if (iand(pw_params%output(1),8)==8) call get_timing(time(14),t_temp)
           if (iter == 1 .and. (.not. lrwfn)) then  ! gen start guess or read
              istartguess = .true. 
              if (iand(pw_params%input, 1) == 1) then      ! read from file
                 istartguess = .false.  
                 call read_band('BAND', bands)  
                 if (bands%nrk <= 0) then                   ! error
                    write(9, 155)  
                    bands%nrk = kpoints%nrk
                    istartguess = .true.
                    goto 20  
                 end if
                 write(wfnname, 120) irk, is  
                 wavefn%data(:, irk, is) = zzero
                 nsp = 1
                 nvecs = bands%nband(irk, is)  
                 call readgsdat(1, ierr, k_gspace(irk), &
                      wavefn%data(1, irk, is), bands%nband(irk, is), nvecs, &
                      nsp, wfnname, 10)
                 if (ierr /= 0) then  
                    write(9, 150) wfnname
                    istartguess = .true.
                    goto 20  
                 end if
                 write(9, 140) nvecs, wfnname  
                 bands%nband(irk, is) = nvecs  
                 if (iand(pw_params%output(1), 8) == 8) &
                             call get_timing(time(12),t_temp)
              end if
20            if (istartguess) then  
                 sub_gspace%gmax = sqrt(pw_params%emaxsub)  
                 sub_gspace%rk = ham%gspace%rk ! inherit shift from big brother
                 sub_gspace%nproc = 1           ! scatter across all processors
                 sub_gspace%name = 'sub'  
                 sub_gspace%myproc = 0
                 sub_gspace%fftsize = (/ 0, 0, 0, 0, 0 /)
                 sub_gspace%istar = .false.
                 sub_gspace%igvec = .true.  
                 sub_gspace%iekin = .true.
                 sub_gspace%iinverse = .false.  
                 sub_gspace%istruc = .false.  
                 sub_gspace%ipackinfo = .false.  
 
                 call generate_gspace(pw_params%output(1), sub_gspace, crys, &
                      syms, 1)
                 if (iter.eq.1 .and. is.eq.1 .and. irk .eq. 1) &
                    call print_mem_use(pot_gspace,ffts,crys,&
                    k_gspace,pspot,pw_params,sub_gspace,ham,kpoints,bands)
                 call create_dense_hamiltonian(densham, sub_gspace)   
                 ! from here Peter
                 if(imode /= 2) then
                    call start_vectors(0, pw_params, ham, densham, &
                         bands%energy(1, irk, is), bands%max, &
                         wavefn%data(1, irk, is), pot_gspace, sub_gspace, &
                         vion%data(1, 1, is), crys)
 
                 elseif(imode == 2) then
                    call start_vectors(0, pw_params, ham, densham, &
                         bands%energy(1, irk, is), bands%max, &
                         wavefn%data(1, 1, is), pot_gspace, sub_gspace, &
                         vion%data(1, 1, is), crys)
                 end if
                 
              if (bands%nband(irk, is) > bands%min(is)) bands%nband(irk, is)= &
                      bands%nband(irk, is) - 1
                 call destroy(densham)  
                 call destroy(sub_gspace) 
               if (pspot%NL_rspace(1)) &
                      call destroy_complex_gspace_array(ham%vnloc)  
              end if
              if (iand(pw_params%output(1),8)==8) &
                                 call get_timing(time(5),t_temp)
           end if
           

           if (lrwfn .and. iter == 1) write(9, 160)  
           if (pw_params%bandgap > dzero) then        ! optimize for insulator
              bands%nband(irk, is) = bands%min(is)  
              bands%energy(bands%min(is) + 1, irk, is) = &
                   bands%energy(bands%min(is), irk, is) + pw_params%bandgap
           end if



           if(iand(pw_params%output(4),2).eq.2) then
              !  
              !    New matrix element bit.
              !
              len=ham%gspace%length
              m=bands%nband(irk, is)
              allocate(f(len,m));
              allocate(ews(ham%vnloc%nvecs,m))
              ews=0
                 
!---------------- find  H(G=0) ---  begin ---------
              allocate(f0(len),hf0(len))
              f0 = 0
              hf0 = 0
              jgb(1) = 0
              jgb(2) = 0
              jgb(3) = 0
              i = findvec (jgb,  k_gspace(irk))     ! find the index of the gvector
              if (i > 0) then                       ! gvector is here, get data
                 f0(i) = zone
              endif
              save_shift = ham%shift
              ham%shift = 0.0
              call apply_ham_wavefn(time2,pw_params%output(1), hf0, ham, f0,  &
                                    1, ews(1, 1), pw_params%nbandsfft)
              ham%shift = save_shift
              if (i > 0) then
                 HG0 = real(hf0(i))
              else
                 HG0 = dzero
              endif         
              deallocate(f0,hf0)
              call all_sum_all_double1(HG0)
!             call all_sum_all(HG0)
              write(9,*)' k = ',irk,' H(G=0) = ',HG0
              ews = 0
!---------------- find  H(G=0) ---   end  ---------
              allocate(ham_overlap(m*m));  ham_overlap=0;
              allocate(wfn_overlap(m*m));  wfn_overlap=0;  ! added by Hyoung Joon Choi
                 
              ! f = H * wfn
              save_shift = ham%shift  ! added by Hyoung Joon Choi
              ham%shift = 0.0         ! added by Hyoung Joon Choi
              call apply_ham_wavefn(time2,pw_params%output(1), f, ham, &
                   wavefn%data(1, irk, is), &
                   m, ews(1, 1), pw_params%nbandsfft)
              ham%shift = save_shift  ! added by Hyoung Joon Choi
              ! eps = wfn^T*f
              call overlap_matrix(ham_overlap(1), wavefn%data(1, irk, is), &
                   f, m, len)
!------------------------------------------------
!      added by hjchoi
!                
!             determine factional translation direction
!                
              if(kpoints%grid(1).le.kpoints%grid(2).and. &
                 kpoints%grid(1).le.kpoints%grid(3)) then
                 idirect = 1
              elseif(kpoints%grid(2).le.kpoints%grid(1).and. &
                 kpoints%grid(2).le.kpoints%grid(3)) then
                 idirect = 2  
              else
                 idirect = 3
              endif      
              do i = 1,m
                 imlen = (i-1)*len
                 do j = 1,len
                    f(j,i) = wavefn%data(j+imlen,irk,is)
                    k=k_gspace(irk)%gvec(idirect,j)
                    if(k.ne.(k/2)*2)f(j,i)=-f(j,i)
                 enddo
              enddo
              call overlap_matrix(wfn_overlap(1), wavefn%data(1, irk, is), &
                   f, m, len)
!------------------------------------------------
36            format('k = ',i5,' is = ',i5,e28.20,'          ')
37            format(2i4,4e18.10)
              write(37,36)irk,is,HG0
              do i=1,m
                 do j=1,m
                    k=(i-1)*m + j
                    write(37,37) i, j,real(ham_overlap(k)), &
                         aimag(ham_overlap(k)), &
                         real(wfn_overlap(k)),&
                         aimag(wfn_overlap(k))
                 enddo
              enddo
              deallocate(wfn_overlap)
              deallocate(ham_overlap)
              deallocate(f)
              deallocate(ews)
              !  
              !    End of new matrix elemet bit.
              !  
           endif 
               
               
           if(iand(pw_params%output(4),1).eq.1) then  ! added by H.J. Choi
              !  
              !    Save H * wfn
              !  
              len=ham%gspace%length
              m=bands%nband(irk, is)
              write(9,*)' save Ham * wfn for k = ',irk,', # of bands = ',m
              allocate(f(len,m));
              allocate(ews(ham%vnloc%nvecs,m))
              ews=0
              ! f = H * wfn
              save_shift = ham%shift
              ham%shift = 0.0
              call apply_ham_wavefn(time2,pw_params%output(1), f, ham, &
                   wavefn%data(1, irk, is), &
                   m, ews(1, 1), pw_params%nbandsfft)
              deallocate(ews)
              ham%shift = save_shift
              write(hwfname,127)irk,is
127           format('HWF',i5.5,'.',i1)
                 
!             call writegsdat(0, k_gspace(irk), wavefn%data(1, irk, is), &
!                  bands%nband(irk, is), bands%nband(irk, is), 1, wfnname, 10)
                 
              call writegsdat(0, k_gspace(irk), f, m, m, 1, hwfname, 10)
              deallocate(f)
              !  
              !    End of save H * wfn
              !  
           endif 

 
! PZ
! if exact diagonalization, perform orthogonalization
           if(pw_params%iexact_diag.eq.1) then 
             call orthogonalization_PZ(pw_params, bands%energy(1, irk, is), bands%max, &
                         wavefn%data(1, irk, is), k_gspace(irk), crys)


           else

!          
! mode 0 and mode 1 are called with different maxitdiag
! PZ
           if(imode == 0) then
             if(pw_params%diag_method == 'Mauri') then
                write(9,*) "diagonalizing kpoint",irk
                call diagonalize_hamiltonian(time,iter,num_HPsi,&
                   pw_params, bands%max, pw_params%maxitdiag_band,   &
                   idiaconv, bands%min(is), bands%nband(irk, is), &
                     pw_params%epsdiag, &
                      bands%energy(1, irk, is), ham, wavefn%data(1, irk, is))


             else if(pw_params%diag_method == 'Grassmann')then
               write(9,*) "diagonalizing kpoint",irk
               call diagonalize_hamiltonian_gcg(time,iter,num_HPsi,pw_params, &
                      bands%max, pw_params%maxitdiag_band,  &
                      idiaconv, bands%min(is), bands%nband(irk, is), &
                      pw_params%epsdiag, &
                      bands%energy(1, irk, is), ham, wavefn%data(1, irk, is))

             else if(pw_params%diag_method == 'Grassmann_metal')then
               write(9,*) "diagonalizing kpoint",irk
               call diag_ham_metal_gcg(pw_params%ilog,time,iter,num_HPsi,&
                 pw_params, pw_params%maxitdiag_band,idiaconv, bands%min(is), &
                     bands%nband(irk,is),pw_params%epsdiag, &
                      bands%energy(1, irk, is), ham, wavefn%data(1, irk, is),&
                      bands%occup(1,irk,is),abs(oldenergy-energs%total),lrwfn)
!              end if
! davegp

              else
                write(9,*) pw_params%diag_method
                call mystop( 'incorrect diag_method' )
              endif

          
           elseif(imode == 1) then
             if(pw_params%diag_method == 'Mauri') then
                write(9,*) "diagonalizing kpoint",irk
                call diagonalize_hamiltonian(time,iter,num_HPsi,&
                   pw_params, bands%max, pw_params%maxitdiag,   &
                   idiaconv, bands%min(is), bands%nband(irk, is), &
                     pw_params%epsdiag, &
                      bands%energy(1, irk, is), ham, wavefn%data(1, irk, is))


             else if(pw_params%diag_method == 'Grassmann')then 
               write(9,*) "diagonalizing kpoint",irk
               call diagonalize_hamiltonian_gcg(time,iter,num_HPsi,pw_params, &
                      bands%max, pw_params%maxitdiag,  &
                      idiaconv, bands%min(is), bands%nband(irk, is), &
                      pw_params%epsdiag, &
                      bands%energy(1, irk, is), ham, wavefn%data(1, irk, is))

             else if(pw_params%diag_method == 'Grassmann_metal')then        
! davegp
!              if (iter .eq. 1 ) then 
!               write(9,*) "diagonalizing kpoint",irk
!               call diagonalize_hamiltonian_gcg(time,iter,num_HPsi,pw_params, &
!                      bands%max, pw_params%maxitdiag,  &
!                      idiaconv, bands%min(is), bands%nband(irk, is), &
!                      pw_params%epsdiag, &
!                      bands%energy(1, irk, is), ham, wavefn%data(1, irk, is))
!              else
!               
               write(9,*) "diagonalizing kpoint",irk
               call diag_ham_metal_gcg(pw_params%ilog,time,iter,num_HPsi,&
                 pw_params, pw_params%maxitdiag,idiaconv, bands%min(is), &
                     bands%nband(irk,is),pw_params%epsdiag, &
                      bands%energy(1, irk, is), ham, wavefn%data(1, irk, is),&
                      bands%occup(1,irk,is),abs(oldenergy-energs%total),lrwfn)
!              end if
! davegp
  
              else
                write(9,*) pw_params%diag_method
                call mystop( 'incorrect diag_method' )
              endif


    if (iand(pw_params%output(1), 8) == 8)  t_temp=gimmetime()
              call kinet(ham, wavefn%data(1, irk, is), bands%nband(irk, is), &
                   bands%ekn(1, irk, is))

           elseif(imode == 2) then

             if(pw_params%diag_method == 'Mauri') then
                write(9,*) "diagonalizing kpoint",irk
                call diagonalize_hamiltonian(time,iter,num_HPsi,&
                   pw_params, bands%max, pw_params%maxitdiag_band,   &
                   idiaconv, bands%min(is), bands%nband(irk, is), &
                     pw_params%epsdiag, &
                      bands%energy(1, irk, is), ham, wavefn%data(1, 1, is))

             else if(pw_params%diag_method == 'Grassmann' .or. &
              pw_params%diag_method == 'Grassmann_metal' )then 
                  write(9,*) "diagonalizing kpoint",irk
                  call diagonalize_hamiltonian_gcg(time,iter,num_HPsi,pw_params, &
                      bands%max, pw_params%maxitdiag_band,  &
                      idiaconv, bands%min(is), bands%nband(irk, is), &
                      pw_params%epsdiag, &
                      bands%energy(1, irk, is), ham, wavefn%data(1, 1, is))
           
              else
                write(9,*) pw_params%diag_method
                call mystop( 'incorrect diag_method' )
              endif

              call kinet(ham, wavefn%data(1, 1, is), bands%nband(irk, is), &
                   bands%ekn(1, irk, is))

              if (crys%nspin == 2 .and. is == 1) write(9, '(1x,a)') 'Spin up'  
              if (crys%nspin == 2 .and. is == 2) write(9, '(1x,a)') 'Spin down'  
              write(9, 500) irk  
              do n = 1, bands%nband(irk, is), 7  
                 write(9, 510) (bands%energy(j, irk, is) + energs%vxc0(is), &
                      j = n, min(bands%nband(irk, is), n + 6))
                 write(9, 515) (bands%occup(j, irk, is) *dhalf / kpoints%w(irk), &
                      j = n, min(bands%nband(irk, is), n + 6))
              end do
500           format(' k-point:',i5)  
510           format(7f11.4)  
515           format(3x,7(' (',f6.4,')  '))  

           end if
     end if ! exact or iterative diagonalization 
           ! up to here Peter
           !
           call destroy(ham)  
           call myflush(9)  

! davegp
! store convergence info
           idiaconv_total = idiaconv_total .and. idiaconv
! davegp

        end do           ! loop over k-points
     end do        ! loop over spins

     if(iand(pw_params%output(4),2).eq.2) then
        close(37)
     endif       
                
     if(iand(pw_params%output(4),1).eq.1) goto 10  ! added by hjchoi

     !
     !        compute occupation numbers and fermi energy
     !  

!-----------------
! LSDA+U. P. Zhang
!-----------------

     if(pw_params%lsdau.eq.1) then

         call lsdau(imode,iter,k_gspace, energs,pw_params, crys, syms,  &
                    kpoints, bands, wavefn,iscfconv,maxiter)
     end if


     if ((imode == 0).or.(imode == 2)) goto 10   ! leave the loop in this mode

! Note: this code previously called flevel even when doing a band structure
! computation.  This doesn't seem correct to me, and I seem to remember it 
! causing some kind of error, but I can't remember what it was.  --David Roundy

!---------------------------------------
! LSDA+U. P. Zhang
! flevel is called within LSDA+U routine 
!---------------------------------------

     if(pw_params%lsdau.eq.0) then

! fix occupation from file

     if(pw_params%occupy.eq.5) write(9,*) "occupy levels: from_file"

     call flevel(imode, 2, pw_params, energs, crys, bands, kpoints)

     if (iand(pw_params%output(1),1073741824 ) == 1073741824) then

     if(pot_gspace%myproc == 0) then
     open(109,file="occ.dat.old")
     write(109,*) MIN(bands%min(1),bands%min(2)),bands%nspin,bands%nrk

! write-out occupation for all k-points

     do is = 1,bands%nspin
     do irk=1,kpoints%grid(1)*kpoints%grid(2)*kpoints%grid(3)
           nkk=abs(kpoints%kmap(irk))
           do ib = 1, MIN(bands%min(1),bands%min(2))
           write(109,*) 0.5d0*bands%occup(ib,nkk,is)/kpoints%w(nkk)
          end do
        end do
     end do
     close(109)
     end if
     end if
     end if



     if (iand(pw_params%output(1),8)==8) call get_timing(time(8),t_temp)
     !
     !        compute new charge
     !

     call charge(ffts, pw_params, pot_gspace, k_gspace, bands, &
          kpoints, wavefn, den%data(1, 1, 1), energs,crys)


!---------------------------------------------------
! see if we need to symmetry the charge density 
! AFM phase
! P. Zhang
!     if(pw_params%lsdau.eq.1) then
!         call lsdau_AFM_symm(pot_gspace,den%data(1,1,1),crys)
!     end if
!---------------------------------------------------

     if (iand(pw_params%output(1),8)==8) call get_timing(time(9),t_temp)
     !
     !        calculate new screening potential
     !

! davegp
     call velect(1, pot_gspace, crys, ffts, energs, pw_params, &
          denc%data(1, 1, 1), vout%data(1, 1, 1), den%data(1, 1, 1), &
          chdr%data(1, 1, 1))
! restore original potential if we are doing non-self-consistent field
     if( pw_params%nscf ) vout%data = vin%data
! davegp

     if (iand(pw_params%output(1),8)==8) call get_timing(time(10),t_temp)
     !
     !        compute total energy
     !
     if (iter .gt. 1) oldenergy=oldenergs%total

!-----------------------------------------------------------------------
! LSDA+U. P. Zhang
!     call etotal(crys, iter, iscfconv, errscfconv, energs, oldenergs, &
!          pw_params, pot_gspace, vion%data(1, 1, 1), vin%data(1, 1, 1), &
!          vout%data(1, 1, 1), den%data(1, 1, 1))


     call etotal(crys, iter, iscfconv, errscfconv, energs, oldenergs, &
          pw_params, pot_gspace, vion%data(1, 1, 1), vin%data(1, 1, 1), &
          vout%data(1, 1, 1), den%data(1, 1, 1),Etot_corr,ELSDA_corr,Etot_corr_old,Eband_diff)

! davegp
! substitute iscfconv with idiaconv_total
     if( pw_params%nscf ) then
       write(9,*) ' hijacked scf for nscf purposes:'
       write(9,*) '  scf convergence = ', iscfconv
       write(9,*) ' nscf convergence = ', idiaconv_total
       iscfconv = idiaconv_total
     endif
! davego

     if(pw_params%lsdau.eq.1) then 
     if (.not. iscfconv2 .and. iter .ne. maxiter) then
        call restore_LSDA(kpoints,k_gspace, bands, wavefn)
     else
        write(9,*) "Last iteration, keep LSDA+U wavefunctions."
     end if
   
     end if
!--------------------------------------------------------------------

     if (iter .eq. 1)  oldenergy=oldenergs%total
     if (iand(pw_params%output(1),8)==8) call get_timing(time(11),t_temp)
     !
     !        write charge density to disk
     !


! Save the CD file, 
!
!     if (.not. iand(pw_params%output(4), 16) == 16) then 
!       if (pw_params%io_scf .or. iscfconv .or. iter .eq. maxiter) &


     write(9,*) 'Is this nscf? ', pw_params%nscf
     call myflush(9)
     if( .not. pw_params%nscf ) then
         write(9,*) 'Dump the charge density file: CD ', den%nvecs
         call myflush(9)
         call writegsdat(0, pot_gspace, den%data(1, 1, 1) , den%nvecs, &
            den%nvecs, den%nspin, 'CD', 2)
         write(9,*) 'Done with CD'
         call myflush(9)
     endif

     if (iand(pw_params%output(1), 8) == 8) call get_timing(time(12),t_temp)
     !
     !        mix with previous iteration
     !

! davegp
! skip mixing
  if( pw_params%nscf ) goto 1234
! davegp

     if (pw_params%mix_method == 'broyden') then  
       call mixer(iter, errscfconv, pw_params, crys, pot_gspace, mix_gspace, &
          vin%data(1, 1, 1), vout%data(1, 1, 1), energs, vion%data(1, 1, 1), &
          ajac(1))

     else if  ((pw_params%mix_method == 'pulay_kerker' .or. &
                pw_params%mix_method == 'pulay_tf') )  then
       !
       ! if this is not the last step then do mix. This is becuase mixing
       ! requires Gzero to equal 0. This would then mess up stress calc.
       !
       if (.not. iscfconv2 .and. iter .ne. maxiter) then
         !
         ! Set Gzero equal to 0
         !
         do is=1,crys%nspin
           if (pot_gspace%ig0 > 0) vin%data(pot_gspace%ig0,1,is) = zzero
           if (pot_gspace%ig0 > 0) vout%data(pot_gspace%ig0,1,is) = zzero
         end do
         !
         ! remove local pseudopotential from Vin. 
         !
         call mzaxpy(pot_gspace%length*crys%nspin,zmone,vin%data(1,1,1),1,&
               vion%data(1,1,1),1)  

         call mch_pulay(vin%data(1,1,1),vout%data(1,1,1),iter,crys, &
               pot_gspace,R0,w_in0,dw,dR,itmaxpul)

         if  (pw_params%mix_method == 'pulay_kerker' .or. &
                        iter .le. pw_params%numpk_befptf ) then 
           do is = 1, crys%nspin 
             call mch_kerk(vin%data(1,1,is),vout%data(1,1,is),&
               vion%data(1, 1, is),iter,pot_gspace,pw_params%alphamix(1))
           end do
        
         else 
           do is = 1, crys%nspin 
             call Thomas_Fermi(vin%data(1,1,is),vout%data(1,1,is), &
              vion%data(1, 1, is), energs%efermi(is),&
               pw_params%tfw_max_cg,    ffts,& 
             crys,pot_gspace, den%data(1, 1, 1), denc%data(1, 1, 1),&
              vnl_ave%data(1,1,1))
           enddo

         end if 

       end if

     else
       write(9,*) pw_params%mix_method
       call mystop( ' unknown mixing type' )
     end if
 
    write(9,*) ' time for potential mixing ', gimmetime()-t_temp

! davegp
 1234 continue
! davegp

     !
     !        write mixed potential to disk
     !
     if (iand(pw_params%output(1), 8) == 8) call get_timing(time(13),t_temp)

     if (.not. iand(pw_params%output(4), 16) == 16) then 
     
!     if (pw_params%io_scf .or. iscfconv .or. iter .eq. maxiter) &
     call writegsdat(0, pot_gspace, vin%data(1, 1, 1) , vin%nvecs, &
          vin%nvecs, vin%nspin, 'VMIX', 4)
     !
     !        write eigenvectors to disk
     !
!     if (iand(pw_params%output(1), 512) == 512) then  
!        write(9, 130)  
!        if (pot_gspace%myproc == 0) then  
!           call mysystem('rm -f BAND')  
!           call write_band('BAND', bands)  
!        end if
!        do is = 1, crys%nspin  
!           do irk = 1, kpoints%nrk  
!              write(wfnname, 120) irk, is  
!              call writegsdat(0, k_gspace(irk), wavefn%data(1, irk, is), &
!                   bands%nband(irk, is), bands%nband(irk, is), 1, wfnname, 10)
!           end do
!        end do
!     end if
     end if

     if (iand(pw_params%output(1), 8) == 8) call get_timing(time(12),t_temp)
     write(9,*)  num_HPsi,' SCFLOOP # H*Psi  AFTER ITERATION',iter 
! need one more iteration for calculating hubbard-U contribution to forces
     if (iscfconv2) goto 10
     if (iscfconv) iscfconv2=.true.
  end do     ! iterative loop

  write(9, 110) iter - 1  

10 continue  
!     write(34,*)  num_HPsi,' SCFLOOP # H*Psi  AFTER ITERATION',iter 
  !
  !     free memory attached to local structures
  !
  deallocate(vloc)
  
  if(pw_params%mix_method == 'broyden') then  
    call destroy(mix_gspace)
    deallocate(ajac)
  else if (pw_params%mix_method == 'pulay_kerker' .or. &
          pw_params%mix_method == 'pulay_tf'  ) then
    deallocate(R0)
    deallocate(w_in0)
    deallocate(dw)
    deallocate(dR)
  end if

  write(9,*)  gimmetime()-t0,' SCFLOOP TIME  AFTER ITERATION',iter 

  if (iand(pw_params%output(1), 8) == 8)&
   write(9,200) time(1),time(2),time(3),time(4),time(5),time(6),time(7),&
                 time(8),time(9),time(10),time(11),time(12),time(13),time(14),&
                time(1)+time(2)+time(3)+time(4)+time(5)+time(6)+time(7)+&
                 time(8)+time(9)+time(10)+time(11)+time(12)+time(13)+time(14)

!-----------------
! LSDA+U. P. Zhang
!-----------------
 
!     if(pw_params%lsdau.eq.1) then
!         call lsdau_deallocate()
!     end if

  return  
  !
100 format(/,20('-'),' BEGIN OF SELF-CONSISTENT LOOP ',20('-'))  
110 format(20('*'),/' WARNING: SCFLOOP NOT CONVERGED IN ', i3, &
       &     ' ITERATIONS')
120 format('WFN',i5.5,'.',i1)  
130 format(/' <<<<<<<< WRITING WAVEFUNCTIONS TO DISK >>>>>>>>>>')  
131 format(/' <<<<<<<< WRITING INVERSE HESSIAN TO DISK >>>>>>>>>>')  
140 format(' <<<<<<<< READ ',i4,' BANDS FROM FILE ',a,'>>>>>>>>')  
150 format(' *** COULD NOT FIND WAVE FUNCTION ',a,'ON DISK!')  
155 format(' *** COULD NOT FIND BAND FILE ON DISK!')  
160 format(' <<<<<<<< REUSING WAVE FUNCTIONS >>>>>>>>')  

200 format(/' DIAGONALIZATION TIME -  BREAKDOWN:',/1x,27('-'), &
       &           /20x,'ACTIONS',15x,'seconds'&
       &     /,11X,'                       =',          &
       &     /,11X,'FFT time               =',2X,F12.4, &
       &     /,11X,'NL time                =',2X,F12.4, &
       &     /,11X,' occ space time        =',2X,F12.4, &
       &     /,11X,' other cg time         =',2X,F12.4, &
       &     /,11X,' init diag             =',2X,F12.4, &
       &     /,11X,' overlap cg time       =',2X,F12.4, &
       &     /,11X,' mzgemm m*m*len        =',2X,F12.4, &
       &     /,11X,' flevel                =',2X,F12.4, &
       &     /,11X,' charge                =',2X,F12.4, &
       &     /,11X,' velect                =',2X,F12.4, &
       &     /,11X,' etotal                =',2X,F12.4, &
       &     /,11X,' io                    =',2X,F12.4, &
       &     /,11X,' mixing                =',2X,F12.4, &
       &     /,11X,' setup                 =',2X,F12.4, &
       &     /,11X,' total                 =',2X,F12.4)


end subroutine scfloop
