! -------------------------------------------------------------------------
!
!   Computes the magnetic susceptibility in the LDA by perturbation theory
!
!   (1996)      Bernd Pfrommer and Francesco Mauri
!               UC Berkeley
!
!   Taken from "magnetic_suscept" -- modified, taken to the molecular limit
!   Using the equation 8 of Gregor et al 1999
!
!   (1999)      Chris Pickard and Francesco Mauri
!               Paris/Kiel
! -------------------------------------------------------------------------
!     @process extchk
  
Subroutine magnetic_suscept_eqn8(qmag,g0mask,pw_params,ffts,bands,&
     pot_gspace,vion,crys,syms,pspot,blop,kpoints)

  Use all_to_all_module
  Include 'use.h'
  Implicit None 
  Include 'interface.h'
  Include 'all_to_all.h'

  Type (kpoint)           :: kpoints     ! the kpoints
  Type (pw_parameter)     :: pw_params ! energy cutoffs etc
  Type (parallel_gspace)  :: pot_gspace
  Type (fft_struc)        :: ffts
  Type (crystal)          :: crys
  Type (band)             :: bands
  Type (symmetry)         :: syms
  Type (pseudo_potential) :: pspot
  Type (blochl_operator)  :: blop
  Type (complex_gspace_array) :: vion ! Screened ionic potential in gspace

  Real(dp) :: g0mask(3,3) ! Tensor to mask G=0 component
  Real(dp) :: qmag

  !
  !     OUTPUT:
  !     ------
  !
  !     DESCRIPTION:
  !     -----------
  !
  !     Computes the diamagnetic susceptibility of insulators.
  !
  !     -------------- Local variables ---------------
  !
  Integer :: alpha,beta,gamma ! Counters for directions of momentum operator
  Integer :: i,j         ! Counter for occupied bands
  Integer :: irkstart    ! Which kpoint to start with
  Integer :: irk         ! kpoint counter
  Integer :: nbands_old  ! Store the old no. of bands - so exit as entry
  Integer ::nt,iat       ! type and ion of type counters
  Integer :: nspin       ! spin check 

  Type (parallel_gspace)   :: k_gspace    ! The gspace for a given k
  Type (parallel_gspace)   :: ksub_gspace ! The subspace for a given k

  Type (hamiltonian)       :: ham     ! The Hamiltonian
  Type (dense_hamiltonian) :: densham ! The subspace Hamiltonian 

  Real(dp) :: t0,tk,tkt,gimmetime
  Real(dp) :: rk(3)               ! Unshifted kpoint
  Real(dp) :: chi_avg             ! Average chi
  Real(dp) :: chi_0(3,3)          ! G=0 contribution
  Real(dp) :: chi_macro(3,3)      ! g(G=0, q.ne.0), in ppm/mol [cgs]
  Real(dp) :: RI(3)               ! position vector

  Complex(dp) :: csum   ! Complex sum

  Complex(dp), Pointer :: u_k(:,:)       ! The wave functions for a given k
  Complex(dp), Pointer :: nloc_deriv_k(:,:,:,:) !Deriv of nl pot at k
  Complex(dp), Pointer :: phi(:),phi1(:) ! Temporary storage
  Complex(dp), Pointer :: vloc(:)        ! Local potential in realspace
  Complex(dp), Pointer :: J_r(:,:,:)     ! J_mag in rspace J_r(r,a0,a1)
  Complex(dp), Pointer :: J2_r(:,:,:)    ! J2_mag in rspace J2_r(r,a0,a1)
  Complex(dp), Pointer :: u_ktil(:,:)    ! u_ktilde in gspace
  Complex(dp), Pointer :: u_ktil2(:,:)   ! u_ktilde2 in gspace
  Complex(dp), Pointer :: graduk(:,:)    ! The gradient of uk
  Complex(dp), Pointer :: u_ktil_r(:)    ! u_ktilde in rspace
  Complex(dp), Pointer :: u_ktil_r2(:)   ! u_ktilde2 in rspace
  Complex(dp), Pointer :: u_k_rspace(:)  ! u_k in real space  
  Complex(dp), Pointer :: rinr(:,:)      ! The position operator  
  Complex(dp), Pointer :: work(:),rJ(:)

  Logical :: idiaconv  ! If diagonalization has converged
  Integer :: num_HPsi              ! op count of H*Psi
  real(dp) time(6) ! 1 fft 2 nl 3occ space 4othercg 5otherscf

  Integer :: indx(-1:+1,3),indx2(3,3),ierr,nvecs

  ! The charge augmentation term
  Integer, parameter :: current=3, magnet=3

  Real(dp),    Allocatable :: dia_corr_shift_tot(:,:,:,:),dia_corr_shift(:,:)

  ! The charge augmentation term

  Real(dp), Allocatable :: para_corr_shift_tot(:,:,:,:),para_corr_shift(:,:,:)
  Real(dp), Allocatable :: para_corr_shift_2(:,:,:,:)

  ! CHI in gspace
  !
  Complex(dp), Allocatable ::     chi_g(:,:,:)


  ! Wavefunction filename

  Character*11 :: wfnname           

  External gimmetime

  ! ----------------------------------------------------------------

  tkt=gimmetime()

  Write (9,*)
  Write (9,*) '***************************'
  Write (9,*) 'MOLECULAR ALGORITHM (Eqn 8)'
  Write (9,*) '***************************'

  ! Allocate the augmentation arrays

  Allocate(dia_corr_shift_tot(magnet,current,crys%ntype,crys%mxdatm))
  Allocate(dia_corr_shift(crys%ntype,crys%mxdatm))

  dia_corr_shift_tot = dzero

  Allocate(para_corr_shift_tot(magnet,current,crys%ntype,crys%mxdatm))
  Allocate(para_corr_shift_2(magnet,current,crys%ntype,crys%mxdatm))
  Allocate(para_corr_shift(current,crys%ntype,crys%mxdatm))

  para_corr_shift_tot = dzero
  para_corr_shift_2   = dzero

  !     Set indexing for vector product

  indx(+1,1) = 3 ; indx(-1,1) = 2
  indx(+1,2) = 1 ; indx(-1,2) = 3
  indx(+1,3) = 2 ; indx(-1,3) = 1

  !     Similar - for application of B x r to rho

  indx2 = 0
  indx2(1,2) = -3 ; indx2(1,3) = +2
  indx2(2,1) = +3 ; indx2(2,3) = -1
  indx2(3,1) = -2 ; indx2(3,2) = +1

  !     Build the position operator

  Allocate(rinr(ffts%r_size,3))

  Call position_in_rspace(ffts,pot_gspace,rinr(1,1),crys) 

  ham%ekinmod=pw_params%ekinmod 

  !
  !     Set up the local potential, which is the same for all kpoints
  !
  Allocate(vloc(ffts%r_size))
  Call setup_local_potential(1,ffts,pot_gspace,vion%data(1,1,1),vloc(1))
  !
  !     Allocate and initialise some other arrays
  !
  Allocate(J_r(ffts%r_size,3,3),J2_r(ffts%r_size,3,3))
  Allocate(work(ffts%r_size),rJ(ffts%r_size))
  J_r  = zzero
  J2_r = zzero

  irk=0
  !
  !     Read checkpoint file if requested
  !
  If(Iand(pw_params%input,4).Eq.4)  &
       Call chkpt_chi_mol2(irk,J_r(1,1,1),J2_r(1,1,1),&
       9*pot_gspace%r_size,pot_gspace%myproc,pw_params%output(1),2)

  irkstart=irk+1

  Do irk = irkstart, kpoints%nrk ! start of k-point loop

     ! Store the number of bands for this k-point so that 
     ! on exit there is no change     

     nbands_old =  bands%nband(irk,1)

     tk=gimmetime()

     rk=kpoints%rk(:,irk)
     Write(9,800) irk, rk
     Call myflush(9)
     !
     !     Set up the sub-gspace and the gspace for a given kpoint irk
     !
     k_gspace%gmax     = Sqrt(pw_params%emax)
     k_gspace%rk       = kpoints%rk(:,irk)
     k_gspace%nproc    = pot_gspace%nproc; ! truly parallel
     k_gspace%name     = 'kmag'
     k_gspace%myproc   = pot_gspace%myproc 
     k_gspace%fftsize  = pot_gspace%fftsize
     k_gspace%istar    = .False. ; k_gspace%igvec = .False.
     k_gspace%iekin    = .True. ; k_gspace%iinverse = .False.
     k_gspace%istruc   = .False.
     k_gspace%ipackinfo= .True.

!     Call generate_gspace(pw_params%output(1),k_gspace,crys,syms,1)
     Call generate_gspace(pw_params%output(1),k_gspace,crys,syms)
     !       
     !     Allocate space for wave function u_k
     !
     Allocate(u_k(k_gspace%length,bands%max))
     Allocate(phi(k_gspace%length),phi1(k_gspace%length))

     !
     !     Now the sub-g-space
     !         
     ksub_gspace%gmax      = Sqrt(pw_params%emaxsub)
     ksub_gspace%rk        = k_gspace%rk ! inherit kpoint from big brother
     ksub_gspace%nproc     = 1;    ! scatter across all processors
     ksub_gspace%name      = 'kmagsub'
     ksub_gspace%myproc    = 0; ksub_gspace%fftsize = (/0,0,0,0,0/)
     ksub_gspace%istar     = .False. ; ksub_gspace%igvec   = .True.
     ksub_gspace%iekin     = .True. ; ksub_gspace%iinverse = .False.
     ksub_gspace%istruc    = .False.
     ksub_gspace%ipackinfo = .False.

!     Call generate_gspace(pw_params%output(1),ksub_gspace,crys,syms,1)
     Call generate_gspace(pw_params%output(1),ksub_gspace,crys,syms)
     !
     !     Set up the nonlocal part of H and the derivative for given k
     !
     Call create_hamiltonian(ham,k_gspace,pspot,vloc,ffts)

     t0=gimmetime()
     Call setup_nonlocal_potential(ham,pspot,crys,1)
     If(Iand(pw_params%output(1),8).Eq.8) Write(9,922) gimmetime()-t0
     !     
     !     Generate starting guess for u_k
     !
     ierr = 1    ! Assume we can't read in any wavefunctions

        nspin=1
     If(Iand(pw_params%input,1).Eq.1) Then ! read from file

        Call read_band('BAND',bands)
        If(bands%nrk.Le.0) Then ! error
           Write(9,155)
           bands%nrk =kpoints%nrk ; Goto 20
        Endif
        

        Write(wfnname,120) irk,nspin
        u_k = zzero
        nvecs = bands%nband(irk,1)
        Call readgsdat(1,ierr,k_gspace,u_k(1,1),bands%nband(irk,1),&
             nvecs,nspin,wfnname,10)

        If(ierr.Ne.0) Then
           Write(9,150) wfnname ; Goto 20
        Endif
        Write(9,140) nvecs, wfnname
        ierr=0
20      Continue
     End If

     !     Looks like we have to calculate the starting guess

     If(ierr.Eq.1) Then

        Call create_dense_hamiltonian(densham,ksub_gspace)

        Call start_vectors(0,pw_params,ham,densham,&
             bands%energy(1,irk,1),bands%max,u_k(1,1),&
             pot_gspace,ksub_gspace,vion%data(1,1,1),crys)

     !     Optimise for insulators

     bands%nband(irk, 1) = bands%nband(irk, 1) - 1  

     End If
     !
     !     Compute the u_k
     !      


     If(pw_params%bandgap.Gt.0) Then
        bands%nband(irk,1)=bands%min(1)
        bands%energy(bands%min(1)+1,irk,1)=&
             bands%energy(bands%min(1),irk,1)+pw_params%bandgap
     Endif

     Call diagonalize_hamiltonian(time,1,num_HPsi,pw_params, &
          bands%max,pw_params%maxitdiag,idiaconv,&
          bands%min(1), bands%nband(irk,1),pw_params%epsdiag,&
          bands%energy(1,irk,1),ham, u_k(1,1))

     !
     !     Calculate the charge augmentation term for this kpt 
     !
     Call diamagnetic_correction(dia_corr_shift,ham,crys,blop, &
          u_k(1,1),bands%ifmax(irk,1))
     do i=1,3
        dia_corr_shift_tot(i,i,:,:) = dia_corr_shift_tot(i,i,:,:)&
             +dia_corr_shift*kpoints%w(irk)
     enddo

     Allocate(u_ktil(k_gspace%length,bands%ifmax(irk,1)))
     Allocate(u_ktil2(k_gspace%length,bands%ifmax(irk,1)))

     !
     !     Allocate and setup the non-local derivative
     !

     Allocate(nloc_deriv_k(k_gspace%length+pspot%nanl,pspot%nanl,2,3)) 
     Call setup_nonlocal_derivative(qmag,rk(1),pspot,&
          k_gspace, nloc_deriv_k(1,1,1,1),crys)

     !     
     !     Loop over the magnetic field directions
     !

     Do beta =1,3            

        !     The first perturbation -> u_ktil

        !                           
        !     Compute [r x p + DVnl_I x R_I]_beta|u_k>
        !     
        Allocate(graduk(k_gspace%length,bands%ifmax(irk,1)))

        Do i=1, bands%ifmax(irk,1)

           graduk(:,i) = zzero

           !     First [r x p]_beta|u_k>

           Do j=1,-1,-2
              Call apply_mom(indx(j,beta),k_gspace,rk(1),u_k(1,i),&
                   phi(1),crys)  
              Call apply_r2psi(indx(-j,beta),k_gspace,ffts,rinr(1,1),&
                   phi(1),phi1(1))
              graduk(:,i) = graduk(:,i) + Real(j,dp)*phi1(:)
           End Do

           !     Now [SUM_I R_I x DVnl_I]_beta|u_k>

           Call apply_rixdvnl(beta,k_gspace,u_k(1,i),ham%vnloc%nvecs,&
                nloc_deriv_k(1,1,1,1),phi(1),crys,pspot)

           graduk(:,i) = graduk(:,i) + phi(:)
           

        End Do

        !     
        !     Compute u_k~ for a given beta
        !       
        ham%shift=dzero      ! no shift here

        Call cg_blocksolve_mol(pw_params%output(1),crys,ham,&
             pw_params%maxitcgsol,pw_params%epsmag,bands%ifmax(irk,1),&
             bands%energy(1,irk,1),u_k(1,1),u_ktil(1,1),graduk(1,1),&
             pw_params%nbandsfft)

        !
        !     The second perturbation -> u_ktil2 : 
        !

        !                           
        !     Compute v_beta|u_k>
        !                    
        Do i=1, bands%ifmax(irk,1)

           graduk(:,i) = zzero

           !     First p_beta|u_k>

           Call apply_mom(beta,k_gspace,rk(1),u_k(1,i),phi(1),crys)  

           graduk(:,i) = graduk(:,i) + phi(:)

           !     Now DVnl_beta|u_k>

           Call apply_dvnl(beta,k_gspace,u_k(1,i),ham%vnloc%nvecs,&
                nloc_deriv_k(1,1,1,1),phi(1),crys)

           graduk(:,i) = graduk(:,i) + phi(:)

        End Do

        !     
        !     Compute u_k~2 for a given beta
        !     
        ham%shift=dzero      ! no shift here

        Call cg_blocksolve_mol(pw_params%output(1),crys,ham,&
             pw_params%maxitcgsol,pw_params%epsmag,bands%ifmax(irk,1),&
             bands%energy(1,irk,1),u_k(1,1),u_ktil2(1,1),graduk(1,1),&
             pw_params%nbandsfft)

        Deallocate(graduk)

        Allocate(u_k_rspace(ffts%r_size))

        Do i=1,bands%ifmax(irk,1) ! Loop over occupied bands

           !     
           !     Fourier transform u_(k,i)
           !     
           Call fourier_transform(-1,ffts,k_gspace,u_k(1,i),u_k_rspace(1),1)


           Do alpha=1,3 
              !     
              !     Apply operator p to u_k~
              !     
              Call apply_mom(alpha,k_gspace,rk(1),u_ktil(1,i),phi(1),crys) 
              !     
              !     Add to the "charge" density
              !     
              Call add_to_j_bare(ffts, k_gspace,u_k_rspace(1),phi(1),&
                   kpoints%w(irk),J_r(1,beta,alpha))
              !     
              !     Apply operator p to u_k~2
              !     
              Call apply_mom(alpha,k_gspace,rk(1),u_ktil2(1,i),phi(1),crys) 
              !     
              !     Add to the "charge" density
              !     
              Call add_to_j_bare(ffts, k_gspace,u_k_rspace(1),phi(1),&
                   kpoints%w(irk),J2_r(1,beta,alpha))

           End Do

           Allocate(u_ktil_r(ffts%r_size),u_ktil_r2(ffts%r_size))

           Call fourier_transform(-1,ffts,k_gspace,u_ktil(1,i),u_ktil_r(1),1)

           Call fourier_transform(-1,ffts,k_gspace,u_ktil2(1,i),u_ktil_r2(1),1)

           Do alpha=1,3 
              !     
              !     Apply operator p to u_k
              !     
              Call apply_mom(alpha,k_gspace,rk(1),u_k(1,i),phi(1),crys) 

              !     
              !     Add to the "charge" density
              !     
              Call add_to_j_bare(ffts,k_gspace,u_ktil_r(1),phi(1),&
                   kpoints%w(irk),J_r(1,beta,alpha)) 

              Call add_to_j_bare(ffts,k_gspace,u_ktil_r2(1),phi(1),&
                   kpoints%w(irk),J2_r(1,beta,alpha)) 

           End Do

           Deallocate(u_ktil_r,u_ktil_r2)

        End Do              ! End of loop over bands, i
        Deallocate(u_k_rspace)    
        Call myflush(9)

        !     
        !     Calculate the current augmetation correction
        !
        Call paramagnetic_correction_mol(para_corr_shift,ham,crys, &
             blop,u_k(1,1),u_ktil(1,1),bands%ifmax(irk,1),(1.d0,0.d0))

        para_corr_shift_tot(beta,:,:,:) = &
             para_corr_shift_tot(beta,:,:,:)+para_corr_shift*kpoints%w(irk)

        Call paramagnetic_correction_mol(para_corr_shift,ham,crys, & 
             blop,u_k(1,1),u_ktil2(1,1),bands%ifmax(irk,1),(1.d0,0.d0))

        para_corr_shift_2(beta,:,:,:) =  &
             para_corr_shift_2(beta,:,:,:)+para_corr_shift*kpoints%w(irk)

     End Do                 ! Loop over beta

     Deallocate(u_k,u_ktil,u_ktil2,phi,phi1,nloc_deriv_k)

     Call destroy(densham)
     Call destroy(ksub_gspace)
     Call destroy(k_gspace)
     Call destroy(ham)

     Write(9,926) gimmetime()-tk
     !
     !     Do checkpointing if required
     !
     If(Iand(pw_params%miscflag,2).Eq.2) Then
        Call chkpt_chi_mol2(irk,J_r(1,1,1),J2_r(1,1,1),9*pot_gspace%r_size,&
             pot_gspace%myproc,pw_params%output(1),1)
     Endif


     !     Set back number of bands, as on entry

     bands%nband(irk,1) = nbands_old

  End Do   ! End of k-point loop 


  !     ---------- Add the "charge" term to the current -----------


  Do beta=1,3
     Do alpha=1,3

        rJ = zzero
        Do j=-1,1,+2
           Call apply_r2rho(indx(-j,beta),ffts,rinr(1,1),&
                J2_r(1,indx(j,beta),alpha),work(1))
           rJ = rJ + real(j,dp)*work
        End Do

        J_r(:,beta,alpha)=J_r(:,beta,alpha) - rJ(:)

     End Do
  End Do


  J_r = -ztwo*J_r ! Factor of two for Hartree <--> Rydberg

  !
  !     ------ Find the macroscopic chi ------
  !

  Do beta=1,3
     Do alpha=1,3

        rJ = zzero
        Do j=1,-1,-2
           Call apply_r2rho(indx(j,alpha),ffts,rinr(1,1),&
                J_r(1,beta,indx(-j,alpha)),work(1))
           rJ = rJ - real(j,dp)*work
        End Do

        !     
        !     ------- Parallel sum --------

        csum = zzero
        Do j=1,ffts%r_size
           csum = csum + rJ(j)
        End Do
        Call all_sum_all(csum)
        j = ffts%r_size
        Call all_sum_all(j)
        csum = csum/real(j,dp)
        chi_0(beta,alpha) = Real(csum)*0.5d0 ! Factor of c it put
        ! in later - nmr_shift

     End Do
  End Do

  ! Convert to units of 10^6 cm**3/mole (note - factor of c)

  chi_macro=chi_0*(.529177d-8)**3.d0*0.6022d24*1.d6/(137.036)**2.d0

  !call symmetrize_tensor(chi_macro(1,1),crys,syms)

  Write(9,190)
  Write(9,200)

  Call printdmatrix(chi_macro(1,1),3,3)

  chi_avg = (chi_macro(1,1)+chi_macro(2,2)+chi_macro(3,3))/3.d0 

  Write(9,400) chi_avg

  chi_avg = (chi_macro(1,1)*g0mask(1,1)+chi_macro(2,2)*g0mask(2,2)+&
       chi_macro(3,3)*g0mask(3,3))/3.d0

  Write(9,415) chi_avg/(crys%vcell*(.529177d-8)**3.d0*0.6022d24)*2*pi2

  !
  !     ------ Construct chi in reciprocal space ------
  !
  Allocate(chi_g(pot_gspace%length,3,3))
  Call j2chi_mol(crys,ffts,pot_gspace,J_r(1,1,1),chi_g(1,1,1),&
       chi_0(1,1),g0mask(1,1))

  !call symmetrize_tensor_global(pw_params%output(1),chi_g(1,1,1),syms,pot_gspace)

  Call writegsdat(1,pot_gspace,chi_g(1,1,1),9,9,nspin,'CHI',3)

  !
  !     Combine the two current-like augmentation terms
  !
  Do nt=1,crys%ntype
     Do iat=1,crys%natom(nt)

        RI = Matmul(crys%avec,crys%rat(:,iat,nt))/pi2

        para_corr_shift_tot(1,:,nt,iat) = para_corr_shift_tot(1,:,nt,iat) + &
             para_corr_shift_2(2,:,nt,iat)*RI(3) &
             -para_corr_shift_2(3,:,nt,iat)*RI(2)
        para_corr_shift_tot(2,:,nt,iat) = para_corr_shift_tot(2,:,nt,iat) + &
             para_corr_shift_2(3,:,nt,iat)*RI(1) &
             -para_corr_shift_2(1,:,nt,iat)*RI(3)
        para_corr_shift_tot(3,:,nt,iat) = para_corr_shift_tot(3,:,nt,iat) + &
             para_corr_shift_2(1,:,nt,iat)*RI(2) &
             -para_corr_shift_2(2,:,nt,iat)*RI(1)

     End Do
  End Do

  !
  !     Write out the nmr chemical shifts
  !     

  Call nmr_shift_new(pw_params,crys,pot_gspace,chi_g(1,1,1),&
       dia_corr_shift_tot,para_corr_shift_tot)

  Deallocate(vloc,J_r,rinr,rJ,dia_corr_shift,dia_corr_shift_tot)
  Deallocate(para_corr_shift_tot,para_corr_shift_2,para_corr_shift,work)
  Deallocate(chi_g) 
  Write (9,927) gimmetime()-tkt

  Return

  !     -------------------------------------------------------------

120 Format('WFN',i5.5,'.',i1)
155 Format(' *** COULD NOT FIND BAND FILE ON DISK!')
140 Format(' <<<<<<<< READ ',i4,' BANDS FROM FILE ',a,'>>>>>>>>')
150 Format(' *** COULD NOT FIND WAVE FUNCTION ',a,'ON DISK!')
190 Format(/' *** WARNING: CHI IS ONLY CORRECT IF SYMMETRY-',&
       /' *** OPERATIONS MAP CARTESIAN AXIS INTO EACH OTHER')
200 Format(/' MACROSCOPIC CHI',/' [10^-6 cm^3/mole]:'/)
400 Format(/' AVERAGE CHI:', f12.4,' *10^-6 cm^3/mole')
415 Format(/' HH SHIFT CONTRIBUTION:', f12.4)
800 Format(//' COMPUTING CONTRIBUTION FROM KPOINT',i4,':',3f9.4)
922 Format(' TIME FOR NONLOCAL POTENTIAL:', f12.3)
923 Format(' TIME FOR BLOCHL OPERATOR:', f12.3)
927 Format(' TOTAL TIME FOR NMR:', f12.3, 'SECONDS')
926 Format(' TIME FOR KPOINT:', f12.3)

End Subroutine magnetic_suscept_eqn8

