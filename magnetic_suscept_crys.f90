! ----------------------------------------------------------------------
!
!   computes the induced current by perturbation theory
!
!   (1996)      Bernd Pfrommer and Francesco Mauri
!               UC Berkeley
!
!   (1999)      Chris Pickard and Francesco Mauri
!               Paris/Kiel
!               Correct treatment of the non-local derivative |kq><k|
!
!   ----------------------------------------------------------------------

Subroutine magnetic_suscept_crys(qmag,g0mask,pw_params,ffts,bands,pot_gspace,&
     vion,crys,syms,pspot,blop,kpoints)

  Use all_to_all_module
  Include 'use.h'
  Implicit None 
  Include 'interface.h'
  Include 'all_to_all.h' 

  Type (kpoint)               :: kpoints     ! the kpoints
  Type (pw_parameter)         :: pw_params   ! energy cutoffs etc
  Type (parallel_gspace)      :: pot_gspace
  Type (fft_struc)            :: ffts
  Type (crystal)              :: crys
  Type (band)                 :: bands
  Type (symmetry)             :: syms
  Type (pseudo_potential)     :: pspot
  Type (blochl_operator)      :: blop
  Type (complex_gspace_array) :: vion ! the screened ionic potential in gspace

  Real(dp) :: g0mask(3,3),qmag  ! tensor to mask G=0 component,magnitude of q

  !
  !     OUTPUT:
  !     ------
  !
  !
  !     DESCRIPTION:
  !     -----------
  !
  !     Computes the current (j_bare) induced in an insulator
  !     due to an applied magnetic field
  !
  !     Calculates corrections to the chemical shielding 
  !     resulting from the use of non-local pseudopotentials
  !     
  !     Comments refer to the PRB paper of Pickard and Mauri (PM)
  !
  !
  !     the loop structure is as follows
  !
  !
  !    loop over k (kpoints)

  !	generate Gspace for k
  !     setup nonlocal derivative dV/dk @ k
  !	compute u_k 
  !     compute diamagnetic correction
  !     compute the q=0 contributions to the G=0 term     

  !	loop over q^  (wavevector of the applied magnetic field)
  !	    loop over +-q
  !		setup H_k+q
  !               setup nonlocal derivative dV/dk @ k+q
  !		compute u_k+q
  !		loop over the two directions b0 perpendicular to q
  !                    p0= q^ x b0 
  !		     solve for u~_k+q,p0 (cg_blocksolve)
  !                  compute paramagnetic correction
  !		     loop over i ( occupied bands)
  !		        Ftrans u_k,i
  !			loop over p1
  !			   apply p1 to u~_k+q,p0,i
  !			   fourier transform u~_k+q,p0,p1,i
  !         j_bare(r,b0,p1)= j_bare+sign(q)*u_k,i(r)*u~_k+q(r),p0,p1,i
  !			end loop over p1
  !		     end loop over i
  !		loop over b0
  !	    end loop over +-q
  !	end loop over q
  !    end loop over k


  !     -------------- local variables --------------------------------
  !
  Integer :: &
       p0, p1,alf, &         ! counters for directions of momentum operator
       iq,         &         ! qpoint counter
       i,j,        &         ! counter for occupied bands
       irkstart,   &         ! which kpoint to start with
       irk,        &         ! kpoint counter
       nbands_old            ! Store the old no. of bands - so exit as entry
  Integer :: nspin,is,spn,&  ! spin check,spin index,a dummy
       isstart,&             ! spin to start with
       num_HPsi              ! op count of H*Psi

  Integer :: dq(3,3)
  Data dq/1,0,0, 0,1,0, 0,0,1/

  Type (parallel_gspace) :: &
       k_gspace,            &! the gspace for a given k
       ksub_gspace           ! the subspace for a given k


  Type (hamiltonian) ::  ham 
  Type (dense_hamiltonian) :: & ! the subspace hamiltonians 
       densham                  ! for k


  Real(dp) :: &
       t0, tk,tkt, gimmetime, chi_avg, &
       rk(3),             &  ! unshifted kpoint
       kq(3),             &  ! shifted kpoint k+q
       chi_macro(3,3),    &  ! g(G=0, q.ne.0), in ppm/mol [cgs]
       binverse(3,3)         ! inverse of the lattice vectors

  Complex(dp), Pointer :: &
       u_kq(:,:),         &  ! the wave functions for a given k+q
       u_k(:,:),          &  ! the wave functions for a given k
       nloc_deriv_k(:,:,:,:),  & !  derivative of nonlocal potential at k
       nloc_deriv_kq(:,:,:,:), & ! same, but at k+q
       phi(:), phi1(:),   &  ! temporary storage
       vloc(:),           &  ! local potential in realspace
       j_bare(:,:,:),     &  ! bare current in realspace j_bare(r,a0,a1)
       u_ktil(:,:,:),     &  ! u_ktilde in gspace
       u_kqtil(:,:),      &  ! u_kqtilde in gspace
       graduk(:,:),       &  ! the gradient of uk
       u_kqtil_r(:),      &  ! u_kqtilde in rspace
       u_k_rspace(:)         ! u_k in real space

  Complex(dp) ::       &
       delchi,         &
       chivv(3,3),     &     ! g(G=0, q.ne.0)
       chihh(3,3),     &     ! g(G=0, q.ne.0)
       kk_termvv(3,3), &     ! g(G=0, q=0)
       kk_termhh(3,3)        ! g(G=0, q=0)

  Real(dp), Allocatable :: e_kq(:) ! the eigenvalues at u_k+q
  real(dp) time(7) ! 1 fft 2 nl 3occ space 4othercg 5otherscf

  Integer, parameter :: current=3, magnet=3

  !     The diamagnetic augmentation term

  Real(dp), Allocatable :: dia_corr_shift_tot(:,:,:,:),dia_corr_shift(:,:)

  !     The paramagnetic augmentation term

  Real(dp), Allocatable :: para_corr_shift_tot(:,:,:,:),para_corr_shift(:,:,:)

  !
  ! CHI in gspace
  !
  Complex(dp), Allocatable ::     chi_g(:,:,:)


  Logical :: idiaconv          ! if diagonalization has converged


  Integer :: b0, b1, iperm0, iperm1, iqsign, iguess, ierr, nvecs


 !     Wavefunction filename

  Character*11 :: wfnname       

  Integer :: corr_size,kp_count

  External gimmetime

  !     ----------------------------------------------------------------
  !*$*OPTIMIZE(0)

  tkt=gimmetime()

  Write (9,*)
  Write (9,*) '*********************'
  Write (9,*) 'NEW CRYSTAL ALGORITHM'
  Write (9,*) '*********************'

  Allocate(dia_corr_shift_tot(magnet,current,crys%ntype,crys%mxdatm))
  Allocate(dia_corr_shift(crys%ntype,crys%mxdatm))
	
  dia_corr_shift_tot = dzero

  Allocate(para_corr_shift_tot(magnet,current,crys%ntype,crys%mxdatm))
  Allocate(para_corr_shift(current,crys%ntype,crys%mxdatm))

  para_corr_shift_tot = dzero

  ham%ekinmod=pw_params%ekinmod 
  iguess = 1

  !
  ! Summery of optimize flags
  ! ^^^^^^^^^^^^^^^^^^^^^^^^^
  !
  !If(Iand(pw_params%optimize,4).Ne.1)
  !   Do not store the nloc_deriv_k, recalculate when needed
  !         This is called with the flag memory

  !If(Iand(pw_params%optimize,4).Ne.4)
  !   Don't use the value at k as a starting guess for k+q in cgsolve
  !         This is called with the flag nocgsolveguess  
  !

  If(Iand(pw_params%optimize,4).Eq.4) iguess =0

  If(pw_params%epsdiag.Gt.1d-12) Write(9,160) pw_params%epsdiag

  binverse=Transpose(crys%avec)/pi2
  !
  !     allocate the local potential, which is the same for all kpoints
  !
  Allocate(vloc(ffts%r_size))
  !
  !     allocate some other arrays
  !
  Allocate(j_bare(ffts%r_size,3,3))
  j_bare=0 
  chivv=0
  chihh=0
  Allocate(e_kq(bands%max))

  irk=0
  is=1
  !
  !     read checkpoint file if requested
  !
  If(Iand(pw_params%input,4).Eq.4)  Then
     corr_size = crys%ntype*crys%mxdatm
     Call chkpt_j_crys(irk,is,j_bare(1,1,1),&
          9*pot_gspace%r_size, pot_gspace%myproc, chivv,chihh, &
          corr_size,dia_corr_shift_tot,para_corr_shift_tot,pw_params%output(1),2)
  End If

  ! trap the situation when we check point at the last k-point
  ! of spin 1.

  if ( (irk-1) .eq. kpoints%nrk) then
     if ( crys%nspin .eq. 2) THEN
        is=2
        irk=0
     endif
  endif

  irkstart=irk+1

  isstart=is

  kp_count = 0 

  Do is=isstart,crys%nspin       !loop over spin

     !
     !The local potential is spin dependent so we have to set it up for
     !each spin component
     !
     Call setup_local_potential(1,ffts,pot_gspace,vion%data(1,1,is),vloc(1))



     Do irk=irkstart, kpoints%nrk     !loop over kpoints

        !Here we need to trap the case where a spin channel is 
        !not occupied eg H Li Na etc
        !this is a slightly experimental part of the code

        If(bands%ifmax(irk,is) .eq. 0)  EXIT



     kp_count = kp_count + 1

     ! Store the number of bands for this k-point so that 
     ! on exit there is no change

     nbands_old =  bands%nband(irk,is)     

     tk=gimmetime()
     rk=kpoints%rk(:,irk)
     Write(9,800) irk, rk
     Write(9,810) is
     Call myflush(9)
     !
     !    Set up the small gspace for the current kpoint

     !    This gspace is used to represent wavefunctions
     !    For a given direction it has half as many points
     !    as the gspace for the potential
     !    ie k_gspace%length=0.125*pot_gspace%length
     !

     k_gspace%gmax=Sqrt(pw_params%emax)
     k_gspace%rk=kpoints%rk(:,irk)
     k_gspace%nproc=pot_gspace%nproc; ! truly parallel
     k_gspace%name='kmag'
     k_gspace%myproc=pot_gspace%myproc; 
     k_gspace%fftsize=pot_gspace%fftsize
     k_gspace%istar=.False. ; k_gspace%igvec=.False.
     k_gspace%iekin=.True. ; k_gspace%iinverse=.False.
     k_gspace%istruc=.False.
     k_gspace%ipackinfo=.True.

!     Call generate_gspace(pw_params%output(1),k_gspace, crys, syms,1)
     Call generate_gspace(pw_params%output(1),k_gspace, crys, syms)

     ! Set up the corresponding sub-gspace (for start guesses)

     ksub_gspace%gmax=Sqrt(pw_params%emaxsub)
     ksub_gspace%rk=k_gspace%rk ! inherit kpoint from big brother
     ksub_gspace%nproc=1;    ! scatter across all processors
     ksub_gspace%name='kmagsub'
     ksub_gspace%myproc=0; ksub_gspace%fftsize=(/0,0,0,0,0/)
     ksub_gspace%istar=.False. ; ksub_gspace%igvec=.True.
     ksub_gspace%iekin=.True. ; ksub_gspace%iinverse=.False.
     ksub_gspace%istruc=.False.
     ksub_gspace%ipackinfo=.False.

!     Call generate_gspace(pw_params%output(1),ksub_gspace, crys,syms, 1)
     Call generate_gspace(pw_params%output(1),ksub_gspace, crys,syms)

     !       
     !     allocate space for wave function u_k
     !

     Allocate(u_k(k_gspace%length,bands%max))
     Allocate(u_kq(k_gspace%length,bands%max))

     Allocate(phi(k_gspace%length))
     Allocate(phi1(k_gspace%length))

     !
     !     set up the nonlocal part of H and the derivative for given k
     !

     Call create_hamiltonian(ham,k_gspace,pspot,vloc,ffts)

     t0=gimmetime()
     Call setup_nonlocal_potential(ham,pspot,crys,1)
     If(Iand(pw_params%output(1),8).Eq.8) Write(9,922) gimmetime()-t0


     !
     !        generate starting guess for u_k
     !
     ! We do a start guess even if we read in the wavefunction
     ! This is because later on we use the approximate eigenvalues
     ! for the cg_blocksolve and I'm not sure if we can replace them
     ! with the exact ones. Must test

     Call create_dense_hamiltonian(densham,ksub_gspace)
     
     Call start_vectors(0,pw_params,ham,densham,bands%energy(1,irk,is),&
          bands%max,u_k(1,is),pot_gspace,ksub_gspace,vion%data(1,1,is),crys) 

     ! nband is greater than the actual number of bands
     ! so we can safely make nbands slightly smaller (a bit arbitary this)
     bands%nband(irk,is) = bands%nband(irk,is) - 1  

     ierr = 0   

     If(Iand(pw_params%input,1).Eq.1) Then ! read from file

        Call read_band('BAND',bands)
        If(bands%nrk.Le.0) Then ! error
           Write(9,155)
           bands%nrk =kpoints%nrk ; Goto 20
        Endif
        spn=1   !this is just a dummy argument. The wavefunctions
                !are stored as a file per spin. Readgsdat expects
                !them to be stored all in one file.

        Write(wfnname,120) irk,is
        u_k = zzero
        nvecs = bands%nband(irk,is)
        Call readgsdat(1,ierr,k_gspace,u_k(1,1),bands%nband(irk,is),&
             nvecs,spn,wfnname,10)

        If(ierr.Ne.0) Then
           Write(9,150) wfnname 
           Call start_vectors(0,pw_params,ham,densham,bands%energy(1,irk,is),&
          bands%max,u_k(1,is),pot_gspace,ksub_gspace,vion%data(1,1,is),crys)
           bands%nband(irk,is) = bands%nband(irk,is) - 1
           Goto 20
        Endif
        Write(9,140) nvecs, wfnname
        ierr=0
20      Continue
     End If




     If(pw_params%bandgap.Gt.0) Then
        bands%nband(irk,is)=bands%min(is)
        bands%energy(bands%min+1,irk,is)=bands%energy(bands%min(is),irk,is)+&
             pw_params%bandgap
     Endif

     ! 
     ! Diagonalize the hamiltonian to get the wvefunctions u_k
     ! at the current k-point
     ! 

     if (pw_params%diag_method == 'Grassmann') then
        call diagonalize_hamiltonian_gcg(time,1,num_HPsi,pw_params, &
             bands%max,pw_params%maxitdiag,idiaconv,&
             bands%min(is),bands%nband(irk,is),pw_params%epsdiag,&
             bands%energy(1,irk,is),ham, u_k(1,1))
     else
        Call diagonalize_hamiltonian(time,1,num_HPsi,pw_params, &
             bands%max,pw_params%maxitdiag,idiaconv,&
             bands%min(is),bands%nband(irk,is),pw_params%epsdiag,&
             bands%energy(1,irk,is),ham, u_k(1,1))
     endif

     !  Calculate the diamagnetic augmentation term for this kpt 
     t0=gimmetime()
     Call diamagnetic_correction(dia_corr_shift,ham,crys,blop, &
          u_k(1,1),bands%ifmax(irk,is))
     do i=1,3
        dia_corr_shift_tot(i,i,:,:) = dia_corr_shift_tot(i,i,:,:) &
             +dia_corr_shift*kpoints%w(irk)/Real(crys%nspin,dp)
     enddo
     If(Iand(pw_params%output(1),8).Eq.8) Write(9,920) gimmetime()-t0
     !
     !     setup nonlocal derivative at k
     !
     If(Iand(pw_params%optimize,1).Ne.1)  Then
        Allocate(nloc_deriv_k(k_gspace%length+pspot%nanl,pspot%nanl,2,3)) 
        Call setup_nonlocal_derivative(qmag,rk(1),pspot, &
             k_gspace, nloc_deriv_k(1,1,1,1), crys)
     Else
        Allocate(nloc_deriv_k(1,1,1,1))
     Endif
     If(Iand(pw_params%optimize,4).Eq.4) Then
        Allocate(u_ktil(ham%gspace%length,bands%ifmax(irk,is),3))
     Else
        Allocate(u_ktil(ham%gspace%length,bands%ifmax(irk,is),1))
     Endif

     !
     ! The following subroutine calculates elements neccessary
     ! for the G=O term (ie the macroscopic susceptibility)
     ! In particular the q=0 term of equation 64 of PM
     !

     Call magnetic_kkterm(ham,nloc_deriv_k(1,1,1,1),qmag,rk(1),pspot, &
          u_k(1,1),ham%gspace%length,bands%energy(1,irk,is),           &
          bands%ifmax(irk,is),pw_params,crys,done,u_ktil(1,1,1),       &
          kk_termvv(1,1),kk_termhh(1,1))

     !
     If(Iand(pw_params%optimize,4).Ne.4) Deallocate(u_ktil)
     !
     !        free to optimize for memory. then recompute below
     !
     If(Iand(pw_params%optimize,1).Eq.1)  Deallocate(nloc_deriv_k)

     !
     !  Now we compute terms at k=k+q
     !



     Do iq=1, 3    ! loop over all q-points
        Do iqsign=1,-1,-2   ! sign of q   
           kq=kpoints%rk(:,irk)+Real(iqsign,dp)*qmag*Matmul(binverse,dq(:,iq))
           !
           !  Generate starting guess for u_k+q. Use same g-space as for u_k,
           !  but modify its kinetic energy appropriately. Since the 
           !  hamiltonian refers to that g-space, it is altered implicitly, too

           k_gspace%rk=kq; 
           Call compute_ekin_gspace(k_gspace,crys)

           t0=gimmetime()
           Call setup_nonlocal_potential(ham,pspot,crys)
           If(Iand(pw_params%output(1),8).Eq.8) Write(9,922) gimmetime()-t0

           !always use the wavefunction at u_k as a start guess   
           u_kq=u_k
           !but use the eigenvalues of the start guess
           ! (I've not tested why we shouldn't use the exact)
           ! though band gap maybe the answer
           e_kq(1:bands%max)=densham%energy(1:bands%max)
           ! e_kq(1:bands%max)=bands%energy(1,irk,is)

           If(pw_params%bandgap.Gt.0) Then
              e_kq(bands%min(is)+1)=e_kq(bands%min(is))+pw_params%bandgap
           Endif

           !
           ! Diagonalize the hamiltonian to get the wvefunctions u_kq
           ! at k+q
           !
           if (pw_params%diag_method == 'Grassmann') then
              call diagonalize_hamiltonian_gcg(time,1,num_HPsi,pw_params,&
                   bands%max,pw_params%maxitdiag,idiaconv, &
                   bands%nband(irk,is), &
                   bands%nband(irk,is),pw_params%epsdiag,e_kq(1),ham,u_kq(1,1))
           else
              Call diagonalize_hamiltonian(time,1,num_HPsi,pw_params,&
                   bands%max,pw_params%maxitdiag,idiaconv, &
                   bands%nband(irk,is), &
                   bands%nband(irk,is),pw_params%epsdiag,e_kq(1),ham,u_kq(1,1))
           endif
           !
           !           setup nonlocal potentials and derivatives
           !
           Allocate(nloc_deriv_kq(k_gspace%length+pspot%nanl,pspot%nanl, 2,3))
           Call setup_nonlocal_derivative(qmag,kq(1),pspot, &
                k_gspace, nloc_deriv_kq(1,1,1,1), crys)
           !
           !   if memory optimization is switched on, we have to recompute
           !   the nonlocal derivative at k
           If(Iand(pw_params%optimize,1).Eq.1) Then
              Allocate(nloc_deriv_k(k_gspace%length+pspot%nanl,&
                   pspot%nanl,2,3)) 
              Call setup_nonlocal_derivative(qmag,rk(1),&
                   pspot,k_gspace,nloc_deriv_k(1,1,1,1), crys)
           Endif

           Allocate(u_kqtil(k_gspace%length,bands%ifmax(irk,is)))

           !
           ! Loop over the two directions of the applied field
           ! perpendicular to the q direction
           !           compute b0, p0 from q
           !               

           Do iperm0 =1,-1,-2
              b0 = Mod(iq+iperm0+2,3) +1 ! gives different b0 for iperm0
              p0 = Mod(iq-iperm0+2,3) +1 ! gives p0 = abs(q x b0)
              !
              !              compute the gradient of uk
              !
              Allocate(graduk(k_gspace%length,bands%ifmax(irk,is)))
              Do i=1, bands%ifmax(irk,is)
                 Call take_nonloc_deriv_kq_k(p0,k_gspace, rk(1),u_k(1,i),&
                      ham%vnloc%nvecs,nloc_deriv_kq(1,1,1,1),&
                      nloc_deriv_k(1,1,1,1),graduk(1,i),crys) 

              End Do

              If(Iand(pw_params%optimize,1).Eq.1) Then
                 Deallocate(nloc_deriv_kq)
                 Deallocate(nloc_deriv_k)
              Endif

              If(Iand(pw_params%optimize,4).Eq.4) u_kqtil(:,:)=u_ktil(:,:,p0)

              !
              !              compute u_kq~ for a given p0
              !

              ham%shift=zzero   ! no shift here
              Call cg_blocksolve(pw_params%output(1),crys,p0,ham,&
                   pw_params%maxitcgsol,pw_params%epsmag, rk(1), &
                   bands%ifmax(irk,is),bands%energy(1,irk,is),e_kq(1), &
                   u_kq(1,1),u_k(1,1),u_kqtil(1,1),iguess,graduk(1,1),&
                   pw_params%nbandsfft)
              Deallocate(graduk)

              !     
              !     Calculate the paramagnetic augmetation correction
              !
              t0=gimmetime()
              Call paramagnetic_correction_crys(para_corr_shift,ham,  &
                   crys,blop,u_k(1,1),u_kqtil(1,1),bands%ifmax(irk,is),&
                   k_gspace,kpoints%rk(:,irk),kq)

              !        1/2q  because derivative is ( f(q)-f(-q) )/2q

              para_corr_shift_tot(b0,:,:,:) = para_corr_shift_tot(b0,:,:,:)-&
                   Real(iperm0,dp)*Real(iqsign,dp)*para_corr_shift*&
                   kpoints%w(irk)/qmag/dtwo/Real(crys%nspin,dp)
              If(Iand(pw_params%output(1),8).Eq.8) Write(9,921) gimmetime()-t0
              !
              !           re-setup nonlocal potentials and derivatives
              !
              If(Iand(pw_params%optimize,1).Eq.1) Then
                 Allocate(nloc_deriv_kq(k_gspace%length+pspot%nanl,&
                      pspot%nanl,2,3))
                 Call setup_nonlocal_derivative(qmag,kq(1),pspot, &
                      k_gspace, nloc_deriv_kq(1,1,1,1), crys)
                 Allocate(nloc_deriv_k(k_gspace%length+pspot%nanl,&
                      pspot%nanl,2,3))
                 Call setup_nonlocal_derivative(qmag,rk(1),pspot, &
                      k_gspace, nloc_deriv_k(1,1,1,1), crys)
              Endif

              Allocate(u_k_rspace(ffts%r_size))

              Do i=1,bands%ifmax(irk,is) ! loop over occupied bands
                 !
                 !                 Fourier transform u_(k,i)
                 !     
                 Call fourier_transform(-1,ffts,k_gspace,u_k(1,i),&
                      u_k_rspace(1),1)

                 Do alf=1,3 
                    !
                    !                    apply operator p1 to u_k+q~
                    !                  
                    Call take_nonloc_deriv(alf,k_gspace, kq(1),u_kqtil(1,i),0,&
                         nloc_deriv_kq(1,1,1,1),phi(1),crys) 
                    !
                    ! add to j_bare. This is equation 52 of PM
                    !
                    Call add_to_j_bare(ffts, k_gspace,u_k_rspace(1), &
                         phi(1),Real(iperm0,dp)*Real(iqsign,dp)*&
                         kpoints%w(irk)/Real(crys%nspin,dp),j_bare(1,b0,alf))
                 End Do

                 Allocate(u_kqtil_r(ffts%r_size))

                 Call fourier_transform(-1,ffts,k_gspace,u_kqtil(1,i), &
                      u_kqtil_r(1),1)

                 Do alf=1,3 
                    !
                    !                    apply operator p1 to u_k
                    !                  
                    Call take_nonloc_deriv(alf,k_gspace, rk(1),u_k(1,i), &
                         0,nloc_deriv_k(1,1,1,1),phi(1), crys) 

                    !
                    ! add to j_bare. This is equation 52 of PM
                    !
                    Call add_to_j_bare(ffts,k_gspace,u_kqtil_r(1),phi(1), &
                         dmone*Real(iperm0,dp)*Real(iqsign,dp)*kpoints%w(irk) &
                         /Real(crys%nspin,dp),j_bare(1,b0,alf)) 
                 End Do

                 Deallocate(u_kqtil_r)

                 !
                 !                 treat the G=0 component
                 ! These are the remaining components of equation 65
                 !
                 Do iperm1 =1,-1,-2
                    b1 = Mod(iq+iperm1+2,3) +1 ! gives different b1 for iperm
                    p1 = Mod(iq-iperm1+2,3) +1 ! gives p1 = abs(q x b1)

                    !
                    !                    --- half-half for G=0 ----
                    !
                    Call take_nonloc_deriv(p1,k_gspace, rk(1),u_kqtil(1,i),0,&
                         nloc_deriv_k(1,1,1,1),phi(1),crys)  

                    Call take_nonloc_deriv(p1,k_gspace, kq(1),u_kqtil(1,i),0,&
                         nloc_deriv_kq(1,1,1,1),phi1(1), crys)  
                    phi(:)=0.5d0*(phi(:)+phi1(:))

                    delchi= Real(iperm1,dp)*Real(iperm0,dp)*kpoints%w(irk)*&
                         (parallel_zdotc2(k_gspace%length,u_k(1,i),1,&
                         phi(1),1)-kk_termhh(p0,p1))/Real(crys%nspin,dp)
                    If(b0.Eq.b1) Then
                       chihh(b0,b1) = chihh(b0,b1) + delchi
                    Else
                       chihh(b0,b1) = chihh(b0,b1) + dtwo*delchi
                    Endif
                    

                    !
                    !                    --- dH/dk  dH/dk for G=0 ----
                    !

                    !
                    ! Note: we swap k with kq wrt definition in the subroutine
                    !  because we wish to apply the operator in the opposite
                    !  order to the case above for which the 
                    ! routine was written
                    Call take_nonloc_deriv_kq_k(p1,k_gspace,rk(1),&
                         u_kqtil(1,i),pspot%nanl,nloc_deriv_k(1,1,1,1),&
                         nloc_deriv_kq(1,1,1,1),phi(1),crys)

                    delchi= Real(iperm1,dp)*Real(iperm0,dp)*kpoints%w(irk)*&
                         (parallel_zdotc2(k_gspace%length,u_k(1,i),1,&
                         phi(1),1)-kk_termvv(p0,p1))/Real(crys%nspin,dp)
                    If(b0.Eq.b1) Then
                       chivv(b0,b1) = chivv(b0,b1) + delchi
                    Else
                       chivv(b0,b1) = chivv(b0,b1) + dtwo*delchi
                    Endif
                    
                 End Do        ! loop over permutation iperm1


              End Do           ! end of loop over bands i
              Deallocate(u_k_rspace)    
              Call myflush(9)
           End Do              ! loop over permutation iperm0
           Deallocate(u_kqtil)
           Deallocate(nloc_deriv_kq)
           If(Iand(pw_params%optimize,1).Eq.1) Then
              Deallocate(nloc_deriv_k)              
           Endif

        End Do                 ! end of iqsign
     End Do                 ! end of loop over q-points

     Deallocate(phi)
     Deallocate(phi1)
     If(Iand(pw_params%optimize,1).Ne.1) Then
        Deallocate(nloc_deriv_k)              
     Endif
     If(Iand(pw_params%optimize,4).Eq.4) Deallocate(u_ktil)


     Call destroy(densham)
     Call destroy(ksub_gspace)
     Call destroy(k_gspace)
     Deallocate(u_k)
     Deallocate(u_kq) 
     Call destroy(ham)



     Write(9,926) gimmetime()-tk
     !
     !        do checkpointing if required
     !
     If(Iand(pw_params%miscflag,2).Eq.2) then
        corr_size = crys%ntype*crys%mxdatm
        Call chkpt_j_crys(irk,is,&
             j_bare(1,1,1),9*pot_gspace%r_size,pot_gspace%myproc, &
             chivv,chihh,corr_size,dia_corr_shift_tot,para_corr_shift_tot,&
             pw_params%output(1),1)
     end if

     !     Set back number of bands, as on entry

     bands%nband(irk,1) = nbands_old

     if((kp_count==pw_params%nmrkpts).and.(irk<kpoints%nrk)) then
        write(9,*) '************************************************'
        write(9,*) 'Stopping now : max k-points for this run reached'
        write(9,*) '************************************************'
        goto 777 
     end if

  End Do                    ! end of loop over kpoints
End Do                    ! end of spin loop


777 continue

  Deallocate(vloc)
  Deallocate(e_kq)

  Allocate(chi_g(pot_gspace%length,3,3)) 
  !
  !          spin    average q     Hartree 
  !
  chivv= dmtwo * dhalf * dtwo * Real(chivv(:,:),dp)/(qmag*qmag) 
  !
  !     first chi with v-v
  !
  chi_macro=Real(chivv,dp)*(.529177d-8)**3.d0*0.6022d24/(137.036)**2.d0*1.d6 

  Call  symmetrize_tensor(chi_macro(1,1),crys, syms)

  Write(9,190);Write(9,200); Call printdmatrix(chi_macro(1,1),3,3)

  chi_avg = (chi_macro(1,1)+chi_macro(2,2)+chi_macro(3,3))/3.d0 

  Write(9,400) chi_avg

  chi_avg = (chi_macro(1,1)*g0mask(1,1)+chi_macro(2,2)*g0mask(2,2)+&
       chi_macro(3,3)*g0mask(3,3))/3.d0 

  Write(9,410) chi_avg/(crys%vcell*(.529177d-8)**3.d0 * 0.6022d24)*2*pi2

  !
  !     now chi computed with half-half
  !
  chihh=-2.d0*0.5d0*2.d0*Real(chihh(:,:),dp)/(qmag*qmag) 
  chi_macro=Real(chihh,dp)*(.529177d-8)**3.d0*0.6022d24/(137.036)**2.d0*1.d6 

  Call  symmetrize_tensor(chi_macro(1,1),crys,syms)

  Write(9,210);   Call printdmatrix(chi_macro(1,1),3,3)

  chi_avg = (chi_macro(1,1)+chi_macro(2,2)+chi_macro(3,3))/3.d0 

  Write(9,400) chi_avg

  chi_avg = (chi_macro(1,1)*g0mask(1,1)+chi_macro(2,2)*g0mask(2,2)+&
       chi_macro(3,3)*g0mask(3,3))/3.d0 

  Write(9,415) chi_avg/(crys%vcell*(.529177d-8)**3.d0 * 0.6022d24)*2*pi2
  !
  !  compute the  magnetic susceptibilty (local magnetic field)
  !  from the bare current. We need to add in the G=0 term
  !  (the macroscopic susceptibility) We use the hh form
  !
  Write(9,350) (g0mask(i,1:3),i=1,3)
 
  if (iand(pw_params%output(2),1) == 1) then
  ! also write out the current density
  Call wrapup_chi(crys,ffts, pot_gspace, j_bare(1,1,1),qmag,g0mask(1,1),&
       chi_g(1,1,1),chihh(1,1),1)
  else
  Call wrapup_chi(crys,ffts, pot_gspace, j_bare(1,1,1),qmag,g0mask(1,1),&
       chi_g(1,1,1),chihh(1,1),0)
  endif
          
  Call symmetrize_tensor_global(pw_params%output(1),chi_g(1,1,1),syms,pot_gspace)

  Call writegsdat(1, pot_gspace, chi_g(1,1,1), 9,9,1,'CHI',3)

  call symmetrize_tensor_j(dia_corr_shift_tot,crys, syms) 
  call symmetrize_tensor_j(para_corr_shift_tot,crys, syms)

  !
  ! Write out the chemical shieldings
  !

  Call nmr_shift_new(pw_params,crys,pot_gspace,chi_g(1,1,1),&
       dia_corr_shift_tot,para_corr_shift_tot)

  Write(9,420) qmag
  Write(9,300)

  Deallocate(j_bare,dia_corr_shift,dia_corr_shift_tot,para_corr_shift,para_corr_shift_tot)
  Deallocate(chi_g)
  Write(9,927) gimmetime()-tkt

  Return

120 Format('WFN',i5.5,'.',i1)
155 Format(' *** COULD NOT FIND BAND FILE ON DISK!')
140 Format(' <<<<<<<< READ ',i4,' BANDS FROM FILE ',a,'>>>>>>>>')
150 Format(' *** COULD NOT FIND WAVE FUNCTION ',a,'ON DISK!')
160 Format(/' *** WARNING: ACCURACY OF DIAGONALIZATION ',&
       'MIGHT NOT BE HIGH ENOUGH:',g12.6)

190 Format(/' *** WARNING: CHI IS ONLY CORRECT IF SYMMETRY-',&
       /' *** OPERATIONS MAP CARTESIAN AXIS INTO EACH OTHER')

200 Format(/' MACROSCOPIC CHI (nonlocal velocity operator)',&
       /' [10^-6 cm^3/mole]:'/)
210 Format(/' MACROSCOPIC CHI (half-half)',/' [10^-6 cm^3/mole]:'/)
220 Format(/' G=0 CORRECTION TO THE NMR SHIFT','(nonlocal velocity operator):')
222 Format(/' G=0 CORRECTION TO THE NMR SHIFT','(half-half):')
300 Format(/' THE NMR SHIFT IS COMPUTED USING THE HALF/HALF ',&
       ' FORM FOR THE G=0 COMPONENT')
350 Format(/' THE G=0 COMPONENT FOR THE SHIFT IS MASKED WITH:'/,3(3f12.6/))
400 Format(/' AVERAGE CHI:', f12.4,' *10^-6 cm^3/mole')

410 Format(/' VV SHIFT CONTRIBUTION:', f12.4)
415 Format(/' HH SHIFT CONTRIBUTION:', f12.4)
420 Format(/' NMR Q:', f10.6)

800 Format(//' COMPUTING CONTRIBUTION FROM KPOINT',i4,':',3f9.4)
810 Format(//' COMPUTING CONTRIBUTION FROM SPIN',i4)

910 Format(' TIME FOR DIAGONALIZATION:', f12.3)
920 Format(' TIME FOR DIAMAGNETIC CORRECTON:', f12.3)
921 Format(' TIME FOR PARAMAGNETIC CORRECTION:', f12.3)
922 Format(' TIME FOR NONLOCAL POTENTIAL:', f12.3)
926 Format(' TIME FOR KPOINT:', f12.3)
927 Format(' TOTAL TIME FOR NMR:', f12.3, 'SECONDS')
930 Format(' TIME FOR ADD_TO_J_BARE:',f12.3)

End Subroutine magnetic_suscept_crys

