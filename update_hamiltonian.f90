!-*-Fortran-*-
!
subroutine update_hamiltonian(kp, ham, ldawfn, len, neig, wfn, &
     pot_gspace, crys, energs, pw_params, denc, rden, den, vion, &
     vbareion,vin,time,t3)
  !
  !     1996 Bernd Pfrommer
  !
  use all_to_all_module  
  include 'use.h'  
  implicit none                 ! never uncomment this line.
  include 'interface.h'  
  include 'all_to_all.h'  
  !
  !     INPUT:
  !     -----
  !
  type(crystal), intent(in) :: crys  
  type(kpoint), intent(in) :: kp  
  integer :: &
       len(kp%nrk), &           ! length of the different gspaces
       neig(crys%nspin), &                  ! number of occupied wave functions
       ldawfn                   ! leading dimension of wavefn array
  type(parallel_gspace), intent(in) :: &
       pot_gspace               ! gspace for charge density and potential
 
  type(pw_parameter), intent(in) :: pw_params   
  complex(dp), intent(in) :: &
       vbareion(pot_gspace%length,crys%nspin), &  ! bare ionic potential
       denc(pot_gspace%length), &  ! the core charge density
       wfn(ldawfn, kp%nrk,crys%nspin)      ! the wave functions in gspace
  !
  !     INPUT/OUTPUT:
  !     -------------
  !
  type(hamiltonian), intent(inout) :: &
       ham(kp%nrk,crys%nspin)              ! the Hamiltonian to be updated
  complex(dp), intent(inout):: &
       vin(pot_gspace%length,crys%nspin)     ! old potenital
  type(energy), intent(inout) :: energs
  real(dp), intent(inout) :: time(num_time) 
  !
  !     OUTPUT:
  !     ------
  !
  complex(dp), intent(out) :: &
       den(pot_gspace%length,crys%nspin)   ! charge density in fourier space
  !
  !     WORK ARRAYS:
  !     -----------
  !
  complex(dp) :: &
       rden(pot_gspace%r_size,crys%nspin), &  ! realspace work array rden(r_size)
       vion(pot_gspace%length,crys%nspin)     ! work array for vion
  !
  !     DESCRIPTION:
  !     -----------
  !
  !     Updates the Hamiltonian to the new wave functions
  !
  !
  !     It is assumed that all Hamiltonians have the same realspace
  !     grid  size for their representation, and it matches the
  !     realspace grid for the potential.
  !
  !     Further, the Hamiltonians are expected to have only a single
  !     local potential, since it doesn't vary with the kpoint. This
  !     means that the vloc pointer of all Hamiltonians points to the
  !     same location.
  !  
  !     --------------- local variables -------------------------
  !
  complex(dp) :: cdtot  
  integer :: n, i, ntot, ids, irk,t3 , iwfn,rlen, nfn, is ,nspin 
  real(dp) errscfconv,spinf 
  real(dp), external :: gimmetime 
  !
  !     ------------- compute new charge -------------------------
  !
  ntot = ham(1,1)%gspace%fftsize(1) * ham(1,1)%gspace%fftsize(2) * &
       ham(1,1)%gspace%fftsize(3)
  
  nspin=crys%nspin
  rden(1:ham(1,1)%fftstruc%r_size,1:nspin) = zzero
  cdtot = zzero

  rlen = ham(1,1)%fftstruc%r_size   ! total number of points in realspace block

  if (crys%nspin >= 2) then
     spinf = dhalf
  else
     spinf = done 
  end if

  do is=1,nspin
   do irk = 1, kp%nrk  

     do n = 1, neig(is),pw_params%nbandsfft        ! loop over occupied bands
        nfn=min(pw_params%nbandsfft, neig(is) + 1 - n)
        !
        !           fourier transform to realspace
        !
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(5),t3)
        call fourier_transform(-1, ham(1,1)%fftstruc, ham(irk,is)%gspace, &
            wfn(1 + (n-1)*len(irk), irk,is),ham(1,1)%fftstruc%rspacebuf(1),nfn)
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(1),t3)
        !
        !           square and add to the charge density
        !
        do iwfn=0,nfn-1
          do i = 1, rlen  
            rden(i,is)=rden(i,is)+ham(1,1)%fftstruc%rspacebuf(i+iwfn*rlen) * &
               conjg(ham(1,1)%fftstruc%rspacebuf(i+iwfn*rlen)) &
                 * dtwo*kp%w(irk)*spinf 
          end do
        end do
        !            call all_sum_all(cdtot)
        !            write(9,*) irk,n,cdtot/dble(ntot)
      end do
   end do
   !
   !     fourier transform back to momentum space
   !
   if (iand(pw_params%output(1), 8) == 8) call get_timing(time(5),t3)
   call fourier_transform(1,ham(1,1)%fftstruc,pot_gspace,den(1,is),&
        rden(1,is),1)
   if (iand(pw_params%output(1), 8) == 8) call get_timing(time(1),t3)
   !
   !     symmetrize charge
   !
   !     symmetrize only if pot_gspace%istar true
   !     because if  pot_gspace%istar false the arrays are not defined
   !     Francesco Mauri
   !
   if (pw_params%isymmetrize) then  
      call symmetrize_scalar_global(den(1,is), pot_gspace)  
   end if
  end do
  !
  !     ------------- compute new electronic potential vion --------
  !
  call velect(0, pot_gspace, crys, ham(1,1)%fftstruc, energs, &
       pw_params, denc, vion, den, rden(1,1))
  !
  !     ------------- add new vhxc to vion -------------------------
  !
  vin=vion

  do is = 1, crys%nspin  
      vion(:, is) = vion(:, is) + vbareion(:,is)  
      if (pot_gspace%ig0 > 0) vion(pot_gspace%ig0,is) = zzero
      !
      !     ------------ setup new local potential for Hamiltonian -------
      !
      call setup_local_potential(1, ham(1,1)%fftstruc, pot_gspace, &
         vion(1,is), ham(1,is)%vloc(1))
  end do  


  if (iand(pw_params%output(1), 8) == 8) call get_timing(time(5),t3)
  return  
  !
end subroutine update_hamiltonian


