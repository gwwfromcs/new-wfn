!-*-Fortran-*-
!
subroutine update_hamiltonian_metal(ipr,kp, ham, ldawfn, len, neig, wfn, &
     pot_gspace, bands, crys, energs, pw_params, denc, rden, den, vion, &
     vbareion,iter,vin,update_pot,oldenergy,time,t3)
  !
  !     2001 David Raczkowski based on update_hamiltonian
  !
  use all_to_all_module  
  include 'use.h'  
  implicit none                 ! never uncomment this line.
  include 'interface.h'  
  include 'all_to_all.h'  
  include 'flibcalls.ph'
  !
  !     INPUT:
  !     -----
  !
  logical update_pot,ipr
  type(kpoint), intent(in) :: kp  
  integer :: &
       len(kp%nrk), &           ! length of the different gspaces
       neig, &                  ! number of occupied wave functions
       ldawfn, &                   ! leading dimension of wavefn array
       iter
  type(parallel_gspace), intent(in) :: &
       pot_gspace               ! gspace for charge density and potential
  type(crystal), intent(in) :: crys  
  real(dp) oldenergy
  type(pw_parameter), intent(in) :: pw_params  
  type(band), intent(inout) :: bands  
  complex(dp), intent(in) :: &
       vbareion(pot_gspace%length,crys%nspin ), &  ! bare ionic potential
       denc(pot_gspace%length), &  ! the core charge density
       wfn(ldawfn, kp%nrk,crys%nspin)      ! the wave functions in gspace
  !
  !     INPUT/OUTPUT:
  !     -------------
  !
  type(hamiltonian), intent(inout) :: &
       ham(kp%nrk,crys%nspin )              ! the Hamiltonian to be updated
  complex(dp), intent(inout):: &
       vin(pot_gspace%length,crys%nspin)     ! old potenital
  type(energy), intent(inout) :: energs 
  real(dp), intent(inout) :: time(num_time),t3
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
     rden(pot_gspace%r_size,crys%nspin ), & ! realspace work array rden(r_size)
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
  real(dp), external :: myderfc,nmerfc 
  complex(dp) :: cdtot  
  integer :: n, i, ntot, ids, irk, iscfconv, ispn, is,lm,lmmax,j, &
          iwfn,rlen, nfn   
  real(dp) errscfconv,FI, eesqh,sq2h, piesqq, zeta, ee, X 
  type(energy), SAVE :: oldenergs
  real(dp)  spinf,delta,dsmear,sqrpi,lambda,num,lambda2
  complex(dp)  &
       vout(pot_gspace%length,crys%nspin)     ! old potenital
    real(dp) bands_old(neig,kp%nrk,crys%nspin),&
     bands_del(neig,kp%nrk,crys%nspin)   
  real(dp) fermi_tmp(2)

  ee=exp(done) 
  eesqh=sqrt(ee)*0.5d0
  sq2h=sqrt(2.0)*0.5d0
  piesqq=sqrt(ee*pi)*0.25d0

  dsmear = pw_params%smearing / ryd  
  ispn=crys%nspin
  sqrpi = sqrt(pi)
  if (bands%nspin >= 2) then
     spinf = dhalf
  else
     spinf = done 
  end if

  if (iter .eq. 1) then 
    oldenergs=dzero
  end if
  !
  !     figure out occupation numbers
  !
  !
  !     ------------- compute new charge -------------------------
  !
  ntot = ham(1,1)%gspace%fftsize(1) * ham(1,1)%gspace%fftsize(2) * &
       ham(1,1)%gspace%fftsize(3)
  
  rden(1:ham(1,1)%fftstruc%r_size,1:crys%nspin) = zzero
  cdtot = zzero

  rlen = ham(1,1)%fftstruc%r_size   ! total number of points in realspace block

  energs%eband = dzero
  energs%ektot = dzero
  num=dzero
  energs%esmear = dzero

  do is=1,ispn
    do irk = 1, kp%nrk  
      call kinet(ham(irk,is), wfn(1, irk,is), bands%nband(irk,is), &
             bands%ekn(1, irk, is))
      do n = 1, neig,pw_params%nbandsfft         ! loop over occupied bands
        nfn=min(pw_params%nbandsfft, neig + 1 - n)
        energs%ektot = energs%ektot + bands%occup(n, irk, is) * &
                                         bands%ekn(n, irk, is)
        energs%eband = energs%eband + &
                bands%occup(n, irk, is) * bands%energy(n, irk, is)

        FI=bands%occup(n, irk, is)/(kp%w(irk)*dtwo)
        if(pw_params%smearing_method .eq. 2) then
          if(abs(FI) .gt. 1.d-06 .and. abs(FI-done) .gt. 1.d-06) then
            energs%esmear=energs%esmear+dtwo*dsmear* kp%w(irk) * &
                   (FI*log(FI)+(done-FI)*log(done-FI))
          end if
         else
          if(abs(FI) .gt. 1.d-10 .and. abs(FI-done) .gt. 1.d-10) then
            if (FI .lt. 0.5) then
              X=sqrt(log(eesqh/FI))-sq2h
            else
              X=sq2h-sqrt(log(eesqh/(done-FI)))
            end if
            X=abs(X)
            zeta=eesqh*X*exp(-(X+sq2h)**2)+piesqq * myderfc(X+sq2h)
            energs%esmear=energs%esmear-dtwo*dsmear*kp%w(irk) *zeta
          end if
        end if
        !
        !           fourier transform to realspace
        !
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(5),t3)
        call fourier_transform(-1, ham(1,1)%fftstruc, ham(irk,is)%gspace, &
        wfn(1 + (n - 1) * len(irk), irk,is), ham(1,1)%fftstruc%rspacebuf(1), 1)
        if (iand(pw_params%output(1), 8) == 8) call get_timing(time(1),t3)
        !
        !           square and add to the charge density
        !
        do iwfn=0,nfn-1
          do i = 1, ham(1,1)%fftstruc%r_size  
            !               cdtot = cdtot+ ham(1)%fftstruc%rspacebuf(i)
            !     $              *conjg(ham(1)%fftstruc%rspacebuf(i))
            !     $              *bands%occup(n,irk,1)
            rden(i,is) = rden(i,is) + ham(1,1)%fftstruc%rspacebuf(i+iwfn*rlen)&
                  * conjg(ham(1,1)%fftstruc%rspacebuf(i+iwfn*rlen)) &
               *bands%occup(n+iwfn, irk, is)*spinf
          end do
        end do
        !            call all_sum_all(cdtot)
        !            write(9,*) irk,n,cdtot/dble(ntot)
     end do  !neig
  end do    !nrk
  !
  !     fourier transform back to momentum space
  !
  if (iand(pw_params%output(1), 8) == 8) call get_timing(time(5),t3)
  call fourier_transform(1, ham(1,1)%fftstruc, pot_gspace, den(1,is),&
                                  rden(1,is), 1)
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
  end do !is
  energs%eband = energs%eband / real(bands%nspin, dp)
  energs%ektot = energs%ektot / real(bands%nspin, dp)
  energs%esmear = energs%esmear / real(bands%nspin, dp)

  !
  !     ------------- compute new electronic potential vion --------
  !
  call velect(0, pot_gspace, crys, ham(1,1)%fftstruc, energs, &
       pw_params, denc, vion, den, rden(1,1))
  !
  !     ------------- add new vhxc to vion -------------------------
  !
  do is = 1, crys%nspin  
     vout(:, is) = vion(:, is) + vbareion(:,is)  
  end do 

  call etotal_dir(ipr,crys, iter, iscfconv, errscfconv, energs, oldenergs, &
          pw_params, pot_gspace, vout(1,1), vin(1,1), vion(1,1), den(1,1))

   if(update_pot)  oldenergs=energs

 10   continue

  if (update_pot) then ! .and. energs%total .lt. oldenergy) then
  
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

  else
    vion=vin
  end if

  if (iand(pw_params%output(1), 8) == 8) call get_timing(time(5),t3)

  return  
  !
end subroutine update_hamiltonian_metal

