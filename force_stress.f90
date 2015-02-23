!     @process extchk
!
subroutine force_stress(crys, syms, bands, kpoints, pspot, &
     energs, pw_params, pot_gspace, k_gspace, vion, vin, vout, den, &
     denc, dvql, ddc, vql, dnc, wavefn, ewaldf, ssum, fsum, ffts)
  !
  !     calls the right subroutines to compute force and stress
  !
  !     1996 by Bernd Pfrommer, while at UCB
  !
  include 'use.h'  
  implicit none                    ! implicit? Just say no!
  include 'interface.h'
  !
  !     INPUT:
  !     -----
  !

  type(fft_struc), intent(in) :: ffts               ! information for the FFTs
  type(symmetry), intent(in) :: syms  
  type(crystal), intent(in) :: crys  
  type(band), intent(in) :: bands  
  type(kpoint), intent(in) :: kpoints  
  type(pseudo_potential), intent(in) :: pspot  
  type(energy), intent(in) :: energs  
  type(force), intent(in) :: ewaldf  
  type(pw_parameter), intent(in) :: pw_params  
  type(parallel_gspace), intent(in) :: &
       pot_gspace, &                                ! gspace for the potential
       k_gspace(*)                         ! gspaces for the different kpoints
  complex(dp), intent(in) :: &
       vion(pot_gspace%length, crys%nspin), &   ! the screened ionic potential
       vin(pot_gspace%length, crys%nspin), & ! electronic V_hxc input potent'l
       vout(pot_gspace%length, crys%nspin), &        ! electronic V_hxc output
       den(pot_gspace%length, crys%nspin), &          ! valence charge density
       denc(pot_gspace%length), &                        ! core charge density
       dvql(pot_gspace%length), &    ! derivative of the local pseudopotential
       ddc(pot_gspace%length)              ! derivative of core charge density
  real(dp), intent(in) :: &
       vql(pot_gspace%length, crys%ntype), &  ! local pspot of different atoms
       dnc(pot_gspace%length, crys%ntype)       ! core charge density of atoms
  type(complex_gspace_array), intent(in) :: wavefn
  !  
  !     INPUT/OUTPUT:
  !     ------------
  !
  real(dp), intent(inout) :: ssum(6), &                     ! sigma * V  [Ryd]
       fsum(3, crys%mxdatm, crys%ntype)                                ! dE/dq
  !
  !     ------------------ local variables ------------------------
  !
  real(dp) :: aadot (6)                                     ! realspace metric
  !
  !     ------------------------------------------------------------------
  !
  fsum = dzero
  ssum = dzero

  !  local part of force and stress

  call forstressloc(1, crys%nspin, fsum, ssum, aadot, energs%xc, &
       energs%alpha, pot_gspace%length, pot_gspace%gvec(1, 1), &
       pot_gspace%ekin(1), vion(1, 1), vin(1, 1), vout(1, 1), den(1, 1), &
       denc(1), dvql(1), ddc(1), vql(1, 1), dnc(1, 1), crys%ntype, &
       crys%natom, crys%rat, crys%bdot, crys%vcell, crys%mxdatm)

  if (pw_params%NL_rspace(2)) then
     call forstressnloc_rsp(1, crys, bands, k_gspace, pw_params%ekinmod, &
       pspot, kpoints%rk, aadot, wavefn,fsum, ssum,pot_gspace,ffts,pw_params) 
  else
    call forstressnloc(1, crys, bands, k_gspace, pw_params%ekinmod, &
       pspot, kpoints%rk, aadot, wavefn, fsum, ssum)
  end if

  if (pw_params%icorr == 'pw') then  
     call stresspw91(ffts, crys, pot_gspace, den, denc, ssum)  
  end if

  if (pw_params%icorr == 'pb') then  
     call stresspbe(ffts, crys, pot_gspace, den, denc, ssum)  
  end if

  call forstresssym(1, 1, fsum, ssum, pw_params%icorr, ewaldf%force, &
       ewaldf%stress, crys%ntype, crys%natom, crys%rat, crys%avec, &
       crys%bvec, crys%bdot, crys%vcell, syms%ntrans, syms%mtrx, &
       syms%tnp, syms%rsymmat, crys%mxdatm)

  return  

end subroutine force_stress
