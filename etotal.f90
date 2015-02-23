!

subroutine etotal(crys, iter, iconv, dvmaxmax, enew, eold, &
     pw_params, gs, vion, vin, vout, den,Etot_corr,ELSDA_corr,Etot_corr_old,Eband_corr)
  !
  !     checks for convergence, and finds the total energy.
  !     adapted from sverre froyen plane wave program
  !     written june 8 1987.jlm
  !
  !     parallel fortran 90 version 1996 by Bernd Pfrommer, UCB
  !
  use all_to_all_module  
  include 'use.h'  
  implicit none  
  include 'interface.h'  
  include 'all_to_all.h'  
  !
  !     INPUT:
  !     ------
  type(crystal), intent(in) :: crys         ! crystal parameters: vcell, ztot
  type(parallel_gspace), intent(in) :: &
       gs                              ! the gspace for potential and density
  type(pw_parameter), intent(in) :: pw_params  
  integer, intent(in) :: iter                              ! iteration number
  complex(dp), intent(in) :: &
       den(gs%length, crys%nspin), &       ! symmetrized charge density in gs
       vion (gs%length, crys%nspin), &   ! ionic potential, screened with vin 
       vin (gs%length, crys%nspin), & ! inp vhxc (with which wavefn are calc)
       vout (gs%length, crys%nspin)      ! output vhxc (calc from new wavefn)
  !
  !     OUTPUT:
  !     -------
  !
  logical, intent(out) :: iconv       ! if true then convergence was achieved
  real(dp), intent(out) :: dvmaxmax    ! max potential deviation of all procs
  real(dp)::Etot_corr,ELSDA_corr,Etot_corr_old,Eband_corr
  !
  !     INPUT/OUTPUT:
  !     -------------
  !
  type(energy), intent(inout) :: enew, &
       eold                                              ! store old energies
  !
  !     ------------------------------------------------------------------
  !                            local variables
  !
  real(dp), parameter :: &
       eps = 1.0d-6, &                    ! maximum loss of charge acceptable
       epscmplx =1.0d-7             ! maximum complex energy to be tolerated

  real(dp) :: t0, dvmax, &
       spindiff, &                          ! differences of spinup/spin down
       diff,&                                   ! deviation of n(G=0) from ztot
       proc,ekinmax
  real(dp), external :: gimmetime  
  type(energy) :: de                                     ! energy differences
  complex(dp) :: dv, e2, e4, cdtot,vinmax,voutmax  
  integer :: nnotconv, &         ! number of unconverged potential components
       i, j, k, idvmax, is, iproc(gs%nproc), iproc_temp  
  logical :: isp  
  logical, SAVE :: iconv1    ! if energy convergence true for last SCF step
  !
  !     ------------------------------------------------------
  !
  t0 = gimmetime()  
  isp = .false.  
  if (crys%nspin == 2) isp = .true.  
  if (iter .eq. 1) iconv1 = .false.
  !
  !     ------------------ check for convergence -----------------------
  !
  nnotconv = 0  
  dvmax = -done  
  idvmax = -1 
  do is = 1, crys%nspin  
     do i = 1, gs%length  
        dv = vout(i, is) - vin(i, is)  
        if (abs(dv) > dvmax) then  
           if (gs%ekin(i) > dzero) then  
              dvmax = abs(dv)  
              idvmax = i  
           end if
        end if
        if (abs(dv) > pw_params%epscv) then
           if (gs%ekin(i) > dzero) then  
              nnotconv = nnotconv + 1  
           end if
        end if
     end do
  end do
  call all_sum_all(nnotconv)  
  if (nnotconv > 0) then  
     iconv = .false.  
  else  
     iconv = .true.  
  end if
  
  dvmaxmax = dvmax  
  ekinmax=gs%ekin(idvmax)
  vinmax=vin(idvmax,1)
  voutmax=vout(idvmax,1)  
  call all_max_all(dvmaxmax)

  iproc=0
  if (dvmaxmax .eq. dvmax) iproc(gs%myproc+1)=gs%myproc+1
  do i=1,gs%nproc
    iproc_temp=iproc(i)
    call all_sum_all(iproc_temp) 
    if (iproc_temp > 0) then
      goto 10
    end if
  end do

10 continue

  call my_broadcast(ekinmax,iproc_temp-1)
  call my_broadcast(vinmax,1,iproc_temp-1)
  call my_broadcast(voutmax,1,iproc_temp-1)

  if (gs%myproc .eq. 0) then
    write(9, 210) nnotconv
    write(9, 220) dvmaxmax 
    write(9, 240) ekinmax  
    call myflush(9)  
  end if
  !
  !     -------------------- initialize energies -------------------------
  !
  enew%hxc = dzero  
  enew%ion = dzero  
  enew%hartree = dzero  
  !
  !     start loop over gspace
  !
  e2 = zzero
  e4 = zzero
  spindiff = dzero  
  diff = dzero
  if (isp) then  
     do i = 1, gs%length  
        !
        !        compute the hxc correction from the potential of the
        !        valence electrons. this will be subtracted from the sum
        !        of the eigenvalues.
        !        then compute the contribution to the total energy from
        !        the interaction with the ions and the electrostatic
        !        (hartree) interaction between the electrons.
        !
        if (gs%ekin(i) > dzero) then  
           e2 = e2 + vin(i, 1) * conjg(den(i, 1)) + &
                vin(i, 2) * conjg(den(i, 2))
           e4 = e4 + (vion(i, 1) - vin(i, 1)) * conjg(den(i, 1)) + &
                (vion(i, 2) - vin(i, 2)) * conjg(den(i, 2))
           cdtot = den(i, 1) + den(i, 2)  
           enew%hartree = enew%hartree + real(cdtot * conjg(cdtot), dp) / &
                gs%ekin(i)
        else  
           diff = real(den(i, 1) + den(i, 2), dp) - crys%ztot  
           !               write(9,*) den(i,1),den(i,2)
           spindiff = real(den(i, 1) - den(i, 2), dp)  
        end if
     end do
  else  
     do i = 1, gs%length  
        if (gs%ekin(i) > dzero) then  
           e2 = e2 + vin(i, 1) * conjg(den(i, 1))  
           e4 = e4 + (vion(i, 1) - vin(i, 1)) * conjg(den(i, 1))  
           enew%hartree = enew%hartree + &
                real(den(i, 1) * conjg(den(i, 1)), dp) / gs%ekin(i)
        else  
           diff = real(den(i, 1)) - crys%ztot  
        end if
     end do
  end if
  call all_sum_all(diff)  
  call all_sum_all(enew%hartree)  
  call all_sum_all(e2)  
  call all_sum_all(e4)  
  call all_sum_all(spindiff)            ! propagate spindiff to all processors
  if (aimag(e2) > epscmplx) then  
     write(9, *) 'etotal: complex energy in e2 encountered:', e2  
     call mystop  
  end if
  if (aimag(e4) > epscmplx) then  
     write(9, *) 'etotal: complex energy in e4 encountered:', e4  
     call mystop  
  end if
  enew%hxc = real(e2, dp)  
  enew%ion = real(e4, dp)  
  enew%hartree = pi4 * enew%hartree / crys%vcell  
  !
  !     find the nonlocal energy.
  !
  enew%nloc = enew%eband - enew%hxc - enew%ektot - enew%ion  
  !
  !     the total energy.
  !
  !     the magnetic contribution is subtracted
  !     the smearing energy (-TS term) is added
  !
  enew%total = enew%eband-enew%hxc + enew%hartree+enew%xc + &
       enew%alpha + enew%ewald + enew%esmear - enew%magnetic+Etot_corr

  enew%totnonv = enew%total - dhalf * enew%esmear
  !
  !     find differences from last iteration.
  !
  de = enew - eold  
  eold = enew  
  if (iter == 1) de = dzero
  !
  ! check for energy convergence
  if(iter .gt. 1) then
    if( abs(de%total) < pw_params%epsce) then
      if(iconv1) then
       iconv = .true.
      else
        iconv1 = .true.
      end if
    else
      iconv1 = .false.
    end if
  end if     
  !
  !     printout
  !
  write(9, 102) iter, enew%eband, de%eband, enew%hxc, de%hxc, &
       enew%ektot, de%ektot, enew%ion, de%ion, enew%nloc, de%nloc, &
       enew%hartree, de%hartree, enew%xc, de%xc, enew%alpha, de%alpha, &
       enew%ewald, de%ewald, enew%esmear, de%esmear, -enew%magnetic, &
       -de%magnetic, Etot_corr,Etot_corr-Etot_corr_old, &
        Eband_corr,0.d0, enew%total+ELSDA_corr,0.d0,&
       enew%totnonv, de%totnonv, enew%total, de%total

  if (abs(diff) > eps) call warn(290, diff, i, 'etotal')  
  !
  if (isp) then  
     write(9, 200) spindiff  
     write(9, 205) enew%afmagmom  
  end if

  if (iand(pw_params%output(1), 8) == 8) write(9, 940) gimmetime() - t0
  call myflush(9)  

  return  

102 format(///,11X,'ITERATION NUMBER',I3,' :',/,11X,21('-'),/ &
       &     /,11X,'BAND ENERGY          =',2X,F17.10,2X,F14.10, &
       &     /,11X,'HXC  CORRECTION      =',2X,F17.10,2X,F14.10, &
       &     /,11X,57('-'), &
       &     /,11X,'KINETIC  ENERGY      =',2X,F17.10,2X,F14.10, &
       &     /,11X,'IONIC    ENERGY      =',2X,F17.10,2X,F14.10, &
       &     /,11X,'NONLOCAL ENERGY      =',2X,F17.10,2X,F14.10, &
       &     /,11X,'HARTREE  ENERGY      =',2X,F17.10,2X,F14.10, &
       &     /,11X,'EXCHANGE CORRELATION =',2X,F17.10,2X,F14.10, &
       &     /,11X,57('-'), &
       &     /,11X,'ALPHA TERM           =',2X,F17.10,2X,F14.10, &
       &     /,11X,'EWALD TERM           =',2X,F17.10,2X,F14.10, &
       &     /,11X,'SMEARING CORRECTION  =',2X,F17.10,2X,F14.10, &
       &     /,11X,'MAGNETIC CORRECTION  =',2X,F17.10,2X,F14.10, &
       &     /,11X,57('-'), &
       &     /,11X,'LSDA+U CORRECTION    =',2X,F17.10,2X,F14.10, &
       &     /,11X,'LSDA+U BAND CORRECT. =',2X,F17.10,2X,F14.10, &
       &     /,11X,57('-'), & 
       &     /,11X,'TOTAL ENERGY(LSDA)   =',2X,F17.10,2X,F14.10, &
       &     /,11X,'NON-VARIATIONAL ETOT =',2X,F17.10,2X,F14.10, &
       &     /,11X,'TOTAL ENERGY(LSDA+U) =',2X,F17.10,2X,F14.10)

103 format(/,' TOTAL ENERGY =',F17.10)  
200 format(/,11x,'TOTAL MAGNETIC MOMENT [MU_BOHR/UCELL]=',f12.8)  

205 format(11x,'TOTAL ANTIFERROMAGNETIC MOMENT [MB/U]=',f12.8)  
210 format(/' NUMBER OF UNCONVERGED POTENTIAL COMPONENTS:',i10)  
220 format(' LARGEST POTENTIAL DIFFERENCE: ',g12.6)
240 format(' WITH KINETIC ENERGY:', f12.6)  
940 format(' TIME FOR ETOTAL:',f12.3)  

end subroutine etotal

!
subroutine etotal_dir(ipr,crys, iter, iconv, dvmaxmax, enew, eold, &
     pw_params, gs, vion, vin, vout, den)
  !
  !     checks for convergence, and finds the total energy.
  !     adapted from sverre froyen plane wave program
  !     written june 8 1987.jlm
  !
  !     parallel fortran 90 version 1996 by Bernd Pfrommer, UCB
  !
  use all_to_all_module  
  include 'use.h'  
  implicit none  
  include 'interface.h'  
  include 'all_to_all.h'  
  !
  !     INPUT:
  !     ------
  integer ipr
  type(crystal), intent(in) :: crys         ! crystal parameters: vcell, ztot
  type(parallel_gspace), intent(in) :: &
       gs                              ! the gspace for potential and density
  type(pw_parameter), intent(in) :: pw_params  
  integer, intent(in) :: iter                              ! iteration number
  complex(dp), intent(in) :: &
       den(gs%length, crys%nspin), &       ! symmetrized charge density in gs
       vion (gs%length, crys%nspin), &   ! ionic potential, screened with vin 
       vin (gs%length, crys%nspin), & ! inp vhxc (with which wavefn are calc)
       vout (gs%length, crys%nspin)      ! output vhxc (calc from new wavefn)
  !
  !     OUTPUT:
  !     -------
  !
  logical, intent(out) :: iconv       ! if true then convergence was achieved
  real(dp), intent(out) :: dvmaxmax    ! max potential deviation of all procs
  !
  !     INPUT/OUTPUT:
  !     -------------
  !
  type(energy), intent(inout) :: enew, &
       eold                                              ! store old energies
  !
  !     ------------------------------------------------------------------
  !                            local variables
  !
  real(dp), parameter :: &
       eps = 1.0d-6, &                    ! maximum loss of charge acceptable
       epscmplx = 1d-4 !1.0d-7             ! maximum complex energy to be tolerated

  real(dp) :: t0, dvmax, &
       spindiff, &                          ! differences of spinup/spin down
       diff,&                                   ! deviation of n(G=0) from ztot
       proc,ekinmax
  real(dp), external :: gimmetime  
  type(energy) :: de                                     ! energy differences
  complex(dp) :: dv, e2, e4, cdtot,vinmax,voutmax  
  integer :: nnotconv, &         ! number of unconverged potential components
       i, j, k, idvmax, is, iproc(gs%nproc), iproc_temp
  logical :: isp  
  !
  !     ------------------------------------------------------
  !
  t0 = gimmetime()  
  isp = .false.  
  if (crys%nspin == 2) isp = .true.  
  !
  !     ------------------ check for convergence -----------------------
  !
  nnotconv = 0  
  dvmax = dzero  
  idvmax = -1  
  do is = 1, crys%nspin  
     do i = 1, gs%length  
        dv = vout(i, is) - vin(i, is)  
        if (abs(dv) > dvmax) then  
           if (gs%ekin(i) > dzero) then  
              dvmax = abs(dv)  
              idvmax = i  
           end if
        end if
        
!        if ((abs(real(dv, dp)) > pw_params%epscv) .or. &
!             (abs(aimag(dv)) > pw_params%epscv)) then
        if (abs(dv) > pw_params%epscv) then
           if (gs%ekin(i) > dzero) then  
              nnotconv = nnotconv + 1  
           end if
        end if
     end do
  end do
  call all_sum_all(nnotconv)  
  if (nnotconv > 0) then  
     iconv = .false.  
  else  
     iconv = .true.  
  end if
  dvmaxmax = dvmax  
  ekinmax=gs%ekin(idvmax)
  vinmax=vin(idvmax,1)
  voutmax=vout(idvmax,1) 
  call all_max_all(dvmaxmax)

  iproc=0
  if (dvmaxmax .eq. dvmax) iproc(gs%myproc+1)=gs%myproc+1
  do i=1,gs%nproc
    iproc_temp=iproc(i)
    call all_sum_all(iproc_temp) 
    if (iproc_temp > 0) then
      goto 10
    end if
  end do

10 continue

  call my_broadcast(ekinmax,iproc_temp-1)
  call my_broadcast(vinmax,1,iproc_temp-1)
  call my_broadcast(voutmax,1,iproc_temp-1)

  if (ipr > 1) then
    if (gs%myproc .eq. 0) then
      write(9, 210) nnotconv 
      write(9, 220) dvmaxmax   
      write(9, 240) ekinmax 
      call myflush(9)  
  !!-------------------------------
  !! George's profiling
  !!-------------------------------
  !call MPI_Pcontrol(2)
  !!-------------------------------
  !! George's profiling
  !!-------------------------------
    end if
  end if
  !
  !     -------------------- initialize energies -------------------------
  !
  enew%hxc = dzero  
  enew%ion = dzero  
  enew%hartree = dzero  
  !
  !     start loop over gspace
  !
  e2 = zzero
  e4 = zzero
  spindiff = dzero  
  diff = dzero
  if (isp) then  
     do i = 1, gs%length  
        !
        !        compute the hxc correction from the potential of the
        !        valence electrons. this will be subtracted from the sum
        !        of the eigenvalues.
        !        then compute the contribution to the total energy from
        !        the interaction with the ions and the electrostatic
        !        (hartree) interaction between the electrons.
        !
        if (gs%ekin(i) > dzero) then  
           e2 = e2 + vin(i, 1) * conjg(den(i, 1)) + &
                vin(i, 2) * conjg(den(i, 2))
           e4 = e4 + (vion(i, 1) - vin(i, 1)) * conjg(den(i, 1)) + &
                (vion(i, 2) - vin(i, 2)) * conjg(den(i, 2))
           cdtot = den(i, 1) + den(i, 2)  
           enew%hartree = enew%hartree + real(cdtot * conjg(cdtot), dp) / &
                gs%ekin(i)
        else  
           diff = real(den(i, 1) + den(i, 2), dp) - crys%ztot  
           !               write(9,*) den(i,1),den(i,2)
           spindiff = real(den(i, 1) - den(i, 2), dp)  
        end if
     end do
  else  
     do i = 1, gs%length  
        if (gs%ekin(i) > dzero) then  
           e2 = e2 + vin(i, 1) * conjg(den(i, 1))  
           e4 = e4 + (vion(i, 1) - vin(i, 1)) * conjg(den(i, 1))  
           enew%hartree = enew%hartree + &
                real(den(i, 1) * conjg(den(i, 1)), dp) / gs%ekin(i)
        else  
           diff = real(den(i, 1)) - crys%ztot  
        end if
     end do
  end if
  call all_sum_all(diff)  
  call all_sum_all(enew%hartree)  
  call all_sum_all(e2)  
  call all_sum_all(e4)  
  call all_sum_all(spindiff)            ! propagate spindiff to all processors
  if (aimag(e2) > epscmplx) then  
     write(9, *) 'etotal: complex energy in e2 encountered:', e2  
     call mystop  
  end if
  if (aimag(e4) > epscmplx) then  
     write(9, *) 'etotal: complex energy in e4 encountered:', e4  
     call mystop  
  end if
  enew%hxc = real(e2, dp)  
  enew%ion = real(e4, dp)  
  enew%hartree = pi4 * enew%hartree / crys%vcell  
  !
  !     find the nonlocal energy.
  !
  enew%nloc = enew%eband - enew%hxc - enew%ektot - enew%ion  
  !
  !     the total energy.
  !
  !     the magnetic contribution is subtracted
  !     the smearing energy (-TS term) is added
  !
  enew%total = enew%eband-enew%hxc + enew%hartree+enew%xc + &
       enew%alpha + enew%ewald + enew%esmear - enew%magnetic

  enew%totnonv = enew%total - dhalf * enew%esmear  
  !
  !     find differences from last iteration.
  !
  de = enew - eold  
!  eold = enew  
  if (iter == 1) de = dzero
  !
  !     printout
  !
  if (ipr .eq. 0) then

  else if(ipr .eq. 1) then
    write(9, 103)  enew%total,iter
  else
    write(9, 102) iter, enew%eband, de%eband, enew%hxc, de%hxc, &
       enew%ektot, de%ektot, enew%ion, de%ion, enew%nloc, de%nloc, &
       enew%hartree, de%hartree, enew%xc, de%xc, enew%alpha, de%alpha, &
       enew%ewald, de%ewald, enew%esmear, de%esmear, -enew%magnetic, &
       -de%magnetic, enew%totnonv, de%totnonv, enew%total, de%total
  end if

  if (abs(diff) > eps) call warn(290, diff, i, 'etotal')  
  !
  if (isp) then  
     write(9, 200) spindiff  
     write(9, 205) enew%afmagmom  
  end if

  if (iand(pw_params%output(1), 8) == 8) write(9, 940) gimmetime() - t0
  call myflush(9)  

  !!-------------------------------
  !! George's profiling
  !!-------------------------------
  !call MPI_Pcontrol(2)
  !!-------------------------------
  !! George's profiling
  !!-------------------------------
  return  

102 format(///,11X,'ITERATION NUMBER',I3,' :',/,11X,21('-'),/ &
       &     /,11X,'BAND ENERGY          =',2X,F17.10,2X,F14.10, &
       &     /,11X,'HXC  CORRECTION      =',2X,F17.10,2X,F14.10, &
       &     /,11X,57('-'), &
       &     /,11X,'KINETIC  ENERGY      =',2X,F17.10,2X,F14.10, &
       &     /,11X,'IONIC    ENERGY      =',2X,F17.10,2X,F14.10, &
       &     /,11X,'NONLOCAL ENERGY      =',2X,F17.10,2X,F14.10, &
       &     /,11X,'HARTREE  ENERGY      =',2X,F17.10,2X,F14.10, &
       &     /,11X,'EXCHANGE CORRELATION =',2X,F17.10,2X,F14.10, &
       &     /,11X,57('-'), &
       &     /,11X,'ALPHA TERM           =',2X,F17.10,2X,F14.10, &
       &     /,11X,'EWALD TERM           =',2X,F17.10,2X,F14.10, &
       &     /,11X,'SMEARING CORRECTION  =',2X,F17.10,2X,F14.10, &
       &     /,11X,'MAGNETIC CORRECTION  =',2X,F17.10,2X,F14.10, &
       &     /,11X,57('-'), &
       &     /,11X,'NON-VARIATIONAL ETOT =',2X,F17.10,2X,F14.10, &
       &     /,11X,'TOTAL ENERGY         =',2X,F17.10,2X,F14.10)

103 format(/,' TOTAL ENERGY =',F17.10,I3)  
200 format(/,11x,'TOTAL MAGNETIC MOMENT [MU_BOHR/UCELL]=',f12.8)  

205 format(11x,'TOTAL ANTIFERROMAGNETIC MOMENT [MB/U]=',f12.8)  
210 format(/' NUMBER OF UNCONVERGED POTENTIAL COMPONENTS:',i10)  
220 format(' LARGEST POTENTIAL DIFFERENCE: ',g12.6)
240 format(' WITH KINETIC ENERGY:', f12.6)  
940 format(' TIME FOR ETOTAL:',f12.3)  

end subroutine etotal_dir
