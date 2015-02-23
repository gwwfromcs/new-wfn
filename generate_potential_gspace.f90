!     @process extchk
!
subroutine generate_potential_gspace(pot_gspace, crys, syms, &
     pw_params, myproc, nproc)
  !
  include 'use.h'  
  implicit none              ! implicit? Just say no!
  include 'interface.h'  
  !
  !     INPUT:
  !     -----
  !
  type(pw_parameter), intent(in) :: pw_params    ! plane wave parameters
  type(crystal), intent(in) :: crys              ! for the lattice vectors
  type(symmetry), intent(in) :: syms  
  integer, intent(in) :: nproc, myproc           ! parallel info
  !
  !     OUTPUT:
  !     ------
  !
  type(parallel_gspace), intent(out) :: pot_gspace    ! potential gspace
  !
  !     DESCRIPTION:
  !     -----------
  !
  !     generates a large potential gspace with all the features:
  !     stars, gvectors, kinetic energy, structure factors
  !
  !
  !     ------------------------------------------------------------------
  !
  real(dp) :: emax    ! energy cutoff for the wave functions!
  integer nfn

  emax = pw_params%emax  
  pot_gspace%gmax = dtwo * sqrt(emax)          ! 4*emax
  pot_gspace%rk = (/ dzero, dzero, dzero /)    ! no shift
  pot_gspace%nproc = nproc                     ! scatter across all processors
  pot_gspace%name = 'potential'  
  pot_gspace%myproc = myproc
  pot_gspace%fftsize = (/ 0, 0, 0, 0, 0 /)
  pot_gspace%istar = .true.
  pot_gspace%igvec = .true.  
  pot_gspace%iekin = .true.
  pot_gspace%iinverse = .false.  
  pot_gspace%istruc = .true.  
  pot_gspace%ipackinfo = .true.  
  !
  !     F. Mauri istar information if there are symmetries
  !
  pot_gspace%istar = .not. (syms%ntrans == 1)  

  call generate_gspace(pw_params%output(1), pot_gspace, crys, syms)  
  write(9, 10) pot_gspace%totlength, pot_gspace%fftsize(1:3)  
  call myflush(9)  

10 format(' SIZE OF POTENTIAL GSPACE = ', i8, &
       &      ' WITH FFT GRID ', 3i6)

end subroutine generate_potential_gspace
