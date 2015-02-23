!     @process extchk
!
subroutine generate_k_gspaces(k_gspace, maxlength, pot_gspace, &
     kpoints, crys, syms, emax, ioutput)
  !
  !
  !
  use all_to_all_module
  include 'use.h'  
  implicit none               ! implicit? Just say no!
  include 'interface.h'  
  include 'all_to_all.h'  
  !
  !     INPUT:
  !     -----
  !
  integer, intent(in) :: ioutput        ! output flag
  type(crystal), intent(in) :: crys     ! for the lattice vectors
  type(symmetry), intent(in) :: syms  
  type(kpoint), intent(in) :: kpoints  
  real(dp), intent(in) :: emax          ! energy cutoff for the wave functions!
  type(parallel_gspace), intent(in) :: pot_gspace  
  !
  !     OUTPUT:
  !     ------
  !
  type(parallel_gspace), intent(out) :: k_gspace(*)                ! gspace
  integer, intent(out) :: maxlength   ! maximum length of all gspaces generated
  integer :: maxlength_all_k   ! maximum length of all gspaces generated
  !
  !     DESCRIPTION:
  !     -----------
  !
  !     generates the small gspaces for the various kpoints
  !
  !
  !     --------  local variables ---------------------------------
  !
  integer :: irk, mtrxdim, i  

  write(9, 10)  
  maxlength = 0  
  maxlength_all_k=0

  do irk = 1, kpoints%nrk  
     k_gspace(irk)%gmax = sqrt(emax)  
     k_gspace(irk)%rk = kpoints%rk(:, irk)   
     k_gspace(irk)%nproc = pot_gspace%nproc
     k_gspace(irk)%myproc = pot_gspace%myproc
     k_gspace(irk)%name = 'k_gspace'  
     k_gspace(irk)%fftsize = pot_gspace%fftsize
     k_gspace(irk)%istar = .false.
     k_gspace(irk)%igvec = .true.
     k_gspace(irk)%iekin = .true.  
     k_gspace(irk)%iinverse = .false.  
     k_gspace(irk)%istruc = .false.  
     k_gspace(irk)%ipackinfo = .true.  
     call generate_gspace(ioutput, k_gspace(irk), crys, syms)  
     if (k_gspace(irk)%length > maxlength) maxlength = k_gspace(irk)%length
     mtrxdim = k_gspace(irk)%length  
     call all_sum_all(mtrxdim)  

     write(9, 20) irk, kpoints%w(irk),(kpoints%rk(i, irk), i = 1, 3), mtrxdim  
     if (mtrxdim > maxlength_all_k) maxlength_all_k=mtrxdim
  end do
     
  do irk = 1, kpoints%nrk  
     k_gspace(irk)%maxlength_all_k=maxlength_all_k
  end do

  call myflush(9)
  
10 format(/' HAMILTON MATRIX DIMENSIONS:',/1x,27('-'), &
       &           /12x,'WEIGHT', 2x,'K-POINT',20x,'DIMENSION')

20 format(i4,4x,4f10.4,12x,i6)  

end subroutine generate_k_gspaces
