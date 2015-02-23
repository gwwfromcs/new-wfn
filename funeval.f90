function funeval(t0, plen, x, p, g, relax_myproc, &
            relaxstress, relaxforce, iwfnreuse_c, adjustpressure,&
          pw_params,altkpoints,bands,altbands,syms,energs,crys,vinit)
  !
  !     ---------------------------------------------------------------
  ! 
  use constants
  use pw_parameter_module
  use kpoint_module
  use band_module
  use symmetry_module
  use energy_module
  use crystal_module
  implicit none
  include 'flibcalls.ph'
  !
  !     INPUT/OUTPUT:
  !     ------------
  !
  type(pw_parameter) :: pw_params ! plane wave parameters
  type(kpoint) :: altkpoints  ! alternative BZ integration kpoints
  type(band) :: bands ! eigenvalues, occupation numbers etc
  type(band) :: altbands    ! another set of eigenvalues etc..
  type(symmetry) :: syms ! symmetry operations
  type(energy) :: energs    ! energies and other information
  type(crystal) :: crys     ! crystal structure
  integer, intent(inout) :: iwfnreuse_c, adjustpressure
  real(dp), intent(inout) :: p
  !
  !     INPUT:
  !     ------
  !
  integer, intent(in) ::& 
        plen,&  
       relax_myproc

  real(dp), intent(in) :: x(plen),&  
      relaxstress, relaxforce, vinit
  !
  !     OUTPUT:
  !     -------
  !
  real(dp) :: funeval
  real(dp), intent(out) :: & 
        g(plen)
  !
  !     -------------------------- local variables -----------------------
  !
  real(dp) ::  at(3, 3)
  real(dp) :: stress(6), fsig(3, 3), feps(3, 3), eps(3, 3), &
       fepsu(3, 3), pst(3, 3), etot
  real(dp) :: v0, dumdum,t0,force(3 * crys%mxdatm * crys%ntype)

  integer :: ierr
  integer :: i, j, k, joblist, njobs

  integer, parameter :: logfile_c = 7
  integer, parameter :: convergefile = 8

  real(dp) :: matinvert

  ! ------------calculate stress and forces --------------------

  if (relax_myproc == 0) then
     write(logfile_c, '(A)') 'Current state:'
     call printpacked(plen, x, crys%lvec)
  end if

  call unpack( plen, x, crys, eps)
  call apply_eps_lat(eps, crys%lvec, at)
  call mattrans(at, crys%avec) 

  if (relax_myproc == 0) then
     call myflush(logfile_c)
     call myflush(convergefile)
     if ( .not. (iand(pw_params%output(1), 268435456 ) == 268435456) &
                  .or. pw_params%ilog <= 1) then
     call mvfile('OUT.old', 'er')
     call mvfile('OUT', '.old')
    end if
  end if

  if (pw_params%correct_stress_cutoff > pw_params%emax + 1.0d-6) &
       pw_params%output(1) = ior(pw_params%output(1), 512)

  if (iand(pw_params%optimize, 16) == 0) iwfnreuse_c = 0  ! no wfn reuse is desired

  if (relaxstress > dzero) iwfnreuse_c = 0       ! no wfn reuse possible here

  joblist = 0
  njobs = 1

  call pwmain( t0, iwfnreuse_c, etot, stress(1), force(1), & 
       pw_params,altkpoints,bands,altbands,syms,energs,crys)

  iwfnreuse_c = 2
  if (iand(pw_params%optimize, 16) == 0) iwfnreuse_c = 0  ! if the user wants it.

  if (relaxstress > dzero) iwfnreuse_c = 0  ! no wfn reuse possible here

  if (pw_params%correct_stress_cutoff > pw_params%emax + 1.0d-6) then
     if ( relax_myproc == 0 ) then
        call mvfile('OUT.old', 'er')
        call mvfile('OUT', '.old')
     end if
     i =  pw_params%iscr ; j = pw_params%maxitscf ; pw_params%input = ior(pw_params%input, 1)
     pw_params%maxitscf=1
     pw_params%iscr=4
     call pwmain( t0, iwfnreuse_c, etot, stress(1), force(1), &
          pw_params,altkpoints,bands,altbands,syms,energs,crys)
       pw_params%maxitscf=j
       pw_params%iscr=i
  end if
  iwfnreuse_c = 2

  call sym2full(stress, fsig)

  do i = 1, 3
     do j = 1, 3
        feps(j, i) = eps(i, j)
     end do
     feps(i, i) = feps(i, i) + done
  end do
  v0 = vinit * matinvert(feps)

  if (relax_myproc == 0) write(logfile_c, '(A,F9.4,3(A,F16.9))') &
       'Enthalpy(',v0,')=', etot, ' + ', p / 14710.75d0 * v0, '=', &
       etot + p / 14710.75d0 * v0

  call matmul(fepsu, fsig, feps)

  if (adjustpressure == 1) then
     p = (stress(1) + stress(2) + stress(3)) * dthird * 14710.75d0 / v0
     adjustpressure = 0
  end if

  fepsu = -fepsu

  call pstress(eps, p, pst, vinit)

  if (relax_myproc == 0) then
     write(logfile_c, '(2(A,F9.4))') &
          'set pressure stress for p= ', p, ' GPa [Ry/a.u**3]:'
     call printpacked(9, pst, crys%lvec)

     write(logfile_c, '(2(A,F9.4))') &
          'actual pressure stress for p= ',&
           (stress(1) + stress(2) + stress(3)) * dthird * 14710.75d0 / v0 ,&
            ' GPa [Ry/a.u**3]:'

     write(logfile_c, '(A)') 'corrected LDA stress * Volume [Ry]:'
     call printpacked(9, fepsu, crys%lvec)
  end if

  fepsu = fepsu + pst

  ! pack it into gradient vector
  call pack(crys%ntype, crys%mxdatm, plen, crys%natom, g, force, fepsu) 

  g(1:9) = g(1:9) * relaxstress  ! correct stress
  g(10:) = -g(10:) * relaxforce  ! correct forces

  ! now zero out forces on individual atoms

  do i = 1, pw_params%nstay_put
     if (pw_params%stay_put(i) > 0 .and. pw_params%stay_put(i) <= (plen - 9) / 3) then
        g(10 + (pw_params%stay_put(i) - 1) * 3:12 + (pw_params%stay_put(i) - 1) * 3) = dzero
     else if (pw_params%stay_put(i) < 0 .and. pw_params%stay_put(i) >= -6) then
        ! Hold a strain fixed
        select case(-pw_params%stay_put(i))
        case(1)
           g(1) = dzero
        case(2)
           g(5) = dzero
        case(3)
           g(9) = dzero
        case(4)
           g(2) = dzero
           g(4) = dzero
        case(5)
           g(6) = dzero
           g(8) = dzero
        case(6)
           g(3) = dzero
           g(7) = dzero
        end select
     else
        if (relax_myproc == 0) write(logfile_c, '(A,I3)') &
             'Stay_put request ignored for atom ', pw_params%stay_put(i)
     end if
  end do

  dumdum = mdnrm2(plen, g(1), 1)
  if (relax_myproc == 0) then
     write(logfile_c, '(A,F9.6,A)') 'gradient has length ', dumdum, ':'
     call printpacked(plen, g, crys%lvec)
     call myflush(logfile_c)
     call myflush(convergefile)
  end if

  funeval = etot + p / 14710.75d0 * v0

  end function funeval


