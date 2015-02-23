!     @process extchk
!
subroutine logfile(ipr, myproc, nproc, ilog, totime, fsum, ssum, &
     crys, pw_params, syms, kpoints, energs)
  !     ----------------------------------------------------------------
  !         print to the logfile unit 9
  !
  include 'use.h'  
  implicit none              ! implicit? No.
  include 'interface.h'

  type(crystal), intent(in) :: crys  
  type(pw_parameter), intent(in) :: pw_params  
  type(symmetry), intent(in) :: syms  
  type(kpoint), intent(in) :: kpoints  
  type(energy), intent(in) :: energs  
  !
  !
  integer, intent(in) :: ipr, &    ! print flag
       nproc, &                    ! number of processors used
       myproc, &
       ilog                        ! control flag/iteration number
  real(dp), intent(in) :: &
       totime, &                   ! total time
       fsum(3, crys%mxdatm, crys%ntype), &  ! symmetrized forces
       ssum(6)                              ! symmetrized stress
  !
  !     ------------------ local variables -----------------------------
  !
  integer :: i, j, k, natomi, n  
  real(dp) :: fac 

  include 'release.h'  

  if (myproc > 0) return  

  if (ilog < 0) return  
  if (ilog <= 0) then  
     write(34, 5) release_string, date_string, author_string  
  end if
  write(34, 10) ilog  
  write(34, 20) syms%ntrans  

  write(34, 30) energs%total  

  fac = 14710.80 / crys%vcell  
  write(34, 9130) crys%avec  

  write(34, 9050) crys%vcell  
  write(34, 9020)  
  !
  do i = 1, crys%ntype  
     natomi = crys%natom(i)  
     write(34, 9029) crys%nameat(i), crys%szet(i)  
     do j = 1, natomi  
        write(34, 9030) (crys%rat(k, j, i) / pi2, k = 1, 3)  
     end do
  end do
9029 format('newtype',1x,a2,1x,f9.5)  

9030 format('coord',5x,3f18.14)  
  if (ilog <= 0) then  
     write(34, 45) nproc  
     write(34, 48) kpoints%grid(1), kpoints%grid(2), kpoints%grid(3)
     write(34, 49) kpoints%shift(1), kpoints%shift(2), kpoints%shift(3)
     write(34, 50) kpoints%nrk  
     write(34, 70) pw_params%emax  
     if (pw_params%ekinmod(1) > dzero) write(34, 75) pw_params%ekinmod(1), &
          pw_params%ekinmod(2), pw_params%ekinmod(3)
     write(34, 60) pw_params%smearing  
     if (pw_params%mfield%gvec(4) /= 0) then  
        write(34, 550) pw_params%mfield%gvec(1), pw_params%mfield%gvec(2), &
             pw_params%mfield%gvec(3), pw_params%mfield%h
     end if
  end if
  write(34, 560) energs%afmagmom  

  write(34, 570) energs%magnetic  
  write(34, 40)  
  n = 0  
  do i = 1, crys%ntype  
     natomi = crys%natom(i)  
     do j = 1, natomi  
        n = n + 1  
        write(34, 41) n, (crys%rat(k, j, i) / pi2, k = 1, 3), &
             ((crys%bvec(k, 1) * fsum(1, j, i) + &
             crys%bvec(k, 2) * fsum(2, j, i) + &
             crys%bvec(k, 3) * fsum(3, j, i)) / pi2, k = 1, 3)
     end do
  end do
  write(34, 42) (ssum(i), i = 1, 6)  

  write(34, 43) (ssum(i) * fac, i = 1, 6)  
  write(34, 80) totime, totime / 3.6d3  

  call myflush(34)

5 format(2x,'paratec release ',a,' ',a,' ',a)  
10 format(//' ------------ Iteration number: ',i3, &
       &     ' ------------'/)
20 format(' number of symmetry operations:',i3)  

30 format(' total energy (not enthalpy) =',f17.10)  
40 format(/, &
       &     ' forces :',15x,'coord',19x, &
       &     'cartesian force [Ry/a.u.]',/ &
       &     ' ',11x,'a1',4x,'a2',4x,'a3',11x, &
       &     '-x-',12x,'-y-',12x,'-z-')
41 format(2x,i3,3x,3f6.3,2x,3f16.12)  
42 format(//,' stress :', &
       &     25x,'sigma * v (ry)',/, &
       &     '(cartes)',2x,'-xx-',8x,'-yy-',8x,'-zz-',8x, &
       &     '-xy-',8x,'-yz-',8x,'-zx-',/, &
       &     5x,6f12.5)
43 format(/,36x,'sigma (gpa)',/, &
       &     10x,'-xx-',8x,'-yy-',8x,'-zz-',8x, &
       &     '-xy-',8x,'-yz-',8x,'-zx-',/, &
       &     5x,6f12.5)

45 format(/' number of processors used: ',i4)  
48 format(/' grid of k-points: ',3i4)  
49 format(/' k-grid shift: ',3f9.4)  
50 format(/' number of irreducible k-points:',i4)  
60 format(' gaussian smearing [eV]:', f10.5)  
70 format(' energy cutoff [Ry]:    ', f10.5)  

75 format(' modification to kinetic energy: a = ',f9.4, &
       &     ', e0 = ',f9.4,', sigma = ',f9.4)
80 format(/' total time:            ', f10.5,'s   =',f10.5,'h')  
550 format(/' MAGNETIC FIELD WITH G= ',3i3,' AND STRENGHT ',f12.6)  
560 format(/' antiferromagnetic moment [Bohr magneton]:', f10.5)  

570 format(' magnetic energy:', f10.5)  
9050 format(/' volume ',f16.10,' au.')  
9020 format(/1x,'atomic positions:')  

9130 format(/' lattice vectors (a.u.)',/3(/'coord',2x,3(1x,f18.14)))  

end subroutine logfile
