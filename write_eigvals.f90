!
subroutine write_eigvals(myproc, pw_params, crys, energs, kp, bands, p)
  !
  use all_to_all_module  
  include 'use.h'  
  implicit none  
  include 'interface.h'  
  include 'all_to_all.h'  
  !
  !     1996 Bernd Pfrommer
  ! CHANGES
  ! 
  !  11/27/00 -  changed p(7, to p(lmax, due to new angular.f90 
  !
  !     INPUT:
  !     -----
  !
  type(crystal), intent(in) ::  crys    ! the crystal structure
  type(band), intent(in) :: bands  
  type(energy), intent(in) :: energs  
  type(kpoint), intent(in) :: kp  
  type(pw_parameter), intent(in) :: pw_params  
  integer, intent(in) :: myproc  
  !
  !     The weights for different angular momenta, different bands,
  !     kpoints, and spins. A fancy interpolation between the function
  !     values at different kpoints is done inside the subroutine
  !     spec_tet.
  !
  real(dp), intent(in) :: p (llmax, crys%mxdatm, crys%ntype, bands%max, &
       bands%nrk, crys%nspin)
  !
  !     DESCRIPTION:
  !     -----------
  !
  !     writes a file with all the eigenvalues, the weights of the
  !     kpoints, and, if applicable, the angular momentum decomposition
  !
  !
  !     -------------------local variables ------------------------------
  !
  integer :: irk, is, n, it, ia  

  write(9, 110)  
  if (myproc /= 0) return  
  !
  !     ------------ print out ---------------------------------------
  !
  open(11, file = 'EIG', status = 'unknown', form = 'formatted')  
!  write(11, '(a,2g20.12)') '# FERMI LEVELS (includes shift):', &
!       energs%efermi(1), energs%efermi(2)
  write(11,*) '# FERMI LEVELS (includes shift):'
  write(11,*)  energs%efermi(1), energs%efermi(2)
  write(11,*) '# SHIFT TO BE ADDED TO EIGENVALUES DUE TO VXC0:'
  write(11,*) energs%evalshift(1), energs%evalshift(2)
!  write(11, '(a,a,2g20.12)') &
!       '# SHIFT TO BE ADDED TO EIGENVALUES DUE', &
!       ' TO VXC0:', energs%evalshift(1), energs%evalshift(2)
  write(11, *) 'number of kpoints:'  
  write(11, *) kp%nrk  
  write(11, *) 'weights of the kpoints:'  
  do irk = 1, kp%nrk  
     write(11, *) irk, kp%w(irk)  
  end do
  write(11, *) 'number of spins:'  
  write(11, *) bands%nspin  
  write(11, *) 'max number of bands:'  
  write(11, *) bands%max  
  write(11, *) 'computed number of bands:'  
  write(11, *) bands%nband  
  write(11, *) 'min number of bands:'  
  write(11, *) bands%min  
  write(11, *) 'eigenvalues in rydbergs:'  
  do is = 1, bands%nspin  
     do irk = 1, bands%nrk  
!        write(11, *) irk, bands%nband(irk, is), &
!             bands%energy(1:bands%nband(irk, is), irk, is)
        write(11, *) irk, bands%nband(irk, is)
        do n=1,bands%nband(irk, is)
          write(11, *)  n,bands%energy(n, irk, is)
        end do
     end do
  end do
  write(11, *) 'angular momentum decomposition info:'  
  write(11, *) (iand(pw_params%output(1), 1024) == 1024)  
  if (iand(pw_params%output(1), 1024) == 1024) then  
     write(11, *) 'number of atomic types:'  
     write(11, *) crys%ntype  
     write(11, *) 'max number of atoms:'  
     write(11, *) crys%mxdatm  
     write(11, *) 'number of atoms for each specie:'  
     write(11, *) crys%natom  
     do is = 1, bands%nspin  
        write(11, *) 'spin:'  
        write(11, *) bands%nspin  
        do irk = 1, bands%nrk  
           write(11, *) 'k-point:'  
           write(11, *) irk  
           do n = 1, bands%nband(irk, is)  
              write(11, *) 'band:'  
              write(11, *) n  
              do it = 1, crys%ntype  
                 write(11, *) 'atom type:'  
                 write(11, *) it  
                 do ia = 1, crys%natom(it)  
                    write(11, *) ia, p(1:7, ia, it, n, irk, is)  
                 end do
              end do
           end do
        end do
     end do
  end if
  close(11)
  
  return  

110 format(/' WRITING EIGENVALUES TO FILE EIG ',/1x,31('-'))  

end subroutine write_eigvals
