!
subroutine dos(myproc, nproc, p, crys, bands, energs, kp, pw_params)

  use all_to_all_module  
  include 'use.h'  
  implicit none  
  include 'interface.h'  
  include 'all_to_all.h'  
  !
  !     1996 Bernd Pfrommer
  !
  !     Based on subroutines by Alberto Garcia, Paul Delaney, and
  !     Michel Cote.
  !
  ! CHANGES
  ! 
  !  11/27/00 - added , stat=ierr to allocate. Bug reported by DJR
  !             and changed p(7, to p(lmax, due to new angular.f90
  !
  !     INPUT:
  !     -----
  !
  type(crystal), intent(in) :: crys                   ! the crystal structure
  type(band), intent(in) :: bands  
  type(energy), intent(in) :: energs  
  type(kpoint), intent(in) :: kp  
  type(pw_parameter), intent(in) :: pw_params  
  integer, intent(in) :: myproc, nproc  
  !
  !     The weights for different angular momenta, different bands,
  !     kpoints, and spins. A fancy interpolation between the function
  !     values at different kpoints is done inside the subroutine
  !     spec_tet.
  !
  real(dp), intent(in) :: &
       p(llmax, crys%mxdatm, crys%ntype, bands%max, bands%nrk, crys%nspin)
  !
  !   DESCRIPTION:
  !   -----------
  !
  !   computes DOS and angular momentum resolved DOS using a
  !   tetrahedron method. The data must be given on a Monkhorst-Pack grid
  !
  !   -------------------local variables ------------------------------
  !
  integer :: n, i, j, k, l, nx, ny, nz, is, ie, it, ilm, ia, b, i1, &
       j1, k1, ibnd, nmesh, kidx(4), itet, ntetra, na, nmytet, &
       itetstart, itetend, neqtet, idum, ierr
  real(dp) :: e(0:3), f(0:3), t0, ecorn(4), fcorn(4), &
       vol, mineigen, maxeigen, delta, fav, t, ebound(2)
  real(dp), allocatable :: et(:,:,:), ados(:,:,:,:,:), edos(:,:), dummy(:)
  integer, allocatable :: itetdef(:,:)  
  real(dp) :: tetrados
  real(dp), external :: gimmetime  
  !
  !     --------------------------------------------------------------
  !
  t0 = gimmetime()  
  write(9, 110)  

  allocate(et(bands%max, bands%nrk, bands%nspin), stat=ierr)  
  call alccheck('et',bands%max*bands%nrk*bands%nspin,ierr)

  do is = 1, bands%nspin
     do i = 1, bands%nrk
        do ibnd = 1, bands%min(is)
           et(ibnd, i, is) = (bands%energy(ibnd, i, is) + &
                energs%evalshift(is) - energs%efermi(is)) * ryd
        end do
     end do
  end do
  !
  !     determine lowest and highest eigenvalues
  !
  mineigen = et(1, 1, 1)
  maxeigen = et(bands%min(1), 1, 1)
  do is = 1, bands%nspin
     do i = 1, bands%nrk
        if (et(1, i, is) < mineigen) mineigen = et(1, i, is)
        if (et(bands%min(is), i, is) > maxeigen) maxeigen = &
             et(bands%min(is), i, is)
     end do
  end do

! P. Zhang
  if(mineigen.lt.(-30.d0)) mineigen=-30.d0

!  nmesh = 2000  

   nmesh=int(100*(maxeigen-mineigen))
   if (nmesh.lt.1000) nmesh=1000

  allocate(edos(nmesh, 2),stat=ierr)  
  call alccheck('edos',nmesh*2,ierr)
  edos = dzero
  allocate(dummy(nmesh), stat =ierr)  
  call alccheck('dummy',nmesh,ierr)
  dummy = dzero

  nx = kp%grid(1)  
  ny = kp%grid(2)  
  nz = kp%grid(3)

  if (.not.associated(kp%kmap)) then  
     write(9, *) 'DOS: kpoint map is not set up!'  
     stop 
  end if
  !
  vol = pi2**3 / crys%vcell / (nx * ny * nz) / dsix  
  !
  !     rewrite the eigenvalues in eV according to fermi level
  !


  allocate(ados(nmesh, 7, crys%mxdatm, crys%ntype, crys%nspin), stat=ierr)  
  call alccheck('ados',0,ierr)

  ados = dzero


  delta = abs(maxeigen - mineigen) / real(nmesh, dp)  
  !
  !     -------- COMPUTE (ANGULAR MOMENTUM RESOLVED) DOS ----------
  !
  ebound(1) = mineigen - 1.0d1 * delta  
  ebound(2) = maxeigen + 1.0d1 * delta  

  delta = abs(ebound(2) - ebound(1)) / real(nmesh, dp)  
  
  ntetra = 6 * nx * ny * nz                            ! number of tetrahedra
  allocate(itetdef(4, ntetra), stat =ierr)  
  call alccheck('itetdef',4*ntetra,ierr)
  call setup_tetrahedra(kp%grid, kp%kmap(1), itetdef, crys%bvec)  
  !
  !     now figure out which tetrahedra are on this processor
  !
  ! minimum number of tetrahedra
  b = ntetra / nproc  
  ! number of procs with one more tet.
  na = mod(ntetra, nproc)  
  nmytet = b  

  if (na > myproc) nmytet = nmytet + 1  
  itetstart = b * myproc + min(myproc, na) + 1  
  itetend = itetstart + nmytet - 1  
  !
  !     run the loop over all tetrahedra
  !
  neqtet = 0  
  idum = 0  
  do itet = itetstart, itetend                  ! loop through all tetrahedra
     kidx = itetdef(:, itet)  
     do is = 1, crys%nspin  
        do ibnd = 1, bands%min(is)  
           ecorn(:) = et(ibnd, kidx(:), is)               ! energy at corners
           fcorn(:) = done
           call spec_tet(vol, ecorn, fcorn, nmesh, ebound, edos(1, is), &
                dummy(1), neqtet)
!           if (iand (pw_params%output(1), 1024) == 1024) then  
!              do it = 1, crys%ntype  
!                 do ia = 1, crys%natom(it)  
!                    do ilm = 1, 7  
!                       ecorn(:) = et(ibnd, kidx(:), is)   ! energy at corners
!                       fcorn(:) = p(ilm, ia, it, ibnd, kidx(:), is)  
!                       call spec_tet(vol, ecorn, fcorn, nmesh, ebound, &
!                            dummy(1), ados(1, ilm, ia, it, is), idum)
!                    end do
!                 end do
!              end do
!           end if
        end do
     end do
  end do
  neqtet = neqtet / (crys%nspin * bands%min(1))  
  call all_sum_all(neqtet)  
  write(9, 210) ntetra  
  write(9, 200) neqtet  
  !
  !     now sum across all processors
  !
  do is = 1, crys%nspin  
     dummy(:) = edos(:, is)  
     call all_sum_all(dummy, nmesh)  
     edos(:, is) = dummy(:)  
  end do
!  if (iand(pw_params%output(1), 1024) == 1024) then  
!     do is = 1, crys%nspin  
!        do it = 1, crys%ntype  
!           do ia = 1, crys%natom(it)  
!              do ilm = 1, 7  
!                 dummy(:) = ados(:, ilm, ia, it, is)  
!                 call all_sum_all(dummy, nmesh)  
!                 ados(:, ilm, ia, it, is) = dummy(:)  
!              end do
!           end do
!        end do
!     end do
!  end if
  if (myproc /= 0) then  
     deallocate(itetdef)  
     deallocate(ados)  
     deallocate(edos)  
     deallocate(et)  
     deallocate(dummy)
     return  
  end if
  !
  !     ------------ print out ---------------------------------------
  !

  if(bands%nspin.eq.1) then
    open(11, file = 'DOS', status = 'unknown', form = 'formatted')  
  end if

  if(bands%nspin.eq.2) then
    open(11, file = 'DOS1', status = 'unknown', form = 'formatted')  
    open(12, file = 'DOS2', status = 'unknown', form = 'formatted')  
    open(13, file = 'DOS', status = 'unknown', form = 'formatted')  
  end if

  write(11, *) '# Units:  states/eV/unit cell vs. eV'
  write(11, '(a15,2f12.6)') '# FERMI LEVELS:', energs%efermi(1) , &
       energs%efermi(2)

  write(11, '(a17)') '# SPIN DENSITIES:'  
  do is = 1, bands%nspin  
     do ie = 1, nmesh  
        write(10+is, * ) ebound(1) + ie * delta, edos(ie, is) * &
             crys%vcell / pi2**3 / bands%nspin
     end do
!     write(11, * ) '&'  
     write(10+is, * )   
  end do
  if (bands%nspin > 1) then  
     write(13, '(a100)') '# sum of spin densities:'  
     do ie = 1, nmesh  
        write(13, * ) ebound(1) + ie * delta, edos(ie, 1) * &
             crys%vcell / pi2**3 / bands%nspin + &
             edos(ie, 2) * crys%vcell / pi2**3 / bands%nspin
     end do
!     write(11, * ) '&'  
     write(13, * ) 
  end if
  !
  !     ------------------------------------------------------------------
  !
  !
  !     print angular resolved DOS if it is requested
  !
!  if (iand(pw_params%output(1), 1024) == 1024) then  
!     do is = 1, bands%nspin  
!        do it = 1, crys%ntype  
!           do ia = 1, crys%natom(it)  
!              do ilm = 1, 7  
!                 if (ilm == 1) then  
!                    write(11, 124) crys%nameat(it), ia, is  
!                 else if (ilm < 4) then  
!                    write(11, 123) crys%nameat(it), ia, ilm - 2, is  
!                 else if (ilm == 7) then  
!                    write(11, 123) crys%nameat(it), ia, 2, is  
!                 else if (ilm < 7) then  
!                    write(11, 125) crys%nameat(it), ia, ilm - 5, is  
!                 end if
!                 do ie = 1, nmesh  
!                    write(11, *) ebound(1) + ie * delta, &
!                         ados(ie, ilm, ia, it, is) * crys%vcell / &
!                         pi2**3 / real(bands%nspin, dp)
!                 end do
!                 write(11, *) '&'  
!              end do
!           end do
!        end do
!     end do
!  end if
  !
  close(11)  

  write(9, 987)  
  if (iand(pw_params%output(1), 8) == 8) then  
     write (9, 988) gimmetime() - t0  
  end if

  deallocate(itetdef)  
  deallocate(ados)  
  deallocate(edos)  
  deallocate(et)
  deallocate(dummy)
  
  return  

110 format(/' GENERATION OF DOS',/1x,17('-'))  
123 format('#',a2,1x,i4,' with angular momentum:',i2,', spin:',i2)  
124 format('#',a2,1x,i4,' local density of states',', spin:',i2)  
125 format('#',a2,1x,i4,' p-state with m-quantum number:', &
       &     i2,', spin:',i2)
200 format(' NUMBER OF TETRAHEDRA WITH EQUAL CORNERS:',i10)  

210 format(' TOTAL NUMBER OF TETRAHEDRA:             ',i10)  
987 format(/' WROTE DENSITY OF STATES TO FILE "DOS"'/)  

988 format(/' TIME FOR DENSITY OF STATES: ',f14.4)  

end subroutine dos
!
!     =================================================================
!
subroutine setup_tetrahedra(grid, k_index, itetra, bvec)  
  !
  !     1996 Bernd Pfrommer
  !
  !     originated from a code by Paul Delaney
  !
  use constants
  implicit none  
  !
  !     INPUT:
  !     -----
  !
  integer, intent(in) :: grid(3), &  ! the Monkhorst pack k-grid size kx,ky,kz
       k_index(0:grid(1)-1, 0:grid(2)-1, 0:grid(3)-1)         ! map to irred k
  real(dp), intent(in) :: bvec(3, 3)  ! the reciprocal lattice vectors, as col
                                                               ! in the matrix
  !
  !     OUTPUT:
  !     ------
  !
  integer, intent(out) :: &
       itetra(4, 6 * grid(1) * grid(2) * grid(3))     ! index into irred. k-pt
  !
  !
  !     DESCRIPTION:
  !     -----------
  !
  !     1) Finds optimal 6-tetrahedral disection of the monkhorst-pack
  !     grid, such that the edgelengths of the tetrahedra are as
  !     regular as possible
  !
  !     2) With this, generates the index array itetra which maps each
  !     corner of the 6*nx*ny*nz tetrahedra to an irreducible kpoint
  !
  !
  !
  !     -------------------------- local variables ------------------
  !
  integer :: ntetra, i, j, k, it, iv, idx, idy, idz, ivert, tset(3, 4, 6), &
       ivert2, ivert1, goodi, goodj, goodk
  real(dp) :: thisedgemin, thisedgemax, cart(3, 4), edge, edgemin, edgemax
  !
  !
  !

!!$  integer, parameter :: initialtetraoffset(3, 4, 6) = (/ &
!!$       0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, &
!!$       0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, &
!!$       1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, &
!!$       1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1 /)

  ! This way in F90 - CJP 16/04/00

  integer, parameter :: initialtetraoffset(3,4,6) = reshape(source=(/ &
       0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, &
       0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, &
       1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, &
       1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1 /),shape=(/3,4,6/))

  edgemin = dzero
  edgemax = 1.0d10  

  !     First set up the matrix of components of the reciprocal lattice
  !     vectors along the Cartesian basis. bcomp(i,j) holds the dot
  !     product of dual basis vector b_i with Cartesian basis vector
  !     e_j.
  !     Now go through all the four long diagonals, specified by the
  !     corner (i,j,k)

  k = 0  
  do i = 0, 1  
     do j = 0, 1  
        do ntetra = 1, 6  
           do ivert = 1, 4  
              call mirrorijk(initialtetraoffset(1, ivert, ntetra), &
                   tset (1, ivert, ntetra), i, j, k)
           end do
        end do
        thisedgemin = 1.0d10  
        thisedgemax = dzero
        do ntetra = 1, 6  
           do ivert = 1, 4  
              cart(:, ivert) = matmul(bvec, real(tset(:, ivert, ntetra), dp))  
           end do
           do ivert1 = 1, 3  
              do ivert2 = ivert1 + 1, 4  
                 edge = (cart(1, ivert1) - cart(1, ivert2))**2 + &
                      (cart(2, ivert1) - cart(2, ivert2))**2 + &
                      (cart(3, ivert1) - cart(3, ivert2))**2
                 edge = sqrt(edge)  
                 thisedgemax = max(thisedgemax, edge)  
                 thisedgemin = min(thisedgemin, edge)  
              end do
           end do
        end do
        if (thisedgemax < edgemax) then  
           goodi = i  
           goodj = j  
           goodk = k  
           edgemax = thisedgemax  
           edgemin = thisedgemin  
        end if
     end do
  end do
  do ntetra = 1, 6  
     do ivert = 1, 4  
        call mirrorijk(initialtetraoffset(1, ivert, ntetra), &
             tset (1, ivert, ntetra), goodi, goodj, goodk)
     end do
  end do
  !      write(0,*)
  !      write(0,*) 'Choosing long-diagonal based at i,j,k=',
  !     $     goodi,goodj,goodk
  !      write(0,*) 'The longest edge length is ',edgemax,' a.u.'
  !      write(0,*) 'The shortest is ',edgemin,' a.u.'
  !     ------------------------------------------------------------------
  !     ----------- now use the optimial tetrahedron to set up the -------
  !     -------------------------- index array ---------------------------
  !
  !      write(9,*) 'TETRAHEDRON GENERATION:'

  ntetra = 0  
  do i = 0, grid(1) - 1  
     do j = 0, grid(2) - 1  
        do k = 0, grid(3) - 1  
           do it = 1, 6  
              !                  write(9,*) 'tetrahedron:',it
              do iv = 1, 4  
                 idx = mod(i + tset(1, iv, it), grid(1))  
                 idy = mod(j + tset(2, iv, it), grid(2))  
                 idz = mod(k + tset(3, iv, it), grid(3))  
                 itetra(iv, ntetra + it) = abs(k_index(idx, idy, idz))  
                 !               write(9,*) idx,idy,idz,itetra(iv,ntetra+it)
              end do
           end do
           ntetra = ntetra + 6  
        end do
     end do
  end do

end subroutine setup_tetrahedra
!
!     =================================================================
!
subroutine mirrorijk(in, out, i, j, k)  
  !
  !     Do mirrors in i around i=.5 if i=1; similarly for j,k
  !
  !     1995/96 by Paul Delaney
  !
  integer, intent(in) :: i, j, k  
  integer, intent(in) :: in(3)
  integer, intent(out) :: out (3)  
  out (1) = in (1)  
  out (2) = in (2)  
  out (3) = in (3)  
  if (i == 1) out(1) = 1 - out(1)  
  if (j == 1) out(2) = 1 - out(2)  
  if (k == 1) out(3) = 1 - out(3)
  
  return  

end subroutine mirrorijk
