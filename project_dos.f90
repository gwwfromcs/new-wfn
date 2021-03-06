subroutine project_dos(myproc, nproc, cnk_all_ibz,syms,  crys, bands, energs, kp,  &
             nband0, nbandmin,nsites,lproj1,mproj,lproj,max_proj,info,nproj, &
                     ind_site_atom,ind_site_isite,ind_isite_site)



  use all_to_all_module  
  include 'use.h'  
  implicit none  
  include 'interface.h'  
  include 'all_to_all.h'  
  !
  !     INPUT:
  !     -----
  !
  type(crystal), intent(in) :: crys                   ! the crystal structure
  type(band), intent(in) :: bands  
  type(energy), intent(in) :: energs  
  type(kpoint), intent(in) :: kp  
  type(symmetry), intent(in) :: syms
  integer, intent(in) :: nband0,nsites,max_proj,lproj1(nsites),mproj(crys%ntype),lproj(max_proj,crys%ntype)
  integer, intent(in) :: myproc, nproc, nbandmin,nproj(max_proj,crys%ntype),info
  integer, intent(in) :: ind_site_atom(500), ind_site_isite(500),ind_isite_site(5,500)

  !
  !     The weights for different angular momenta, different bands,
  !     kpoints, and spins. A fancy interpolation between the function
  !     values at different kpoints is done inside the subroutine
  !     spec_tet.
  !
!
! projection up to L=2, each projection is considered a site
! for example, if two projects are done on an atomic site,
! they are considered two projection sites.
!
  double complex, intent(in) :: cnk_all_ibz(5, nbandmin-nband0,bands%nrk,crys%nspin,nsites)
  !
  !   DESCRIPTION:
  !   -----------
  !
  !   computes angular momentum resolved DOS using a
  !   tetrahedron method. The data must be given on a Monkhorst-Pack grid
  !
  !   -------------------local variables ------------------------------
  !
  character*2 symbol
  character*1 symbol_tem

  integer :: n, i, j, k, l, nx, ny, nz, is, ie, it, ilm, ia, b, i1, &
       j1, k1, ibnd, nmesh, kidx(4),kidx_ibz(4), itet, ntetra, na, nmytet, &
       itetstart, itetend, neqtet, idum, ierr,isite,max_L, max_LM

  parameter(max_L=2,max_LM=max_L*2+1)

  real(dp) :: e(0:3), f(0:3), t0, ecorn(4), fcorn(4), &
       vol, mineigen, maxeigen, delta, fav, t, ebound(2)
  integer, allocatable :: itetdef(:,:) , itetdef_fbz(:,:) 
  real(dp), allocatable :: et(:,:,:), ados(:,:,:,:), edos(:,:), dummy(:),dnk_all(:,:,:,:,:),qtrand(:,:,:),vtrand(:,:,:)
  real(dp), allocatable :: int_ados(:,:,:,:)
  real(dp) :: tetrados
  real(dp), external :: gimmetime  
  double complex,allocatable::cnk_all_fbz(:,:,:,:,:)

  !
  !     --------------------------------------------------------------
  !
  t0 = gimmetime()  
  write(9, 110)  

!----------------------------------------------
! rotate projection if needed
!
  allocate(dnk_all(max_LM,nbandmin-nband0,kp%nrk_fbz,crys%nspin,nsites))

  if(kp%reduced) then
    allocate(cnk_all_fbz(max_LM,nbandmin-nband0,kp%nrk_fbz,crys%nspin,nsites))

!    do i=1,kp%nrk_fbz 
!    write(9,*) i,kp%ind_symtr(i),kp%kmap(i)
!    end do

    allocate(qtrand(5,5,syms%ntrans))
    allocate(vtrand(3,3,syms%ntrans))

    vtrand(:,:,1:syms%ntrans)=syms%rsymmat(:,:,1:syms%ntrans)


    call findtransd(syms%ntrans,vtrand,qtrand)

     call rotate_cnk(lproj1,nband0,nbandmin,vtrand,qtrand,cnk_all_ibz, &
                     cnk_all_fbz,syms,kp, nsites, bands%nspin,ind_site_atom,ind_site_isite,ind_isite_site)

!    do it = 1, nsites
!    do is = 1, bands%nspin
!       call rotate_cnk(lproj1(it),nband0,nbandmin,vtrand,qtrand,cnk_all_ibz(1,1,1,is,it), &
!                       cnk_all_fbz(1,1,1,is,it),syms,kp,is,it, &
!                     ind_site_atom,ind_site_isite,ind_isite_site)
!    end do
!    end do

    dnk_all=cnk_all_fbz*DCONJG(cnk_all_fbz)

  else

  dnk_all=cnk_all_ibz*DCONJG(cnk_all_ibz)

  end if


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


  mineigen = minval(et(nband0+1, :,:))
  maxeigen = maxval(et(bands%min(1), :,:))
  if(mineigen.lt.(-30.d0)) mineigen=-30.d0

!  do is = 1, bands%nspin
!     do i = 1, bands%nrk
!        if (et(1, i, is) < mineigen) mineigen = et(1, i, is)
!        if (et(bands%min(is), i, is) > maxeigen) maxeigen = &
!             et(bands%min(is), i, is)
!     end do
!  end do

!  nmesh = 2000  

   nmesh=int(1000*(maxeigen-mineigen)/15)
   if (nmesh.lt.1000) nmesh=1000

!  allocate(edos(nmesh, 2),stat=ierr)  
!  call alccheck('edos',nmesh*2,ierr)
!  edos = dzero

  allocate(dummy(nmesh), stat =ierr)  
  call alccheck('dummy',nmesh,ierr)
  dummy = dzero

  nx = kp%grid(1)  
  ny = kp%grid(2)  
  nz = kp%grid(3)

  if (.not.associated(kp%kmap)) then  
     write(9, *) 'DOS: kpoint map is not set up!'  
     call mystop()
  end if
  !
  vol = pi2**3 / crys%vcell / (nx * ny * nz) / dsix  
  !
  !     rewrite the eigenvalues in eV according to fermi level
  !


  allocate(ados(nmesh, 5, nsites, crys%nspin), stat=ierr)  
  allocate(int_ados(nmesh, 5, nsites, crys%nspin), stat=ierr)  
  call alccheck('ados',0,ierr)
  call alccheck('int_ados',0,ierr)

  ados = dzero
  int_ados = dzero


  delta = abs(maxeigen - mineigen) / real(nmesh, dp)  
  !
  !     -------- COMPUTE (ANGULAR MOMENTUM RESOLVED) DOS ----------
  !
  ebound(1) = mineigen - 1.0d1 * delta  
  ebound(2) = maxeigen + 1.0d1 * delta  

  delta = abs(ebound(2) - ebound(1)) / real(nmesh, dp)  
  
  ntetra = 6 * nx * ny * nz                            ! number of tetrahedra
  allocate(itetdef(4, ntetra), stat =ierr)  
  allocate(itetdef_fbz(4, ntetra), stat =ierr)  
  call alccheck('itetdef',4*ntetra,ierr)
  call setup_tetrahedra_pdos(kp%grid, kp%kmap(1), itetdef, itetdef_fbz,crys%bvec)  
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
     kidx(:) = itetdef_fbz(:, itet)
     kidx_ibz(:) = itetdef(:, itet)
     do is = 1, crys%nspin  
        do ibnd = nband0+1, nbandmin+nband0
           ecorn(:) = et(ibnd, kidx_ibz(:), is)               ! energy at corners
           do it = 1, nsites
              do ilm = 1, 2*lproj1(it)+1
                 fcorn(:) = dnk_all(ilm,ibnd-nband0,kidx(:),is,it)
!                 fcorn(:) = 1.d0
                 call spec_tet(vol, ecorn, fcorn, nmesh, ebound, &
                         dummy(1), ados(1, ilm, it, is), neqtet)
              end do
           end do

        end do
     end do
  end do
  neqtet = neqtet / (crys%nspin * nbandmin)  
  call all_sum_all(neqtet)  
  !
  !     now sum across all processors
  !

     do is = 1, crys%nspin  
        do it = 1, nsites
           do ilm = 1, 2*lproj1(it)+1
              dummy(:) = ados(:, ilm, it, is)  
              call all_sum_all(dummy, nmesh)  
              ados(:, ilm, it, is) = dummy(:)  
           end do
        end do
     end do

  if (myproc /= 0) then  
     deallocate(itetdef)  
     deallocate(ados)  
     deallocate(int_ados)  
!     deallocate(edos)  
     deallocate(et)  
     deallocate(dummy)
     return  
  end if
  !
  !
  !     ------------------------------------------------------------------
  !
     ados=ados* crys%vcell/pi2**3 / real(bands%nspin, dp)

! integrated pdos 08/05/2012 shanghai

     do is = 1, bands%nspin  
        isite=0

! mproj(itype): how many projectors for each atomic type
! lproj(i,itype): i=1,mproj(itype): angular momentum of each projector
! nproj(i,itype): principle quantum number


        do it = 1, crys%ntype

           symbol=crys%nameat(it)(1:2)
           if(mproj(it).ge.1) then
              do i=1,mproj(it)


                 if(lproj(i,it).eq.0) symbol_tem="s"
                 if(lproj(i,it).eq.1) symbol_tem="p"
                 if(lproj(i,it).eq.2) symbol_tem="d"
                 if(lproj(i,it).eq.3) symbol_tem="f"

                 if(symbol(2:2)==' ') then
                 if(info.eq.1) then
                    open(339,file=symbol(1:1)//CHAR(it+48)//"_"//CHAR(nproj(i,it)+48)//symbol_tem//"_dos_int"//CHAR(is+48)//".dat")
                 else
                    open(339,file=symbol(1:1)//CHAR(it+48)//"_"//symbol_tem//"_dos_int"//CHAR(is+48)//".dat")
                 end if

                 else
                 if(info.eq.1) then
                    open(339,file=symbol(1:1)//symbol(2:2)//CHAR(it+48)//"_"//CHAR(nproj(i,it)+48)//symbol_tem//"_dos_int"//CHAR(is+48)//".dat")
                 else
                    open(339,file=symbol(1:1)//symbol(2:2)//CHAR(it+48)//"_"//symbol_tem//"_dos_int"//CHAR(is+48)//".dat")
                 end if
                 end if


! now integrate pdos
! not done yet
                 do j=isite+1,isite+crys%natom(it)
                 do ilm=1,lproj(i,it)*2+1

! first two points
                    int_ados(1,ilm,j,is)=0.d0
                    int_ados(2,ilm,j,is)=(ados(1,ilm,j,is)+ados(2,ilm,j,is))*delta/2.d0

                    do ie=1,nmesh-2
                    int_ados(ie+2,ilm,j,is)=int_ados(ie,ilm,j,is)+ (ados(ie,ilm,j,is)+4.d0*ados(ie+1,ilm,j,is)+ ados(ie+2,ilm,j,is))*delta/3.d0
                    end do
                 end do
                 end do

                 do ie = 1, nmesh  
                    write(339, 110) ebound(1) + ie * delta, &
                         ((int_ados(ie, ilm, j,is),ilm=1,lproj(i,it)*2+1),j=isite+1,isite+crys%natom(it)), &
                         (SUM(int_ados(ie, 1:lproj(i,it)*2+1,j,is)),j=isite+1,isite+crys%natom(it))
                 end do

                 isite=isite+crys%natom(it) 
                 close(339)  

              end do
           end if
        end do
     end do
  !

     do is = 1, bands%nspin  
        isite=0

! mproj(itype): how many projectors for each atomic type
! lproj(i,itype): i=1,mproj(itype): angular momentum of each projector
! nproj(i,itype): principle quantum number


        do it = 1, crys%ntype

           symbol=crys%nameat(it)(1:2)
           if(mproj(it).ge.1) then
              do i=1,mproj(it)


                 if(lproj(i,it).eq.0) symbol_tem="s"
                 if(lproj(i,it).eq.1) symbol_tem="p"
                 if(lproj(i,it).eq.2) symbol_tem="d"
                 if(lproj(i,it).eq.3) symbol_tem="f"

                 if(symbol(2:2)==' ') then
                 if(info.eq.1) then
                    open(239,file=symbol(1:1)//CHAR(it+48)//"_"//CHAR(nproj(i,it)+48)//symbol_tem//"_dos"//CHAR(is+48)//".dat")
                 else
                    open(239,file=symbol(1:1)//CHAR(it+48)//"_"//symbol_tem//"_dos"//CHAR(is+48)//".dat")
                 end if

                 else
                 if(info.eq.1) then
                    open(239,file=symbol(1:1)//symbol(2:2)//CHAR(it+48)//"_"//CHAR(nproj(i,it)+48)//symbol_tem//"_dos"//CHAR(is+48)//".dat")
                 else
                    open(239,file=symbol(1:1)//symbol(2:2)//CHAR(it+48)//"_"//symbol_tem//"_dos"//CHAR(is+48)//".dat")
                 end if
                 end if


                 do ie = 1, nmesh  
                    write(239, 110) ebound(1) + ie * delta, &
                         ((ados(ie, ilm, j,is),ilm=1,lproj(i,it)*2+1),j=isite+1,isite+crys%natom(it)), &
                         (SUM(ados(ie, 1:lproj(i,it)*2+1,j,is)),j=isite+1,isite+crys%natom(it))
                 end do
                 isite=isite+crys%natom(it) 
                 close(239)  

              end do
           end if
        end do
     end do
  !


  deallocate(itetdef)  
  deallocate(ados)  
!  deallocate(edos)  
  deallocate(et)
  deallocate(dummy)
  
  return  

110 format(100f12.6)

200 format(' NUMBER OF TETRAHEDRA WITH EQUAL CORNERS:',i10)  
210 format(' TOTAL NUMBER OF TETRAHEDRA:             ',i10)  


end subroutine project_dos

subroutine setup_tetrahedra_pdos(grid, k_index, itetra, itetra_fbz,bvec)  

!----------------------------------------------------------
! 2007 P. Zhang
! Modified for calculating projected DOS
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
  integer, intent(out) :: &
       itetra_fbz(4, 6 * grid(1) * grid(2) * grid(3))     ! index into irred. k-pt
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
                 itetra_fbz(iv, ntetra + it) = idx+1+idy*grid(1)+idz*grid(1)*grid(2)
!                                write(9,*) itetra_fbz(iv,ntetra+it),itetra(iv,ntetra+it)
              end do
           end do
           ntetra = ntetra + 6  
        end do
     end do
  end do

end subroutine setup_tetrahedra_pdos

   subroutine  rotate_cnk(lproj,nband0,nbandmin,vtrand,qtrand,cnk_all_ibz, &
                     cnk_all_fbz,syms,kp, nsites, nspin,ind_site_atom,ind_site_isite,ind_isite_site)


! subroutine rotate_cnk(lproj,nband0,nbandmin,vtrand,qtrand,cnk_all_ibz,cnk_all_fbz,syms,kp,is,it, &
!                     ind_site_atom,ind_site_isite,ind_isite_site)



  include 'use.h'
  implicit none

  type(symmetry), intent(in) :: syms
  type(kpoint) kp

  integer::lproj(nsites),nbandmin,nband0,nsites,nspin
  integer, intent(in) :: ind_site_atom(500), ind_site_isite(500),ind_isite_site(5,500)

  double precision::vtrand(3,3,syms%ntrans),qtrand(5,5,syms%ntrans)
  double complex::cnk_all_ibz(5,nbandmin-nband0,kp%nrk,nspin,nsites)
  double complex::cnk_all_fbz(5,nbandmin-nband0,kp%nrk_fbz,nspin,nsites)

! local variables
 
  integer::ik,ikk,ib,ikt,it,is,itt,iit,imap
  double complex::cnk_tem(5)

  do it = 1, nsites
  do is = 1, nspin


  do ik=1,kp%nrk_fbz
     cnk_all_fbz(1:2*lproj(it)+1,1:nbandmin-nband0,ik,is,it)=cnk_all_ibz(1:2*lproj(it)+1,1:nbandmin-nband0,abs(kp%kmap(ik)),is,it) 
  end do

! now rotate cnk if necessary

  do ik=1,kp%nrk_fbz
     ikk=kp%kmap(ik)

     ikt=abs(kp%ind_symtr(ik))

! this sym op maps atom it to imap
     imap=syms%ind_rot(ikt,ind_site_atom(it))

     iit=ind_site_isite(it)
     itt=ind_isite_site(iit,imap)
!     write(9,*) "it,ik,itt,ikk",it,ik,itt,ikk



!     write(9,*) ik,kp%kmap(ik),kp%ind_symtr(ik),ikt
     if(ikk.gt.0) then

     do ib=1,nbandmin-nband0
        cnk_tem(1:2*lproj(it)+1)=cnk_all_ibz(1:2*lproj(it)+1,ib,ikk,is,itt)
        if(lproj(it).eq.2) then
           cnk_all_fbz(1:2*lproj(it)+1,ib,ik,is,it)=MATMUL(qtrand(:,:,ikt),cnk_tem(1:2*lproj(it)+1))
        end if
        if(lproj(it).eq.1) then
           cnk_all_fbz(1:2*lproj(it)+1,ib,ik,is,it)=MATMUL(vtrand(:,:,ikt),cnk_tem(1:2*lproj(it)+1))
        end if
        if(lproj(it).eq.0) then
           cnk_all_fbz(1,ib,ik,is,it)=cnk_tem(1)
        end if
     end do
     end if
  end do
  end do
  end do

  end subroutine rotate_cnk
