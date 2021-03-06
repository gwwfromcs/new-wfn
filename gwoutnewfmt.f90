! Weiwei Gao: adapted from gwout.f90
! Assume that bands%num_gwout == -1, which means output all bands
! further improvements is needed


!
subroutine gwoutnewfmt (fname, crys, syms, bands, gs, kgs, wfn, kpoints, &
     iflag, pw)
  !
  use all_to_all_module  
  include 'use.h'
  implicit none                    ! implicit? Just say no!
  include 'interface.h'  
  include 'all_to_all.h'
  !  
  !     include 'mpif.h'
  !
  !     writes a data array to disk such that it can be read and used in
  !     later GW calculations
  !
  !     INPUT:
  !     -----
  !
  type(crystal), intent(in) :: crys  
  type(symmetry), intent(in) :: syms  
  type(complex_gspace_array), intent(in) :: &
       wfn                                ! the wave functions
  type(band), intent(in) :: bands  
  type(parallel_gspace), intent(in) :: &
       kgs(bands%nrk), &                  ! the gspaces at the different kpoints
       gs                                 ! the big gspaces ---> pot_gspace
  type(kpoint), intent(in) :: kpoints     ! all the kpoints
  integer, intent(in) :: iflag      ! iflag =1 real version, iflag =2, complex
  character(len=*), intent(in) :: fname   ! file name
  !     ---------------- local variables ---------------------------------
  integer :: i, info, iproc, irk, is, j, k, n, ndim, ng
  integer :: idum(1),n_bands_out
  integer, allocatable :: ifmin(:,:), ifmax(:,:), nkpt(:)    
  complex(dp), allocatable :: zc(:,:)

! Weiwei Gao
  !     ---------------- wfnselect.inp -----------------
  ! selective output wavefunctions, for speed-up GW calculation with energy
  ! integrations
!WG: parameters for energy integration
  type(pw_parameter), intent(in) :: pw
  real(DP) :: startingenergy1, startingenergy2
  real(DP) :: energystep1, energystep2, highestenergy
  real(DP) :: Vxc0, cellvol, e1, e2, total_w
  integer  :: nsp, nsp1, nsp2, index_startEI
  integer,allocatable  :: bandls(:,:)
  integer  :: bandflag, acc, ncount
  integer  :: write_flag, weight_flag
  integer  :: irecord
  integer,allocatable :: allgvec(:,:), allkgvec(:,:)
  integer  :: ngcount, nkgcount

!WG: following definitions are copied from
!WG: BGW-1.0.4/Common/wfn_rho_vxc_io.p.f/write_binary_header
  integer  :: ns, ntran, cell_symmetry, nat, nk, nbands, ngkmax !< numbers of spins, G-vectors, symmetries,
                             !! cell type (0 = cubic, 1 = hexagonal), numbers of atoms, k-points, bands, max(ngk)
  real(DP) :: ecutrho, ecutwfc !< charge-density and wave-function cutoffs, in Ry
  real(DP) :: celvol, alat, avec(3, 3), adot(3, 3) !< cell volume, lattice constant,
                              !! lattice vectors, metric tensor in real space (in a.u., avec in units of alat)
  real(DP) :: recvol, blat, bvec(3, 3), bdot(3, 3) !< cell volume, lattice constant,
                              !! lattice vectors, metric tensor in reciprocal space (in a.u., bvec in units of blat)
  character(len=32) :: sdate, stime !< if read, result from file is returned; if write, current is returned
  character :: sheader*3, sdate_*32, stime_*32, stitle*32, sflavor*7
  integer :: ii, jj, ib, igg, ik, itran, iat, ierr
  character :: adate*11, atime*14
  real(DP), allocatable :: apos(:,:), dummypos(:)
  integer, allocatable :: atype(:) !atomic species
  character :: in_cell_symmetry*10
  integer :: in_nat 
  integer :: kgrid(3)
!WG

!------------------------------------------------------------
! switch band index, modify this part if needed
! switch the second and the fifth band at the Gamma point 
!
! PZ
!
! bands%gwout_index(n,irk,is)

!  if (bands%gw_band_reordering.eq.1) then
!      bands%gwout_index(2,1,:)=5
!      bands%gwout_index(5,1,:)=2
!  end if
!------------------------------------------------------

  if (gs%myproc == 0) then
!WG: open and read wfnselect.inp
    if (pw%gw_sum == 1) then
      open(unit=400078, file = 'wfnselect.inp', &
           status = 'old', form = 'formatted') 
      read(400078,*) index_startEI
      read(400078,*) energystep1
      read(400078,*) startingenergy2
      read(400078,*) energystep2
      read(400078,*) weight_flag
      read(400078,*) Vxc0
      read(400078,*) cellvol
      read(400078,*) in_cell_symmetry
      if (trim(in_cell_symmetry).eq.'hexagonal') then
        cell_symmetry = 1
      else if(trim(in_cell_symmetry).eq.'cubic') then
        cell_symmetry = 0
      else
        write(*,*) 'wrong cell_symmetry in wfnselect.inp, set to default(cubic).'
        cell_symmetry = 0
      endif
      read(400078,*) in_nat
      nat = sum(crys%natom)
      if (nat.ne.in_nat) then
          write(*,*) 'nat=',nat,',in_nat=',in_nat
          !call mystop('nat is not equal to in_nat')
          write(*,*) 'nat is not equal to in_nat'
          write(*,*) 
      endif
      allocate(apos(3,nat))
      allocate(atype(nat))
      do iat = 1,nat
        read(400078,*) (apos(ii,iat),ii=1,3),atype(iat) 
        !atom positions in bohr!!! and atom type
      enddo
      read(400078,*) (kgrid(ii),ii=1,3)
      !WG: when there is gw shift, kpoints%grid is 0 0 0, so I need to
      !specify kgrid in wfnselect.inp
      close(400078)
      
      !WG: now, we calculate startingenergy1
      startingenergy1 = 0
      highestenergy=0
      total_w = 0
      do ik = 1,kpoints%nrk
       do is = 1,wfn%nspin
       ! startingenergy1 = startingenergy1 + kp%w(ik) * &
       !      kp%el(index_startEI,ik,is)   !in Ry
        startingenergy1 = startingenergy1 + kpoints%w(ik) * &
             (bands%energy(index_startEI,ik,is)+Vxc0)
        highestenergy   = highestenergy + kpoints%w(ik) * &
            (bands%energy(bands%min(1),ik,is)+Vxc0)       !in Ry
        write(6,*) 'ik=',ik,' highest band energy', &
        bands%energy(bands%min(1),ik,is)
        total_w = total_w + kpoints%w(ik)
       enddo
       write(6,*) 'ik=',ik,' weight is', kpoints%w(ik)
      enddo
      startingenergy1 = startingenergy1 / total_w
      highestenergy = highestenergy / total_w
!WG: 
      write(6,*) 'total k points weight', total_w
      write(6,*) 'average energy of band ',index_startEI,' = ',startingenergy1
      !write(6,*) 'average energy of highest band ',kp%mnband,' = ',highestenergy
      write(6,*) 'average energy of highest band ',bands%min(1),' = ',highestenergy
!WG:  recalculating startingenergy1 and highest energy
      startingenergy1 = Vxc0 + ( (6.0*pi**2 * index_startEI)/cellvol )**(2.0/3.0)
!WG: 
      if ((highestenergy-Vxc0)**1.5 * cellvol / (6 * pi**2) .gt. bands%min(1)-1) then
        highestenergy = Vxc0 + ( (6.0*pi**2 * bands%min(1)) / cellvol )**(2.0/3.0)
      endif
      
!WG:  recalculate startingenergy2
      startingenergy2 = startingenergy1 + & 
        int((startingenergy2-startingenergy1)/energystep2) * energystep2
      
!WG:  print out what we got from wfnselect.inp
      write(6,*) 'recalculated startingenergy1',startingenergy1
      write(6,*) 'startingband for energy integration:',index_startEI
      write(6,*) 'energystep1:',energystep1
      write(6,*) 'recalculated startingenergy2:',startingenergy2
      write(6,*) 'energystep1:',energystep2
      write(6,*) 'recalculated highest energy:',highestenergy
      
!WG:  calculating the weight of each energy step
      nsp1=int( (startingenergy2-startingenergy1)/energystep1 )
      nsp2=int( (highestenergy-startingenergy2)/energystep2 )
        
      nsp = index_startEI + nsp1 + nsp2
      ALLOCATE(bandls(2, nsp))
      total_w=0
      do ib = 1,index_startEI
        bandls(2,ib) = 1
        bandls(1,ib) = ib
        total_w = total_w + bandls(2,ib)
      enddo
      do ib = 1,nsp1
        e1 = startingenergy1 + (ib-1)*energystep1
        e2 = startingenergy1 + (ib)  *energystep1
        bandls(2,ib+index_startEI) = nint( &
         ((e2 - Vxc0)**1.5 * cellvol / ( 6 * pi**2 )) - &
         ((e1 - Vxc0)**1.5 * cellvol / ( 6 * pi**2 )) )
        bandls(1,ib+index_startEI) = &
         (e2 - Vxc0)**1.5 * cellvol / ( 6 * pi**2 ) - &
         bandls(2,ib+index_startEI)/2 
        total_w = total_w + bandls(2,ib+index_startEI)
      enddo
      do ib = 1,nsp2
        e1 = startingenergy2 + (ib-1)*energystep2
        if (ib.lt.nsp2) then
           e2 = startingenergy2 + (ib)  *energystep2
        else
           e2 = highestenergy
        endif
        bandls(2,index_startEI+nsp1+ib) = nint( &
         ((e2 - Vxc0)**1.5 * cellvol / ( 6 * pi**2 )) - &
         ((e1 - Vxc0)**1.5 * cellvol / ( 6 * pi**2 )) )
        bandls(1,index_startEI+nsp1+ib) = &
         (e2 - Vxc0)**1.5 * cellvol / ( 6 * pi**2 ) - &
         bandls(2,ib+index_startEI+nsp1)/2 
        total_w = total_w + bandls(2,ib+index_startEI+nsp1)
      enddo
      print *, 'total_w=',total_w
      
!WG:  calculating the indices of sampled band 
!WG:  printout the band list
      write(6,'( "List the indices of bands we sampled" )') 
      write(6,'( " number     real_band_index        weight      energy ")')
      do ib = 1,nsp
        write(6,'(3(1x,i12),1x,f18.9)') ib, bandls(1,ib), bandls(2,ib), &
          bands%energy(bandls(1,ib),1,1)+Vxc0
      enddo 
!WG:  write the band list to weight.dat
      open(unit=400079,file='weight.dat',form='formatted',status='replace')
      write(400079,'(4(1x,f10.5))') &
        startingenergy1,energystep1,startingenergy2,energystep2
      write(400079,'(2(1x,i10))') nsp,index_startEI
      do ib = 1,nsp
        write(400079,'(3(1x,i12),1x,f18.9,f18.9)') ib, bandls(1,ib), bandls(2,ib), &
          bands%energy(bandls(1,ib),1,1)+Vxc0
      enddo
      close(400079)
    else !pw%gw_sum = 2
!WG:  if pw%gw_sum = 2, then write all the wavefunctions in newformat
!WG:  write the band list to weight.dat
!WG:  even though we don't do band sampling here, we still need some info for
!WG:  new format wavefunctions
      open(unit=400078, file = 'wfnselect.inp', &
           status = 'old', form = 'formatted') 
      read(400078,*) ! index_startEI
      read(400078,*) ! energystep1
      read(400078,*) ! startingenergy2
      read(400078,*) ! energystep2
      read(400078,*) ! weight_flag
      read(400078,*) Vxc0
      read(400078,*) cellvol
      read(400078,*) in_cell_symmetry
      if (trim(in_cell_symmetry).eq.'hexagonal') then
        cell_symmetry = 1
      else if(trim(in_cell_symmetry).eq.'cubic') then
        cell_symmetry = 0
      else
        write(*,*) 'wrong cell_symmetry in wfnselect.inp, set to default(cubic).'
        cell_symmetry = 0
      endif
      read(400078,*) in_nat
      nat = sum(crys%natom)
      if (nat.ne.in_nat) then
          write(*,*) 'nat=',nat,',in_nat=',in_nat
          !call mystop('nat is not equal to in_nat')
          write(*,*) 'nat is not equal to in_nat'
          write(*,*) 
      endif
      allocate(apos(3,nat))
      allocate(atype(nat))
      do iat = 1,nat
        read(400078,*) (apos(ii,iat),ii=1,3),atype(iat) 
        !atom positions in bohr!!! and atom type
      enddo
      read(400078,*) (kgrid(ii),ii=1,3)
      close(400078)

      nsp = bands%min(1)
      Allocate(bandls(2,nsp))
      do ib = 1,nsp
        bandls(2,ib) = 1 !weight of the sampled band
        bandls(1,ib) = ib 
      enddo
      open(unit=400079,file='weight.dat',form='formatted',status='replace')
      write(400079,*) 'full band sampling'
      do ib = 1,nsp
        write(400079,'(3(1x,i12),1x,f18.9,f18.9)') ib, bandls(1,ib), bandls(2,ib), &
          bands%energy(bandls(1,ib),1,1)+Vxc0
      enddo
      close(400079)
    endif
  endif !if gs%myproc

!WG: calculate ng at each k points
  allocate(nkpt(kpoints%nrk))
  idum(1) = gs%length  
  call all_sum_all_integer(idum, 1)  
  ng = idum(1)
  ngkmax = 0
  do irk = 1, kpoints%nrk  
     ndim = kgs(irk)%length  
     idum(1) = ndim  
     call all_sum_all_integer(idum, 1)  
     nkpt(irk) = idum(1)
     if (nkpt(irk) .gt. ngkmax) then
       ngkmax = nkpt(irk)
     endif
  enddo

! WG: Headers!!!
  if (gs%myproc == 0) then  
     !
     !     erase old file
     !
     write(*,*) 'open file ',fname
     open(21, file = fname, status = 'replace', form = 'unformatted')

     !------old code------
     !write(21) ((crys%bdot(i, j), i = 1, 3), j = 1, 3)  
     !write(21) crys%vcell  
     !
     !write(21) syms%ntrans  
     !do n = 1, syms%ntrans  
     !   write(21) ((syms%mtrx(i, j, n), i = 1, 3), j = 1, 3)  
     !end do
     !do n = 1, syms%ntrans  
     !   write(21) (syms%tnp(k, n), k = 1, 3)  
     !end do
     !
     !write(21) crys%nspin  
     !write(21) kpoints%nrk  
     !write(21) (kpoints%grid(i), i = 1, 3)  
     !write(21) (kpoints%shift(i), i = 1, 3)  
     !write(21) (kpoints%w(i), i = 1, kpoints%nrk)  
     !do j = 1, kpoints%nrk  
     !   write(21) (kpoints%rk(i, j), i = 1, 3)  
     !end do
     !
     !if (bands%num_gwout .eq. -1) then     
     !  write(21) bands%min(1)  
     !  n_bands_out=bands%min(1)
     !else
     !  write(21) 2*bands%num_gwout       
     !  n_bands_out=2*bands%num_gwout
     !end if
     !------old code-------

     !WG: time info in new format
     call date_time(adate, atime)
     sheader = 'WFN'
     if (iflag == 1) then
         sflavor = "Real"
     else 
         sflavor = "Complex"
     endif
     sdate_ = adate
     stime_ = atime
     stitle = sheader // "-" // sflavor
     write(*,*) 'write date and time to ',fname
     write(21, iostat=ierr) stitle, sdate_, stime_
     
     
     write(*,*) 'start to calculate some data'
     !get some data
     ns = crys%nspin
     ntran = syms%ntrans
     cell_symmetry = 0 
     !nat = sum(crys%natom)
     ecutrho = gs%gmax*gs%gmax
     nk = kpoints%nrk
     !nbands = bands%min(1)
     nbands = nsp
     ecutwfc = pw%emax
     write(*,*) 'print out for debug'
     write(*,*) ns, ng, ntran, cell_symmetry, nat, ecutrho, nk, nbands, ngkmax, ecutwfc
     !
     write(21) ns, ng, ntran, cell_symmetry, nat, ecutrho, nk, nbands, ngkmax, ecutwfc
     write(21) (gs%fftsize(ii), ii = 1, 3), (kgrid(ii), ii = 1, 3), & 
               (kpoints%shift(ii), ii = 1, 3) 
               !WG: when there is gw shift, kpoints%grid is 0 0 0, so I need to
               !specify kgrid in wfnselect.inp
                        

     write(*,*) (gs%fftsize(ii), ii = 1, 3), (kgrid(ii), ii = 1, 3), & 
               (kpoints%shift(ii), ii = 1, 3)
     
     !allocate(apos(3,nat))
     allocate(dummypos(3))
     !convert the atomic postion from crystal coord to cartesian(bohr)
     do iat = 1, nat
       dummypos(1:3) = apos(1:3, iat)
       apos(1:3, iat) = 0
       do ii = 1, 3
         write(*,*) (crys%avec(jj,ii),jj=1,3)
         do jj = 1, 3
           apos(jj, iat) = apos(jj, iat) + &
             dummypos(ii) * crys%avec(jj, ii)
         enddo
       enddo
     enddo

     !calculate some crystal info from old wavefunctions
     alat = sqrt(crys%avec(1, 1)**2 + crys%avec(2, 1)**2 + crys%avec(3,1)**2)
     avec(1:3, 1:3) = crys%avec(1:3, 1:3) / alat
     apos(1:3, 1:nat) = apos(1:3, 1:nat) / alat
     celvol = avec(1, 1) * (avec(2, 2) * avec(3, 3) - avec(2, 3) * avec(3, 2)) - &
              avec(2, 1) * (avec(1, 2) * avec(3, 3) - avec(1, 3) * avec(3, 2)) + &
              avec(3, 1) * (avec(1, 2) * avec(2, 3) - avec(1, 3) * avec(2, 2))

     blat = 2.0d0 * pi / alat
     bvec(1,1) = (avec(2,2) * avec(3,3) - avec(3,2) * avec(2,3)) / celvol
     bvec(2,1) = (avec(3,2) * avec(1,3) - avec(1,2) * avec(3,3)) / celvol
     bvec(3,1) = (avec(1,2) * avec(2,3) - avec(2,2) * avec(1,3)) / celvol
     bvec(1,2) = (avec(2,3) * avec(3,1) - avec(3,3) * avec(2,1)) / celvol
     bvec(2,2) = (avec(3,3) * avec(1,1) - avec(1,3) * avec(3,1)) / celvol
     bvec(3,2) = (avec(1,3) * avec(2,1) - avec(2,3) * avec(1,1)) / celvol
     bvec(1,3) = (avec(2,1) * avec(3,2) - avec(3,1) * avec(2,2)) / celvol
     bvec(2,3) = (avec(3,1) * avec(1,2) - avec(1,1) * avec(3,2)) / celvol
     bvec(3,3) = (avec(1,1) * avec(2,2) - avec(2,1) * avec(1,2)) / celvol
   
     celvol = abs(celvol) * alat**3
     recvol = (2.0d0 * pi)**3 / celvol

     write(*,*) 'finished calculating some basic parameters'
   
     do ii=1,3
       do jj=1,3
         adot(jj,ii) = dot_product(avec(1:3,jj), avec(1:3,ii)) * alat**2
       enddo
     enddo
   
     do ii=1,3
       do jj=1,3
         bdot(jj,ii) = dot_product(bvec(1:3,jj), bvec(1:3,ii)) * blat**2
       enddo
     enddo

     write(21) celvol, alat, ((avec(jj, ii), jj = 1, 3), ii = 1, 3), &
     ((adot(jj, ii), jj = 1, 3), ii = 1, 3)
     write(21) recvol, blat, ((bvec(jj, ii), jj = 1, 3), ii = 1, 3), &
     ((bdot(jj, ii), jj = 1, 3), ii = 1, 3)
     write(21) (((syms%mtrx(ii, jj, itran), ii = 1, 3), jj = 1, 3), itran = 1, syms%ntrans)
     write(21) ((syms%tnp(jj, itran), jj = 1, 3), itran = 1, syms%ntrans)
     write(21) ((apos(ii, iat), ii = 1, 3), atype(iat), iat = 1, nat)

     write(*,*) celvol, alat, ((avec(jj, ii), jj = 1, 3), ii = 1, 3), &
     ((adot(jj, ii), jj = 1, 3), ii = 1, 3)
     write(*,*) recvol, blat, ((bvec(jj, ii), jj = 1, 3), ii = 1, 3), &
     ((bdot(jj, ii), jj = 1, 3), ii = 1, 3)
     write(*,*) (((syms%mtrx(ii, jj, itran), ii = 1, 3), jj = 1, 3), itran = 1, syms%ntrans)
     write(*,*) ((syms%tnp(jj, itran), jj = 1, 3), itran = 1, syms%ntrans)
     write(*,*) ((apos(ii, iat), ii = 1, 3), atype(iat), iat = 1, nat)
     !write(21) (ngk(ik), ik = 1, nk)
     !write(21) (kw(ik), ik = 1, nk)
     !write(21) ((kpt(ii, ik), ii = 1, 3), ik = 1, nk)
     !write(21) ((ifmin(ik, is), ik = 1, nk), is = 1, ns)
     !write(21) ((ifmax(ik, is), ik = 1, nk), is = 1, ns)
     !write(21) (((energies(ib, ik, is), ib = 1, nbands), ik = 1,nk), is = 1, ns)
     !write(21) (((occupations(ib, ik, is), ib = 1, nbands), ik = 1, nk), is = 1, ns)

     write(21) (nkpt(irk), irk = 1,kpoints%nrk)
     write(21) (kpoints%w(irk), irk = 1,kpoints%nrk)
     write(21) ((kpoints%rk(i, irk), i = 1, 3), irk = 1, kpoints%nrk)

     write(*,*) 'starting determine ifmin and ifmax'
     allocate(ifmin(bands%nrk, crys%nspin))  
     allocate(ifmax(bands%nrk, crys%nspin))  
     ifmin = 0  
     ifmax = 0  
     do is = 1, wfn%nspin  
        do irk = 1, bands%nrk  
           do n = 1, bands%nband(irk, is)
!
! SIB:  1/1/2001
! changed occupancy criterion...bands%occup is actually the band occupancy
! times 2*weight, so for many k-points, it can get very small!
! Original code:
!              if(bands%occup(n,irk,is).gt.1.0d-3) then
!
              if(bands%occup(n,irk,is)/(2.0d0*kpoints%w(irk)) &
                                         .gt.1.0d-3) then
                 if (ifmin(irk, is) == 0) ifmin(irk, is) = n  
                 ifmax(irk, is) = n  
              end if
           end do
        end do
     end do
     write(21) ((ifmin(irk, is), irk = 1, bands%nrk), is = 1, crys%nspin)
     write(21) ((ifmax(irk, is), irk = 1, bands%nrk), is = 1, crys%nspin)
     write(*,*) ((ifmin(irk, is), irk = 1, bands%nrk), is = 1, crys%nspin)
     write(*,*) ((ifmax(irk, is), irk = 1, bands%nrk), is = 1, crys%nspin)
!     do is = 1, wfn%nspin  
!        do irk = 1, bands%nrk  
!           write(21) (bands%energy(n, irk, is), n = 1, bands%min(is))  
!           write(21) (bands%energy( bands%gwout_index(n,irk,is),&
!                      irk, is), n = 1,n_bands_out )  
!        end do
!     end do
     write(21) (((bands%energy(bandls(1,ib), irk, is),ib = 1,nsp), &
             irk = 1,kpoints%nrk), is = 1, wfn%nspin)
     write(21) (((bands%occup(bandls(1,ib), irk, is)/(2.0d0*kpoints%w(irk)),ib = 1,nsp), &
             irk = 1,kpoints%nrk), is = 1, wfn%nspin)
!     write(*,*) (((bands%energy(ib, irk, is),ib = 1,bands%min(is)), &
!             irk = 1,kpoints%nrk), is = 1, wfn%nspin)
!     write(*,*) (((bands%occup(ib, irk, is)/(2.0d0*kpoints%w(irk)),ib = 1,bands%min(is)), &
!             irk = 1,kpoints%nrk), is = 1, wfn%nspin)

     close(21)  
     write(*,*) 'finish writing header'
  end if
  call my_broadcast(n_bands_out,0)
  !

  !------old code------
  !idum(1) = gs%length  
  !call all_sum_all_integer(idum, 1)  
  !ng = idum(1)  
  !if (gs%myproc == 0) then  
  !   open(unit = 21, file = fname, position = 'append', status = &
  !        'unknown', form = 'unformatted')
  !   write(21) (gs%fftsize(i), i = 1, 3)  
  !   write(21) ng  
  !   close(21)  
  !end if
  !------old code------

  if (gs%myproc == 0) then
     open(unit = 21, file = fname, position = 'append', status = &
               'unknown', form = 'unformatted')
     write(21) 1
     write(21) ng
     close(21)
  end if

  allocate(allgvec(3,ng))
  allgvec = 0
  ngcount = 0
  do iproc = 0, gs%nproc - 1  
     if (iproc == gs%myproc) then  
        do i = 1, gs%length  
          !write(21) gs%gvec(1, i), gs%gvec(2, i), gs%gvec(3, i)  
          allgvec(1:3,ngcount+i) = gs%gvec(1:3,i)
        end do
        !open(unit = 21, file = fname, position = 'append', status = &
        !       'unknown', form = 'unformatted')
        !write(21) ((gs%gvec(ii,igg),ii=1,3),igg=1,gs%length)
        ngcount = ngcount + gs%length
        !close(21)
     end if
     call my_broadcast(ngcount,iproc)
     if (gs%nproc > 1) call parallel_barrier()  
  end do
  !call mpi_allreduce(allgvec, allgvec, 3*ngcount, MPI_INTEGER, &
  !   MPI_SUM, MPI_COMM_WORLD, info)
  call all_sum_all_integer(allgvec(1,1),3*ngcount)

  if (gs%myproc == 0) then
    write(*,*) 'write the big gvector'
  !write all gvectors for the big gspace
    open(unit = 21, file = fname, position = 'append', status = &
               'old', form = 'unformatted')
    write(21) ((allgvec(ii, igg), ii = 1, 3), igg = 1, ng)
    close(21)
  endif
  deallocate(allgvec)

  !
  if (iflag == 1) call gwreal(bands, kgs, wfn, kpoints) !  
  !
  do irk = 1, kpoints%nrk  
     ndim = kgs(irk)%length
     if (gs%myproc == 0) then  
        write(*,*) 'write nkpt for irk ',irk,' nkpt(irk):',nkpt(irk)
        open(unit = 21, file = fname, position = 'append', status = &
             'unknown', form = 'unformatted')
        write(21) 1
        write(21) nkpt(irk)  
        close(21)  
     end if

     !---old code here---
!    do iproc = 0, gs%nproc - 1  
!       if (iproc == gs%myproc) then  
!          open(unit = 21, file = fname, position = 'append', status = &
!               'old', form = 'unformatted')
!          do i = 1, kgs(irk)%length  
!             write(21) kgs(irk)%gvec(1, i), kgs(irk)%gvec(2, i), &
!                  kgs(irk)%gvec(3, i)
!          end do
!          close(21)  
!       end if
!       if (gs%nproc > 1) call parallel_barrier()  
!    end do
     !---old code end---

     allocate(allkgvec(3,nkpt(irk)))
     allkgvec = 0
     nkgcount = 0
     do iproc = 0, gs%nproc - 1  
        if (iproc == gs%myproc) then  
           do i = 1, kgs(irk)%length  
             !write(21) gs%gvec(1, i), gs%gvec(2, i), gs%gvec(3, i)  
             allkgvec(1:3,nkgcount+i) = kgs(irk)%gvec(1:3,i)
           end do
           !open(unit = 21, file = fname, position = 'append', status = &
           !  'unknown', form = 'unformatted')
           !write(21) ((kgs(irk)%gvec(ii,igg),ii=1,3),igg=1,kgs(irk)%length)
           nkgcount = nkgcount + kgs(irk)%length
           !close(21)
        end if
        call my_broadcast(nkgcount,iproc)
        if (gs%nproc > 1) call parallel_barrier()  
     end do
     if(nkpt(irk).ne.nkgcount) then
         write(*,*) 'nkpt=',nkpt,',nkgcount=',nkgcount
         call mystop('nkpt not equal to nkgcount')
     endif
     !call mpi_allreduce(allkgvec, allkgvec, 3*nkgcount, MPI_INTEGER, &
     !    MPI_SUM, MPI_COMM_WORLD, info)
     call all_sum_all_integer(allkgvec,3*nkgcount)

     if (gs%myproc == 0) then
     ! write all gvectors for the big gspace
       write(*,*) 'write kgvec space.'
       open(unit = 21, file = fname, position = 'append', status = &
                  'old', form = 'unformatted')
       write(21) ((allkgvec(ii, igg), ii = 1, 3), igg = 1, nkpt(irk))
       close(21)
     endif
     deallocate(allkgvec)

     if (gs%myproc == 0) open(unit = 21, file = fname, position = &
          'append', status = 'unknown', form = 'unformatted')
     allocate(zc(nkpt(irk), wfn%nspin))  
!     do n = 1, bands%min(1)  
     ncount = 1

     write(*,*) 'before band loop myproc=',gs%myproc,',bands%min(1)=',bands%min(1)
     call my_broadcast(nsp,0)
     !call my_broadcast_int(bandls,2*nsp,0)
     do n = 1, bands%min(1)
        do is = 1, wfn%nspin  
!           call data_gather(1, gs%myproc, gs%nproc, &
!                wfn%data(1 + (n - 1) * ndim, irk, is), zc(1, is), ndim)
           call data_gather(1, gs%myproc, gs%nproc, &
                wfn%data(1 + (bands%gwout_index(n,irk,is) - 1) * ndim, irk, is), zc(1, is), ndim)
        end do

        if (iflag == 1) then  
           if (gs%myproc == 0) then
             if(bandls(1,ncount).eq.n) then
               write(*,*) 'start to write wavefunction for band',n
               write(21) 1
               write(21) nkpt(irk)
               write(21) ((real(zc(j, k), dp), j = 1, &
                  nkpt(irk)), k = 1, wfn%nspin)
  !            write(*,*) ((real(zc(j, k), dp), j = 1, &
  !               nkpt(irk)), k = 1, wfn%nspin)
               ncount = ncount + 1
             endif
           endif
        else  
           if (gs%myproc == 0) then
             if(bandls(1,ncount).eq.n) then
               write(*,*) 'start to write wavefunction for band',n
               write(21) 1
               write(21) nkpt(irk)
               write(21) ((zc(j, k), j = 1, nkpt(irk)), &
                  k = 1, wfn%nspin)
               ncount = ncount + 1
             endif
           endif
        end if

        call my_broadcast(ncount,0)      
        if (ncount.gt.nsp) then
          exit
        endif
          
     end do
     deallocate(zc)  
     if (gs%myproc == 0) close(21)  
     write(*,*) 'irk=',irk,',gs%myproc=',gs%myproc
  end do

  if (gs%myproc == 0) then  
    deallocate(ifmin)
    deallocate(ifmax)
  end if

  return  

end subroutine gwoutnewfmt

  subroutine date_time(bdate,btime)
    character, intent(out) :: bdate*11,btime*14
    
    integer :: lmonth
    integer :: idate (8)
    character :: day*2,year*4
    character :: adate*8,atime*10,azone*5
    character :: hour*2,min*2,sec*2
    character*3 :: month(12)
      
    DATA month/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep', &
      'Oct','Nov','Dec'/


    call date_and_time(adate,atime,azone,idate)
    
    read(adate,101) year,lmonth,day
101 format(a4,i2,a2)
    write(bdate,102) day,month(lmonth),year
102 format(a2,'-',a3,'-',a4)
    read(atime,201) hour,min,sec
201 format(a2,a2,a2,a4)
    write(btime,202) hour,min,sec,azone
202 format(a2,':',a2,':',a2,1x,a5)

    return    
  end subroutine date_time
