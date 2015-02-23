! m4undef.m4
!
! resets various m4 commands to have a prefix m4_
! this can be achieved using option -P but only on some architectures
! any file to be preprocessed with m4 should use this set of macro definitions
!
! David Prendergast, June 6, 2006



! fft_macros.m4
!
!     fft_aux_space(aux,naux,aux2,naux2,dec1,dec2)
!     fft_local_init(N,A,LDA,LOT,DIREC)
!     fft_perm_init(N,AUX)
!
!     fft_multiple_backward(N,A,LOTSTRIDE,LOT,LOOPDUMMY,VECOUT,OFFSET)
!     fft_multiple_forward (N,A,LOTSTRIDE,LOT,LOOPDUMMY,VECOUT,OFFSET)
!
!     fft_backward(N,A,LOTSTRIDE)
!     fft_forward (N,A,LOTSTRIDE)
!
!     fft_convol(N,A,B,POS)
!
!     fft_multiple_scale(N,A,LOTSTRIDE,LOT,LOOPDUMMY,SCALE)
!
!     timeget (t0)
!     timediff(t1,t0)
!
!
!     fft_fcblock     defines blocking factor for fast convolute
!     fft_local_free  release memory
!

!-*-F90-*-
!
!     Fortran 90 structures.
!
!     (1996) Bernd Pfrommer
!
!     Just include the modules you need
!

!
!     ========== parallel gspace structure ==========================
!

module parallel_gspace_module  

  use constants
  implicit none  
  !
  !     first some public data which is global to ALL GSPACES
  !
  ! the blacs context (processor grid layout)

  integer :: blacs_context

  type parallel_gspace  
     sequence
     !
     !     ------ true information --------------
     !
     real(dp) :: gmax             ! cutoff length
     real(dp) :: rk(3)            ! the k-vector by which it is shifted
     integer :: length            ! number of g-space vectors on this proc
     integer :: totlength         ! total number of g-space vectors (all procs)
     integer :: nstar             ! total number of stars in gspace
     integer :: ig0               ! index of the g=0. If not on this proc, = -1
     integer :: ntype             ! number of types of atoms
     logical :: igvec             ! if gvectors are available
     logical :: iekin             ! if kinetic energy is available
     logical :: iinverse          ! if inverse index array is available
     logical :: istruc            ! if structure factors are there
     logical :: istar             ! if star arrays are available
     integer, pointer :: &
          gvec(:,:)         ! the gspace vectors themselves
     integer, pointer :: &
          inverse(:)        ! index to inverse of G. only for local.
     real(dp), pointer :: &
          ekin(:)           ! kinetic energy |k+G|^2
     complex(dp), pointer :: &
          struc(:,:)        ! struc factors struc(ntype, length)
     integer, pointer :: inds(:)   ! maps gvecs to stars inds(length)
     integer, pointer :: mstar(:)    ! size of stars mstar(nstar)
     complex(dp), pointer :: &
          phase(:)          ! phases(length)
     !
     !     ------ parallel information ----------
     !

     integer, pointer :: &
          order(:,:)        ! orders the gspace layout order(3, lorder)
     integer, pointer :: &
          fastfind(:,:)     ! help array for fastfind(3, fftsize(1))
     integer :: lorder            ! length of order array
     integer :: nproc             ! number of procs for gspace
     integer :: myproc            ! myproc for this gspace
     integer :: maxlength_all_k   ! max g-vector length of all k points

     !
     !     ------ FFT information ---------------
     !
     integer :: fftsize(5)        ! FFT grid boundaries for gspace
     integer :: naux(4)           ! the size of the four workspaces
     integer :: r_size            ! the size of the complementary realspace
     integer :: maxpacksize_t1    ! sizeof 1st dim of packinfo (1st transpose)
     integer :: maxpacksize_t2    ! maximum packet size for second transpose
     integer :: totpacksize_t2    ! total size of all packets sent during t2
     integer :: maxnchunk         ! dimension of 2nd index of chunk array
     integer :: ipad(1)  
     integer, pointer :: &
          packinfo(:,:,:)   ! packing/unpacking information
     integer, pointer :: &
          chunk(:,:,:,:)    ! chunk(4, pack/unpack, maxnchunk, nprocs)
     integer, pointer :: &
          packsize(:,:,:)   ! see setup_packinfo routine
     integer, pointer :: &
          pickup(:,:)       ! pickup(nproc, inverse/forward)
     integer, pointer :: &
          packtype(:)       ! mpi data type array
     character(len=16) :: name    ! a string identifying the gspace
     logical :: ipackinfo         ! if packinfo is there or not
     integer :: ipad2(1)  

  end type parallel_gspace

contains  

  subroutine destroy_parallel_gspace(gs)  

    implicit none
    type(parallel_gspace), intent(inout) :: gs  

    if (associated(gs%order)) deallocate(gs%order)  
    if (associated(gs%gvec)) deallocate(gs%gvec)  
    if (associated(gs%ekin)) deallocate(gs%ekin)  
    if (associated(gs%inds)) deallocate(gs%inds)  
    if (associated(gs%mstar)) deallocate(gs%mstar)  
    if (associated(gs%phase)) deallocate(gs%phase)  
    if (associated(gs%inverse)) deallocate(gs%inverse)  
    if (associated(gs%struc)) deallocate(gs%struc)  
    if (associated(gs%packinfo)) deallocate(gs%packinfo)  
    if (associated(gs%packsize)) deallocate(gs%packsize)  
    if (associated(gs%fastfind)) deallocate(gs%fastfind)  
    if (gs%ipackinfo .and. associated(gs%packsize)) deallocate(gs%chunk)
    if (gs%ipackinfo .and. associated(gs%pickup)) deallocate(gs%pickup)
    if (gs%ipackinfo .and. associated(gs%packtype)) deallocate(gs%packtype)

  endsubroutine destroy_parallel_gspace
  !
  !     assign two gspaces to each other
  !
  subroutine parallel_gspace_assignment(gt, gs)

    implicit none 
    type(parallel_gspace), intent(out) :: gt  
    type(parallel_gspace), intent(in) :: gs

    gt%gmax = gs%gmax  
    gt%rk = gs%rk  
    gt%length = gs%length  
    gt%totlength = gs%totlength  
    gt%nstar = gs%nstar  
    gt%igvec = gs%igvec  
    if (gs%igvec .and. associated(gs%gvec)) then  
       allocate(gt%gvec(3, gt%length))  
       gt%gvec = gs%gvec  
    end if
    gt%iekin = gs%iekin
    if (gs%iekin .and. associated(gs%ekin)) then  
       allocate(gt%ekin(gt%length))  
       gt%ekin(:) = gs%ekin(:)  
    end if
    gt%iinverse = gs%iinverse  
    if (gs%iinverse .and. associated(gs%inverse)) then  
       allocate(gt%inverse(gt%length))  
       gt%inverse(:) = gs%inverse(:)  
    end if
    gt%istruc = gs%istruc  
    gt%ntype = gs%ntype  
    if (gs%istruc .and. associated(gs%struc)) then  
       allocate(gt%struc(gt%ntype, gt%length))  
       gt%struc(:,:) = gs%struc(:,:)  
    end if
    gt%istar = gs%istar  
    if (gs%istar) then  
       allocate(gt%inds(gt%length))  
       gt%inds(:) = gs%inds(:)  
       allocate(gt%mstar(gt%nstar))  
       gt%mstar(1:gs%nstar) = gs%mstar(1:gs%nstar)  
       allocate(gt%phase(gt%length))  
       gt%phase(:) = gs%phase(:)  
    end if
    gt%lorder = gs%lorder  
    allocate(gt%order(3, gt%lorder))  
    gt%order(:,:) = gs%order(:,:)  
    gt%fastfind(:,:) = gt%fastfind(:,:)  
    gt%nproc = gs%nproc  
    gt%myproc = gs%myproc  
    gt%fftsize = gs%fftsize  
    gt%naux = gs%naux  
    gt%r_size = gs%r_size  
    gt%name = gs%name  
    gt%maxpacksize_t1 = gs%maxpacksize_t1  
    gt%maxpacksize_t2 = gs%maxpacksize_t2
    gt%totpacksize_t2 = gs%totpacksize_t2
    gt%maxnchunk = gs%maxnchunk
    gt%ipackinfo = gs%ipackinfo  
    if (gs%ipackinfo .and. associated(gs%packinfo)) then  
       allocate(gt%packinfo(gt%maxpacksize_t1, gt%nproc, 2))  
       gt%packinfo = gs%packinfo  
    end if
    if (gs%ipackinfo .and. associated(gs%packsize)) then  
       allocate(gt%packsize(gt%nproc, 2, 2))  
       gt%packsize = gs%packsize  
    end if
    if (gs%ipackinfo .and. associated(gs%chunk)) then  
       allocate(gt%chunk(4, 2, gt%maxnchunk, gt%nproc))  
       gt%chunk = gs%chunk  
    end if
    if (gs%ipackinfo .and. associated(gs%pickup)) then  
       allocate(gt%pickup(gt%nproc, 2))  
       gt%pickup = gs%pickup  
    end if
    if (gs%ipackinfo .and. associated(gs%packtype)) then  
       allocate(gt%packtype(gt%nproc))  
       gt%packtype = gs%packtype  
    end if

  end subroutine parallel_gspace_assignment

end module parallel_gspace_module
!
!     ==================== FFT structure =============================
!
module fftstruc_module  

  use constants
  use parallel_gspace_module  

  
  ! edit davegp
  !
  ! include the fftw_hitlist module for creation and destruction of
  !  1D plans
  !
  use fftw_hitlist, only : fftw_hitlist_create, fftw_hitlist_destroy
  

  implicit none

  type fft_struc  
     sequence

     integer :: fftsize(3)        ! the three fftgrid sizes
     integer :: r_size            ! the size of the realspace array
     integer :: naux(4)           ! the size of the four workspaces
     complex(dp), pointer :: &
          unfoldbuf(:)      ! unfolding buffer
     complex(dp), pointer :: &
          rspacebuf(:)      ! for convolutions etc
     complex(dp), pointer :: &
          t1buf(:)          ! buffer for intermediate fft
     complex(dp), pointer :: &
          sendbuf(:)        ! send buffer for mpi
     complex(dp), pointer :: &
          recbuf(:)         ! receive buffer for mpi
     complex(dp), pointer :: &
          workspace(:)      ! space for call to local
     integer, pointer :: &
          mpireq(:)         ! space for mpi requests
     integer, pointer :: &
          mpistat(:)        ! space for mpi status
     integer :: unfoldbufsize     ! sizes of the various buffers
     integer :: t1bufsize  
     integer :: sendbufsize  
     integer :: recbufsize  

  end type fft_struc

contains

  subroutine destroy_fftstruc(ffts)

    implicit none
    type(fft_struc), intent(inout) :: ffts

    if (associated(ffts%unfoldbuf)) deallocate(ffts%unfoldbuf)  
    if (associated(ffts%rspacebuf)) deallocate(ffts%rspacebuf)  
    if (associated(ffts%t1buf)) deallocate(ffts%t1buf)  
    if (associated(ffts%sendbuf)) deallocate(ffts%sendbuf)  
    if (associated(ffts%recbuf)) deallocate(ffts%recbuf)  
    if (associated(ffts%workspace)) deallocate(ffts%workspace)  
    if (associated(ffts%mpireq)) deallocate(ffts%mpireq)  
    if (associated(ffts%mpistat)) deallocate(ffts%mpistat)  

    
    ! edit davegp
    !
    ! if there are  plans created for these sizes destroy them
    call fftw_hitlist_destroy( ffts%fftsize(1) )
    call fftw_hitlist_destroy( ffts%fftsize(2) )
    call fftw_hitlist_destroy( ffts%fftsize(3) )
    

  end subroutine destroy_fftstruc

  ! allocate memory for fft
  subroutine create_fftstruc(ffts, gs, nfn)

    implicit none
    integer, intent(in) :: nfn
    type(fft_struc), intent(inout) :: ffts  
    type(parallel_gspace), intent(in) :: gs  
    integer :: mpistatsize, alcstat

    ffts%fftsize = gs%fftsize(1:3)  
    ffts%r_size = (gs%fftsize(2) * gs%fftsize(3)) / gs%nproc

    if (mod(gs%fftsize(2) * gs%fftsize(3), gs%nproc) > gs%myproc) &
         ffts%r_size = ffts%r_size + 1

    ffts%r_size = ffts%r_size * gs%fftsize(1)           ! realspace
    allocate(ffts%rspacebuf(ffts%r_size * nfn), stat = alcstat)  
    call alccheck('rspacebuf', ffts%r_size * nfn, alcstat)

    ffts%sendbufsize = max(gs%maxpacksize_t1, gs%maxpacksize_t2, &
         gs%totpacksize_t2)
    nullify(ffts%sendbuf)  
    if (ffts%sendbufsize > 0) allocate(ffts%sendbuf(ffts%sendbufsize * nfn), &
         stat = alcstat)
    call alccheck('sendbuf', ffts%sendbufsize * nfn, alcstat)

    ffts%recbufsize = max(gs%maxpacksize_t1, gs%maxpacksize_t2)  
    nullify(ffts%recbuf)  
    if (ffts%recbufsize > 0) allocate(ffts%recbuf(ffts%recbufsize * nfn), &
         stat = alcstat)
    call alccheck('recbuf', ffts%recbufsize * nfn, alcstat)  

    call fft_workspace(gs%nproc, ffts%fftsize, ffts%naux, ffts%unfoldbufsize, &
         ffts%t1bufsize, mpistatsize)                    ! size of fft
    allocate(ffts%mpireq(gs%nproc))  
    allocate(ffts%mpistat(gs%nproc * mpistatsize))  

    allocate(ffts%t1buf(ffts%t1bufsize * nfn), stat = alcstat)  
    call alccheck('t1buf', ffts%t1bufsize * nfn, alcstat)  

    allocate(ffts%workspace(ffts%naux(1) + ffts%naux(2) + ffts%naux(3) + &
         ffts%naux(4)))

    allocate(ffts%unfoldbuf(ffts%unfoldbufsize * nfn), stat = alcstat)
    call alccheck('unfoldbuf', ffts%unfoldbufsize * nfn, alcstat)  

    call fft_iworkspace(ffts%workspace(1), ffts%naux(1), ffts%fftsize(1))

    
    ! edit davegp
    !
    ! create 1D  plans for these grid sizes
    call fftw_hitlist_create( ffts%fftsize(1) )
    call fftw_hitlist_create( ffts%fftsize(2) )
    call fftw_hitlist_create( ffts%fftsize(3) )
    

  end subroutine create_fftstruc

end module fftstruc_module
!
!     ==================== crystal structure =========================
!
module crystal_module  

  use constants
  use iolib_module  
  implicit none

  type crystal  
     sequence  
     
     real(dp) :: avec(3, 3)       ! lattice vectors in cartesian coordinates
     real(dp) :: lvec(3, 3)       ! transpose of avec - cols. are cart.
     real(dp) :: adot(3, 3)       ! adot(i,j) = a_i dot a_j
     real(dp) :: bvec(3, 3)       ! reciprocal lattice vectors
     real(dp) :: bdot(3, 3)       ! b_i dot b_k for reciprocal lattice vectors
     real(dp), pointer :: szet(:)    ! spin polarization szet(ntype)
     real(dp), pointer :: &
          dipole(:)                  ! dipole strength dipole(ntype)
     real(dp), pointer :: coord(:,:,:)  !  positions of the atoms w.r.t. lvec
     !
     !     positions of the atoms rat(3,mxdatm, ntype), with respect to the
     !     lattice vectors, multiplied by a factor of 2 pi
     real(dp), pointer :: rat(:,:,:)     
     real(dp) :: ztot             ! total charge of all atoms
     real(dp), pointer :: zv(:)      ! valence charge zv(ntype)
     real(dp), pointer :: &
          mtsphere(:)                ! size of mt sphere in a
     real(dp), pointer :: mass(:)
     real(dp) :: vcell            ! unit cell volume
     real(dp) :: totalmass        ! the total mass of all the atoms in u
     real(dp) :: net_charge
     integer, pointer :: natom(:)    ! # atoms of each type natom(ntype)
     integer :: nspin             ! number of spins
     integer :: nalpha            ! number of alpha electrons
     integer :: nbeta             ! number of beta electrons
     integer :: ntype             ! number of types
     integer :: mxdatm            ! maximum number of atoms of one type
     integer :: ilsdau
     integer :: igwout
     integer :: ivc
     logical :: icomplex          ! false if crystal has inversion symmetry
     character(len=2), pointer :: &
          nameat(:)                  ! names of atoms nameat(ntype)

  end type crystal

contains  

  subroutine destroy_crystal(crys)

    implicit none
    type(crystal), intent(inout) :: crys

    if (associated(crys%rat)) deallocate(crys%rat)  
    if (associated(crys%natom)) deallocate(crys%natom)  
    if (associated(crys%nameat)) deallocate(crys%nameat)  
    if (associated(crys%szet)) deallocate(crys%szet)  
    if (associated(crys%dipole)) deallocate(crys%dipole)  
    if (associated(crys%zv)) deallocate(crys%zv)  
    if (associated(crys%mtsphere)) deallocate(crys%mtsphere)  

  end subroutine destroy_crystal

  subroutine create_crystal(crys, mxdatm, ntype, rat, natom, nameat, avec, &
       szet, mtsphere, dipole,lvec) 

    implicit none

    type(crystal), intent(out) :: crys  
    integer, intent(in) :: mxdatm, ntype, natom(ntype)
    real(dp), intent(in) :: rat(3, mxdatm, ntype), avec(3, 3), &
         szet(ntype), mtsphere(ntype), dipole(ntype), lvec(3,3)
!    character(len=2), intent(in) :: nameat(ntype)
  character, intent(in) :: &
       nameat(ntype * 2) ! transfer string for name of atoms

  integer i,j

    allocate(crys%rat(3, mxdatm, ntype))  
    allocate(crys%natom(ntype))  
    allocate(crys%nameat(ntype))  
    allocate(crys%szet(ntype))  
    allocate(crys%zv(ntype))  
    allocate(crys%mtsphere(ntype))  
    allocate(crys%dipole(ntype))

    crys%ntype = ntype
    crys%mxdatm = mxdatm
    crys%avec = avec  
    crys%lvec = lvec  
    crys%natom = natom  
    crys%szet = szet  
    crys%mtsphere = mtsphere  
    crys%dipole = dipole  
!    crys%nameat = nameat  
    crys%rat = rat  

  do i = 1, ntype
     do j = 1, 2
        crys%nameat(i)(j:j) = nameat((i - 1) * 4 + j)
     end do
  end do

  end subroutine create_crystal

  ! ------- iolib write crystal

  subroutine write_crystal(fname, crys)  

    implicit none  

    character(len=*), intent(in) :: fname      ! the filename
    type(crystal), intent(in) :: crys          ! the crystal structure to write

    integer :: info, i, j, k, l  
    integer, parameter :: mymarker = 7  

    open(unit = 21, file = fname, position = 'append', status = 'unknown', &
         form = 'formatted')
    write(21, *) mymarker, ' crystal'  
    write(21, *) crys%avec  
    write(21, *) crys%adot  
    write(21, *) crys%bvec  
    write(21, *) crys%bdot  
    write(21, *) crys%ztot, crys%vcell, crys%nspin, crys%ntype, crys%mxdatm, &
         crys%icomplex
    write(21, *) crys%nameat  
    write(21, *) crys%natom  
    write(21, *) crys%szet  
    write(21, *) crys%zv  
    write(21, *) crys%mtsphere  
    write(21, *) crys%rat  
    close(21)

  end subroutine write_crystal
  !
  !
  !
  subroutine read_crystal(fname, crys)  

    implicit none      ! implicit? no!

    character(len=*), intent(in) :: fname      ! the filename
    type(crystal), intent(out) :: crys
  
    integer :: info, i, j, k, l  
    integer, parameter :: mymarker = 7

    open(unit = 21, file = fname, status = 'unknown', form = 'formatted')
    call iolib_find_marker(mymarker, 21, info)  
    if (info == 0) then  
       write(0, *) 'iolib: could not find marker:', mymarker  
       call mystop  
    end if

    read(21, *, end = 102) crys%avec  
    read(21, *, end = 102) crys%adot  
    read(21, *, end = 102) crys%bvec  
    read(21, *, end = 102) crys%bdot  
    read(21, *, end = 102) crys%ztot, crys%vcell, crys%nspin, crys%ntype, &
         crys%mxdatm, crys%icomplex

    allocate(crys%rat(3, crys%mxdatm, crys%ntype))  
    allocate(crys%natom(crys%ntype))  
    allocate(crys%nameat(crys%ntype))  
    allocate(crys%szet(crys%ntype))  
    allocate(crys%zv(crys%ntype))  
    allocate(crys%mtsphere(crys%ntype))  

    read(21, *, end = 102) crys%nameat  
    read(21, *, end = 102) crys%natom  
    read(21, *, end = 102) crys%szet  
    read(21, *, end = 102) crys%zv  
    read(21, *, end = 102) crys%mtsphere  
    read(21, *, end = 102) crys%rat  
    close(21)

    return

102 continue  
    write(0, *) 'iolib: reached eof before reading object:', mymarker
    call mystop  

  end subroutine read_crystal

end module crystal_module
!
!     ==================== symmetry structure =========================
!
module symmetry_module

  use constants
  use iolib_module
  implicit none

  type symmetry  
     sequence  

     real(dp) :: tnp(3, 48)       ! the nonprimitive translations
     real(dp) :: rsymmat(3, 3, 48)! the cartesian symmetry transformations
     integer :: mtrx(3, 3, 48)    ! the transformation matrices
     integer :: rmtrx(3, 3, 48)   ! the transformation matrices in realspace
     integer :: ntrans
     integer :: mtrans            ! symmetry operations used to relate kps in FBZ to those in IBZ
     integer :: ind_mtrans(48)
     integer :: ind_ktrans(48)
     integer :: ind_inv(48)       ! keep track of inverse symmetry operations
     integer :: ind_rot(48,5000)  ! rotation mapping table 

!--------------------
! LSDA+U
! P. Zhang

     integer :: nfrac             ! number of fractional translation operations
     integer :: indfrac(48)       ! index of fractional translation operations

  end type symmetry
  !
  !     Assume the matrix h has the lattice vectors in its COLUMNS.
  !     Suppose a symmetry operation S takes a vector from r to r':
  !
  !         r' = Sr = alpha r + b
  !
  !     where alpha = mirror-rotation,    b = nonprimitive translation.
  !
  !     then    mtrx = (h^(-1) alpha^-1 h)^T
  !
  !     such that mtrx gives the transformation for the G vectors
  !     in reciprocal lattice vector coordinates
  !
  !        G' = alpha G       (in cartesian coordinates)
  !        G' =  mtrx^T * G   (in lattice coordinates)
  !
  !     It turns out that:
  !
  !       rsymmat = alpha
  !
  !
  !     the translations tnp are actually 2pi* b, and given relative
  !     to the realspace lattice coordinates

contains  

  ! ----- iolib write symmetry

  subroutine write_symmetry(fname, syms)

    implicit none

    character(len=*), intent(in) :: fname      ! the filename
    type(symmetry), intent(in) :: syms         ! the symmetry info to write

    integer :: info  
    integer, parameter :: mymarker = 8

    open(unit = 21, file = fname, position = 'append', status = 'unknown', &
         form = 'formatted')
    write(21, *) mymarker, ' symmetry'  
    write(21, *) syms%ntrans  
    write(21, *) syms%mtrx  
    write(21, *) syms%tnp  
    write(21, *) syms%rmtrx  
    write(21, *) syms%rsymmat  
    close(21)

  end subroutine write_symmetry
  !
  !
  !
  subroutine read_symmetry(fname, syms)

    implicit none      ! implicit? no!

    character(len=*), intent(in) :: fname      ! the filename
    type(symmetry), intent(out) :: syms  

    integer :: info  
    integer, parameter :: mymarker = 8

    open(unit = 21, file = fname, status = 'unknown', form = 'formatted')
    call iolib_find_marker(mymarker, 21, info)  
    if (info == 0) then  
       write(0, *) 'iolib: could not find marker:', mymarker  
       call mystop  
    end if

    read(21, *, end = 102) syms%ntrans  
    read(21, *, end = 102) syms%mtrx  
    read(21, *, end = 102) syms%tnp  
    read(21, *, end = 102) syms%rmtrx  
    read(21, *, end = 102) syms%rsymmat  
    close(21)  

    return

102 continue  
    write(0, *) 'iolib: reached eof before reading object:', mymarker
    call mystop  

  end subroutine read_symmetry

end module symmetry_module
!
!     ==================== magnetic field structure ====================
!
module magnetic_field_module

  use constants
  implicit none

  type magnetic_field  
     sequence

     real(dp) :: h             ! strength of magnetic field
     integer :: gvec(4)        ! direction of magnetic field, switch

  end type magnetic_field

end module magnetic_field_module
!
!     ==================== pw parameter structure ======================
!
module pw_parameter_module

  use constants
  use magnetic_field_module
  implicit none

  type pw_parameter  
     sequence  

     type(magnetic_field) :: &
          mfield,&                  ! magnetic field
          epot                    ! electric potential

     real(dp) :: &
          q(3), &                 ! gw use
          alphamix(4), &          ! linear mixing parameters
          smearing, &             ! gaussian smearing parameter eV
          emaxsub, &              ! cutoff energy for submatrix
          emaxmix, &              ! mixing energy cutoff
          emax, &                 ! energy cutoff
          epsdiag, &              ! accuracy of diagonalization
          shiftsafety, &          ! extra shift to avoid bombout in matrix
          random_startvec, &      ! how much random startvector to mix in0, 1
          bandgap, &              ! estimated bandgap in Ryd
          epsmag, &               ! accuracy of cg inversion for mag suscep
          epscv, &                ! convergence precision for potential
          epsce, &                ! convergence precision for energy
          chi_q,&                 ! |q|-vector (relative to reciprocal lattice)
          g0mask(9),&             ! mask for g=0
          ekinmod(3),&            ! modification parameters for kinetic energy
          bd_lgth, &              ! maximal length size for angle calculation
          bra(3),&                ! vector from the left for tensor visual
          ket(3),&                ! vector from the right for tensor visual
          correct_stress_cutoff,&
          MD_alpha,MD_beta,&      ! coefficients for mixing of past WF for MD
          init_dir_econv


     real(dp), pointer :: energywindow(:)        ! energy windows for 
                                                       ! various plots, in eV
     real(dp), pointer :: lineplot(:)   ,&       ! lineplot
                           sliceplot(:)   ,&     ! sliceplot
             pp_pwfocc(:,:)    ! atomic pseudo-wavefunctions occupations

     integer :: &
          checkpoint_wfn, &       ! how frequently wfn is chkpted
          occupy, &               ! how to occupy the bands with electrons
          maxitdiag, &            ! maximum number of iterations for diagonal.
          maxitdiag_band, &       ! max # of iter. for diag. in bandstructure
          maxitcgsol, &           ! maximum number of iterations for cg solver
          maxitscf, &             ! maximum number of iterations for scf loop
          optimize, &             ! flag to trigger optimizations
          miscflag, &             ! miscellaneous flags
          nbandsfft, &            ! maximum number of bands to FFT simultan.
          input, &                ! flag to trigger additional input
          output(4),&                ! flag to trigger additional output
          njobs, &                ! number of jobs to be performed
          nband_tddft,&           ! number of bands for TDDFT
          iscr,&                  ! type of screening used on first iteration
          nplotwv,&               ! number of wave functions to plot
          nlineplot,&             ! number of wave functions to plot
          nsliceplot,&            ! number of sliceplots to do
          ilog,&                  ! flag to trigger log output by pw code
          nenergywindow,&         ! number of energy windows
          nrg,&                   ! # grid points for radial grid integrations
          nmrkpts,&               ! # of k-points calculated in one go (NMR)
          nstay_put,&             ! number of atoms to stay puts
          nprocpvm,&              ! number of procs for PVM
          pp_format,&             ! pseudopotenital format 
          numpk_befptf,&          ! # of pulay_kerker steps before pulay_tf
          tfw_max_cg,&            ! max # of cg steps for solving TFW .eq 
          smearing_method,&       ! choice of Fermi-Dirac, Gaussian, and others
          extrapolation,&         ! extrap method for MD or relaxation
          itdiag_add,&            ! additional iterations for diag
          itdiag_add_metal, &     ! additional iterations for metal diag
          lsdau, &                ! do LSDA+U calculation
          wannier, &              ! output matrix elements for wannier code
          iexact_diag             ! exact diagonalization       

                                  ! (l=0) - 1 (l=1) - 2  (l=2) - 3 
     integer, pointer :: plotwv(:)    ! memory location plotwv is located
     integer, pointer :: joblist(:)    ! the job codes
     integer, pointer :: stay_put(:)   ! list of atoms which stay put
     integer, pointer :: pp_iloc(:)    ! PP to be used for local
     real(dp), pointer :: NLPP_rcut(:)    ! cutoff radius for real space

     logical :: isymmetrize       ! true if k-point symmetrization is necessary
     logical NL_rspace(2)
     logical io_scf
     character, pointer :: bslabels(:)    ! a list of character labels
                                                ! for the bstruc plot
     character(len=2) :: icorr    ! the correlation type
     character(len=20) :: diag_method    ! method use for diagonalization
     character(len=20) :: mix_method      ! method use for diagonalization

! weiweigao
     integer :: gw_sum = 0        ! for energy integration of GW calculation, default is 0 
! weiweigao
       
! davegp
     logical :: nscf              ! flag indicating non-self-consistent field
                                  ! no update of potential each scf step
! davegp

  end type pw_parameter

end module pw_parameter_module
!
!     ==================== kpoint structure =========================
!
module kpoint_module  

  use constants
  use iolib_module                ! because here already the iolib functions
                                  ! are available as regular member functions
  type kpoint
     sequence
  
     integer :: nrk               ! number of k-points
     integer :: grid(3)           ! the generation grid
     integer :: nrk_fbz       ! number of k-points

     real(dp) :: shift(3)         ! shift vector
     logical :: generated         ! true if generated, false if read
     logical :: reduced           ! whether it is a reduced k-point grid
     integer, pointer :: &
          kmap(:)   , &     ! map: the irred kpoint of each reducib.
          ind(:,:), &       ! index of six nearest neighbors: ind(6, nrk)
          ind_symtr(:)      ! index of symmetry operation which connects a k-point in
                            ! full BZ to a k-point in IBZ; if the i-th k-point in the FBZ
                            ! connects to -k (k is in IBZ) through a symmetery op j,
                            ! then ind_symtr(i)=-j

     real(dp), pointer :: &
          w(:)   , &        ! symmetry weights w(nrk)
          label(:), &       ! label for bandstructure plot
          ind_G0(:,:)       ! Umklapp vectors

     real(dp), pointer :: &
          rk(:,:),rk_fbz(:,:)           ! k-points rk(3, nrk)


  end type kpoint

contains  

  subroutine destroy_kpoint(kp)

    implicit none
    type(kpoint), intent(inout) :: kp  

    if (associated(kp%w)) deallocate(kp%w)  
    if (associated(kp%rk)) deallocate(kp%rk)  
    if (associated(kp%ind)) deallocate(kp%ind)  
    if (associated(kp%label)) deallocate(kp%label)  

  end subroutine destroy_kpoint

  subroutine create_kpoint(kp, nrk)

    implicit none
    type(kpoint), intent(out) :: kp  
    integer, intent(in) :: nrk

    kp%nrk = nrk  
    allocate(kp%w(nrk))
    allocate(kp%label(nrk))  
    allocate(kp%rk(3, nrk))
    allocate(kp%ind(6, nrk))

    kp%label = dzero
    kp%ind = 0 
    kp%w = dzero
    kp%rk = dzero

  end subroutine create_kpoint

  subroutine write_kpoint(fname, struc)  

    implicit none      ! implicit? no!

    character(len=*), intent(in) :: fname      ! the filename
    type(kpoint), intent(in) :: struc
  
    !     ----------- local variables
    integer :: info, i, j, k, l  
    integer, parameter :: mymarker = 4

    open(unit = 21, file = fname, position = 'append', status = 'unknown', &
         form = 'formatted')
    write(21, *) mymarker, ' kpoint'  
    write(21, *) struc%nrk  
    write(21, *) struc%grid  
    write(21, *) struc%shift  
    write(21, *) struc%generated  
    write(21, *) struc%reduced  
    write(21, *) struc%kmap  
    write(21, *) struc%ind  
    write(21, *) struc%w  
    write(21, *) struc%label  
    write(21, *) struc%rk  
    close(21)

    return  

  end subroutine write_kpoint

  subroutine read_kpoint(fname, struc)  

    implicit none      ! implicit? no!

    character(len=*), intent(in) :: fname      ! the filename
    type(kpoint), intent(out) ::  struc

    !     ----------- local variables
    integer :: info, i, j, k, l  
    integer, parameter :: mymarker = 4

    open(unit = 21, file = fname, status = 'unknown', form = 'formatted')
    call iolib_find_marker(mymarker, 21, info)  
    if (info == 0) then  
       write(0, *) 'iolib: could not find marker:', mymarker  
       call mystop  
    end if

    read(21, *, end = 102) struc%nrk  
    allocate(struc%w(struc%nrk))  
    allocate(struc%rk(3, struc%nrk))  
    allocate(struc%ind(6, struc%nrk))  
    allocate(struc%label(struc%nrk))  
    read(21, *, end = 102) struc%grid  
    read(21, *, end = 102) struc%shift  
    read(21, *, end = 102) struc%generated  
    read(21, *, end = 102) struc%reduced  
    allocate(struc%kmap(struc%grid(1) * struc%grid(2) * struc%grid(3)))
    read(21, *, end = 102) struc%kmap  
    read(21, *, end = 102) struc%ind  
    read(21, *, end = 102) struc%w  
    read(21, *, end = 102) struc%label  
    read(21, *, end = 102) struc%rk  
    close(21)  

    return

102 continue  
    write(0, *) 'iolib: reached eof before reading object:', mymarker

    call mystop  

  end subroutine read_kpoint

end module kpoint_module
!
!     ==================== pseudopotential  structure ==================
!
module pseudo_potential_module  

  use constants
  implicit none

  type pseudo_potential  
     sequence  

     integer, pointer :: nqnl(:)    ! mesh points of nonlocal pot (ntype)
     integer, pointer :: nkb(:,:)   ! norm for nonlocal pot nkb(3, ntype)
     integer, pointer :: lo(:,:)    ! l quantum number for NL projector
     real(dp), pointer :: &
          delqnl(:)                 ! nloc pot mesh dist delqnl(ntype)
     real(dp), pointer :: &
          vkb(:,:,:)                ! nloc pspot vkb(mxdlqp, 3, ntype)
     real(dp), pointer :: &
          d2vkbdq2(:,:,:)           ! d^2 vkb/ dq^2
     integer :: nanl              ! number of nonlocal projectors
     integer :: ntype             ! number of atomic types
     integer :: mxdlqp            ! max number of mesh points for nonloc pot
     integer :: natom_tot         ! # of total atoms
     integer, pointer::&
          nrr(:)            ! # of points for real space mesh
     real(dp), pointer :: &       ! radial points
          r(:,:)         
     real(dp), pointer :: &       ! real space potential for each projector 
          vr_pp(:,:,:)                
     real(dp), pointer :: &       ! real space pseudo-wfs for each projector
          wr_pp(:,:,:)                     
     ! Following are used for NON-local real space implementation
     real(dp), pointer :: &       ! real space potential for each projector 
          vkb_mask(:,:,:)      
     real(dp), pointer ::  kb_rsp(:)  ,& ! real space projectors
               xyzmap(:)  , &    ! xyz dist. from atom of mesh point
               qi(:,:)          ! 1-d q-space interpolation positions
     real(dp) amr(201),ri(201)
     integer :: NLPP_rsp_pts      ! # of total used mesh pts for all atoms
     integer :: NLPP_rsp_nrefpts  ! # of total mesh pts for all projectors
     integer :: mrb2_matom_node   ! # of pts per PE for a given cutoff radius
                  
     integer,pointer ::  nmap(:)  ,& !# mesh pts. within rcut for an atom
            indm(:)  ,&      ! memory pos. of a mesh pt for WF
            numref(:)        ! # of proj. for atom
     logical NL_rspace(2)

  end type pseudo_potential

contains  

  subroutine destroy_pseudo_potential(pspots)

    implicit none
    type(pseudo_potential), intent(inout) :: pspots  

    if (associated(pspots%nqnl)) deallocate(pspots%nqnl)  
    if (associated(pspots%nkb)) deallocate(pspots%nkb)  
    if (associated(pspots%lo)) deallocate(pspots%lo)  
    if (associated(pspots%delqnl)) deallocate(pspots%delqnl)  
    if (associated(pspots%vkb)) deallocate(pspots%vkb)  
    if (associated(pspots%d2vkbdq2)) deallocate(pspots%d2vkbdq2)  
    if (associated(pspots%nrr)) deallocate(pspots%nrr)
    if (associated(pspots%r)) deallocate(pspots%r)
    if (associated(pspots%vr_pp)) deallocate(pspots%vr_pp)
    if (associated(pspots%wr_pp)) deallocate(pspots%wr_pp)
    if (associated(pspots%qi)) deallocate(pspots%qi)
    if (associated(pspots%vkb_mask)) deallocate(pspots%vkb_mask)

    if(pspots%NL_rspace(1)) then
      if (associated(pspots%kb_rsp)) deallocate(pspots%kb_rsp)
      if (associated(pspots%xyzmap)) deallocate(pspots%xyzmap)
      if (associated(pspots%nmap)) deallocate(pspots%nmap)
      if (associated(pspots%indm)) deallocate(pspots%indm)
      if (associated(pspots%numref)) deallocate(pspots%numref)
    end if

  end subroutine destroy_pseudo_potential
  
end module pseudo_potential_module
!
!     ==================== gspace_array structure ==================
!
module gspace_array_module  

  use constants
  use parallel_gspace_module
  implicit none

  type complex_gspace_array  
     sequence  
 
     complex(dp), pointer :: &
          data(:,:,:)       ! data(length, nvecs, nspin)
     integer :: nspin             ! number of spin components
     integer :: nvecs             ! number of one-dimensional arrays
     integer :: length            ! length of a single vector

  end type complex_gspace_array

  type double_gspace_array  
     sequence  


     real(dp), pointer :: &
          data(:,:,:)       ! data(length, nvecs ,nspin)
     integer :: nspin             ! number of spin components
     integer :: nvecs             ! number of one-dimensional arrays
     integer :: length            ! length of a single vector

  end type double_gspace_array

contains  

  subroutine destroy_complex_gspace_array(cgsarrays)

    implicit none
    type(complex_gspace_array), intent(inout) :: cgsarrays

    if (associated(cgsarrays%data)) deallocate(cgsarrays%data)  

  end subroutine destroy_complex_gspace_array

  subroutine create_complex_gspace_array(cgsarrays, len, nvecs, nspin)

    implicit none
    type(complex_gspace_array), intent(out) :: cgsarrays  
    integer, intent(in) :: len, nvecs, nspin
    integer :: alcstat
  
    if (len > 0 .and. nvecs > 0 .and. nspin > 0) then  
       allocate(cgsarrays%data(len, nvecs, nspin), stat = alcstat)  
       call alccheck('gsarray', len * nvecs * nspin, alcstat)  
    else  
       nullify(cgsarrays%data)  
    end if
    cgsarrays%nspin = nspin  
    cgsarrays%length = len  
    cgsarrays%nvecs = nvecs  

  end subroutine create_complex_gspace_array

  subroutine print_complex_gspace_array(gs, cgsa)

    implicit none
    type(complex_gspace_array), intent(in) :: cgsa  
    type(parallel_gspace), intent(in) :: gs  
    integer :: i, is, n

    do is = 1, cgsa%nspin  
       write(9, *) 'spin: ', is  
       do n = 1, cgsa%nvecs  
          write(9, * ) 'vector: ', n  
          do i = 1, gs%length  
             write(9, *) i, cgsa%data(i, n, is)  
          end do
       end do
    end do

  end subroutine print_complex_gspace_array

  subroutine destroy_double_gspace_array(dgsarrays)

    implicit none
    type (double_gspace_array), intent(inout) :: dgsarrays

    if (associated(dgsarrays%data)) deallocate(dgsarrays%data)  

  end subroutine destroy_double_gspace_array

  subroutine create_double_gspace_array(dgsarrays, len, nvecs, nspin)

    implicit none
    type(double_gspace_array), intent(out) :: dgsarrays  
    integer, intent(in) :: len, nvecs, nspin
    integer :: alcstat

    if (len > 0 .and. nvecs > 0 .and. nspin > 0) then  
       allocate(dgsarrays%data(len, nvecs, nspin), stat = alcstat)  
       call alccheck('gsarray', len * nvecs * nspin, alcstat)  
    else  
       nullify(dgsarrays%data)  
    end if
    dgsarrays%nspin = nspin  
    dgsarrays%length = len  
    dgsarrays%nvecs = nvecs  

  end subroutine create_double_gspace_array

end module gspace_array_module
!
!     ==================== rspace_array structure ==================
!
module rspace_array_module  

  use constants
  implicit none

  type complex_rspace_array  
     sequence  

     complex(dp), pointer :: &
          data(:,:,:)                
     integer :: nspin             ! number of spin components
     integer :: nvecs             ! number of one-dimensional arrays
     integer :: length            ! length of a single vector

  end type complex_rspace_array

contains  

  subroutine destroy_complex_rspace_array(crsarrays)

    implicit none
    type(complex_rspace_array), intent(inout) :: crsarrays

    if (associated(crsarrays%data)) deallocate(crsarrays%data)  

  end subroutine destroy_complex_rspace_array

  subroutine create_complex_rspace_array(crsarrays, len, nvecs, nspin)

    implicit none
    type(complex_rspace_array), intent(out) :: crsarrays  
    integer, intent(in) :: len, nvecs, nspin
    integer :: alcstat

    allocate(crsarrays%data(len, nvecs, nspin), stat = alcstat)  
    call alccheck('rs_array', nvecs * nspin * len, alcstat)  
    crsarrays%nspin = nspin  
    crsarrays%length = len  
    crsarrays%nvecs = nvecs  

  end subroutine create_complex_rspace_array

end module rspace_array_module
!
!     ==================== blochl operator structure ==================
!
module blochl_operator_module
  use constants
  
  type blochl_operator
     sequence
     integer, pointer :: nqnl(:)   !mesh pts of nonlocal proj nqnl(ntype)
     integer, pointer :: nkb(:,:)  ! norm of NL proj nkb(nlprjmx,ntype)
     integer, pointer :: lo(:,:)   ! l quantum number for NL projector
     integer, pointer :: numl(:,:)      ! Number of prj in each channel
     integer, pointer :: ln2np(:,:,:)   ! l,n -> num. of projector
     
     real(dp), pointer :: delqnl(:)    !nloc proj mesh dist delqnl(ntype)
     real(dp), pointer :: r_rec(:,:)   ! Core radii for reconstruction 
     real(dp), pointer :: vkb(:,:,:)   ! nl prj vkb(mxdlqp,nlprjmx,ntype)
     real(dp), pointer :: d2vkbdq2(:,:,:)    ! d^2 vkb/ dq^2
     
     integer nanl              ! number of nonlocal projectors
     integer ntype             ! number of atomic types
     integer mxdlqp            ! max number of mesh points for nonloc pot
     integer nlprjmx           ! max number of l-projs     
 
     real(dp), pointer :: delta_dia(:,:,:,:)   ! The recon dia  mat. ele
     real(dp), pointer :: delta_para(:,:,:,:)  ! The recon para mat. ele
    
  end type blochl_operator
  
contains
  subroutine destroy_blochl_operator(blop)

    implicit none
    type(blochl_operator), intent(inout) :: blop

    if(associated(blop%nqnl)) deallocate(blop%nqnl)
    if(associated(blop%nkb)) deallocate(blop%nkb)
    if(associated(blop%lo)) deallocate(blop%lo)
    if(associated(blop%r_rec)) deallocate(blop%r_rec)
    if(associated(blop%numl)) deallocate(blop%numl)
    if(associated(blop%ln2np)) deallocate(blop%ln2np)
    if(associated(blop%delqnl)) deallocate(blop%delqnl)
    if(associated(blop%vkb)) deallocate(blop%vkb)
    if(associated(blop%d2vkbdq2)) deallocate(blop%d2vkbdq2)

  end subroutine destroy_blochl_operator
  
end module blochl_operator_module


!
!     ==================== energy structure =========================
!
module energy_module  

  use constants
  implicit none

  type energy  
     sequence

     real(dp) :: &                ! all energies in Rydbergs
          inputfermi, &           ! fermi energy as given by input
          afmagmom, &             ! antiferromagnetic sublattice moment in Bohr
          vxc0(2), &              ! spin up and down of the xc potential at G=0
          magnetic, &             ! magnetic energy due to external field
          hxc, &                  ! hartree exchange-correlation correction
          xc, &                   ! the exchange-correlation energy
          ewald, &                ! ewald energy, with compensating background
          eband, &                ! band energy
          ion, &                  ! ion-electron energy
          nloc, &                 ! nonlocal energy contribution
          hartree, &              ! hartree electron-electron energy
          ektot, &                ! total kinetic energy
          esmear, &               ! smearing energy contribution
          efermi(2), &            ! fermi energy for both spins
          evalshift(2), &         ! shift of the energy eigenvalues
          total, &                ! the total energy
          totnonv, &              ! etot - 1/2 TS term, is not variational
                                  !  it is the estimate of the zero T energy
          alpha                   ! the alpha term which contains divergences
       

  end type energy

contains  

  subroutine energy_assignment(et, es)

    implicit none
    type(energy), intent(out) :: et  
    type(energy), intent(in) :: es
  
    et%inputfermi = es%inputfermi  
    et%afmagmom = es%afmagmom  
    et%vxc0 = es%vxc0  
    et%magnetic = es%magnetic  
    et%hxc = es%hxc  
    et%xc = es%xc  
    et%ewald = es%ewald  
    et%eband = es%eband  
    et%ion = es%ion  
    et%nloc = es%nloc  
    et%hartree = es%hartree  
    et%ektot = es%ektot  
    et%esmear = es%esmear  
    et%efermi = es%efermi  
    et%total = es%total  
    et%alpha = es%alpha  
    et%totnonv = es%totnonv  
    et%evalshift = es%evalshift  

  end subroutine energy_assignment

  subroutine energy_assignment0(et, e0)

    implicit none
    type(energy), intent(out) :: et  
    real(dp), intent(in) :: e0  

    et%inputfermi = e0  
    et%afmagmom = e0  
    et%vxc0 = e0  
    et%magnetic = e0  
    et%hxc = e0  
    et%xc = e0  
    et%ewald = e0  
    et%eband = e0  
    et%ion = e0  
    et%nloc = e0  
    et%hartree = e0  
    et%ektot = e0  
    et%esmear = e0  
    et%efermi = e0  
    et%total = e0  
    et%totnonv = e0  
    et%alpha = e0  
    et%evalshift = e0  

  end subroutine energy_assignment0

  function energy_subtract(et, es) result(en_s)

    type(energy) en_s  
    type(energy), intent(in) :: et  
    type(energy), intent(in) :: es  

    en_s%inputfermi = et%inputfermi  
    en_s%afmagmom = et%afmagmom - es%afmagmom  
    en_s%vxc0 = et%vxc0 - es%vxc0  
    en_s%magnetic = et%magnetic - es%magnetic  
    en_s%hxc = et%hxc - es%hxc  
    en_s%xc = et%xc - es%xc  
    en_s%ewald = et%ewald-es%ewald  
    en_s%eband = et%eband-es%eband  
    en_s%ion = et%ion - es%ion  
    en_s%nloc = et%nloc - es%nloc  
    en_s%hartree = et%hartree-es%hartree  
    en_s%ektot = et%ektot - es%ektot  
    en_s%esmear = et%esmear - es%esmear  
    en_s%efermi = et%efermi - es%efermi  
    en_s%total = et%total - es%total  
    en_s%totnonv = et%totnonv - es%totnonv  
    en_s%alpha = et%alpha - es%alpha  
    en_s%evalshift = et%evalshift - es%evalshift  

  end function energy_subtract

end module energy_module
!
!     ==================== force structure =========================
!
module force_module  

  use constants
  implicit none

  type force  
     sequence  

     integer :: mxdatm  
     integer :: ntype  
     real(dp) :: stress(6)        ! the symmetric stress xx,yy,zz,xy,yz,xz
     real(dp), pointer :: &
          force(:,:,:)      ! forces on the atoms force(3, mxdatm, ntype)

  end type force

contains  

  subroutine destroy_force(f)

    implicit none
    type(force), intent(inout) :: f

    if (associated(f%force)) deallocate(f%force)  

  end subroutine destroy_force

end module force_module
!
!     ==================== band structure =========================
!
module band_module  

  use constants
  use iolib_module
  implicit none

  type band  
     sequence  

     integer :: nrk               ! number of kpoints
     integer :: min(2)               ! minimum number of bands computed
     integer :: max               ! maximum number of bands computed
     integer :: nband1            ! for used in flevel, P. Zhang
     integer :: nband0            ! for used in flevel, P. Zhang
     integer :: nspin             ! number of spins
     integer :: num_gwout         ! number of bands above
!
! the detail of band index reordering is still to be implemented.
! please modify gwout.f90
     integer :: gw_band_reordering      ! reordering band index (in case the LDA screws up band ordering)

     real(dp) :: gw_mid_energy
     real(dp) :: gw_low_energy
     real(dp) :: gw_high_energy

     integer, pointer :: &
          nband(:,:)   , &  ! number of bands computed nband(nrk, nspin)
          ifmax(:,:),&        ! highest occupied band ifmax(nrk, nspin)
          gwout_index(:,:,:) ! index of outputting certain bands for gw code

     real(dp), pointer :: &
          energy(:,:,:)  ,& ! eigenvalues energy(max, nrk, nspin)
          ekn(:,:,:)  ,&    ! kinetic energy(max, nrk, nspin)
          occup(:,:,:)      ! occupation numbers occup(max, nrk, nspin)

  end type band

contains  

  subroutine destroy_band(b)

    implicit none
    type(band), intent(inout) :: b

    if (associated(b%ifmax)) deallocate(b%ifmax)  
    if (associated(b%energy)) deallocate(b%energy)  
    if (associated(b%ekn)) deallocate(b%ekn)  
    if (associated(b%occup)) deallocate(b%occup)  
    if (associated(b%nband)) deallocate(b%nband)  
  end subroutine destroy_band

  subroutine create_band(b, nrk, max, nspin)

    implicit none
    type(band), intent(out) :: b  
    integer, intent(in) :: nrk, max, nspin

    b%nrk = nrk  
    b%max = max  
    b%nspin = nspin     
    allocate(b%nband(nrk, nspin)) ; b%nband = 0  
    allocate(b%occup(max, nrk, nspin)) ;b%occup = dzero
    allocate(b%energy(max, nrk, nspin)) ;b%energy = dzero
    allocate(b%ekn(max, nrk, nspin)) ;b%ekn = dzero
    allocate(b%ifmax(nrk, nspin)) ;b%ifmax = 0  
    allocate(b%gwout_index(max, nrk, nspin)) ;b%gwout_index = 0
  end subroutine create_band

  subroutine write_band(fname, struc)  

    implicit none      ! implicit? no!

    character(len=*), intent(in) :: fname      ! the filename
    type(band), intent(in) :: struc

    !     ----------- local variables
    integer :: info, i, j, k, l  
    integer, parameter :: mymarker = 6

    open(unit = 21, file = fname, position = 'append', status = 'unknown', &
         form = 'formatted')
    write(21, *) mymarker, ' band info'  
    write(21, *) struc%nrk, struc%min, struc%max, struc%nspin  
    write(21, *) ((struc%nband(k, j), k = 1, struc%nrk), j = 1, struc%nspin)
    write(21, *) ((struc%ifmax(k, j), k = 1, struc%nrk), j = 1, struc%nspin)

    write(21, *) (((struc%energy(l, k, j), l = 1, struc%max), &
         k = 1, struc%nrk), j = 1, struc%nspin)
    write(21, *) (((struc%ekn(l, k, j), l = 1, struc%max), &
         k = 1, struc%nrk), j = 1, struc%nspin)
    write(21, *) (((struc%occup(l, k, j), l = 1, struc%max), &
         k = 1, struc%nrk), j = 1, struc%nspin)
    close(21)  
    return  

  end subroutine write_band

  subroutine read_band(fname, bands)  

    implicit none      ! implicit? no!

    character(len=*), intent(in) :: fname      ! the filename
    type(band), intent(inout) :: bands
  
    !     ----------- local variables

    integer :: info, i, j, k, l, nrk, mmin(2), mmax, nspin  
    integer, parameter :: mymarker = 6 

    open(unit = 21, file = fname, status = 'unknown', form = 'formatted')
    call iolib_find_marker(mymarker, 21, info)  
    if (info == 0) then  
       bands%nrk = -1  
       close(21)  
       return  
    end if
    read(21, *, end = 102) nrk, mmin, mmax, nspin  
    if (nrk /= bands%nrk) then  
       write(9, *) 'read_band: kpoint mismatch: kpoints in file:', nrk
       call mystop  
    end if
    if (nspin /= bands%nspin) then  
       write(9, *) 'read_band: spin mismatch, spins in file:', nspin
       call mystop  
    end if
    if (mmax > bands%max) then  
       write(9, *) 'read_band: too many bands in file:', mmax  
       call mystop  
    end if
    if ((.not. associated(bands%nband)) .or. &
         (.not. associated(bands%ifmax)) .or. &
         (.not. associated(bands%energy)) .or. &
         (.not. associated(bands%ekn)) .or. &
         (.not. associated(bands%occup))) then
       write(9, *) 'read_band: structure not allocated!'
       call mystop  
    end if
    read(21, *, end = 102) ((bands%nband(k, j), k = 1, nrk), j = 1, nspin)
    read(21, *, end = 102) ((bands%ifmax(k, j), k = 1, nrk), j = 1, nspin)

    read(21, *, end = 102) (((bands%energy(l, k, j), l = 1, mmax), &
         k = 1, nrk), j = 1, nspin)
    read(21, *, end = 102) (((bands%ekn(l, k, j), l = 1, mmax), &
         k = 1, nrk), j = 1, nspin)
    read(21, *, end = 102) (((bands%occup(l, k, j), l = 1, mmax), &
         k = 1, nrk), j = 1, nspin)
    close(21)  

    return

102 continue  
    write(0, *) 'iolib: reached eof when reading bands'  

    call mystop  

  end subroutine read_band


end module band_module
!
!     ==================== hamiltonian structure =======================
!
!     by convention, only the vnonloc part of the Hamiltonian gets
!     allocated and deallocated upon creation and destruction.
!     For the dense Hamiltonian, different rules apply
module hamiltonian_module  

  use constants
  use parallel_gspace_module  
  use rspace_array_module  
  use gspace_array_module  
  use pseudo_potential_module  
  use fftstruc_module
  implicit none

  type hamiltonian  
     sequence


     real(dp) :: shift            ! global energy shift for the hamiltonian
     real(dp) :: ekinmod(3)       ! modification parameters for kinetic energy
     type (parallel_gspace), pointer :: gspace     
      complex(dp), pointer :: &
          vloc(:)           ! local potential
     real(dp), pointer :: &
          xnorm(:)          ! normalization xnorm(nanl)
     type(fft_struc), pointer :: fftstruc     
     type(pseudo_potential), pointer :: &
          pspot             ! non-local potential pseudopotential
     type(complex_gspace_array) :: &
          vnloc                   ! nonlocal potential projectors
     integer :: dim               ! local number plane waves for this hamilt.
     real(dp), pointer :: &
          rsp_norm(:)       ! normalization xnorm(nanl)
     complex,pointer :: cphase(:)   ! phase for mesh pt in real space 
                                          ! for the current k-point
  end type hamiltonian

  type dense_hamiltonian  
     sequence  

     type(parallel_gspace), pointer :: gspace     
     complex(dp), pointer :: &
          matrix(:)                  ! the dense matrix itself
     real(dp), pointer :: &
          energy(:)                  ! the spectrum
     integer :: neig                       ! number of eigenvalues

  end type dense_hamiltonian

contains  

  subroutine destroy_hamiltonian(h)

    implicit none
    type(hamiltonian), intent(inout) :: h
 
    call destroy_complex_gspace_array(h%vnloc)  
    if (associated(h%xnorm)) deallocate(h%xnorm)  

    if (h%pspot%NL_rspace(1)) then 
      if (associated(h%rsp_norm)) deallocate(h%rsp_norm)
      if (associated(h%cphase)) deallocate(h%cphase)
    end if

  end subroutine destroy_hamiltonian

  subroutine create_hamiltonian(h, gs, pspot, vloc, ffts)

    type(hamiltonian), intent(out) :: h  
    type(parallel_gspace), intent(in), target :: gs  
    type(pseudo_potential), intent(in), target :: pspot  
    type(fft_struc), intent(in), target :: ffts  
    complex(dp), pointer :: vloc(:)              ! local potential

    nullify(h%gspace)  
    h%gspace => gs  
    h%dim = gs%length  
    nullify(h%vloc)  
    h%vloc => vloc  
    nullify(h%fftstruc)  
    h%fftstruc => ffts  
    h%pspot=>pspot
    if (pspot%NL_rspace(1))then
       allocate(h%cphase(pspot%NLPP_rsp_pts))
       allocate(h%rsp_norm(pspot%nanl))
    end if
    if (pspot%nanl > 0) then  
       call create_complex_gspace_array(h%vnloc, gs%length, pspot%nanl, 1)
       allocate(h%xnorm(pspot%nanl))  
    else  
       h%vnloc%nspin = 0  
       h%vnloc%nvecs = 0  
       h%vnloc%length = 0  
       nullify(h%vnloc%data)  
       nullify(h%xnorm)  
    end if

  end subroutine create_hamiltonian

  subroutine destroy_dense_hamiltonian(h)

    implicit none
    type(dense_hamiltonian), intent(inout) :: h  

    nullify(h%gspace)  
    if (associated(h%matrix)) deallocate(h%matrix)  
    if (associated(h%energy)) deallocate(h%energy)  
    h%neig = 0  

  end subroutine destroy_dense_hamiltonian

  subroutine create_dense_hamiltonian(h, gs)

    implicit none
    type(dense_hamiltonian), intent(out) :: h  
    type(parallel_gspace), intent(in), target :: gs

    h%gspace => gs  
    nullify(h%matrix)  
    allocate(h%energy(gs%length))  
    h%neig = 0  

  end subroutine create_dense_hamiltonian

end module hamiltonian_module
!
!     =========== vqmc_gspace_module ==========================
!
module vqmc_gspace_module  

  use constants
  use iolib_module  
  implicit none

  type gspace  
     sequence                     ! want sequential memory layout.

     complex(dp), pointer :: &
          phase(:)                   ! phases(length)
     integer, pointer :: &
          gvec(:,:)                  ! the gspace vectors themselves
     integer, pointer :: inds(:)     ! maps gvecs to stars inds(length)
     integer, pointer :: mstar(:)    ! size of stars mstar(nstar)
     integer :: fftsize(3)        ! FFT grid boundaries for gspace
     integer :: length            ! the total number of g-space vectors
     integer :: nstar             ! total number of stars in gspace
     logical :: istar             ! if star arrays are available

  end type gspace

contains  

  subroutine destroy_gspace(gs)

    implicit none
    type(gspace), intent(inout) :: gs

    if (associated(gs%gvec)) deallocate(gs%gvec)  
    if (associated(gs%inds)) deallocate(gs%inds)  
    if (associated(gs%mstar)) deallocate(gs%mstar)  
    if (associated(gs%phase)) deallocate(gs%phase)  

  end subroutine destroy_gspace

  subroutine write_gspace(fname, gsstruc)  

    implicit none      ! implicit? no!
    character(len=*), intent(in) :: fname   ! the filename
    type(gspace), intent(in) :: gsstruc     ! the gspace structure

    integer :: info, i, j  
    integer, parameter :: mymarker = 1      ! the marker for the gspace object

    open(unit = 21, file = fname, position = 'append', status = 'unknown', &
         form = 'formatted')
    write(21, *) mymarker, ' gspace'  
    write(21, *) (gsstruc%fftsize(i), i = 1, 3)  
    write(21, *) gsstruc%length  
    write(21, *) ((gsstruc%gvec(j, i), j = 1, 3), i = 1, gsstruc%length)
    write(21, *) gsstruc%istar  
    if (gsstruc%istar) then  
       write(21, *) (gsstruc%inds(i), i = 1, gsstruc%length)  
       write(21, *) (gsstruc%phase(i), i = 1, gsstruc%length)  
       write(21, *) gsstruc%nstar  
       write(21, *) (gsstruc%mstar(i), i = 1, gsstruc%nstar)  
    else  
       write(21, *)
       write(21, *)
       write(21, *)
       write(21, *)
    end if
    close(21)  
    return  

  end subroutine write_gspace

  subroutine read_gspace(fname, gsstruc)  

    implicit none      ! implicit? no!
    character(len=*), intent(in) :: fname    ! the filename
    type(gspace), intent(out) :: gsstruc     ! the gspace structure

    integer :: info, i, j  
    integer, parameter :: mymarker = 1       ! the marker for the gspace object

    open(unit = 21, file = fname, status = 'unknown', form = 'formatted')
    call iolib_find_marker(mymarker, 21, info)  
    if (info == 0) then  
       write(0, *) 'iolib: could not find marker:', mymarker  
       call mystop  
    end if
    read(21, *, end = 102) (gsstruc%fftsize(i), i = 1, 3)  
    read(21, *, end = 102) gsstruc%length  
    allocate(gsstruc%gvec(3, gsstruc%length))  
    read(21, *, end = 102) ((gsstruc%gvec(j, i), j = 1, 3), &
         i = 1, gsstruc%length)
    read(21, *) gsstruc%istar  
    if (gsstruc%istar) then  
       allocate(gsstruc%inds(gsstruc%length))  
       read(21, *) (gsstruc%inds(i), i = 1, gsstruc%length)  
       allocate(gsstruc%phase(gsstruc%length))  
       read(21, *) (gsstruc%phase(i), i = 1, gsstruc%length)  
       read(21, *) gsstruc%nstar  
       allocate(gsstruc%mstar(gsstruc%nstar))  
       read(21, *) (gsstruc%mstar(i), i = 1, gsstruc%nstar)  
    else  
       gsstruc%nstar = 0  
       nullify(gsstruc%phase)
       nullify(gsstruc%mstar)
       nullify(gsstruc%inds)
       read(21, *)
       read(21, *)
       read(21, *)
       read(21, *)
    end if
    close(21)  
    return
  
102 continue  
    write(0, *) 'iolib: reached eof before reading object:', mymarker
    call mystop  

  end subroutine read_gspace

end module vqmc_gspace_module
!
!     ==================== isrt structure ===============
!
module vqmc_isrt_module  

  use constants
  use iolib_module
  implicit none

  type isrt  
     sequence  

     integer, pointer :: kndim(:)    ! # of plane-waves at each k-point
     integer, pointer :: &
          isrt(:,:)         ! the actual isrt array: isrt(lda, nrk)
     integer :: nrk               ! number of k-points
     integer :: lda               ! leading dimension of isrt array

  end type isrt

contains  

  subroutine destroy_isrt(is)

    implicit none
    type(isrt), intent(inout) :: is

    if (associated(is%kndim)) deallocate(is%kndim)  
    if (associated(is%isrt)) deallocate(is%isrt)  

  end subroutine destroy_isrt

  subroutine write_isrt(fname, struc)  

    implicit none      ! implicit? no!
    character(len=*), intent(in) :: fname      ! the filename
    type(isrt), intent(in) :: struc  

    integer :: info, i, j  
    integer, parameter :: mymarker = 2

    open(unit = 21, file = fname, position = 'append', status = 'unknown', &
         form = 'formatted')
    write(21, *) mymarker, ' isrt'  
    write(21, *) struc%nrk  
    write(21, *) (struc%kndim(i), i = 1, struc%nrk)  
    write(21, *) struc%lda  
    write(21, *) ((struc%isrt(i, j), i = 1, struc%lda), j = 1, struc%nrk)
    close(21)  
    return  

  end subroutine write_isrt

  subroutine read_isrt(fname, struc)  

    implicit none      ! implicit? no!

    character(len=*), intent(in) :: fname      ! the filename
    type(isrt), intent(out) :: struc  

    integer :: info, i, j  
    integer, parameter :: mymarker = 2

    open(unit = 21, file = fname, status = 'unknown', form = 'formatted')
    call iolib_find_marker(mymarker, 21, info)  
    if (info == 0) then  
       write(0, *) 'iolib: could not find marker:', mymarker  
       call mystop  
    end if
    read(21, *, end = 102) struc%nrk  
    allocate(struc%kndim(struc%nrk))  
    read(21, *, end = 102) (struc%kndim(i), i = 1, struc%nrk)  
    read(21, *, end = 102) struc%lda  
    allocate(struc%isrt(struc%lda, struc%nrk))  
    read(21, *, end = 102) ((struc%isrt(i, j), i = 1, struc%lda), &
         j = 1, struc%nrk)
    close(21)  

    return  

102 continue  
    write(0, *) 'iolib: reached eof before reading object:', mymarker
    call mystop  

  end subroutine read_isrt

end module vqmc_isrt_module
!
!     ==================== wavefn arrays structure ===============
!
module vqmc_wavefn_module  

  use constants
  use iolib_module
  implicit none

  type wavefn  
     sequence  

     integer :: nrk               ! number of k-points
     integer :: nspin             ! number of spins, 1 or 2
     integer :: mnband            ! maximum number of bands at all k-point
     integer :: lda               ! leading dimension of wfns array
     integer, pointer :: kndim(:)    ! # of plane-waves at each k-point
     integer, pointer :: &
          nband(:,:)                 ! # of bands computed (kpoint, spin)
     complex(dp), pointer :: &
          wfn(:,:,:,:)               ! the actual wfn array:
                                           ! wfn(lda, mnband, nrk, nspin)
  end type wavefn

contains

  subroutine destroy_wavefn(wf)

    implicit none
    type(wavefn), intent(inout) :: wf

    if (associated(wf%kndim)) deallocate(wf%kndim)  
    if (associated(wf%nband)) deallocate(wf%nband)  
    if (associated(wf%wfn)) deallocate(wf%wfn)  

  end subroutine destroy_wavefn

  subroutine write_wavefn(fname, struc)  

    implicit none      ! implicit? no!
    character(len=*), intent(in) :: fname      ! the filename
    type(wavefn), intent(in) :: struc

    !     ----------- local variables
    integer :: info, i, j, k, l  
    integer, parameter :: mymarker = 3

    open(unit = 21, file = fname, position = 'append', status = 'unknown', &
         form = 'formatted')
    write(21, *) mymarker, ' wavefn'  
    write(21, *) struc%nrk  
    write(21, *) (struc%kndim(i), i = 1, struc%nrk)  
    write(21, *) struc%nspin  
    write(21, *) struc%mnband  
    write(21, *) struc%lda  
    write(21, *) ((struc%nband(i, j), i = 1, struc%nrk), j = 1, struc%nspin)
    write(21, *) ((((struc%wfn(i, j, k, l), i = 1, struc%lda), &
         j = 1, struc%mnband), k = 1, struc%nrk), l = 1, struc%nspin)
    close(21)  
    return  

  end subroutine write_wavefn

  subroutine read_wavefn(fname, struc)  

    implicit none      ! implicit? no!
    character(len=*), intent(in) :: fname      ! the filename
    type(wavefn), intent(out) :: struc

    !     ----------- local variables
    integer :: info, i, j, k, l  
    integer, parameter :: mymarker = 3

    open(unit = 21, file = fname, status = 'unknown', form = 'formatted')
    call iolib_find_marker(mymarker, 21, info)  
    if (info == 0) then  
       write(0, *) 'iolib: could not find marker:', mymarker  
       call mystop  
    end if
    read(21, *, end = 102) struc%nrk  
    allocate(struc%kndim(struc%nrk))  
    read(21, *, end = 102) (struc%kndim(i), i = 1, struc%nrk)  
    read(21, *, end = 102) struc%nspin  
    read(21, *, end = 102) struc%mnband  
    read(21, *, end = 102) struc%lda  
    allocate(struc%nband(struc%nrk, struc%nspin))  
    read(21, *, end = 102) ((struc%nband(i, j), i = 1, struc%nrk), &
         j = 1, struc%nspin)
    allocate(struc%wfn(struc%lda, struc%mnband, struc%nrk, struc%nspin))
    read(21, *, end = 102) ((((struc%wfn(i, j, k, l), i = 1, struc%lda), &
         j = 1, struc%mnband), k = 1, struc%nrk), l = 1, struc%nspin)
    close(21)  

    return  

102 continue  
    write(0, *) 'iolib: reached eof before reading object:', mymarker
    call mystop

  end subroutine read_wavefn

end module vqmc_wavefn_module
!
!     ==================== level array structure ===============
!
module vqmc_level_module  

  use constants
  use iolib_module
  implicit none

  type level  
     sequence  

     real(dp) :: smearing,&         ! the smearing parameter used Ryd
                 fermilevel       ! the fermi level Ryd

     integer :: nrk               ! number of k-points
     integer :: mnband            ! maximum number of bands
     integer :: nspin             ! number of spins
     integer, pointer :: &
          nband(:,:)           ! number of bands computed (kpoint, spin)
     integer, pointer :: &
          ifmin(:,:)           ! lowest occupied band (kpoint, spin)
     integer, pointer :: &
          ifmax(:,:)           ! highest occupied band (kpoint, spin)
     real(dp), pointer :: &
          occup(:,:,:)         !occup number (band, kpoint, spin)
     real(dp), pointer :: &
          el(:,:,:)            ! energy Ryd (band, kpoint, spin)


  end type level

contains  

  subroutine destroy_level(lv)

    implicit none
    type(level), intent(inout) :: lv

    if (associated(lv%nband)) deallocate(lv%nband)  
    if (associated(lv%ifmin)) deallocate(lv%ifmin)  
    if (associated(lv%ifmax)) deallocate(lv%ifmax)  
    if (associated(lv%occup)) deallocate(lv%occup)  
    if (associated(lv%el)) deallocate(lv%el)  

  end subroutine destroy_level

  subroutine write_level(fname, struc)  

    implicit none      ! implicit? no!
    character(len=*), intent(in) :: fname      ! the filename
    type(level), intent(in) :: struc
  
    !     ----------- local variables
    integer :: info, i, j, k, l  
    integer, parameter :: mymarker = 5

    open(unit = 21, file = fname, position = 'append', status = 'unknown', &
         form = 'formatted')
    write(21, *) mymarker, ' level'  
    write(21, *) struc%nrk  
    write(21, *) struc%mnband  
    write(21, *) struc%nspin  
    write(21, *) ((struc%nband(i, j), i = 1, struc%nrk), j = 1, struc%nspin)
    write(21, *) ((struc%ifmin(i, j), i = 1, struc%nrk), j = 1, struc%nspin)
    write(21, *) ((struc%ifmax(i, j), i = 1, struc%nrk), j = 1, struc%nspin)
    write(21, *) (((struc%occup(k, i, j), k = 1, struc%mnband), &
         i = 1, struc%nrk), j = 1, struc%nspin)
    write(21, *) (((struc%el(k, i, j), k = 1, struc%mnband), &
         i = 1, struc%nrk), j = 1, struc%nspin)
    write(21, *) struc%smearing  
    write(21, *) struc%fermilevel  
    close(21)  
    return  

  end subroutine write_level

  subroutine read_level(fname, struc)  

    implicit none      ! implicit? no!
    character(len=*), intent(in) :: fname      ! the filename
    type(level), intent(out) :: struc
  
    !     ----------- local variables
    integer :: info, i, j, k, l  
    integer, parameter :: mymarker = 5

    open(unit = 21, file = fname, status = 'unknown', form = 'formatted')
    call iolib_find_marker(mymarker, 21, info)  
    if (info == 0) then  
       write(0, *) 'iolib: could not find marker:', mymarker  
       call mystop  
    end if
    read(21, *, end = 102) struc%nrk  
    read(21, *, end = 102) struc%mnband  
    read(21, *, end = 102) struc%nspin  
    allocate(struc%nband(struc%nrk, struc%nspin))  
    read(21, *, end = 102) ((struc%nband(i, j), i = 1, struc%nrk), &
         j = 1, struc%nspin)
    allocate(struc%ifmin(struc%nrk, struc%nspin))  
    allocate(struc%ifmax(struc%nrk, struc%nspin))  
    allocate(struc%occup(struc%mnband, struc%nrk, struc%nspin))  
    allocate(struc%el (struc%mnband, struc%nrk, struc%nspin) )  
    read(21, *, end = 102) ((struc%ifmin(i, j), i = 1, struc%nrk), &
         j = 1, struc%nspin)
    read(21, *, end = 102) ((struc%ifmax(i, j), i = 1, struc%nrk), &
         j = 1, struc%nspin)
    read(21, *, end = 102) (((struc%occup(k, i, j), k = 1, struc%mnband), &
         i = 1, struc%nrk), j = 1, struc%nspin)
    read(21, *, end = 102) (((struc%el(k, i, j), k = 1, struc%mnband), &
         i = 1, struc%nrk), j = 1, struc%nspin)
    read(21, *, end = 102) struc%smearing  
    read(21, *, end = 102) struc%fermilevel  
    close(21)  

    return  

102 continue  
    write(0, *) 'iolib: reached eof before reading object:', mymarker
    call mystop  

  end subroutine read_level

end module vqmc_level_module

!
!     ==================== molecular_dynamics structure =======================
!
module molecular_dynamics_module

  use constants
  implicit none

  type  molecular_dynamics 
 
    integer  itmax,&           ! # of MD steps
             ensemble,&        ! 1 NVE, 2 NVT
             n_deg_freedom,&   ! 3*natom_tot - 3
             extrapolation     ! extrapolation method

    real(dp) temp,&            ! inout temp at begin or for entire run
             time_step,&       ! time step for integration of eq. of motion
             Q_mass            ! mass of extended system for connstant T MD
             
    real(dp), pointer :: velocity(:) ! velocity of atoms 

  end type molecular_dynamics 


end module molecular_dynamics_module
