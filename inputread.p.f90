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
! Originally inputread.c
!
!       Service routines for reading the input file 
!       1995 Bernd Pfrommer
!
!       Parallel implementation by Greg McMullan 1999
!       Merged into existing code by Peter Haynes 1999
!
!       Converted to Fortran90 by Peter Haynes March 2000
!       and to use ESDF of Chris Pickard (see header of esdf_mod.f90)
!  
subroutine crystal_struc(myproc,mycomm,convergefile,switchlv,crys)
       
  use constants  
  use esdf   
  use crystal_module
  implicit none
include 'mpif.h'
  integer, intent(in) :: myproc, mycomm,convergefile

  type(crystal),intent(out) :: crys           ! crystal structure
  logical, intent(out) :: switchlv
  !
  !     -------------------------- local variables -----------------------
  !
  integer :: i,j,k,ntyp,mxatm,alen
  integer :: ierr
  integer :: natoms, nvec, nvol,nch
  real(dp) :: gvec(3, 3), volume, dummy, vec(3),vinit, matinvert,lvec(3,3)
  integer :: ip, ic
  integer :: ibuf(10)
  integer :: nlines, iline, symlen
  logical :: found
  character(len=llength) :: str  
  real(dp) :: getmass,temp,net_charge
  character, allocatable :: atm(:)
  real(dp), allocatable :: szet(:),mts(:),dipole(:),coord(:),mass(:)
  integer, allocatable :: natom(:)

  integer :: isupcell(3),isup,isupcell_tot(3),  info,scmax,scvmax
  !
  !     -------------------------- BEGIN  -----------------------
  !
  alen = 4
  ierr = 0  
  net_charge=0.d0

  if (myproc == 0) then

     ! check for supercell calculation  

     str = esdf_string('super_cell','1 1 1')
     read(str, *) isupcell(1:3)
     str = esdf_string('super_cell_vac', '0 0 0')
     read(str, *) isupcell_tot(1:3)
     isupcell_tot=isupcell_tot+isupcell
     scmax=max(isupcell_tot(1),isupcell_tot(2),isupcell_tot(3))

     if (esdf_defined('coordinates_absolute') .and. scmax .gt. 1) &
       call mystop( ' ''coordinates_absolute'' and ''super_cell'' are incompatible  ' )

     ! count the number of atomic types and coordinates per type

     found = esdf_block('coordinates',nlines)
     if (.not. found) call mystop( 'inputread:66 cannot find ''begin coordinates'' in input file'//trim(str) )
     ntyp = 0 ; mxatm = 0 ; natoms = 0
     do iline = 1, nlines
        if (index(block_data(iline), 'newtype') > 0) then  ! found new type
           ntyp = ntyp + 1
           natoms = 0
        else if (index(block_data(iline), 'coord') > 0) then ! found coordinate
          natoms = natoms +isupcell(1)*isupcell(2)*isupcell(3)
           mxatm = max(mxatm, natoms)
        else
           write(*, '(A/A)') 'found funny keyword in coordinate structure:', &
                block_data(iline)
        end if
     end do

     if (ntyp <= 0) call mystop( 'cannot find atomic types in input file' )
     if (mxatm <= 0) call mystop( 'cannot find coordinates in input file' )

     ibuf(1) = ntyp
     ibuf(2) = mxatm

  end if


  call MPI_Bcast(ibuf, 2, MPI_INTEGER, 0, mycomm, ierr)
  ntyp = ibuf(1)
  mxatm = ibuf(2)


  ! allocate arrays
  allocate(natom(ntyp),coord(3* mxatm*ntyp),szet(ntyp), mts(ntyp), &
       dipole(ntyp), mass(ntyp) )
  coord=dzero;

  ! NOTE: allow at least the size of an integer for each label
  allocate(atm(ntyp * alen))
  atm = repeat(' ', ntyp * alen)

  if (myproc == 0) then

     ! reposition to begin coordinates
     found = esdf_block('coordinates', nlines)
     if (.not. found) call mystop( 'inputread:108 cannot find ''begin coordinates'' in input file' )

     ! read the coordinates
     ic = 0
     do iline = 1, nlines
        if (index(block_data(iline), 'newtype') > 0) then    ! found new type
           ic = ic + 1
           natom(ic) = 0
           str = adjustl(block_data(iline) &
                (index(block_data(iline), 'newtype') + 7:))
           symlen = scan(str, ' ') - 1
           if (symlen > alen) write(*, '(A/A)') &
                'found funny symbol in coordinate structure:', &
                block_data(iline)
           symlen = min(symlen, alen)
           do i = 1, symlen	
              atm((ic - 1) * alen + i) = str(i:i)
           end do
           str = adjustl(str(scan(str, ' '):))
           mass(ic) = getmass(alen, atm((ic - 1) * alen + 1:ic * alen))
           szet(ic) = 0.2d0
           mts(ic) = done
           dipole(ic) = dzero
           if (len_trim(str) > 0) then                          ! found spin
              read(str, *) szet(ic)
              str = adjustl(str(scan(str, ' '):))
              if (len_trim(str) > 0) then          ! found muffin tin radius
                 read(str, *) mts(ic)
                 str = adjustl(str(scan(str, ' '):))
                 if (len_trim(str) > 0) then           ! found dipole moment
                    read(str, *) dipole(ic)
                 end if
              end if
           end if
        else if(index(block_data(iline), 'coord') > 0) then ! found coordinate
          str = adjustl(block_data(iline) &
                (index(block_data(iline), 'coord') + 5:))
           natom(ic) = natom(ic) +  1
           ip = 3 * (mxatm * (ic - 1) + natom(ic) - 1) + 1
           read(str, *) coord(ip), coord(ip+1), coord(ip+2)           
            coord(ip)= coord(ip)/isupcell_tot(1);
            coord(ip+1)= coord(ip+1)/isupcell_tot(2);
            coord(ip+2)= coord(ip+2)/isupcell_tot(3);
           do i=1,isupcell(1)
            do j=1,isupcell(2)
             do k=1,isupcell(3)
               if (i /=1 .or. j /=1 .or. k/=1 ) then
                 natom(ic) = natom(ic) +  1
                 isup = 3 * (mxatm * (ic - 1) + natom(ic) - 1) + 1    
                 coord(isup)=coord(ip)+real(i-1,dp)/isupcell_tot(1)
                 coord(isup+1)=coord(ip+1)+real(j-1,dp)/isupcell_tot(2)    
                 coord(isup+2)=coord(ip+2)+real(k-1,dp)/isupcell_tot(3)
               end if
             end do
            end do
          end do
        else
           write(*, '(A/A)') 'found funny keyword in coordinate structure:', &
                block_data(iline)
        end if
     end do

     if (ic /= ntyp) then
        write(*, '(2(A,I3))') 'ERROR: ic= ', ic, ' ntyp= ', ntyp
        call mystop
     end if

     do i=1,ntyp
        if (natom(i) < 1) then
           write(*, '(2A)') 'found no coordinates for type ', &
                atm((i - 1) * alen + 1:i * alen)
           call mystop
        end if

     end do

     ! --------------- read lattice vectors ----------------------

     found = esdf_block('latticevecs', nlines)
     if (.not. found) call mystop( 'cannot find ''begin latticevecs'' in input file' )
     nvec = 0 ; nvol = 0 ; nch=0
     do iline = 1, nlines
        if (index(block_data(iline), 'coord') > 0) then     ! found coordinate
           nvec = nvec + 1
        else if (index(block_data(iline), 'volume') > 0) then   ! found volume
           nvol = nvol + 1
        else if (index(block_data(iline), 'net_charge') > 0) then   ! found net_charge
           nch=nch+1
        else
           write(*, '(A/A)') 'found funny keyword in latticevecs structure:', &
                block_data(iline)
        end if
     end do

     if (nvec < 3) call mystop( 'cannot find all coordinates for lattice vectors' )
     if (nvec > 3) call mystop( 'too many coordinates for lattice vectors' )
     if (nvol > 1) call mystop( 'too many volumes for lattice vectors' )

     nvec = 0
     do iline = 1, nlines
        if (index(block_data(iline), 'coord') > 0) then    ! found coordinate
           str = adjustl(block_data(iline) &
                (index(block_data(iline), 'coord') + 5:))
           nvec = nvec + 1
           read(str, *) lvec(nvec, 1), lvec(nvec, 2), lvec(nvec, 3)
        else if (index(block_data(iline), 'volume') > 0) then  ! found volume
           str = adjustl(block_data(iline) &
                (index(block_data(iline), 'volume') + 6:))
           read(str, *) volume
        else if (index(block_data(iline), 'net_charge') > 0) then  ! found net_charge
           str = adjustl(block_data(iline) &
                (index(block_data(iline), 'net_charge') + 10:))
           read(str, *) net_charge
        else
           write(*, '(A/A)') 'found funny keyword in latticevecs structure:', &
                block_data(iline)
        end if
     end do

     gvec = lvec
     if (nvol == 0) then
        volume = abs(matinvert(gvec))
     else                     ! rescale coordinates
        volume = (volume / abs(matinvert(gvec)))**dthird
        lvec = lvec * volume
     end if

     lvec(1,1)=lvec(1,1)*isupcell_tot(1); lvec(1,2)=lvec(1,2)*isupcell_tot(1);
     lvec(1,3)=lvec(1,3)*isupcell_tot(1);
     lvec(2,1)=lvec(2,1)*isupcell_tot(2); lvec(2,2)=lvec(2,2)*isupcell_tot(2);
     lvec(2,3)=lvec(2,3)*isupcell_tot(2);
     lvec(3,1)=lvec(3,1)*isupcell_tot(3); lvec(3,2)=lvec(3,2)*isupcell_tot(3);
     lvec(3,3)=lvec(3,3)*isupcell_tot(3);

     ! -----calculate relative coordinates if given absolute -------

     if (esdf_defined('coordinates_absolute')) then
        gvec = transpose(lvec)
        dummy = matinvert(gvec)
        do i = 1, ntyp                         ! ----- loop over types -----
           do j = 1, natom(i)       ! ----- loop over atoms of this type---
              ip = 3 * (j - 1 + (i - 1) * mxatm) + 1
              vec(1:3) = coord(ip:ip + 2)
              coord(ip) = vec(1) * gvec(1, 1) + vec(2) * gvec(1, 2) + &
                   vec(3) * gvec(1, 3)
              coord(ip + 1) = vec(1) * gvec(2, 1) + vec(2) * gvec(2, 2) + &
                   vec(3) * gvec(2, 3)
              coord(ip + 2) = vec(1) * gvec(3, 1) + vec(2) * gvec(3, 2) + &
                   vec(3) * gvec(3, 3)
           end do
        end do
     end if

     gvec = lvec

     vinit = matinvert(gvec(1, 1))

     if (vinit < dzero) then

       write(convergefile, '(A)') &
           'Alas, your lattice vectors are not right handed'  
       write(convergefile, '(A)') &
            'I will swap two of them for you'

       switchlv=.true.
       do i=1,3
         temp=lvec(3,i)
         lvec(3,i)=lvec(2,i)
         lvec(2,i)=temp
       end do
       do i = 1, ntyp                     ! ----- loop over types -----
         do j = 1, natom(i)       ! ----- loop over atoms of this type---
           ip = 3 * (j - 1 + (i - 1) * mxatm) + 1
           temp=coord(ip+2)
           coord(ip+2)=coord(ip+1)
           coord(ip+1)=temp
         end do
       end do

       gvec=lvec  
       vinit=-vinit

    end if


  end if
  !
  ! Broadcast variables from proc 0 if in parallel
  !

  call MPI_Bcast(vinit, 1, MPI_DOUBLE_PRECISION, 0, mycomm, ierr)
  call MPI_Bcast(natom, ntyp, MPI_INTEGER, 0, mycomm, ierr)
  call MPI_Bcast(coord, 3 * mxatm * ntyp, MPI_DOUBLE_PRECISION, 0, &
       mycomm, ierr)
  call MPI_Bcast(szet, ntyp, MPI_DOUBLE_PRECISION, 0, mycomm, ierr)
  call MPI_Bcast(mts, ntyp, MPI_DOUBLE_PRECISION, 0, mycomm, ierr)
  call MPI_Bcast(dipole, ntyp, MPI_DOUBLE_PRECISION, 0, mycomm, ierr)
  call MPI_Bcast(mass, ntyp, MPI_DOUBLE_PRECISION, 0, mycomm, ierr)
  call MPI_Bcast(net_charge, 1, MPI_DOUBLE_PRECISION, 0, mycomm, ierr)
  call MPI_Bcast(lvec, 9, MPI_DOUBLE_PRECISION, 0, mycomm, ierr)
  call MPI_Bcast(atm, alen * ntyp, MPI_CHARACTER, 0, mycomm, ierr)

  ! 
  ! Setup crystal variable
  !
  allocate(crys%natom(ntyp),crys%coord(3, mxatm,ntyp),crys%szet(ntyp), &
           crys%mtsphere(ntyp),crys%dipole(ntyp), crys%mass(ntyp), &
           crys%nameat(ntyp),crys%rat(3, mxatm,ntyp),crys%zv(ntyp) )
  crys%coord=dzero; crys%rat=dzero;

  do i = 1, ntyp
     do j = 1, 2
        crys%nameat(i)(j:j) =atm((i - 1) * 4 + j)
     end do
  end do

  ip=0
  do i = 1, ntyp                     ! ----- loop over types -----
    do j = 1, natom(i)  
      do k=1,3
!        ip=ip+1
        ip = 3 * (j - 1 + (i - 1) * mxatm) + k !1
        crys%rat(k,j,i) = coord(ip)
        crys%coord(k,j,i) = coord(ip)
      end do
    end do
  end do       

  crys%ntype = ntyp
  crys%mxdatm = mxatm
  crys%vcell=vinit
  crys%avec = transpose(lvec)  
  crys%lvec = lvec  
  crys%natom = natom  
  crys%szet = szet  
  crys%mtsphere = mts  
  crys%dipole = dipole  
  crys%mass=mass
  crys%net_charge = net_charge


  deallocate(natom,coord,szet, mts, dipole, mass )

  return

end subroutine crystal_struc
!
! READ parameters for mostly relaxation and MD
!
subroutine read_param(myproc, mycomm, eps, fixed_vector, vector_plane,&
     decoup_rc, relax_method, MD,&
     relaxepsfac, relaxcoordfac, accuracy, relaxstring, p, &
     adjustpressure, jobstring, lambdalimit, relaxstress, relaxforce, &
     iupdateh, lambda_tolerance, lambda2_tolerance, itmax, avmass, vinit)
  !
  !     ---------------------------------------------------------------
  !
  use constants
  use esdf
  use molecular_dynamics_module
  implicit none
include 'mpif.h'
  !
  !     INPUT:
  !     ------
  !
  integer, intent(in) :: myproc, mycomm
  real(dp), intent(in) :: avmass, vinit
  !
  !     OUTPUT:
  !     ------
  !
  type(molecular_dynamics), intent(out) :: MD  ! molecular dynamics structure
  integer, intent(out) :: decoup_rc, adjustpressure, iupdateh, itmax
  real(dp), intent(out) :: eps(9), fixed_vector(3), vector_plane(3), &
         relaxepsfac, &
       relaxcoordfac, accuracy, p, lambdalimit, relaxstress, relaxforce, &
       lambda_tolerance, lambda2_tolerance
  integer, intent(out) :: relax_method
  character(len=llength), intent(out) :: relaxstring, jobstring
  !
  !     -------------------------- local variables -----------------------
  !
  integer :: ierr, icount, len, i, nlines
  character(len=llength) :: buf
  logical :: found
  real(dp) :: vector_length,bulkmodulus, optfreq
  integer, parameter :: logfile_c = 7

  integer :: ibuf(7)
  real(dp) :: dbuf(30)

  !
  !     -------------------------- BEGIN  -----------------------
  !

  if (myproc == 0) then

      ! --------------- read deformation tensor ----------------------
 
      found = esdf_block('deformation', nlines)
      if (found) then
         
         if (nlines .ne. 3) call mystop( 'cannot find all components of deformation tensor' )
         
         read(block_data(1),*) eps(1), eps(2), eps(3)
         read(block_data(2),*) eps(4), eps(5), eps(6)
         read(block_data(3),*) eps(7), eps(8), eps(9)
      else 
         eps = 0.0;
      end if
 
      buf = esdf_string('fixed_vector', '0.0 0.0 0.0')
      read(buf, *) fixed_vector(1:3)
      buf = esdf_string('other_vector_in_fixed_plane', '0.0 0.0 0.0')
      read(buf, *) vector_plane(1:3)
 
      if (vector_length(fixed_vector) .ne. dzero) then
         fixed_vector = fixed_vector/vector_length(fixed_vector)
      end if
      if (vector_length(vector_plane) .ne. dzero) then
         vector_plane = vector_plane/vector_length(vector_plane)
      end if
      
      ! End of funky new stuff by weidong and david

     decoup_rc = 0
     if (index(esdf_string('decouple_coordinates', 'no'), 'yes') > 0) &
          decoup_rc = 1

     bulkmodulus = esdf_physical('estimated_bulk_modulus', 1.0d0, 'GPa')

     if (bulkmodulus <= 1.0d-6) then
        write(logfile_c, '(A,F6.2,A)') &
             'estimated_bulk_modulus too small: ', bulkmodulus, ' MBar'
     end if

     relaxepsfac = esdf_double('relax_eps_factor', &
          147.1076d0 / (dthree * bulkmodulus * vinit))

     optfreq = esdf_physical('estimated_optical_phonon_frequency', &
          15.0d0, 'THz')
     relaxcoordfac = esdf_double('relax_coord_factor', &
          done / (avmass * optfreq * optfreq * 8.4157d-5))

     accuracy = esdf_double('relax_accuracy', 0.01d0)

     buf = esdf_string('relax_pressure', '0.0')
     if (index(buf,'adjust') > 0) then
        p = 0.0d0
        adjustpressure = 1
     else
        read(buf, *) p
        adjustpressure = 0
     end if

     jobstring = esdf_string('job', 'relax')

     lambdalimit = esdf_double('lambda_limit', 10.d0)

     relax_method = esdf_integer('relax_method', 1)

     relaxforce = dzero ; relaxstress = dzero
     relaxstring = esdf_string('relax_what', 'force_and_stress')
     if (index(relaxstring, 'force') > 0) relaxforce = done
     if (index(relaxstring, 'stress') > 0) relaxstress = done

     if (relax_method .eq. 1) then
       relaxstring = esdf_string('relax_how', 'normal')
     else
       relaxstring = esdf_string('relax_how', 'fast')
     end if

     lambda_tolerance  = 0.1d0 ; lambda2_tolerance = 0.1d0
     iupdateh =	1
     if (index(relaxstring, 'gradient') > 0) iupdateh = 0
     if (index(relaxstring, 'slow') > 0) then
        lambda_tolerance = 1.0d-6 ; lambda2_tolerance = 1.0d-6
     end if
     if (index(relaxstring, 'normal') > 0) then
        lambda_tolerance = 0.6d0 ; lambda2_tolerance = 0.3d0
     end if
     if (index(relaxstring, 'fast') > 0) then
        lambda_tolerance = 10.0d0 ; lambda2_tolerance = dtwo
     end if

     itmax = esdf_integer('relax_max_iter', 0)
   
     MD%itmax = esdf_integer('MD_max_iter', 100)    
     MD%itmax = max(itmax,MD%itmax)

     MD%ensemble = esdf_integer('ensemble_type', 2)
     MD%Q_mass =  esdf_double('MD_Q_mass', 0d0)
     MD%time_step = esdf_double('MD_time_step', 3d0)
     MD%temp = esdf_double('MD_temp', 500d0)

    if (trim(jobstring) == 'MD') then
      relaxstress=dzero
      MD%extrapolation = esdf_integer('extrapolation_method', 2)
    else 
      MD%extrapolation = esdf_integer('extrapolation_method', 1)
    end if

    if (relaxstress .eq. done) then
      MD%extrapolation=-1
    end if

  end if



  if (myproc == 0) then
     ibuf(1) = decoup_rc
     ibuf(2) = iupdateh
     ibuf(3) = adjustpressure
     ibuf(4) = itmax
     ibuf(5) = MD%itmax
     ibuf(6) = MD%extrapolation
     ibuf(7) = MD%ensemble
  end if
  call MPI_Bcast(ibuf(1), 7, MPI_INTEGER, 0, mycomm, ierr)
  if (myproc /= 0) then
     decoup_rc = ibuf(1)
     iupdateh = ibuf(2)
     adjustpressure = ibuf(3)
     itmax = ibuf(4)
     MD%itmax = ibuf(5)
     MD%extrapolation = ibuf(6)
     MD%ensemble = ibuf(7)
  end if

  call MPI_Bcast(relaxstring, llength, MPI_CHARACTER, 0, mycomm, ierr)
  call MPI_Bcast(jobstring, llength, MPI_CHARACTER, 0, mycomm, ierr)

  if (myproc == 0) then
     dbuf(1:9) = eps(1:9)
     dbuf(10) = bulkmodulus
     dbuf(11) = relaxepsfac
     dbuf(12) = optfreq
     dbuf(13) = relaxcoordfac
     dbuf(14) = accuracy
     dbuf(15) = p
     dbuf(16) = lambdalimit
     dbuf(17) = relaxstress
     dbuf(18) = relaxforce
     dbuf(19) = lambda_tolerance
     dbuf(20) = lambda2_tolerance
     dbuf(21:23) = fixed_vector(1:3)
     dbuf(24:26) = vector_plane(1:3)
     dbuf(27) = relax_method
     dbuf(28) = MD%time_step
     dbuf(29) = MD%temp
     dbuf(30) = MD%Q_mass
  end if
  call MPI_Bcast(dbuf(1), 30, MPI_DOUBLE_PRECISION, 0, mycomm, ierr)
  if (myproc /= 0) then
     eps(1:9) = dbuf(1:9)
     bulkmodulus = dbuf(10)
     relaxepsfac = dbuf(11)
     optfreq = dbuf(12)
     relaxcoordfac = dbuf(13)
     accuracy = dbuf(14)
     p = dbuf(15)
     lambdalimit = dbuf(16)
     relaxstress = dbuf(17)
     relaxforce = dbuf(18)
     lambda_tolerance = dbuf(19)
     lambda2_tolerance = dbuf(20)
     fixed_vector(1:3) = dbuf(21:23)
     vector_plane(1:3) = dbuf(24:26)
     relax_method = dbuf(27) 
     MD%time_step = dbuf(28) 
     MD%temp = dbuf(29) 
     MD%Q_mass =  dbuf(30)
  end if



end subroutine read_param
!
! read all the parameters for the planewave code
!
subroutine read_pwparam(myproc, mycomm,pw_params,altkpoints,bands,&
        altbands,syms,energs,crys,lvec)
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
  use esdf
  implicit none
include 'mpif.h'
  !
  !     INPUT:
  !     ------
  !
  integer, intent(in) :: myproc, mycomm
  !
  !     INPUT/OUTPUT:
  !     ------------
  !
  type(crystal) :: crys     ! crystal structure
  !
  !     OUTPUT:
  !     ------
  !
  type(pw_parameter) :: pw_params ! plane wave parameters
  type(kpoint) :: altkpoints  ! alternative BZ integration kpoints
  type(band) :: bands ! eigenvalues, occupation numbers etc
  type(band) :: altbands    ! another set of eigenvalues etc..
  type(symmetry) :: syms ! symmetry operations
  type(energy) :: energs    ! energies and other information
  !
  !     -------------------------- local variables -----------------------
  !
  real(dp) :: lvec(3, 3)               ! lattice vectors
  integer :: ierr, istop, nwv, nline, ifound, i, njob,nvec
  integer :: itmp, icount, nitems, j
  integer :: nlines, iline
  character(len=llength) :: buf
  integer :: generate_kpts
  logical found
  character(len=llength) :: str

  ! always update pwmain string array also. only add at the end !

! davegp
!  integer, parameter :: NUM_PW_JOBS = 6
  integer, parameter :: NUM_PW_JOBS = 7

!  character(len=14), parameter :: pw_jobstrings(NUM_PW_JOBS) = (/ &
!       'scf           ', &
!       'nmr_shift     ', &
!       'nmr_plot      ', &
!       'band_structure', &
!       'pot_plot      ', &
!       'tddft         ' /)
  character(len=14), parameter :: pw_jobstrings(NUM_PW_JOBS) = (/ &
       'scf           ', &
       'nmr_shift     ', &
       'nmr_plot      ', &
       'band_structure', &
       'pot_plot      ', &
       'tddft         ', &
       'nonselfcon    ' /)
! davegp



! LSDA+U. P. Zhang
!  integer :: ibuf(53)
  integer :: ibuf(57)
! Weiwei Gao for gw_sum
  integer :: ibuf_gw_sum
  real(dp) :: dbuf(48)
  character(42) :: cbuf

  !
  !     -------------------------- BEGIN  -----------------------
  !
  allocate(pw_params%pp_iloc(crys%ntype),pw_params%pp_pwfocc(3,crys%ntype))
  allocate(pw_params%NLPP_rcut(crys%ntype))

  if (myproc == 0) then

     altkpoints%nrk = esdf_integer('number_kpoints', 0)
     if (altkpoints%nrk <= 0) then
        buf = esdf_string('k_grid', '1 1 1')
        read(buf, *) altkpoints%grid(1:3)
        buf = esdf_string('k_grid_shift', '0.0 0.0 0.0')
        read(buf, *) altkpoints%shift(1:3)
     end if

     buf = esdf_string('gw_shift', '0.0 0.0 0.0')
     read(buf, *) pw_params%q(1:3)

     pw_params%smearing = esdf_physical('gaussian_smearing', 0.05d0, 'eV')

     if ( pw_params%smearing .eq. 0.05d0) &
       pw_params%smearing = esdf_physical('smearing_energy', 0.05d0, 'eV') 

     pw_params%smearing_method = esdf_integer('smearing_method', 1)

     pw_params%epsdiag = esdf_double('accuracy_diag', 1.0d-12)

     pw_params%init_dir_econv = esdf_double('init_dir_econv', 1.0d-1) 

     pw_params%shiftsafety = esdf_physical('diagsafety', 1.0d0, 'eV')

     pw_params%random_startvec =&
                         esdf_double('randomize_diag_start_guess', 0.0d0)

     pw_params%bandgap = esdf_physical('bandgap', 0.0d0, 'eV')/ryd

     pw_params%emax = esdf_physical('energy_cutoff', 20.0d0, 'Ry')

     pw_params%emaxsub =  esdf_physical('submatrix_energy_cutoff', 5.0d0, 'Ry')

! PZ  full exact diagonalization
     pw_params%iexact_diag=0
     if(abs((pw_params%emax-pw_params%emaxsub)/pw_params%emax).le.0.01d0) then
     pw_params%iexact_diag=1
     write(9,*) "Exact diagonalization is employed."
     end if

     pw_params%emaxmix = esdf_physical('mixing_energy_cutoff', 5.0d0, 'Ry')

     ! Determine the method of mixing
     pw_params%mix_method = esdf_string('mix_method','pulay_kerker')

     pw_params%numpk_befptf = esdf_integer('PTF_num_init_PK', 2) 

     pw_params%tfw_max_cg = esdf_integer('PTF_max_iter_cg', 40) 

     pw_params%pp_format = esdf_integer('pp_format',1)

     pw_params%NL_rspace(1) = esdf_boolean('NLPP_rspace',.false.)
     pw_params%NL_rspace(2) = esdf_boolean('NLPP_rspace_force',.false.)

     pw_params%io_scf = esdf_boolean('io_scf',.false.)     

     if (.not. pw_params%NL_rspace(1) .and. pw_params%NL_rspace(2) ) then
       write(9,*) ' cannot do real space force without real-space energy'
       call mystop
     end if        

     if(pw_params%NL_rspace(1) .and. pw_params%pp_format .eq. 1) then
       write(9,*) ' cannot do real space NL with pp_format 1'
       call mystop
     end if   

     if(pw_params%NL_rspace(1)) then
       str = esdf_string('NLPP_rcut',' 4.7')
       read(str,*) pw_params%NLPP_rcut(1)
 
       do j=1,3
         if ( 2*pw_params%NLPP_rcut(1) .gt. sqrt(lvec(j,1)**2 &
            + lvec(j,2)**2 + lvec(j,3)**2) ) then
           write(6,*) '2*rcut=',2*pw_params%NLPP_rcut(1),' > the',j,&
                 '  a lattive vector: cannot use real space non-local',&
             sqrt(lvec(j,1)**2 + lvec(j,2)**2 + lvec(j,3)**2)      
           call mystop 
         end if
       end do
 
       do i=2,crys%ntype
         pw_params%NLPP_rcut(i)=pw_params%NLPP_rcut(1)
       end do
       do i=2,crys%ntype
         str = adjustl(str(scan(str, ' '):))
         if (len_trim(str) > 0) then 
           read(str,*) pw_params%NLPP_rcut(i)
           do j=1,3
             if ( 2*pw_params%NLPP_rcut(i) .gt. sqrt(lvec(j,1)**2 &
                  + lvec(j,2)**2 + lvec(j,3)**2) ) then
               write(6,*) '2*rcut=',2*pw_params%NLPP_rcut(i),' > the',j,&
                 '  a lattive vector: cannot use real space non-local'   
               call mystop 
             end if
           end do
         else
           go to 100
         end if
       end do
     end if
 
100  continue

     if (pw_params%mix_method .eq. 'pulay_tf' .and. &
          pw_params%pp_format .eq. 1) &
          call mystop( ' pulay_tf mixing incompatible with pp_format 1' )

     if(pw_params%pp_format .eq. 3) then
       found = esdf_block('pseudopotential', nlines)
       if (.not. found) &
              call mystop( 'cannot find ''begin pseudopotential'' in input file' )

       nvec=0
       do iline = 1, crys%ntype
         if (index(block_data(iline), 'pp_data') > 0) then  ! found coordinate
          nvec = nvec + 1
        else
          write(*, '(A/A)') 'funny keyword in pseudopotenital structure:', &
          block_data(iline)
        end if
     end do
     
     if (nvec < crys%ntype) call mystop( 'cannot find all data for pseudopotentials' )
     if (nvec > crys%ntype) call mystop( 'too much data for pseudopotenitals' )

     nvec = 0
     do iline = 1,crys%ntype
        if (index(block_data(iline), 'pp_data') > 0) then    ! found coordinate
           str = adjustl(block_data(iline) &
                (index(block_data(iline), 'pp_data') + 7:))
           nvec = nvec + 1
           pw_params%pp_pwfocc(:,nvec)=-1d0 ! for checking when reading pseudopotential file
           if (len_trim(str) > 0) then                          ! found local index
              read(str, *) pw_params%pp_iloc(nvec)
              str = adjustl(str(scan(str, ' '):))
              if (len_trim(str) > 0) then          ! found occupation of s projector
                 read(str, *) pw_params%pp_pwfocc(1,nvec)
                 itmp=1
                 str = adjustl(str(scan(str, ' '):))
                 if (len_trim(str) > 0) then           ! found occupation of p projector
                    read(str, *) pw_params%pp_pwfocc(2,nvec)
                    itmp=2
                    str = adjustl(str(scan(str, ' '):))
                    if (len_trim(str) > 0) then           ! found occupation of d projector
                       read(str, *) pw_params%pp_pwfocc(3,nvec)
                       itmp=3
                    end if
                 end if
              end if
           end if
           if (pw_params%pp_iloc(nvec) .le. 0 .or. &
                pw_params%pp_iloc(nvec) .gt. itmp ) call mystop( ' iloc for atom nvec is bad' )
        else
           write(*, '(A/A)') 'funny keyword in pseudopotential structure:', &
                block_data(iline)
        end if
     end do

     end if

     ! corrected kinetic energy functional
  
     buf = esdf_string('modify_kinetic_energy', 'off')
     if (index(buf, 'on') > 0) then
        pw_params%ekinmod(1) = pw_params%emax * dtwo
        pw_params%ekinmod(2) =  pw_params%emax
        pw_params%ekinmod(3) =  pw_params%emax * 0.1d0
     else if (index(buf, 'off') > 0) then
        pw_params%ekinmod(1) = dzero
        pw_params%ekinmod(2) = 1.0d3
        pw_params%ekinmod(3)  = done
     else
        read(buf, *) pw_params%ekinmod(1:3)
     end if

     pw_params%nrg = esdf_integer('number_radial_gridpoints', 10)

     pw_params%nbandsfft = esdf_integer('number_bands_fft', 1)


     pw_params%correct_stress_cutoff = &
          esdf_physical('polished_energy_cutoff', -1.0d0, 'Ry')

! LSDA+U option; default: do L(S)DA. P. Zhang

     pw_params%lsdau= esdf_integer('lsda_u', 0)
     pw_params%wannier= esdf_integer('wannier', 0)
!     pw_params%ipr = esdf_integer('ipr', 0) !bshih20100409

! GW band index reordering, P. Zhang

     bands%gw_band_reordering = esdf_integer('gw_band_reordering', 0)

! LSDA+U must be done with LSDA. P. Zhang
     crys%nspin = esdf_integer('number_of_spins', 1)

!Weiwei Gao: this flag is for wavefunctions for GW with energy integration
     pw_params%gw_sum = esdf_integer('gw_sum', 0) 


! allow LDA+U
!     if (pw_params%lsdau.eq.1) crys%nspin=2

     crys%nalpha = esdf_integer('number_of_alpha', -1)

     crys%nbeta = esdf_integer('number_of_beta', -1)  !F.Mauri

     pw_params%maxitdiag = esdf_integer('max_iter_diag', 3)
     if (pw_params%maxitdiag .eq. 3) then
        pw_params%maxitdiag = esdf_integer('min_iter_diag', 3)        
     end if

     pw_params%maxitdiag_band = esdf_integer('max_iter_diag_band', 137)
     if (pw_params%maxitdiag > pw_params%maxitdiag_band) &
          pw_params%maxitdiag_band = pw_params%maxitdiag
     pw_params%maxitscf = esdf_integer('max_iter_scfloop', 40)

     pw_params%itdiag_add = esdf_integer('iter_diag_add', 5)
     pw_params%itdiag_add_metal = esdf_integer('iter_diag_add_metal', 3)


     pw_params%nprocpvm = esdf_integer('nproc_pvm', 1)

     buf = esdf_string('screening_type', 'previous')
     if (index(buf, 'previous') > 0) then
        pw_params%iscr = 4
     else if (index(buf, 'atomic') > 0) then
        pw_params%iscr = 5
     else if (index(buf, 'vmix') > 0) then
        pw_params%iscr = 6
     else
        write(*,'(3A)') 'unknown screening type: ', trim(buf), &
             '; using ''previous'''
        pw_params%iscr = 4
     end if

     pw_params%output(:) = 0
     buf = esdf_string('output_flags', '')
     if (index(buf, 'diagperf') > 0) pw_params%output(1) = pw_params%output(1) + 1
     if (index(buf, 'vqmc') > 0) pw_params%output(1) = pw_params%output(1) + 2
     if (index(buf, 'nmrshift') > 0) pw_params%output(1) = pw_params%output(1) + 4
     if (index(buf, 'timing') > 0) pw_params%output(1) = pw_params%output(1) + 8
     if (index(buf, 'cdplot') > 0) pw_params%output(1) = pw_params%output(1) + 16
     if (index(buf, 'ballnstick') > 0) pw_params%output(1) = pw_params%output(1) + 32
     if (index(buf, 'neighbor') > 0) pw_params%output(1) = pw_params%output(1) + 64
     if (index(buf, 'dos') > 0) pw_params%output(1) = pw_params%output(1) + 128
     if (index(buf, 'waveplot') > 0) pw_params%output(1) = pw_params%output(1) + 256
     if (index(buf, 'wavefn') > 0) pw_params%output(1) = pw_params%output(1) + 512
     if (index(buf, 'angdos') > 0) pw_params%output(1) = pw_params%output(1) + 1024
     if (index(buf, 'eigval') > 0) pw_params%output(1) = pw_params%output(1) + 2048
     if (index(buf, 'momden') > 0) pw_params%output(1) = pw_params%output(1) + 4096
     if (index(buf, 'angle') > 0) pw_params%output(1) = pw_params%output(1) + 8192
     if (index(buf, 'ballswrapped') > 0) pw_params%output(1) = pw_params%output(1) + 16384
     if (index(buf, 'khoros') > 0) pw_params%output(1) = pw_params%output(1) + 32768
     if (index(buf, 'memusage') > 0) pw_params%output(1) = pw_params%output(1) + 65536
     if (index(buf, 'efield') > 0) pw_params%output(1) = pw_params%output(1) + 131072
     if (index(buf, 'potplot') > 0) pw_params%output(1) = pw_params%output(1) + 262144
     if (index(buf, 'gwr') > 0) pw_params%output(1) = pw_params%output(1) + 524288
     if (index(buf, 'gwc') > 0) pw_params%output(1) = pw_params%output(1) + 1048576
     if (index(buf, 'gwscreening') > 0) pw_params%output(1) = pw_params%output(1) + 2097152
     if (index(buf, 'project') > 0) pw_params%output(1) = pw_params%output(1) + 4194304
     if (index(buf, 'eqn3') > 0) pw_params%output(1) = pw_params%output(1) + 8388608
     if (index(buf, 'eqn8') > 0) pw_params%output(1) = pw_params%output(1) + 16777216
     if (index(buf, 'crystal') > 0) pw_params%output(1) = pw_params%output(1) + 33554432
     if (index(buf, 'original') > 0) pw_params%output(1) = pw_params%output(1) + 67108864
     if (index(buf, 'xyz') > 0) pw_params%output(1) = pw_params%output(1) + 134217728
     if (index(buf, 'appendout') > 0) pw_params%output(1) = pw_params%output(1) + 268435456
     if (index(buf, 'diagconv') > 0) pw_params%output(1) = pw_params%output(1) + 536870912
     if (index(buf, 'occ') > 0) pw_params%output(1) = pw_params%output(1) + 1073741824


     if (index(buf, 'electronphonon') > 0) &          ! added by hjchoi
          pw_params%output(4) = pw_params%output(4) + 2 ! added by hjchoi
     if (index(buf, 'hamwavefn') > 0) &               ! added by hjchoi
          pw_params%output(4) = pw_params%output(4) + 1 ! added by hjchoi
     if (index(buf, 'fermisurf') > 0) &
          pw_params%output(4) = pw_params%output(4) + 4
     if (index(buf, 'nofrcstr') > 0) &
          pw_params%output(4) = pw_params%output(4) + 8
     if (index(buf, 'noioscf') > 0) &
          pw_params%output(4) = pw_params%output(4) + 16
     !
     ! check if smearing methods are compatible with optimize metal
     !
     if (iand(pw_params%optimize,64) == 64 ) then
       if (pw_params%smearing_method .ne. 2 ) then
         write(9,*) 'OPTIMIZE METAL can only be used with smearing 2 (FD)' 
         write(9,*) 'switching to smearing 2'
         pw_params%smearing_method=2
       end if
     end if

     pw_params%miscflag = 0
     buf = esdf_string('no', '')
     if (index(buf, 'restore_inversion') > 0) pw_params%miscflag = &
                                                pw_params%miscflag + 1
     if (index(buf, 'submatrix_diag') > 0) pw_params%miscflag = &
                                                pw_params%miscflag + 4
     buf = esdf_string('checkpoint', '')
     if (index(buf, 'nmr') > 0) pw_params%miscflag = pw_params%miscflag + 2



     pw_params%nmrkpts = esdf_integer('nmrkpts',huge(1))

     pw_params%checkpoint_wfn  = esdf_integer('checkpoint_emin', 1000000)

     buf = esdf_string('exchange_correlation', 'ceperley_alder')

     if (index(buf, 'ceperley_alder') > 0 .or. index(buf, 'ca') > 0) then
        pw_params%icorr = 'ca'
     else if (index(buf,'perdew_wang_91') > 0 .or. index(buf, 'pw91') > 0) then
        pw_params%icorr = 'pw'
     else if (index(buf, 'perdew_burke_ernzerhof') > 0 .or. &
          index(buf, 'pbe') > 0) then
        pw_params%icorr = 'pb'
     else
        write(*, '(3A)') 'unknown exchange-correlation type: ', trim(buf), &
             '; using ceperley-alder'
        pw_params%icorr = 'ca'
     end if

     pw_params%input = 0
     buf = esdf_string('input_flags', '')
     if (index(buf, 'wavefn') > 0) pw_params%input = pw_params%input + 1
     if (index(buf, 'lukman') > 0) pw_params%input = pw_params%input + 2
     if (index(buf, 'chkptchi') > 0) pw_params%input = pw_params%input + 4

     pw_params%optimize = 2 + 4 + 16 + 32 + 256
     buf = esdf_string('optimize', '')
     if (index(buf, 'memory') > 0)  pw_params%optimize = pw_params%optimize + 1
     if (index(buf, 'nostartguess') > 0)  pw_params%optimize = pw_params%optimize - 2
     if (index(buf, 'nocgsolveguess') > 0)  pw_params%optimize = pw_params%optimize - 4
     if (index(buf, 'insulator') > 0)  pw_params%optimize = pw_params%optimize + 8
     if (index(buf, 'nowfnreuse') > 0)  pw_params%optimize = pw_params%optimize - 16
     if (index(buf, 'noprecond') > 0)  pw_params%optimize = pw_params%optimize - 32
     if (index(buf, 'metal') > 0)  pw_params%optimize = pw_params%optimize + 64
     if (index(buf, 'descend') > 0)  pw_params%optimize = pw_params%optimize + 128
     if (index(buf, 'noadjustshift') > 0)  pw_params%optimize = pw_params%optimize - 256

 
     if (iand(pw_params%optimize,8) == 8 ) then
       pw_params%diag_method=esdf_string('diagonalization_method','Grassmann') 
     else
       pw_params%diag_method=esdf_string&
                         ('diagonalization_method','Grassmann_metal')
     end if

     bands%min(1) = esdf_integer('number_bands', -1)
     if (bands%min(1) < 0) call mystop( 'cannot find ''number_bands'' in input file' )
     bands%min(2)=bands%min(1)

     if (esdf_defined('accuracy_diag_per_band')) then
         pw_params%epsdiag =  pw_params%epsdiag*bands%min(1) 
     end if

     altbands%min(1) = esdf_integer('number_alt_bands', bands%min(1)   )
     altbands%min(2) = altbands%min(1) 

     pw_params%nband_tddft = esdf_integer('number_bands_tddft', bands%min(1) )


     if (iand(pw_params%output(1), 1048576) == 1048576  .or. &
          (iand(pw_params%output(1), 524288) == 524288)) then
         bands%num_gwout = esdf_integer('num_gwout', -1)
         bands%gw_mid_energy= esdf_double('gwout_mid_energy', 1d4)
         bands%gw_low_energy=esdf_double('gwout_low_energy', 1d4)
         bands%gw_high_energy=esdf_double('gwout_high_energy', 1d4)
     end if
  
     buf = esdf_string('occupy_levels', 'normal')
     if (index(buf, 'normal') > 0) then
       pw_params%occupy = 0
     else if (index(buf, 'fixed') > 0) then
       pw_params%occupy = 1
     else if (index(buf, 'antiferro') > 0) then
       pw_params%occupy = 2
     else if (index(buf, 'excited') > 0) then
       if ( crys%nspin  > 1) &
             call mystop( 'Cannot use excited occupations with spin polarization' )
       pw_params%occupy = 3
     else if (index(buf, 'from_input') > 0) then
       if ( crys%nalpha < 0 ) &
             call mystop( 'Please specify number_of_alpha' )
       if ( crys%nbeta  < 0 ) &
             call mystop( 'Please specify number_of_beta' )

       bands%min(1)= esdf_integer('number_of_bands_alpha', crys%nalpha)
       bands%min(2)= esdf_integer('number_of_bands_beta', crys%nbeta)  !F.Mauri

       if (iand(pw_params%optimize,8) == 8 ) then
         bands%min(1)=crys%nalpha
         bands%min(2)=crys%nbeta
       end if
       
        pw_params%smearing = 0.d0 ! NO SMEARING FOR FIXED OCCUPATIONS
        pw_params%occupy = 4  !F. Mauri

! P. Zhang, for calculating bare phonon energy
     else if (index(buf, 'from_file') > 0) then
       write(9,*) "occupy levels: fixed_from_input_file"
       pw_params%occupy = 5
     else
        write(*, '(3A)') 'unknown occupy_levels type: ', trim(buf), &
             '; using ''normal'''
        pw_params%occupy = 0
     end if

     energs%inputfermi = esdf_physical('fermi_level', 0.0d0, 'Ry')
     energs%evalshift = dzero

     pw_params%alphamix(1) = esdf_double('mixing_parameter', 0.77d0)
     buf = esdf_string('linear_mixing', '0.0 0.6 0.05')
     read(buf, *) pw_params%alphamix(2:4)

     pw_params%epscv = esdf_double('potential_convergence_criterion', 1.0d-6)

     pw_params%epsce = esdf_double('energy_convergence_criterion', 1.0d-9)

     pw_params%chi_q = esdf_double('nmr_q', 0.01d0)

     pw_params%bd_lgth = esdf_double('bond_length', 4.0d0) 

     buf = esdf_string('nmr_g0mask', '0.6666666666666667 0.0 0.0  &
          & 0.0 0.6666666666666667 0.0  &
          & 0.0 0.0 0.6666666666666667')
     read(buf, *) pw_params%g0mask(1:9)

     syms%ntrans = esdf_integer('number_symmetry_ops', -1)

     buf = esdf_string('magnetic_field_kvec', '0 0 0 0')
     read(buf, *) pw_params%mfield%gvec(1:4)

     pw_params%mfield%h = esdf_double('magnetic_field_strength', 0.0d0)

     buf = esdf_string('electric_potential_kvec', '0 0 0 0')
     read(buf, *) pw_params%epot%gvec(1:4)

     pw_params%epot%h = esdf_double('electric_potential_strength', 0.0d0)

     ! --------------- read the plot_wave_function information ---------

     if (esdf_block('plot_wave_function', nlines)) then
        nwv = 0
        do iline = 1, nlines
           if (index(block_data(iline), 'wavefunction') > 0) then
              nwv = nwv + 1
           else
              write(*, '(A/A)') &
                   'found funny keyword in plot_wave_function structure:', &
                   block_data(iline)
           end if
        end do
        if (nwv == 0) &
             call mystop( 'found no wavefunctions in plot_wave_function structure' )
        pw_params%nplotwv = nwv
        allocate(pw_params%plotwv(2 * nwv))
        nwv = 0
        do iline = 1, nlines
           if (index(block_data(iline), 'wavefunction') > 0) then
              nwv = nwv + 1
              buf = adjustl(block_data(iline) &
                   (index(block_data(iline), 'wavefunction') + 12:))
              read(buf, *) pw_params%plotwv((nwv - 1) * 2 + 1:nwv * 2)
           else
              write(*, '(A/A)') &
                   'found funny keyword in plot_wave_function structure:', &
                   block_data(iline)
           end if
        end do
     else
        pw_params%nplotwv = 0
        allocate(pw_params%plotwv(1))
     end if

     ! --------------- read the line_plot information ---------

     buf = esdf_string('tensor_bra', '1.0 0.0 0.0')
     read(buf, *) pw_params%bra(1:3)

     buf = esdf_string('tensor_ket', '1.0 0.0 0.0')
     read(buf, *) pw_params%ket(1:3)

     if (esdf_block('line_plot', nlines)) then
        nline = 0
        do iline = 1, nlines
           if (index(block_data(iline), 'line') > 0) then
              nline = nline + 1
           else
              write(*, '(A/A)') &
                   'found funny keyword in line_plot structure:', &
                   block_data(iline)
           end if
        end do
        if (nline == 0) call mystop( 'found no lines in line_plot structure' )
        pw_params%nlineplot = nline
        allocate(pw_params%lineplot(7 * nline))
        nline = 0
        do iline = 1, nlines
           if (index(block_data(iline), 'line') > 0) then
              nline = nline + 1
              buf = adjustl(block_data(iline) &
                   (index(block_data(iline), 'line') + 4:))
              read(buf, *) pw_params%lineplot((nline - 1) * 7 + 1:nline * 7)
           else
              write(*, '(A/A)') &
                   'found funny keyword in line_plot structure:', &
                   block_data(iline)
           end if
        end do
     else
        pw_params%nlineplot = 0
        allocate(pw_params%lineplot(1))
     end if

      !        -------- read the slice plot information
   
     if (esdf_block('slice_plot', nlines)) then
         nline = 0
         do iline = 1, nlines
            if (index(block_data(iline), 'slice') > 0) then
               nline = nline + 1
            else
               write(*, '(A/A)') &
                    'found funny keyword in slice_plot structure:', &
                    block_data(iline)
            end if
         end do
         if (nline == 0) call mystop( 'found no lines in slice_plot structure' )
         pw_params%nsliceplot = nline
         allocate(pw_params%sliceplot(11 * nline))
         nline = 0
         do iline = 1, nlines
            if (index(block_data(iline), 'slice') > 0) then
               nline = nline + 1
               buf = adjustl(block_data(iline) &
                    (index(block_data(iline), 'slice') + 5:))
               read(buf, *) pw_params%sliceplot((nline - 1) * 11 + 1:nline * 11)
            else
               write(*, '(A/A)') &
                    'found funny keyword in slice_plot structure:', &
                    block_data(iline)
            end if
         end do
      else
         pw_params%nsliceplot = 0
         allocate(pw_params%sliceplot(1))
      end if
 
   

     !         -------- read the joblist for the pw code --------

     if (esdf_block('pw_jobs', nlines)) then
        njob = 0
        do iline = 1, nlines
           if (index(block_data(iline), 'pw_job') > 0) then
              njob = njob + 1
           else
              write(*, '(A/A)') &
                   'found funny keyword in pw_jobs structure:', &
                   block_data(iline)
           end if
        end do
        if (njob == 0) call mystop( 'pw_jobs list is empty!' )
        pw_params%njobs = njob
        allocate(pw_params%joblist(njob))
        njob = 0
        do iline = 1, nlines
           if (index(block_data(iline), 'pw_job') > 0) then
              njob = njob + 1
              buf = adjustl(block_data(iline) &
                   (index(block_data(iline), 'pw_job') + 6:))
              ifound = 0
              do i = 1, NUM_PW_JOBS
                 if (index(buf, trim(pw_jobstrings(i))) > 0) then
                    ifound = 1
                    pw_params%joblist(njob) = i - 1
                    exit
                 end if
              end do
              if (ifound == 0) then
                 write(*, '(A/A)') 'illegal pw_job:', buf
                 call mystop
              end if
           else
              write(*, '(A/A)') &
                   'found funny keyword in pw_jobs structure:', &
                   block_data(iline)
           end if
        end do
     else                                ! default: do just the scf loop
        pw_params%njobs = 1
        allocate(pw_params%joblist(1))
        pw_params%joblist(1) = 0
     end if

     ! ---------- read the energy window information ------------------

     if (esdf_block('energy_window', nlines)) then
        nline = 0
        do iline = 1, nlines
           if (index(block_data(iline), 'window') > 0) then
              nline = nline + 1
           else
              write(*, '(A/A)') &
                   'found funny keyword in energy_window structure:', &
                   block_data(iline)
           end if
        end do
        if (nline == 0) &
             call mystop( 'found no energy windows in energy_window structure' )
        pw_params%nenergywindow = nline
        allocate(pw_params%energywindow(2 * nline))
        nline = 0
        do iline = 1, nlines
           if (index(block_data(iline), 'window') > 0) then
              nline = nline + 1
              buf = adjustl(block_data(iline) &
                   (index(block_data(iline), 'window') + 6:))
              read(buf, *) pw_params%energywindow((nline - 1) * 2 + 1:nline * 2)
           else
              write(*, '(A/A)') &
                   'found funny keyword in energy_window structure:', &
                   block_data(iline)
           end if
        end do
     else
        pw_params%nenergywindow = 0
        allocate(pw_params%energywindow(1))
     end if

     ! ---------- read the stay_put information ------------------

     if (esdf_block('stay_put', nlines)) then
        nline = 0
        do iline = 1, nlines
           if (index(block_data(iline), 'stay_put') > 0) then
              nline = nline + 1
           else
              write(*, '(A/A)') &
                   'found funny keyword in stay_put structure:', &
                   block_data(iline)
           end if
        end do
        if (nline == 0) &
             call mystop( 'found no constraints in stay_put structure' )
        pw_params%nstay_put = nline
        allocate(pw_params%stay_put(nline))
        nline = 0
        do iline = 1, nlines
           if (index(block_data(iline), 'stay_put_lattice') > 0) then
              nline = nline + 1
              buf = adjustl(block_data(iline) &
                   (index(block_data(iline), 'stay_put_lattice') + 16:))
              if (index(buf, 'xx') > 0) then
                 pw_params%stay_put(nline) = -1
              else if (index(buf, 'yy') > 0) then
                 pw_params%stay_put(nline) = -2
              else if (index(buf, 'zz') > 0) then
                 pw_params%stay_put(nline) = -3
              else if (index(buf, 'xy') > 0 .or. index(buf, 'yx') > 0) then
                 pw_params%stay_put(nline) = -4
              else if (index(buf, 'yz') > 0 .or. index(buf, 'zy') > 0) then
                 pw_params%stay_put(nline) = -5
              else if (index(buf, 'xz') > 0 .or. index(buf, 'zx') > 0) then
                 pw_params%stay_put(nline) = -6
              else
                 write(*, '(A/A)') &
                      'illegal option to stay_put_lattice keyword:', buf
                 call mystop
              end if
           else if (index(block_data(iline), 'stay_put') > 0) then
              nline = nline + 1
              buf = adjustl(block_data(iline) &
                   (index(block_data(iline), 'stay_put') + 8:))
              read(buf, *) pw_params%stay_put(nline)
           else
              write(*, '(A/A)') &
                   'found funny keyword in stay_put structure:', &
                   block_data(iline)
           end if
        end do
     else
        pw_params%nstay_put = 0
        allocate(pw_params%stay_put(1))
     end if

     !         ---- generate the kpoints for a band structure plot --------

     istop = 0  ! are any of the jobs a band structure calculation ?
     do i = 1, pw_params%njobs
        if (pw_params%joblist(i) == 3) istop = 1
     end do

     nitems = 0
     if (istop == 1) then
        nitems = generate_kpts(pw_params,altbands%min(1))
     else
        allocate(pw_params%bslabels(1))
     end if

  end if  ! END of myproc == 0 test


  call MPI_Barrier(mycomm, ierr)
  if (myproc == 0) then
    ibuf(1) = altkpoints%nrk
    ibuf(2:3) = bands%min(1:2) 
    ibuf(4:6) = altkpoints%grid(1:3)
    ibuf(7) = crys%nspin
    ibuf(8) = pw_params%input
    ibuf(9) = pw_params%maxitdiag
    ibuf(10) = pw_params%maxitscf
    ibuf(11) = pw_params%ilog
    ibuf(12) = pw_params%nprocpvm
    ibuf(13) = pw_params%iscr
    ibuf(14) = pw_params%occupy
    ibuf(15) = syms%ntrans
    ibuf(16) = pw_params%tfw_max_cg
    ibuf(17) = pw_params%maxitdiag_band
    ibuf(18) = pw_params%smearing_method
    ibuf(19:22) = pw_params%mfield%gvec(1:4)
    ibuf(23:26) = pw_params%epot%gvec(1:4)
    ibuf(27) = pw_params%nplotwv
    ibuf(28) = pw_params%nlineplot
    ibuf(29) = pw_params%nenergywindow
    ibuf(30) = pw_params%nstay_put
    ibuf(31) = pw_params%miscflag
    ibuf(32:35) = pw_params%output(1:4)
    ibuf(36) = pw_params%optimize
    ibuf(37) = pw_params%checkpoint_wfn
    ibuf(38) = pw_params%njobs
    ibuf(39) = pw_params%nrg
    ibuf(40) = pw_params%nbandsfft
    ibuf(41) = nitems
    ibuf(42:43) = altbands%min(1:2)
    ibuf(44) = pw_params%nband_tddft
    ibuf(45) = pw_params%nmrkpts
    ibuf(46) = pw_params%pp_format
    ibuf(47) = pw_params%numpk_befptf
    ibuf(48) = pw_params%nsliceplot
    ibuf(49) = crys%nalpha
    ibuf(50) = crys%nbeta  !F. Mauri
    ibuf(51) = bands%num_gwout
    ibuf(52) = pw_params%itdiag_add
    ibuf(53) = pw_params%itdiag_add_metal
! LSDA+U. P. Zhang
    ibuf(54) = pw_params%lsdau
    ibuf(55) = pw_params%wannier
! GW band re-ordering. P. Zhang
    ibuf(56) = bands%gw_band_reordering
    ibuf(57) = pw_params%iexact_diag
! Weiwei Gao: for GW energy integration
    ibuf_gw_sum = pw_params%gw_sum
  end if

  call MPI_Bcast(ibuf(1), 57, MPI_INTEGER, 0, mycomm, ierr)

! Weiwei Gao: for GW energy integration gw_sum
  call MPI_Bcast(ibuf_gw_sum, 1, MPI_INTEGER, 0, mycomm, ierr)

  if (myproc /= 0) then

! Weiwei Gao:
    pw_params%gw_sum = ibuf_gw_sum
    altkpoints%nrk = ibuf(1) 
    bands%min(1:2)   = ibuf(2:3)
    altkpoints%grid(1:3) = ibuf(4:6)
    crys%nspin = ibuf(7)
    pw_params%input = ibuf(8)
    pw_params%maxitdiag = ibuf(9)
    pw_params%maxitscf = ibuf(10)
    pw_params%ilog = ibuf(11)
    pw_params%nprocpvm = ibuf(12)
    pw_params%iscr = ibuf(13)
    pw_params%occupy = ibuf(14)
    syms%ntrans = ibuf(15)
    pw_params%tfw_max_cg = ibuf(16) 
    pw_params%maxitdiag_band = ibuf(17)
    pw_params%smearing_method = ibuf(18)
    pw_params%mfield%gvec(1:4) = ibuf(19:22)
    pw_params%epot%gvec(1:4) = ibuf(23:26)
    pw_params%nplotwv = ibuf(27)
    pw_params%nlineplot = ibuf(28)
    pw_params%nenergywindow = ibuf(29)
    pw_params%nstay_put = ibuf(30)
    pw_params%miscflag = ibuf(31)
    pw_params%output(1:4) = ibuf(32:35)
    pw_params%optimize = ibuf(36)
    pw_params%checkpoint_wfn = ibuf(37)
    pw_params%njobs = ibuf(38)
    pw_params%nrg = ibuf(39)
    pw_params%nbandsfft = ibuf(40)
    nitems = ibuf(41)
    altbands%min(1:2) = ibuf(42:43)
    pw_params%nband_tddft = ibuf(44)
    pw_params%nmrkpts = ibuf(45)
    pw_params%pp_format = ibuf(46)
    pw_params%numpk_befptf = ibuf(47) 
    pw_params%nsliceplot =  ibuf(48) 
    crys%nalpha = ibuf(49)
    crys%nbeta = ibuf(50)
    bands%num_gwout = ibuf(51) 
    pw_params%itdiag_add = ibuf(52) 
    pw_params%itdiag_add_metal = ibuf(53) 

! LSDA+U. P. Zhang
    pw_params%lsdau=ibuf(54)
    pw_params%wannier=ibuf(55)
    bands%gw_band_reordering=ibuf(56)
    pw_params%iexact_diag=ibuf(57)



  end if


  call MPI_Barrier(mycomm, ierr)
  if (myproc == 0) then
    cbuf(1:2)=pw_params%icorr(1:2)
    cbuf(3:22)=pw_params%mix_method  
    cbuf(23:42)=pw_params%diag_method  
  end if
  call MPI_Bcast(cbuf, 42, MPI_CHARACTER, 0, mycomm, ierr)
  if (myproc /= 0) then
    pw_params%icorr(1:2)=cbuf(1:2)
    pw_params%mix_method=cbuf(3:22)
    pw_params%diag_method =cbuf(23:42)
  end if


 if (myproc == 0) then
     dbuf(1:3) = pw_params%q(1:3)
     dbuf(4) = pw_params%smearing
     dbuf(5:7) = altkpoints%shift(1:3)
     dbuf(8:10) = pw_params%bra(1:3)
     dbuf(11:13) = pw_params%ket(1:3)
     dbuf(14) = pw_params%epsdiag
     dbuf(15) = pw_params%shiftsafety
     dbuf(16) = pw_params%random_startvec
     dbuf(17) = pw_params%bandgap
     dbuf(18) = pw_params%emax
     dbuf(19) = pw_params%emaxsub
     dbuf(20) = pw_params%emaxmix
     dbuf(21) = pw_params%correct_stress_cutoff
     dbuf(22) = pw_params%chi_q
     dbuf(23:31) = pw_params%g0mask(1:9)
     dbuf(32) = pw_params%mfield%h
     dbuf(33) = pw_params%epot%h
     dbuf(34) = pw_params%epscv
     dbuf(35:38) = pw_params%alphamix(1:4)
     dbuf(39:41) = pw_params%ekinmod(1:3)
     dbuf(42) = energs%inputfermi
     dbuf(43) = pw_params%epsce
     dbuf(44) = pw_params%bd_lgth
     dbuf(45) = bands%gw_mid_energy
     dbuf(46) = bands%gw_low_energy
     dbuf(47) = bands%gw_high_energy
     dbuf(48) = pw_params%init_dir_econv    
  end if
  call MPI_Bcast(dbuf(1), 48, MPI_DOUBLE_PRECISION, 0, mycomm, ierr)
  if (myproc /= 0) then
     pw_params%q(1:3) = dbuf(1:3)
     pw_params%smearing = dbuf(4)
     altkpoints%shift(1:3) = dbuf(5:7)
     pw_params%bra(1:3) = dbuf(8:10)
     pw_params%ket(1:3) = dbuf(11:13)
     pw_params%epsdiag = dbuf(14)
     pw_params%shiftsafety = dbuf(15)
     pw_params%random_startvec = dbuf(16)
     pw_params%bandgap = dbuf(17)
     pw_params%emax = dbuf(18)
     pw_params%emaxsub = dbuf(19)
     pw_params%emaxmix = dbuf(20)
     pw_params%correct_stress_cutoff = dbuf(21)
     pw_params%chi_q = dbuf(22)
     pw_params%g0mask(1:9) = dbuf(23:31)
     pw_params%mfield%h = dbuf(32)
     pw_params%epot%h = dbuf(33)
     pw_params%epscv = dbuf(34)
     pw_params%alphamix(1:4) = dbuf(35:38)
     pw_params%ekinmod(1:3) = dbuf(39:41)
     energs%inputfermi = dbuf(42)
     pw_params%epsce =  dbuf(43) 
     pw_params%bd_lgth = dbuf(44)
     bands%gw_mid_energy = dbuf(45) 
     bands%gw_low_energy = dbuf(46) 
     bands%gw_high_energy = dbuf(47) 
     pw_params%init_dir_econv =  dbuf(48)  


     if (pw_params%nplotwv > 0) then
        allocate(pw_params%plotwv(2 * pw_params%nplotwv))
     else
        allocate(pw_params%plotwv(1))
     end if

     if (pw_params%njobs > 0) then
        allocate(pw_params%joblist(pw_params%njobs))
     else
        allocate(pw_params%joblist(1))
     end if

     if (pw_params%nstay_put > 0) then
        allocate(pw_params%stay_put(pw_params%nstay_put))
     else
        allocate(pw_params%stay_put(1))
     end if

     if (pw_params%nlineplot > 0) then
        allocate(pw_params%lineplot(7 * pw_params%nlineplot))
     else
        allocate(pw_params%lineplot(1))
     end if

     if (pw_params%nsliceplot > 0) then
        allocate(pw_params%sliceplot(11 * pw_params%nsliceplot))
     else
        allocate(pw_params%sliceplot(1))
     end if

     if (pw_params%nenergywindow > 0) then
        allocate(pw_params%energywindow(2 * pw_params%nenergywindow))
     else
        allocate(pw_params%energywindow(1))
     end if

     if (nitems > 0) then
        allocate(pw_params%bslabels(100 * nitems))
     else
        allocate(pw_params%bslabels(1))
     end if

  end if

  call MPI_Bcast(pw_params%NL_rspace(1), 2, MPI_LOGICAL, 0, mycomm, ierr)
  call MPI_Bcast(pw_params%io_scf, 1, MPI_LOGICAL, 0, mycomm, ierr)

 if (pw_params%pp_format .eq. 3 ) call MPI_Bcast(pw_params%pp_iloc(1),& 
          crys%ntype, MPI_INTEGER, 0, mycomm, ierr)
 if (pw_params%pp_format .eq. 3 ) call MPI_Bcast(pw_params%pp_pwfocc,& 
          3*crys%ntype, MPI_DOUBLE_PRECISION, 0, mycomm, ierr)
 if (pw_params%NL_rspace(1) ) call MPI_Bcast(pw_params%NLPP_rcut(1),& 
          crys%ntype, MPI_DOUBLE_PRECISION, 0, mycomm, ierr)

  if (pw_params%nplotwv > 0) call MPI_Bcast(pw_params%plotwv(1),& 
          2 * pw_params%nplotwv, MPI_INTEGER, 0, mycomm, ierr)
  if (pw_params%njobs > 0) call MPI_Bcast(pw_params%joblist(1), &
          pw_params%njobs,  MPI_INTEGER, 0, mycomm, ierr)
  if (pw_params%nstay_put > 0) call MPI_Bcast(pw_params%stay_put(1), &
          pw_params%nstay_put, MPI_INTEGER, 0, mycomm, ierr)
  if (pw_params%nlineplot > 0) call MPI_Bcast(pw_params%lineplot(1), &
          7 * pw_params%nlineplot, MPI_DOUBLE_PRECISION, 0, mycomm, ierr)
   if (pw_params%nsliceplot > 0) call MPI_Bcast(pw_params%sliceplot(1), &
          11 *pw_params%nsliceplot, MPI_DOUBLE_PRECISION, 0, mycomm, ierr)
  if (pw_params%nenergywindow > 0) call MPI_Bcast(pw_params%energywindow, &
          2 * pw_params%nenergywindow, MPI_DOUBLE_PRECISION, 0, mycomm, ierr)
  if (nitems > 0) call MPI_Bcast(pw_params%bslabels(1), 100 * nitems, &
          MPI_CHARACTER, 0, mycomm, ierr)


end subroutine read_pwparam
!
! ============================ support routines ===================
!
subroutine mvfile(name1, name2)

  implicit none
  character(len=*) :: name1, name2

  call mysystem(' if test -f '//trim(name1)//' ; then mv -f '// &
       trim(name1)//' '//trim(name1)//trim(name2)//' ; fi')

end subroutine mvfile
! ---------------------------------------------------------------
function matinvert(mat)

  use constants
  implicit none

  real(dp) :: matinvert
  real(dp), intent(inout) :: mat(9)
  real(dp) :: w(9), vol
  vol=(mat(1)*(mat(5)*mat(9)-mat(8)*mat(6)) &
	 -mat(2)*(mat(4)*mat(9)-mat(7)*mat(6)) &
	 +mat(3)*(mat(4)*mat(8)-mat(7)*mat(5)))
  w = mat

  mat(1) = (w(5)*w(9)-w(6)*w(8))/vol
  mat(4) = (w(6)*w(7)-w(4)*w(9))/vol
  mat(7) = (w(4)*w(8)-w(5)*w(7))/vol
  
  mat(2) = (w(8)*w(3)-w(9)*w(2))/vol
  mat(5) = (w(9)*w(1)-w(7)*w(3))/vol
  mat(8) = (w(7)*w(2)-w(8)*w(1))/vol

  mat(3) = (w(2)*w(6)-w(3)*w(5))/vol
  mat(6) = (w(3)*w(4)-w(1)*w(6))/vol
  mat(9) = (w(1)*w(5)-w(2)*w(4))/vol

  matinvert = vol

end function matinvert
! ---------------------------------------------------------------
function getmass(alen, sym)

  use constants
  implicit none
  real(dp) :: getmass
  integer, intent(in) :: alen
  character(len=alen), intent(in) :: sym

  ! a table with all the masses of the atoms in u

  type mass_table
     character(len=2) :: name
     real(dp) :: mass
  end type mass_table

  integer, parameter :: AT_TB_SIZE = 107

  type(mass_table), parameter :: table_of_masses(AT_TB_SIZE) = (/ &
       mass_table('Ac',  227.0d0 ), mass_table('Al',  26.981539d0 ), &
       mass_table('Am',  243.0d0 ), mass_table('Sb',  121.760d0 ), &
       mass_table('Ar',  39.948d0 ), mass_table('As',  74.92159d0 ), &
       mass_table('At',  210d0 ), mass_table('Ba',  137.327d0 ), &
       mass_table('Bk',  247.0d0 ), mass_table('Be',  9.012182d0 ), &
       mass_table('Bi',  208.98037d0 ), mass_table(' B', 10.811d0 ), &
       mass_table('Br',  79.904d0 ), mass_table('Cd',  112.411d0 ), &
       mass_table('Cs',  132.90543d0 ), mass_table('Ca',  40.078d0 ), &
       mass_table('Cf',  251.0d0 ), mass_table(' C', 12.011d0 ), &
       mass_table('Ce',  140.115d0 ), mass_table('Cl',  35.4527d0 ), &
       mass_table('Cr',  51.9961d0 ), mass_table('Co',  58.93320d0 ), &
       mass_table('Cu',  63.546d0 ), mass_table('Cm',  247.0d0 ), &
       mass_table('Dy',  162.50d0 ), mass_table('Es',  252.0d0 ), &
       mass_table('Er',  167.26d0 ), mass_table('Eu',  151.965d0 ), &
       mass_table('Fm',  257d0 ), mass_table(' F', 18.9984032d0 ), &
       mass_table('Fr',  223.0d0 ), mass_table('Gd',  157.25d0 ), &
       mass_table('Ga',  69.723d0 ), mass_table('Ge',  72.61d0 ), &
       mass_table('Au',  196.96654d0 ), mass_table('Hf',  178.49d0 ), &
       mass_table('Jl',  262.0d0 ), mass_table('He',  4.002602d0 ), &
       mass_table('Ho',  164.93032d0 ), mass_table(' H', 1.00794d0 ), &
       mass_table('H0',  0.50032d0 ), mass_table('H1', 1.50794d0 ), &
       mass_table('In',  114.818d0 ), mass_table(' I', 126.90447d0 ), &
       mass_table('Ir',  192.217d0 ), mass_table('Fe',  55.845d0 ), &
       mass_table('Kr',  83.80d0 ), mass_table('La',  138.9055d0 ), &
       mass_table('Lr',  262.0d0 ), mass_table('Pb',  207.2d0 ), &
       mass_table('Li',  6.941d0 ), mass_table('Lu',  174.967d0 ), &
       mass_table('Mg',  24.3050d0 ), mass_table('Mn',  54.93805d0 ), &
       mass_table('Md',  258.0d0 ), mass_table('Hg',  200.59d0 ), &
       mass_table('Mo',  95.94d0 ), mass_table('Nd',  144.24d0 ), &
       mass_table('Ne',  20.1797d0 ), mass_table('Np',  237.0d0 ), &
       mass_table('Ni',  58.6934d0 ), mass_table('Nb',  92.90638d0 ), &
       mass_table(' N', 14.00674d0 ), mass_table('No',  259.0d0 ), &
       mass_table('Os',  190.23d0 ), mass_table(' O', 15.9994d0 ), &
       mass_table('Pd',  106.42d0 ), mass_table(' P', 30.973762d0 ), &
       mass_table('Pt',  195.08d0 ), mass_table('Pu',  244.0d0 ), &
       mass_table('Po',  209.0d0 ), mass_table(' K', 39.0983d0 ), &
       mass_table('Pr',  140.90765d0 ), mass_table('Pm',  145.0d0 ), &
       mass_table('Pa',  231.03588d0 ), mass_table('Ra',  226.0d0 ), &
       mass_table('Rn',  222.0d0 ), mass_table('Re',  186.207d0 ), &
       mass_table('Rh',  102.90550d0 ), mass_table('Rb',  85.4678d0 ), &
       mass_table('Ru',  101.07d0 ), mass_table('Db',  261.0d0 ), &
       mass_table('Sm',  150.36d0 ), mass_table('Sc',  44.955910d0 ), &
       mass_table('Se',  78.96d0 ), mass_table('Si',  28.0855d0 ), &
       mass_table('Ag',  107.8682d0 ), mass_table('Na',  22.989768d0 ), &
       mass_table('Sr',  87.62d0 ), mass_table(' S', 32.066d0 ), &
       mass_table('Ta',  180.9479d0 ), mass_table('Tc',  98.0d0 ), &
       mass_table('Te',  127.60d0 ), mass_table('Tb',  158.92534d0 ), &
       mass_table('Tl',  204.3833d0 ), mass_table('Th',  232.0381d0 ), &
       mass_table('Tm',  168.93421d0 ), mass_table('Sn',  118.710d0 ), &
       mass_table('Ti',  47.867d0 ), mass_table(' W', 183.84d0 ), &
       mass_table(' U', 238.0289d0 ), mass_table(' V', 50.9415d0 ), &
       mass_table('Xe',  131.29d0 ), mass_table('Yb',  173.04d0 ), &
       mass_table(' Y', 88.90585d0 ), mass_table('Zn',  65.39d0 ), &
       mass_table('Zr',  91.224d0 ) /)

  integer :: i

  do i = 1, AT_TB_SIZE
     if (adjustr(sym(1:2)) == table_of_masses(i)%name) then
        getmass = table_of_masses(i)%mass
        return
     end if
  end do

  write(*,'(3A)') 'Could not find the element ', trim(sym), &
       ' in the table of masses!'
  call mystop

end function getmass
