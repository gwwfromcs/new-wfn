!     @process extchk opt
!
subroutine readgsdat(ipr, ierr, gs, data, ldnvecs, nvecs, nspin, &
     filename, namelen)
  !
  use all_to_all_module  
  include 'use.h'  
  implicit none              ! implicit? Just say no!
  include 'interface.h'  
  include 'all_to_all.h'  
  !
  !     INPUT:
  !     -----
  !
  integer, intent(in) :: &
       namelen, &  ! length of filename
       ldnvecs, &  ! leading dimension of nvecs
       ipr         ! print flag
  type(parallel_gspace), intent(in) :: &
       gs          ! the gspace on which the data is set up
  character(len=namelen), intent(in) :: &
       filename    ! name of file to be written
  !
  !     INPUT/OUTPUT:
  !     ------------
  !
  integer, intent(inout) :: &
       nvecs, &    ! number of vectors to be read/ found
       nspin       ! number of spins to be read/ found
  !
  !     OUTPUT:
  !     -------
  !
  complex(dp), intent(out) :: &
       data(gs%length, ldnvecs, nspin)    ! data that should be read from
  integer, intent(out) :: &
       ierr                   ! error code in case something goes wrong
  !
  !     reads a g-space data array from disk
  !
  !     1995 Bernd Pfrommer while at UCB
  !     1997 Bernd Pfrommer improved for better parallel I/O
  !
  !          Now, the data is read only by processor 0, and distributed
  !          across all processors via broadcast
  !
  !     2006 David Prendergast - welcome to the 21st Century
  !          Increased the block size for reading to the maximum and
  !          allowed for decrease if too large
  !          Removed major bug where root process could jump to a point
  !          beyond mpi_bcast statements which would paralyze slave processes
  !
  !     ---------------- local variables ---------------------------------
  !
  integer :: i, j, k, ig1, ig2, ig3, ig3start, ig3end, &
       iord, iproc, i1, i2, i3, nfound, idat, is, n, nspind, nvecsd, ib
  ! integer, parameter :: nb = 64  ! block size for parallel io
  integer :: nb ! block size for parallel io
  real(dp) :: t0
  real(dp), external :: gimmetime
  integer, external :: findvec  
  complex(dp) :: rdata  
  !integer :: igb(3, nb)  
  integer,allocatable :: igb(:,:)  
  complex(dp), allocatable :: dummy(:,:,:)  
  character(len=24) :: bdate  
  character(len=8) :: btime  
  !     ------------------------------------------------------------------
  !
  t0 = gimmetime()  
  ierr = 0  
  data = zzero

  nfound = 0  
  if (gs%myproc == 0) then  
     write(9,*) 'inside readgsdat'
     call myflush(9)
     open(unit = 21, file = filename, status = 'unknown', &
          form = 'unformatted')
     rewind(21)  
     read(21, end = 912) bdate, btime  
     read(21, end = 912) nvecsd, nspind  
     goto 913  
912  ierr = -1  
     call myflush(21)
     close(21)  
913  continue  
  end if
  call my_broadcast(ierr, 0)  
  call my_broadcast(bdate, 24, 0)  
  call my_broadcast(btime, 8, 0)  
  call my_broadcast(nvecsd, 0)  
  call my_broadcast(nspind, 0)  

  if (ierr == -1) goto 911  
  if (ipr >= 1) write(9, 100) filename, bdate  
  if (nvecsd > nvecs .or. nspind > nspin) then  
     write(9, 210) nvecsd, nvecs, nspind, nspin  
     call mystop  
  end if

  if( gs%myproc==0 ) then
    write(9,*) 'about to read'
    call myflush(9)
  endif

  nb=gs%length

  ! determine size of dummy array that fits in memory on all processes

121 continue

  allocate(dummy(nvecsd, nspind, nb),igb(3,nb),stat=ierr)  
  if( ierr/=0 ) then
    if( nb < 2 ) &
      call mystop('readgsdat: problem allocating space dummy')
    nb=nb/2
    if( allocated(dummy) ) deallocate(dummy)
    if( allocated(igb)   ) deallocate(igb)
    goto 121
  endif

  if( gs%myproc==0 ) then
    write(9,*) ' readgsdat:       root block size = ', nb
    call myflush(9)
  endif

  
  ! find the minimum size over processes
  call all_min_all( nb )

  if( allocated(dummy) ) deallocate(dummy)
  if( allocated(igb)   ) deallocate(igb)
  allocate(dummy(nvecsd, nspind, nb),igb(3,nb),stat=ierr)
  if( ierr/=0 ) call mystop('readgsdat: problem allocating space dummy')

  if( gs%myproc==0 ) then
    write(9,*) ' readgsdat: global min block size = ', nb
    call myflush(9)
  endif

  idat = 0  

  do while (.true.)  
     if (gs%myproc == 0) then  
        do ib = 1, nb  
           ! this may be a problem line
           ! read(21, end = 123) igb(1, ib), igb(2, ib), igb(3, ib)  
           read(21, end = 122) igb(1, ib), igb(2, ib), igb(3, ib)  
           if (igb(1, ib) /= -1234567) then  
              read(21,end=922) ((dummy(n, is, ib), n = 1, nvecsd), is = 1, nspind)
           else  
              exit  
           end if
        end do
     end if
     ! new jump point
     
 122 continue
     if( gs%myproc==0 ) then
     ib = min(ib, nb)  
     end if

     ! if( gs%myproc==0 ) then
     !   write(9,*) 'do broadcasting'
     !   call myflush(9)
     ! endif
     call my_broadcast(ib, 0)  
     call my_broadcast(igb(1, 1), 3 * ib, 0)  
     call my_broadcast(dummy(1, 1, 1), ib * nvecsd * nspind, 0)  
     do k = 1, ib  
        if (igb(1, k) == -1234567) goto 123          ! done
        idat = idat + 1  
        i = findvec (igb (1, k), gs)          ! find the index of the gvector
        if (i > 0) then                       ! gvector is here, get data
           data(i, 1:nvecsd, 1:nspind) = dummy(1:nvecsd, 1:nspind, k)  
           nfound = nfound + 1  
        end if
     end do
  end do
123 continue    ! continue here if read reaches end of file

  if (ipr >= 1) write(9, 200) nfound, idat, gimmetime() - t0  

  deallocate(dummy)  
  if (gs%myproc == 0) then
    call myflush(21)
    close(21)  
  end if
  nvecs = nvecsd  
  nspin = nspind  

  ! write(9,*) 'about to leave readgsdat'
  ! call myflush(9)
  return

911 continue    ! in case it couldn't find the header
  ierr = -1  
  if (ipr >= 1) write(9, 220) filename  

  nvecs = 0
  nspin = 0  

  return

922 continue   ! in case incomplete charge data
  ierr = -1
  if( ipr >= 1) write(9,*) 'problem reading block ', ib, ' from ', trim(filename)
  call mystop

  return

100 format(/' READING DATA FROM FILE: ',A, &
       &     /' WRITTEN ON THE ',a24)
200 format(' READ ',i7,' OF ',i7,' COMPONENTS FOR THIS PROCESSOR ', &
       &     'IN ',f12.3,' SECONDS')
210 format(' *** READING ERROR: ', i4,' >',i4,' or ', i4, ' > ',i4)  
220 format(' *** COULD NOT READ VALID DATA FROM FILE: ',a)  

end subroutine readgsdat
