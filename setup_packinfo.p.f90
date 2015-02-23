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

!-*-Fortran-*-
!     
!     @process extchk
!
subroutine setup_packinfo(gs)







  use all_to_all_module
  include 'use.h'
  implicit none             ! implicit? Just say no!
include 'mpif.h'
  include 'interface.h'
  include 'all_to_all.h'
  !     1996 by Bernd Pfrommer while at UC Berkeley
  !
  !     INPUT/OUTPUT:
  !     ------------
  !
  type(parallel_gspace), intent(inout) :: &
       gs                   ! the gspace for which the data is set up
  !     
  !     DESCRIPTION:
  !     -----------
  !
  !
  !     Sets up the necessary arrays to perform the 3d fft in parallel.
  !     This involves attaching the following memory blocks to the
  !     gspace structure:
  !
  !                                      ,---- pack/unpack
  !     gs%packinfo(maxpacksize_t1,nproc,2)
  !                     \              
  !                      `-- which element of the packet
  !
  !                        ,---- pack/unpack
  !     gs%packsize(nproc,2,2)                 = number of dc variables sent
  !                          \              
  !                           `-- which phase of the transpose (first or second)
  !
  !                     ,---- inverse=1, forward =2 fft
  !     gs%pickup(nproc,2)  = where on the other proc the packed data
  !                           for the t2 phase resides (only used for shmem)
  !          WARNING: altered for simultaneous FFT of all wavefunctions
  !    
  !
  !     packing index:
  !
  !     the index=1 refers to "packing    for proc p" for inverse FFT
  !                       and "unpacking from proc p"     forward
  !     the index=2 refers to "unpacking from proc p" for inverse FFT
  !                       and "packing    for proc p"     forward
  !
  !     transposition index:
  !
  !     the index=1 refers to  1st transpose of inverse FFT
  !                        and 2nd              forward
  !               2 refers to  2nd transpose of inverse FFT
  !                        and 1st              forward
  !
  !     Also, gs%maxpacksize_t1/t2 is set to the maximum package size 
  !     for the first/second transpose
  !
  !
  !     packinfo and packsize contain the necessary information to 
  !     pack and unpack the messages exchanged between the processors
  !     during the transposition phase of the inverse and forward FFT.
  !
  !     For the second transpose of the FFT grid, the following array
  !     is used:
  !                    1    2                 
  !     gs%chunk(1..4,pack/unpack,whichchunk,proc)
  !
  !     Each chunk is described by 4 numbers in the first index:
  !
  !     1 = starting rod, as indexed on the local processor
  !     2 = starting y offset for packing, x offset for unpacking
  !     3 = length of chunk along y 
  !     4 =  width of chunk along x 
  !
  !     For the last chunk, the width is negative. This serves as a 
  !     flag to stop packing/unpacking.
  !
  !
  !
  !
  !     WARNING: do not call and ask to set up packinfo with gs%nproc=1
  !     unless there is really only one processor up and running.
  !
  !     This will crash, since the gather routine
  !     will try to collect data into the arrays, which are dimensioned
  !     too small then.
  !
  !
  !     ------------------ local variables ------------------------------
  !
  integer :: &
       nnx, nny, nnz, maxpacksize, maxpacksize_t2, ptarget, &
       ipack, iord, irod, totpacksize_t2, &
       rodrange(2), nx, ny, nz, nxp, nrmxz, nrmyz, &
       nxs, nys, irods, k, ipack1, ldx, ldy, ldz, &
       j, i, ierror, nrodst1, nrodst2, ilen, nchnk, &
       nrodst2t, b2t, na2t, dnx,b,b2,na,na2, ip
  logical :: inewchunk
  integer, allocatable :: &
       mypack(:), helper(:), destinfo(:,:), myunpack(:)
  integer, allocatable :: aofb(:), packdum(:), &
       recvbuf(:,:,:)

  real(dp) :: t0
  real(dp), external :: gimmetime 

  t0 = gimmetime()  

  if (.not. gs%ipackinfo) then
     gs%maxpacksize_t1 = -1
     gs%maxpacksize_t2 = -1
     nullify(gs%packinfo)
     nullify(gs%packsize)
     nullify(gs%pickup)
     return
  end if

  nnx = gs%fftsize(1) ; ldx = nnx
  nny = gs%fftsize(2) ; ldy = nny
  nnz = gs%fftsize(3) ; ldz = nnz




  !
  !     x-dimension size of the slab obtained after the first FFT
  !     and the first transpose
  !
  dnx = mod(gs%fftsize(4) - gs%fftsize(5) + nnx, nnx) + 1
  !
  !     minimum number of rows each processor holds after the first FFT
  !     and the first transpose
  !
  b = (dnx * nnz) / gs%nproc 
  !
  !     number of processors with an additional rod. The lower-number
  !     processors hold the additional rods.
  !
  na = mod(dnx * nnz, gs%nproc)
  !
  !     number of rods on this processor after first fft
  !
  nrodst1 = b
  if (na> gs%myproc) nrodst1 = nrodst1 + 1
  !
  !     minimum number of rows each processor holds after the 2nd
  !     transpose
  !
  b2 = (nny * nnz) / gs%nproc 
  !
  !     number of processors with an additional rod after the 2nd FFT. 
  !     The lower-number processors hold the additional rods.
  !
  na2 = mod(nny * nnz, gs%nproc)
  !
  !     number of rods this processor holds
  !
  nrodst2 = b2
  if (na2 > gs%myproc) nrodst2 = nrodst2 + 1
  !
  !     figure out the maximum packet size,
  !     including the packet to itself
  !
  nrmxz = (dnx * nnz) / gs%nproc + 1
  nrmyz = (nny * nnz) / gs%nproc + 1

  maxpacksize = 0
  do ptarget = 0, gs%nproc - 1   ! where the packet goes
     rodrange(1) = b * ptarget + min(ptarget, na)
     rodrange(2) = rodrange(1) + b
     if (na <= ptarget) rodrange(2) = rodrange(2) - 1
     ipack = 0
     do iord = 1, gs%lorder
        irod = gs%order(1, iord)
        nx = irod / nny
        ny = mod(irod, nny)
        nxp = mod(nx - gs%fftsize(5) + nnx, nnx)
        do nz = 0, nnz - 1
           i = nz * dnx + nxp     ! rod index
           if (i >= rodrange(1) .and. i <= rodrange(2)) then
              ipack = ipack + 1
           end if
        end do
     end do

     if (ipack > maxpacksize) maxpacksize = ipack
  end do

  call all_max_all(maxpacksize)

  gs%maxpacksize_t1 = maxpacksize

  allocate(gs%packinfo(maxpacksize, gs%nproc, 2))
  gs%packinfo = 0
  allocate(gs%packsize(gs%nproc, 2, 2))
  gs%packsize = 0
  allocate(gs%pickup(gs%nproc, 2))
  gs%pickup = 0

  allocate(mypack(maxpacksize))
  allocate(helper(max(nny * nrmxz, nnx * nrmyz)))



  !
  !     ---- set up for the 1st transposition phase of FFT ---------
  !


  !      write(9,*) '------ first transposition packing info----- '
  !      call myflush(9)

  !
  !     
  !
  do ptarget = 0, gs%nproc - 1   ! where the packet goes
     !         write(9,*) ' === target proc:',ptarget
     !
     !        those are the rod indices held by the target processor
     !
     rodrange(1) = b * ptarget + min(ptarget, na)
     rodrange(2) = rodrange(1) + b
     if (na <= ptarget) rodrange(2) = rodrange(2) - 1


     !
     !        for better cache behavior, loop in stride one through the
     !        target array.
     !
     !
     !        first set up a helper index array to cut down operation count
     !
     helper(:) = -1
     do iord = 1, gs%lorder
        irod = gs%order(1, iord)
        nx = irod / nny
        ny = mod(irod, nny)
        nxp = mod(nx - gs%fftsize(5) + nnx, nnx)
        do nz = 0, nnz - 1
           i = nz * dnx + nxp     ! rod index
           if (i >= rodrange(1) .and. i <= rodrange(2)) then
              helper((i - rodrange(1)) * nny + ny + 1) = (iord - 1) * nnz + &
                   nz + 1
           end if
        end do
     end do
     !
     !        now pick up the data points
     !
     ipack = 0
     do irod=rodrange(1), rodrange(2)
        nz = irod / dnx
        nx = mod(irod, dnx)
        do ny = 0, nny - 1
           if (helper((irod - rodrange(1)) * nny + ny + 1) >= 0) then
              ipack = ipack + 1
              !
              !   set up unpacking(inverse) and packing(forward) info. 
              !
              mypack(ipack) = (irod - rodrange(1)) * nny + ny + 1
              k = ipack

              gs%packinfo(ipack, ptarget + 1, 1) = &
                   helper((irod - rodrange(1)) * nny + ny + 1)

              if (k > gs%maxpacksize_t1) then
                 write(9, *)'bombout in setup_packinfo, 1st transp:', &
                      k, gs%maxpacksize_t1
                 call mystop
              end if
           end if
        end do
     end do

     !
     !        communicate the unpacking info and packet size 
     !        to the target processor. This is done by specifying
     !        ptarget as root of the gather call.
     !

     call mpi_gather(mypack, maxpacksize, MPI_INTEGER, &
          gs%packinfo(1, 1, 2), maxpacksize, MPI_INTEGER, &
          ptarget, MPI_COMM_WORLD, ierror) 
     call mpi_gather(ipack, 1, MPI_INTEGER, &
          gs%packsize(1, 2, 1), 1, MPI_INTEGER, ptarget, &
          MPI_COMM_WORLD, ierror) 

     gs%packsize(ptarget + 1, 1, 1) = ipack


  end do


  nullify(gs%packtype)



  !
  !     ---- set up for the 2nd transposition phase of  FFT ---------
  !

  !     
  !     startrod (not index)
  !     start x/y
  !     len
  !     width       1       2
  !     chunk(1..4,pack/unpack,whichchunk,proc)
  !     

  !      write(9,*) '------ second transposition packing info----- '
  !      call myflush(9)

  gs%maxnchunk = gs%nproc * (gs%fftsize(3) / gs%nproc + 3)
  allocate(gs%chunk(4, 2, gs%maxnchunk, gs%nproc))
  allocate(destinfo(4, gs%maxnchunk))

  allocate(recvbuf(4, gs%maxnchunk, gs%nproc))


!  write(9,*)' mid', gimmetime() - t0  
!  call myflush(9)  

  destinfo = 0
  gs%chunk = 0
  maxpacksize_t2 = 0
  do ptarget = 0, gs%nproc - 1   ! where the packet goes
     rodrange(1) = b2 * ptarget + min(ptarget, na2)
     rodrange(2) = rodrange(1) + b2
     if (na2 <= ptarget) rodrange(2) = rodrange(2) - 1

     ipack = 0              ! initial packet length is 0
     nchnk = 0              ! assume no chunks from this proc to ptarget

     inewchunk = .true.
     do irod = 0, nrodst1 - 1    ! loop through all rods on this proc
        i = b * gs%myproc + min(gs%myproc, na) + irod ! source rod index
        nz = i / dnx
        nxp = mod(i, dnx)
        nx = mod(nxp + gs%fftsize(5), nnx)
        !     
        !           see if a new chunk starts here
        !     
        if (nx == gs%fftsize(5)) inewchunk = .true.

        !     
        !  if yes, then find out about its length, and whether it goes to
        !  ptarget at all.
        !            
        if (inewchunk) then  ! this could be a new chunk. lets see...
           ilen = 0
           do ny = 0, nny - 1  ! run along y to see if it really is a new chunk
              j = ny + nz * nny ! target rod index
              if (j >= rodrange(1) .and. j <= rodrange(2)) then
                 !     
                 !     found start of new chunk
                 !
                 if (ilen == 0) then
                    if (nchnk + 1 > gs%maxnchunk) then
                       write(9, *) 'nchnk too big in setup_packinfo'
                       call mystop
                    end if
                    gs%chunk(2, 1, nchnk + 1, ptarget + 1) = ny
                    destinfo(1, nchnk + 1) = j - rodrange(1) ! target rod index
                 end if
                 ilen = ilen + 1
              end if
           end do
           if (ilen > 0) then ! there really is a new chunk for this processor
              !                  write(9,*) ' found new chunk!'
              nchnk = nchnk + 1
              gs%chunk(3, 1, nchnk, ptarget + 1) = ilen ! length of chunk
              destinfo(3, nchnk) = ilen ! length of chunk is same for dest.
              !     
              !     determine width of chunk
              !     
              gs%chunk(4, 1, nchnk, ptarget + 1) = &
                   min(dnx - nxp, nrodst1 - irod)
              destinfo(4, nchnk) = gs%chunk(4, 1, nchnk, ptarget + 1)
              !     
              ! this is the rod on the source processor where the chunk starts
              !
              gs%chunk(1, 1, nchnk, ptarget + 1) = irod
              !     
              !     this is the x offset for the destination processor
              !
              destinfo(2, nchnk) = nx
              !
              !                 increase the total packet size
              !
              ipack = ipack + ilen * gs%chunk(4, 1, nchnk, ptarget + 1)
           end if
           inewchunk = .false. ! go on with search for chunk start
        end if
     end do

     if (nchnk > 0) then
        !     
        !           mark last chunk by setting the width field negative!
        !
        gs%chunk(4, 1, nchnk, ptarget + 1) = &
             -gs%chunk(4, 1, nchnk, ptarget + 1) 
        destinfo(4, nchnk) = -destinfo(4, nchnk)
     end if
     !
     !        communicate the unpacking info in destinfo to the target
     !        processor. 
     !        This is done by specifying   ptarget as root of the gather call.
     !
     gs%packsize(ptarget + 1, 1, 2) = ipack 

     call mpi_gather(destinfo(1, 1), gs%maxnchunk * 4, MPI_INTEGER, &
          recvbuf(1, 1, 1), gs%maxnchunk * 4, MPI_INTEGER, &
          ptarget, MPI_COMM_WORLD, ierror) 
!     gs%chunk(:, 2, :,:) = recvbuf(:,:,:) ! massage into shape...     

     call mpi_gather(ipack, 1, MPI_INTEGER, &
          gs%packsize(1, 2, 2), 1, MPI_INTEGER, ptarget, &
          MPI_COMM_WORLD, ierror)

     if (gs%myproc .eq. ptarget) then
     do j=1,gs%nproc
       do i=1,gs%maxnchunk
         gs%chunk(1, 2, i,j) = recvbuf(1,i,j) ! massage into shape...
         gs%chunk(2, 2, i,j) = recvbuf(2,i,j)
         gs%chunk(3, 2, i,j) = recvbuf(3,i,j)
         gs%chunk(4, 2, i,j) = recvbuf(4,i,j) 
       end do
     end do
     end if

 

     if (ipack > maxpacksize_t2) maxpacksize_t2 = ipack
  end do                    ! end of loop over target processor

!  write(9,*)' pre max 2', gimmetime() - t0, gs%maxnchunk * 4  
!  call myflush(9)

  call all_max_all(maxpacksize_t2)
  gs%maxpacksize_t2 = maxpacksize_t2

  deallocate(recvbuf)

  deallocate(destinfo)
  deallocate(helper)
  deallocate(mypack)
  gs%totpacksize_t2 = gs%maxpacksize_t2


  !
  !     Now set up the pickup info for shmem. This small array tells
  !     shmem_get where to pick up the packed data from the remote
  !     node during the t2 stage.
  !
  !     First, we fill the helper array with the location of the
  !     packets to the target processors. If there is no packet, the  
  !     entry will be zero.
  !
  allocate(helper(gs%nproc))
  helper = 0
  ipack = 1 
  do ip = 1, gs%nproc - 1
     ptarget = mod(ip + gs%myproc, gs%nproc) ! the target processor 
     if (gs%packsize(ptarget + 1, 1, 2) > 0) then
        helper(ptarget + 1) = ipack
        ipack = ipack + gs%packsize(ptarget + 1, 1, 2)  
     end if
  end do
  totpacksize_t2 = ipack
  !
  !     now communicate the pickup information to all the other procs
  !
  do ip = 1, gs%nproc
     call mpi_gather(helper(ip), 1, MPI_INTEGER, &
          gs%pickup(1, 1), 1, MPI_INTEGER, ip - 1, &
          MPI_COMM_WORLD, ierror) 
  end do
  !
  !     Now the same thing for the reverse direction, i.e. the t2 stage
  !     during the forward FFT. This is also a good occasion to 
  !
  helper = 0
  ipack = 1 
  do ip = 1, gs%nproc - 1
     ptarget = mod(ip + gs%myproc, gs%nproc) ! the target processor
     if (gs%packsize(ptarget + 1, 2, 2) > 0) then
        helper(ptarget + 1) = ipack
        ipack = ipack + gs%packsize(ptarget + 1, 2, 2) 
     end if
  end do

  gs%totpacksize_t2 = max(totpacksize_t2, ipack)
  !
  !     now communicate the pickup information to all the other procs
  !
  do ip = 1, gs%nproc
     call mpi_gather(helper(ip), 1, MPI_INTEGER, &
          gs%pickup(1, 2), 1, MPI_INTEGER, ip - 1, &
          MPI_COMM_WORLD, ierror) 
  end do

  deallocate(helper)

  !
end subroutine setup_packinfo
