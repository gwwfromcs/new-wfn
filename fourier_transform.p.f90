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
subroutine fourier_transform(idir, ffts, gs, gsdat, rspace, nfn)
  !
  !





  !
  include 'use.h' 
  implicit none             ! implicit? Just say no!
include 'mpif.h'
  include 'interface.h'

  !
  !     INPUT:
  !     -----
  !
  integer, intent(in) :: idir                      ! if >=0  rspace ==> gspace
                                                   ! if <0   rspace <== gspace
  integer, intent(in) :: nfn                      ! number of functions to FFT
  type(fft_struc), intent(in) :: ffts                ! work arrays for the FFT
  type(parallel_gspace), intent(in) :: &
       gs                            ! the gspace for which the data is set up
  !
  !     INPUT/OUTPUT:
  !     ------------
  !
  !     WARNING: UPON TRANSFORM REALSPACE => GSPACE, THE REALSPACE
  !     ARRAY IS DESTROYED.
  !
  complex(dp), intent(inout) :: &
       rspace(gs%r_size * nfn)     ! the data in realspace rspace(ffts%r_size)
  complex(dp), intent(inout) :: &
       gsdat(gs%length * nfn)            ! the data in gspace gsdat(gs%length)
  !
  !     --------------------------------------------------------------------
  !     summer 1996 by Bernd Pfrommer while at UC Berkeley
  !
  !     
  !     DESCRIPTION:
  !     -----------
  !
  !     Takes an array gsdat in fourier space, and transforms it to
  !     realspace, or the other way round. This is done in parallel
  !     using the  message passing  library.
  !
  !
  !     Returns in rspace a bunch of rods along the x direction which
  !     belong to this processor. The data layout is as follows:
  !
  !
  !           ------------ Nz
  !           |          |
  !           |          |        
  !           |          |
  !           |   2222222|
  !           |2211111111|
  !           |1000000000|
  !           ------------
  !           Ny
  !
  !
  !     All the information about which processor receives what from
  !     where is contained in the packinfo array, which is attached to
  !     the gspace structure.
  !
  !
  
      !
      ! edit davegp
      ! I changed this parameter setting so that you can
      ! link to the correct header associated with the library
      ! that you link. Requires an include directory to be
      ! added to conf.mk
      !
      include 'fftw_f77.i'
      !
      ! INTEGER FFTW_FORWARD,FFTW_BACKWARD
      ! PARAMETER (FFTW_FORWARD=-1,FFTW_BACKWARD=1)
      ! INTEGER FFTW_ESTIMATE,FFTW_MEASURE
      ! PARAMETER (FFTW_ESTIMATE=0,FFTW_MEASURE=1)
      ! INTEGER FFTW_IN_PLACE,FFTW_USE_WISDOM
      ! PARAMETER (FFTW_IN_PLACE=8,FFTW_USE_WISDOM=16)

! for 64-bit machines (integer is not enough)
! it should also work on 32-bit machines

      ! the plan tags from create_plan
      INTEGER*8 FFTW_PLAN_FWD, FFTW_PLAN_BWD

  !
  real(dp) :: t0, ttot, normfac
  integer :: ix, iy,iz, icx, icy, icz, i
  real(dp), external :: gimmetime





  if (idir >= 0) then        ! rspace ->gspace
     !
     !        do the first stage of the 3d FFT
     !
     
     call fourier_transform_3(idir, ffts, gs, rspace(1), nfn)
     !
     !        do the last two stages of the 3d FFT
     !
     call fourier_transform_12(idir, ffts, gs, gsdat(1), rspace(1), nfn)


  else                      ! gspace -> rspace
     !
     !        do the last two stages of the 3d FFT
     !
     
     call fourier_transform_12(idir, ffts, gs, gsdat(1), rspace(1), nfn)
     !
     !        do the first stage of the 3d FFT
     !
     call fourier_transform_3(idir, ffts, gs, rspace(1), nfn)


  end if

  return

123 format('tFFTI:',3f8.4,', tUFOLD:',f8.4,', tTRODS:',f8.4, &
       ', tZBUF:', f8.4,', tZRSPC:', f8.4,', tTLT1:', f8.4, &
       ', tTLT2:', f8.4,', tTC1:', f8.4,', tTC2:', f8.4)
124 format('%FFTI:',3f8.1,', %UFOLD:',f8.1,', %TRODS:',f8.1, &
       ', %ZBUF:', f8.1,', %ZRS%C:', f8.1,', %TLT1:', f8.1, &
       ', %TLT2:', f8.1,', %TC1:', f8.1,', %TC2:', f8.1)
223 format('tFFTF:',3f8.4,', tUFOLD:',f8.4,', tTRODS:',f8.4, &
       ', tZBUF:', f8.4,', tZRSPC:', f8.4,', tTLT1:', f8.4, &
       ', tTLT2:', f8.4,', tTC1:', f8.4,', tTC2:', f8.4)
224 format('%FFTF:',3f8.1,', %UFOLD:',f8.1,', %TRODS:',f8.1, &
       ', %ZBUF:', f8.1,', %ZRSPC:', f8.1,', %TLT1:', f8.1, &
       ', %TLT2:', f8.1,', %TC1:', f8.1,', %TC2:', f8.1)

end subroutine fourier_transform
!
!     ================================================================
!
subroutine fourier_transform_12(idir, ffts, gs, gsdat, rspace, nfn)
  !
  use all_to_all_module
  !
  ! edit davegp
  use fftw_hitlist, only: fftw_hitlist_plan
  !
  include 'use.h'
  implicit none             ! implicit? Just say no!
include 'mpif.h'
  include 'interface.h'
  include 'all_to_all.h'
  include 'flibcalls.ph'

  !
  !     INPUT:
  !     -----
  !
  integer, intent(in) :: idir                      ! if >=0  rspace ==> gspace
                                                   ! if <0   rspace <== gspace
  integer, intent(in) :: nfn                      ! number of functions to FFT
  type(fft_struc), intent(inout) :: ffts             ! work arrays for the FFT
  type(parallel_gspace), intent(in) :: &
       gs                            ! the gspace for which the data is set up
  !
  !     INPUT/OUTPUT:
  !     ------------
  !
  complex(dp), intent(inout), target :: &
       rspace(gs%r_size * nfn)     ! the data in realspace rspace(ffts%r_size)
  complex(dp), intent(inout) :: &
       gsdat(gs%length * nfn)            ! the data in gspace gsdat(gs%length)
  !
  !     DESCRIPTION:
  !     -----------
  !
  !     performs the first two stages of a forward or an inverse 3dFFT
  !
  !
  !     ---------------------- local variables -------------------------
  !
  integer :: &
       i, j, k, iord, nrods, b, b2, na, na2, rodrange(2), &
       nx, ny, nz, nxp, ireqr, ireqs, ip, ierr, irod, ipack, iofy, &
       psource, ptarget, packinfo, &
       ldx, ldy, ldz, &     ! leading dimensions in x, y, z
       nnx, nny, nnz, &     ! sizes of the fourier grid
       nrodst1, &           ! number of rods after first transposition
       nrodst2, &           ! number of rods after second transposition
       dnx, &               ! number of nx slabs after t1
       ioffset, &           ! rod offset for myproc on nz/ny grid
       srod, sy, slen, swidth, drod, dx, dlen, ix, ichunk, &
       poffset, dwidth, &   !
       npacksend, &         ! counter for num pack send
       dwx, dwy, dwz, dw2,& ! pointers into the work array
       ifn                  ! labels function
 
  integer :: istatus(MPI_STATUS_SIZE)

  real(dp) :: t0, t1, t2
  complex(dp) :: dummy

  
      !
      ! edit davegp
      ! I changed this parameter setting so that you can
      ! link to the correct header associated with the library
      ! that you link. Requires an include directory to be
      ! added to conf.mk
      !
      include 'fftw_f77.i'
      !
      ! INTEGER FFTW_FORWARD,FFTW_BACKWARD
      ! PARAMETER (FFTW_FORWARD=-1,FFTW_BACKWARD=1)
      ! INTEGER FFTW_ESTIMATE,FFTW_MEASURE
      ! PARAMETER (FFTW_ESTIMATE=0,FFTW_MEASURE=1)
      ! INTEGER FFTW_IN_PLACE,FFTW_USE_WISDOM
      ! PARAMETER (FFTW_IN_PLACE=8,FFTW_USE_WISDOM=16)

! for 64-bit machines (integer is not enough)
! it should also work on 32-bit machines

      ! the plan tags from create_plan
      INTEGER*8 FFTW_PLAN_FWD, FFTW_PLAN_BWD


  real(dp), external :: gimmetime


  complex(dp), pointer :: t1buf(:), ufb(:), sendbuf(:)
  t1buf => ffts%t1buf
  ufb => rspace
  sendbuf => ffts%sendbuf


  nnx = gs%fftsize(1) ; nny = gs%fftsize(2) ; nnz = gs%fftsize(3)
  dnx = mod(gs%fftsize(4) - gs%fftsize(5) + nnx, nnx) + 1

  nrodst1 = (dnx * nnz) / gs%nproc 
  if (mod(dnx * nnz, gs%nproc) > gs%myproc) nrodst1 = nrodst1 + 1

  nrodst2 = (nny * nnz) / gs%nproc 
  if (mod(nny * nnz, gs%nproc) > gs%myproc) nrodst2 = nrodst2 + 1
  !
  !     variables for local T2 transpose. See setup_packinfo
  !
  b = (dnx * nnz) / gs%nproc ; na = mod(dnx * nnz, gs%nproc)
  b2 = (nny * nnz) / gs%nproc ; na2 = mod(nny * nnz, gs%nproc)
  !
  rodrange(1) = b2 * gs%myproc + min(gs%myproc, na2)
  rodrange(2) = rodrange(1) + b2
  if (na2 <= gs%myproc) rodrange(2) = rodrange(2) - 1
  !
  ioffset = b * gs%myproc + min(gs%myproc, na)
  !
  !     do some sanity checks
  !
  if (ffts%fftsize(1) /= gs%fftsize(1) .or. &
       ffts%fftsize(2) /= gs%fftsize(2) .or. &
       ffts%fftsize(3) /= gs%fftsize(3))then 
     write(9, *) 'fftsize mismatch between gspace:', &
          gs%fftsize(1:3), 'and ffts:', ffts%fftsize
     call mystop
  end if
  !
  !     figure out leading dimensions for the cray optimization
  !
  ldx = nnx ; ldy = nny ; ldz = nnz



  if (ldz * gs%lorder > gs%r_size) then
     write(9, *) ' FFT: strange geometry used...'
     write(9, *) ' change code not to abuse rspace as buffer!'
     call mystop
  end if
  !
  !     compute pointers into the work array 
  !
  dwx = 1 ; dwy = ffts%naux(1) + 1 ; dwz = ffts%naux(1) + ffts%naux(2) + 1
  dw2 = ffts%naux(1) + ffts%naux(2) + ffts%naux(3) + 1

  if (idir < 0) then ! from gspace to realspace
     !
     !     ================   FROM GSPACE TO REALSPACE ===================
     !                        (inverse or backward FFT)
     !
     !        Set the data up for the first fourier transform.
     !        All the rods on this processor are put next to each other
     !        into the unfold buffer, and then fourier transformed with 
     !        a block routine for optimum performance.
     ! 

     ! pdh print out input data (gspace)
     !      do ifn=0,nfn-1
     !      write(*,*) 'WFN:',ifn+1,' of ',nfn
     !      write(*,*) gsdat(1+ifn*gs%length:(ifn+1)*gs%length)
     !      end do

     

  

     call fftunfold(gsdat(1), gs, ufb(1), ldz, nfn)
  

     
     ! pdh print out unfold buffer
     !
     !      do ifn=0,nfn-1
     !      write(*,*) 'WFN:',ifn+1,' of ',nfn
     !      write(*,*) 'ldz=',ldz
     !      write(*,*) 'gs%lorder=',gs%lorder
     !      write(*,*) 'gs%length=',gs%length
     !      write(*,*) ufb(1+ifn*ldz*gs%lorder:ldz*gs%lorder*(ifn+1))
     !      end do

     !
     !     ---- inverse fourier transform in z direction locally ----
     !
     


     







     
      !
      ! edit davegp
      ! new function call which provides existing plans
      ! if these dimensions/directions have been used before
      !
      FFTW_PLAN_FWD = fftw_hitlist_plan( nnz,  1 )
      FFTW_PLAN_BWD = fftw_hitlist_plan( nnz, -1 )
      !
      ! This was done before
      !
      !call fftw_f77_create_plan(FFTW_PLAN_FWD, &
      !   nnz,FFTW_FORWARD,FFTW_MEASURE+FFTW_IN_PLACE+FFTW_USE_WISDOM)
      !call fftw_f77_create_plan(FFTW_PLAN_BWD, &
      !   nnz,FFTW_BACKWARD,FFTW_MEASURE+FFTW_IN_PLACE+FFTW_USE_WISDOM)
      ! I just changed the above from FFTW_ESTIMATE to FFTW_MEASURE
      !

     
     
      !
      call fftw_f77(FFTW_PLAN_BWD, gs%lorder * nfn, ufb(1),1, ldz, ffts%workspace(dwz),1, ldz)
      !


     
      !
      ! edit davegp
      ! I do not destroy old plans anymore since I keep them
      ! for future use by FFTs on arrays of the same dimension
      !
      ! This was done before
      !
      !call fftw_f77_destroy_plan(FFTW_PLAN_FWD)
      !call fftw_f77_destroy_plan(FFTW_PLAN_BWD)
      !

     

     

     ! pdh print out unfold buffer after zFFT
     !
     !      do ifn=0,nfn-1
     !      write(*,*) 'WFN:',ifn+1,' of ',nfn
     !      write(*,*) 'ldz=',ldz
     !      write(*,*) 'gs%lorder=',gs%lorder
     !      write(*,*) 'gs%length=',gs%length
     !      write(*,*) ufb(1+ifn*ldz*gs%lorder:ldz*gs%lorder*(ifn+1))
     !      end do


  
     t1buf(1:nrodst1 * ldy * nfn) = zzero
  


     
     !
     !        do local transpose without any message passing
     !
     

     do ifn = 0, nfn - 1
        call frepack(gs%packsize(gs%myproc + 1, 2, 1), &
             t1buf(ifn * ldy * nrodst1 + 1), &
             gs%packinfo(1, gs%myproc + 1, 2), &
             ldy * nrodst1, ufb(ifn * ldz * gs%lorder + 1), &
             gs%packinfo(1, gs%myproc + 1, 1), ldz * gs%lorder)
     end do

     

     ! pdh print out t1buf buffer after local transpose t1
     !
     !      do ifn=0,nfn-1
     !      write(*,*) 'WFN:',ifn+1,' of ',nfn
     !      write(*,*) 'ldy=',ldy
     !      write(*,*) 'nrodst1=',nrodst1
     !      write(*,*) t1buf(1+ifn*ldy*nrodst1:ldy*nrodst1*(ifn+1))
     !      end do

     !
     !        communication phase to transpose. New grid will be Nx/Nz
     !

     
  
     do ip = 1, gs%nproc - 1
        ptarget = mod(ip + gs%myproc, gs%nproc)         ! the target processor
        psource = mod(gs%myproc + gs%nproc - ip, gs%nproc)  ! source processor

  
         !
         !        post a non-blocking receive request to the source processor
         !
        if (gs%packsize(psource + 1, 2, 1) > 0) &
             call mpi_irecv(ffts%recbuf(1), ffts%recbufsize * nfn, &
             MPI_DOUBLE_COMPLEX, psource, psource, MPI_COMM_WORLD, ireqr, ierr)
        !
        !           generate and send off the package for the target processor
        !
    
        do ifn = 0, nfn - 1
           do ipack = 1, gs%packsize(ptarget + 1, 1, 1)
              sendbuf(ipack + ifn * gs%packsize(ptarget + 1, 1, 1)) = &
                   ufb(gs%packinfo(ipack, ptarget + 1, 1) + &
                   ifn * ldz * gs%lorder)
           end do
        end do

        if (gs%packsize(ptarget + 1, 1, 1) > 0) &
             call mpi_isend(sendbuf(1), &
             gs%packsize(ptarget + 1, 1, 1) * nfn, MPI_DOUBLE_COMPLEX, &
             ptarget, gs%myproc, MPI_COMM_WORLD, ireqs, ierr)
    
        !     
        !           unpack the package from the neighbor. 
        !
        if (gs%packsize(psource + 1, 2, 1) > 0) &
             call mpi_wait(ireqr, istatus, ierr)      ! wait until packet here
        do ifn = 0, nfn - 1
           do ipack = 1, gs%packsize(psource + 1, 2, 1)
              t1buf(gs%packinfo(ipack, psource + 1, 2) + &
                   ifn * ldy * nrodst1) = &
                   ffts%recbuf(ipack + ifn * gs%packsize(psource + 1, 2, 1))
           end do
        end do
        !
        !              wait until send has completed
        !
        if (gs%packsize(ptarget + 1, 1, 1) > 0) &
             call mpi_wait(ireqs, istatus, ierr)
  
     end do
  
     

     
     !
     !        inverse fourier transform in y direction locally 
     !


     







     
      !
      ! edit davegp
      ! new function call which provides existing plans
      ! if these dimensions/directions have been used before
      !
      FFTW_PLAN_FWD = fftw_hitlist_plan( nny,  1 )
      FFTW_PLAN_BWD = fftw_hitlist_plan( nny, -1 )
      !
      ! This was done before
      !
      !call fftw_f77_create_plan(FFTW_PLAN_FWD, &
      !   nny,FFTW_FORWARD,FFTW_MEASURE+FFTW_IN_PLACE+FFTW_USE_WISDOM)
      !call fftw_f77_create_plan(FFTW_PLAN_BWD, &
      !   nny,FFTW_BACKWARD,FFTW_MEASURE+FFTW_IN_PLACE+FFTW_USE_WISDOM)
      ! I just changed the above from FFTW_ESTIMATE to FFTW_MEASURE
      !

     
      !
      call fftw_f77(FFTW_PLAN_BWD, nrodst1 * nfn, t1buf(1),1, ldy, ffts%workspace(dwy),1, ldy)
      !

     
      !
      ! edit davegp
      ! I do not destroy old plans anymore since I keep them
      ! for future use by FFTs on arrays of the same dimension
      !
      ! This was done before
      !
      !call fftw_f77_destroy_plan(FFTW_PLAN_FWD)
      !call fftw_f77_destroy_plan(FFTW_PLAN_BWD)
      !

     

     ! pdh print out t1buf buffer after yFFT
     !
     !      do ifn=0,nfn-1
     !      write(*,*) 'WFN:',ifn+1,' of ',nfn
     !      write(*,*) 'ldy=',ldy
     !      write(*,*) 'nrodst1=',nrodst1
     !      write(*,*) t1buf(1+ifn*ldy*nrodst1:ldy*nrodst1*(ifn+1))
     !      end do

     !
     !     -------- second transposition phase. grid will be Ny/Nz  ------------
     !                                                       Ny  running faster
     !

     !
     !        clear the array if no fast convolution is wanted
     !
     
     iofy = gs%fftsize(5) - gs%fftsize(4)



     
     
     !
     !        do local transpose
     !
     if (gs%packsize(gs%myproc + 1, 1, 2) > 0) then
        do ifn = 0, nfn - 1

           call fftlt2(-1, gs%chunk(1, 1, 1, gs%myproc + 1), &
                8 * gs%maxnchunk, &
  rspace(1 + ifn * ldx * nrodst2), &
                ldx * nrodst2, t1buf(1 + ifn * ldy * nrodst1), ldy * nrodst1, &
                nnx, nny, ldx, ldy)

        end do
     endif
     

     ! pdh print out rspace buffer after local transpose t2
     !
     !      do ifn=0,nfn-1
     !      write(*,*) 'WFN:',ifn+1,' of ',nfn
     !      write(*,*) 'ldx=',ldx
     !      write(*,*) 'nrodst2=',nrodst2
     !      write(*,*) rspace(1+ifn*ldx*nrodst2:ldx*nrodst2*(ifn+1))
     !      end do

     !
     !        communicate with the other processors.
     !


     
     !
     !        First pack all the packets to all processors
     !        In case mpi is used, also send out immediately
     !
     poffset = 1 ; npacksend = 0
     do ip = 1, gs%nproc - 1
        ptarget = mod(ip + gs%myproc, gs%nproc)         ! the target processor
        if (gs%packsize(ptarget + 1, 1, 2) > 0) then
           ichunk = 1
           ipack = poffset
           do while (.true.)
              srod    = gs%chunk(1, 1, ichunk, ptarget + 1)
              sy      = gs%chunk(2, 1, ichunk, ptarget + 1) + 1
              slen    = gs%chunk(3, 1, ichunk, ptarget + 1) - 1
              swidth  = abs(gs%chunk(4, 1, ichunk, ptarget + 1))
              do ifn = 0, nfn - 1
                 do irod = srod, srod + swidth - 1
                   call mzcopy(slen + 1, t1buf(irod * ldy + sy), 1, &
                    sendbuf(poffset), 1 ) 
!                    sendbuf(poffset:poffset + slen) = &
!                         t1buf(irod * ldy + sy:irod * ldy + sy + slen)
                    poffset = poffset + slen + 1
                 end do
                 sy = sy + ldy * nrodst1
              end do
              if (gs%chunk(4, 1, ichunk, ptarget + 1) < 0) exit
              ichunk = ichunk + 1
           end do
  
           npacksend = npacksend + 1
           call mpi_isend(sendbuf(ipack), &
                nfn * gs%packsize(ptarget + 1, 1, 2), &
                MPI_DOUBLE_COMPLEX, ptarget, gs%myproc, MPI_COMM_WORLD, &
                ffts%mpireq(npacksend), ierr)
  
        end if
     end do                 ! loop over all packets to other processors


  
     !
     !        now get the data from the other processor
     !
     !
     !        The  code has room for improvement here. One should just
     !        pull in the messages, read the tag, and dispatch accordingly.
     !
     !
     do ip = 1, gs%nproc - 1
        psource = mod(gs%myproc + gs%nproc - ip, gs%nproc) ! source processor
        if (gs%packsize(psource + 1, 2, 2) > 0) then
  
           !              for , do a blocking receive on the source proc
           call mpi_recv(ffts%recbuf(1), ffts%recbufsize * nfn, &
                MPI_DOUBLE_COMPLEX, psource, psource, &
                MPI_COMM_WORLD, istatus, ierr)
  
           !
           !              unpack the received data
           !
           ichunk = 1
           poffset = 1
           do while (.true.)
              drod    = gs%chunk(1, 2, ichunk, psource +1 )
              dx      = gs%chunk(2, 2, ichunk, psource +1 )
              dlen    = gs%chunk(3, 2, ichunk, psource +1 ) - 1
              dwidth  = abs(gs%chunk(4, 2, ichunk, psource + 1))
              do ifn = 0 , nfn - 1
                 do j = 0, dwidth - 1
                    ix = mod(dx + j, nnx) + 1 + ifn * ldx * nrodst2
  
                     call mzcopy(dlen + 1, ffts%recbuf(poffset), 1, &
                        rspace(drod * ldx + ix), ldx ) 
!                    rspace(drod * ldx + ix:(drod + dlen) * ldx + ix: &
!                         ldx) = ffts%recbuf(poffset:poffset + dlen)
  
                    poffset = poffset + dlen + 1
                 end do
              end do
              if (gs%chunk(4, 2, ichunk, psource + 1) < 0) exit
              ichunk = ichunk + 1
           end do
        end if
     end do                 ! loop over processors
  
     if (npacksend > 0) call mpi_waitall(npacksend, &
          ffts%mpireq(1), ffts%mpistat(1), ierr) 
  
     

  else
     !
     !     ================   FROM REALSPACE TO GSPACE ===================
     !
     !                        (forward FFT)
     !
     !
     !     ---------------- second transposition phase -----------------------
     !                 (really the first transpose for the forward FFT)

     !
     !     do local transpose without any message passing
     !

     ! pdh print out rspace buffer before local transpose t2
     !
     !      do ifn=0,nfn-1
     !      write(*,*) 'WFN:',ifn+1,' of ',nfn
     !      write(*,*) 'ldx=',ldx
     !      write(*,*) 'nrodst2=',nrodst2
     !      write(*,*) rspace(1+ifn*ldx*nrodst2:ldx*nrodst2*(ifn+1))
     !      end do

     
     if (gs%packsize(gs%myproc + 1, 1, 2) > 0) then
        do ifn = 0, nfn - 1

           call fftlt2(1, gs%chunk(1, 1, 1, gs%myproc + 1), &
                8*gs%maxnchunk, &
  rspace(1 + ifn * ldx * nrodst2), &
                ldx * nrodst2, ffts%t1buf(1 + ifn * ldy * nrodst1), &
                ldy * nrodst1, nnx, nny, ldx, ldy)

        end do
     end if
     
     !
     !        now communicate (t2 forward FFT)
     !
     

     !
     !        First pack all the packets to all processors
     !        In case mpi is used, also send out immediately
     !
     poffset = 1 ; npacksend = 0
     do ip = 1, gs%nproc - 1
        ptarget = mod(ip + gs%myproc, gs%nproc)         ! the target processor
        if (gs%packsize(ptarget + 1, 2, 2) > 0) then
           ichunk = 1
           ipack = poffset
           do while(.true.)
              drod    = gs%chunk(1, 2, ichunk, ptarget + 1)
              dx      = gs%chunk(2, 2, ichunk, ptarget + 1)
              dlen    = gs%chunk(3, 2, ichunk, ptarget + 1) - 1
              dwidth  = abs(gs%chunk(4, 2, ichunk, ptarget + 1))
              do ifn = 0, nfn - 1
                 do j = 0, dwidth - 1
                    ix = mod(dx + j, nnx) + 1 + ifn * ldx * nrodst2
  
                     call mzcopy(dlen + 1, rspace(drod * ldx + ix), ldx, &
                            sendbuf(poffset), 1 ) 
  
                    poffset = poffset + dlen + 1
                 end do
              end do
              if (gs%chunk(4, 2, ichunk, ptarget + 1) < 0) exit
              ichunk = ichunk + 1
           end do

  
           npacksend = npacksend + 1
           call mpi_isend(sendbuf(ipack), &
                nfn * gs%packsize(ptarget + 1, 2, 2), &
                MPI_DOUBLE_COMPLEX, ptarget, gs%myproc, &
                MPI_COMM_WORLD, ffts%mpireq(npacksend), ierr)
  
        end if
     end do                 ! loop over all packets to other procs
     !
     !        now get the data from the other processor
     !
     !        The  code has room for improvement here. One should just
     !        pull in the messages, read the tag, and dispatch accordingly.
     !
     !
  
     do ip = 1, gs%nproc - 1
        psource = mod(gs%myproc + gs%nproc - ip, gs%nproc)   ! the source proc
        if (gs%packsize(psource + 1, 1, 2) > 0) then
  
           !              for , do a blocking receive on the source proc
           call mpi_recv(ffts%recbuf(1), ffts%recbufsize * nfn, &
                MPI_DOUBLE_COMPLEX, psource, psource, &
                MPI_COMM_WORLD, istatus, ierr)
  
           !
           !         unpack the package from the neighbor. The packinfo array
           !         gives the index for the incoming data.
           !
           ichunk = 1
           poffset = 1
           do while (.true.)
              srod    = gs%chunk(1, 1, ichunk, psource + 1)
              sy      = gs%chunk(2, 1, ichunk, psource + 1) + 1
              slen    = gs%chunk(3, 1, ichunk, psource + 1) - 1
              swidth  = abs(gs%chunk(4, 1, ichunk, psource + 1))
              do ifn = 0, nfn - 1
                 do irod = srod, srod + swidth - 1
                    call mzcopy(slen + 1, ffts%recbuf(poffset), 1, &
                      ffts%t1buf(irod * ldy + sy), 1 ) 
                    poffset = poffset + slen + 1
                 end do
                 sy = sy + ldy * nrodst1
              end do
              if (gs%chunk(4, 1, ichunk, psource + 1) < 0) exit
              ichunk = ichunk + 1
           end do
        end if
     end do                 ! loop over processors

     !        make sure the send buffer can be deallocated
  
     if (npacksend > 0) call mpi_waitall(npacksend, &
          ffts%mpireq(1), ffts%mpistat(1), ierr) 
  
     


     ! pdh print out t1buf buffer after yFFT
     !
     !      do ifn=0,nfn-1
     !      write(*,*) 'WFN:',ifn+1,' of ',nfn
     !      write(*,*) 'ldy=',ldy
     !      write(*,*) 'nrodst1=',nrodst1
     !      write(*,*) t1buf(1+ifn*ldy*nrodst1:ldy*nrodst1*(ifn+1))
     !      end do

     !
     !        fourier transform in y direction locally 
     !
     


     

     







     
      !
      ! edit davegp
      ! new function call which provides existing plans
      ! if these dimensions/directions have been used before
      !
      FFTW_PLAN_FWD = fftw_hitlist_plan( nny,  1 )
      FFTW_PLAN_BWD = fftw_hitlist_plan( nny, -1 )
      !
      ! This was done before
      !
      !call fftw_f77_create_plan(FFTW_PLAN_FWD, &
      !   nny,FFTW_FORWARD,FFTW_MEASURE+FFTW_IN_PLACE+FFTW_USE_WISDOM)
      !call fftw_f77_create_plan(FFTW_PLAN_BWD, &
      !   nny,FFTW_BACKWARD,FFTW_MEASURE+FFTW_IN_PLACE+FFTW_USE_WISDOM)
      ! I just changed the above from FFTW_ESTIMATE to FFTW_MEASURE
      !

     
      !
      call fftw_f77(FFTW_PLAN_FWD, nrodst1 * nfn, ffts%t1buf(1),1, ldy, ffts%workspace(dwy),1, ldy)
      !

     
      !
      ! edit davegp
      ! I do not destroy old plans anymore since I keep them
      ! for future use by FFTs on arrays of the same dimension
      !
      ! This was done before
      !
      !call fftw_f77_destroy_plan(FFTW_PLAN_FWD)
      !call fftw_f77_destroy_plan(FFTW_PLAN_BWD)
      !

     
     !
     !     --------------------- first transpose --------------------
     !
     !
     !        do local transpose without any message passing
     !
     

     do ifn = 0, nfn - 1
        call frepack(gs%packsize(gs%myproc + 1, 2, 1), &
             ufb(ifn * ldz * gs%lorder + 1), &
             gs%packinfo(1, gs%myproc+1, 1), ldz * gs%lorder, &
             ffts%t1buf(ifn * ldy * nrodst1 + 1), &
             gs%packinfo(1, gs%myproc + 1, 2), ldy * nrodst1)
     end do

     


     
  
     do ip = 1, gs%nproc - 1
        ptarget = mod(ip + gs%myproc, gs%nproc)         ! the target processor
        psource = mod(gs%myproc + gs%nproc - ip, gs%nproc)  ! source processor
  
        !
        !           post a non-blocking receive request to the source processor
        !
        if (gs%packsize(psource + 1, 1, 1) > 0) &
             call mpi_irecv(ffts%recbuf(1), ffts%recbufsize * nfn, &
             MPI_DOUBLE_COMPLEX, psource, psource, &
             MPI_COMM_WORLD, ireqr, ierr)
        !
        !           generate and send off the package for the target processor
        !
        if (gs%packsize(ptarget + 1, 2, 1) > 0) then
           do ifn = 0, nfn - 1
              do ipack = 1, gs%packsize(ptarget + 1, 2, 1)
                 sendbuf(ipack + ifn * gs%packsize(ptarget + 1, 2, 1)) = &
                      ffts%t1buf(gs%packinfo(ipack, ptarget + 1, 2) + &
                      ifn * ldy * nrodst1)
              end do
           end do
           call mpi_isend(sendbuf(1), gs%packsize(ptarget + 1, 2, 1) * nfn, &
                MPI_DOUBLE_COMPLEX, ptarget, gs%myproc, &
                MPI_COMM_WORLD, ireqs, ierr)
        end if
        !
        !           unpack the package from the neighbor. 
        !
        if (gs%packsize(psource + 1, 1, 1) > 0) then
           call mpi_wait(ireqr, istatus, ierr)        ! wait until packet here
           do ifn = 0, nfn - 1
              do ipack = 1, gs%packsize(psource + 1, 1, 1)
                 ufb(gs%packinfo(ipack, psource + 1, 1) + ifn * ldz * &
                      gs%lorder) = ffts%recbuf(ipack + ifn * &
                      gs%packsize(psource + 1, 1, 1))
              end do
           end do
        end if
        !
        !           wait until send has completed
        !
        if (gs%packsize(ptarget + 1, 2, 1) > 0) &
             call mpi_wait(ireqs, istatus, ierr)
  
     end do
  
     

     !
     !     ---- fourier transform in z direction locally ----
     !
     
     
     







     
      !
      ! edit davegp
      ! new function call which provides existing plans
      ! if these dimensions/directions have been used before
      !
      FFTW_PLAN_FWD = fftw_hitlist_plan( nnz,  1 )
      FFTW_PLAN_BWD = fftw_hitlist_plan( nnz, -1 )
      !
      ! This was done before
      !
      !call fftw_f77_create_plan(FFTW_PLAN_FWD, &
      !   nnz,FFTW_FORWARD,FFTW_MEASURE+FFTW_IN_PLACE+FFTW_USE_WISDOM)
      !call fftw_f77_create_plan(FFTW_PLAN_BWD, &
      !   nnz,FFTW_BACKWARD,FFTW_MEASURE+FFTW_IN_PLACE+FFTW_USE_WISDOM)
      ! I just changed the above from FFTW_ESTIMATE to FFTW_MEASURE
      !

     
      !
      call fftw_f77(FFTW_PLAN_FWD, gs%lorder * nfn, ufb(1),1, ldz, ffts%workspace(dwz),1, ldz)
      !

     
      !
      ! edit davegp
      ! I do not destroy old plans anymore since I keep them
      ! for future use by FFTs on arrays of the same dimension
      !
      ! This was done before
      !
      !call fftw_f77_destroy_plan(FFTW_PLAN_FWD)
      !call fftw_f77_destroy_plan(FFTW_PLAN_BWD)
      !


     

     !        Set the data up for the first fourier transform.
     !        All the rods on this processor are put next to each other
     !        into the unfold buffer, and then fourier transformed with 
     !        a block routine for optimum performance.
     ! 
     
     call fftfold(gsdat(1), gs, ufb(1), ldz, nfn)
     
  end if




  return

100 format('TIME FOR ',a30,':',f12.6)

end subroutine fourier_transform_12
!
!     ======================================================================
!
subroutine fourier_transform_3(idir, ffts, gs, rspace, nfn)
  !
  ! edit davegp
  use fftw_hitlist, only: fftw_hitlist_plan
  !
  !
  include 'use.h'
  implicit none             ! implicit? Just say no!
include 'mpif.h'
  include 'interface.h'

  !
  !     INPUT:
  !     -----
  !
  integer, intent(in) :: idir                     ! if >=0  rspace ==> gspace
                                                  ! if <0   rspace <== gspace
  integer, intent(in) :: nfn                     ! number of functions to FFT
  type(fft_struc), intent(in) :: ffts               ! work arrays for the FFT
  type(parallel_gspace), intent(in) :: &
       gs                           ! the gspace for which the data is set up
  !
  !     INPUT/OUTPUT:
  !     ------------
  !
  complex(dp), intent(inout) :: &
       rspace(gs%r_size * nfn)! the data in realspace rspace(ffts%r_size)
  !                                 on the cray pvp, this is untouched, and the
  !                                 output is in ffts%unfoldbuf
  !
  !     DESCRIPTION:
  !     -----------
  !
  !     performs the first stage of a forward or the last stage of an
  !     inverse 3dFFT.
  !
  !     ON THE VECTOR MACHINES, THE OUTPUT IS ACTUALLY RETURNED IN 
  !     THE FFTS%UNFOLDBUF ARRAY. This is because on the vector
  !     machines, the stride should not be a power of two, but that is
  !     likely to be the case for a nice FFT grid.
  !
  !
  !
  integer :: i, j, nnx, nny, nnz, nrodst2, dwx, dw2, ldx, ldy, ldz, &
       nz, ny, k, ierr, ifn
  real(dp) :: t0
  real(dp), external :: gimmetime

  
      !
      ! edit davegp
      ! I changed this parameter setting so that you can
      ! link to the correct header associated with the library
      ! that you link. Requires an include directory to be
      ! added to conf.mk
      !
      include 'fftw_f77.i'
      !
      ! INTEGER FFTW_FORWARD,FFTW_BACKWARD
      ! PARAMETER (FFTW_FORWARD=-1,FFTW_BACKWARD=1)
      ! INTEGER FFTW_ESTIMATE,FFTW_MEASURE
      ! PARAMETER (FFTW_ESTIMATE=0,FFTW_MEASURE=1)
      ! INTEGER FFTW_IN_PLACE,FFTW_USE_WISDOM
      ! PARAMETER (FFTW_IN_PLACE=8,FFTW_USE_WISDOM=16)

! for 64-bit machines (integer is not enough)
! it should also work on 32-bit machines

      ! the plan tags from create_plan
      INTEGER*8 FFTW_PLAN_FWD, FFTW_PLAN_BWD



  nnx = gs%fftsize(1) ; nny = gs%fftsize(2) ; nnz = gs%fftsize(3)
  nrodst2 = (nny * nnz) / gs%nproc 
  if (mod(nny * nnz, gs%nproc) > gs%myproc) nrodst2 = nrodst2 + 1
  dwx = 1 ; dw2 = ffts%naux(1) + ffts%naux(2) + ffts%naux(3) + 1
  ldx = nnx ; ldy = nny ; ldz = nnz




  if (idir < 0) then ! from gspace to realspace
     !
     !           inverse fourier transform in x direction locally
     !
     !     
     
     








     
      !
      ! edit davegp
      ! new function call which provides existing plans
      ! if these dimensions/directions have been used before
      !
      FFTW_PLAN_FWD = fftw_hitlist_plan( nnx,  1 )
      FFTW_PLAN_BWD = fftw_hitlist_plan( nnx, -1 )
      !
      ! This was done before
      !
      !call fftw_f77_create_plan(FFTW_PLAN_FWD, &
      !   nnx,FFTW_FORWARD,FFTW_MEASURE+FFTW_IN_PLACE+FFTW_USE_WISDOM)
      !call fftw_f77_create_plan(FFTW_PLAN_BWD, &
      !   nnx,FFTW_BACKWARD,FFTW_MEASURE+FFTW_IN_PLACE+FFTW_USE_WISDOM)
      ! I just changed the above from FFTW_ESTIMATE to FFTW_MEASURE
      !

     do ifn = 0, nfn - 1
        do i = 0, nrodst2 - 1
           j = i * nnx + 1 + ifn * nrodst2 * nnx
           if (gs%fftsize(5) > gs%fftsize(4) + 1) &
                rspace(j + gs%fftsize(4) + 1:j + gs%fftsize(5) - 1) = zzero
           
      !
      call fftw_f77(FFTW_PLAN_BWD,1,rspace(j),1,1,ffts%workspace(dwx),1,1)
      !

        end do
     end do

     
      !
      ! edit davegp
      ! I do not destroy old plans anymore since I keep them
      ! for future use by FFTs on arrays of the same dimension
      !
      ! This was done before
      !
      !call fftw_f77_destroy_plan(FFTW_PLAN_FWD)
      !call fftw_f77_destroy_plan(FFTW_PLAN_BWD)
      !




     

  else                      ! -------------- forward FFT ------------

     !
     !        fourier transform in x direction locally
     !
     
     







     
      !
      ! edit davegp
      ! new function call which provides existing plans
      ! if these dimensions/directions have been used before
      !
      FFTW_PLAN_FWD = fftw_hitlist_plan( nnx,  1 )
      FFTW_PLAN_BWD = fftw_hitlist_plan( nnx, -1 )
      !
      ! This was done before
      !
      !call fftw_f77_create_plan(FFTW_PLAN_FWD, &
      !   nnx,FFTW_FORWARD,FFTW_MEASURE+FFTW_IN_PLACE+FFTW_USE_WISDOM)
      !call fftw_f77_create_plan(FFTW_PLAN_BWD, &
      !   nnx,FFTW_BACKWARD,FFTW_MEASURE+FFTW_IN_PLACE+FFTW_USE_WISDOM)
      ! I just changed the above from FFTW_ESTIMATE to FFTW_MEASURE
      !


     
      !
      call fftw_f77(FFTW_PLAN_FWD, nrodst2 * nfn, rspace(1),1, ldx, ffts%workspace(dwx),1, ldx)
      !


     
      !
      ! edit davegp
      ! I do not destroy old plans anymore since I keep them
      ! for future use by FFTs on arrays of the same dimension
      !
      ! This was done before
      !
      !call fftw_f77_destroy_plan(FFTW_PLAN_FWD)
      !call fftw_f77_destroy_plan(FFTW_PLAN_BWD)
      !




     

  end if


  return

end subroutine fourier_transform_3
!
!     ======================================================================
!
subroutine fast_convolute(ffts, gs, gsdat, vgsdat, rspace, vloc)
  !
  ! edit davegp
  use fftw_hitlist, only: fftw_hitlist_plan
  !
  !
  include 'use.h'
  implicit none             ! implicit? Just say no!
include 'mpif.h'
  include 'interface.h'

  !
  !     INPUT:
  !     -----
  !
  type(fft_struc), intent(in) :: ffts                ! work arrays for the FFT
  type(parallel_gspace), intent(in) :: &
       gs                            ! the gspace for which the data is set up
  complex(dp), intent(in) :: vloc(ffts%r_size)           ! realspace potential
  !
  !     INPUT/OUTPUT:
  !     ------------
  !
  complex(dp), intent(inout) :: &
       rspace(gs%r_size)           ! the data in realspace rspace(ffts%r_size)
  complex(dp), intent(inout) :: &
       gsdat(gs%length)                  ! the data in gspace gsdat(gs%length)
  !
  !     OUTPUT:
  !     ------
  !
  complex(dp), intent(inout) :: &
       vgsdat(gs%length)            ! the data in gspace multiplied with vloc
  !                                                         vgsdat(gs%length)
  !
  !     DESCRIPTION:
  !     -----------
  !
  !     performs the last stage of the inverse FFT, multiplies with
  !     vloc, and does the first stage of the forward FFT
  !
  integer :: &
       i, j, k, nnx, nny, nnz, nrodst2, dwx, dw2, ldx, ldy, ldz, nblock, n, &
       ifn
  real(dp) :: t0
  complex(dp), allocatable :: work(:)
  real(dp), external :: gimmetime

  
      !
      ! edit davegp
      ! I changed this parameter setting so that you can
      ! link to the correct header associated with the library
      ! that you link. Requires an include directory to be
      ! added to conf.mk
      !
      include 'fftw_f77.i'
      !
      ! INTEGER FFTW_FORWARD,FFTW_BACKWARD
      ! PARAMETER (FFTW_FORWARD=-1,FFTW_BACKWARD=1)
      ! INTEGER FFTW_ESTIMATE,FFTW_MEASURE
      ! PARAMETER (FFTW_ESTIMATE=0,FFTW_MEASURE=1)
      ! INTEGER FFTW_IN_PLACE,FFTW_USE_WISDOM
      ! PARAMETER (FFTW_IN_PLACE=8,FFTW_USE_WISDOM=16)

! for 64-bit machines (integer is not enough)
! it should also work on 32-bit machines

      ! the plan tags from create_plan
      INTEGER*8 FFTW_PLAN_FWD, FFTW_PLAN_BWD



  nnx = gs%fftsize(1) ; nny = gs%fftsize(2) ; nnz = gs%fftsize(3)
  nrodst2 = (nny * nnz) / gs%nproc 
  if (mod(nny * nnz, gs%nproc) > gs%myproc) nrodst2 = nrodst2 + 1
  dwx = 1 ; dw2 = ffts%naux(1) + ffts%naux(2) + ffts%naux(3) + 1
  nblock = 1
  ldx = nnx ; ldy = nny ; ldz = nnz



  

  nblock = 2



  







  
      !
      ! edit davegp
      ! new function call which provides existing plans
      ! if these dimensions/directions have been used before
      !
      FFTW_PLAN_FWD = fftw_hitlist_plan( nnx,  1 )
      FFTW_PLAN_BWD = fftw_hitlist_plan( nnx, -1 )
      !
      ! This was done before
      !
      !call fftw_f77_create_plan(FFTW_PLAN_FWD, &
      !   nnx,FFTW_FORWARD,FFTW_MEASURE+FFTW_IN_PLACE+FFTW_USE_WISDOM)
      !call fftw_f77_create_plan(FFTW_PLAN_BWD, &
      !   nnx,FFTW_BACKWARD,FFTW_MEASURE+FFTW_IN_PLACE+FFTW_USE_WISDOM)
      ! I just changed the above from FFTW_ESTIMATE to FFTW_MEASURE
      !

  !
  !     now do the blocked loop
  !
  do i = 0, nrodst2 - nblock, nblock
     j = i * ldx + 1
     if (gs%fftsize(5) > gs%fftsize(4) + 1) then
        do n = 0, nblock - 1     ! clear the array
           rspace(n * ldx + j + gs%fftsize(4) + 1: &
                n * ldx + j + gs%fftsize(5) - 1) = zzero
        end do
     end if
     









     
      !
      call fftw_f77(FFTW_PLAN_BWD, nblock, rspace(1 + i * nnx),1, ldx, ffts%workspace(dwx),1, ldx)
      !


     !
     !        mult with realspace potential
     !
     
      do ifn=0,1-1
         rspace(j+(nnx * nblock*ifn): &
     	     j+(nnx * nblock*(ifn+1))-1)= &
             rspace(j+(nnx * nblock*ifn): &
             j+(nnx * nblock*(ifn+1))-1)* &
             vloc(j:j+(nnx * nblock)-1) 
      end do

     !
     !        generic loop for arbitrary block size
     !         do n=0,nblock-1    
     !            rspace(j+n*ldx:j+n*ldx+nnx-1)= &
     !             rspace(j+n*ldx:j+n*ldx+nnx-1)*vloc(j+n*ldx:j+n*ldx+nnx-1)
     !         end do
     !



     
      !
      call fftw_f77(FFTW_PLAN_FWD, nblock, rspace(1 + i * nnx),1, ldx, ffts%workspace(dwx),1, ldx)
      !


  end do
  
      !
      ! edit davegp
      ! I do not destroy old plans anymore since I keep them
      ! for future use by FFTs on arrays of the same dimension
      !
      ! This was done before
      !
      !call fftw_f77_destroy_plan(FFTW_PLAN_FWD)
      !call fftw_f77_destroy_plan(FFTW_PLAN_BWD)
      !




  







  
      !
      ! edit davegp
      ! new function call which provides existing plans
      ! if these dimensions/directions have been used before
      !
      FFTW_PLAN_FWD = fftw_hitlist_plan( nnx,  1 )
      FFTW_PLAN_BWD = fftw_hitlist_plan( nnx, -1 )
      !
      ! This was done before
      !
      !call fftw_f77_create_plan(FFTW_PLAN_FWD, &
      !   nnx,FFTW_FORWARD,FFTW_MEASURE+FFTW_IN_PLACE+FFTW_USE_WISDOM)
      !call fftw_f77_create_plan(FFTW_PLAN_BWD, &
      !   nnx,FFTW_BACKWARD,FFTW_MEASURE+FFTW_IN_PLACE+FFTW_USE_WISDOM)
      ! I just changed the above from FFTW_ESTIMATE to FFTW_MEASURE
      !

  !
  !     now finish off what is left
  !
  k = i
  do i = k, nrodst2 - 1
     j = i * nnx + 1
     if (gs%fftsize(5) > gs%fftsize(4) + 1) &
          rspace(j + gs%fftsize(4) + 1:j + gs%fftsize(5) - 1) = zzero
     







     
      !
      call fftw_f77(FFTW_PLAN_BWD,1,rspace(1 + nnx * i),1,1,ffts%workspace(dwx),1,1)
      !


     
      do ifn=0,1-1
         rspace(j+(nnx*ifn): &
     	     j+(nnx*(ifn+1))-1)= &
             rspace(j+(nnx*ifn): &
             j+(nnx*(ifn+1))-1)* &
             vloc(j:j+(nnx)-1) 
      end do



     
      !
      call fftw_f77(FFTW_PLAN_FWD,1,rspace(1 + nnx * i),1,1,ffts%workspace(dwx),1,1)
      !

  end do

  
      !
      ! edit davegp
      ! I do not destroy old plans anymore since I keep them
      ! for future use by FFTs on arrays of the same dimension
      !
      ! This was done before
      !
      !call fftw_f77_destroy_plan(FFTW_PLAN_FWD)
      !call fftw_f77_destroy_plan(FFTW_PLAN_BWD)
      !




  


  return

end subroutine fast_convolute
!
!     ======================== SERVICE ROUTINES ==========================
!
subroutine fftunfold(folded, gs, unfolded, ldz, nfn)

  include 'use.h'
  implicit none             ! implicit? Just say no!
  include 'interface.h'
  !
  !     INPUT:
  !     -----
  !
  type(parallel_gspace), intent(in) :: gs            ! the gspace for the data
  integer, intent(in) :: &
       ldz, &                                              ! leading dimension
       nfn                         ! number of functions to FFT simultaneously
  complex(dp), intent(in) :: &
       folded(gs%length, nfn)                             ! the data in gspace
  !
  !     OUTPUT:
  !     ------
  !
  complex(dp), intent(out) :: &
       unfolded(ldz, gs%lorder, nfn)
  !
  !     the unfold array must have the size Nz* number of rods held by
  !     this proc, i.e. nfn * gs%lorder * Nz
  !
  !
  !     DESCRIPTION:
  !     -----------
  !
  !     Transfers the rods from the packed gspace arrangement onto
  !     a linear grid of length gs%fftsize(3), such that a subsequent
  !     1d FFT can be performed
  !
  !     ----------------- local variables -----------------------
  !
  integer :: i, ii, nz, iord, istart, iend, iid, it, ifn

  unfolded = zzero                ! blank out target array

  !
  !     do not remove the it variable. It helps to circumvent a
  !     compiler bug on the SGI PowerChallenge Fortran 90 v.7.0
  !     
  !

  nz = gs%fftsize(3)

  do ifn = 1, nfn

     ii = 1
     it = 0
     do iord = 1, gs%lorder
        istart = gs%order(2, iord)
        iend   = gs%order(3, iord)
        it = it + iend
        if (istart <= iend) then         ! only one contiguous block to copy:
           iid    = ii + mod(iend - istart + nz, nz)
           unfolded(istart + 1:iend + 1, iord, ifn) = &
                folded(ii:iid, ifn)
        else                             ! two blocks to copy
           iid    = ii + nz - istart - 1
           unfolded(istart + 1:nz, iord, ifn) = &
                folded(ii:iid, ifn)
           ii     = iid + 1
           iid    = ii + iend
           unfolded(1:iend + 1, iord, ifn) = &
                folded(ii:iid, ifn)
        end if

        ii = iid + 1             ! new start  =    old end +1

     end do

  end do

  if (it < 0) write(9, *) it

end subroutine fftunfold


!     -----------------------------------------------------------------------

subroutine fftfold(folded, gs, unfolded, ldz, nfn)
  !
  include 'use.h'
  implicit none             ! implicit? Just say no!
  include 'interface.h'
  !
  !     1996 Bernd Pfrommer
  !
  !     INPUT:
  !     -----

  type(parallel_gspace), intent(in) :: gs           ! the gspace for the data
  integer, intent(in) :: ldz, nfn
  complex(dp), intent(in) :: &
       unfolded(ldz, gs%lorder, nfn)
  !
  !     OUTPUT:
  !     ------
  !
  complex(dp), intent(out) :: &
       folded(gs%length, nfn)                           ! the data in gspace

  !     DESCRIPTION:
  !     -----------
  !
  !     Collects the data from a 1d grid into the packed gspace
  !     arrangement.
  !
  !     Scales the result by   1/nx*ny*nz, except for the alpha and convex
  !
  !     the unfold array must have the size Nz* number of rods held by
  !     this proc, i.e. nfn * gs%lorder * Nz.
  !
  !     the fold array must have the length of the gspace
  !
  !
  !     ----------------- local variables -----------------------
  !
  integer :: i, ii, nz, iord, istart, iend, iid, it, ifn
  !
  !     do not remove the it variable. It helps to circumvent a
  !     compiler bug on the SGI PowerChallenge Fortran 90 v.7.0
  !
  real(dp) :: sc

  sc = done / real(gs%fftsize(1) * gs%fftsize(2) * gs%fftsize(3), dp)



  nz = gs%fftsize(3)

  do ifn = 1, nfn

     ii = 1
     it = 0
     do iord = 1, gs%lorder
        istart = gs%order(2, iord)
        iend   = gs%order(3, iord)
        it = it + iend
        if(istart <= iend) then           ! only one contiguous block to copy:
           iid    = ii + mod(iend - istart + nz, nz)
           folded(ii:iid, ifn) = unfolded(istart + 1:iend + 1, iord, ifn) * sc
        else                                    ! there are two blocks to copy
           iid    = ii + nz - istart - 1
           folded(ii:iid, ifn) = unfolded(istart + 1:nz, iord, ifn) * sc
           ii     = iid + 1
           iid    = ii + iend
           folded(ii:iid, ifn) = unfolded(1:iend + 1, iord, ifn) * sc
        end if
        ii = iid + 1          ! new start  =    old end +1
     end do
  end do

  if (it < 0) write(9, *) it

end subroutine fftfold

!
!     ============================================================

subroutine fft_convolute(gs, ffts, vpsi, vloc, psi, nwfn)
  !
  include 'use.h'
  implicit none             ! implicit? Just say no!
  include 'interface.h'
  include 'timing.h'
  !
  !     1996 Bernd Pfrommer
  !
  !     INPUT:
  !     -----
  !
  type(parallel_gspace), intent(in) :: gs            ! the gspace for the data
  type(fft_struc), intent(in) :: ffts                ! work arrays for the fft 
  integer, intent(in) :: nwfn            ! number of wavefunctions to convolve
  complex(dp), intent(in) :: &
       vloc(gs%r_size), &                   ! the local potential in realspace
       psi(gs%length * nwfn)                   ! the original wave function(s)
  !
  !     OUTPUT:
  !     ------
  !
  complex(dp), intent(out) :: &                            ! vpsi = vloc * psi
       vpsi(gs%length * nwfn)
  !
  !     DESCRIPTION:
  !     -----------
  !
  !     Applies the local potential vloc to a wave function psi 
  !     by 
  !            1) inverse FFT to realspace
  !            2) Multiply with vloc
  !            3) FFT to gspace
  !
  !
  !     ---------------------- local variables ----------------------
  !
  integer :: i, ifn
  real(dp) :: t0, ttot, t12i, t12f, tfc
  real(dp), external :: gimmetime



  
      !
      ! edit davegp
      ! I changed this parameter setting so that you can
      ! link to the correct header associated with the library
      ! that you link. Requires an include directory to be
      ! added to conf.mk
      !
      include 'fftw_f77.i'
      !
      ! INTEGER FFTW_FORWARD,FFTW_BACKWARD
      ! PARAMETER (FFTW_FORWARD=-1,FFTW_BACKWARD=1)
      ! INTEGER FFTW_ESTIMATE,FFTW_MEASURE
      ! PARAMETER (FFTW_ESTIMATE=0,FFTW_MEASURE=1)
      ! INTEGER FFTW_IN_PLACE,FFTW_USE_WISDOM
      ! PARAMETER (FFTW_IN_PLACE=8,FFTW_USE_WISDOM=16)

! for 64-bit machines (integer is not enough)
! it should also work on 32-bit machines

      ! the plan tags from create_plan
      INTEGER*8 FFTW_PLAN_FWD, FFTW_PLAN_BWD


  !
  !     now convolute with local potential
  !


  !
  !     this is the faster convolute technique
  !
  

  call fourier_transform_12(-1, ffts, gs, psi(1), ffts%rspacebuf(1), nwfn)

  
  

!------------------------------------------------------------------
! something is wrong here, gs%r_size is much greater than gs%length
! vpsi upper bound is over-reached.
! shouldn't this
!          vpsi(ifn * gs%r_size + 1)
! be 
!          vpsi(ifn * gs%length + 1)
! here?
! PZ
!------------------------------------------------------------------

  do ifn = 0, nwfn - 1

     call fast_convolute(ffts, gs, psi(ifn * gs%length + 1), &
!          vpsi(ifn * gs%r_size + 1), &
          vpsi(ifn * gs%length + 1), &
          ffts%rspacebuf(ifn * gs%r_size + 1), vloc(1))

  end do
  



  call fourier_transform_12(1, ffts, gs,vpsi(1), ffts%rspacebuf(1), nwfn)

  


  return

100 format('TIME FOR ',a30,f12.6)
110 format('FFT:',f8.4,'   FFTI:',f8.1,'  CONV:',f8.1,'  FFTF:',f8.1)
123 format('FFTI:',2f8.4,', UFOLD:',f8.4,', TRODS:',f8.4, &
       ', ZBUF:', f8.4,', ZRSPC:', f8.4,', TLT1:', f8.4, &
       ', TLT2:', f8.4,', TC1:', f8.4,', TC2:', f8.4)
124 format('FFTI:',2f8.1,', UFOLD:',f8.1,', TRODS:',f8.1, &
       ', ZBUF:', f8.1,', ZRSPC:', f8.1,', TLT1:', f8.1, &
       ', TLT2:', f8.1,', TC1:', f8.1,', TC2:', f8.1)
223 format('FFTF:',2f8.4,', UFOLD:',f8.4,', TRODS:',f8.4, &
       ', ZBUF:', f8.4,', ZRSPC:', f8.4,', TLT1:', f8.4, &
       ', TLT2:', f8.4,', TC1:', f8.4,', TC2:', f8.4)
224 format('FFTF:',2f8.1,', UFOLD:',f8.1,', TRODS:',f8.1, &
       ', ZBUF:', f8.1,', ZRSPC:', f8.1,', TLT1:', f8.1, &
       ', TLT2:', f8.1,', TC1:', f8.1,', TC2:', f8.1)

end subroutine fft_convolute
