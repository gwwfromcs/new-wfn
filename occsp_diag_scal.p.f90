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

  subroutine occsp_diag_scal(myproc,h,nmat,nproc,eval,work,lwork,rwork)
  use constants
  use all_to_all_module 
  include 'use.h'
  IMPLICIT NONE
  include 'interface.h'  
  include 'all_to_all.h'  
  include 'flibcalls.ph'        
!
!      distributed solver (block cyclic blacs layout):
!
!           nbl = blocksize
!           nprow = processor grid row
!           npcol = processor grid column

  include "mpif.h"

  !
  !     INPUT:
  !     -----
  !
  integer, intent(in) :: myproc,nproc,nmat,lwork
  !
  !     INPUT/OUTPUT:
  !     ------------
  !
  complex(dp) h(nmat*nmat)   ! matrix to diagonalize 
  !
  !     OUTPUT:
  !     -------
  !
  real(dp) eval(nmat)
  !
  !     WORK:
  !     -------
  !
  complex(dp) work(lwork)
  real(dp) rwork(max(1, 3*nmat-2))
  !
  !     --------------------- local variables ----------------------
  !
  integer nb,info
  complex(dp), allocatable :: h_n(:),cbuf(:)
! scalapack and blacks arrays
  integer nbl,nprow,npcol,myprow,mypcol,icntxt
  integer nbc,nbr,ngc,nbc_max,nbr_max,num_max,ir_low,ilen,numbl_max,length
  integer idiff,ibc,ngr,ibr,ic
  integer sendprow,sendpcol,i,nbc_t,nbr_t
  integer color,comm_diag,ierror
  logical, SAVE :: print_occ = .false.
  integer iproc


 !
 !  choose sclapack layout. Block cyclic. Block size as
 !  close to 32 as possible for large matrices
 !  processor layout as close as possible to square
 !
 if (nmat .gt. 150) then 

!! ============
!! davegp debug
!
!  write(100+myproc,*) 'process: ', myproc, 'about to call layout_occ_scalapack'
!  write(100+myproc,*) 'before: ', nmat, nbl, nproc, nprow, npcol

  call layout_occ_scalapack(nmat, nbl, nproc, nprow, npcol)

!  write(100+myproc,*) ' after: ', nmat, nbl, nproc, nprow, npcol
!!
!  !write(100+myproc,*) 'process: ', myproc, 'instead to call layout_scalapack'
!  !write(100+myproc,*) 'before: ', nmat, nbl, nproc, nprow, npcol
!  !call layout_scalapack(nmat, nbl, nproc, nprow, npcol)
!  !write(100+myproc,*) ' after: ', nmat, nbl, nproc, nprow, npcol
!
!  call myflush(100+myproc)

 ! figure where myproc is in proc grid
  call blacs_get(-1,0,icntxt)
  call blacs_gridinit(icntxt,'c',nprow, npcol)
  call blacs_gridinfo(icntxt, nprow, npcol, myprow, mypcol)

  color=0
  if (myproc .le. nprow*npcol-1)  color=1
  call mpi_comm_split(MPI_COMM_WORLD,color,myproc, comm_diag,ierror)
  if (color .eq. 0) go to 20

! ============
! davegp debug

  !if (.not. print_occ) then
  write(9,*) 'DOING OCC. SPACE DIAG in parallel'
  write(9,400) nmat,nprow, npcol,nbl
400  format (' OCC. SPACE SIZE ',i5, &
         /' PROCESSOR GRID ',i3,' x ',i3,' BLOCKSIZE ',i3)
   print_occ =.true.
  !end if

!
! figure out number of blocks per processor in the column/row
! and array sizes and maximums as well for allocations
!
  nbc = nmat/(nbl*npcol)
  if(mod(nmat,(nbl*npcol)).gt.mypcol*nbl) nbc=nbc+1  ! for partial matrix
  nbr = nmat/(nbl*nprow)
  if(mod(nmat,(nbl*nprow)).gt.myprow*nbl) nbr=nbr+1
  nbc_max = nbc+1   !zero contribution
  nbr_max = nbr+1
  numbl_max = nbc_max*nbr_max
  num_max=numbl_max

  call all_max_all(num_max,comm_diag)

  allocate (h_n( num_max*num_max*nbl*nbl))
  allocate (cbuf( num_max*num_max*nbl*nbl)) 

  idiff=0; ngc=0;
  ngr=0;
  do ibc=0, nbc-1           ! loop over column blocks
    do ic=(ibc*npcol+mypcol)*nbl+1, &
                     min((ibc*npcol+mypcol)*nbl+nbl,nmat) ! loop over cols
      ngr = 0
      ngc=ngc+1
      do ibr =0, nbr-1    ! loop over row blocks
        ir_low = (ibr*nprow+myprow)*nbl+1
        ilen = min((ibr*nprow+myprow)*nbl+nbl,nmat)-ir_low+1
        ngr = ngr + ilen

        if (ilen .gt. 0) &
         call mzcopy(ilen,h(ir_low+(ic-1)*nmat),1, h_n(idiff+1),1)
        idiff = idiff + ilen
      end do
    end do
  end do





     call zdiag_occ_scalapack(icntxt, nmat, nbl, nprow, npcol, &
          h_n(1), ngc, ngr, eval(1),comm_diag)

  do iproc=0,nprow*npcol-1
   
    if (myproc.eq.iproc) then
      cbuf=h_n
      length=ngc*ngr
      nbc_t=nbc
      nbr_t=nbr
      sendprow=myprow
      sendpcol=mypcol
    end if
       
      call my_broadcast_small(length,iproc,comm_diag)
      call my_broadcast_small(nbc_t,iproc,comm_diag)
      call my_broadcast_small(nbr_t,iproc,comm_diag)
      call my_broadcast_small(sendprow,iproc,comm_diag)
      call my_broadcast_small(sendpcol,iproc,comm_diag)
      call my_broadcast_small(cbuf,length,iproc,comm_diag)

    idiff=0;
    do ibc=0, nbc_t-1           ! loop over column blocks
      do ic=(ibc*npcol+sendpcol)*nbl+1, &
                     min((ibc*npcol+sendpcol)*nbl+nbl,nmat) ! loop over cols
        do ibr =0, nbr_t-1    ! loop over row blocks
          ir_low = (ibr*nprow+sendprow)*nbl+1
          ilen = min((ibr*nprow+sendprow)*nbl+nbl,nmat)-ir_low+1

        if (ilen .gt. 0) &
          call mzcopy(ilen,cbuf(idiff+1),1,h(ir_low+(ic-1)*nmat),1)

          idiff = idiff + ilen
        end do
      end do
    end do

  end do

  deallocate (h_n)
  deallocate (cbuf)
  call BLACS_GRIDEXIT(icntxt)

   20  continue

  call parallel_barrier()

   call mpi_comm_free(comm_diag,ierror)

   call my_broadcast(h,nmat*nmat,0)
   call my_broadcast(eval,nmat,0)

  else

      if(myproc.eq.0) then
        call mzheev('V','L',nmat,h(1),nmat, eval(1),work(1),lwork, &
             rwork(1),info) 
        if(info.ne.0) then
          if(myproc.eq.0)&
            write(9,*) 'diagonalization failed with info=',info
            call mystop
        endif
      endif
      call myflush(9)
      call my_broadcast(h,nmat*nmat,0)
      call my_broadcast(eval,nmat,0)

  end if



  return
  end subroutine occsp_diag_scal
