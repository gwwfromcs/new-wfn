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

!
module all_to_all_module
  !
  !     first some global variables. 
  !
  use constants

  integer :: ata_nproc         ! number of processors

  
contains

  subroutine allocate_shmem_work(nproc,len)

    implicit none
    integer, intent(in) :: nproc, len
    integer :: ierr

   
  end subroutine allocate_shmem_work

  subroutine deallocate_shmem_work()

    implicit none
    integer :: ierr


  end subroutine deallocate_shmem_work
  !
  !     allocs and initializes the global reduction work arrays
  !     required  for shmem
  !
  subroutine fast_all_sum_all_alloc_dc(len)

    implicit none
    integer, intent(in) :: len        ! number of elements to be reduced
    integer :: ierr


  end subroutine fast_all_sum_all_alloc_dc



  !
  !     ================ summation subroutines ==============================
  !
  !     this subroutine uses the fast shmem routines,
  !     but needs a call to the allocation routine before
  !     it is called the first time
  !
  subroutine fast_all_sum_all_complex(data, len)

    implicit none
include 'mpif.h'
    integer :: ierr,i,mblock,istart,iend,mblock2,nblock
    integer, intent(in) :: len
    complex(dp), intent(inout) :: data(*)
!    complex(dp), intent(inout) :: data(len)


    !
    !     just the regular thing for 
    !
      
    complex(dp), allocatable :: mpibuf(:)
    if (len <= 0) return
  

    if(len.le.(4*1024*1024)) then
    allocate(mpibuf(len))
    call mpi_allreduce(data(1), mpibuf(1), len, &
         MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
    data(1:len) = mpibuf(1:len)
    else
    write(9,*) "Dangeous! in fast_all_sum_all_complex, len",len



    mblock=1024*1024
    nblock=(len+mblock-1)/mblock
    allocate(mpibuf(mblock))
    
    do i=1,nblock
    istart=1+(i-1)*mblock 
    iend=min(i*mblock,len)
    mblock2=iend-istart+1
!    write(9,*) "i,istart,iend",i,istart,iend
!    call mpi_allreduce(data(1), mpibuf(1), len, &
!         MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
    call mpi_allreduce(data(istart), mpibuf(1), mblock2, &
         MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
!    data(1:len) = mpibuf(1:len)
    data(istart:iend) = mpibuf(1:mblock2)
    end do
    end if

    deallocate(mpibuf)



  end subroutine fast_all_sum_all_complex

  subroutine all_sum_all_complex(data, len)

    implicit none
include 'mpif.h'
    integer :: ierr
    integer, intent(in) :: len
    complex(dp), intent(inout) :: data(*)
    complex(dp), allocatable :: mpibuf(:)

      
    if(len.ge.(8*1024*1024)) write(9,*) "Dangeous! in all_sum_all_complex, len",len
    if (len <= 0) return
    allocate(mpibuf(len))
    call mpi_allreduce(data(1), mpibuf(1), len, &
         MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
    data(1:len) = mpibuf(1:len)
    deallocate(mpibuf)


  end subroutine all_sum_all_complex

  subroutine all_sum_all_complex1(data)

    implicit none
include 'mpif.h'
    integer :: ierr
    complex(dp), intent(inout) :: data
    complex(dp) :: mpibuf

      
    call mpi_allreduce(data,mpibuf,1, &
         MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
    data = mpibuf


  end subroutine all_sum_all_complex1

  subroutine all_sum_all_complex2(data, len)

    implicit none
include 'mpif.h'
    integer :: ierr
    integer, intent(in) :: len
    complex(dp), intent(inout) :: data(:,:)
    complex(dp), allocatable :: mpibuf(:,:)

      
    if(len.ge.(8*1024*1024)) write(9,*) "Dangeous! in all_sum_all_complex2, len",len
    if (len <= 0) return
    allocate(mpibuf(size(data,1),size(data,2)))
    call mpi_allreduce(data(1,1), mpibuf(1,1), len, &
         MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
    data = mpibuf
    deallocate(mpibuf)


  end subroutine all_sum_all_complex2

  subroutine all_sum_all_complex3(data, len)

    implicit none
include 'mpif.h'
    integer :: ierr
    integer, intent(in) :: len
    complex(dp), intent(inout) :: data(:,:,:)
    complex(dp), allocatable :: mpibuf(:,:,:)

      
    if(len.ge.(8*1024*1024)) write(9,*) "Dangeous! in all_sum_all_complex3, len",len
    if (len <= 0) return
    allocate(mpibuf(size(data,1),size(data,2),size(data,3)))
    call mpi_allreduce(data(1,1,1), mpibuf(1,1,1), len, &
         MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
    data = mpibuf
    deallocate(mpibuf)


  end subroutine all_sum_all_complex3

  subroutine all_sum_all_complex4(data, len)

    implicit none
include 'mpif.h'
    integer :: ierr
    integer, intent(in) :: len
    complex(dp), intent(inout) :: data(:,:,:,:)
    complex(dp), allocatable :: mpibuf(:,:,:,:)

      
    if(len.ge.(8*1024*1024)) write(9,*) "Dangeous! in all_sum_all_complex4, len",len
    if (len <= 0) return
    allocate(mpibuf(size(data,1),size(data,2),size(data,3), &
         size(data,4)))
    call mpi_allreduce(data(1,1,1,1), mpibuf(1,1,1,1), len, &
         MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
    data = mpibuf
    deallocate(mpibuf)


  end subroutine all_sum_all_complex4

  subroutine all_sum_all_complex_sub(comm, data, len)

    implicit none
include 'mpif.h'
    integer :: ierr
    integer, intent(in) :: len, comm
    complex(dp), intent(inout) :: data(*)
    complex(dp), allocatable :: mpibuf(:)

      
    if(len.ge.(8*1024*1024)) write(9,*) "Dangeous! in all_sum_all_complex_sub, len",len
    if(len <= 0) return
    allocate(mpibuf(len))
    call mpi_allreduce(data(1), mpibuf(1), len, &
         MPI_DOUBLE_COMPLEX, MPI_SUM, comm, ierr)
    data(1:len) = mpibuf(1:len)
    deallocate(mpibuf)


  end subroutine all_sum_all_complex_sub

  subroutine all_sum_all_complex1_sub(comm, data)

    implicit none
include 'mpif.h'
    integer :: ierr
    integer, intent(in) :: comm
    complex(dp), intent(inout) :: data
    complex(dp) :: mpibuf

      
    call mpi_allreduce(data, mpibuf, 1, &
         MPI_DOUBLE_COMPLEX, MPI_SUM, comm, ierr)
    data = mpibuf


  end subroutine all_sum_all_complex1_sub

  subroutine all_sum_all_double(data, len)

    implicit none
include 'mpif.h'
    integer :: ierr
    integer, intent(in) :: len
    real(dp), intent(inout) :: data(len)
    real(dp), allocatable :: mpibuf(:)

      
    if(len.ge.(8*1024*1024)) write(9,*) "Dangeous! in all_sum_all_double, len",len
    if (len <= 0) return
    allocate(mpibuf(len))
    call mpi_allreduce(data(1), mpibuf(1), len, &
         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    data = mpibuf
    deallocate(mpibuf)


  end subroutine all_sum_all_double

  subroutine all_sum_all_double1(data)

    implicit none
include 'mpif.h'
    integer :: ierr
    real(dp), intent(inout) :: data
    real(dp) :: databuf

      
    call mpi_allreduce(data, databuf, 1, &
         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    data = databuf


  end subroutine all_sum_all_double1

  subroutine all_sum_all_double2(data, len)

    implicit none
include 'mpif.h'
    integer :: ierr
    integer, intent(in) :: len
    real(dp), intent(inout) :: data(:,:)
    real(dp), allocatable :: mpibuf(:,:)

      
    if(len.ge.(8*1024*1024)) write(9,*) "Dangeous! in all_sum_all_double2, len",len
    if (len <= 0) return
    allocate(mpibuf(size(data,1),size(data,2)))
    call mpi_allreduce(data(1,1), mpibuf(1,1), len, &
         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    data = mpibuf
    deallocate(mpibuf)


  end subroutine all_sum_all_double2

  subroutine all_sum_all_double3(data, len)

    implicit none
include 'mpif.h'
    integer :: ierr
    integer, intent(in) :: len
    real(dp), intent(inout) :: data(:,:,:)
    real(dp), allocatable :: mpibuf(:,:,:)

      
    if(len.ge.(8*1024*1024)) write(9,*) "Dangeous! in all_sum_all_double3, len",len
    if (len <= 0) return
    allocate(mpibuf(size(data,1),size(data,2),size(data,3)))
    call mpi_allreduce(data(1,1,1), mpibuf(1,1,1), len, &
         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    data = mpibuf
    deallocate(mpibuf)


  end subroutine all_sum_all_double3

  subroutine all_sum_all_integer(data, len)

    implicit none
include 'mpif.h'
    integer :: ierr
    integer, intent(in) :: len
    integer, intent(inout) :: data(*)
    integer, allocatable :: mpibuf(:)

      
    if(len.ge.(8*1024*1024)) write(9,*) "Dangeous! in all_sum_all_integer, len",len
    if (len <= 0) return
    allocate(mpibuf(len))
    call mpi_allreduce(data, mpibuf, len, &
         MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    data(1:len) = mpibuf(1:len)
    deallocate(mpibuf)


  end subroutine all_sum_all_integer

  subroutine all_sum_all_integer1(data)

    implicit none
include 'mpif.h'
    integer :: ierr
    integer, intent(inout) :: data
    integer :: databuf

      
    call mpi_allreduce(data, databuf, 1, &
         MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    data = databuf


  end subroutine all_sum_all_integer1
  !
  !     ================= maximum subroutines ===========================
  !
  subroutine all_max_all_double1(data)

    implicit none
include 'mpif.h'
    integer :: ierr
    real(dp), intent(inout) :: data
    real(dp) :: mpibuf

      
    call mpi_allreduce(data, mpibuf, 1, &
         MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    data = mpibuf


  end subroutine all_max_all_double1

  subroutine all_min_all_double1(data)

    implicit none
include 'mpif.h'
    integer :: ierr
    real(dp), intent(inout) :: data
    real(dp) :: mpibuf

      
    call mpi_allreduce(data, mpibuf, 1, &
         MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    data = mpibuf


  end subroutine all_min_all_double1

  subroutine all_max_all_int1(data)

    implicit none
include 'mpif.h'
    integer :: ierr
    integer, intent(inout) :: data
    integer :: mpibuf

      
    call mpi_allreduce(data, mpibuf, 1, &
         MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    data = mpibuf


  end subroutine all_max_all_int1

  subroutine all_min_all_int1(data)

    implicit none
include 'mpif.h'
    integer :: ierr
    integer, intent(inout) :: data
    integer :: mpibuf

      
    call mpi_allreduce(data, mpibuf, 1, &
         MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
    data = mpibuf


  end subroutine all_min_all_int1

  subroutine all_max_all_int2(data,comm)

    implicit none
include 'mpif.h'
    integer :: ierr
    integer, intent(inout) :: data
    integer :: mpibuf,comm

      
    call mpi_allreduce(data, mpibuf, 1, &
         MPI_INTEGER, MPI_MAX, comm, ierr)
    data = mpibuf


  end subroutine all_max_all_int2

  subroutine all_min_all_int2(data,comm)

    implicit none
include 'mpif.h'
    integer :: ierr
    integer, intent(inout) :: data
    integer :: mpibuf,comm

      
    call mpi_allreduce(data, mpibuf, 1, &
         MPI_INTEGER, MPI_MIN, comm, ierr)
    data = mpibuf


  end subroutine all_min_all_int2

  !
  !     ================= broadcast subroutines ===========================
  !
  subroutine my_broadcast_int1(data, root)

    implicit none
include 'mpif.h'
    integer :: ierr
    integer :: data
    integer :: root

      
    call mpi_bcast(data, 1, MPI_INTEGER, &
         root, MPI_COMM_WORLD, ierr)


  end subroutine my_broadcast_int1

  subroutine my_broadcast_int1_smcomm(data, root,comm)

    implicit none
include 'mpif.h'
    integer :: ierr
    integer :: data
    integer :: root,comm

      
    call mpi_bcast(data, 1, MPI_INTEGER, &
         root, comm, ierr)


  end subroutine my_broadcast_int1_smcomm

  subroutine my_broadcast_int(data, count, root)

    implicit none
include 'mpif.h'
    integer :: ierr
    integer, intent(in) :: count
    integer, intent(inout) :: data
    integer, intent(in) :: root

      
    if(count.ge.(8*1024*1024)) write(9,*) "Dangeous! in broadcast_int, count",count
    call mpi_bcast(data, count, MPI_INTEGER, &
         root, MPI_COMM_WORLD, ierr)


  end subroutine my_broadcast_int

  subroutine my_broadcast_dcmplx(data, count, root)

    implicit none
include 'mpif.h'
    integer :: ierr,mblock,mblock2,nblock,i,istart,iend
    integer, intent(in) :: count
    complex(dp), intent(inout) :: data(count)
    integer, intent(in) :: root

      


    if(count.le.(4*1024*1024)) then
    call mpi_bcast(data(1), count, MPI_DOUBLE_COMPLEX, &
         root, MPI_COMM_WORLD, ierr)

    else
    if(count.ge.(4*1024*1024)) write(9,*) "Dangeous! in my_broadcast_dcmplx, count",count

    mblock=1024*1024
    nblock=(count+mblock-1)/mblock

    do i=1,nblock
    istart=1+(i-1)*mblock
    iend=min(i*mblock,count)
    mblock2=iend-istart+1

!    write(9,*) "i,istart,iend",i,istart,iend
    call mpi_bcast(data(istart), mblock2, MPI_DOUBLE_COMPLEX, &
         root, MPI_COMM_WORLD, ierr)

    end do
    end if

!    call mpi_bcast(data(1), count, MPI_DOUBLE_COMPLEX, &
!         root, MPI_COMM_WORLD, ierr)


  end subroutine my_broadcast_dcmplx

  subroutine my_broadcast_dcmplx_small(data, count, root,comm)
 
    implicit none
include 'mpif.h'
    integer :: ierr
    integer, intent(in) :: count
    complex(dp), intent(inout) :: data(count)
    integer, intent(in) :: root,comm

      
    if(count.ge.(8*1024*1024)) write(9,*) "Dangeous! in my_broadcast_dcmplx_small, count",count
    call mpi_bcast(data(1), count, MPI_DOUBLE_COMPLEX, &
         root, comm, ierr)


  end subroutine my_broadcast_dcmplx_small


  subroutine my_broadcast_dcmplx2(data, count, dim, root)

    implicit none
include 'mpif.h'
    integer :: ierr
    integer, intent(in) :: count
    integer, intent(in) ::  dim
    complex(dp), intent(inout) :: data(count,dim)
    integer, intent(in) :: root

      
    if((count*dim).ge.(8*1024*1024)) write(9,*) "Dangeous! in my_broadcast_dcmplx2, count*dim",count*dim
    call mpi_bcast(data(1,1), count*dim, MPI_DOUBLE_COMPLEX, &
         root, MPI_COMM_WORLD, ierr)


  end subroutine my_broadcast_dcmplx2

  
  subroutine my_broadcast_dcmplx3(data, count,root)

    implicit none
include 'mpif.h'
    integer :: ierr
    integer, intent(in) :: count
    complex(dp), intent(inout) :: data
    integer, intent(in) :: root

      
    if(count.ge.(8*1024*1024)) write(9,*) "Dangeous! in my_broadcast_dcmplx3, count",count
    call mpi_bcast(data, count, MPI_DOUBLE_COMPLEX, &
         root, MPI_COMM_WORLD, ierr)


  end subroutine my_broadcast_dcmplx3

  subroutine my_broadcast_double(data, count, root)

    implicit none
include 'mpif.h'
    integer :: ierr
    integer, intent(in) :: count
    real(dp), intent(inout) :: data(count)
    integer, intent(in) :: root
      
    if(count.ge.(8*1024*1024)) write(9,*) "Dangeous! in my_broadcast_double, count",count
    call mpi_bcast(data(1), count, MPI_DOUBLE_PRECISION, &
         root, MPI_COMM_WORLD, ierr)


  end subroutine my_broadcast_double

  subroutine my_broadcast_double_smcomm(data, count, root,comm)

    implicit none
include 'mpif.h'
    integer :: ierr
    integer, intent(in) :: count
    real(dp), intent(inout) :: data(count)
    integer, intent(in) :: root,comm
      
    if(count.ge.(8*1024*1024)) write(9,*) "Dangeous! in my_broadcast_double_smcomm, count",count
    call mpi_bcast(data(1), count, MPI_DOUBLE_PRECISION, &
         root, comm, ierr)


  end subroutine my_broadcast_double_smcomm

  subroutine my_broadcast_double1(data, root)

    implicit none
include 'mpif.h'
    integer :: ierr
    real(dp), intent(inout) :: data
    integer, intent(in) :: root
      
    call mpi_bcast(data, 1, MPI_DOUBLE_PRECISION, &
         root, MPI_COMM_WORLD, ierr)


  end subroutine my_broadcast_double1

  subroutine my_broadcast_double3(data, count, root)

    implicit none
include 'mpif.h'
    integer :: ierr
    integer, intent(in) :: count
    real(dp), intent(inout) :: data
    integer, intent(in) :: root
      
    if(count.ge.(8*1024*1024)) write(9,*) "Dangeous! in my_broadcast_double3, count",count
    call mpi_bcast(data, count, MPI_DOUBLE_PRECISION, &
         root, MPI_COMM_WORLD, ierr)


  end subroutine my_broadcast_double3

  subroutine my_broadcast_char(data, count, root)

    implicit none
include 'mpif.h'
    integer :: ierr
    integer, intent(in) :: count
    character(len=count), intent(inout) :: data
    integer, intent(in) :: root
      
    call mpi_bcast(data, count, MPI_CHARACTER, &
         root, MPI_COMM_WORLD, ierr)


  end subroutine my_broadcast_char
  !
  !     ================= dot product subroutines =======================
  !
  complex(dp) function parallel_zdotc(n, x, incx, y, incy) 

    implicit none
    include 'flibcalls.ph'
    integer, intent(in) :: n, incx, incy
    complex(dp) :: x(*), y(*)

    parallel_zdotc = mzdotc(n, x(1), incx, y(1), incy)
    call all_sum_all_complex1(parallel_zdotc)
    return

  end function parallel_zdotc

  complex(dp) function parallel_zdotc2(n, x, incx, y, incy) 

    implicit none
    include 'flibcalls.ph'
    integer n, incx, incy
    complex(dp) :: x, y

    parallel_zdotc2 = mzdotc(n, x, incx, y, incy)
    call all_sum_all_complex1(parallel_zdotc2)
    return

  end function parallel_zdotc2
  !
  !     ================ parallel barrier routine =======================
  !
  subroutine parallel_barrier()

    implicit none
include 'mpif.h'
    integer :: info

      
    call mpi_barrier(MPI_COMM_WORLD, info)


  end subroutine parallel_barrier
  !
  !
  !     ================ subgroup communication routines =======================
  !
  !     creates communicators and corresponding groups for a list of ranks
  !
  subroutine all_split_all(key, rank, subcomm)

    implicit none
include 'mpif.h'
    integer, intent(in) :: rank, key
    integer, intent(out) :: subcomm        ! output
      
    integer :: info

    call mpi_comm_split(MPI_COMM_WORLD, key, rank, subcomm, info)


  end subroutine all_split_all

end module all_to_all_module
