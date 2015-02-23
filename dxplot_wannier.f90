!
subroutine dxplot_wannier(imode, ioutput, data, gs, ffts, crys, irk, is)
  !
  use all_to_all_module  
  include 'use.h'  
  implicit none                      ! implicit? Just say no!
  include 'interface.h'  
  include 'all_to_all.h'
  include 'mpif.h'
  !
  !     prepare a realspace  plot of a complex data array data.
  !     writes a UNK0000k.s file which can be imported by the wannier90 code
  !     for plotting wannier orbits.
  !     
  !     Revised from dxplot.f90                        Mar 14, 2008 
  !     
  !	
  !     INPUT
  !     -----
  !
  type(parallel_gspace), intent(in) :: gs ! gspace on which the data is defined
  type(fft_struc), intent(in) :: ffts  
  type(crystal), intent(in) :: crys  
!  integer, intent(in) ::  numout 
  integer, intent(in) :: ioutput, &    ! flag to determine format (dx/khoros)
       imode                           ! mode flag
                                       ! imode = 1 print real part of data
                                       !       = 2 print imaginary part of data
                                       !       = 3 print absolute squared
                                       !       = 0 print absolute of data
  integer, intent(in) :: irk, is
  complex(dp), intent(in) :: data(gs%length)            ! data to be visualized
  !
  !     ------------- local variables
  !
  integer :: iproc, i, j, k, ic, nbnd
  Integer :: ffts_total, info, myproc
  Integer :: Gdisp(0:gs%nproc-1)
  Integer :: Rcounts(0:gs%nproc-1)
  complex(dp), allocatable :: chd(:),chd1(:)            ! for realspace data
  character(len=11) :: wanname
  
  !
  !     -------------------------------------------------------------
  !
  call parallel_barrier()

  allocate(chd(ffts%r_size))
  chd = (0.0D0, 0.0D0)

  call fourier_transform(-1,ffts, gs, data, chd, 1)

  call parallel_barrier()
  ic = 0
  Call MPI_GATHER(ffts%r_size, 1, MPI_INTEGER, Rcounts, 1, MPI_INTEGER, &
                   0, MPI_COMM_WORLD, info)
  if (gs%myproc .eq. 0) then
     Gdisp(0) = 0
     do iproc = 1, gs%nproc-1
        Gdisp(iproc) = Gdisp(iproc - 1) +Rcounts(iproc -1)
     end do

     ffts_total = ffts%fftsize(1)*ffts%fftsize(2)*ffts%fftsize(3)
     allocate(chd1(ffts_total))
     chd1 = (0.0D0, 0.0D0)
  end if

  Call MPI_GATHERV(chd, ffts%r_size, MPI_DOUBLE_COMPLEX, chd1, Rcounts , &
                Gdisp ,MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, info)


  if(gs%myproc .eq. 0) then

    write(18) (chd1(i), i=1,ffts_total)
    deallocate(chd1)
  end if

  deallocate(chd)
  return    
  !------ wannier out put (UNK file) 
180 format(6f10.6)
200 format ('UNK',i5.5,'.',i1)
201 format (5i5)
280 format(F10.4,4X,F10.4)
end subroutine dxplot_wannier
