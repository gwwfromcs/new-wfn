!
subroutine dxplot(imode, ioutput, data, gs, ffts, crys, outfilename,numout)
  !
  use all_to_all_module  
  include 'use.h'  
  implicit none                      ! implicit? Just say no!
  include 'interface.h'  
  include 'all_to_all.h'  
  !
  !     prepare a realspace  plot of a complex data array data.
  !     writes a .dx file which can be imported by the IBM data explorer,
  !     a visualization software package. If ioutput has the khoros flag
  !     set, it produces a file suitable for the khoros graphics
  !     package instead of a dx file.
  !
  !     1996, 1997 Bernd Pfrommer
  !
  !     INPUT
  !     -----
  !
  type(parallel_gspace), intent(in) :: gs ! gspace on which the data is defined
  type(fft_struc), intent(in) :: ffts  
  type(crystal), intent(in) :: crys  
  integer, intent(in) ::  numout 
  integer, intent(in) :: ioutput, &    ! flag to determine format (dx/khoros)
       imode                           ! mode flag
                                       ! imode = 1 print real part of data
                                       !       = 2 print imaginary part of data
                                       !       = 3 print absolute squared
                                       !       = 0 print absolute of data
  character(len=numout) :: &
       outfilename             ! name of output file. Suffix .dx will be added
  complex(dp), intent(in) :: data(gs%length)            ! data to be visualized
  !
  !     ------------- local variables
  !
  integer :: iproc, i, j, k, ic  
  logical :: ikhoros  
  real(dp) :: delta(3, 3), sc1, sc2, sc3, ri(3), x(3)  
  complex(dp), allocatable :: chd(:)                      ! for realspace data
  !
  !     -------------------------------------------------------------
  !
  ikhoros = (iand(ioutput, 32768) == 32768)                  ! if khoros or dx

  call dx_latticevec(crys, gs%myproc)  

  call parallel_barrier()  

  allocate(chd(ffts%r_size))  
  !
  !     inverse fourier transform the data to real space
  !
  call fourier_transform(-1, ffts, gs, data, chd, 1)  
  !
  !     --------- open files and write the header
  !
  if (gs%myproc == 0) then  
     if (ikhoros) then  
        open(unit = 18, file = 'vis/'//outfilename//'.avs', &
             status = 'unknown', form = 'formatted')
        write(18, 105) gs%fftsize(1), gs%fftsize(2), gs%fftsize(3) ! dx header
     else  
        open(unit = 18, file = 'vis/'//outfilename//'.dx', &
             status = 'unknown', form = 'formatted')
        write(18, 5)                          ! write the header for the array
        write(18, 10) gs%fftsize(1) * gs%fftsize(2) * gs%fftsize(3)
     end if
     close(18)  
  end if

  call parallel_barrier()  
  !
  !     --------- now write the data
  !
  !     Due to the specific data layout, the a1 axis runs the fastest,
  !     then the a2 axis, and finally the a3 axis runs slowest.
  !
  sc1 = done / real(gs%fftsize(1), dp)  
  sc2 = done / real(gs%fftsize(2), dp)  
  sc3 = done / real(gs%fftsize(3), dp)  

  ic = 0  

  do iproc = gs%nproc - 1, 0, -1  
     if (iproc == gs%myproc) then  
        if (ikhoros) then  
           open(unit = 18, file = 'vis/'//outfilename//'.avs', &
                position = 'append', status = 'unknown', form = 'formatted')
           do i = 1, ffts%r_size  
              ri(1) = sc1 * mod(ic, gs%fftsize(1))  
              ri(2) = sc2 * mod(ic / (gs%fftsize(1)), gs%fftsize(2))  
              ri(3) = sc3 * (ic / (gs%fftsize(1) * gs%fftsize(2)))  
              x = matmul(crys%avec, ri)  
              ic = ic + 1  
              if (imode == 1) then  
                 write(18, 180) x, real(chd(ffts%r_size - i + 1), dp)  
              else if (imode == 2) then  
                 write(18, 180) x, aimag(chd(ffts%r_size - i + 1))  
              else if (imode == 3) then  
                 write(18, 180) x, abs(chd(ffts%r_size - i + 1)) * &
                      abs(chd(ffts%r_size - i + 1))
              else  
                 write(18, 180) x, abs(chd(ffts%r_size - i + 1))  
              end if
           end do
        else  
           open(unit = 18, file = 'vis/'//outfilename//'.dx', &
                position = 'append', status = 'unknown', form = 'formatted')
           if (imode == 1) then  
              write(18, 80) (real(chd(i), dp), i = ffts%r_size, 1, -1)
           else if (imode == 2) then  
              write(18, 80) (aimag(chd(i)), i = ffts%r_size, 1, -1)
           else if (imode == 3) then  
              write(18, 80) (abs(chd(i)) * abs(chd(i)), i = ffts%r_size, 1, -1)
           else  
              write(18, 80) (abs(chd(i)), i = ffts%r_size, 1, -1)
           end if
        end if
        call myflush(18)  
        close(18)  
     end if
     call my_broadcast(ic, iproc)  
  end do

  if (.not.ikhoros) then                             !     now the dx trailer
     if (gs%myproc == 0) then  
        open (unit = 18, file = 'vis/'//outfilename//'.dx', &
             position = 'append', status = 'unknown', form = 'formatted')
        !        grid connections object
        write(18, 20) gs%fftsize(3), gs%fftsize(2), gs%fftsize (1)
        !        calculate deltas
        delta(1:3, 1) = crys%avec(1:3, 3) / gs%fftsize(3)  
        delta(1:3, 2) = crys%avec(1:3, 2) / gs%fftsize(2)  
        delta(1:3, 3) = crys%avec(1:3, 1) / gs%fftsize(1)  
        !     write header for gridpositions
        write(18, 30) gs%fftsize(3), gs%fftsize(2), gs%fftsize(1)
        !     compute and write origin
        write(18, 40) dzero, dzero, dzero
        !     write deltas
        write(18, 50) ((delta(i, j), i = 1, 3), j = 1, 3)  
        !     write header for field object
        write(18, 60)  
        close(18)  
        ! i am proc 0
     end if
     
  end if
  if (ikhoros) then                                              ! khoros
     write(9, 199) outfilename  
  else  
     write(9, 99) outfilename  
  end if

  deallocate(chd)  

  return  

  !     ----------------- dx output -----------------------------
5 format(/'# this is the object defining the data values')  
10 format('object 1 class array type float rank 0 items ',i8, &
       &     ' data follows')

20 format(/'# this is the object defining the grid connections', &
       &     /'object 2 class gridconnections counts ',3i4)

30 format(//'# this is the object defining the grid ', &
       &     /'object 3 class gridpositions counts ',3i4 )
40 format('origin', 3f20.10)  

50 format('delta ', 3f20.10)  


60 format(/'# this is the collective object, one for each grid', &
       &     /'object 4 class field', &
       &     /'component "positions"   value 3', &
       &     /'component "connections" value 2', &
       &     /'component "data"        value 1')
80 format(5g15.7)  


99 format(//' GRAPHICAL OUTPUT:', &
       &     /' ----------------', &
       &     //' INPUT DATA FOR IBM DATA EXPLORER WRITTEN TO FILE', &
       &     /'    - ',a13,'.dx')
  !     ----------------- khoros output -----------------------------
105 format('3',/'3',3(/i6),/'1')  

180 format(3f10.6,5x, g15.7)  
199 format(//' GRAPHICAL OUTPUT:', &
       &     /' ----------------', &
       &     //' INPUT DATA FOR AVS/KHOROS WRITTEN TO FILE', &
       &     /'    - ',a13,'.avs')

end subroutine dxplot
