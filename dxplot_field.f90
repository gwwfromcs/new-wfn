!     @process extchk
!
subroutine dxplot_field(imode, ioutput, data, gs, ffts, crys, &
     outfilename)

  use all_to_all_module  
  include 'use.h'  
  implicit none                  ! implicit? Just say no!
  include 'interface.h'
  include 'all_to_all.h'  
  !
  !     1996 Bernd Pfrommer
  !
  !     INPUT
  !     -----
  !
  type(parallel_gspace), intent(in) :: gs ! gspace on which the data is defined
  type(fft_struc), intent(in) :: ffts  
  type(crystal), intent(in) :: crys  
  integer, intent(in) :: ioutput, &                  ! output flags (khoros/dx)
       imode                           ! mode flag
                                       ! imode = 1 print real part of data
                                       !       = 2 print imaginary part of data
                                       !       = 0 print absolute of data
  character(len=7), intent(in) :: &
       outfilename              ! name of output file. Suffix .dx will be added
  complex(dp), intent(in) :: &
       data(gs%length, 3)                  ! vector field data to be visualized
  !
  !     DESCRIPTION:
  !     -----------
  !
  !
  !     prepare a realspace  plot of a complex vector field data array dat
  !     writes a .dx file which can be imported by the IBM data explorer,
  !     a visualization software package.
  !
  !     ------------- local variables --------------------------------
  !
  integer :: n1, n2, n3, iproc, i, j, k, l, ic  
  logical :: ikhoros  
  real(dp) :: delta(3, 3), sc1, sc2, sc3, ri(3), x(3)  
  complex(dp), allocatable :: chd(:,:)                 ! for realspace data
  !
  !     ------------------------------------------------------------------
  !
  ikhoros = (iand(ioutput, 32768) == 32768)               ! if khoros or dx
  allocate(chd(ffts%r_size, 3))  
  !
  !     write the lattice vector file first
  !
  call dx_latticevec(crys, gs%myproc)  
  call parallel_barrier()  
  !
  !     get dimensions of my block
  !
  n1 = gs%fftsize(1)  
  n2 = gs%fftsize(2)  
  n3 = gs%fftsize(3)  
  !
  !     inverse fourier transform the data to real space
  !
  do i = 1, 3  
     call fourier_transform(-1, ffts, gs, data(1, i), chd(1, i), 1)
  enddo
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
  !     now write the data
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
                 write(18, 180) x, (real(chd(ffts%r_size - i + 1, j), dp), &
                      j = 1, 3)
              else if (imode == 2) then  
                 write(18, 180) x, (aimag(chd(ffts%r_size - i + 1, j)), &
                      j = 1, 3)
              else if (imode == 3) then  
                 write(18, 180) x, (abs(chd(ffts%r_size - i + 1, j)), &
                      j = 1, 3)
              else  
                 write(18, 180) x, (abs(chd(ffts%r_size - i + 1, j) ), &
                      j = 1, 3)
              end if
           end do
        else  
           open(unit = 18, file = 'vis/'//outfilename//'.dx', &
                position = 'append', status = 'unknown', form = 'formatted')
           if (imode == 1) then  
!              write(18, 100) ((real(chd(i, j), dp), j = 1, 3), &
!                   i = ffts%r_size, 1, -1)
              write(18, 100) ((real(chd(i, j), dp), j = 1, 3), &
                   i = 1, ffts%r_size, 1)
           else if (imode == 2) then  
!              write(18, 100) ((aimag(chd(i, j)), j = 1, 3), &
!                   i = ffts%r_size, 1, -1)
              write(18, 100) ((aimag(chd(i, j)), j = 1, 3), &
                   i = 1, ffts%r_size, 1)
           else if (imode == 3) then  
!              write(18, 100) ((abs(chd (i, j) ), j = 1, 3), &
!                   i = ffts%r_size, 1, -1)
              write(18, 100) ((abs(chd(i, j) ), j = 1, 3), &
                   i = 1, ffts%r_size, 1)
           else  
!              write(18, 100) ((abs(chd(i, j)), j = 1, 3), &
!                   i = ffts%r_size, 1, -1)
              write(18, 100) ((abs(chd(i, j)), j = 1, 3), &
                   i = 1, ffts%r_size, 1)
           end if
        end if
        call myflush(18)  
        close(18)  
     end if
     call my_broadcast(ic, iproc)  
  end do

  if (.not.ikhoros) then                             !     now the dx trailer
     if (gs%myproc == 0) then  
        open(unit = 18, file = 'vis/'//outfilename//'.dx', &
             position = 'append', status = 'unknown', form = 'formatted')
        !        grid connections object
        write(18, 20) gs%fftsize(3), gs%fftsize(2), gs%fftsize(1)
        !        calculate deltas
        delta(1:3, 1) = crys%avec(1:3, 3) / gs%fftsize(3)  
        delta(1:3, 2) = crys%avec(1:3, 2) / gs%fftsize(2)  
        delta(1:3, 3) = crys%avec(1:3, 1) / gs%fftsize(1)  
        !     write header for gridpositions
        write(18, 30) gs%fftsize(3), gs%fftsize(2), gs%fftsize (1)
        !     compute and write origin
        write(18, 40) dzero, dzero, dzero  
        !     write deltas
        write(18, 50) ((delta(i, j), i = 1, 3), j = 1, 3)  
        !     write header for field object
        write(18, 60)  
        close(18)  
     end if
     write(9, 999) outfilename  
  else  
     write(9, 199) outfilename  
  endif

  deallocate (chd)  

  return  

  !     ----------------- volumetric data -----------------------------
5 format(/'# this is the object defining the data values')  
10 format('object "vdata" class array type ', &
       &     'float rank 1 shape 3 items ', &
       &     i8, ' data follows')

20 format(/'# this is the object defining the grid connections', &
       &     /'object "grd" class gridconnections counts ',3i4)

30 format(//'# this is the object defining the grid ' &
       &     /'object "pos" class gridpositions counts ',3i4 )
40 format('origin', 3f20.10)  

50 format('delta ', 3f20.10)  

60 format(/'# this is the collective object, one for each grid', &
       &     /'object "vecfield" class field', &
       &     /'component "positions"   value "pos"', &
       &     /'component "connections" value "grd"', &
       &     /'component "data"        value "vdata"')
100 format(3g15.7)  
105 format('3',/'3',3(/i6),/'3')  

180 format(3f10.6,5x, 3g15.7)  

199 format(//' GRAPHICAL OUTPUT:', &
       &     /' ----------------', &
       &     //' VECTOR FIELD FOR AVS/KHOROS WRITTEN TO FILE', &
       &     /'    - ',a7,'.avs')
999 format(//' GRAPHICAL OUTPUT:', &
       &     /' ----------------', &
       &     //' VECTOR FIELD FOR IBM DATA EXPLORER WRITTEN TO FILE', &
       &     /'    - ',a7,'.dx')

end subroutine dxplot_field
