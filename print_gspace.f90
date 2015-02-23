!
subroutine print_gspace(ipr, gs)  
  !
  !     prints out information about the gspace
  !
  !     (1996) Bernd Pfrommer
  !
  include 'use.h'  
  implicit none           ! never
  include 'interface.h'  
  !
  !     INPUT:
  !     ------
  !
  integer, intent(in) :: ipr                 ! print level flag
  type(parallel_gspace), intent(in) :: gs    ! the gspace to be printed out
  !
  !     ----------   local variables  --------------------------
  !
  integer :: i, j, k, p

  if (ipr > 0) then  
     write(9, 1000) gs%name  
     write(9, 1005) gs%length  
     write(9, 1007) (gs%fftsize(i), i = 1, 3)  
  end if

  if (ipr > 1) then  
     if (gs%igvec) write(9, 1070)  
     if (gs%iekin) write(9, 1071)  
     write(9, 1072)  
     if (gs%istar) write(9, 1073)  
     write(9, 1002) gs%gmax**2  
     write(9, 1050) gs%lorder  
     write(9, 1013) gs%nproc  
     write(9, 1012) gs%myproc  
     write(9, 1015) (gs%fftsize(i), i = 1, 5)  
     write(9, 1080) gs%rk(1), gs%rk(2), gs%rk(3)  
  end if

  if (ipr > 2) then  
     write(9, 1060)  
     write(9, '(4(i6,1x,3i3,1x))') (j, (gs%order(k, j), k = 1, 3), &
          j = 1, gs%lorder)
     if (gs%istar) then  
        write(9, 1095)  
        if (gs%igvec) then  
           write(9, 1098) (i, gs%ekin(i), gs%inds(i), gs%mstar(gs%inds(i)), &
                gs%phase(i), (gs%gvec(j, i), j = 1, 3), i = 1, gs%length)
        else  
           write(9, 1097) (i, gs%ekin(i), gs%inds(i), gs%mstar(gs%inds(i)), &
                gs%phase(i), i = 1, gs%length)
        end if
     end if
  end if

  if (ipr > 3 .and. gs%ipackinfo) then  
     write(9, 1062)  
     do j = 1, 2                    ! packing unpacking
        write(9, 1064) 1, j  
        do p = 0, gs%nproc - 1  
           write(9, 1063) p  
           write(9, 1065) (gs%packinfo(k, p + 1, j), &
                k = 1, gs%packsize(p + 1, j, 1))
        end do
     end do
     do j = 1, 2  
        write(9, 1064) 2, j  
        do p = 0, gs%nproc - 1  
           write(9, 1063) p  
           if (j == 1) write(9, 1066)  
           if (j == 2) write(9, 1068)  
           k = 0  
           if (gs%packsize(p + 1, j, 2) > 0) then  
              do while (.true.)  
                 k = k + 1  
                 write(9, '(4i15)') gs%chunk(1:4, k, p + 1, j)  
                 if (gs%chunk(4, k, p + 1, j) < 0) exit  
              end do
           end if
        end do
     end do
  end if

  call myflush(9)
  
1000 format(/1x,a10,' gspace information:', &
       &        /' -----------------------------',/)
1002 format (' sphere radius squared (gmax^2):', f10.2)  
1005 format (' number of gvectors on this processor:', i6)  
1007 format (' FFT grid:',3i4)  
1012 format (' processor number:',i3)  
1013 format (' total number of processors:',i3)  
1015 format ('                      ',3x,'N 1',2x,'N 2',2x,'N 3', &
       &        /' gspace FFT grid:     ',5i5)
1050 format (/' order array size:',i5)  
1060 format (/' order array:')  
1062 format (/' first stage parallel packing info:')  
1063 format (/' from/to processor:',i4)  
1064 format (/' phase:',i2,', packing:',i2)  
1065 format (2i12)  
1066 format (/7x,' startrod',10x,'starty',9x,'length',10x,'width')  
1067 format (/' second transpose info:')  
1068 format (/7x,' startrod',10x,'start x',9x,'length',10x,'width')  
1070 format (' set up gspace vectors.')  
1071 format (' set up kinetic energy.')  
1072 format (' set up order array.')  
1073 format (' set up phase factors and stars.')  
1080 format (' gspace shift vector:', 3f12.6)  
1095 format (' gspace:',/' gvec #',5x, '|G|^2',7x,'star  size', &
       &      8x,'phase',16x,'gvec')
1097 format (1x,i6,3x,f10.4,3x,i4,3x,i3,2f10.4)  
1098 format (1x,i6,3x,f10.4,3x,i4,3x,i3,2f10.4,1x,3i6)  

end subroutine print_gspace
