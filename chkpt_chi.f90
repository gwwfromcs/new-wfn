!
subroutine chkpt_chi(irk, j_bare, r_size, myproc, chivv, chihh, &
     ioutput, ireadwrite)

  use constants
  implicit none                 ! implicit? Just say no!
  !
  !     1997 Bernd Pfrommer while at UCB
  !
  !     INPUT:
  !     -----
  !
  integer, intent(in) :: ireadwrite, &            ! if 1 then write, else read
       myproc, &                                             ! which proc I am
       ioutput                                      ! flag for printing output
  !
  !     INPUT/OUTPUT:
  !     -------------
  !
  integer, intent(inout) :: &
       irk, &                       ! at which kpoint the last chkpointing was
       r_size                                             ! size of rhom array
  complex(dp), intent(inout) :: chivv(3, 3), &          ! the data to be saved
       chihh(3, 3), &
       j_bare(r_size)  
  !
  !     DESCRIPTION:
  !     -----------
  !
  !     For each processor, a checkpoint file is written/read which
  !     saves/retreives the necessary state for the chi calculation
  !     ---------------- local variables ---------------------------------
  !
  real(dp) :: t0
  real(dp), external :: gimmetime  
  character(len=255) :: filename                 ! name of file to be written

  t0 = gimmetime()  

  write(filename, '(''CHKPTCHI.'',i3.3)') myproc  

  open(75, file = filename, form = 'unformatted')  
  if (ireadwrite == 1) then  
     write(75) irk, r_size  
     write(75) chivv, chihh  
     write(75) j_bare  
     close(75)  
     if (iand(ioutput, 8) == 8) write(9, 100) gimmetime() - t0
  else  
     read(75) irk, r_size  
     read(75) chivv, chihh  
     read(75) j_bare  
     close(75)  
     if (iand(ioutput, 8) == 8) write(9, 110) gimmetime() - t0
  end if

  call myflush(9)  

  return  

100 format(' TIME TO WRITE CHKPTCHI:',f12.3)  
110 format(' TIME TO READ CHKPTCHI:',f12.3)  

end subroutine chkpt_chi
