!
subroutine tpage(nproc)  
  !
  implicit none    ! no implicit. even in the simplest case
  !
  character(len=24) :: bdate  
  character(len=8) :: btime  
  character(len=20) :: title  
  integer, intent(in) :: nproc  
  include 'release.h'  
  !
  call zedate(bdate)  
  call zetime(btime)  
  !
  write(9, 100) release_string, date_string
  write(9,103) author_string
  if (nproc == 1) write(9, 150) bdate, nproc  

  if (nproc > 1) write(9, 160) bdate, nproc  

  return

100  FORMAT(2X,'PARATEC RELEASE ',A,' ',A)
103  FORMAT(2X,' ',A)
150 format(/15X,'RUN ',A24,' ON ',i4,' PROCESSOR'/)  
160 format(/15X,'RUN ',A24,' ON ',i4,' PROCESSORS'/)  
101 FORMAT(A20)  
102 FORMAT(//5X,A20//)

end subroutine tpage
