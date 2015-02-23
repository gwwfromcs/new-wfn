!
!     ===============  the primitive iolib modules =====================
!
!
!     The following integer markers have been used:
!
!     Int      Object:
!
!       1      gspace (not the parallel one)
!       2      isrt arrays
!       3      wave functions
!       4      kpoint
!       5      level
!       6      bands
!       7      crystal structure
!       8      symmetry info
!
!     If you add a new one, make sure you update iolib_skip()
!     in the iolib_module
!
!
!     ******************************************************************
!
module iolib_module  

contains  
  !
  !     ======== close an iofile =============
  !
  subroutine iolib_close(fname)  

    implicit none  

    character(len=*) :: fname  
    close(21)  
    !
    !     ! hah ! didnot use fname at all.
    !
  end subroutine iolib_close
  !
  !     =============== generic routine to find certain marker ===========
  !
  subroutine iolib_find_marker(imark, iunit, info)

    implicit none  
    !
    !     INPUT:
    !
    integer, intent(in) :: imark, &    ! integer marker
         iunit                         ! the unit to be browsed
    !
    !     OUTPUT:
    !
    integer, intent(out) :: info       ! 1 if marker found, 0 otherwise
    !
    !     ------------- local variables
    !
    integer :: marker, ierr  

    rewind(iunit)  
    !      write(0,*) 'searching for marker:', imark
    do while (.true.)  
       read(iunit, *, end = 100) marker  
       !         write(0,*) 'found marker:', marker
       ! found the right object
       if (marker == imark) then  
          info = 1  
          return  
       else  
          ! found other object. skip
          call iolib_skip(iunit, marker, ierr)  
          if (ierr /= 0) goto 100  
       end if
    end do

    return  

100 continue  
    info = 0  

    return  

  end subroutine iolib_find_marker
  !
  !     ============= routine to skip an object. has to be updated =======
  !
  !
  !
  !
  subroutine iolib_skip(iunit, imark, ierr)

    implicit none  
    !
    !     INPUT:
    !
    integer, intent(in) :: iunit, &    ! the io unit to be used
         imark                   ! the marker code of the object to be skipped
    !
    !     OUTPUT:
    !
    integer, intent(out) :: ierr  ! error code. if error, then 1. Zero otherwise
    !
    !     ------- local variables
    !
    integer :: i, ngtot, j, k, l  
    !
    !     take whatever action is necessary to skip the object
    !
    select case(imark)  
    case default  
       write(0, *) 'iolib: illegal marker: ', imark  
       write(0, *) '       is the linked iolib up to date?'  
       call mystop  
    case(1)  
       ! ------ how to skip object 1
       read(iunit, *, err = 101)  
       read(iunit, *, err = 101)  
       read(iunit, *, err = 101)  
       read(iunit, *, err = 101)  
       read(iunit, *, err = 101)  
       read(iunit, *, err = 101)  
       read(iunit, *, err = 101)  
       read(iunit, *, err = 101)  
       ierr = 0  
       return  
    case(2)  
       ! ------ how to skip object 2
       do i = 1, 4  
          read(iunit, *, err = 101)  
       end do
       ierr = 0  
       return  
    case(3)  
       ! ------ how to skip object 3
       do i = 1, 7  
          read(iunit, *, err = 101)  
       end do
       ierr = 0  
       return  
    case(4)  
       ! ------ how to skip object 4
       do i = 1, 10  
          read(iunit, *, err = 101)  
       end do
       ierr = 0  
       return  
    case(5)  
       ! ------ how to skip object 5
       do i = 1, 10  
          read(iunit, *, err = 101)  
       end do
       ierr = 0  
       return  
    case(6)  
       ! ------ how to skip object 6
       do i = 1, 6  
          read(iunit, *, err = 101)  
       end do
       ierr = 0  
       return  
    case(7)  
       ! ------ how to skip object 7
       do i = 1, 11  
          read(iunit, *, err = 101)  
       end do
       ierr = 0  
       return  
    case(8)  
       ! ------ how to skip object 8
       do i = 1, 5  
          read(iunit, *, err = 101)  
       end do
       ierr = 0  
       return  
    end select
101 continue  
    ierr = 1  

    return  

  end subroutine iolib_skip

end module iolib_module
