!
!     @process extchk opt(2)
!
integer function findvec(gvec, gs)  

  include 'use.h'            
  implicit none  ! implicit? Just say no!
  include 'interface.h'  
  !
  !     function searches for a vector gvec in the gspace. If successful,
  !     it returns the index of gvec. If not, it returns -1
  !
  !     1995 Bernd Pfrommer, while at UC Berkeley
  !
  !     optimized Dec 26, 1995      Bernd
  !     the search is best performed in the translated space
  !
  !
  !     INPUT:
  !     -----
  !
  type(parallel_gspace), intent(in) :: &
       gs                                   ! the parallel gspace to search in
  integer, intent(in) :: gvec(3)                    ! the vector to search for
  !
  !     --------------------- local variables ---------------------------
  !
  integer :: i, iord, i3end, i3start, nhalf(3), irod, gv(3), &
       iordstart, iordend, iidxstart, &
       gvmod(3)                                          ! wrapped-around gvec

  findvec = -1  

  nhalf(1:3) = (gs%fftsize(1:3) - 1) / 2  
  !
  !     check if outside of the gspace anyways
  !
  if (gvec(1) > nhalf(1) .or. gvec(1) < -nhalf(1) .or. &
       gvec(2) > nhalf(2) .or. gvec(2) < -nhalf(2) .or. &
       gvec(3) > nhalf(3) .or. gvec(3) < -nhalf(3)) return
  gvmod(1) = mod(gvec(1) + gs%fftsize(1), gs%fftsize(1))  
  gvmod(2) = mod(gvec(2) + gs%fftsize(2), gs%fftsize(2))  
  gvmod(3) = mod(gvec(3) + gs%fftsize(3), gs%fftsize(3))  
  iordstart = gs%fastfind(1, gvmod(1) + 1)  
  if (iordstart < 0) return  
  iordend = gs%fastfind(2, gvmod(1) + 1)  
  iidxstart = gs%fastfind(3, gvmod(1) + 1)  
  i = iidxstart  
  do iord = iordstart, iordend  
     i3start = gs%order(2, iord)  
     i3end = gs%order(3, iord)  
     irod = gs%order(1, iord)  
     gv(1) = irod / gs%fftsize(2)  
     if (gvmod(1) == gv(1)) then  
        gv(2) = mod(irod, gs%fftsize(2))  
        if (gvmod(2) == gv(2)) then  
           if (i3start > i3end) then  
              if (gvmod(3) <= i3end .or. gvmod(3) >= i3start) then
                 findvec = i + 1 + mod(gvmod(3) + gs%fftsize(3) - i3start, &
                      gs%fftsize(3))
                 return  
              end if
           else  
              if (gvmod(3) <= i3end .and. gvmod(3) >= i3start) then
                 findvec = i + 1 + mod(gvmod(3) + gs%fftsize(3) - i3start, &
                      gs%fftsize(3))
                 return  
              end if
           end if
        end if
     end if
     i = i + mod(i3end - i3start + gs%fftsize(3), gs%fftsize(3)) + 1
  end do

  return  

end function findvec
