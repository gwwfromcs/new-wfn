!
subroutine gspaceinverse(inv, blockdesct, lordert, ordert, nfftt)  
  !
  !     set up inversion information array inv. ONLY WORKS FOR LOCAL GSPAC
  !
  !     1995 by Bernd Pfrommer, while at UCB
  !
  implicit none             ! implicit? Just say no!
  !
  !     INPUT:
  !     -----
  !
  integer, intent(in) :: &
       blockdesct(*), &           ! block descriptor, for all procs
       nfftt(5), &                ! fft grid dimensions
       lordert, &                 ! length of order array
       ordert(4, *)               ! order array, describes the source data
  !
  !     OUTPUT:
  !     -------
  !
  integer, intent(out) :: &
       inv(*)                     ! the inversion array, pairing (G,-G)
  !
  !     ------------------ local variables ------------------------
  !
  integer :: i, ii, iords, is1, is2, is3, is3start, is3end, is, isi, gv(3)
  integer, external :: findvec
  !     ------------------ first symmetrize with symmetry operations -----

  is = 0  
  do iords = 1, lordert                     ! loop through gspace
     is1 = blockdesct(ordert(1, iords) + 1)  
     is1 = is1 - ((2 * is1) / nfftt(1)) * nfftt(1)  
     is2 = ordert(2, iords)  
     is2 = is2 - ((2 * is2) / nfftt(2)) * nfftt(2)  
     is3start = ordert(3, iords)  
     is3start = is3start - ((2 * is3start) / nfftt(3)) * nfftt(3)  
     is3end = ordert(4, iords)  
     is3end = is3end - ((2 * is3end) / nfftt(3)) * nfftt(3)  
     do is3 = is3start, is3end  
        is = is + 1  
        gv(1) = -is1  
        gv(2) = -is2  
        gv(3) = -is3  
        isi = findvec(0, blockdesct, lordert, ordert, gv, nfftt)  
        if (isi <= 0) then  
           write(9, *) 'gspaceinverse: cannot find inverse!'  
           call mystop  
        end if
        inv(is) = isi  
     end do
  end do

  return  

end subroutine gspaceinverse
