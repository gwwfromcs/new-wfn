!
! $Id: string.f90,v 1.1.1.1 2001/03/08 16:53:57 uid655 Exp $
!
! $Log: string.f90,v $
! Revision 1.1.1.1  2001/03/08 16:53:57  uid655
! initial import of version v5.0 of paratec
! M.Profeta tag v5
!
! Revision 1.1  1996/03/07  01:49:12  pfrommer
! Initial revision
!
! Revision 1.1  1991/12/16  00:07:11  alberto
! Initial revision
!
!!$subroutine chrlen(string, nchar, lchar)  
!!$  !
!!$  !***********************************************************************
!!$  !
!!$  !  CHRLEN accepts a STRING of NCHAR characters and returns LCHAR,
!!$  !  the length of the string up to the last nonblank, nonnull.
!!$  !
!!$  implicit none
!!$  !
!!$  character(len=1) :: char
!!$  character(len=*) :: string
!!$  integer, intent(in) :: nchar
!!$  integer, intent(out) :: lchar  
!!$  !
!!$  integer :: i, ncopy  
!!$  !
!!$  ncopy = nchar
!!$  if (ncopy <= 0) ncopy = len(string)  
!!$  !
!!$  do i = 1, ncopy
!!$     lchar = ncopy + 1 - i  
!!$     if (string(lchar:lchar) /= ' ' .and. &
!!$          string(lchar:lchar) /= char(0)) return
!!$  end do
!!$  lchar = 0
!!$
!!$  return
!!$  
!!$end subroutine chrlen
!
subroutine chrcap(string, nchar)  
  !
  !***********************************************************************
  !
  !  CHRCAP accepts a STRING of NCHAR characters and replaces
  !  any lowercase letters by uppercase ones.
  !
  !
  !     .. Scalar Arguments ..
  integer, intent(in) :: nchar  
  character(len=*), intent(inout) :: string
  !     ..
  !     .. Local Scalars ..
  integer :: i, itemp, ncopy  
  !     ..
  ncopy = nchar  
  if (ncopy <= 0) ncopy = len(string)  
  do i = 1, ncopy  
     !
     if (lge(string(i:i), 'a') .and. lle(string(i:i), 'z')) then
        itemp = ichar(string(i:i)) + ichar('A') - ichar('a')  
        string(i:i) = char(itemp)  
     end if
  end do
  !
  return  
  !
end subroutine chrcap
!
logical function leqi(strng1, strng2)  
  !
  !***********************************************************************
  !
  !  Case-insensitive lexical equal-to comparison
  !
  !
  !     .. Scalar Arguments ..
  character(len=*), intent(in) :: strng1, strng2
  !     ..
  !     .. Local Scalars ..
  integer :: i, len1, len2, lenc  
  character(len=1) :: s1, s2
  !     ..
  !     .. External Subroutines ..
  external chrcap  
  !     ..
  len1 = len(strng1)  
  len2 = len(strng2)  
  lenc = min(len1, len2)  
  !
  leqi = .false.  
  do i = 1, lenc  
     s1 = strng1(i:i)  
     s2 = strng2(i:i)  
     call chrcap(s1, 1)  
     call chrcap(s2, 1)  
     if (s1 /= s2) return  
  end do
  !
  if (len1 > lenc .and. strng1(lenc + 1:len1) /= ' ') return  
  if (len2 > lenc .and. strng2(lenc + 1:len2) /= ' ') return  
  leqi = .true.  
  !
  return  
  !
end function leqi
!
subroutine loc_des(message)  
  !
  !     Processes message to locate the delimiters % and $ used
  !     in the warnp routine.
  !
  !     Alberto Garcia, Feb 1, 1991
  !
  !     .. Scalar Arguments ..
  character(len=*), intent(inout) :: message
  !     ..
  !     .. Scalars in Common ..
  integer :: form_length  
  character(len=200) :: form_spec
  !     ..
  !     .. Local Scalars ..
  integer :: dol_pos, pct_pos  
  character(len=200) :: work
  !     ..
  !     .. Common blocks ..
  common / form_des / form_spec  
  common / form_len / form_length  
  !
  save / form_des /, / form_len /  
  !     ..
  pct_pos = index(message, '%')  
  !
  form_spec = form_spec(1:form_length)//message(1:(pct_pos - 1))//''','
  form_length = form_length + (pct_pos - 1) + 2  
  work = message(pct_pos + 1:)  
  !
  dol_pos = index(work, '$')  
  form_spec = form_spec(1:form_length)//work(1:dol_pos - 1)//','''
  form_length = form_length + (dol_pos - 1) + 2  
  !
  !        Return the rest of message
  !
  message = ' '  
  message = work(dol_pos + 1:)  
  !
  return  
  !
end subroutine loc_des
