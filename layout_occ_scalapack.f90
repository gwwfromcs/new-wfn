!
subroutine layout_occ_scalapack(matsize, b, nproc, p, q)  
  !
  use constants
  implicit none
  !
  !     INPUT:
  !     -----
  !
  integer, intent(in) :: nproc, &        ! number of processors
       matsize                           ! size of matrix
  !
  !     OUTPUT:
  !     ------
  !
  integer, intent(out) :: b, &           ! block size
       p, &                              ! processor grid row
       q                                 ! processor grid column
  !
  !     DESCRIPTION:
  !     -----------
  !
  !     Figures out a  p by q processor grid layout for the scalapack
  !     library. This p*q grid is used to partition the matrix with
  !     a block size b. For more, see the scalapack documentation
  !
  !     The goal is to get a processor grid which is as close to
  !     "square" as possible.
  !
  !     1997 Andrew Canning, Bernd Pfrommer
  !     2001 modified for occupied space by DBR 
  !     
  !     ----------------------- local variables ----------------------
  !
  integer :: i, b_t, p_max,q_max  
  !
  !     Find maximum processor grid size
  !
  p = int(sqrt(real(nproc, dp) + 1.0d-6))  
  do i = p, 1, -1  
     if (mod(nproc, i) == 0) exit  
  end do

  p_max = i  
  q_max = nproc / p_max
  !
  !  if matsize < 250 use only 1 processor, otherwise ensure block size is
  !   between 40 and 60 with approximately equal amounts on each processor
  !
  if (matsize .gt. 250) then

    do i=q_max,1,-1
      b_t=matsize/i
      if (mod(matsize, i) .ne. 0) b_t=b_t+1
      if (b_t .gt. 40) exit
    end do

    q_max=i
    p_max=min(q_max,p_max)

!   ensure that b_t is large enough such that the p_max processors can handle
!   all of the cols given the block size
    if (b_t .lt. matsize/p_max) b_t=matsize/p_max  
  
    do i=1,b_t
      if ( (b_t)/i .le. 60) then
        b= (b_t)/i
        go to 10
      end if
    end do

  else
    b=matsize
  end if 

10 continue
  !
  !  ensure enough processors are used.
  !
  q= matsize/b
  if (q .eq. 0) then 
    q=1
  else if (mod(matsize, b) .ne. 0) then
    q=q+1  
  end if 
  q=min(q,q_max)

  p= matsize/b
  if (p .eq. 0) then 
    p=1
  else if (mod(matsize, b) .ne. 0) then 
    p=p+1
  end if

  p=min(p,p_max)

  return  

end subroutine layout_occ_scalapack
