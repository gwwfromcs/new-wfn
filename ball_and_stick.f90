!
subroutine ball_and_stick(ioutput, crys, myproc)  
  !
  include 'use.h'  
  implicit none                     ! implicit? Just say no!
  include 'interface.h'  
  !
  !     ------------------------------------------------------------------
  !            writes a data file for the ball and stick model
  !     to be read by the appropriate module from the IBM data explorer
  !     a visualization software package.
  !
  !
  !
  !     1996 Bernd Pfrommer
  !
  !     INPUT
  !     -----
  !
  type(crystal), intent(in) :: crys  
  integer, intent(in) :: myproc, ioutput  
  !
  !     ------------- local variables
  !
  integer :: i, j, k, ia, it, irat(3)  
  real(dp) :: rat(3)  
  logical :: ikhoros  

  ! want only proc 0 to do this part

  if (myproc /= 0) return  
  if ((iand(ioutput, 32) /= 32) .and. (iand(ioutput, 16384) /=16384)) return  

  call dx_latticevec(crys, myproc)  

  ikhoros = (iand(ioutput, 32768) == 32768)       ! if khoros or dx

  if (ikhoros) then
     open(19, file = 'vis/BALLS.avs', status = 'unknown', form = &
          'formatted')
  else  
     open(19, file = 'vis/BALLS.dx', status = 'unknown', form = &
          'formatted')
     write(19, 10)  
  end if
  !     write position array
  k = 0  
  do it = 1, crys%ntype  
     do ia = 1, crys%natom(it)  
        k = k + 1  
     end do
  end do
  if (ikhoros) then  
     write(19, 120) 1, 3, k, 1  
  else  
     write(19, 20) k  
  end if

  do it = 1, crys%ntype  
     do ia = 1, crys%natom(it)  
        rat = crys%rat(:, ia, it) / pi2  
        if (ikhoros) then                   ! translate to first octant
           rat(:) = rat(:) + 1.0d3 
           irat(:) = int(rat(:))  
           rat(:) = rat(:) - real(irat(:), dp)  
        end if
        if (ikhoros) then  
           write(19, '(3f17.10,f6.1)') matmul(crys%avec, rat), real(it, dp)
        else  
           write(19, '(3f17.10)') matmul(crys%avec, rat)  
        end if
     end do
  end do
  !     write the data field with the types of atoms
  if (.not.ikhoros) then  
     write(19, 50) k  
     do it = 1, crys%ntype  
        do ia = 1, crys%natom(it)  
           write(19, '(f12.4)') real(it, dp)  
        end do
     end do
     write(19, 60)  
  end if

  close(19)  

  return  
  !
  !     --------- dx formats --------------------
  !
10 format(2('#',/),'#    BALL AND STICK INFO:',i4,2(/'#'))  

20 format('object "ballcoord" array type float rank 1 shape 3', &
       &     ' items ',i6,' data follows',/, &
       &     '# the positions of the atoms')
50 format('object "data" array type float rank 0 items',i6, &
       &     ' data follows')
60 format(' object "molecule" field', &
       &     /' component "positions" value "ballcoord"' &
       &     /' component "data" value "data"', &
       &     /'end')
  !
  !     --------- khoros formats --------------------
  !
120 format(i6/,i6/,i6/,i6)  

end subroutine ball_and_stick
