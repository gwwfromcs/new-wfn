!*
!     @process extchk
!
subroutine dx_latticevec(crys, myproc)
  
  include 'use.h'  
  implicit none        ! implicit? Just say no!
  include 'interface.h'  
  !
  !     ------------------------------------------------------------------
  !     writes the file LATTICE_VEC.dx for the IBM data explorer,
  !     a visualization software package.
  !
  !
  !     1996 Bernd Pfrommer
  !
  !     added the UCELL_FRAME.dx file 08/10/96
  !
  !
  !     INPUT
  !     -----
  !
  type(crystal), intent(in) :: crys  
  integer, intent(in) :: myproc  
  !
  !     ------------- local variables
  !
  integer :: i, j  
  logical :: ex  
  
  if (myproc /= 0) return                 ! want only proc 0 to do this part

  inquire(file = 'vis', exist = ex)

  if (.not.ex) call mysystem('mkdir vis')  

  open(19, file = 'vis/LATTICE_VEC.dx', status = 'unknown', &
       form = 'formatted')
  write(19, 10)  
  write(19, 20) 1, 3  
  write(19, 30) ((crys%avec(i, j), i = 1, 3), j = 1, 3)  
  write(19, 20) 2, 3  
  write(19, 35)  
  write(19, 40) 3, 1, 2  
  write(19, 50)  
  close(19)
  
  open(19, file = 'vis/UCELL_FRAME.dx', status = 'unknown', &
       form = 'formatted')
  write(19, 60)  
  write(19, 15) 3, 12  
  write(19, '(2i3)') 0, 1, 0, 2, 0, 3, 1, 4, 1, 5, 3, 5, 3, 6, 2, &
       6, 2, 4, 7, 5, 7, 6, 7, 4
  write(19, *) 'attribute "element type" string "lines"'  
  write(19, 20) 4, 8  
  write(19, 30) dzero, dzero, dzero  
  write(19, 30) (crys%avec(i, 1), i = 1, 3)  
  write(19, 30) (crys%avec(i, 2), i = 1, 3)  
  write(19, 30) (crys%avec(i, 3), i = 1, 3)  
  write(19, 30) (crys%avec(i, 1) + crys%avec(i, 2), i = 1, 3)  
  write(19, 30) (crys%avec(i, 1) + crys%avec(i, 3), i = 1, 3)  
  write(19, 30) (crys%avec(i, 2) + crys%avec(i, 3), i = 1, 3)  
  write(19, 30) (crys%avec(i, 1) + crys%avec(i, 2) + crys%avec(i, 3), i = 1, 3)
  write(19, 70) 5, 12  
  write(19, '(12(''1.0''/))')  
  write(19, *) 'attribute "dep" string "connections"'  
  write(19, 45) 6, 5, 4, 3  
  write(19, 50)  
  close(19)  

  return
  
10 format(2('#',/),'#    LATTICE VECTOR INFO:',2(/'#'))  
15 format('object ',i2,' class array type ', &
       &     'int rank 1 shape 2 items ',i4,' data follows')
20 format('object ',i2,' class array type ', &
       &     'float rank 1 shape 3 items ',i4,' data follows')
30 format(3(f15.8))  
35 format(3('0 0 0'/))  
40 format('object ',i3,' class field',/'component "data" value ',i2, &
       &     /'component "positions" value ',i2)
45 format('object ',i3,' class field',/'component "data" value ' &
       &     ,i2,/'component "positions" value ',i2, &
       &     /'component "connections" value ',i2)

50 format('end')  
60 format(2('#',/),'#   UNIT CELL FRAME:',2(/'#'))  

70 format('object ',i4,' array type float rank 0 items ',i4 &
       &     ,' data follows')

end subroutine dx_latticevec
