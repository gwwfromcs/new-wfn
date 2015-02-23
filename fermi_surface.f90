!
subroutine fermi_surface(myproc, pw_params, crys, energs, kp, bands)
  !
  use all_to_all_module  
  include 'use.h'  
  implicit none  
  include 'interface.h'  
  include 'all_to_all.h'  
  !
  !     2001 David Roundy
  !
  !     INPUT:
  !     -----
  !
  type(crystal), intent(in) ::  crys    ! the crystal structure
  type(band), intent(in) :: bands  
  type(energy), intent(in) :: energs  
  type(kpoint), intent(in) :: kp  
  type(pw_parameter), intent(in) :: pw_params  
  integer, intent(in) :: myproc  
  !
  !     DESCRIPTION:
  !     -----------
  !
  !  writes a few files which have the fermi surface and some reciprocal
  !  space stuff.
  !
  !     -------------------local variables ------------------------------
  !
  integer :: irk, is, n, it, ia 

! David's variables:
  integer :: myband, bandnum, i, j, k 
  real(dp) :: delta(3, 3)

  write(9, 110)  
  if (myproc /= 0) return  

! First output the reciprocal lattice vectors...
  call dx_reciprocal_lvec(crys, myproc)

! Then output all the eigenvalues as a function of k (so an isosurface
! at zero gives you the fermi surface).
  open(13, file = 'vis/fermi.dx', status = 'unknown', form = 'formatted') 
  is = 1
  ! Write dx trailer
  write(13, 20) kp%grid(1), kp%grid(2), kp%grid (3)
  !        calculate deltas
  delta(1:3, 1) = crys%bvec(1:3, 1) / kp%grid(1)  
  delta(1:3, 2) = crys%bvec(1:3, 2) / kp%grid(2)  
  delta(1:3, 3) = crys%bvec(1:3, 3) / kp%grid(3)  
  !     write header for gridpositions
  write(13, 30) kp%grid(1), kp%grid(2), kp%grid(3)
  !     compute and write origin
  write(13, 40) dzero, dzero, dzero
  !     write deltas
  write(13, 50) ((delta(i, j), i = 1, 3), j = 1, 3)  
  write(9,*) 'plotting' , bands%min(1), 'bands'
  call myflush(9)
  do myband=1,bands%min(1)
     n = 0
   ! Write dx header
     write(13, 5)                          ! write the header for the array
     write(13, 10) 2+myband, kp%grid(1) * kp%grid(2) * kp%grid(3)
     ! End of dx header
     do i = 1, kp%grid(1)
        do j = 1, kp%grid(2) 
           do k = 1, kp%grid(3) 
              n = n + 1  
              ! negative indices indicate irred points! (see generate_kpoints)
              irk = abs(kp%kmap(n))
              write(13, *) (bands%energy(myband, irk, is)-energs%efermi(1) &
                   +energs%evalshift(1))
           end do
        end do
     end do
     !     write header for field object
     write(13, 60) myband, 2+myband
     call myflush(13)
  end do
  close(13)
  
  return  

110 format(/' Outputting the Fermi surface ',/1x,31('-'))  

  !     ----------------- dx output -----------------------------
5 format(/'# this is the object defining the data values')  
10 format('object ',i4,' class array type float rank 0 items ',i8, &
       &     ' data follows')

20 format(/'# this is the object defining the grid connections', &
       &     /'object 1 class gridconnections counts ',3i4)

30 format(//'# this is the object defining the grid ', &
       &     /'object 2 class gridpositions counts ',3i4 )
40 format('origin', 3f20.10)  

50 format('delta ', 3f20.10)  


60 format(/'# this is the collective object, one for each grid', &
       &     /'object "',i3,'" class field', &
       &     /'component "positions"   value 2', &
       &     /'component "connections" value 1', &
       &     /'component "data"        value ',i4)
80 format(5g15.7)  

end subroutine fermi_surface


