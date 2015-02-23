!     @process extchk
!
subroutine direct_current_process(crys, ffts, gs,pw_params, &
  nlineplot,line_plot,nsliceplot,slice_plot,bra,ket,qmag)
  !
  include 'use.h'  
  implicit none        ! implicit? Just say no!
  include 'interface.h'  
  !
  !     INPUT:
  !     -----
  !
  Type (pw_parameter)         :: pw_params   ! energy cutoffs etc
  type(crystal), intent(in) :: crys           ! for the lattice vectors
  type(fft_struc), intent(in) :: ffts  
  type(parallel_gspace), intent(in) :: gs     ! potential gspace
  complex(dp) :: &
       rho(ffts%r_size, 3, 3)                 ! raw f(r, p0, p1)
  real(dp), intent(in) :: &
       qmag, &                                ! magnitude of q
   bra(3),ket(3),    &
   line_plot(7,*),   &
   slice_plot(11,*)  
   integer, intent(in) :: &
      nlineplot,     &
      nsliceplot
  !
  !     OUTPUT:
  !     ------
  !   
  !    Jonathan Yates Paris Dec 2000
  !
  !     This plots the induced current obtained directly from the
  !     magnetic_sus_crys subroutine. 
  !     Formally this current differs from that obtained from
  !     CHI using induced_current.f90 (see Gregor Mauri Car paper)
  !     (basically because our current is not divergentless
  !     However in testing I found there was little difference between
  !     the two. You should check this, especially if you are working
  !     in the core area....
  !
  !     It's not very clean. Sorry
  !
  !     ---------------- local variables ---------------------------------
  !
  integer :: i, j, ig, p0, p1, info,ierr,nspin,nvecs
  real(dp) :: pi4vcell
  complex(dp), allocatable :: r_g(:,:,:)         ! rho_magnetic in gspace
  complex(dp), allocatable :: r_g2(:,:)         ! rho_magnetic in gspace
  complex(dp), allocatable :: r_g3(:,:)         ! rho_magnetic in gspace


  !
  !
  !
  allocate(r_g(gs%length, 3, 3), stat = info)  
  allocate(r_g2(gs%length, 9), stat = info)  
  allocate(r_g3(gs%length, 3), stat = info)  
  if (info /= 0) then  
     write(9, *) '*** ERROR: ALLOCATION FAILED IN direct_current!'  
     write(9, *) '*** CURE:  REDUCE MEMORY CONSUMPTION:'  
     call mystop  
  end if

  pi4vcell = dfour * pi / crys%vcell  
  nspin=1
  nvecs=9
  Call readgsdat(1,ierr, gs, r_g(1,1,1), 9,nvecs,nspin,'JGS',3)


  do p0 =  1,3
  do p1 =  1,3
  r_g2(:,3*(p0-1)+p1)=r_g(:,p1,p0)
  enddo
  enddo


  r_g2=-1.0*r_g2* pi4vcell* (1/qmag)  / (137.036**2 * pi4) * 2.99792458d10 

  do p0=1,gs%length                   
  if (gs%ekin(p0) <= dzero) then      !make sure G=0 term is zero
  write(9,*) 'correcting G=0 term'
!   write(9,*)   r_g2(p0,:)
  r_g2(p0,:)=0.0
  endif
  enddo


 do ig = 1, gs%length  
     do i = 1, 3  
        r_g3(ig,i) = r_g2(ig,(3*i)-2) * ket(1) + r_g2(ig,(3*i)-1) * ket(2) + &
             r_g2(ig,(3*i)) * ket(3)
     end do
  end do


      if (nlineplot >= 1) then
     call lineplot(2,crys,1.0,r_g3(1,1),gs,nlineplot,  &
         line_plot,bra,ket,3,'dirlin1')
           endif
           if (nsliceplot >= 1) then
           call sliceplot(2,crys,1.0,r_g3(1,1),gs,nsliceplot,   &
                slice_plot,bra,ket,3,0,'dirslice')
           endif 


  deallocate(r_g,r_g2)

end subroutine direct_current_process
