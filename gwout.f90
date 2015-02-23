!
subroutine gwout (fname, crys, syms, bands, gs, kgs, wfn, kpoints, &
     iflag)
  !
  use all_to_all_module  
  include 'use.h'
  implicit none                    ! implicit? Just say no!
  include 'interface.h'  
  include 'all_to_all.h'
  !  
  !     include 'mpif.h'
  !
  !     writes a data array to disk such that it can be read and used in
  !     later GW calculations
  !
  !     INPUT:
  !     -----
  !
  type(crystal), intent(in) :: crys  
  type(symmetry), intent(in) :: syms  
  type(complex_gspace_array), intent(in) :: &
       wfn                                ! the wave functions
  type(band), intent(in) :: bands  
  type(parallel_gspace), intent(in) :: &
       kgs(bands%nrk), &                  ! the gspaces at the different kpoints
       gs                                 ! the big gspaces
  type(kpoint), intent(in) :: kpoints     ! all the kpoints
  integer, intent(in) :: iflag      ! iflag =1 real version, iflag =2, complex
  character(len=*), intent(in) :: fname   ! file name
  !     ---------------- local variables ---------------------------------
  integer :: i, info, iproc, irk, is, j, k, n, ndim, ng, nkpt  
  integer :: idum(1),n_bands_out
  integer, allocatable :: ifmin(:,:), ifmax(:,:)  
  complex(dp), allocatable :: zc(:,:)

!------------------------------------------------------------
! switch band index, modify this part if needed
! switch the second and the fifth band at the Gamma point 
!
! PZ
!
! bands%gwout_index(n,irk,is)

!  if (bands%gw_band_reordering.eq.1) then
!      bands%gwout_index(2,1,:)=5
!      bands%gwout_index(5,1,:)=2
!  end if
!------------------------------------------------------

  
  if (gs%myproc == 0) then  
     !
     !     erase old file
     !
     open(21, file = fname, status = 'replace', form = 'unformatted')
     write(21) ((crys%bdot(i, j), i = 1, 3), j = 1, 3)  
     write(21) crys%vcell  
     !
     write(21) syms%ntrans  
     do n = 1, syms%ntrans  
        write(21) ((syms%mtrx(i, j, n), i = 1, 3), j = 1, 3)  
     end do
     do n = 1, syms%ntrans  
        write(21) (syms%tnp(k, n), k = 1, 3)  
     end do
     !
     write(21) crys%nspin  
     write(21) kpoints%nrk  
     write(21) (kpoints%grid(i), i = 1, 3)  
     write(21) (kpoints%shift(i), i = 1, 3)  
     write(21) (kpoints%w(i), i = 1, kpoints%nrk)  
     do j = 1, kpoints%nrk  
        write(21) (kpoints%rk(i, j), i = 1, 3)  
     end do
     !
     if (bands%num_gwout .eq. -1) then     
       write(21) bands%min(1)  
       n_bands_out=bands%min(1)
     else
       write(21) 2*bands%num_gwout       
       n_bands_out=2*bands%num_gwout
     end if

     allocate(ifmin(bands%nrk, crys%nspin))  
     allocate(ifmax(bands%nrk, crys%nspin))  
     ifmin = 0  
     ifmax = 0  
     do is = 1, wfn%nspin  
        do irk = 1, bands%nrk  
           do n = 1, bands%nband(irk, is)
!
! SIB:  1/1/2001
! changed occupancy criterion...bands%occup is actually the band occupancy
! times 2*weight, so for many k-points, it can get very small!
! Original code:
!              if(bands%occup(n,irk,is).gt.1.0d-3) then
!
              if(bands%occup(n,irk,is)/(2.0d0*kpoints%w(irk)) &
                                         .gt.1.0d-3) then
                 if (ifmin(irk, is) == 0) ifmin(irk, is) = n  
                 ifmax(irk, is) = n  
              end if
           end do
        end do
     end do

     write(21) ((ifmin(irk, is), irk = 1, bands%nrk), is = 1, crys%nspin)
     write(21) ((ifmax(irk, is), irk = 1, bands%nrk), is = 1, crys%nspin)
     do is = 1, wfn%nspin  
        do irk = 1, bands%nrk  
!           write(21) (bands%energy(n, irk, is), n = 1, bands%min(is))  
           write(21) (bands%energy( bands%gwout_index(n,irk,is),&
                      irk, is), n = 1,n_bands_out )  
        end do
     end do
     !
     close(21)  
  end if
  call my_broadcast(n_bands_out,0)
  !
  idum(1) = gs%length  
  call all_sum_all_integer(idum, 1)  
  ng = idum(1)  
  if (gs%myproc == 0) then  
     open(unit = 21, file = fname, position = 'append', status = &
          'unknown', form = 'unformatted')
     write(21) (gs%fftsize(i), i = 1, 3)  
     write(21) ng  
     close(21)  
  end if
  do iproc = 0, gs%nproc - 1  
     if (iproc == gs%myproc) then  
        open(unit = 21, file = fname, position = 'append', status = &
             'old', form = 'unformatted')
        do i = 1, gs%length  
           write(21) gs%gvec(1, i), gs%gvec(2, i), gs%gvec(3, i)  
        end do
        close(21)  
     end if
     if (gs%nproc > 1) call parallel_barrier()  
  end do
  !
  if (iflag == 1) call gwreal(bands, kgs, wfn, kpoints)  
  !
  do irk = 1, kpoints%nrk  
     ndim = kgs(irk)%length  
     idum(1) = ndim  
     call all_sum_all_integer(idum, 1)  
     nkpt = idum(1)  
     if (gs%myproc == 0) then  
        open(unit = 21, file = fname, position = 'append', status = &
             'unknown', form = 'unformatted')
        write(21) nkpt  
        close(21)  
     end if
     do iproc = 0, gs%nproc - 1  
        if (iproc == gs%myproc) then  
           open(unit = 21, file = fname, position = 'append', status = &
                'old', form = 'unformatted')
           do i = 1, kgs(irk)%length  
              write(21) kgs(irk)%gvec(1, i), kgs(irk)%gvec(2, i), &
                   kgs(irk)%gvec(3, i)
           end do
           close(21)  
        end if
        if (gs%nproc > 1) call parallel_barrier()  
     end do
     if (gs%myproc == 0) open(unit = 21, file = fname, position = &
          'append', status = 'unknown', form = 'unformatted')
     allocate(zc(nkpt, wfn%nspin))  
!     do n = 1, bands%min(1)  
     do n = 1, n_bands_out
        do is = 1, wfn%nspin  
!           call data_gather(1, gs%myproc, gs%nproc, &
!                wfn%data(1 + (n - 1) * ndim, irk, is), zc(1, is), ndim)
           call data_gather(1, gs%myproc, gs%nproc, &
                wfn%data(1 + (bands%gwout_index(n,irk,is) - 1) * ndim, irk, is), zc(1, is), ndim)
        end do
        if (iflag == 1) then  
           if (gs%myproc == 0) write(21) ((real(zc(j, k), dp), j = 1, &
                nkpt), k = 1, wfn%nspin)
        else  
           if (gs%myproc == 0) write(21) ((zc(j, k), j = 1, nkpt), &
                k = 1, wfn%nspin)
        end if
     end do
     deallocate(zc)  
     if (gs%myproc == 0) close(21)  
  end do

  if (gs%myproc == 0) then  
    deallocate(ifmin)
    deallocate(ifmax)
  end if

  return  

end subroutine gwout
