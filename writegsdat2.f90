!
subroutine writegsdat2(ipr, gs, data, ldnvecs, nvecs, nspin, filename, namelen)
  !
  use all_to_all_module  
  include 'use.h'  
  implicit none             ! implicit? Just say no!
  include 'interface.h'
  !
  !---------------------------------------------------------------------
  !J    for VXC and CD95 output in pwmain.fp, always the case ipr=0,
  !J    data(1,1,1:crys%spin),ldnvecs=nvecs=1. nspin should be 1 or 2
  !     writes a data array to disk such that it can be read and used in
  !     later calculations
  !
  !     .. HEADER ..
  !
  !      g1 g2 g3 ((data(i1,n,is), n=1,nvecs), is=1,nspin)
  !      g1 g2 g3 ((data(i2,n,is), n=1,nvecs), is=1,nspin)
  !      ...
  !      g1 g2 g3 ((data(il,n,is), n=1,nvecs), is=1,nspin)
  !
  !     INPUT:
  !     -----
  !
  integer, intent(in) :: &
       ipr, &              ! print flag
       namelen, &          ! length of filename
       nvecs, &            ! number of vectors
       ldnvecs, &          ! leading dimension
       nspin               ! number of spins
  type(parallel_gspace), intent(in) :: &
       gs                  ! the gspace on which the data is set up
  complex(dp), intent(in) :: &
       data(gs%length, ldnvecs, nspin)  ! data that should be written to file
  character(len=*), intent(in) :: &
       filename            ! name of file to be written
  !
  !     ---------------- local variables ---------------------------------
  !
  integer :: i, j, k, iproc, wios, n, is, igs2  
  character(len=24) :: bdate  
  character(len=8) :: btime  
  character(len=255) :: buffer  
  !
  !     ----- variables for the gspace loop ---------
  !
  integer :: igv(4), fftn(4), ffth(4), igv3, irod, iord, igs  
  !
  buffer = filename(1:namelen)  
  fftn(1:3) = gs%fftsize(1:3)  
  fftn(4) = gs%fftsize(3)  
  ffth(:) = fftn(:) / 2  
  !     -------- write header -----------------
  if (gs%myproc == 0) then  
     open(unit = 21, file = buffer, err = 99, iostat = wios, &
          status = 'unknown', form = 'unformatted')
     rewind(21)  
     call zedate (bdate)       ! get time and date
     call zetime (btime)  
     write(21) bdate, btime  
     write(21) nvecs, nspin  
     write(21) gs%nproc  
     close(21)  
  end if
  !     ------------------------------------------------------------------
  !
  !     Now the processors write gspace and gsarr values out on disk
  !
  !      do iproc = 0, gs%nproc-1     ! loop over processors
  !         if(iproc.eq.gs%myproc) then
  !            open(unit=21,file=buffer,
  !     $           err=99, iostat=wios,
  !     $           position = 'append',
  !     $           status='old',form='formatted')
  !           loop over g-vectors
  !            write(21) gs%length
  igs = 0  
  do iord = 1, gs%lorder               ! loop through x/y gspace
     irod = gs%order(1, iord)  
     igv(1) = irod / gs%fftsize(2)  
     igv(2) = mod(irod, gs%fftsize(2))  
     igv(3) = gs%order(2, iord)  
     igv(4) = gs%order(3, iord)  
     igv(:) = mod(igv(:) + ffth(:), fftn(:)) - ffth(:)  
     do igv3 = igv(3), igv(4)          ! loop over z axis
        igs = igs + 1  
        !                  write(21) igv(1),igv(2),igv3
        !                  write(9,*) igv(1),igv(2),igv3
        !                  write(21) ((data(igs,n,is),
        !     $                 n=1,nvecs), is=1,nspin)
        !                  write(9,*) ((data(igs,n,is),
        !     $                 n=1,nvecs), is=1,nspin)
     end do
  end do

  write(9, *) 'wrote', igs, 'components'  
  !            if(gs%nproc.gt.1)  close(21)
  !          end if
  !         if(gs%nproc.gt.1)  call parallel_barrier()
  ! figured out how many g-vectors on this processor
  igs2 = igs  
  do iproc = 0, gs%nproc - 1                    ! loop over processors
     if (iproc == gs%myproc) then  
        open(unit = 21, file = buffer, err = 99, iostat = wios, &
             position = 'append', status = 'unknown', form = 'unformatted')
        !           loop over g-vectors
        write(21) igs2  
        igs = 0  
        do iord = 1, gs%lorder          ! loop through x/y gspace
           irod = gs%order(1, iord)  
           igv(1) = irod / gs%fftsize(2)  
           igv(2) = mod(irod, gs%fftsize(2))  
           igv(3) = gs%order(2, iord)  
           igv(4) = gs%order(3, iord)  
           igv(:) = mod(igv(:) + ffth(:), fftn(:)) - ffth(:)  
           do igv3 = igv(3), igv(4)             ! loop over z axis
              igs = igs + 1  
              write(21) igv(1), igv(2), igv3  
              !                 write(9,*) igv(1),igv(2),igv3
              write(21) ((data(igs, n, is), n = 1, nvecs), is = 1, nspin)  
!              write(9, *) ((data(igs, n, is), n = 1, nvecs), is = 1, nspin)
           end do
        end do
        if (gs%nproc > 1) close(21)  
     end if
     if (gs%nproc > 1) call parallel_barrier()  
  end do
  if (gs%myproc == 0) then  
     if (gs%nproc > 1) then  
        open(unit = 21, file = buffer, position = 'append', &
             status = 'unknown', form = 'unformatted')
     end if
     write(21) -1234567, 0, 0       ! mark end of record
     close(21)  
  end if
  if (gs%nproc > 1) then  
     call parallel_barrier()  
  end if

  if (ipr >= 1) write(9, 100) filename  

  return

99 write(0, *) 'open failed in writegsdat!'  
  write(0, *) 'ios error code:', wios  

  call mystop  

100 format(/' SAVED DATA TO FILE ',a)  

end subroutine writegsdat2
