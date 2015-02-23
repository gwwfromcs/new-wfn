!
subroutine expand(gsub, gfull, input, output, unfolded)  
  !
  !     1996 Bernd Pfrommer while at  UC Berkeley
  !
  include 'use.h'  
  implicit none  
  include 'interface.h'
  !
  !     INPUT:
  !     -----
  !
  type(parallel_gspace), intent(in) :: gsub, &                    ! sub gspace
       gfull                                                     ! full gspace
  complex(dp), intent(in) :: input(gsub%length)           ! input data, folded
  !
  !     OUTPUT:
  !     ------
  !
  complex(dp), intent(out) :: output(gfull%length)   ! the output data, folded
  !
  !     WORK:
  !     -----
  !
  complex(dp) :: unfolded(gsub%fftsize(3), gsub%fftsize(2), gsub%fftsize(1))
  !
  !     DESCRIPTION:
  !     -----------
  !
  !     Takes data on a small, replicated gspace, and sets it on a
  !     large, possibly distributed gspace.
  !
  !     To avoid any searches, this is done as follows:
  !
  !      1) Set the data from the small gspace on a fourier grid
  !      2) Pick up the data from the fourier grid into the large gspace
  !
  !
  !     As a consequence, the expansion scales only linearly with the
  !     size of the gspace.
  !
  !
  !     ==================================================================
  !
  !     local variables:
  integer :: ii, istart, iend, iid, nz, i, ix, iy, fftmargin(3, 2), &
       ift(3), &                               ! pointer into the fourier grid
       ifg(3)
  !
  !     ----- variables for the gspace loop ---------
  !
  integer :: igv(4), fftn(4), ffth(4), igv3, irod, iord, igs

  if (gsub%nproc > 1) then  
     write(9, *) 'expand: works only for local subspaces'  
     call mystop  
  end if

  unfolded = zzero                                      ! clear unfolded array
  !
  !     first put data from small space to fourier grid
  !
  nz = gsub%fftsize(3)  
  ii = 1  
  do iord = 1, gsub%lorder  
     irod = gsub%order(1, iord)  
     ix = irod / gsub%fftsize(2) + 1  
     iy = mod(irod, gsub%fftsize(2)) + 1  
     istart = gsub%order(2, iord)  
     iend = gsub%order(3, iord)  
     if (istart <= iend) then     ! there is only one contiguous block to copy
        iid = ii + mod(iend - istart + nz, nz)  
        unfolded(istart + 1:iend + 1, iy, ix) = input(ii:iid)  
     else                                       ! there are two blocks to copy
        iid = ii + nz - istart - 1  
        unfolded(istart + 1:nz, iy, ix) = input(ii:iid)  
        ii = iid + 1  
        iid = ii + iend  
        unfolded(1:iend + 1, iy, ix) = input(ii:iid)  
     end if
     !         write(9,*) iord, unfolded(:,iy,ix)
     ! new start  =    old end +1
     ii = iid + 1  
  end do
  !
  !     now pick up the stuff for the large gspace if
  !
  !            -(fftsub-1)/2 < index < (fftsub+1)/2
  fftn(1:3) = gfull%fftsize(1:3)  
  fftn(4) = gfull%fftsize(3)  

  ffth(:) = fftn(:) / 2  
  fftmargin(:, 1) = -(gsub%fftsize(1:3) - 1) / 2  
  fftmargin(:, 2) = (gsub%fftsize(1:3) - 1) / 2  
  !      write(9,*)'sub fftsize:',gsub%fftsize(1:3)
  !      write(9,*)'fftmargin:',fftmargin
  igs = 0  
  do iord = 1, gfull%lorder                         ! loop through x/y gspace
     irod = gfull%order(1, iord)  
     igv(1) = irod / gfull%fftsize(2)  
     igv(2) = mod(irod, gfull%fftsize(2))  
     igv(3) = gfull%order(2, iord)  
     igv(4) = gfull%order(3, iord)  
     igv(:) = mod(igv(:) + ffth(:), fftn(:)) - ffth(:)  
     do igv3 = igv(3), igv(4)                              ! loop over z axis
        igs = igs + 1  
        ! within range of sub gspace
        if (fftmargin(1, 1) <= igv(1) .and. igv(1) <= fftmargin(1, 2) .and. &
             fftmargin(2, 1) <= igv(2) .and. igv(2) <= fftmargin(2, 2) .and. &
             fftmargin(3, 1) <= igv3 .and. igv3 <= fftmargin(3, 2)) then
           ift(1:2) = mod(igv(1:2) + gsub%fftsize(1:2), gsub%fftsize(1:2))
           ift(3) = mod(igv3 + gsub%fftsize(3), gsub%fftsize(3))  
           output(igs) = unfolded(ift(3) + 1, ift(2) + 1, ift(1) + 1)
           !               write(9,*) igs,ift(:), igv(1:2),igv3,output(igs)
        else  
           output(igs) = zzero
        end if
     end do
  end do

  return

end subroutine expand
