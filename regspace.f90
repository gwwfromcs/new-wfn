!
subroutine regspace(gs, datas, gt, datat)

  use all_to_all_module  
  include 'use.h'  
  implicit none              ! just say no!
  include 'interface.h'  
  include 'all_to_all.h'
  ! 
  !     1995 Bernd Pfrommer, UC Berkeley
  !
  !
  !     INPUT:
  !     -----
  !
  type(parallel_gspace), intent(in) :: &
       gs, &     ! source gspace
       gt        ! target gspace
  complex(dp), intent(in) :: &
       datas(gs%length)    ! the source data
  !
  !     OUTPUT:
  !     ------
  !
  complex(dp), intent(out) :: &
       datat(gt%length)    ! the target data
  !
  !     DESCRIPTION:
  !     -----------
  !
  !     Transfers data from source gspace to target gspace.
  !
  !     The source gspace can be distributed across multiple processors,
  !     whereas the target gspace must be replicated identically on
  !     each processor. This is the case for instance with the mixing
  !     gspace. If the source gspace is distributed, a lot of
  !     communication has to be done, so use this routine frugally.
  !
  !     Meanwhile, with get_gset we have a better routine. It scales
  !     like order(n*log(n)) with the size of the gspace and is clearly
  !     superior. Please try to use that routine instead!
  !
  !
  !
  !     ==================================================================
  !
  !     local variables:
  !
  complex(dp) :: dcdum, rdat  
  logical :: packit  
  integer :: i, j, k, svec(3)  
  !
  !     ----- variables for the gspace loop ---------
  !
  integer :: igvt(4), fftnt(4), fftht(4), igv3t, irodt, iordt, &
       igt, igs
  integer, external :: findvec    ! finds index of gvec in gspace

  if (gt%nproc > 1) then  
     write(9, *) 'regspace: target gspace must be replicated!'  
     call mystop  
  end if
  fftnt(1:3) = gt%fftsize(1:3)  
  fftnt(4) = gt%fftsize(3)  
  fftht(:) = fftnt(:) / 2  
  !
  ! clear the target array
  !
  datat = zzero
  !     ----------------------------------------------------------------
  !     Pick up the desired components while looping through
  !     the target gspace. when the target  gspace is smaller, this is
  !     more efficient than looping through the source space, because
  !     searching a gspace scales like n^(2/3).
  !
  igt = 0  
  do iordt = 1, gt%lorder                    ! loop through x/y gspace
     irodt = gt%order(1, iordt)  
     igvt(1) = irodt / gt%fftsize(2)  
     igvt(2) = mod(irodt, gt%fftsize(2))  
     igvt(3) = gt%order(2, iordt)  
     igvt(4) = gt%order(3, iordt)  
     igvt(:) = mod(igvt(:) + fftht(:), fftnt(:)) - fftht(:)  
     do igv3t = igvt(3), igvt(4)             ! loop over z axis
        igt = igt + 1  
        !
        !              search in the source gspace
        !
        svec(1:2) = igvt(1:2)  
        svec(3) = igv3t  
        igs = findvec(svec, gs)  
        if (igs > 0) datat(igt) = datas(igs)  
        !            write(9,*) igs,datas(igs),igt,datat(igt)
     end do
  end do
  !     In case the source gspace is distributed, the target array is
  !     completed by a global sum. This is very brute force, but the
  !     simplest to program. Check if this causes an unacceptable
  !     amount of communication.
  !
  if (gs%nproc > 1) call all_sum_all(datat, gt%length)  

end subroutine regspace
