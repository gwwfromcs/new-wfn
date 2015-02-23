!
subroutine get_gset(ngsmlist, gsmlist, nvec, vsub, gmaxsub, gs, v)
  !
  use all_to_all_module  
  include 'use.h'  
  implicit none                         ! implicit? Just say no!
  include 'interface.h'  
  include 'all_to_all.h'  
  !
  !
  !     1997 Bernd Pfrommer, Andrew Canning
  !
  !     INPUT:
  !     -----
  !
  integer, intent(in) :: nvec, &     ! number of vectors in v and vsub
       ngsmlist, &                   ! number of G vectors in gsmlist
       gsmlist(ngsmlist)             ! list of G vectors for which
                                     ! v(G) is requested. The G vectors are
                                     ! specified by a global ordering index
  type(parallel_gspace), intent(in) :: &
       gs                            ! large parallel gspace for source data
  real(dp), intent(in) :: gmaxsub    ! max g length for the small gspace
  complex(dp), intent(in) :: &
       v(gs%length, nvec)            ! the source data, given on the large gs
  !
  !     OUTPUT:
  !     ------
  !
  complex(dp), intent(out) :: &
       vsub(ngsmlist, nvec)          ! the output data for the g vectors
                                     ! requested in the gsmlist
  !
  !     DESCRIPTION:
  !     -----------
  !
  !     Given a complex data field with nvec vectors on a large
  !     distributed g-space,
  !     this routine retrieves the data on a subset of G-vectors inside
  !     a smaller sphere of size emaxsub. From those, it only returns
  !     the ones specified in gsmlist, which are specified by a global
  !     index common to both the large and small gspace.
  !
  !     ==================================================================
  !
  !     -------------------- local variables -----------------------------
  !
  complex(dp), allocatable :: &
       temp(:), &                    ! temporary var. to receive vsubloc
       vsubloc(:,:)            ! potential for gvecs<esub residing on myproc
  real(dp) :: emaxsub  
  integer :: igssub, i, j, &
       ivec, &                       ! index for vector number
       iproc, nsize, il, igidx, &
       ngllist, &                    ! number of gvecs<esub residing on myproc
       ngllistmax                    ! largest ngllist across all processors
  integer, allocatable :: &
       itemp(:), &                   ! temp to receive gllist from other procs
       gsmsort(:), &                 ! sorting index array for output gvec list
       glsort(:), &              ! sorting index array for gvecs<esub on myproc
       gllist(:)                 ! global order index for those gvectors<esub
                                 ! residing on this processor
  !
  !     ----- variables for the gspace loop ---------
  !
  integer :: igv(4), fftn(4), ffth(4), igv3, irod, iord, igs  

  emaxsub = gmaxsub * gmaxsub  
  !
  !     run through my portion of the large, distributed  gspace,
  !     and find all g-vectors within the smaller gspace.
  !
  fftn(1:3) = gs%fftsize(1:3)  
  fftn(4) = gs%fftsize(3)  
  ffth(:) = fftn(:) / 2  

  igs = 0 ; igssub = 0
  do igs = 1, gs%length  
     if (gs%ekin(igs) <= emaxsub) igssub = igssub + 1  
  end do

  ngllist = igssub  
  allocate(vsubloc(ngllist, nvec))  
  allocate(gllist(ngllist))  
  igs = 0  
  igssub = 0  
  do iord = 1, gs%lorder                        ! loop through x/y gspace
     irod = gs%order(1, iord)  
     igv(1) = irod / gs%fftsize(2)  
     igv(2) = mod(irod, gs%fftsize(2))  
     igv(3) = gs%order(2, iord)  
     igv(4) = gs%order(3, iord)  
     igv(:) = mod(igv(:) + ffth(:), fftn(:)) - ffth(:)  
     do igv3 = igv(3), igv(4)                   ! loop over z axis
        igs = igs + 1                           ! index into large gspace
        if (gs%ekin(igs) <= emaxsub) then  
           igssub = igssub + 1  
           gllist(igssub) = mod(igv3 + fftn(3), fftn(3)) + &
                mod(igv(2) + fftn(2), fftn(2)) * fftn(3) + &
                mod(igv(1) + fftn(1), fftn(1)) * fftn(3) * fftn(2)
           vsubloc(igssub, 1:nvec) = v(igs, 1:nvec)  
        end if
     end do
  end do
  !
  !     sort the target gsmlist according to its global ordering values
  !
  allocate(gsmsort(ngsmlist))  
  do i = 1, ngsmlist
     gsmsort(i) = i
  end do
  call iheapsort(ngsmlist, gsmlist(1), gsmsort(1))  
  !
  !     sort the list according to the gllist, and rearrange vsubloc and
  !     gllist.
  !
  allocate(glsort(ngllist))  
  do i = 1, ngllist
     glsort(i) = i
  end do

  if (ngllist > 0) call iheapsort(ngllist, gllist(1), glsort(1))
  ngllistmax = ngllist  

  call all_max_all(ngllistmax)  
!!!! make temp bigger
  allocate(temp(ngllistmax * nvec))  
  do i = 1, nvec  
     do j = 1, ngllist  
        temp(ngllist * (i - 1) + j) = vsubloc(glsort(j), i)  
     end do
     do j = 1, ngllist  
        vsubloc(j, i) = temp(ngllist * (i - 1) + j)  
     end do
  end do
  allocate(itemp(ngllistmax))  
  if (ngllist > 0) then  
     itemp(1:ngllist) = gllist(glsort(1:ngllist))  
     gllist(1:ngllist) = itemp(1:ngllist)  
  end if
  deallocate(glsort)  
  !
  !     now loop through all processors, and get their vsubloc
  !
  vsub = zzero
  do iproc = 0, gs%nproc - 1  
     nsize = ngllist  
     call my_broadcast(nsize, iproc)       ! broadcast the size of the list
     if (nsize == 0) cycle                 ! nothing from this guy
     do i = 1, nvec  
        temp(ngllist * (i - 1) + 1:ngllist * i) = vsubloc(1:ngllist, i)  
     end do

     itemp(1:ngllist) = gllist(1:ngllist)  
     call my_broadcast(temp, nsize * nvec, iproc)     !  broadcast vsubloc
     call my_broadcast(itemp(1), nsize, iproc)        !  broadcast index list
     !
     !        now pick up components matching gsmlist
     !
     il = 1  
     do i = 1, ngsmlist  
        igidx = gsmlist(gsmsort(i))  
        !
        !           search for match and grab vsub if it finds it
        !
        do while (itemp(il) < igidx .and. il < nsize)  
           il = il + 1  
        end do
        if (itemp(il) == igidx) then  
           do ivec = 1, nvec  
              vsub(gsmsort(i), ivec) = temp((ivec - 1) * nsize + il)  
           end do
        end if
        if (il == nsize .and. igidx > itemp(il)) exit  
     end do
  end do
  deallocate(temp)  
  deallocate(itemp)  
  deallocate(vsubloc)  
  deallocate(gllist)  
  deallocate(gsmsort)  

  return

end subroutine get_gset
