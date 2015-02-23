!
subroutine put_gset(ngsmlist, gsmlist, nvec, veclist, vsub, &
     gmaxsub, gs, v, nvvec)

  use all_to_all_module  
  include 'use.h'  
  implicit none             ! implicit? Just say no!
  include 'interface.h'  
  include 'all_to_all.h'  
  !
  !
  !     1997 Andrew Canning, Bernd Pfrommer
  !
  !     INPUT:
  !     -----
  !
  integer, intent(in) :: &
       nvvec, &                 ! number of vectors allowed in v
       nvec, &                  ! number of vectors hold by myproc
       veclist(nvec), &         ! list of vectors on myproc
       ngsmlist                 ! number of vectors in gsmlist
  type(parallel_gspace), intent(in) :: &
       gs                       ! large parallel gspace for the target
  real(dp), intent(in) :: &
       gmaxsub                  ! max g length for the small gspace
  complex(dp), intent(in) :: &
       vsub (ngsmlist, nvec)    ! the data for the g vectors
                                ! listed in the gsmlist.
  integer, intent(in) :: &
       gsmlist(ngsmlist)        ! list of G vectors for which
                                ! vsub(G) is given. The G vectors are
                                ! specified by a global ordering index.
  !     OUTPUT:
  !     ------
  !
  complex(dp), intent(out) :: &
       v(gs%length, nvvec)      ! the source data, given on the large gs
  !
  !     DESCRIPTION:
  !     -----------
  !
  !     Given a complex data field on a small, distributed g-space,
  !     with a maximum g length of gmaxsub, this routine puts the data
  !     on a large, distributed gspace. The input data field vsub can
  !     consist of multiple (neig) vectors (e.g. the different bands at a
  !     a given k-point), and the veclist describes which ones
  !     this processor holds. Finally, gsmlist contains the list of
  !     global g-vector indices for which vsub is given.
  !
  !     While the (distributed) input data consisting of nvec vectors
  !     can refer to a list of arbitrary vectors (given by veclist),
  !     those which have larger index than nvvec will be ignored.
  !
  !
  !
  !     ==================================================================
  !
  !     -------------------- local variables -----------------------------
  !
  complex(dp), allocatable :: temp(:)     ! temporary var. to receive vsub
  real(dp) :: emaxsub  
  integer :: igssub, i, ivec, j, &
       nvecr, &                   ! number of vectors on remote processor
       nvecmax, &                 ! maximum number of vectors on all processors
       iproc, nsize, il, igidx, &
       ngllist, &                 ! number of gvecs<esub residing on myproc
       ngsmlistmax                ! largest ngllist across all processors
  integer, allocatable :: &
       veclisttemp(:), &          ! temp to receive veclist from other procs
       glidx(:), &                ! maps subspace to large gspace. This is
                                  ! used to finally place the data.
       itemp(:), &                ! temp to receive gllist from other procs
       gsmsort(:), &              ! sorting index array for output gvec list
       glsort(:), &               ! sorting index array for gvecs<esub on myproc
       gllist(:)                  ! global order index for those gvectors<esub
                                  ! residing on this processor
  !
  !     ----- variables for the gspace loop ---------
  !
  integer :: igv(4), fftn(4), ffth(4), igv3, irod, iord, igs

  emaxsub = gmaxsub * gmaxsub  
  nvecmax = nvec  

  call all_max_all(nvecmax)  
  !
  !     run through my portion of the large, distributed  gspace,
  !     and find all g-vectors within the smaller gspace.
  !
  !
  !     to dimension arrays, first just count how many there are.
  !
  igssub = 0
  do igs = 1, gs%length  
     if (gs%ekin(igs) <= emaxsub) igssub = igssub + 1  
  end do

  ngllist = igssub  
  allocate(gllist(ngllist))  
  allocate(glidx(ngllist))
 
  fftn(1:3) = gs%fftsize(1:3)  
  fftn(4) = gs%fftsize(3)  
  ffth(:) = fftn(:) / 2  
  igs = 0  
  igssub = 0  
  do iord = 1, gs%lorder                    ! loop through x/y gspace
     irod = gs%order(1, iord)  
     igv(1) = irod / gs%fftsize(2)  
     igv(2) = mod(irod, gs%fftsize(2))  
     igv(3) = gs%order(2, iord)  
     igv(4) = gs%order(3, iord)  
     igv(:) = mod(igv(:) + ffth(:), fftn(:)) - ffth(:)  
     do igv3 = igv(3), igv(4)               ! loop over z axis
        igs = igs + 1                       ! index into large gspace
        if (gs%ekin(igs) <= emaxsub) then  
           igssub = igssub + 1  
           gllist(igssub) = mod(igv3 + fftn(3), fftn(3)) + &
                mod(igv(2) + fftn(2), fftn(2)) * fftn(3) + &
                mod(igv(1) + fftn(1), fftn(1)) * fftn(3) * fftn(2)
           glidx (igssub) = igs  
        end if
     end do
  end do
  !
  !     sort the source gsmlist according to its global ordering values
  !     and reorder vsub/gsmlist accordingly.
  !
  allocate(gsmsort(ngsmlist))  
  do i = 1, ngsmlist
     gsmsort(i) = i
  end do
  call iheapsort(ngsmlist, gsmlist(1), gsmsort(1))  
  ngsmlistmax = ngsmlist  
  call all_max_all(ngsmlistmax)  
  allocate(temp(ngsmlistmax * nvecmax))  
  allocate(itemp(ngsmlistmax))  
  !
  !     sort the target list according to gllist
  !
  allocate(glsort(ngllist))  
  do i = 1, ngllist
     glsort(i) = i
  end do

  if (ngllist > 0) call iheapsort(ngllist, gllist(1), glsort(1))

  allocate(veclisttemp(nvecmax))  
  !
  !     now loop through all processors, and get their v
  !
  do iproc = 0, gs%nproc - 1  
     nsize = ngsmlist  
     nvecr = nvec  
     if (iproc == gs%myproc) then  
        do i = 1, nvec  
           do j = 1, ngsmlist  
              temp(ngsmlist * (i - 1) + j) = vsub(gsmsort(j), i)  
           end do
        end do
        do j = 1, ngsmlist  
           itemp(j) = gsmlist(gsmsort(j))  
        end do
        veclisttemp(1:nvec) = veclist(1:nvec)  
     end if
     call my_broadcast(nvecr, iproc)       ! bcast the number of vectors
     call my_broadcast(nsize, iproc)       ! broadcast the size of the list
     call my_broadcast(veclisttemp(1), nvecr, iproc) ! the list of v
     call my_broadcast(itemp(1), nsize, iproc)       !  broadcast index list
     call my_broadcast(temp, nsize * nvecr, iproc)   !  broadcast vsub
     !
     !        now pick up components matching gllist
     !
     il = 1  
     do i = 1, ngllist  
        igidx = gllist(glsort(i))  
        !
        !           search for match and grab vsub if it finds it
        !
        do while (itemp(il) < igidx .and. il < nsize)  
           il = il + 1  
        end do
        if (itemp(il) == igidx) then  
           do ivec = 1, nvecr  
              if (veclisttemp(ivec) <= nvvec) v(glidx(glsort(i)), &
                   veclisttemp(ivec)) = temp((ivec - 1) * nsize + il)
           end do
        end if
        if (il == nsize .and. igidx > itemp(il)) exit  
     end do
  end do

  deallocate(gsmsort)  
  deallocate(veclisttemp)  
  deallocate(glidx)  
  deallocate(temp)  
  deallocate(itemp)  
  deallocate(gllist)  
  deallocate(glsort)  

  return

end subroutine put_gset
