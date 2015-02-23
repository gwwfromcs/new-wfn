!*
subroutine generate_gspace(ioutput, gs, crys, symms)  

  use all_to_all_module
  include 'use.h'  
  implicit none            ! implicit? Just say no!
  include 'interface.h' 
  include 'all_to_all.h'
  !
  !      INPUT:
  !      ------
  !
  integer, intent(in) :: ioutput                               ! output flag
  type(crystal), intent(inout) :: crys  
  type(symmetry), intent(in) :: symms  
  !
  !      INPUT/OUTPUT:
  !      -------------
  !
  type(parallel_gspace), intent(inout) :: gs  
  !
  !     ------------------------------------------------------------------
  !
  !     Multiprocessor version of the traditional gspace routine.
  !     The gspace is now maintained differently: Each processor only
  !     holds the particular area of the g-space, which is assigned to it.
  !
  !     The layout of the data is dictated by the FFT: it is the
  !     only part which requires a substantial amount of communication
  !     between the processors:
  !
  !       * each processor holds a certain number of rods,
  !         which are lines along the z-axis of the FFT grid.
  !         To get good load balancing, the rods are handed out
  !         to processors in descending length:
  !
  !
  !         +
  !         ++
  !         +++
  !         ++++
  !         +++++
  !         ++++++
  !         012210....
  !
  !         This idea is due to Andrew Canning.
  !
  !       * each processor holds only the data on those gvectors
  !         that fall on one of the rods it holds. In other words, the
  !         ball of gvectors is pierced by the rods which point along the
  !         axis, and each processor holds those G-vectors that fall on th
  !         rods which it owns. The map between gvectors and rods is done
  !         in the array  gs%order. It is used to unfold the data
  !         when FFTs are performed, and is also necessary whenever a loop
  !         over gspace is performed.
  !
  !
  !     The routine to collect the gvectors into stars was written by JLM.
  !
  !     All the other stuff was written by  Bernd Pfrommer, UC Berkeley
  !
  !
  !
  !     ------------------------- actions --------------------------------
  !
  !     *  generates the gspace according to the variable
  !                gs%gmax
  !                gs%rk(3)  (the k-point)
  !                gs%nproc
  !                gs%myproc
  !     *  if the fourier grid dimensions gs%fftsize  are zero, they
  !        get automatically set to snugly hold the sphere
  !        corresponding FFT in the fftstruc
  !
  !     *  if
  !           gs%igvec    generate the gvectors (h,k,l)
  !           gs%iekin    compute the kinetic energy of the g-vectors
  !           gs%istar    generate the star information: inds, mstar, phas
  !           gs%iinverse find the inverse of each gvec (local gspace only
  !           gs%struc    compute the structure factors, and determine if
  !                       crystal has inversion symmetry
  !
  !     * the gs%order array is of great importance. It contains
  !       the mapping between the packed gspace and the FFT grid
  !
  !     gs%order(1,i)=irod   gives the rod number of the FFT rod owned by
  !     gs%myproc, on which the data i belongs. The numbering starts
  !     from 0. If you know the rod number, the fft grid position
  !     is given by:
  !
  !           nx = irod/gs%fftsize(2)
  !           ny = mod(irod, gs%fftsize(2))
  !
  !
  !     ^ y
  !     |
  !     ++++++++++++++++++++++++++++++++++
  !     |  |  |  |  |  |  |  |  |  |  |  |
  !     ++++++++++++++++++++++++++++++++++
  !     |  |  |  |  |  |  |  |  |  |  |  |
  !     ++++++++++++++++++++++++++++++++++
  !     |  |  |  |  |  |  |  |  |  |  |  |
  !     ++++++++++++++++++++++++++++++++++
  !     |  |  |  |  |  |  |  |  |  |  |  |
  !     ++++++++++++++++++++++++++++++++++
  !     |  |  |  |  |  |  |  |  |  |  |  |
  !     ++++++++++++++++++++++++++++++++++
  !     |  |  |  |  |  |  |  |  |  |  |  |
  !     G++G++++++++++++++++++++++++++++++
  !     |  |  |  |  |  |  |  |  |  |  |  |
  !     G++G++G+++++++++++++++++++++++++++
  !     |  |  |  |  |  |  |  |  |  |  |  |
  !     G++G++G++G++++++++++++++++++++++++
  !     |  |  |  |  |  |  |  |  |  |  |  |
  !     G++G++G++G++++++++++++++++++++++++ ---> x
  !
  !
  !
  !
  !     ------------------------ memory allocation -----------------------
  !
  !     this subroutine allocates memory, and attaches it to various
  !     structures:
  !
  !     * always
  !
  !     gs%order       (the order array)
  !     gs%fastfind    (a small index array for faster finding)
  !
  !
  !     * additionally, if
  !
  !           gs%igvec    ==>  gs%gvec
  !           gs%iekin    ==>  gs%ekin
  !           gs%istar    ==>  gs%inds, gs%phase, gs%mstar
  !           gs%iinverse ==>  gs%inverse
  !           gs%istruc   ==>  gs%struc
  !
  !     ------------------ local variables ------------------------------
  !
  real(dp) :: qk(3), gmod, gmax2, rgm, rkmod, fi, gmr, diff, t0
  complex(dp) :: ph, ex(48), str
  integer :: &
       ns, &            ! number of stars
       i, j, k, &       ! g-space integers
       l, ir, ih, &
       ncount, &        ! counts number of symmops leaving G invariant
       ktran(3, 48), &  ! symmetry transformed vector
       ngt, &           ! number of vectors in total gspace
       ng, &            ! number of vectors on myproc
       myblocks, kstart, nrods, &
       ntotrod, &       ! total number of rods everywhere
       ib, &            ! counts blocks
       ip, &            ! counts planes within one block
       iprot, iavail, itest, &
       kk(3), kmax(3), nx, ny, p, iproc, lorder, nt, ja, jmax, &
       igv(3), ili, &
       sflag            ! detect if inside sphere
  !
  !     variables for the generic gspace loop
  !
  integer :: il, fftn(4), ffth(4), iorder, irod, gv(4), gv3,icomp1  
  real(dp), parameter :: delta = 1.0d-6
  integer, allocatable :: gtvec(:,:), ihelp(:), rodmap(:,:), &
       rodlen(:), idxarr(:)
  real(dp), allocatable :: ekt(:)  
  logical :: imyrod, iflag
  integer, external :: findvec
  real(dp), external :: gimmetime  
  !     -----------------------------------------------------------------
  !     calculate length of k-vector

  gs%ntype = crys%ntype  
  t0 = gimmetime()  

  rkmod = sqrt(dot_product(gs%rk, matmul(crys%bdot, gs%rk)))  
  !     determine boundaries kmax such that grid holds the sphere snugly
  call findbound(gs%gmax + rkmod, crys%bvec, kmax(1), kmax(2), kmax(3))
  !
  !     Pick nice FFT grid if it has not already been picked
  !
  !
  if (gs%fftsize(1) <= 0 .or. gs%fftsize(2) <= 0 .or. gs%fftsize(3) <= 0) then
     gs%fftsize(1:3) = 2 * kmax(1:3) + 1  
     call adjustfft (gs%fftsize, symms)  
  end if
  gs%fftsize(4) = kmax(1)                                  ! right limit as is
  gs%fftsize(5) = mod(gs%fftsize(1) - kmax(1), gs%fftsize(1))  
  !     sanity checks: sphere should not be larger than grid limits:
  if (2 * kmax(1) + 1 > gs%fftsize(1) .or. &
       2 * kmax(2) + 1 > gs%fftsize(2) .or. 2 * kmax(3) + 1 > gs%fftsize(3)) &
       then
     write(9, *) ' sphere does not fit into FFT grid!'  
     write(9, *) ' gmax:', gs%gmax  
     write(9, *) ' nfft:', gs%fftsize  
     write(9, *) ' kmax:', kmax  
     write(9, *) ' presumably, the cutoff is too low, and |k| is ', &
          &'comparable to gmax'
     call mystop  
  end if
  if (iand(ioutput, 8) == 8) write(9,*)' adjustfft', gimmetime() - t0  
  call myflush(9)  
  !
  !     compute the size of the FFT grid that this processor holds
  !     in realspace
  !
  gs%r_size = (gs%fftsize(2) * gs%fftsize(3)) / gs%nproc  
  if (mod(gs%fftsize(2) * gs%fftsize(3), gs%nproc) > gs%myproc) &
       gs%r_size = gs%r_size + 1
  gs%r_size = gs%r_size * gs%fftsize(1)  
  gs%ig0 = -1                         ! will be set to >0 if myproc holds G=0
  gmax2 = gs%gmax * gs%gmax  
  !
  !     ------- count the number of elements in total and small gspace. --
  !                     set up the array with rod lengths
  !
  gs%lorder = 0  
  gs%length = 0  
  gs%totlength = 0  
  ntotrod = 0  
  allocate(rodlen(gs%fftsize(1) * gs%fftsize(2)))  
  allocate(rodmap(gs%fftsize(2), gs%fftsize(1)))  
  rodmap = 0  
  allocate(idxarr(gs%fftsize(1) * gs%fftsize(2)))  
  do i = -kmax(1), kmax(1)  
     qk(1) = real(i, dp) + gs%rk(1)  
     nx = mod(i + gs%fftsize(1), gs%fftsize(1))  
     do j = -kmax(2), kmax(2)  
        qk(2) = real(j,dp) + gs%rk(2)  
        ny = mod(j + gs%fftsize(2), gs%fftsize(2))  
        sflag = 0
        do k = -kmax(3), kmax(3)  
           qk(3) = real(k, dp) + gs%rk(3)  
           !              calculate norm with special metric
           gmod = (qk(1) * crys%bdot(1, 1) + qk(2) * crys%bdot(2, 1) + &
                qk(3) * crys%bdot(3, 1)) * qk(1) + &
                (qk(1) * crys%bdot(1, 2) + qk(2) * crys%bdot(2, 2) + &
                qk(3) * crys%bdot(3, 2)) * qk(2) + &
                (qk(1) * crys%bdot(1, 3) + qk(2) * crys%bdot(2, 3) + &
                qk(3) * crys%bdot(3, 3)) * qk(3)
           if (gmod <= gmax2) then                  ! we are inside the sphere
              if (sflag == 0) then             ! mark where we intruded sphere
                 !                    write(9,*) 'intruded sphere'
                 kstart = k  
                 sflag = 1  
                 ntotrod = ntotrod + 1  
                 rodmap(ny + 1, nx + 1) = ntotrod  
              end if
              gs%totlength = gs%totlength + 1  
           else                                    ! we are outside the sphere
              if (sflag == 1) then                   ! we have left the sphere
                 sflag = 0  
                 rodlen(ntotrod) = k - kstart  
              end if
           end if
        end do                                      ! loop along 3rd dimension
        if (sflag == 1) then                        ! reached border of sphere
           sflag = 0  
           rodlen(ntotrod) = k - kstart  
        end if
     end do      ! 2nd dim
  end do         ! 1st dim
  !
  !     sort array with rod lengths in ascending order
  !

  if (iand(ioutput, 8) == 8) write(9,*)' compute size', gimmetime() - t0  
  call myflush(9)  

  do i = 1, ntotrod  
     idxarr(i) = i  
  end do
  call iheapsort (ntotrod, rodlen, idxarr)  
  !
  !     figure out how many rods this processor holds, and allocate array
  !
  nrods = ntotrod / gs%nproc  
  if (mod(ntotrod / gs%nproc, 2) == 0) then  
     if (mod(ntotrod, gs%nproc) > gs%myproc) nrods = nrods + 1  
  else  
     if (mod(ntotrod, gs%nproc) >= gs%nproc - gs%myproc) nrods = nrods + 1
  end if
  gs%lorder = nrods  
  !      write(9,*) 'total number of rods:',gs%lorder
  !
  !     assign rods to processors. modify rodmap.
  !
  !     this scales like N^4/3, so watch out for large systems
  !
  gs%length = 0  
  do i = 0, ntotrod - 1  
     if ((mod(i / gs%nproc, 2) == 0 .and. mod(i, gs%nproc) == gs%myproc) .or. &
          (mod(i / gs%nproc, 2) /= 0 .and. mod(i, gs%nproc) == &
          (gs%nproc - gs%myproc - 1))) then
        gs%length = gs%length + rodlen(idxarr(ntotrod - i))  
        !            write(9,*) 'hold:', i+1,rodlen(idxarr(ntotrod-i))
        do j = 1, gs%fftsize(1)  
           do k = 1, gs%fftsize(2)  
              if (rodmap(k, j) == idxarr(ntotrod - i)) then  
                 rodmap(k, j) = -1             ! mark as belonging to this proc
                 !                     write(9,*) 'marked rod as mine'
              end if
           end do
        end do
     end if
  end do  
  !
  !     allocate arrays
  !
  allocate(gs%order(3, gs%lorder))  
  allocate(gs%fastfind(3, gs%fftsize(1)))  
  gs%fastfind(:,:) = -1                                           ! invalidate
  if (gs%iekin) then
     allocate(gs%ekin(gs%length))  
  else
     nullify(gs%ekin)
  end if
  if (gs%igvec) then
     allocate(gs%gvec(3, gs%length))  
  else
     nullify(gs%gvec)
  end if
  if (gs%istruc) then
     allocate(gs%struc(gs%ntype, gs%length))  
  else
     nullify(gs%struc)
  end if
  if (gs%istar) then  
     !        permanent ones, attached to gspace structure
     allocate(gs%inds(gs%length))
     allocate(gs%phase(gs%length)) 
     allocate(gs%mstar(gs%totlength))  
     !        volatile work arrays
     allocate(gtvec(3, gs%totlength))  
     allocate(ihelp(gs%totlength))  
     allocate(ekt(gs%totlength))  
  else  
     nullify(gs%mstar)
     nullify(gs%inds)
     nullify(gs%phase)
  end if

  if (iand(ioutput, 8) == 8) write(9,*)' allocate', gimmetime() - t0  
  call myflush(9)  
  !
  !     ---------------- now really generate the gspace--------------
  !
  if (gs%istruc) crys%icomplex = .false.  
  icomp1=0
  lorder = 0  
  ng = 0  
  ngt = 0  
  do i = -kmax(1), kmax(1)  
     qk(1) = real(i, dp) + gs%rk(1)  
     nx = mod(i + gs%fftsize(1), gs%fftsize(1))  
     do j = -kmax(2), kmax(2)  
        qk(2) = real(j, dp) + gs%rk(2)  
        ny = mod(j + gs%fftsize(2), gs%fftsize(2))  
        sflag = 0  
        !            write(9,*) nx,ny,'rodmap:', rodmap(ny+1,nx+1)
        imyrod = .false.  
        if (rodmap(ny + 1, nx + 1) < 0) then     ! aha, it is me who holds it
           imyrod = .true.                       ! set flag
           if (gs%fastfind(1, nx + 1) < 0) then  ! started new xy plane
              gs%fastfind(1, nx + 1) = lorder + 1  
              gs%fastfind(3, nx + 1) = ng  
           end if
           gs%fastfind(2, nx + 1) = lorder + 1   ! catch end of xy plane
        end if
        do k = -kmax(3), kmax(3)  
           qk(3) = real(k, dp) + gs%rk(3)  
           !              calculate norm with special metric
           gmod = (qk(1) * crys%bdot(1, 1) + qk(2) * crys%bdot(2, 1) + &
                qk(3) * crys%bdot(3, 1)) * qk(1) + &
                (qk(1) * crys%bdot(1, 2) + qk(2) * crys%bdot(2, 2) + &
                qk(3) * crys%bdot(3, 2)) * qk(2) + &
                (qk(1) * crys%bdot(1, 3) + qk(2) * crys%bdot(2, 3) + &
                qk(3) * crys%bdot(3, 3)) * qk(3)
           if (gmod <= gmax2) then  !     -------- we are inside the sphere
              if (sflag == 0) then  ! mark where we intruded sphere
                 !                    write(9,*) 'intruded sphere'
                 sflag = 1  
                 if (imyrod) then  
                    lorder = lorder + 1  
                    gs%order(1, lorder) = nx * gs%fftsize(2) + ny  
                    gs%order(2, lorder) = mod(k + gs%fftsize(3), &
                         gs%fftsize (3) )
                 end if
              end if
              if (imyrod) then  
                 ng = ng + 1  
                 if (i == 0 .and. j == 0 .and. k == 0) gs%ig0 = ng  
                 ! i hold the G=0 vector
              end if
              if (gs%istar) ngt = ngt + 1  
              if (gs%igvec .and. imyrod) then  
                 gs%gvec(1, ng) = i  
                 gs%gvec(2, ng) = j  
                 gs%gvec(3, ng) = k  
                 !                     write(9,*) ng, nx,ny,gs%gvec(:,ng)
              end if
              if (gs%istar) then  
                 gtvec(1, ngt) = i  
                 gtvec(2, ngt) = j  
                 gtvec(3, ngt) = k  
                 ekt(ngt) = gmod  
                 if (imyrod) then  
                    ihelp(ngt) = ng  
                 else  
                    ihelp(ngt) = 0  
                 end if
              end if
              if (gs%iekin .and. imyrod) then  
                 gs%ekin(ng) = gmod  
              end if
              if (gs%istruc .and. imyrod) then  
                 do nt = 1, crys%ntype  ! loop over atomic types
                    str = zzero
                    jmax = crys%natom(nt)  
                    do ja = 1, jmax  
                       fi = real(i, dp) * crys%rat(1, ja, nt) + &
                            real(j, dp) * crys%rat(2, ja, nt) + &
                            real(k, dp) * crys%rat(3, ja, nt)
                       str = str + cmplx(cos(fi), -sin(fi), dp)  
                    end do
                    if (abs(aimag(str)) > 1.0d-6)  icomp1 =1  !crys%icomplex = .true.  
                    gs%struc(nt, ng) = str  
                 end do
              end if
           else  !     -------- we are outside the sphere
              if (sflag == 1) then  
                 if (imyrod) then  
                    gs%order(3, lorder) = mod(k - 1 + gs%fftsize(3), &
                         gs%fftsize(3))
                 end if
                 !                    write(9,*) 'left sphere'
                 sflag = 0  
              end if
           end if
        end do    ! loop along 3rd dimension
        if (sflag == 1) then                       ! didn't leave the sphere
           sflag = 0 
           if (imyrod) gs%order(3, lorder) = mod(k - 1 + gs%fftsize(3), &
                gs%fftsize(3))
        end if
     end do     ! 2nd dim
  end do        ! 1st dim

  if (iand(ioutput, 8) == 8) write(9,*)' generate', gimmetime() - t0  
  call myflush(9)  

  call all_sum_all(icomp1)
  if (icomp1 .gt. 0) crys%icomplex = .true.  
  if (iand(ioutput, 8) == 8) write(9,*)' all sum', gimmetime() - t0  
  call myflush(9)  

  !     ------------------------------------------------------------------
  !                          set up stars if required
  !
  if (gs%istar) then  
     !
     !        sort the g-vectors according to length and collect
     !        heapsort algorithm from numerical recipes,
     !        original coding by JLM

     if (ngt > 1) then  
        l = ngt / 2 + 1  
        ir = ngt  
20      continue  
        if (l > 1) then  
           l = l - 1  
           rgm = ekt(l)  
           kk(1:3) = gtvec(1:3, l)  
           ih = ihelp(l)  
        else  
           rgm = ekt(ir)  
           kk(1:3) = gtvec(1:3, ir)  
           ih = ihelp(ir)  
           ekt(ir) = ekt(1)  
           gtvec(1:3, ir) = gtvec(1:3, 1)  
           ihelp(ir) = ihelp(1)  
           ir = ir - 1  
           if (ir == 1) then  
              ekt(1) = rgm  
              gtvec(1:3, 1) = kk(1:3)  
              ihelp(1) = ih  
              !
              goto 25  
              !
           end if
        end if
        i = l  
        j = l + l  
21      if (j <= ir) then  
           if (j < ir) then  
              if (ekt(j) < ekt(j + 1)) j = j + 1  
           end if
           if (rgm < ekt(j)) then  
              ekt(i) = ekt(j)  
              gtvec(1:3, i) = gtvec(1:3, j)  
              ihelp(i) = ihelp(j)  
              i = j  
              j = j + j  
           else  
              j = ir + 1  
           end if
           goto 21  
        end if
        ekt(i) = rgm  
        gtvec(1:3, i) = kk(1:3)  
        ihelp(i) = ih  
        goto 20  
25      continue  
        ! end of if(ngt.gt.1)...
     end if

  if (iand(ioutput, 8) == 8) write(9,*)' stars 1', gimmetime() - t0  
  call myflush(9)  

     !
     !        Collect G-vectors into stars and compute "phase" factors
     !                   First star is (0,0,0)
     !
     ns = 1  
     gs%mstar(1) = 1  
     if (ihelp(1) > 0) then  
        gs%inds(ihelp(1)) = 1  
        gs%phase(ihelp(1)) = zone  
     end if
     !
     iprot = 1  
     iavail = 2  
     itest = 2  
     !
30   continue  
     !
     if (iavail > ngt) goto 50  
     !
     if (itest > ngt) then  
        diff = done  
     else  
        diff = abs(ekt(itest) - ekt(iprot))  
     end if
     if (diff > delta) then  
        !
        !        new star. calculate all gvectors in star
        !
        iprot = iavail  
        iavail = iprot + 1  
        itest = iprot + 1  
        ns = ns + 1  
        gs%mstar(ns) = 1  
        do i = 1, symms%ntrans  
           ktran(1, i) = symms%mtrx(1, 1, i) * gtvec(1, iprot) + &
                symms%mtrx(1, 2, i) * gtvec(2, iprot) + &
                symms%mtrx(1, 3, i) * gtvec(3, iprot)
           ktran(2, i) = symms%mtrx(2, 1, i) * gtvec(1, iprot) + &
                symms%mtrx(2, 2, i) * gtvec(2, iprot) + &
                symms%mtrx(2, 3, i) * gtvec(3, iprot)
           ktran(3, i) = symms%mtrx(3, 1, i) * gtvec(1, iprot) + &
                symms%mtrx(3, 2, i) * gtvec(2, iprot) + &
                symms%mtrx(3, 3, i) * gtvec(3, iprot)
        end do
        do i = 1, symms%ntrans  
           fi = symms%tnp(1, i) * real(ktran(1, i), dp) + &
                symms%tnp(2, i) * real(ktran(2, i), dp) + &
                symms%tnp(3, i) * real(ktran(3, i), dp)
           ex(i) = cmplx(cos(fi), sin(fi), dp)  
        end do
        ph = zzero  
        ncount = 0  
        do i = 1, symms%ntrans  
           if (ktran(1, i) - gtvec(1, iprot) == 0 .and. &
                ktran(2, i) - gtvec(2, iprot) == 0 .and. &
                ktran(3, i) - gtvec(3, iprot) == 0) then
              ph = ph + ex(i)  
              ncount = ncount + 1  
           end if
        end do
        if (ihelp(iprot) > 0) then  
           gs%phase(ihelp(iprot)) = ph / real(ncount, dp)  
           gs%inds(ihelp(iprot)) = ns  
        end if
        !
        goto 30  
        !
     end if
     !
     iflag = .false.  
     ncount = 0  
     ph = zzero  
     !
     !        calculate phase factor for gvec(itest)
     !
     do i = 1, symms%ntrans  
        if (ktran(1, i) - gtvec(1, itest) == 0 .and. &
             ktran(2, i) - gtvec(2, itest) == 0 .and. &
             ktran(3, i) - gtvec(3, itest) == 0) then
           iflag = .true.  
           ph = ph + ex(i)  
           ncount = ncount + 1  
        end if
     end do
     ! gvec(itest) is related to prototype
     if (iflag) then  
        gs%mstar(ns) = gs%mstar(ns) + 1  
        if (itest /= iavail) then  
           kk(1) = gtvec(1, iavail)  
           kk(2) = gtvec(2, iavail)  
           kk(3) = gtvec(3, iavail)  
           gtvec(1, iavail) = gtvec(1, itest)  
           gtvec(2, iavail) = gtvec(2, itest)  
           gtvec(3, iavail) = gtvec(3, itest)  
           gtvec(1, itest) = kk(1)  
           gtvec(2, itest) = kk(2)  
           gtvec(3, itest) = kk(3)  
           gmr = ekt(iavail)  
           ekt(iavail) = ekt(itest)  
           ekt(itest) = gmr  
           ih = ihelp(iavail)  
           ihelp(iavail) = ihelp(itest)  
           ihelp(itest) = ih  
        end if
        if (ihelp(iavail) > 0) then  
           gs%phase(ihelp(iavail)) = ph / real(ncount, dp)  
           gs%inds(ihelp(iavail)) = ns  
           !               write(9,*) gtvec(:,iavail), ph/dble(ncount)
        end if
        iavail = iavail + 1  
     end if
     itest = itest + 1  
     !
     goto 30  
50   continue  

     gs%nstar = ns  

  end if
  if (gs%istar) then  
     deallocate(gtvec)  
     deallocate(ihelp)  
     deallocate(ekt)  
  end if

  if (iand(ioutput, 8) == 8) write(9,*)' stars ', gimmetime() - t0  
  call myflush(9)  

  !
  !     find gspace inversion information if requested
  !
  if (gs%iinverse) then  
     if (gs%nproc > 1) then  
        write(9, *) 'inverse array only for LOCAL gspaces!'  
        call mystop  
     end if
     
     allocate(gs%inverse(gs%length))  
     il = 0  
     fftn(1:3) = gs%fftsize(1:3)  
     fftn(4) = gs%fftsize(3)  
     ffth(:) = fftn(:) / 2  
     do iorder = 1, gs%lorder                      ! loop through x/y gspace
        irod = gs%order(1, iorder)  
        gv(1) = irod / gs%fftsize(2)  
        gv(2) = mod(irod, gs%fftsize(2))  
        gv(3) = gs%order(2, iorder)  
        gv(4) = gs%order(3, iorder)  
        gv(:) = mod(gv(:) + ffth(:), fftn(:)) - ffth(:)  
        do gv3 = gv(3), gv(4)                             ! loop over z axis
           il = il + 1  
           igv(1) = -gv(1)  
           igv(2) = -gv(2)  
           igv(3) = -gv3  
           ili = findvec(igv, gs)  
           if (ili <= 0) then  
              write(9, *) 'generate_gspace: cannot find inverse!'  
              call mystop  
           end if
           gs%inverse(il) = ili  
        end do
     end do
  else  
     nullify(gs%inverse)  
  end if
  !
  !     set up the packing information arrays for the parallel 3dfft
  !

  if (iand(ioutput, 8) == 8) write(9,*)' inverses ', gimmetime() - t0  
  call myflush(9)  

  call setup_packinfo(gs)   

  if (iand(ioutput, 8) == 8) write(9,*)' setup pack ', gimmetime() - t0  
  call myflush(9)  

  call print_gspace(0, gs)  
  deallocate(rodlen)  
  deallocate(rodmap)  
  deallocate(idxarr)  

  if (iand(ioutput, 8) == 8) write(9,*)' end ', gimmetime() - t0  
  call myflush(9)  

  if (iand(ioutput, 8) == 8) write(9, 940) gimmetime() - t0  
  call myflush(9)  

  return   

940 format(' TIME FOR GSPACE SETUP:',f12.3)  

end subroutine generate_gspace
