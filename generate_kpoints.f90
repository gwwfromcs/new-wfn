!
subroutine generate_kpoints(ipr, kp, crys, symms, bands, iflag, &
     qshift, nsafety)

  include 'use.h'  
  implicit none  
  include 'interface.h'  
  !
  !     INPUT:
  !     -----
  !
  integer, intent(in) :: ipr        ! print flag
  integer, intent(in) :: nsafety    ! number of extra states to start diag with
  integer, intent(in) :: iflag  
  real(dp), intent(in) :: qshift(3)  
  type(crystal), intent(in) :: crys  
  type(symmetry)  :: symms  
  !
  !     OUTPUT:
  !     ------
  !
  type(kpoint), intent(out) :: kp  
  type(band), intent(out) :: bands  
  !
  !     computes k-points for integration
  !     over the brillouin zone. periodic mesh.
  !     adapted from sverre froyen plane wave program
  !     written june 1 1987. jlm
  !     modified july 31 1987. jlm
  !
  !     fortran 90 version 1996 Bernd Pfrommer
  !
  !     DESCRIPTION:
  !     -----------
  !     * computes kpoints in irreducible zone
  !     * computes kmap, which maps reducible to irreducible kpoints
  !
  !     ----------------- local variables ------------------------------
  !
  integer, parameter :: kptfile = 13
  real(dp), allocatable :: xk(:,:), ekk(:)  
  real(dp) :: fk2(3, 100000)  
  integer :: nfk2, indsubq(48), ntranq, nrk2, indrq2(100000),identity(3,3)  
  integer :: nx, ny, nz, ntrans, neg, km, kl, i, j, k, l, n, maxn, &
       nrk, jp, irk, im, jm, ip, jcar, jabs, is
  real(dp) :: wsum, &                                        ! sum of weights
       sx, sy, sz, ek1, ek, ek2, dx, dy, dz, dw, &
       diff, rkt(3), rktran(3), rkmod(3, 4), rkcar(3), rktran_inv(3),rktran_save(3)
  real(dp) :: reduce_01, x_dummy, qnorm
  logical :: lshift  
  integer :: nband1, ifmax1  
  real(dp), parameter :: eps = 1.0d-6
  !  
  !     -----------------------------------------------------------------
  !
  reduce_01(x_dummy) = mod(x_dummy + 5.0d2, done)

  ntrans = symms%ntrans  
  lshift = .true.  
  qnorm = sqrt(qshift(1)**2 + qshift(2)**2 + qshift(3)**2)


  if (qnorm < eps) lshift = .false.  
  !
  !     switch off reduction
  !
  if (.not. kp%reduced) ntrans = 1  
  if (.not. kp%reduced) symms%mtrans=0

  sx = kp%shift(1)  
  sy = kp%shift(2)  
  sz = kp%shift(3)  
  nx = kp%grid(1)  
  ny = kp%grid(2)  
  nz = kp%grid(3)  
  nrk = kp%nrk  
  !
  maxn = nx * ny * nz  
  nullify(kp%kmap)  

  if (kp%generated) then
     !     ------- count the number of irreducible k-points -------------
     !     allocate dynamic arrays

     allocate(xk(3, maxn))  
     allocate(kp%kmap(maxn))  
     !
     !     set up uniform array of k points
     !     kmap is used to mark reducible k points and also to
     !     map reducible to irreducible k points
     !
     dx = done / real(nx, dp)  
     dy = done / real(ny, dp)  
     dz = done / real(nz, dp)  
     n = 0  
     do i = 1, nx  
        do j = 1, ny  
           do k = 1, nz  
              n = n + 1  
              xk(1, n) = (real(i - 1, dp) + sx) * dx  
              xk(2, n) = (real(j - 1, dp) + sy) * dy  
              xk(3, n) = (real(k - 1, dp) + sz) * dz  
              kp%kmap(n) = n  
           end do
        end do
     end do
     !
     !     reduce to irreducible zone
     !
     nrk = 0  
     dw = done / real(n, dp)  
     do i = 1, n  
        if (kp%kmap(i) == i) then  
           !
           !     new irreducible point
           !     mark with negative kmap
           !
           nrk = nrk + 1  
           rkt(1) = xk(1, i)  
           rkt(2) = xk(2, i)  
           rkt(3) = xk(3, i)  
           kp%kmap(i) = -nrk  
           if (i /= n) then  
              !
              !     operate on irreducible rk with the symmetry operations
              !
              do j = 1, ntrans  

! restore rotations with fractional translation
!
!                 if((crys%ilsdau.ne.1) &
!                     .or.(SUM(abs(symms%tnp(:,j))).lt.1.d-6)) then

                 rktran(1) = real(symms%mtrx(1, 1, j), dp) * rkt(1) + &
                      real(symms%mtrx(1, 2, j), dp) * rkt(2) + &
                      real(symms%mtrx(1, 3, j), dp) * rkt (3)
                 rktran(2) = real(symms%mtrx(2, 1, j), dp) * rkt(1) + &
                      real(symms%mtrx(2, 2, j), dp) * rkt(2) + &
                      real(symms%mtrx(2, 3, j), dp) * rkt(3)
                 rktran(3) = real(symms%mtrx(3, 1, j), dp) * rkt(1) + &
                      real(symms%mtrx(3, 2, j), dp) * rkt(2) + &
                      real(symms%mtrx(3, 3, j), dp) * rkt(3)
                 !
                 !     translate to interval 0-1.
                 !
                 do k = 1, 3  
                    rktran(k) = reduce_01(rktran(k))  
                    rktran_inv(k) = reduce_01(-rktran(k))  
                 end do
                 !
                 ! remove (mark) k points related to irreducible rk by symmetry
                 !
                 ip = i + 1  
                 do k = ip, n  
                    if (kp%kmap(k) == k) then  
                       !
                       !     both the transformed rk ...
                       !
                       do l = 1, 3  
                          if (abs(rktran(l) - xk(l, k)) > eps) goto 930  
                       end do
                       kp%kmap(k) = nrk  
                       goto 940  
                       !
                       !     ... and its inverse (time reversal)
                       !
930                    continue  
                       if (iand(iflag, 1048576) == 1048576) goto 940  
                       do l = 1, 3  
                          diff = abs(rktran_inv(l) - xk(l, k))  
                          if ((.not. kp%reduced) .or. diff > eps) goto 940  
                       end do
                       kp%kmap(k) = nrk  
940                    continue  
                    end if
                 end do
!              end if
              end do
           end if
        end if
     end do
  end if
  !     -------- read kpt ------------------------------------------------
  if (.not. kp%generated) then  
 

     !
     !     read in number of eigenvalues to be computed, the lowest and
     !     highest filled band, the integration weight and the reduced
     !     k-vector.
     !
     open(unit = kptfile, file = 'KPOINTS', status = 'unknown', &
          form = 'formatted')
     read(kptfile, *) nrk  
     kp%nrk = nrk  
     call create_kpoint(kp, nrk)     
     ! give 5 bands safety
     call create_band(bands, nrk, max(bands%min(1),bands%min(2)) + nsafety, &
                  crys%nspin)
     allocate(ekk(nrk))  
     do i = 1, nrk  
        read(kptfile, *) kp%label(i), bands%nband(i, 1), k, &
             bands%ifmax(i, 1), kp%w(i), (kp%rk(j, i), j = 1, 3)
        !            read(kptfile,9000) kp%label(i), bands%nband(i,1),
        !     $           k, bands%ifmax(i,1), kp%w(i), (kp%rk(j,i),j=1,3)
        ! 9000       format(f8.3,3i4,f7.2,3f16.8)
     end do
     bands%nband = bands%nband + nsafety  
     bands%ifmax(:, bands%nspin) = bands%ifmax(:, 1)  
     bands%nband(:, bands%nspin) = bands%nband(:, 1)  

     write(9, 100) nrk  

     close(unit = kptfile)  
     ! calculate length of vectors
     do i = 1, nrk  
        ekk(i) = dzero  
        do j = 1, 3  
           do k = 1, 3  
              ekk(i) = ekk(i) + kp%rk(j, i) * crys%bdot(j, k) * kp%rk(k, i)
           end do
        end do
     end do
     !
     !     --------------------- generate kpts ---------------------------
     !
  else  
     !
     !        allocate space for kpoints and other work arrays
     !
     call create_kpoint(kp, nrk)  
     call create_band(bands, nrk, max(bands%min(1),bands%min(2)) + nsafety, & ! give 5 bands safety
          crys%nspin)  
  
     allocate(ekk(nrk))  
     !
     !     the points will be computed from parameters
     !     nx, ny, and nz are the number of points in the three
     !     directions dermined by the lattice wave vectors. sx, sy, and
     !     sz shifts the grid of integration points from the origin.
     !
     maxn = nx * ny * nz  

     allocate(kp%rk_fbz(3,maxn))
     allocate(kp%ind_symtr(maxn))
     kp%nrk_fbz=maxn

     kp%ind_symtr(:)=1   ! identity
     allocate(kp%ind_G0(3,maxn))
     kp%ind_G0=0.d0

     !
     !     set up uniform array of k points
     !     kmap is used to mark reducible k points and also to
     !     map reducible to irreducible k points
     !
     dx = done / real(nx, dp)  
     dy = done / real(ny, dp)  
     dz = done / real(nz, dp)  
     n = 0  
     do i = 1, nx  
        do j = 1, ny  
           do k = 1, nz  
              n = n + 1  
              xk(1, n) = (real(i - 1, dp) + sx) * dx  
              xk(2, n) = (real(j - 1, dp) + sy) * dy  
              xk(3, n) = (real(k - 1, dp) + sz) * dz  
              kp%kmap(n) = n  
           end do
        end do
     end do

     kp%rk_fbz(:,:)=xk(:,:)
     !
     !      reduce to irreducible zone
     !
     nrk = 0  
     dw = done / real(n, dp)  
     do i = 1, n  
        if (kp%kmap(i) == i) then  
           !
           !     new irreducible point
           !     mark with negative kmap
           !
           nrk = nrk + 1  
           kp%rk(1, nrk) = xk(1, i)  
           kp%rk(2, nrk) = xk(2, i)  
           kp%rk(3, nrk) = xk(3, i)  
           kp%kmap(i) = -nrk  
           kp%w(nrk) = dw  
           if (i /= n) then  
              !
              !     operate on irreducible rk with the symmetry operations
              !
              do j = 1, ntrans  

!                 if((crys%ilsdau.ne.1) &
!                     .or.(SUM(abs(symms%tnp(:,j))).lt.1.d-6)) then

                 rktran(1) = real(symms%mtrx(1, 1, j), dp) * kp%rk(1, nrk) + &
                      real(symms%mtrx(1, 2, j), dp) * kp%rk(2, nrk) + &
                      real(symms%mtrx(1, 3, j), dp) * kp%rk(3, nrk)
                 rktran(2) = real(symms%mtrx(2, 1, j), dp) * kp%rk(1, nrk) + &
                      real(symms%mtrx(2, 2, j), dp) * kp%rk(2, nrk) + &
                      real(symms%mtrx(2, 3, j), dp) * kp%rk(3, nrk)
                 rktran(3) = real(symms%mtrx(3, 1, j), dp) * kp%rk(1, nrk) + &
                      real(symms%mtrx(3, 2, j), dp) * kp%rk(2, nrk) + &
                      real(symms%mtrx(3, 3, j), dp) * kp%rk(3, nrk)

                 rktran_save=rktran

                 !
                 !     translate to interval 0-1.
                 !


                 do k = 1, 3  
                    rktran(k) = reduce_01(rktran(k))  
                    rktran_inv(k) = reduce_01(-rktran(k))  
                 end do
                 !
                 ! remove (mark) k points related to irreducible rk by symmetry
                 !
                 ip = i + 1  
                 do k = ip, n  
                    if (kp%kmap(k) == k) then  
                       !
                       !     both the transformed rk ...
                       !
                       do l = 1, 3  
                          if (abs(rktran(l) - xk(l, k)) > eps) goto 130
                       end do
                       kp%w(nrk) = kp%w(nrk) + dw  
                       kp%kmap(k) = nrk  
                       kp%ind_symtr(k)=j
                       kp%ind_G0(:,k)=xk(:,k)-rktran_save(:)

                       goto 140  
                       !
                       !     ... and its inverse (time reversal)
                       !
130                    continue  
                       if (iand(iflag, 1048576) == 1048576) goto 140
                       do l = 1, 3  
                          diff = abs(rktran_inv(l) - xk(l, k))  
                          if ((.not. kp%reduced) .or. diff > eps) goto 140
                       end do
                       kp%w(nrk) = kp%w(nrk) + dw  
                       kp%kmap(k) = nrk  
                       kp%ind_symtr(k)=-j

!?
!                       kp%ind_G0(:,k)=xk(:,k)-rktran_save(:)
                       kp%ind_G0(:,k)=xk(:,k)+rktran_save(:)

140                    continue  
                    end if
                 end do
!                 end if
              end do
           end if
        end if
     end do


! count how many symmetry operations are used to relate k-points in
! FBZ to those in IBZ

     symms%mtrans=0

     do i = 1, maxn

        if(abs(kp%ind_symtr(i)).ne.1) then
        do j=1,symms%mtrans
           if(abs(kp%ind_symtr(i)).eq.abs(symms%ind_mtrans(j))) goto 111
        end do

        symms%mtrans=symms%mtrans+1
        symms%ind_mtrans(symms%mtrans)=abs(kp%ind_symtr(i))
        symms%ind_ktrans(abs(kp%ind_symtr(i)))=symms%mtrans

  111   continue
        end if
      end do

! now keep track of inverse symmetry operations

      identity=0
      do i=1,3
         identity(i,i)=1
      end do

      do i=1,symms%mtrans
         k=abs(symms%ind_mtrans(i))

         do j=1,symms%ntrans

            if(SUM(ABS(MATMUL(symms%mtrx(:,:,k),symms%mtrx(:,:,j))-identity)).eq.0) symms%ind_inv(i)=j

         end do
      end do


     !
     !     set up array indk
     !     indk(1-6,i) gives the index for irreducible points that
     !     are neighbors of point i
     !
     n = 0  
     do i = 1, nx  
        do j = 1, ny  
           do k = 1, nz  
              n = n + 1  
              if (kp%kmap(n) <= 0) then  
                 irk = -kp%kmap(n)  
                 ip = i + 1  
                 jp = j + 1  
                 kl = k + 1  
                 if (i == nx) ip = 1  
                 if (j == ny) jp = 1  
                 if (k == nz) kl = 1  
                 im = i - 1  
                 jm = j - 1  
                 km = k - 1  
                 if (i == 1) im = nx  
                 if (j == 1) jm = ny  
                 if (k == 1) km = nz  
                 ip = ((ip - 1) * ny + j - 1) * nz + k  
                 kp%ind(1, irk) = abs(kp%kmap(ip))  
                 im = ((im - 1) * ny + j - 1) * nz + k  
                 kp%ind(2, irk) = abs(kp%kmap(im))  
                 jp = ((i - 1) * ny + jp - 1) * nz + k  
                 kp%ind(3, irk) = abs(kp%kmap(jp))  
                 jm = ((i - 1) * ny + jm - 1) * nz + k  
                 kp%ind(4, irk) = abs(kp%kmap(jm))  
                 kl = ((i - 1) * ny + j - 1) * nz + kl  
                 kp%ind(5, irk) = abs(kp%kmap(kl))  
                 km = ((i - 1) * ny + j - 1) * nz + km  
                 kp%ind(6, irk) = abs(kp%kmap(km))  
              end if
           end do
        end do
     end do
     wsum = dzero  
     do is=1, bands%nspin
       do i = 1, nrk  
!         bands%nband(i, bands%nspin) = bands%max  
!         bands%ifmax(i, bands%nspin) = bands%min  

         bands%nband(i, is) = bands%max  
         bands%ifmax(i, is) = bands%min(is)  
         wsum = wsum + kp%w(i)  
       end do
     end do
     !
     !     find rk**2, etc por the printout
     !
     write(9, 101) nrk, nx, ny, nz, sx, sy, sz, bands%min(1)  
     write(9, 102)  
     !
     do i = 1, nrk  
        do j = 1, 4  
           do k = 1, 3  
              rkmod(k, j) = kp%rk(k, i)  
           end do
        end do
        rkmod(1, 2) = rkmod(1, 2) - done  
        rkmod(2, 3) = rkmod(2, 3) - done
        rkmod(3, 4) = rkmod(3, 4) - done
        ek = 2.0d3
        do j = 1, 4  
           ek1 = dzero  
           ek2 = dzero  
           do k = 1, 3  
              do l = 1, 3  
                 ek1 = ek1 + rkmod(k, j) * crys%bdot(k, l) * rkmod(l, j)  
                 ek2 = ek2 + (rkmod(k, j) - done) * crys%bdot(k, l) * &
                      (rkmod(l, j) - done)
              end do
           end do
           if (ek1 < ek) then  
              ek = ek1  
              jcar = j  
           end if
           if (ek2 < ek) then  
              ek = ek2  
              jcar = -j  
           end if
        end do
        !
        !        rk in cartesian coordinates
        !
        neg = 0  
        do k = 1, 3  
           jabs = abs(jcar)  
           rkcar(k) = crys%bvec(k, 1) * rkmod(1, jabs) + &
                crys%bvec(k, 2) * rkmod(2, jabs) + &
                crys%bvec(k, 3) * rkmod(3, jabs)
           if (jcar < 0) rkcar(k) = rkcar(k) - crys%bvec(1, k) - &
                crys%bvec(2, k) - crys%bvec(3, k)
           if (rkcar(k) < dzero) neg = neg + 1  
        end do
        if (neg <= 1) neg = 1  
        if (neg > 1) neg = -1  
        !
        !     printout
        !

! if too many k points, I guess nobody cares to look at the
! details

        if(nrk.lt.200)  then 
        write(9, 103) i, kp%w(i), (kp%rk(j, i), j = 1, 3), &
             (real(neg, dp) * rkcar(j), j = 1, 3), ek
        end if
     end do
     !
     if (abs(wsum) < eps) wsum = dzero  
     !     ------ calculate |k|^2 for all kpoints
     do i = 1, nrk  
        ekk(i) = dzero
        do j = 1, 3  
           do k = 1, 3  
              ekk(i) = ekk(i) + kp%rk(j, i) * crys%bdot(j, k) * kp%rk(k, i)
           end do
        end do
     end do
     deallocate(xk)  
     !     -------- write kpt -----------------
     if (kp%generated) then
        !
        !  write out the number of eigenvalues to be computed, the lowest
        !  and highest filled band, the integration weight and the reduced
        !  k-vector.
        !
        open(unit = kptfile, file = 'SCF_KPOINTS', form = 'formatted')
        open(unit = 111, file = 'kpoints.dat', form = 'formatted')
        rewind(kptfile)
        write(kptfile, *) nrk
        write(111, *) nrk
        do i = 1, nrk
        !   write(kptfile, '(F8.3,X,3(I3,X),F6.2,3(X,F15.8))') &
           write(kptfile, '(F5.2,X,3(I3,X),F12.9,3(X,F13.8))') &
        !   write(kptfile, *) &
                kp%label(i), bands%min(1), 1, &
                bands%min(1), kp%w(i), (kp%rk(j, i), j = 1, 3)
           write(111, '(3f20.14,2i3)') (kp%rk(j, i), j = 1, 3),1,0
        end do
        close(unit = kptfile)   
     end if            ! for write KPOINTS file
  end if               ! for generate kpoints
  deallocate(ekk)  

  call myflush(9)  

  if (lshift) then    !xav
     call fullbzx(nrk, kp%rk, nfk2, fk2, symms)  
     call subgrpqx(qshift, ntranq, indsubq, symms)  
     call irrbzqx(nfk2, fk2, ntranq, indsubq, nrk2, indrq2, symms)  
     deallocate(kp%rk)  
     deallocate(kp%w)  
     deallocate(kp%label)  
     deallocate(kp%ind)  
     allocate(kp%rk(3, nrk2))  
     allocate(kp%w(nrk2))  
     allocate(kp%label(nrk2))  
     allocate(kp%ind(6, nrk2))  

     kp%label = 0
     nband1 = bands%nband(1, 1)  
     !      ifmin1 = bands%ifmin(1,1)
     !      ifmin1 = 1

     ifmax1 = bands%ifmax(1, 1)  
     !      deallocate(bands%ifmin)

     deallocate(bands%nband)  
     deallocate(bands%occup)  
     deallocate(bands%energy)  
     deallocate(bands%ekn)  
     deallocate(bands%ifmax)  
     deallocate(bands%gwout_index)

     !      bands%occup = 0.0
     !      bands%energy = 0.0
     !      bands%ekn = 0.0
     !      kp%ind = 0
     !      kp%label = 0
     allocate(bands%nband(nrk2, bands%nspin))  
     !      allocate(kp%ifmin(nrk2,crys%nspin))
     allocate(bands%ifmax(nrk2, bands%nspin))  
     allocate(bands%ekn(bands%max, nrk2, bands%nspin))  
     allocate(bands%occup(bands%max, nrk2, bands%nspin))  
     allocate(bands%energy(bands%max, nrk2, bands%nspin))  

     allocate(bands%gwout_index(bands%max, nrk2, bands%nspin))  
     bands%occup = dzero
     bands%energy = dzero 
     bands%ekn = dzero
     kp%ind = 0  
     kp%label = 0  
     kp%shift = dzero
     kp%grid = 0  
     kp%nrk = nrk2       !xav
     bands%nrk = nrk2  
     kp%generated = .false.  
     write(9, 104) nrk2  
     do i = 1, nrk2  
        kp%rk(1, i) = fk2(1, indrq2(i)) + qshift(1)  
        kp%rk(2, i) = fk2(2, indrq2(i)) + qshift(2)  
        kp%rk(3, i) = fk2(3, indrq2(i)) + qshift(3)  
        kp%w(i) = done / real(nrk2, dp)  
        bands%nband(i, 1) = nband1  
        !        bands%ifmin(i,1) = 1
        bands%ifmax(i, 1) = bands%min(1)  
        !        bands%ifmax(i,1) = 4
     end do
!     bands%ifmax(:, bands%nspin) = bands%ifmax(:, 1)  
!
! what is this following doing?

     bands%ifmax(:, bands%nspin) = bands%min(bands%nspin)
     bands%nband(:, bands%nspin) = bands%nband(:, 1)  
  end if
  !      call qualitykp(crys, symms, kp)

  return
  
100 FORMAT(//' K POINTS READ FROM FILE: ',i5)  
101 FORMAT(//' K POINT GENERATION:', &
       &     /' ------------------', &
       &     /1X,I5,' K POINTS GENERATED BY PROGRAM FROM', &
       &     ' PARAMETERS :',/,4X,'N =',3I5,6X,'S =',3F6.2,6X, &
       &     'NB =',I5,/)
102 FORMAT(4X,'I',7X,'W',17X,'RK',22X,'RKCAR',12X,'EKIN',/)  
103 FORMAT(1X,I4,3X,F11.7,3X,3F8.4,3X,3F7.3,3X,F7.3)  

104 FORMAT(/,1X,I5,' K POINTS GENERATED FOR', &
       &     ' THE GW SLIGHTLY SHIFTED GRID')
120 format(' READ ',i5,' K-POINTS FORM FILE KPOINTS')  

end subroutine generate_kpoints
!
!
!
! ADDITIONAL SERVICE ROUTINES, X. Blase, E. Chang
!
!-------------------------------------
subroutine subgrpqx(qi, ntranq, indsub, symms)  
  !
  include 'use.h'  
  implicit none
  !
  !       include 'symmetry.h'
  !       include 'use.h'
  type(symmetry), intent(in) :: symms  
  !
  !......comunication scalar and arrays
  integer, intent(out) :: ntranq, indsub(48)  
  real(dp), intent(in) :: qi(3)  
  !
  !......local scalars
  integer :: itran, i, j
  real(dp) :: diff
  !......local arrays
  real(dp) :: qk(3), g0(3)  
  !
  ! we look for operations such that  R(q) = q + g0q
  !
  ntranq = 0  
  do itran = 1, symms%ntrans  
     !
     do i = 1, 3  
        qk(i) = dzero  
        do j = 1, 3  
           qk(i) = qk(i) + symms%mtrx(i, j, itran) * qi(j)  
        end do
     end do
     !
     call bring_to_01(qk, g0)  
     !
     diff = abs(qk(1) - qi(1)) + abs(qk(2) - qi(2)) + abs(qk(3) - qi(3))
     !
     if (diff > 1.0d-6) cycle  
     !
     ! store index element of subroup
     !
     ntranq = ntranq + 1  
     indsub(ntranq) = itran  
     !
  end do
  !
  ! print out
  !
  !      write(9,9000) ntranq,(indsub(i),i=1,ntranq)
  !9000  format(/,i3,' operations in small group of q:',48i4)
  !
  return  

end subroutine subgrpqx

!cccccccccccccccccccccccccccccccccccccccccccc

subroutine fullbzx(nrk, rk, nfk, fk, symms)  
  !
  include 'use.h'  
  implicit none
  integer, parameter :: nfkp = 100000  
  !       include 'bands.i'
  !
  !      SUBROUTINE TO DETERMINE POINTS IN FBZ FROM
  !      K-POINTS IN IRR BZ
  !
  integer, intent(in) :: nrk
  real(dp), intent(in) :: rk(3, 100000)  
  integer, intent(out) :: nfk
  real(dp), intent(out) :: fk(3, nfkp)
  !
  !       include 'symmetry.h'
  !       include 'use.h'
  type(symmetry), intent(in) :: symms  
  !
  !
  !
  !  O   NFK,FK             K-POINTS IN FULL BZ
  !                         ONE POINT IN THE SET FK
  !
  !.....local scalars
  integer :: irk, it, i, j, iflag, ifk
  !.....local array
  real(dp) :: qk(3)  
  integer :: mati(3, 3)  
  real(dp) :: xmat(3, 3), g0(3)
  real(dp), parameter :: eps = 1.0d-5
  !
  !      LOOP OVER Q-POINTS
  !
  nfk = 0  
  do irk = 1, nrk  
     !
     !      LOOP OVER TRANSFORMATIONS
     !
     do it = 1, symms%ntrans  
        !
        !      ROTATE RK
        !
        do i = 1, 3  
           qk(i) = dzero  
           do j = 1, 3  
              qk(i) = qk(i) + symms%mtrx(i, j, it) * rk(j, irk)  
           end do
        end do
        !
        !      TRANSLATE TO INTERVAL 0-1
        !
        call bring_to_01(qk, g0)  
        !
        !      COMPARE TO OTHER POINTS IN FULL ZONE
        !
        if (nfk == 0) goto 40  
        iflag = 0  
        do ifk = 1, nfk  
           if (abs(qk(1) - fk(1, ifk)) < eps .and. &
                abs(qk(2) - fk(2, ifk)) < eps .and. &
                abs(qk(3) - fk(3, ifk)) < eps) iflag = 1
        end do
        if (iflag == 1) cycle
        !
        !      STORE NEW KPOINT IN FBZ
        !
40      nfk = nfk + 1  
        !
        if (nfk > nfkp) then  
           write(9, *) 'nfk ,nfkp=', nfk, nfkp  
           call mystop(210)  
        end if
        fk(1, nfk) = qk(1)  
        fk(2, nfk) = qk(2)  
        fk(3, nfk) = qk(3)  
     end do
  end do
  !
  !      WRITE(9,70) NFK
  !70    FORMAT(5X,4HNFK=,I5)

  return
  
end subroutine fullbzx
!---------------------------------------------------------
subroutine irrbzqx(nfk, fk, ntranq, indsub, nrq, indrq, symms)  
  !
  include 'use.h'  
  implicit none
  !
  integer, parameter :: nfkp = 100000  
  !
  !      SUBROUTINE REDUCES FULL BZ (FK) TO IRREDUCIBLE PART (RQ)
  !      USING OPERATIONS IN THE SUBGROUP OF KN
  !
  real(dp), intent(in) :: fk(3, *)
  integer, intent(out) :: indrq(*)
  integer, intent(in) :: ntranq, nfk, indsub(48)
  integer, intent(out) :: nrq
  !
  !       include 'symmetry.h'
  !       include 'use.h'
  type(symmetry), intent(in) :: symms  
  !
  !
  !
  !  I   NFK,FK             K-POINTS IN FULL BZ
  !  O   NRQ,INDRQ          Q-POINTS IN IRR BZ WITH RESPECT TO SUBGROUP,
  !                         INDEX OF RQ STORED IN INDRQ REFERING TO FK
  !
  !      UTILITY ARRAYS
  !
  real(dp) :: qk(3), g0(3), xmati(3, 3)
  integer :: ik, it, i, j, k, l, irq, iflag
  real(dp) :: det, matinvert
  !
  !      INITIALIZE NUMBER OF POINTS IN IRR. BZ
  !
  nrq = 0  
  !
  !
  !      LOOP OVER K-POINTS IN FULL ZONE
  !
  do ik = 1, nfk
     if (ik == 1) goto 60  
     !
     !      LOOP OVER TRANSFORMATION MATRICES
     !
     do it = 1, ntranq  
        !
        do k = 1, 3  
           do l = 1, 3  
              xmati(k, l) = real(symms%mtrx(k, l, indsub(it)), dp)  
           end do
        end do
        det=matinvert(xmati)
        !
        do i = 1, 3  
           qk(i) = dzero  
           do j = 1, 3  
              qk(i) = qk(i) + xmati(i, j) * fk(j, ik)  
           end do
        end do
        !
        call bring_to_01(qk, g0)  
        !
        !      COMPARE TO OTHER K-POINTS IN THE IRR. BZ WITH RESPECT TO QVEC
        !
        do irq = 1, nrq
           iflag = 0  
           do i = 1, 3  
              if (abs(fk(i, indrq(irq)) - qk(i)) > 1.0d-6) iflag = 1
           end do
           if (iflag == 0) goto 70  
        end do
     end do
     !
     !      END LOOP OVER TRANSFORMATION MATRICES
     !
     !      STORE INDEX OF NEW K-POINT IN IRR. BZ WITH RESPECT TO Q
     !
60   nrq = nrq + 1  
     indrq(nrq) = ik  
70 end do
  !
  !      END LOOP OVER FULL BZ
  !
  ! print out
  !
  !      write(9,9000) NRQ
  !9000  format(/,' irr k-points according to subgroup of q ',i5)
  !      do IRQ=1,NRQ
  !      write(9,9010) INDRQ(IRQ), (FK(I,INDRQ(IRQ)),I=1,3)
  !9010  format(i5,3x,3f8.5)
  !      enddo
  !
  return
  
end subroutine irrbzqx
!---------------------------------------------------------
subroutine bring_to_01(x, g)  
  !
  use constants
  implicit none  
  real(dp), intent(inout) :: x(3), g(3)  
  real(dp) :: one_eps, y(3)
  !
  integer :: i  
  !
  one_eps = done - 1.0d-6  
  !
  do i = 1, 3  
     y(i) = mod(x(i) + 5.0d2, done)  
     if (y(i) > one_eps) y(i) = dzero  
     g(i) = x(i) - y(i)  
     g(i) = -g(i)  
     x(i) = y(i)  
  end do

  return  

end subroutine bring_to_01


