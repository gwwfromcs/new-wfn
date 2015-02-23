
subroutine angular(p,gs,ffts, crys,bands,kpoints, wavefn, syms, pw_params,nrg)
!
  use all_to_all_module
  include 'use.h'
  implicit none             
  include 'interface.h'
  include 'all_to_all.h'

!  integer llmax
!  parameter (llmax = (2+1)**2+1) ! maximum lm quantum number
!
!     (c) 1996 Bernd Pfrommer
!     
!     Based on the decomp subroutine from Albertos code
!     I admit, this is a real mess, but I wont clean it up.
!
!     INPUT:
!     -----

  type(symmetry), intent(in) :: syms  
  type(crystal), intent(in) :: crys  
  type(fft_struc), intent(in) :: ffts         ! storage for FFT
  type(parallel_gspace), intent(in) :: gs(*)  ! gspaces where wavefn computed
  type(band), intent(in) :: bands             ! for nrk
  type(kpoint), intent(in) :: kpoints         ! for the weights kpoints%w
  type(complex_gspace_array), intent(in) :: wavefn  
  type(pw_parameter), intent(in) :: pw_params ! to know if timing should be done
  integer, intent(in) :: nrg                  ! number of radial grid points


!     OUTPUT:
!     ------
!
!     the symmetrized weights for the whole sphere, s, p, m=-1,0,1 and d
!     (first index). 
!
!     p index                   meaning
!     
!     1                         whole sphere
!     2                         s component
!     3                         px component
!     4                         py component
!     5                         pz component
!
!     6                         xy component
!     7                         yz component
!     8                         zx component
!     9                         x2-y2
!    10                         z2-r2


  real(dp), intent(out) ::    &
        p(llmax,crys%mxdatm,crys%ntype,bands%max,bands%nrk,crys%nspin)


!     DESCRIPTION:
!     -----------
!
!     computes the angular momentum decomposition with respect to
!     each site.
!
! CHANGES
! 
!  11/27/00 -  
!     included d components and corrected symmetrization of the decomposed 
!     DOS. Also got rid of symmetrization before squaring P component
!       from Seung-Hoon Jhi
!
!
!
!     ---------------------- local variables ----------------------------

  integer :: ndim, itype, ia, m, irk, is, n, l, i, idif, kk, nnx, &
       nny, nnz, ntot, icount, ir, idx, idy, idz, ic, b2, na2, nb
  complex(dp) :: psi, psq, epfi, epfip, y11, y1m1, dw(nrg)  

  real(dp), parameter :: eps = 1.0d-8, delta = 1.0d-5  
  real(dp) :: t0, cro(3), xdif, cs, cp, cd1, cd2, cd3, rcut, rcut2, rb(3), &
       ra(3), qmod, fi, ys, yx, yy, yz, yxy, yyz, yzx, yxy2, yz2, r(3), &
       qkcart(3), qr, sj0, sj1, sj2, prjc(10), pos(3), rc(3), dr(3), drl, &
       dintot, fip, gpk(3), gcart(3), gmag, gtr, cpc
  real(dp), external :: gimmetime  
  !
  !   dynamic arrays 
  !
  complex(dp), allocatable :: wpr(:,:,:,:), dhd(:), psisq(:,:,:), &
       sumtot(:,:)
! wprxyz(:,:,:,:,:,:)
  real(dp), allocatable :: psym(:,:,:,:,:,:), rga(:,:,:), wga(:,:,:), myslab(:)
  !
  !     ----- variables for the gspace loop ---------
  !
  integer :: igv(4), fftn(4), ffth(4), igv3, irod, iord, igs  
  real(dp) :: gv(3)  
  !
  !     ------------------------------------------------------------------
  !
  real(dp) & 
           vtran(48,3,3), &   !vector transformation matrices
           qtran(48,5,5)     !rank 2 tensor transformation matrices

!
!     -------------------------------------------------------------------
!
  t0=gimmetime()

  nnx=ffts%fftsize(1)
  nny=ffts%fftsize(2)
  nnz=ffts%fftsize(3)
  ntot = nnx*nny*nnz
  dintot = 1.d0/real(ntot,dp)
  !
  !     initial coordinates, where the slab of this processor starts
  !
  b2 = (nny*nnz)/gs(1)%nproc 
  na2 = mod(nny*nnz,gs(1)%nproc)
  icount=(b2*gs(1)%myproc+min(gs(1)%myproc,na2))*nnx


  write(9,200)
  write(9,100) nrg
!  write(9,105) 

  cs = done / sqrt(pi4)  
  cp = sqrt(dthree) * cs  
  cpc = sqrt(dthree / (8.0d0 * pi))  
  cd1 = sqrt(3.75d0) * cs  
  cd2 = dtwo * cd1  
  cd3 = sqrt(1.25d0) * cs  
  !
  !     compute the weights and abscissas for a Gauss-Legendre interpolation
  !
  !
  allocate(rga(nrg,crys%mxdatm,crys%ntype))
  allocate(wga(nrg,crys%mxdatm,crys%ntype))
  do itype=1, crys%ntype
     do ia=1, crys%natom(itype)
        call mygauleg(0,crys%mtsphere(itype),rga(1,ia,itype), &
                wga(1,ia,itype),nrg)

     end do
  end do

  allocate(wpr(nrg,9,bands%max,crys%nspin))
!     allocate(wprxyz(nrg,3,bands%max,crys%mxdatm,
!    $     crys%ntype,crys%nspin))
!     wprxyz=0
  allocate(psym(llmax-1,crys%mxdatm,crys%ntype,bands%max,bands%nrk,crys%nspin))
  psym=0
      
  allocate(sumtot(bands%max,crys%nspin))
  call findtransf(crys%bvec,crys%avec,syms%ntrans, syms%rsymmat,vtran,qtran)

  do irk = 1, bands%nrk     ! ------- loop over different kpoints -----
     !
     !     fourier transform the wave functions at this kpoint to
     !     realspace, square, and fourier transform back
     !
     nb=bands%nband(irk,1)
     if(crys%nspin.gt.1) nb=max(bands%nband(irk,1), bands%nband(irk,2))


     allocate(psisq(gs(irk)%length, nb, crys%nspin))
     psisq=0
     allocate(dhd(gs(irk)%r_size))

     do is = 1, crys%nspin                      ! -------- loop over spins --
        do n = 1, bands%nband(irk, is)            ! ------ loop over bands --
           call fourier_transform(-1, ffts, gs(irk), &
                wavefn%data((n-1)*gs(irk)%length+1, irk, is), dhd(1), 1)
           do i = 1, gs(irk)%r_size  
              dhd(i) = dhd(i) * conjg(dhd(i))  
           end do
           call fourier_transform(1, ffts, gs(irk), psisq(1, n, is), &
                dhd(1), 1)
        end do
     end do
     deallocate (dhd)  


     do itype = 1, crys%ntype          ! ----- loop over different types ----
        do ia = 1, crys%natom(itype)     ! ---- loop over different atom ----
           !
           !-------- now start the ugly loop over G-space -----------
           !
           wpr = zzero  

           sumtot = zzero  
           fftn(1:3) = gs(irk)%fftsize(1:3)  
           fftn(4) = gs(irk)%fftsize(3)   
           ffth(:) = fftn(:) / 2  
           ndim = gs(irk)%length  
           igs = 0  
           ! loop through x/y gspace
           do iord = 1, gs(irk)%lorder  
              irod = gs(irk)%order(1, iord)  
              igv(1) = irod / gs(irk)%fftsize (2)  
              igv(2) = mod(irod, gs(irk)%fftsize(2))  
              igv(3) = gs(irk)%order(2, iord)  
              igv(4) = gs(irk)%order(3, iord)  
              igv(:) = mod(igv(:) + ffth(:), fftn(:)) - ffth(:)  
              gv(1:2) = real(igv(1:2), dp)  
              gpk(1:2) = gv(1:2) + gs(irk)%rk(1:2)  
              ! loop over z axis
              do igv3 = igv(3), igv(4)  
                 gv(3) = real(igv3, dp)  
                 gpk(3) = gv(3) + gs(irk)%rk(3)  

                 igs = igs + 1  
                 qmod = sqrt(gs(irk)%ekin(igs))  
                 fip = dot_product(gpk, crys%rat(:, ia, itype))

                 epfip = exp(cmplx(0.d0, fip, dp))  
                 gcart = matmul(crys%bvec, gv)  
                 ! length |G|
                 gmag = sqrt(dot_product(gcart, gcart)) 

                 fi = dot_product(gv, crys%rat(:, ia, itype))  
                 if (gmag > delta) then  
                    ! |G|*|R|
                    gtr = crys%mtsphere(itype) * gmag  
                    epfi = exp(cmplx(0.d0, fi, dp)) * pi4 / gmag**3 * &
                         (sin(gtr) - gtr * cos(gtr))
                 else  
                    epfi = pi43 * crys%mtsphere(itype)**3  
                 end if

                 do is = 1, crys%nspin               ! -------- loop over spins
                    do n = 1, bands%nband(irk, is)      ! ------ loop over band
                       psq = psisq(igs, n, is) * epfi  
                       sumtot(n, is) = sumtot(n, is) + psisq(igs, n, is) * epfi 
                    end do

                 end do
                 !        Express k+G in cartesian coordinates and compute
                 !        the relevant spherical harmonics. Note that by
                 !        dividing by qmod we only need to multiply components
                 !        (See any book dealing with Ylms)
                 !
                 qkcart = matmul(crys%bvec, gpk)  
                 if (qmod > delta) qkcart = qkcart / qmod  
                 !
                 !                 s, x, y, z, xy, yz, zx, x2-y2, 3z2-r2
                 !
                 Ys = cs
                 !                 yx = cp*qkcart(1) ! this for px decomposition
                 !                 yy = cp*qkcart(2) ! this for py decomposition
                ! this for Y_1,1 
                 Yx = cp*qkcart(1) ! this for px decomposition
                 Yy = cp*qkcart(2) ! this for py decomposition

! ...  not necessary JSH ...
!                 y11 = -cpc * cmplx(qkcart(1), -qkcart(2), dp)  
                 ! this for Y_1,-1
!                 y1m1 = cpc * cmplx(qkcart(1), qkcart(2), dp)
                 ! Y_1,0

                 yz = cp * qkcart(3)  
                 yxy = cd2 * qkcart(1) * qkcart(2)  
                 yyz = cd2 * qkcart(2) * qkcart(3)  
                 yzx = cd2 * qkcart(3) * qkcart(1)  
                 yxy2 = cd1 * (qkcart(1) * qkcart(1) - qkcart(2) * qkcart(2))  
                 yz2 = cd3 * (dthree * qkcart(3) * qkcart(3) - done) 
                 !
                 !        Compute spherical Bessel functions at the radial grid
                 !
                 do m = 1, nrg               ! ---------- loop over radial grid
                    qr = qmod * rga(m, ia, itype)  
                    sj0 = done  
                    sj1 = dzero  
                    sj2 = dzero  
                    if (qr >= delta) then  
                       sj0 = sin(qr) / qr  
                       sj1 = (sj0 - cos(qr)) / qr  
                       sj2 = 3 * sj1 / qr - sj0  
                    end if
                    do is = 1, crys%nspin            ! -------- loop over spins
                       do n = 1, bands%nband(irk, is)  ! ------ loop over bands
                          psi = wavefn%data((n-1)*ndim+igs, irk, is)  

                          psi = psi * epfip  
                          wpr(m, 1, n, is) = wpr(m, 1, n, is) + psi*sj0*ys  
                          !         wpr(m,2,n,is) = wpr(m,2,n,is) + psi*sj1*yx
                          !         wpr(m,3,n,is) = wpr(m,3,n,is) + psi*sj1*yy
                          wpr(m,2,n,is) = wpr(m,2,n,is) + psi*sj1*Yx 
                          wpr(m,3,n,is) = wpr(m,3,n,is) + psi*sj1*Yy
                          wpr(m,4,n,is) = wpr(m,4,n,is) + psi*sj1*Yz

!                          wpr(m, 2, n, is) = wpr(m, 2, n, is) + psi*sj1*y1m1  
!                          wpr(m, 3, n, is) = wpr(m, 3, n, is) + psi*sj1*yz  
!                          wpr(m, 4, n, is) = wpr(m, 4, n, is) + psi*sj1*y11  

                          wpr(m, 5, n, is) = wpr(m, 5, n, is) + psi*sj2*yxy  
                          wpr(m, 6, n, is) = wpr(m, 6, n, is) + psi*sj2*yyz  
                          wpr(m, 7, n, is) = wpr(m, 7, n, is) + psi*sj2*yzx  
                          wpr(m, 8, n, is) = wpr(m, 8, n, is) + psi*sj2*yxy2  
                          wpr(m, 9, n, is) = wpr(m, 9, n, is) + psi*sj2*yz2  

                       end do
                    end do
                 end do                         ! ------- loop over radial grid

              end do
           end do                             ! ------- end of loop over Gspace
           call all_sum_all(wpr, nrg*9*bands%max*crys%nspin)  
           call all_sum_all(sumtot, bands%max*crys%nspin)  

!            wprxyz(1:nrg,1:3,:,ia,itype,:)=wpr(1:nrg,2:4,:,:)
           !
           !              square and integrate the radial wave function
           !              except for the px,py,pz quantum numbers
           !
           do is = 1, crys%nspin  
              do n = 1, bands%nband(irk, is)  
                 do l = 1, 9  
                    prjc(l) = dzero           
                    do m = 1, nrg                   ! loop over radial grid
                       prjc(l) = prjc(l) + wga(m, ia, itype) * &
                            (rga(m, ia, itype)*abs(wpr(m, l, n, is)))**2
                    end do
                    prjc(l) = prjc(l) * pi4**2 / crys%vcell  
                 end do
                 p(1, ia, itype, n, irk, is) = abs(sumtot(n, is)) / &
                      crys%vcell
!                     p(2:10,ia,itype,n,irk,is) = prjc(1:9)
                 p(2, ia, itype, n, irk, is) = prjc(1)  
                 p(3, ia, itype, n, irk, is) = prjc(2) + prjc(3) + prjc(4)  
                 p(7, ia, itype, n, irk, is) = prjc(5) + prjc(6) + &
                      prjc(7) + prjc(8) + prjc(9)
              end do
           end do
        end do                    ! -------- end loop over atoms of one type
     end do                            ! -------- end loop over atomic types
     deallocate (psisq)  
     !
     ! ...    symmetrize
     !
     do is = 1, crys%nspin  

        do n = 1, bands%nband(irk, is)  
           !               write(9,*) '------- band ',n,' -----------'
           ! --- loop over all symmetry operation
           do i = 1, syms%ntrans  
              !                  write(9,*) 'symmat:'
              !                  write(9,*) syms%rsymmat(:,:,i)
              do itype = 1, crys%ntype   
                do ia = 1, crys%natom (itype)  
                    !              write(9,*) '#### destination atom:',ia
                    !
                    !              -1
                    !     find mtrx    * (rat - tnp)
                    !
                    do l = 1, 3  
                       cro(l) = dzero  
                       do m = 1, 3  
                          cro(l) = cro(l) + real(syms%mtrx(m, l, i), dp) * &
                               (crys%rat (m, ia, itype) - syms%tnp(m, i))
                       end do
                    end do
                    !
                    !     search for matching atom
                    !
                    do l = 1, crys%natom(itype)  
                       do m = 1, 3  
                          xdif = abs(cro(m) - crys%rat(m, l, itype)) * &
                               ootwopi  
                          idif = xdif + eps  
                          if (abs(xdif - real(idif,dp)) > eps) goto 26  
                       end do
                       kk = l  
                       goto 27  
26                     continue  
                    end do
                    write(9,*) 'subroutine angular:', &
                         'unable to find equivalent atom for k=', ia
                    call mystop
27                  continue
  
                    !
                    ! ... sclar
                    !
                    psym(1,ia,itype,n,irk,is)= &
                        psym(1,ia,itype,n,irk,is)+p(2,kk,itype,n,irk,is)
                    !
                    ! ... vector 
                    !
                    do l=1,3
                          psym(l+1,ia,itype,n,irk,is)= &
                           psym(l+1,ia,itype,n,irk,is) &
                          +vtran(i,l,1)*p(3,kk,itype,n,irk,is) &
                          +vtran(i,l,2)*p(4,kk,itype,n,irk,is) &
                          +vtran(i,l,3)*p(5,kk,itype,n,irk,is)
                    enddo
                    !
                    ! ... 2nd rank tensor 
                    !
                    do l=1,5
                      psym(l+4,ia,itype,n,irk,is)= &
                         psym(l+4,ia,itype,n,irk,is) &
                          +qtran(i,l,1)*p(6,kk,itype,n,irk,is) &
                          +qtran(i,l,2)*p(7,kk,itype,n,irk,is) &
                          +qtran(i,l,3)*p(8,kk,itype,n,irk,is) &
                          +qtran(i,l,4)*p(9,kk,itype,n,irk,is) &
                          +qtran(i,l,5)*p(10,kk,itype,n,irk,is)
                    enddo

!                    !                write(9,*) '*** from source:',kk
!                    do l = 1, 3  
!                       dw = zzero  
!                       do m = 1, 3  
!                          dw(:) = dw(:) + syms%rsymmat(l, m, i) * &
!                               wprxyz(:, m, n, kk, itype, is)
!                       end do
!                       prjc(l) = dzero  
!                       do m = 1, nrg               ! loop over radial grid
!                          prjc(l) = prjc(l) + wga(m, ia, itype) * &
!                               (rga(m, ia, itype) * abs(dw(m)))**2
!                       end do
!                       prjc(l) = prjc(l) * pi4**2 / crys%vcell  
!                       ! write(9,*) 'l=',l-2,':',prjc(l)
!                       psym(3+l, ia, itype, n, irk, is) = &
!                            psym(3+l, ia, itype, n, irk, is) + prjc(l)

                 end do         ! -------- end loop over atoms of one type
              end do                 ! -------- end loop over atomic types
           end do             ! -------- end loop over symmetry operations

           !
           !              normalize
           !
           psym(:,:,:,n,irk,is)=psym(:,:,:,n,irk,is)/  dble(syms%ntrans)
           if(mod(is,2).eq.0) then
             write(9,60)irk,n
           else
             write(9,61)irk,n
           endif

           do itype=1, crys%ntype
             do ia=1, crys%natom(itype)
               write(9,51) ia,itype,p(1,ia,itype,n,irk,is), &
                   psym(:,ia,itype,n,irk,is)
             end do
           end do

        end do                              ! -------- end loop over bands
     end do                                 ! -------- end loop over spins
  end do                          ! -------- end of loop over kpoints ----

!  deallocate(wprxyz)  
  deallocate(sumtot)  
  deallocate(wpr)  
  
  deallocate(wga) ; deallocate(rga)


  p(2:llmax,:,:,:,:,:)=psym(1:llmax-1,:,:,:,:,:)
  deallocate(psym)
  !
  !-----------------Now compute the SLICES in  realspace ------------
  !

  if (gs(1)%myproc == 0) open(71, file = 'SLICES')  

  allocate(myslab(nnz))  
  do is = 1, crys%nspin  
     do irk = 1, bands%nrk  
        do n = 1, bands%nband(irk, is)  
           !
           !              fourier transform to realspace
           !
           ndim = gs(irk)%length  
           allocate(dhd(gs (irk)%r_size))  
           call fourier_transform(-1, ffts, gs(irk), &
                wavefn%data((n-1)*ndim+1, irk, is), dhd, 1)
           myslab = dzero
           do ir = 1, gs(irk)%r_size             ! loop over all points
              !
              !                 compute the current position vector rc
              !
              !                 x runs fastest, then y, then z
              !
              ! position in total grid
              ic = icount + ir - 1  
              !
              idz = (ic / (nnx * nny))  
              if ((idz >= nnz) .or. (idz < 0)) then  
                 write (9,*) 'BUG IN ANGULAR!'  
                 call myflush(9)  
              else  
                 myslab(idz + 1) = myslab(idz + 1) + &
                      abs(dhd(ir))**2  
              end if
           end do           ! loop over all points of this wave function
           deallocate(dhd)  
           !
           call all_sum_all(myslab, nnz)  
           if (gs(1)%myproc == 0) then  
              myslab = myslab * dintot  
              write(71, '(a,i4,1x,a,i4)') '#kpoint:', irk, 'band:', n  
              write(71, '(i5,g16.8)') (i, myslab(i) , i = 1, nnz)  
              write(71, '(''&'')')  
           end if

        end do                                         ! loop over bands
     end do                                          ! loop over kpoints
  end do                                               ! loop over spins
  deallocate(myslab)  
  if (gs(1)%myproc == 0) close(71)  
  if (iand(pw_params%output(1), 8) == 8) write(9, 110) gimmetime() - t0

  return  
50  format(i1,1x,i4,1x,i3,1x,i3,1x,i2,1x,7f9.5)  
 ! 50   format(5i4,f12.6)
 51   format(i3,1x,i3,1x,10f6.3)
 60   format('k-point :',i4,1x,'band :',i3,1x,'spin up')
 61   format('k-point :',i4,1x,'band :',i3,1x,'spin down')
100 format(' NUMBER OF RADIAL GRID POINTS:',i4)  
!105  format('#ATOM TYPE total    s     px     py     pz     xy     yz   &
!       zx     x2-y2     3z2-r2')
110 format(' TIME FOR ANGULAR MOMENTUM DECOMPOSITION:',f12.3)  
200 format(/' PERFORMING ANGULAR MOMENTUM DECOMPOSITION', &
              &     /1x,41('-'))

end subroutine angular
!
!     ==================================================================
!
subroutine findtransf(bvec,avec,ntrans,mtrx,vtran,qtran)

  include 'use.h'
  implicit none

!
!     from SNU code -JSH 
!     mtrx convention is different ..
!
!     OUTPUT
!     vtran(k,i,j) vector transformation matrix
!                  for the k-th symmetry operation
!     qtran(k,i,j) 2nd rank tensor transformation matrix
!                  for the k-th symmetry operation
!     avec(i,j)    i-th comp. of j-th primitive vector
!     bvec(i,j)    i-th comp. of j-th reciprocal vector
!
  real(dp) zero,um,dois,six
  parameter (zero=0.0D0, um=1.0D0, dois=2.0D0, six=6.0D0)
!
  real(dp) vtran(48,3,3),qtran(48,5,5),coef(5,3,3)
  real(dp) avec(3,3),bvec(3,3)
  real(dp) rt2i,rt6i,cjmn,delta
  real(dp) mtrx(3,3,48)
  integer ntrans
  integer i,j,k,m,n
!
  rt2i = um/sqrt(dois)
  rt6i = um/sqrt(six)
  delta=1.0D-7
!
!      compose the vector transformation matrix
!      find rotation matrices in real space
!
           
!          T      -1
!      ai=b *(2pi)
!
!      do i=1,3
! ai(:,i)=bvec(i,:)/(2*pi)
!      enddo
  vtran=0
  do i=1,3
    do j=1,3
      do k=1,ntrans
        vtran(k,i,j)=mtrx(i,j,k)
!     do m=1,3
!            do n=1,3
!              vtran(k,i,j) = vtran(k,i,j)
!    $         + bvec(m,i)*mtrx(m,n,k)*avec(n,j)
!            enddo
!            enddo
      enddo
    enddo
  enddo
!
!      compose the 2nd rank tensor transformation matrix
!
  coef = zero
  qtran = zero
  coef(1,1,2) = rt2i
  coef(1,2,1) = rt2i
  coef(2,2,3) = rt2i
  coef(2,3,2) = rt2i
  coef(3,1,3) = rt2i
  coef(3,3,1) = rt2i
  coef(4,1,1) = rt2i
  coef(4,2,2) = -rt2i
  coef(5,1,1) = -rt6i
  coef(5,2,2) = -rt6i
  coef(5,3,3) = dois*rt6i
  do i=1,5
   do j=1,5
    do k=1,ntrans
     do m=1,3
      do n=1,3
           cjmn = vtran(k,1,m)*(coef(j,1,1)*vtran(k,1,n) &
                             + coef(j,1,2)*vtran(k,2,n) &
                             + coef(j,1,3)*vtran(k,3,n)) &
               + vtran(k,2,m)*(coef(j,2,1)*vtran(k,1,n) &
                             + coef(j,2,2)*vtran(k,2,n) &
                             + coef(j,2,3)*vtran(k,3,n)) &
               + vtran(k,3,m)*(coef(j,3,1)*vtran(k,1,n) &
                             + coef(j,3,2)*vtran(k,2,n) &
                             + coef(j,3,3)*vtran(k,3,n)) 
           qtran(k,i,j) = qtran(k,i,j) + coef(i,m,n)*cjmn
      end do
     end do    
    end do
   end do
  end do
!
  do k=1,ntrans
     vtran(k,:,:) = abs(vtran(k,:,:))
     qtran(k,:,:) = abs(qtran(k,:,:))
  enddo
  return

end subroutine findtransf
!
!     ==================================================================
!

subroutine mygauleg(x1, x2, x, w, n)  

  use constants
  implicit none  
  !
  !     from Numerical Recipes
  !
  !     .. Parameters ..
  real(dp), parameter :: eps = 3.0d-14  
  !     ..
  !     .. Scalar Arguments ..
  real(dp), intent(in) :: x1, x2  
  integer, intent(in) :: n  
  !     ..
  !     .. Array Arguments ..
  real(dp), intent(out) :: w(n), x(n)  
  !     ..
  !     .. Local Scalars ..
  real(dp) :: p1, p2, p3, pp, xl, xm, z, z1, dj
  integer :: i, j, m  
  !     ..
  m = (n + 1) / 2  
  xm = 0.5d0 * (x2 + x1)  
  xl = 0.5d0 * (x2 - x1)  
  do i = 1, m  
     z = cos(pi * (real(i, dp) - 0.25d0) / (real(n, dp) + 0.5d0))  
     do
        p1 = done  
        p2 = dzero  
        do j = 1, n  
           p3 = p2  
           p2 = p1
           dj = real(j, dp)
           p1 = ((dtwo*dj - done)*z*p2 - (dj - done)*p3) / dj
        end do
        pp = real(n, dp) * (z * p1 - p2) / (z * z - done)  
        z1 = z  
        z = z1 - p1 / pp  
        if (abs(z - z1) <= eps) exit
     end do
     x(i) = xm - xl * z  
     x(n + 1 - i) = xm + xl * z  
     w(i) = dtwo * xl / ((done - z * z) * pp * pp)  
     w(n + 1 - i) = w(i)  
  end do

  return  

end subroutine mygauleg
