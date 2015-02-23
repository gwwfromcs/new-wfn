!     @process extchk
!
subroutine forstresssym(ipr, mode, fsum, ssum, icorr, force, &
     stress, ntype, natom, rat, avec, bvec, bdot, vcell, ntrans, mtrx, &
     tnp, rsymmat, mxdatm)
  !
  !     symmetrizes and computes the hellman-feynman forces and stress
  !     adapted from sverre froyen plane wave program
  !     written january 18 1988. jlmo
  !     explicit typing and other cosmetics 1995 Bernd Pfrommer
  !
  use constants
  use lsdau_shared

  implicit none            ! implicit? Just say no!
  !
  !     INPUT:
  !     -----
  integer, intent(in) :: &
       ipr, &              ! print flag
       natom(*), &         ! number of atoms of each type
       mtrx(3, 3, 48), &   ! symmetry operation: matrices
       mode, &             ! if mode=1 then symmetrize stress/force at the end
       ntype, &            ! number of atomic types
       ntrans, &           ! number of symmetry operations
       mxdatm              ! array dimension
  real(dp), intent(in) :: &
       force(3, mxdatm, ntype), &  ! forces as from ewald subroutine (lattice)
       stress(6), &             ! stress tensor as from ewald subroutine (ryd)
       bdot(3, 3), &                              ! metric in reciprocal space
       bvec(3, 3), &                              ! reciprocal lattice vectors
       avec(3, 3), &                               ! realspace lattice vectors
       vcell, &
       tnp(3, 48), &                               ! nonprimitive translations
       rsymmat(3, 3, 48), &            ! the realspace point symmetry matrices
       rat(3, mxdatm, ntype)   ! they are needed for the stress-symmetrization
  character(len=2), intent(in) :: &
       icorr                             ! string for the exchange-correlation
  !
  !     OUTPUT:  fsum and ssum will return symmetrized force/stress
  !     ------
  !
  real(dp), intent(inout) :: &
       fsum(3, mxdatm, ntype), &                        ! unsymmetrized forces
       ssum(6)                                          ! unsymmetrized stress
  !
  !     ------------------------- local variables ------------------------
  !
  real(dp), parameter :: eps = 1.0d-6
  real(dp) :: fcart(3), select(6), ssym(6), cro(3), fac, xdif, rotfor(3),deltaf
  real(dp),allocatable::fcart_all(:,:)
  !
  !     WORK ARRAYS:
  !     -----------
  !
  real(dp) :: fsym(3, mxdatm, ntype)  
  !
  integer :: i, j, k, natomi, natomj, kk, l, m, n, idif,nn  
  !
  !      symmetrize
  !
  n=0
  do i = 1, ntype  
     natomi = natom(i)  
     do j = 1, natomi  
        n=n+1
        do k = 1, 3  
           fsym(k, j, i) = dzero  
        end do
     end do
  end do

  allocate(fcart_all(3,n))

  ssym = dzero  
  !
  !     symmetrizes forces
  !
  do i = 1, ntrans  
     do j = 1, ntype  
        natomj = natom(j)  
        do k = 1, natomj  
           !
           !                     -1
           !            find mtrx    * (rat - tnp)
           !
           do l = 1, 3  
              cro(l) = dzero  
              do m = 1, 3  
                 cro(l) = cro(l) + real(mtrx(m, l, i), dp) * &
                      (rat(m, k, j) - tnp(m, i))
              end do
           end do
           do l = 1, natomj  
              do m = 1, 3  
                 xdif = abs(cro(m) - rat(m, l, j)) / pi2  
                 idif = int(xdif + eps)  
                 if (abs(xdif - real(idif, dp)) > eps) goto 16  
              end do
              kk = l  
              goto 17  
16            continue  
           end do
           write(9, *) 'subroutine forstresssym:'  
           write(9, *) 'unable to find equivalent atom for k=', k  
           stop  
17         continue  
           !
           !              rotate force and add
           !
           do l = 1, 3  
              do m = 1, 3  
                 fsym(l, k, j) = fsym(l, k, j) + real(mtrx(l, m, i), dp) * &
                      fsum(m, kk, j)
              end do
           end do
        end do
     end do
  end do
  !
  !      symmetrize stress
  !
  do i = 1, 6  
     j = i  
     k = i  
     if (i > 3) j = i - 3  
     if (i > 3) k = j + 1  
     if (k > 3) k = 1  
     do l = 1, ntrans  
        ssym(i) = ssym(i) + &
             ssum(1) * real(mtrx(j, 1, l) * mtrx(k, 1, l), dp) + &
             ssum(2) * real(mtrx(j, 2, l) * mtrx(k, 2, l), dp) + &
             ssum(3) * real(mtrx(j, 3, l) * mtrx(k, 3, l), dp) + &
             ssum(4) * real(mtrx(j, 1, l) * mtrx(k, 2, l), dp) + &
             ssum(5) * real(mtrx(j, 2, l) * mtrx(k, 3, l), dp) + &
             ssum(6) * real(mtrx(j, 3, l) * mtrx(k, 1, l), dp) + &
             ssum(4) * real(mtrx(j, 2, l) * mtrx(k, 1, l), dp) + &
             ssum(5) * real(mtrx(j, 3, l) * mtrx(k, 2, l), dp) + &
             ssum(6) * real(mtrx(j, 1, l) * mtrx(k, 3, l), dp)
     end do
     ssym(i) = ssym(i) / real(ntrans, dp)  
  end do
  !
  !     transform to cartesian coordinates and printout
  !     forces
  !
  if (ipr /= 0) then  
     write(9, 100)  
  end if
  n = 0  
  do i = 1, ntype  
     natomi = natom(i)  
     do j = 1, natomi  
        n = n + 1  
        !
        !           compute force in cartesian coordinates: f_cart = a*f_q
        !           where a has the lattice vectors as columns
        !
        do k = 1, 3  
           fcart(k) = avec(k, 1) * force(1, j, i) + &
                avec(k, 2) * force(2, j, i) + &
                avec(k, 3) * force(3, j, i) + &
                (bvec(k, 1) * fsym(1, j, i) + &
                bvec(k, 2) * fsym(2, j, i) + &
                bvec(k, 3) * fsym(3, j, i)) / real(ntrans, dp)
        end do



       
        if (ipr /= 0) write(9, 101) n, (matmul(avec, rat(:, j, i)) / pi2), &
             (fcart(k), k = 1, 3)
        fcart_all(:,n)=fcart(:)
        !
        !      compute the dE/dq, which is NOT the force in lattice coordinates
        !      but  dE/dq = aT* f_cart
        !
        !          f_sum = aT * f_cart
        !
        do k = 1, 3  
           fsum(k, j, i) = dzero  
           do l = 1, 3  
              fsum(k, j, i) = fsum(k, j, i) + avec(l, k) * fcart(l)  
           end do
        end do
     end do
  end do

!------------------------------------------------------------
!                     !!!!!!CAUTION!!!!!!!
!------------------------------------------------------------
! There could be bugs in the code. The error in calculated
! forces contributed by the Hubbard U term is larger
! than expected.  However, it could well be due to numerical error.
! After days of debugging, I have to leave it this way.

! PZ
        if(iLSDAU.eq.1) then

        write(9,*)
        write(9,*) "LSDA+U force correction - as calculated"
        write(9,*)

  n=0
  nn=0
  do i = 1, ntype  
     natomi = natom(i)  
     do j = 1, natomi  
        n=n+1
        if(lproj_tem(i).ge.0) then
        nn=nn+1
        fcart_all(:,n)= fcart_all(:,n)+force_U(nn,:)
        end if

        if (ipr /= 0) write(9, 101) n, (matmul(avec, rat(:, j, i)) / pi2), fcart_all(:,n)

     end do
  end do


! If the hubbard contribution to the forces has a large numerical
! error, this could be a quick fix. I am afraid that the errors
! come from bugs instead of numerical. 
!

!        do i=1,3
!           deltaf=SUM(fcart_all(i,:))/n
!           fcart_all(i,:)=fcart_all(i,:)-deltaf
!        end do
! 
!        write(9,*)
!        write(9,*) "LSDA+U force correction - shifted"
!        write(9,*)
!
!------------------------ 
!  deallocate(lproj_tem)
!  deallocate(force_U)
!------------------------ 

!  n=0
!  do i = 1, ntype  
!     natomi = natom(i)  
!     do j = 1, natomi  
!        n=n+1
!        if (ipr /= 0) write(9, 101) n, (matmul(avec, rat(:, j, i)) / pi2), fcart_all(:,n)
!     end do
!  end do
!        write(9,*)
!        write(9,*) "NOTE: LSDA+U stress correction has not bee implemented"
!        write(9,*)

       end if
!
!------------------------------------------------------------------------
  !
  !     stress:      select(j,k) = b_j * ssym * b_k
  !
  fac = 14710.80d0 / vcell  
    ! j,k are corresponding indices in full
  do i = 1, 6  
     j = i  
     k = i  
     if (i > 3) j = i - 3  
     if (i > 3) k = j + 1  
     if (k > 3) k = 1  
     select(i) = bvec(j, 1) * ssym(1) * bvec(k, 1) + &
          bvec(j, 2) * ssym(2) * bvec(k, 2) + &
          bvec(j, 3) * ssym(3) * bvec(k, 3) + &
          bvec(j, 1) * ssym(4) * bvec(k, 2) + &
          bvec(j, 2) * ssym(5) * bvec(k, 3) + &
          bvec(j, 3) * ssym(6) * bvec(k, 1) + &
          bvec(k, 1) * ssym(4) * bvec(j, 2) + &
          bvec(k, 2) * ssym(5) * bvec(j, 3) + &
          bvec(k, 3) * ssym(6) * bvec(j, 1)
  end do
  !     add electronic stress to ewald contribution
  do i = 1, 6  
     ssum(i) = stress(i) + select(i)  
  end do
  !      write(9,*) 'forces on atoms in dE/dq'
  !      do i=1,ntype
  !         natomi = natom(i)
  !         do j=1,natomi
  !            write(9,*)  fsum(1,j,i), fsum(2,j,i), fsum(3,j,i)
  !         end do
  !      end do
  !
  !
  !     explicitly symmetrize stress and forces if desired.
  !
  !     this has to be done to avoid small, but significant errors:
  !     the lattice vectors might not be perfectly symmetric. this
  !     in turn produces an error in the stress, which is used to
  !     move the lattice vectors. after about 10 iterations, the error
  !     gets large enough to cause trouble.
  !
  if (mode >= 1) then  
     !         write(9,*) 'mtrx:'
     !         do l=1,ntrans
     !            write(9,'(i3/,3(3i5/))') l,
     !     $           mtrx(1,1,l),mtrx(2,1,l),mtrx(3,1,l),
     !     $           mtrx(1,2,l),mtrx(2,2,l),mtrx(3,2,l),
     !     $           mtrx(1,3,l),mtrx(2,3,l),mtrx(3,3,l)
     !         end do
     !
     !        first the stress. This is the simple part. Ideally:
     !                         -1
     !             sigma = rmat    *  sigma   * rmat
     !
     do i = 1, 6  
        j = i  
        k = i  
        if (i > 3) j = i - 3  
        if (i > 3) k = j + 1  
        if (k > 3) k = 1  
        ssym(i) = dzero  
        do l = 1, ntrans              !  loop through all symmetry operations
           !
           !                 do sigma = rmat * sigma *rmat
           !
           ssym(i) = ssym(i) + &
                ssum(1) * rsymmat(1, j, l) * rsymmat(1, k, l) + &
                ssum(2) * rsymmat(2, j, l) * rsymmat(2, k, l) + &
                ssum(3) * rsymmat(3, j, l) * rsymmat(3, k, l) + &
                ssum(4) * rsymmat(1, j, l) * rsymmat(2, k, l) + &
                ssum(5) * rsymmat(2, j, l) * rsymmat(3, k, l) + &
                ssum(6) * rsymmat(3, j, l) * rsymmat(1, k, l) + &
                ssum(4) * rsymmat(2, j, l) * rsymmat(1, k, l) + &
                ssum(5) * rsymmat(3, j, l) * rsymmat(2, k, l) + &
                ssum(6) * rsymmat(1, j, l) * rsymmat(3, k, l)
        end do
        ssym(i) = ssym(i) / real(ntrans, dp)  
     end do
     do i = 1, 6  
        ssum(i) = ssym(i)  
     end do
     !
     !        now symmetrize the forces again.
     !
     !        Assume in realspace, a point-symmetry rotation matrix
     !        looks like:
     !
     !          alpha = (cos(phi)  -sin(phi)    0)
     !                  (sin(phi)   cos(phi)    0)
     !                  (0          0           1)
     !
     !      *  the forces f_q in lattice coordinates transform like
     !
     !          f_q' = beta * f_q
     !
     !          where  beta := a^-1 * alpha * a
     !
     !      *  the quantity F_q:= dE/dq is the projection of the force
     !         onto the lattice vectors, and transforms differently:
     !
     !                  T
     !          F_q := a * a  *  f_q
     !
     !                          T
     !          F_q' = (beta^-1)  * F_q
     !
     !          now:       beta     = rmtrx
     !                     beta^-1  = mtrx
     !
     !
     do i = 1, ntype  
        natomi = natom(i)  
        do j = 1, natomi  
           do k = 1, 3  
              fsym(k, j, i) = dzero  
           end do
        end do
     end do
     do i = 1, ntrans  
        do j = 1, ntype  
           natomj = natom(j)  
           do k = 1, natomj  
              !
              !                           -1
              !                  find mtrx    * (rat - tnp)
              !
              do l = 1, 3  
                 cro(l) = dzero  
                 do m = 1, 3  
                    cro(l) = cro(l) + real(mtrx(m, l, i), dp) * &
                         (rat(m, k, j) - tnp(m, i))
                 end do
              end do
              do l = 1, natomj  
                 do m = 1, 3  
                    xdif = abs(cro(m) - rat(m, l, j)) / pi2  
                    idif = int(xdif + eps)  
                    if (abs(xdif - real(idif, dp)) > eps) goto 116  
                 end do
                 kk = l  
                 goto 117  
116              continue  
              end do
              write(9, *) 'subroutine forstresssym:'  
              write(9, *) 'unable to find equivalent atom for k=', k  
              stop  
117           continue  
              !
              !                 now we know that mtrx maps atom k to atom kk
              !                 and mtrx maps atom kk to atom k.
              !
              !
              !                 rotate force and add
              !
              !                  write(9,*) 'atom kk=',kk,'mapped to k=',k
              !                  write(9,*) 'force on atom kk:',
              !     $                 fsum(1,kk,j),fsum(2,kk,j),fsum(3,kk,j)
              !
              !                 apply transpose of mtrx to fsum
              !
              do l = 1, 3  
                 rotfor(l) = dzero  
                 do m = 1, 3  
                    fsym(l, k, j) = fsym(l, k, j) + &
                         real(mtrx(l, m, i), dp) * fsum(m, kk, j)
                    rotfor(l) = rotfor(l) + real(mtrx(l, m, i), dp) * &
                         fsum(m, kk, j)
                 end do
              end do
              !                  write(9,*) 'rotated force:',
              !     $                 rotfor(1),rotfor(2),rotfor(3)
              !                  write(9,'(3i4,6f10.5)') i,j,k,
              !     $                 fsum(1,k,j),fsum(2,k,j),fsum(3,k,j),
              !     $                 fsym(1,k,j),fsym(2,k,j),fsym(3,k,j)
           end do
        end do
        !            write(9,*)
     end do
     !
     !        normalize coordinates and write them into the output array
     !
     do i = 1, ntype  
        natomi = natom(i)  
        do j = 1, natomi  
           do k = 1, 3  
              fsum(k, j, i) = fsym(k, j, i) / real(ntrans, dp)  
           end do
        end do
     end do
  end if

  deallocate(fcart_all)

  !
  if (ipr /= 0) write(9, 102) (ssum(i), i = 1, 6)  
  if (ipr /= 0) write(9, 103) (ssum(i) * fac, i = 1, 6)  
  !
  return  

100   FORMAT(//' FORCE/STRESS ANALYSIS',/,1X,21('*'),///, &
           &     ' FORCES :',9X,'COORD',25X, &
           &     'FORCE (Ry/a.u.)',/ &
           &     ' ',11X,'A1',4X,'A2',4X,'A3',11X, &
           &     '-X-',12X,'-Y-',12X,'-Z-')
101   FORMAT(I3,3X,3F8.3,2X,3F14.10)  
102   FORMAT(//,' STRESS :', &
           &     25X,'SIGMA * V (Ry)',/, &
           &     '(cartes)',2x,'-XX-',8X,'-YY-',8X,'-ZZ-',8X, &
           &     '-XY-',8X,'-YZ-',8X,'-ZX-',/, &
           &     5X,6F12.8)
103   FORMAT(/,36X,'SIGMA (GPa)',/, &
           &     9X,'-XX-',8X,'-YY-',8X,'-ZZ-',8X, &
           &     '-XY-',8X,'-YZ-',8X,'-ZX-',/, &
           &     5X,6F12.5)


end subroutine forstresssym
