!-*-Fortran-*-
!
!     @process extchk
!
subroutine crstl(ioutput, crys, syms, miscflag,iter_relax)  
  !
  use flibcalls_module  
  include 'use.h'  
  implicit none  
  include 'interface.h'  
  !
  !     INPUT:
  !     -----
  !
  integer iter_relax
  integer, intent(in) :: miscflag, &                       ! for inv symm rest
       ioutput                                           ! for timing printout
  !
  !     INPUT/OUTPUT:
  !     ------------
  !
  type(crystal), intent(inout) :: crys  
  type(symmetry), intent(inout) :: syms  
  !
  !     DESCRIPTION:
  !     -----------
  !
  !     Subroutine computes various crystal parameters and the symmetry
  !     Based on Alberto's crystal subroutine
  !     1996 Bernd Pfrommer
  !
  !     input:
  !     -----
  !
  !     crys%natom,     crys%ntype,     crys%mxdatm,
  !     crys%rat,       crys%avec,      crys%szet,
  !     crys%nameat
  !
  !     syms%ntrans   (if <0, compute symmetry operations)
  !
  !     action:
  !     ------
  !
  !     computes
  !
  !     crys%vcell,    crys%adot,     crys%bdot,
  !     crys%bvec
  !
  !
  !     syms%ntrans,   syms%mtrx,     syms%rmtrx,
  !     syms%rsymmat,  syms%tnp
  !
  !
  !     --------------- local variables ----------------------------------
  !
  real(dp) :: t0, matinvert  
  integer :: i, j, k, nops,icar, itest, lmtrx(48, 3, 3), &
       lrmtrx(48, 3, 3), invers_no, ipr
  real(dp) :: a1(3), a2(3), a3(3), am(3), ba(3), b1(3), b2(3), b3(3), &
       ltnp(48, 3)
  character(len=10) :: group
  real(dp), allocatable :: lcoorat(:)  
  logical :: rest_inv_flag
  !     ..
  !     .. External Functions ..
  logical, external :: rest_inv  
  real(dp), external :: gimmetime
  integer, SAVE ::nops_orig
  !     ..
  !
  !     ----------------------------------------------------------------
  group = 'unknown'  

  ipr = 1  
  lmtrx = 0

  ltnp = dzero  
  t0 = gimmetime()  

  allocate(lcoorat(3 * crys%mxdatm * crys%ntype))  

  if (ipr == 1) then  
     write(9, 9000)
     write(9, 9020)
     do i = 1, crys%ntype  
        do j = 1, crys%natom(i)
           write(9, 9030) crys%nameat(i), j, crys%rat(1:3, j, i), &
                crys%szet(i)
        end do
     end do
     write(9, 9021)  
     do i = 1, crys%ntype  
        do j = 1, crys%natom(i)
           write(9, 9030) crys%nameat(i), j, &
                matmul(crys%avec, crys%rat(1:3, j, i)), crys%szet(i)
        end do
     end do
  end if

9529 format('newtype',1x,a2,1x,f9.5)
9530 format('coord',5x,3f18.14)

  !
  !      Compute reciprocal lattice and reset scales.
  !
  !               T                       T        -1          T -1
  !              B A  = 2pi (1), so b =  B  = 2pi(A  )  = 2pi(a )
  !
  !     and we can compute the determinant (celvol) and the inverse.
  !     But note that the 2pi factor is NOT needed here yet. (It
  !     will be put in below)
  !
  crys%bvec = transpose(crys%avec)  
    crys%vcell=matinvert(crys%bvec)
  if (crys%vcell < dzero) then  
     if (ipr == 1) then  
        write(9, '(/,x,a,/)') 'Note: The lattice vectors as given'// &
             ' do not form a direct triad'
     end if
     crys%vcell = -crys%vcell  
  end if
  !
  !
  if (ipr == 1) then  
     write(9, 9050) crys%vcell  
  end if
  !
  !     Compute the inner products ba, put in the 2pi factor and clean up.
  !
  a1(:) = crys%avec(:, 1)  
  a2(:) = crys%avec(:, 2)  
  a3(:) = crys%avec(:, 3)  
  crys%bvec = pi2 * crys%bvec  
  !
  am(1) = sqrt(a1(1) * a1(1) + a1(2) * a1(2) + a1(3) * a1(3))
  am(2) = sqrt(a2(1) * a2(1) + a2(2) * a2(2) + a2(3) * a2(3))
  am(3) = sqrt(a3(1) * a3(1) + a3(2) * a3(2) + a3(3) * a3(3))
  !
  ba(:) = pi2 / am(:)  
  b1(:) = crys%bvec(:, 1)  
  b2(:) = crys%bvec(:, 2)  
  b3(:) = crys%bvec(:, 3)  
  !
  !      printout
  !
  if (ipr == 1) then  
     write(9, 9130) a1(:), a2(:), a3(:)  
     write(9, 9140) b1(:), b2(:), b3(:)  
  end if
     write(9,*)
     write(9, 9020)
  
  do i = 1, crys%ntype
     write(9, 9529) crys%nameat(i), crys%szet(i)
     do j = 1, crys%natom(i)
        write(9, 9530) (crys%rat(k, j, i), k = 1, 3)
     end do
  end do
  !
  !     Compute metrics crys%bdot(i,j) and crys%adot(i,j)
  !
  crys%bdot = matmul(transpose(crys%bvec), crys%bvec)  
  crys%adot = matmul(transpose(crys%avec), crys%avec)  
  !     Symmetry operations for the symmetry group of the crystal.
  !
  !
  !      If nops is negative,  calculate the symmetry operations.
  !      If icar=1 the matrices and primitive translations should
  !      be given in cartesian coordinates and the translation vectors
  !      should be scaled the same way as the basis vectors (atomic units
  !      divided by scale). Otherwise the matrices must be given in compon
  !      of the lattice wave vectors and the translations in units of the
  !      lattice basis vectors.
  !
  !      if ipr=1 the matrices are printed.
  !
  !      itest is ignored. The test is always performed.
  !
  !      group is some group identifier.
  !
  !
  !     transform the shit to weird coordinates
  !
  do k = 1, crys%ntype  
     do i = 1, crys%natom(k)  
        do j = 1, 3  
           lcoorat((j-1)*crys%mxdatm*crys%ntype + (i-1)*crys%ntype + k) = &
                crys%rat(j, i, k)
        end do
     end do
  end do
  
  if (iter_relax .eq. 0) then
    nops_orig=syms%ntrans
  end if

  nops=nops_orig

  if (nops <= 0) then    
     call symgen(ipr, crys%avec, lcoorat, crys%natom, crys%ntype, &
          invers_no, nops,syms%ntrans, lmtrx, ltnp, crys%mxdatm, lrmtrx, &
          syms%rsymmat,syms%ind_rot)
     !
     !        Restore manifest inversion symmetry if possible and desired
     !
     if (iand(miscflag, 1) == 0) then  
        ! Do not remove this reference to rest_inv_flag: it overcomes a 
        !  bug in the SGI compiler
        rest_inv_flag = rest_inv(ipr, invers_no, crys%ntype, crys%natom, &
             lcoorat(1), crys%nameat, ltnp, crys%mxdatm)
        if (rest_inv_flag) then
           call symgen(ipr, crys%avec, lcoorat, crys%natom, crys%ntype, &
                invers_no, nops,syms%ntrans, lmtrx, ltnp, crys%mxdatm, lrmtrx, &
                syms%rsymmat,syms%ind_rot)
        end if
     end if
     do k = 1, 48
        do i = 1, 3  
           syms%mtrx(i,:, k) = lmtrx(k, i,:)  
           syms%rmtrx(i,:, k) = lrmtrx(k, i,:)  
        end do
     end do
     syms%tnp = pi2 * transpose(ltnp)  
  else  
     call inpmat(crys, syms)  
  end if

  if (syms%ntrans <= 0) then  
     write(9, *) ' *** SYMMETRY  COULD NOT FIND IDENTITY.'  
     write(9, *) ' ***  ERROR    CHECK INPUT FILE COORDINATES'  
     write(9, *) ' ***           ARE TWO ATOMS AT SAME COORD?'  

     call mystop  
  end if
  !     ------------- rescale coordinates and clean up -------------------
  !
  do k = 1, crys%ntype  
     do i = 1, crys%natom(k)  
        do j = 1, 3  
           crys%rat(j, i, k) = pi2 * &
                lcoorat((j-1)*crys%mxdatm*crys%ntype + (i-1)*crys%ntype + k)
        end do
     end do
  end do
  !
  deallocate (lcoorat)  
  !
  !     --------- debugging printout ----------------------
  !
  !      do i=1, syms%ntrans
  !         write(9,9010) i, syms%mtrx(:,:,i), syms%tnp(:,i)
  !         write(9,9010) i, syms%rmtrx(:,:,i)
  !         write(9,9022) i, syms%rsymmat(:,:,i)
  !      end do
  if (iand(ioutput, 8) == 8) then  
     write(9, 910) gimmetime() - t0  
     call myflush(9)
  end if

  return  

 910  format(' TIME FOR CRYSTAL SYMMETRY DETERMINATION: ',f12.3,' SECONDS')
 9000 format(1x,'Lattice information',/1x,19('-'))  
 9010 format(i5,3(3x,3i3),4x,3f8.5)  
 9020 format(/1x,'atomic positions relative to lattice vectors:'/)  
 9021 format(/1x,'atomic positions in cartesian coordinates:'/)  
 9022 format(i5,3(3x,3f10.6))  
 9030 format(1x,a2,' #',i4,':',3f20.14,4x,f5.2)  
 9050 format(/' Cell volume =',f14.4,' au.')  
 9130 format(/' lattice vectors (a.u.)',/3(/1x,3f20.14))
 9140 format(/' reciprocal basis (a.u.)',/3(/1x,3f20.14))

end subroutine crstl
