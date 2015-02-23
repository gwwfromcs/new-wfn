! m4undef.m4
!
! resets various m4 commands to have a prefix m4_
! this can be achieved using option -P but only on some architectures
! any file to be preprocessed with m4 should use this set of macro definitions
!
! David Prendergast, June 6, 2006



! fft_macros.m4
!
!     fft_aux_space(aux,naux,aux2,naux2,dec1,dec2)
!     fft_local_init(N,A,LDA,LOT,DIREC)
!     fft_perm_init(N,AUX)
!
!     fft_multiple_backward(N,A,LOTSTRIDE,LOT,LOOPDUMMY,VECOUT,OFFSET)
!     fft_multiple_forward (N,A,LOTSTRIDE,LOT,LOOPDUMMY,VECOUT,OFFSET)
!
!     fft_backward(N,A,LOTSTRIDE)
!     fft_forward (N,A,LOTSTRIDE)
!
!     fft_convol(N,A,B,POS)
!
!     fft_multiple_scale(N,A,LOTSTRIDE,LOT,LOOPDUMMY,SCALE)
!
!     timeget (t0)
!     timediff(t1,t0)
!
!
!     fft_fcblock     defines blocking factor for fast convolute
!     fft_local_free  release memory
!

!
subroutine start_vectors(ipr, pw_params, ham, densham, eval, neig, xvec, &
     pot_gspace, sub_gspace, vion, crys)
  !
  include 'use.h'
  use flibcalls_module
  use all_to_all_module
  implicit none             ! implicit? Just say no!
  include 'interface.h'
  include 'flibcalls.ph' 
  include 'all_to_all.h'
  !
  !     --------------------  arguments ------------------------------------
  !
  !     INPUT:
  !     -----
  !
  integer, intent(in) :: &
       ipr, &               ! print flag
       neig                 ! number of eigenvectors to search for.
                            ! this is the number of approximate eigenvecs
                            ! returned in the array xvecs
  type(pw_parameter), intent(in) :: &
       pw_params            ! several plane wave parameters 
  type(parallel_gspace), intent(in) :: &
       pot_gspace, &        ! the gspace structure for the potential
       sub_gspace           ! the gspace structure for the submatrix
  complex(dp), intent(in) :: &
       vion(*)              ! screened ionic potential in fourier space
  type(hamiltonian), intent(in) :: &
       ham                  ! all information necessary for the hamiltonian
  type(crystal), intent(in) :: &
       crys                 ! crystal structure info,e.g. inv symmetry
  !
  !     OUTPUT:
  !     ------
  !
  type(dense_hamiltonian), intent(inout) :: &
       densham              ! a structure containing the spectrum etc.
  real(dp), intent(out) :: &
       eval(neig)           ! the approximate eigenvalues on return
  complex(dp), intent(out) :: &
       xvec(*)              ! approx. eigenvectors, on the full gspace 
                            ! xvec(neig * ham%gspace%length)
  !
  !     DESCRIPTION:
  !     -----------
  !
  !     Generates approximate eigenvectors in the full G-space by
  !     diagonalizing a Hamiltonian set up only on a smaller gspace.
  !     The approximate eigenvectors are then expanded onto the full G
  !     -space by either filling in zeroes, or throwing in random
  !     numbers. 
  !
  !     The full spectrum of the approximate eigenvectors is returned 
  !     in densham%matrix. If scalapack is used, only those components
  !     and eigenvectors that belong to this processor are stored in
  !     densham%matrix. The component densham%energy will have the full
  !     spectrum in the subspace.
  !      
  !     1997 Bernd Pfrommer, Andrew Canning 
  !
  !
  !     --------------------- local variables ------------------------------
  !
  integer :: &
       gr(3), gc(3), &
       ngc, ngr, &          ! number of gvecs in cols/rows on this proc
       b, na, ngv, &
       n, &                 ! submatrix size
       ibc, ibr, &          ! block counters
       ic, ir, &            ! col/row counters
       igc, igr, &          ! col/row gvec counters
       nbc, &               ! number of blocks in one column
       nbr, &               ! number of blocks in one row
       nb, nprow, npcol
  integer :: &
       nprocblacs, &        ! number of procs as blacs sees it
       myprow, mypcol, &    ! myprocs place in blacs proc grid
       myproc, nproc, &
       idiff, &             ! index into G-G' array
       ngdiff, &            ! number of G-G' vectors needed on this proc
       ngdiffmax, &         ! max num of G-G' vectors needed on this proc
       ndsqh, ndsq, worksubtot, &
       idum, i1, i2, i, j, k, l, subcomm, info, &
       dg(3), nfft(3), kk, ngmax, &
       ioutput              ! output flags
  real(dp) :: &
       t0, tdiag, x, y, alf, bet, zbuff_real, zbuff_imag
  complex(dp) :: &
       hsum, vnonloc, dcdum, vloc, phfac
  real(dp), allocatable :: &
       rworksub(:)          ! size 3*mxdsdim, for matrix diagonalization
  complex(dp), allocatable :: &
       vnlsubr(:,:), &      ! the nonlocal projectors for the subspace
       vnlsubc(:,:), &      ! the nonlocal projectors for the subspace
       worksub(:)
  integer, allocatable :: &
       wfnlist(:), &        ! which wfns are on this proc
       glist(:), &          ! ordering list for G-G' vectors in subspace
       glistr(:), &         ! ordering list for G vectors in subspace
       glistc(:)            ! ordering list for G' vectors in subspace
  real(dp), external :: gimmetime
  !
  !     ---------  sanity checks etc  --------------------------------
  !
  ioutput = pw_params%output(1)
  myproc = pot_gspace%myproc
  nproc = pot_gspace%nproc
  n = sub_gspace%length

  t0 = gimmetime()

  if (neig > sub_gspace%length) then
     write(9, 920) sub_gspace%length, neig
     call mystop
  end if
  if (sub_gspace%length > densham%gspace%length) then
     write(9, *)  'dense hamiltonian is too small for start vectors!'
  end if

  !     ---------  do the parallel data layout  -----------------------
  !
  !
  !     generate list of G-G' vectors required to set up those
  !     components of H_sub(G,G') residing on this processor
  !

  call layout_scalapack(n, nb, nproc, nprow, npcol)
  ! figure where myproc is in proc grid
  call blacs_gridinfo(blacs_context, nprow, npcol, myprow, mypcol) 
  write(9, 400) n, nprow, npcol, nb

  if (myprow == -1) then
     write(9, *) '*** start_vectors: BLACS ERROR on proc', myproc
     call mystop
  end if

  !
  !
  !     figure out number of blocks per processor in the column/row
  !
  nbc = n / (nb *npcol)
  if (mod(n, (nb * npcol)) > mypcol * nb) nbc = nbc + 1
  nbr = n / (nb *nprow)
  if (mod(n, (nb * nprow)) > myprow * nb) nbr = nbr + 1
  !
  !     ------- set up the local potential from large to small gspace --------
  !
  !
  ngmax = nb * max(nbc, nbr) ! upper bound for g list
  nfft = pot_gspace%fftsize(1:3)
  !
  !     set up the list of column gvectors
  !
  allocate(glistc(ngmax))
  ngc = 0
  do ibc = 0, nbc - 1           ! loop over column blocks
     do ic = (ibc * npcol + mypcol) * nb + 1, &
          min((ibc * npcol + mypcol) * nb + nb, n) 
        ngc = ngc + 1
        gc(1:3) = sub_gspace%gvec(1:3, ic)
        glistc(ngc) = mod(gc(3) + nfft(3), nfft(3)) + &
             mod(gc(2) + nfft(2), nfft(2)) * nfft(3) + &
             mod(gc(1) + nfft(1), nfft(1)) * nfft(3) * nfft(2)
     end do
  end do
  !
  !     set up the list of row gvectors as well
  !
  allocate(glistr(ngmax))
  ngr = 0
  do ibr = 0, nbr - 1          ! loop over row blocks
     do ir = (ibr * nprow + myprow) * nb + 1, &
          min((ibr * nprow + myprow) * nb + nb, n) ! loop over row
        ngr = ngr + 1
        gr(1:3) = sub_gspace%gvec(1:3, ir)
        glistr(ngr) = mod(gr(3) + nfft(3), nfft(3)) + &
             mod(gr(2) + nfft(2), nfft(2)) * nfft(3) + &
             mod(gr(1) + nfft(1), nfft(1)) * nfft(3) * nfft(2)
     end do
  end do
  !
  !     now set up list of g-g' vectors
  !
  ngdiffmax = ngc * ngr         ! should be exact
  allocate(glist(ngdiffmax))

  ngdiff = 0
  do ibc = 0, nbc - 1           ! loop over column blocks
     do ic = (ibc * npcol + mypcol) * nb + 1, &
          min((ibc * npcol + mypcol) * nb + nb, n) ! loop over column
        do ibr = 0, nbr - 1    ! loop over row blocks
           do ir = (ibr * nprow + myprow) * nb + 1, &
                min((ibr * nprow + myprow) * nb + nb, n) ! loop over row
              dg = sub_gspace%gvec(:, ir) - sub_gspace%gvec(:, ic) 
              ngdiff = ngdiff + 1
              if (ngdiff > ngdiffmax) then
                 write(9, *) 'start_vectors: ngdiff too large:', ngdiff
                 call mystop
              end if
              glist(ngdiff) = mod(dg(3) + nfft(3), nfft(3)) + &
                   mod(dg(2) + nfft(2), nfft(2)) * nfft(3) + &
                   mod(dg(1) + nfft(1), nfft(1)) * nfft(3) * nfft(2)
           end do
        end do
     end do
  end do

  if(associated(densham%matrix)) &
       deallocate(densham%matrix) ! make sure we have no leak...
  allocate(densham%matrix(ngr * ngc))
  !
  !     reset local potential from pot_gspace to 2*sub_gspace
  !
  call get_gset(ngdiff, glist(1), 1, densham%matrix(1), &
       2 * sub_gspace%gmax, pot_gspace, vion(1)) 
  !
  !     reset the nonlocal potential projectors from the normal gspace
  !     onto the sub gspace
  !
  if (ham%vnloc%nvecs > 0) then
     allocate(vnlsubr(ngr, ham%vnloc%nvecs))
     call get_gset(ngr, glistr(1), ham%vnloc%nvecs, vnlsubr(1, 1), &
          sub_gspace%gmax, ham%gspace, ham%vnloc%data(1, 1, 1)) 
     allocate(vnlsubc(ngc, ham%vnloc%nvecs))
     call get_gset(ngc, glistc(1), ham%vnloc%nvecs, vnlsubc(1, 1), &
          sub_gspace%gmax, ham%gspace, ham%vnloc%data(1, 1, 1)) 
  end if
  !
  !     --------- construct Hamiltonian matrix in subspace -----------------
  !
  allocate(wfnlist(neig))

  igc = 0 ; idiff = 0
  do ibc = 0, nbc - 1           ! loop over column blocks
     do ic = (ibc * npcol + mypcol) * nb + 1, &
          min((ibc * npcol + mypcol) * nb + nb, n) 
        igc = igc + 1
        if (igc <= neig) wfnlist(igc) = ic
        igr = 0
        do ibr = 0, nbr - 1    ! loop over row blocks
           do ir = (ibr * nprow + myprow) * nb + 1, &
                min((ibr * nprow + myprow) * nb + nb, n) ! loop over row
              igr = igr + 1
              idiff = idiff +1
              vloc = densham%matrix(idiff) ! local part: V(G_col-G_row)
              vnonloc = zzero
              do kk = 1, ham%vnloc%nvecs
                 vnonloc = vnonloc + &
                      vnlsubr(igr, kk) * conjg(vnlsubc(igc, kk)) * &
                      ham%xnorm(kk)
              end do
              if (ir == ic) then ! add kinetic energy for the diagonal
                 hsum = sub_gspace%ekin(ir) + vnonloc
              else
                 hsum = vloc + vnonloc ! loc pot if off-diagonal
              end if
              densham%matrix(igr + (igc - 1) * ngr) = hsum
              !                  write(9,410) ir,ic, hsum, vloc, vnonloc
           end do
        end do
     end do
  end do
  call myflush(9)

  !     ------------ now diagonalize the submatrix -------------------------

  tdiag = gimmetime() 

  if (.not. crys%icomplex) then ! if system has inversion symm, pack hsub -
     ndsq = ngc * ngr
     ndsqh = ndsq / 2

     if (ipr >= 2) then
        write(9, *)  'hsub:'
        do i = 1, ngr
           write(9, '(2000f20.14)') &
                (real(densham%matrix(i + (j - 1) * ngr), dp), j = 1, ngc)
        end do
     end if
     !
     !        pack hsub: Remove all the imaginary parts to prepare call to
     !                   real symmetric eigensolver
     !
     do i = 0, ndsqh - 1
        densham%matrix(i + 1) = cmplx(real(densham%matrix(2 * i + 1), dp), &
             real(densham%matrix(2 * i + 2), dp), dp)
     end do

     if (mod(ndsq, 2) == 1) then
        densham%matrix(ndsqh + 1) = &
             cmplx(real(densham%matrix(ndsq), dp), dzero, dp)
     end if

     !        diagonalize it

     call rdiag_scalapack(blacs_context, neig, n, nb, nprow, npcol, &
          densham%matrix(1), ngc, ngr, densham%energy(1))

     !        now we have to unpack it from real to complex again

     !        ! the if statement must go first!

     if (mod(ndsq, 2) == 1) then
        densham%matrix(ndsq) = &
             cmplx(real(densham%matrix(ndsqh + 1), dp), dzero, dp)
     end if

     do i = ndsqh - 1, 0, -1
        zbuff_real = real(densham%matrix(i + 1), dp)
        zbuff_imag = aimag(densham%matrix(i + 1))

        densham%matrix(2 * i + 1) = zbuff_real
        densham%matrix(2 * i + 2) = zbuff_imag
     end do

     if (ipr >= 2) then
        write(9, *)  'packed eigenvectors hsub (column vectors):'
        do i = 1, ngr
           write(9, '(3f20.14)') &
                (real(densham%matrix(i + (j - 1) * ngr), dp), j = 1, ngc)
        end do
     end if

  else                      ! the complex case is easier. 

     call zdiag_scalapack(blacs_context, neig, n, nb, nprow, npcol, &
          densham%matrix(1), ngc, ngr, densham%energy(1))
     !        
     !        enforce the phase
     !
     !         call all_split_all(mypcol,iam,subcomm)
     !
     !         do igc =1,ngc
     !            if(myprow.eq.0) then
     !               x = dble(densham%matrix((igc-1)*ngr+1))
     !               y = aimag(densham%matrix((igc-1)*ngr+1))
     !               if(abs(y).gt.1d-10) then
     !                  bet = 1/sqrt(1+(x/y)**2)
     !                  alf = -x*bet/y
     !                  phfac = cmplx(alf,bet,kind=8)
     !               else
     !                  phfac = cmplx(1.d0,0.d0,kind=8)
     !               endif
     !            else
     !               phfac =0
     !            endif
     !            call all_sum_all_sub(subcomm,phfac)
     !            densham%matrix((igc-1)*ngr+1:igc*ngr)= &
     !                densham%matrix((igc-1)*ngr+1:igc*ngr)*phfac
     !
     !         end do
     !         call mpi_comm_free(subcomm, info)
     !


  end if

  tdiag = gimmetime() - tdiag

  if (ipr >= 1) then
     write(9, 100)
     write(9, 110) (densham%energy(l), l = 1, neig)
     call myflush(9)
  end if
  densham%neig = sub_gspace%length ! have all of them

  !     ---------- expand the packed eigenvectors into arrays ------------
  !
  !      do j=0, neig-1 ! check norm
  !         hsum =0
  !         do i=1,sub_gspace%length         
  !            hsum=hsum+densham%matrix(j*sub_gspace%length+i)*&
  !                conjg(densham%matrix(j*sub_gspace%length+i))
  !         end do
  !         call all_sum_all(hsum)
  !         write(9,*) 'unexpanded norm',j+1,hsum
  !      end do

  xvec(1:neig * ham%gspace%length) = zzero

  call put_gset(ngr, glistr, min(neig, ngc), wfnlist, &
       densham%matrix(1), sub_gspace%gmax, ham%gspace, xvec(1), neig)

  !      do j=0, neig-1 ! check norm
  !         hsum =0
  !         do i=1,ham%gspace%length         
  !            hsum=hsum+xvec(j*ham%gspace%length+i)*&
  !                conjg(xvec(j*ham%gspace%length+i))
  !         end do
  !         call all_sum_all(hsum)
  !         write(9,*) 'expanded norm',j+1,hsum
  !      end do

  deallocate(glistr)
  deallocate(glistc)
  deallocate(glist)
  deallocate(wfnlist)

  eval(1:neig) = densham%energy(1:neig)

  if (pw_params%random_startvec > dzero) &
       call randomize_startguess(crys, ham%gspace, &
       sub_gspace%gmax, pw_params%random_startvec, neig, xvec(1))

  if (ipr >= 3) then
     write(9, 300)
     do j = 1, ham%gspace%length
        if (.not. crys%icomplex) then
           write(9, 310) (real(xvec((l - 1) * ham%gspace%length + j), dp), &
                l = 1, neig)
        else
           write(9, 310) &
                (xvec((l - 1) * ham%gspace%length + j), l = 1, neig)
        end if
     end do
  end if



  if (ham%vnloc%nvecs > 0) then
     deallocate(vnlsubr)
     deallocate(vnlsubc)
  end if

  if (iand(ioutput, 65536) == 65536) then
     write(9, 930) 2 * (ngmax + ngdiffmax) * 8 + 16 * ngdiff + &
          (ngc + ngr) * ham%vnloc%nvecs * 16
  end if
  if (iand(ioutput, 8) == 8) then
     write(9, 910) sub_gspace%length, tdiag, gimmetime() - t0
     call myflush(9)
  end if

  return

100 format (/' Eigenvalues from submatrix diagonalization:',&
       /' -------------------------------------------')
110 format (4f19.14)
200 format ('Eigenvectors from submatrix diagonalization:',&
       /'-------------------------------------------')
210 format (1008f9.3)

300 format ('Startvectors after unpacking:',&
       /'-------------------------------------------')
310 format (400f20.14)

400 format (' START GUESS SIZE ',i5,&
       /' PROCESSOR GRID ',i3,' x ',i3,' BLOCKSIZE ',i3)
410 format (2i5,8f20.12)

910 format(' START GUESS SIZE ',i4,' TOOK ',&
       f8.2,' SECONDS FOR DIAG, AND ',f8.2,' TOTAL')

920 format(' *** ERROR: the subspace for the starting guess',&
       /' ***        has size ',i5,', but you have asked for',i5,&
       /' ***        bands to be computed. Increase the energy',&
       /' ***        cutoff of the subspace.')

930 format(' START GUESS MEMORY USAGE IS ',i10,' BYTES')

end subroutine start_vectors


 subroutine mdsyev_wrapper(v, l, m, a1, n, a2, a3, k, i)
  include 'flibcalls.ph' 
    character :: v, l
    real(kind(1.0d0)) :: a1(n*m)
    real(kind(1.0d0)):: a2(n), a3(3*n)
    integer :: m, n, k, i


   call mdsyev(v, l, m, a1(1), n, a2(1), a3(1), k, i)

  return
  end subroutine 

