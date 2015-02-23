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

!-*-Fortran-*- 
!
!     Angular momentum projection of wavefunctions onto states
!     confined in a spherical box
!     -----------------------------------------------------------------------
!
subroutine project(myproc, emax, nbands, crys, wavefn, nkpts, k_gspace, &
     filename)

  include 'use.h'
  implicit none             ! never uncomment this line. 
  !
  !     INPUT:
  !     ------
  !
  integer, intent(in) :: myproc, nbands, nkpts
  real(dp), intent(in) :: emax
  type(crystal), intent(in) :: crys
  type(complex_gspace_array), intent(in) :: &
       wavefn               ! wave functions for all kpoints, bands
  type(parallel_gspace), intent(in) :: k_gspace(nkpts)
  character(len=*), intent(in) :: filename
  !     
  !     -------------------------- local variables --------------------
  !     
  ! variables for angular momentum projection of wavefunctions
  !
  real(dp), allocatable :: bess(:,:)                  ! more parameters
  real(dp) :: alpha(4,7)
  real(dp) :: gmax, step, rc, tau(3)
  integer, allocatable :: kx(:), ky(:), kz(:)
  integer :: ll, ii, jj, kk, il, im, nrgrid  ! parameters for angular momentum

  !     angular momentum projection of wavefunctions. 

  if (myproc == 0) then
     open(unit=39, file=filename)
  end if
  nrgrid = 3000
  gmax = sqrt(emax)
  allocate(bess(nrgrid, 4))
  rc = dfour
  tau = dzero
  do ii = 1, nbands
     do jj = 1, wavefn%nspin
        do kk = 1, nkpts
           allocate(kx(k_gspace(kk)%length))
           allocate(ky(k_gspace(kk)%length))
           allocate(kz(k_gspace(kk)%length))
           do ll = 1, k_gspace(kk)%length
              kx(ll) = k_gspace(kk)%gvec(1,ll)
              ky(ll) = k_gspace(kk)%gvec(2,ll)
              kz(ll) = k_gspace(kk)%gvec(3,ll)
           end do
           call angular_wavefn(wavefn%data((ii-1)*k_gspace(kk)%length+1, &
                kk, jj), k_gspace(kk)%length, kx, ky, kz, bess, nrgrid, &
                alpha, tau, ii*jj*kk, gmax, rc, step, crys%bdot, crys%bvec)
           deallocate(kx)
           deallocate(ky)
           deallocate(kz)
           if (myproc == 0) then
              write(39,'(3i4)') ii, jj, kk
              do il = 1, 4
                 do im = 1, 2*il-1
                    write(39,'(2i5,2x,f5.3)') il, im, alpha(il, im)
                 end do
              end do
           end if
        end do
     end do
  end do
  deallocate(bess)

  if (myproc == 0) close(39)
  !     end of angular momentum projection of wavefunction


end subroutine project


subroutine angular_wavefn(wavefn, nkpt, kx, ky, kz, bess, ngrid, &
     result, tau, istart, gmax, rc, step, bdot, cart)

  use constants
  implicit none
include 'mpif.h'
  !
  !     INPUT:
  !     ------
  !
  integer, intent(in) :: nkpt, ngrid, istart
  complex(dp), intent(in) :: wavefn(nkpt)
  integer, intent(in) :: kx(nkpt), ky(nkpt), kz(nkpt)
  real(dp), intent(in) :: tau(3), gmax, rc, bdot(3, 3), cart(3, 3)
  !
  !     OUTPUT:
  !     -------
  !
  real(dp), intent(out) :: bess(ngrid, 4), result(4, 7), step
  !     
  !     -------------------------- local variables --------------------
  !     
  integer :: ierr
  integer :: ii, jj, kk, ll, index
  complex(dp) :: alpha(4,7), phase
  real(dp) :: k1, k2, k3, kc(3), norm
  real(dp) :: vec(3), ylm(4, 7), ff(4), gval, angle, sum
  real(dp) :: cellvol, j0, j1, j2, j3, rval, oogval
  real(dp) :: ekin(nkpt)
  !
  ! DESCRIPTION: computes the projected wavefunctions onto
  !              angular momentum eigenstates with a square-well
  !              radial part with 
  !
  !             nkpt = number of g-vectors on this processor
  !             wavefn = complex wavefunction            
  !             rc = cutoff radius of square-well function
  !             alpha = value of scalar product 
  !             tau = relative coordinate of atom in question
  !                     (in lattice coordinates) 
  !             bdot = metric used to compute scalar products of
  !                     vectors in reciprocal space
  !             cart = transformation matrix from lattice coordinates
  !                     to cartesian coordinates
  !             istart = if istart =1 generate bessel function integral
  !
  !
  ! COMPUTE RELEVANT PARAMETERS
  !
  cellvol = abs((pi2**3) / ( &
       cart(1, 1)*cart(2, 2)*cart(3, 3) + &
       cart(1, 2)*cart(2, 3)*cart(3, 1) + &
       cart(1, 3)*cart(2, 1)*cart(3, 2) - &
       cart(3, 1)*cart(2, 2)*cart(1, 3) - &
       cart(3, 2)*cart(2, 3)*cart(1, 1) - &
       cart(3, 3)*cart(2, 1)*cart(1, 2)))
  step = gmax * rc / real(ngrid, dp) 
  if (istart == 1) then
     ! GENERATE BESSEL FUNCTION INTEGRAL
     bess(:,:) = dzero
     do ii = 1, ngrid
        do jj = 1, ii
           rval = jj * step 
           !              
           ! BESSEL FUNCTIONS TIMES rval**2*step
           !
           j0 = sin(rval)*rval*step
           j1 = j0/rval - rval*cos(rval)*step
           j2 = dthree/rval*j1 - j0
           j3 = dfive / rval * j2 - j1
           bess(ii,1) = bess(ii,1) + j0
           bess(ii,2) = bess(ii,2) + j1
           bess(ii,3) = bess(ii,3) + j2
           bess(ii,4) = bess(ii,4) + j3
        end do
     end do
  end if  ! end of startup if then
  !
  ! GENERATE KINETIC ENERGY
  !
  ekin(:) = dzero
  do ii = 1, nkpt
     vec(1) = kx(ii)
     vec(2) = ky(ii)
     vec(3) = kz(ii)
     do jj = 1, 3
        do kk = 1, 3
           ekin(ii) = ekin(ii) + vec(jj) * bdot(jj,kk) * vec(kk)
        end do
     end do
  end do
  !
  ! COMPUTE alpha
  !
  alpha(:,:) = zzero
  do ii = 1, nkpt   ! sum over gvectors
     gval = sqrt(ekin(ii))
     index = floor(gval * rc / step) 
     k1 = kx(ii)
     k2 = ky(ii)
     k3 = kz(ii)
     angle = k1*tau(1) + k2*tau(2) + k3*tau(3)
     phase = exp(pi2*cmplx(dzero, angle))

     if (index > ngrid) index = ngrid

     if (index == 0) then
        alpha(1,1) = alpha(1,1) + phase * wavefn(ii) * rc**3 / dthree
     else
        oogval = done / gval
        kc(1) = cart(1, 1)*k1 + cart(1, 2)*k2 + cart(1, 3)*k3
        kc(2) = cart(2, 1)*k1 + cart(2, 2)*k2 + cart(2, 3)*k3
        kc(3) = cart(3, 1)*k1 + cart(3, 2)*k2 + cart(3, 3)*k3
        k1 = kc(1) * oogval
        k2 = kc(2) * oogval
        k3 = kc(3) * oogval
        do ll = 1, 4
           ff(ll) = bess(index, ll)
        end do
        ylm(1, 1) = done

        ylm(2, 1) = oort3 * k1
        ylm(2, 2) = oort3 * k2
        ylm(2, 3) = oort3 * k3

        ylm(3, 1) = oort3 * k1 * k2
        ylm(3, 2) = oort3 * k2 * k3
        ylm(3, 3) = oort3 * k1 * k3
        ylm(3, 4) = drt3 * 0.5d0 * (k1**2 - k2**2)
        ylm(3, 5) = 1.50d0 * k3**2 - 0.5d0

        ylm(4, 1) = 2.5d0 * k1**3 - 1.5d0 * k1
        ylm(4, 2) = 2.5d0 * k2**3 - 1.5d0 * k2
        ylm(4, 3) = 2.5d0 * k3**3 - 1.5d0 * k3
        ylm(4, 4) = 0.5d0 * drt15 * k1 * (k2**2 - k3**2)
        ylm(4, 5) = 0.5d0 * drt15 * k3 * (k1**2 - k2**2)
        ylm(4, 6) = 0.5d0 * drt15 * k2 * (k3**2 - k1**2)
        ylm(4, 7) = drt15 * k1 * k2 * k3

        do kk = 1, 4
           do ll = 1, 2*kk-1
              alpha(kk, ll) = alpha(kk, ll) + &
                   phase * wavefn(ii) * ff(kk) * ylm(kk, ll) * oogval**3  
           end do
        end do

     end if                      ! end of if  statement for g=0 component
  end do                         ! end of sum over g-vectors
 
  call mpi_allreduce(alpha, alpha, 28, MPI_DOUBLE_COMPLEX, MPI_SUM, &
       MPI_COMM_WORLD, ierr)

  norm = dzero
  do kk = 1, 4
     do ll = 1, 2*kk-1
        norm = norm + alpha(kk, ll)*conjg(alpha(kk, ll))
     end do
  end do

  result(:,:) = alpha(:,:)*conjg(alpha(:,:)) / norm

  return

end subroutine angular_wavefn
