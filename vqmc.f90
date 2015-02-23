!
subroutine vqmc(fname, crys, syms, bands, gs, kgs, wfn, kpoints, &
     efermi, smearing)
  !
  include 'use.h'  
  use all_to_all_module  
  use vqmc_gspace_module  
  use vqmc_isrt_module  
  use vqmc_wavefn_module  
  use vqmc_level_module  
  implicit none                 ! implicit? Just say no!
  include 'interface.h'  
  include 'all_to_all.h'  
  !
  !     writes a data array to disk such that it can be read and used in
  !     later vqmc calculations
  !
  !     INPUT:
  !     -----
  !
  type(crystal), intent(in) :: crys  
  type(symmetry), intent(in) :: syms  
  type(complex_gspace_array), intent(in) :: &
       wfn                      ! the wave functions
  type(band), intent(in) :: bands  
  type(parallel_gspace), intent(in) :: &
       kgs (bands%nrk), &       ! the gspaces at the different kpoints
       gs                       ! the big gspaces
  type(kpoint), intent(in) :: &
       kpoints                  ! all the kpoints
  real(dp), intent(in) :: &
       efermi, &                ! the fermi energy
       smearing                 ! the gaussian smearing parameter [Ryd]
  character(len=*), intent(in) :: &
       fname                    ! file name
  !
  !     ---------------- local variables ---------------------------------
  !
  type(gspace) :: gsstruc       ! the structure containing the whole gspace
  type(isrt) :: isortstruc      ! complete isort array
  type(wavefn) :: wavefnstruc   ! complete wave functions
  type(level) :: levelstruc     ! energy levels and occupation numbers
  integer :: mpdim, &           ! largest number of plane waves at all kpoints
       ngsm, ngtot, mpdim_all, irk, is, n, i, ndim  
  real(dp), parameter :: occup_epsilon = 1.0d-3

  mpdim = wfn%length / bands%max    ! should work.
  !
  !     erase old file
  !
  open(21, file = 'VQMC', status = 'replace')  
  close(21)  
  !
  !     generate k-point structure, and write it to disk
  !
  call write_crystal('VQMC', crys)  
  call write_symmetry('VQMC', syms)  
  call write_kpoint('VQMC', kpoints)  
  !
  !     generate level structure, and write it to disk
  !
  levelstruc%nrk = bands%nrk  
  levelstruc%mnband = bands%max  
  levelstruc%nspin = wfn%nspin  
  allocate(levelstruc%nband(bands%nrk, wfn%nspin))  
  do is = 1, wfn%nspin  
     levelstruc%nband(:, is) = bands%nband(:, is)  
  end do
  allocate(levelstruc%ifmin(bands%nrk, wfn%nspin))  
  allocate(levelstruc%ifmax(bands%nrk, wfn%nspin))  
  allocate(levelstruc%occup(bands%max, bands%nrk, wfn%nspin))  
  allocate(levelstruc%el (bands%max, bands%nrk, wfn%nspin))  
  levelstruc%ifmin = 0
  levelstruc%ifmax = 0
  levelstruc%occup = dzero
  levelstruc%el = dzero

  do is = 1, wfn%nspin  
     do irk = 1, bands%nrk  
        do n = 1, bands%nband(irk, is)  
           levelstruc%occup(n, irk, is) = bands%occup(n, irk, is) / &
                (kpoints%w(irk) * real(wfn%nspin, dp))
           levelstruc%el(n, irk, is) = bands%energy(n, irk, is)  
           if (bands%occup(n, irk, is) > occup_epsilon) then  
              levelstruc%ifmax(irk, is) = n  
           end if
        end do
        levelstruc%ifmin(irk, is) = 1  
     end do
  end do
  levelstruc%smearing = smearing  
  levelstruc%fermilevel = efermi  

  call write_level('VQMC', levelstruc)  
  !
  !     collect gvecs from all the processors into one array
  !
  call gspace_gather(gs%myproc, gs%nproc, gsstruc, gs%fftsize, &
       gs%gvec, gs%length, gs%istar, gs%nstar, gs%mstar(1), gs%phase(1), &
       gs%inds(1))
  !
  !     write it to the file
  !
  call write_gspace('VQMC', gsstruc)  
  !
  !     Allocate space for wave functions and isort array. Initialize.
  !
  mpdim_all = mpdim  

  if (gs%nproc > 1) call all_sum_all(mpdim_all)  
  isortstruc%nrk = bands%nrk  
  allocate(isortstruc%kndim(bands%nrk))  
  isortstruc%lda = mpdim_all  
  allocate(isortstruc%isrt(mpdim_all, bands%nrk))  
  isortstruc%isrt = 0    ! set array to zero
  wavefnstruc%nrk = bands%nrk  
  wavefnstruc%nspin = wfn%nspin  
  wavefnstruc%mnband = maxval(bands%nband)  
  wavefnstruc%lda = mpdim_all  
  allocate(wavefnstruc%kndim(bands%nrk))  
  allocate(wavefnstruc%nband(bands%nrk, wfn%nspin))  
  allocate(wavefnstruc%wfn(mpdim_all, wavefnstruc%mnband, bands%nrk, &
       wfn%nspin))
  do is = 1, wfn%nspin  
     !         wavefnstruc%nband(:,is) = bands%nband(:,is)
     ! to save space
     wavefnstruc%nband(:, is) = bands%min(is)  
  end do

  wavefnstruc%wfn = zzero    ! set array to zero

  do i = 1, bands%nrk  
     wavefnstruc%kndim(i) = kgs(i)%length  
     isortstruc%kndim(i) = kgs(i)%length  
  end do
  !
  !     now fill the data part of the wave function and isort structures
  !
  do irk = 1, bands%nrk  
     ndim = kgs(irk)%length       ! on this processor
     call isort_gather(1, isortstruc%isrt(1, irk), gsstruc, kgs(irk))
     do is = 1, wfn%nspin  
        !            do n=1,bands%nband(irk,is)
        do n = 1, wavefnstruc%nband(irk, is)  
           call data_gather(1, kgs(irk)%myproc, kgs(irk)%nproc, &
                wfn%data(1 + (n - 1) * ndim, irk, is), &
                wavefnstruc%wfn(1, n, irk, is), ndim)
        end do
     end do
  end do
  call write_isrt('VQMC', isortstruc)  
  call write_wavefn('VQMC', wavefnstruc)  
  call iolib_close('VQMC')  
  call destroy_isrt(isortstruc)  
  call destroy_wavefn(wavefnstruc)  
  call destroy_gspace(gsstruc)  
  call destroy_level(levelstruc)  

  return

end subroutine vqmc
