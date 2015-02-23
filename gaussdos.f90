!
program gaussdos  
  !
  use constants
  implicit none  
  !
  !     computes Gaussian-broadened DOS from the EIG file
  !
  !
  !
  !     1997 Bernd Pfrommer
  !
  !     DESCRIPTION:
  !     -----------
  !
  !     reads the file EIG, and writes out a gaussian-broadened DOS
  !     The DOS is in states/(unit cell and eV), on the axis, it shows
  !     energy in eV with respect to the fermi level.
  !
  !
  !
  !     The weights for different angular momenta, different bands,
  !     kpoints, and spins. A fancy interpolation between the function
  !     values at different kpoints is done inside the subroutine
  !     spec_tet.
  !
  ! CHANGES
  ! 
  !  11/27/00 - changed p(7, to p(lmax, due to new angular.f90
  !
  !     INPUT:
  !     -----
  !     reads from stdin the gaussian smearing and the number of bins
  !
  !
  !     -------------------local variables ------------------------------
  !
  real(dp) :: invsig, fac, &
       eshift(2), &         ! energy shift
       efermi(2), &         ! fermi level
       e, &                 ! energy
       deltae, &            ! energy difference between two mesh points
       d, &                 ! dos
       nos, &               ! number of states
       sigma, &             ! gaussian smearing in eV
       window(2)            ! the energy window for plot (in Rydb)
  real(dp), allocatable :: &
       dos(:,:), &          ! dos
       energy(:,:,:), &     ! energy of the different states
       p(:,:,:,:,:,:), &    ! the weights
       w(:)                 ! symmetry weights of kpoints
  integer :: irk, is, n, it, ia, nrk, i, j, nspin, maxb, minb, &
       ntype, mxdatm, l, &
       nbin                 ! number of bins for the DOS plot
  logical :: ip  
  integer, allocatable :: nband(:,:), natom(:), iord(:)  

  character(len=1024) :: dummystring  
  !
  !     ------------ read in ---------------------------------------
  !
  read(5, *) sigma  
  sigma = sigma / ryd  

  read(5, *) nbin  
  write(0, *) 'gaussian smearing parameter [eV]=', sigma * ryd  
  write(0, *) 'number of mesh points=', nbin  

  invsig = done / sigma  

  open(11, file = 'EIG', status = 'old', form = 'formatted')  
  read(11, '(a32,2g20.12)') dummystring, efermi (1), efermi (2)  
  read(11, '(a47,2g20.12)') dummystring, eshift (1), eshift (2)  
  read(11, *)  
  read(11, *) nrk  

  allocate(w(nrk))  
  read(11, *)  
  do irk = 1, nrk  
     read(11, *) j, w(irk)  
  end do
  read(11, *)  
  read(11, *) nspin  
  read(11, *)  
  read(11, *) maxb  
  read(11, *)  
  allocate(nband(nrk, nspin))  
  read(11, *) nband  
  read(11, *)  

  read(11, *) minb  

  read(11, *)  
  allocate(energy(maxb, nrk, nspin))  
  energy = dzero
  do is = 1, nspin  
     do irk = 1, nrk  
        read(11, *) i, j, energy(1:j, irk, is)  
     end do
  end do

  energy = energy + eshift(1) - efermi(1)  
  read(11, *)  

  read(11, *) ip  
  if (.not. ip) then  
     ntype = 1  
     allocate(natom(ntype))  
     natom(1) = 1  
     mxdatm = 1  
  else  
     read(11, *)  
     read(11, *) ntype  
     read(11, *)  
     read(11, *) mxdatm  
     allocate(natom(ntype))  
     read(11, *)  
     read(11, *) natom  
     allocate(p(llmax, mxdatm, ntype, maxb, nrk, nspin))
     do is = 1, nspin  
        read(11, *)
        read(11, *) j  
        do irk = 1, nrk  
           read(11, *)  
           read(11, *) j  
           do n = 1, nband(irk, is)  
              read(11, *)  
              read(11, *) j  
              do it = 1, ntype  
                 read(11, *)  
                 read(11, *) j  
                 do ia = 1, natom(it)  
                    read(11, *) j, p(1:7, ia, it, n, irk, is)  
                 end do
              end do
           end do
        end do
     end do

  end if
  !
  !     ----------- read in ordering info --------------------
  !
  open(7, file = 'ORD')  

  allocate(iord(natom(1)))  
  write(0, *) 'ordering info:'  
  do i = 1, natom(1)  
     read(7, *) sigma, iord(i)  
     write(0, *) i, iord(i)  
  end do

  close(7)  
  !
  !     ---------- set energy window ---------------------------
  !
  window(1) = minval(energy) - 0.1d0
  window(2) = maxval(energy) + 0.1d0

  deltae = (window(2) - window(1)) / real(nbin, dp)  

  fac = dtwo * invsig / sqrt(pi) / real(nspin, dp)  
  allocate(dos(nbin + 1, mxdatm))  

  open(12, file = 'DOSMATRIX')  
  do is = 1, nspin  
     do it = 1, ntype  
        do ia = 1, natom(it)  
           write(9, '(''# spin:'',i2,''type:'',i3,''atom:'',i3)') is, it, ia
           nos = 0
           do i = 0, nbin  
              e = window(1) + real(i, dp) * deltae  
              d = dzero
              do irk = 1, nrk  
                 do n = 1, nband(nrk, nspin)  
                    d = d + w(irk) * fac * p(1, ia, it, n, irk, is) * &
                         exp(-(invsig * (e - energy(n, irk, is)))**2)
                 end do
              end do
              nos = nos + d * deltae  
              if (d < 1.0d-10) d = dzero
              dos(i + 1, ia) = d / ryd  
              write(9, *) e * ryd, d / ryd, nos  
           end do
           if (ia < natom(it)) write(9, '(''&'')')  
        end do
        do ia = 1, natom(it)  
           l = iord(ia)  
           write(12, *) (dos(i, l), i = 1, nbin + 1)  
        end do
        if (it < ntype) write(9, '(''&'')')  
     end do
     if (is < nspin) write(9, '(''&'')')  
  end do

  close(12)  
  close(11)  

end program gaussdos
