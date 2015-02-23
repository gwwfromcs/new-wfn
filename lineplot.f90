!     @process extchk
!
subroutine lineplot(mode, crys, dscale, data, gs, nlineplot, &
     line_plot, bra, ket, ndim, outfilename, numout)

  use all_to_all_module  
  include 'use.h'  
  implicit none                    ! implicit? Just say no!
  include 'interface.h'  
  include 'all_to_all.h'  
  !
  !     do line plots.
  !
  !     if the number of bins is >0, then a regular line plot is done
  !     if <0, then the average on the plane perpendicular to the line
  !     is computed
  !
  !
  !     1996 Bernd Pfrommer
  !
  !     INPUT
  !     -----
  !
  type(parallel_gspace), intent(in) :: &
       gs                               ! gspace on which the data is defined
  type(crystal), intent(in) :: crys     ! for the metric in realspace
  integer, intent(in) :: ndim, &        ! dimensionality of the tensor data
       nlineplot                        ! number of line plots to do

  integer, intent(in) ::  numout 
  integer, intent(in) ::  mode       ! mode flag 
                                        !  mode = 1 print real part of data
                                        !       = 2 print imaginary part of data
                                        !       = 3 print abs square of data
                                        !       = 0 print absolute of data
                                        ! mode>3 imode=imode-4
                                        !       and project  onto bra
                       !eg ndim=3 mode=1 plot real part of 3 vector components
                       !   ndim=3 mode=5 plot real part of projection of vector
                       !                  in the direction of bra
  character(len=numout), intent(in) :: &   
       outfilename           ! name of output file. Suffix .dx will be appended
  real(dp), intent(in) :: dscale, &     ! scaling factor
       bra(3), &                        ! for tensor data left vector
       ket(3), &                        ! for tensor data right vector
       line_plot(7, nlineplot)          ! line plot ((start,end,bins), lines)
  complex(dp), intent(in) :: &
       data(gs%length, ndim)            ! data to be visualized
  !
  !     ------------- local variables
  !
  integer :: iline, nbins, n, i, j, ig1, ig2, ig3, i3start, i3end, &
       igi, igtot, igmax, iproc ,sproj,imode
  real(dp) :: rp, rip, startp(3), endp(3), d(3), r(3), dr, ovl, nvec(3), nlen
  real(dp), allocatable :: lgvec(:)  
  complex(dp), allocatable :: ldata(:)  
  character(len=3) :: numstring  
  complex(dp) :: ft(ndim), ftv  
  !
  !     ----- variables for the gspace loop ---------
  !
  integer :: igv(4), fftn(4), ffth(4), igv3, irod, iord, igs  
  real(dp) :: gv(3)    ! if required

  fftn(1:3) = gs%fftsize(1:3)  
  fftn(4) = gs%fftsize(3)  
  ffth(:) = fftn(:) / 2  

  if (ndim /= 1 .and. ndim /= 3 .and. ndim /= 9) then  
     write(9, *) '*** NDIM =', ndim, ' NOT IMPLEMENTED YET'  
     call mystop  
  end if

  !
  ! if imode is gt 3 we set a flag (sproj) to project onto the direction
  ! bra. (done this way to avoid changing the calling of lineplot)
  ! Jonathan Yates Paris Dec 2000
 
  sproj=0
  if(mode .gt.  3) then
     imode=mode-4
     sproj=1
  else
     imode=mode
  endif
               
  write(9, *)  
  do iline = 1, nlineplot  
     write(numstring, '(i3.3)') iline  
     if (gs%myproc == 0) then  
        open(18, file = outfilename//'.line.'//numstring, form = &
             'formatted', status = 'unknown')
     end if
     write(9, 10) iline, outfilename//'.line.'//numstring  
     call myflush(9)  
     startp(:) = line_plot(1:3, iline)  
     endp(:) = line_plot(4:6, iline)  
     nbins = int(line_plot(7, iline))  
     if (nbins > 0) then  
        d(:) = (endp(:) - startp(:)) / nbins  
        do n = 0, nbins  
           r(:) = startp(:) + real(n, dp) * d(:)  
           dr = sqrt(dot_product(d, matmul(crys%adot, d)))  
           ft = 0  
           igs = 0  
           do iord = 1, gs%lorder             ! loop through x/y gspace
              irod = gs%order(1, iord)  
              igv(1) = irod / gs%fftsize(2)  
              igv(2) = mod(irod, gs%fftsize(2))  
              igv(3) = gs%order(2, iord)  
              igv(4) = gs%order(3, iord)  
              igv(:) = mod(igv(:) + ffth(:), fftn(:)) - ffth(:)  
              gv(1:2) = real(igv(1:2), dp) + gs%rk(1:2)  
              do igv3 = igv(3), igv(4)                ! loop over z axis
                 gv(3) = real(igv3, dp) + gs%rk(3)  
                 igs = igs + 1  
                 do i = 1, ndim  
                    ft(i) = ft(i) + data(igs, i) * exp(cmplx(dzero, &
                         pi2 * (r(1) * gv(1) + r(2) * gv(2) + r(3) * gv(3)), &
                         dp))
                 end do
              end do
           end do
           call all_sum_all(ft, ndim)  
           if (gs%myproc == 0) then  
1768          format            (8g20.10)  
              ft = ft * dscale  
              if (ndim == 3) then  
                 if (sproj ==1) then  ! project
                    ftv = ft(1)*bra(1)+ft(2)*bra(2)+ft(3)*bra(3)
                    if (imode == 1) write(18, 1768) real(n, dp) * dr, &
                         real(ftv, dp)
                    if (imode == 2) write(18, 1768) real(n, dp) * dr, &
                         aimag(ftv)
                    if (imode == 3) write(18, 1768) real(n, dp) * dr, &
                         abs(ftv) * abs(ftv)
                    if (imode == 0) write(18, 1768) real(n, dp) * dr, &
                         abs(ftv)
                 else 
                    if (imode == 1) write(18, 1768) real(n, dp) * dr, &
                         (real(ft(i), dp), i = 1, ndim)
                    if (imode == 2) write(18, 1768) real(n, dp) * dr, &
                         (aimag(ft(i)), i = 1, ndim)
                    if (imode == 3) write(18, 1768) real(n, dp) * dr, &
                         (abs(ft(i)) * abs(ft(i)), i = 1, ndim)
                    if (imode == 0) write(18, 1768) real(n, dp) * dr, &
                         (abs(ft(i)), i = 1, ndim)
                 end if
              else if (ndim == 9) then  

                 if (sproj ==1) then   ! project
                    ftv = zzero
                    do j = 0, 2
                       ftv = ftv + (bra(1) * ft(j * 3 + 1) + &
                            bra(2) * ft(j * 3 + 2) + bra(3) * ft(j * 3 + 3)) * &
                            ket(j + 1)
                    end do
                    if (imode == 1) write(18, 1768) real(n, dp) * dr, &
                         real(ftv, dp)
                    if (imode == 2) write(18, 1768) real(n, dp) * dr, &
                         aimag(ftv)
                    if (imode == 3) write(18, 1768) real(n, dp) * dr, &
                         abs(ftv) * abs(ftv)
                    if (imode == 0) write(18, 1768) real(n, dp) * dr, abs(ftv)
                 else             
                    if (imode == 1) write(18, 1768) real(n, dp) * dr, &
                         (real(ft(i), dp), i = 1, ndim)
                    if (imode == 2) write(18, 1768) real(n, dp) * dr, &
                         (aimag(ft(i)), i = 1, ndim)
                    if (imode == 3) write(18, 1768) real(n, dp) * dr, &
                         (abs(ft(i)) * abs(ft(i)), i = 1, ndim)
                    if (imode == 0) write(18, 1768) real(n, dp) * dr, &
                         (abs(ft(i)), i = 1, ndim)            
                 end if
              else  
                 ftv = ft(1)  
                 if (imode == 1) write(18, 1768) real(n, dp) * dr, real(ftv, dp)
                 if (imode == 2) write(18, 1768) real(n, dp) * dr, aimag(ftv)
                 if (imode == 3) write(18, 1768) real(n, dp) * dr, &
                      abs(ftv) * abs(ftv)
                 if (imode == 0) write(18, 1768) real(n, dp) * dr, abs(ftv)  
              end if
           end if
        end do
        !
        !     ------------------- compute plane average ----------------
        !
     else if (nbins < 0) then  
        igmax = max(gs%fftsize(1), gs%fftsize(2), gs%fftsize(3))  
        allocate(ldata(igmax))  
        allocate(lgvec(igmax))  
        nvec = endp - startp  
        nlen = dot_product(nvec, matmul(crys%adot, nvec))  
        rip = dot_product(startp, matmul(crys%adot, nvec)) / sqrt(nlen)
        igs = 0  
        igi = 0  
        do iord = 1, gs%lorder                   ! loop through x/y gspace
           irod = gs%order(1, iord)  
           igv(1) = irod / gs%fftsize(2)  
           igv(2) = mod(irod, gs%fftsize(2))  
           igv(3) = gs%order(2, iord)  
           igv(4) = gs%order(3, iord)  
           igv(:) = mod(igv(:) + ffth(:), fftn(:)) - ffth(:)  
           gv(1:2) = real(igv(1:2), dp) + gs%rk(1:2)  
           do igv3 = igv(3), igv(4)             ! loop over z axis
              gv(3) = real(igv3, dp) + gs%rk(3)  
              igs = igs + 1  
              ovl = pi2 * dot_product(nvec, gv)  
              !                  write(9,*)  igs, igv(1), igv(2), igv3
              if (abs(ovl * ovl - nlen * gs%ekin(igs)) < 1.0d-8) then  
                 igi = igi + 1  
                 if (igi > igmax) then  
                    write(9, *) '*** lineplot: fatal error'  
                    call mystop  
                 end if
                 lgvec(igi) = ovl / sqrt(nlen)  
                 if (ndim == 3) then  
                    ftv = bra(1) * data(igs, 1) + bra(2) * data(igs, 2) + &
                         bra(3) * data(igs, 3)
                 else if (ndim == 9) then  
                    ftv = zzero  
                    do j = 0, 2  
                       ftv = ftv + (bra(1) * data(igs, j * 3 + 1) + &
                            bra(2) * data(igs, j * 3 + 2) + &
                            bra(3) * data(igs, j * 3 + 3)) * ket(j + 1)
                    end do
                 else  
                    ftv = data(igs, 1)  
                 end if
                 ldata(igi) = ftv  
                 !                     write(9,*) igi,  lgvec(igi), ldata(igi)
                 !                     write(9,*) data(igs,1:ndim)
              end if
           end do
        end do
        !
        !           now write the data to the file
        !
        close(18)  
        igtot = igi  
        call all_sum_all(igtot)  
        do iproc = 0, gs%nproc - 1  
           if (gs%myproc == iproc) then  
              if (gs%myproc == 0) then  
                 open(18, file = outfilename//'.line.'//numstring, form = &
                      'formatted', status = 'unknown')
              else  
                 open(18, file = outfilename//'.line.'//numstring, form = &
                      'formatted', position = 'append', status = 'unknown')
              end if
              if (gs%myproc == 0) then  
                 write(18, '(a,g12.6)') 'r initial parallel:', rip  
                 write(18, '(a,g12.6)') '|rf-ri| initial parallel:', &
                      sqrt(nlen)
                 write(18, '(a,3g12.6)') 'direction (cartesian):', &
                      matmul(crys%adot, nvec) / sqrt(nlen)
                 write(18, *) igtot  
              end if
              do i = 1, igi  
                 write(18, *) i, lgvec(i), dscale * ldata(i)  
              end do
              call myflush(18)  
              close(18)  
           end if
           call parallel_barrier()  
        end do
        deallocate(lgvec)  
        deallocate(ldata)  
        if (gs%myproc == 0) then  
           allocate(lgvec(igtot))  
           allocate(ldata(igtot))  
           open(18, file = outfilename//'.line.'//numstring, form = &
                'formatted', status = 'unknown')
           read(18, *) ; read(18, *) ; read(18, *)  
           read(18, *) igtot  
           do i = 1, igtot  
              read(18, *) j, lgvec(i), ldata(i)  
           end do
           close(18)  
           open(18, file = outfilename//'.fft.'//numstring, form = &
                'formatted', status = 'unknown')
           do i = 1, abs(nbins)  
              ftv = zzero
              rp = rip + real(i - 1, dp) * sqrt(nlen) / real(abs(nbins) - 1, dp)
              do j = 1, igtot  
                 ftv = ftv + ldata(j) * exp(cmplx(dzero, lgvec(j) * rp, dp))  
              end do
              write(18, *) rp, real(ftv, dp)  
           end do
           deallocate(lgvec)  
           deallocate(ldata)  
        end if
     end if
     if (gs%myproc == 0) close(18)  
     write(9, 20) bra, ket, line_plot(1:7, iline)  
     call myflush(9)  
  end do

10 format(' WRITING LINEPLOT #',i4,' TO FILE: ',a)  

20 format(' BRA:',3f6.2,' KET:',3f6.2,/' PATH: ', 7f8.3)  

end subroutine lineplot
