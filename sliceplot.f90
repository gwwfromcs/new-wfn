!     @process extchk
!
subroutine sliceplot(mode, crys, dscale, data, gs, nsliceplot, &
     slice_plot, bra, ket, ndim,outfilename)

  use all_to_all_module  
  include 'use.h'  
  implicit none                    ! implicit? Just say no!
  include 'interface.h'  
  include 'all_to_all.h'  
  !
  !     Will plot 2D slices of data
  !    
  !     Dec2000 Jonathan Yates ,Paris (jry20@phy.cam.ac.uk) 
  !     basically a hacked version of Bernd Prommer's lineplot
  !     
  !
  !
  !     
  !
  !     INPUT
  !     -----
  !
  type(parallel_gspace), intent(in) :: &
       gs                               ! gspace on which the data is defined
  type(crystal), intent(in) :: crys     ! for the metric in realspace
  integer, intent(in) :: ndim, &        ! dimensionality of the tensor data
       nsliceplot                       ! number of slices to do
  integer, intent(in) :: mode           ! mode flag
                                        !  mode = 1 print real part of data
                                        !       = 2 print imaginary part of data
                                        !       = 3 print abs square of data 
                                        !       = 0 print absolute of data
                                        !  mode>3 imode=imode-4      
                                        !       and project  onto bra
                       !eg ndim=3  mode=1 plot real part of 3 vector components
                       !   ndim=3  mode=5 plot real part of projection of vector
                       !                  in the direction of bra 


  character(len=7), intent(in) :: &
       outfilename                      ! name of output file. 
  real(dp), intent(in) :: dscale, &     ! scaling factor
       bra(3), &                        ! for tensor data left vector
       ket(3), &                        ! for tensor data right vector
       slice_plot(11, nsliceplot)       ! slice plot ((start,end,top,
                                        !  lengthbins, upbins), slices)
  complex(dp), intent(in) :: &
       data(gs%length, ndim)            ! data to be visualized
  !
  !     ------------- local variables
  !
  integer :: islice, nbins, n, i, j , upn, upnbins,sproj,imode
  real(dp) ::  startp(3), endp(3), topp(3), d(3),&
     r(3), dr,  nvec(3), upd(3),upr(3), updr
  character(len=3) :: numstring  
  complex(dp) :: ft(ndim), ftv ! ,sum
  !
  !     ----- variables for the gspace loop ---------
  !
  integer :: igv(4), fftn(4), ffth(4), igv3, irod, iord, igs  
  real(dp) :: gv(3)    ! if required

  fftn(1:3) = gs%fftsize(1:3)  
  fftn(4) = gs%fftsize(3)  
  ffth(:) = fftn(:) / 2  

  
  sproj=0
  if(mode .gt. 3) then
     imode=mode-4
     sproj=1
  else 
  imode=mode
  endif

  if (ndim /= 1 .and. ndim /= 3 .and. ndim /= 9) then  
     write(9, *) '*** NDIM =', ndim, ' NOT IMPLEMENTED YET'  
     call mystop  
  end if

  !
  !Normalise bra and ket
  ! hmmmm we can't do that cos bra and ket are intent in
  ! must change this sometime.....
!  bra(:)=bra(:)/sqrt(bra(1)**2+bra(2)**2+bra(3)**2) 
!  ket(:)=ket(:)/sqrt(ket(1)**2+ket(2)**2+ket(3)**2)


  write(9, *)  
  do islice = 1, nsliceplot  
     write(numstring, '(i3.3)') islice  
     if (gs%myproc == 0) then  
        open(18, file = outfilename//'.slice.'//numstring, form = &
             'formatted', status = 'unknown')
     end if
     write(9, 10) islice, outfilename//'.slice.'//numstring  
     call myflush(9)  
     startp(:) = slice_plot(1:3, islice)  
     endp(:) = slice_plot(4:6, islice)  
     topp(:) =   slice_plot(7:9, islice)
       nbins = int(slice_plot(10, islice))
     upnbins = int(slice_plot(11, islice))
       upd(:) = (topp(:) - startp(:)) / upnbins
        updr = sqrt(dot_product(upd, matmul(crys%adot, upd)))  
        do upn = 0,upnbins              !loop over lines
        d(:) = ( endp(:) - startp(:)  ) / nbins  
        do n = 0, nbins                 ! loop along lines
           r(:) = startp(:) + real(n, dp)   * d(:) &
                            + real(upn, dp) * upd(:)  
           dr = sqrt(dot_product(d, matmul(crys%adot, d)))  
           ft = (0.0,0.0)  
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
            call myflush(9)
           if (gs%myproc == 0) then  
1768          format            (8g20.10)  
              ft = ft * dscale  
              if (ndim == 3) then   !for a vector quantity
                 if (sproj ==1 ) then !project onto bra
                    ftv=ft(1)*bra(1)+ ft(2)*bra(2)+ ft(3)*bra(3)               
                    if (imode == 1) write(18, 1768) real(n, dp) * dr, &
                         real(upn, dp) * updr, real(ftv, dp)
                    if (imode == 2) write(18, 1768) real(n, dp) * dr, &
                         real(upn, dp) * updr, aimag(ftv)
                    if (imode == 3) write(18, 1768) real(n, dp) * dr, &
                         real(upn, dp) * updr, abs(ftv) * abs(ftv)
                    if (imode == 0) write(18, 1768) real(n, dp) * dr, &
                         real(upn, dp) * updr, abs(ftv)
                 else              !plot out all the vector
                    if (imode == 1) write(18, 1768) real(n, dp) * dr, &
                         real(upn, dp) * updr, &
                         (real(ft(i), dp), i = 1, ndim)
                    if (imode == 2) write(18, 1768) real(n, dp) * dr, &
                         real(upn, dp) * updr, &
                         (aimag(ft(i)), i = 1, ndim)
                    if (imode == 3) write(18, 1768) real(n, dp) * dr, &
                         real(upn, dp) * updr, &
                         (abs(ft(i)) * abs(ft(i)), i = 1, ndim)
                    if (imode == 0) write(18, 1768) real(n, dp) * dr, &
                         real(upn, dp) * updr, &
                         (abs(ft(i)), i = 1, ndim)
                 endif
              else if (ndim == 9) then  
                 if (sproj ==1) then
                    ftv = zzero  
                    do j = 0, 2  
                       ftv = ftv + (bra(1) * ft(j * 3 + 1) + &
                            bra(2) * ft(j * 3 + 2) + bra(3) *ft(j * 3 + 3)) * &
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
                 endif
              else  
                 ftv = ft(1)  
                 if (imode == 1) write(18, 1768) real(n, dp) * dr, real(ftv, dp)
                 if (imode == 2) write(18, 1768) real(n, dp) * dr, aimag(ftv)
                 if (imode == 3) write(18, 1768) real(n, dp) * dr, &
                      abs(ftv) * abs(ftv)
                 if (imode == 0) write(18, 1768) real(n, dp) * dr, abs(ftv)  
              end if
           end if
        end do    !loop over line
        write(18,*) 
       enddo      !loop pver slice

     if (gs%myproc == 0) close(18)  
     write(9, 20) bra, ket, slice_plot(1:7, islice)  
     call myflush(9)  
  end do

10 format(' WRITING SLICEPLOT #',i4,' TO FILE: ',a)  

20 format(' BRA:',3f6.2,' KET:',3f6.2,/' PATH: ', 11f8.3)  

end subroutine sliceplot
