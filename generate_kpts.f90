! Originally generate_kpts.f90
!
! generates a file KPOINTS with the points for the bandstructure */
!
! Converted to Fortran90 by Peter Haynes March 2000

integer function generate_kpts(pw_params,nband_alt)

  use constants
  use pw_parameter_module
  use esdf
  implicit none

  type(pw_parameter), intent(inout) :: pw_params
  integer nband_alt

  double precision :: kpoint(7), kp(3)
  integer, parameter :: kptfile = 13
  integer :: i, j, k, l, ierr, nrk, n, nitems

  integer :: nlines, iline, ios
  character(len=llength) :: buf

  open(unit = kptfile, file = 'KPOINTS', iostat = ios, form = 'formatted', &
       status = 'replace')

  if (esdf_block('bandstructure', nlines)) then
     nitems = 0 ; nrk = 0
     do iline = 1, nlines
        if (index(block_data(iline), 'kpoint') > 0) then
           buf = block_data(iline)(index(block_data(iline), 'kpoint') + 6:)
           read(buf, *) kpoint(1:7)
           n = int(kpoint(7))
           if (n <= 0) then
              write(*, '(A/A)') 'number of bins in error:', buf
              call mystop
           end if
           nitems = nitems + 1
           nrk = nrk + n + 1
        else if (index(block_data(iline), 'label') <= 0) then
           write(*, '(A/A)') &
                'found funny keyword in bandstructure structure:', &
                block_data(iline)
        end if
     end do
  else
     call mystop( 'cannot find ''begin bandstructure'' in input file' )
  end if

  if (nitems == 0) call mystop( 'found no kpoints in bandstructure structure' )

  write(kptfile, '(I4)') nrk
  allocate(pw_params%bslabels(nitems * 100))
  pw_params%bslabels = repeat(' ', 100 * nitems)

  if (esdf_block('bandstructure', nlines)) then
     nitems = 0 ; nrk = 0
     do iline = 1, nlines
        if (index(block_data(iline), 'kpoint') > 0) then
           buf = block_data(iline)(index(block_data(iline), 'kpoint') + 6:)
           read(buf, *) kpoint(1:7)
           n = int(kpoint(7))
           if (n <= 0) then
              write(*, '(A/A)') 'number of bins in error:', buf
              call mystop
           end if
           do i = 1, n + 1
              do j = 1, 3
                 kp(j) = (i - 1) * (kpoint(j + 3) - kpoint(j)) / n
              end do
              l = 1
              if (i == n + 1) l = 0
              write(kptfile, '(F8.3,X,3(I4,X),F6.2,3(X,F15.8))') &
                   real(l, dp), nband_alt, 1, nband_alt, done, &
                   kp(1) + kpoint(1), kp(2) + kpoint(2), &
                   kp(3) + kpoint(3)
              nrk = nrk + 1
           end do
           nitems = nitems + 1
        else if (index(block_data(iline), 'label') > 0) then
           buf = adjustl(block_data(iline) &
                (index(block_data(iline), 'label') + 5:))
!           pw_params%bslabels(nitems * 100 + 1:nitems * 100 + 1 + len_trim(buf)) = &
!                trim(buf)
            do i = 1, len_trim(buf)
               pw_params%bslabels(nitems*100 + i) = buf(i:i)
            end do
        else
           write(*, '(A/A)') &
                'found funny keyword in bandstructure structure:', &
                block_data(iline)
        end if
     end do
  else
     call mystop( 'cannot find ''begin bandstructure'' in input file' )
  end if

  close(unit = kptfile)

  generate_kpts = nitems

end function generate_kpts
