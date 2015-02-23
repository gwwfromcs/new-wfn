!-*-Fortran-*-
!
!     contains timings for the FFT routine
!
  real(dp) :: &
       tunfold, txfft, tyfft, tzfft, &
       t0rspace, ttrods, t0t1buf, tlt1, tt1, tlt2, tt2
  common /comtime/ tunfold, txfft, tyfft, tzfft, &
       t0rspace, ttrods,  t0t1buf, tlt1, tt1, tlt2, tt2 
